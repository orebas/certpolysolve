#include "numeric/arithmetic/polynomial_arithmetic.hpp"
#include <certpolysolve/certpolysolve.hpp>
#include "api_test_support.hpp"
#include "exact_winding_api.hpp"
#include "internal.hpp"
#include <charconv>
#include <chrono>
#include <iostream>
#include <random>
#include <sstream>

// A replayable property campaign, not a timing benchmark. Intentionally uses a
// slow independent exact checker and permits legitimate unmet solver goals.
namespace {
using namespace reference;
using Clock = std::chrono::steady_clock;
struct Configuration {
    uint64_t seed = 2026090601, start = 0, cases = 100, seconds = 0, work_units = 100000, portfolio = 0;
    // Explicit methods require exact coefficients and a complex approximation
    // goal: ball families become random-coefficient cases, cluster goals are off.
    cps_method method = CPS_METHOD_AUTO;
    // Also solve exact transforms of every planted-root case (see run_case).
    bool transforms = false;
    // Planted-root degrees are drawn from 1..max_degree. Values of 64 and above
    // reach the explicit methods' centering (heuristics::centering_min_degree).
    uint64_t max_degree = 12;
};
// Exact transform constants, as in tests/test_equivalence.cpp: 64-bit dyadic
// stand-ins for pi and sqrt(2), and a unit at an irrational angle.
const Rational pi64 = Rational(Integer("57952155664616982739")) / power2(64);
const Rational sqrt2_64 = Rational(Integer("26087635650665564424")) / power2(64);
const Complex unit{Rational(3, 5), Rational(4, 5)};
Complex conjugate(const Complex& z) { return {z.real, -z.imag}; }
struct Statistics {
    uint64_t cases = 0, regions = 0, known_counts = 0, winding_counts = 0;
    uint64_t boundary_checks_skipped = 0, goals_met = 0, local_regions = 0;
    uint64_t largest_outer_bits = 0, transformed = 0, unfinished = 0;
};
uint64_t integer(const char *text) {
    const std::string value(text);
    uint64_t result = 0;
    const auto parsed = std::from_chars(value.data(), value.data() + value.size(), result);
    require(parsed.ec == std::errc{} && parsed.ptr == value.data() + value.size(), "expected unsigned integer option");
    return result;
}
uint64_t mix(uint64_t x) {
    x += UINT64_C(0x9e3779b97f4a7c15);
    x = (x ^ (x >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    x = (x ^ (x >> 27)) * UINT64_C(0x94d049bb133111eb);
    return x ^ (x >> 31);
}
int signed_small(std::mt19937_64& generator, int magnitude = 16) {
    return static_cast<int>(generator() % static_cast<uint64_t>(2 * magnitude + 1)) - magnitude;
}
bool cancelled(void *data) { return Clock::now() >= *static_cast<Clock::time_point *>(data); }
int cancel_callback(void *data) { return cancelled(data) ? 1 : 0; }
// The 0.5 calls return codes where the old ones threw: an option the library refuses fails the case.
void accept(cps_error code, const char *what) { require(code == CPS_OK, what); }

struct Box {
    Rational left, right, bottom, top;
    bool contains(const Complex& z) const {
        return left <= z.real && z.real <= right && bottom <= z.imag && z.imag <= top;
    }
    Rational error() const { return (right - left + top - bottom) / 2; }
};
// Auto must always cover a fixed-degree finite polynomial. An explicit method
// may legitimately stop incomplete at the work-bit cap; such runs are counted,
// not checked.
void verify(const certpolysolve::Result& result, const ComplexPolynomial& polynomial,
            const std::vector<Complex> *known_roots, int bits, int cluster_bits,
            Statistics& statistics, bool allow_incomplete = false) {
    if (!result.complete()) {
        require(allow_incomplete, "a fixed-degree finite polynomial lost complete coverage");
        ++statistics.unfinished;
        return;
    }
    const auto roots = result.roots();
    check_disjoint(roots);
    uint64_t count = 0;
    for (const auto& root : roots) {
        const Box box{endpoint(root.box.real.lo), endpoint(root.box.real.hi), endpoint(root.box.imag.lo), endpoint(root.box.imag.hi)};
        require(root.count > 0, "empty region was published as a root enclosure");
        count += root.count;
        if (known_roots) {
            uint64_t exact_count = 0;
            std::optional<Complex> first;
            bool unique = true;
            for (const auto& z : *known_roots) if (box.contains(z)) {
                ++exact_count;
                if (first) unique = unique && z.real == first->real && z.imag == first->imag;
                else first = z;
                if (root.realness == certpolysolve::Realness::real) require(z.imag == 0, "false real-root classification");
                if (root.realness == certpolysolve::Realness::not_real) require(z.imag != 0, "false nonreal-root classification");
            }
            require(exact_count == root.count, "known roots disagree with an exported rectangle count");
            if (root.multiplicity) require(unique && root.multiplicity == exact_count, "distinct roots mislabeled as one multiple root");
            ++statistics.known_counts;
        } else {
            const auto independent_count = rectangle_count(polynomial, box.left, box.right, box.bottom, box.top);
            if (independent_count) {
                require(static_cast<uint64_t>(*independent_count) == root.count, "exact Sturm--Cauchy count disagrees with solver");
                ++statistics.winding_counts;
            } else ++statistics.boundary_checks_skipped;
        }
        if (result.goal_met()) {
            if (cps_result_goal(result.get()) == CPS_CLUSTER_RADIUS)
                require(box.error() <= power2(-cluster_bits), "cluster-radius goal is false after export");
            else if (cps_result_goal(result.get()) == CPS_ROOT_APPROXIMATION)
                // Counted clusters are allowed. max(|re|,|im|) >= |z|/sqrt(2),
                // so allow that factor, as tests/test_equivalence.cpp does.
                check_accuracy(root, bits - 1);
            else {
                require(root.multiplicity == root.count, "individual-root goal contains an unresolved cluster");
                check_accuracy(root, bits);
            }
        }
        ++statistics.regions;
        statistics.local_regions += root.certificate.find("pellet") == 0 ? 1u : 0u;
    }
    require(count == polynomial.size() - 1, "complete partition counts do not sum to polynomial degree");
}

bool run_case(const Configuration& configuration, uint64_t index, Clock::time_point& deadline, Statistics& statistics) {
    const uint64_t seed = configuration.seed, work_units = configuration.work_units, portfolio = configuration.portfolio;
    const bool explicit_method = configuration.method != CPS_METHOD_AUTO;
    std::mt19937_64 generator(mix(seed) ^ mix(index));
    const int kind = explicit_method && index % 10 == 9 ? 8 : static_cast<int>(index % 10);
    const int bits_choices[] = {16, 40, 53, 80, 128};
    const int bits = bits_choices[generator() % 5];
    const uint64_t maximum = UINT64_C(64) << (generator() % 4);
    const int cluster_bits = static_cast<int>(generator() % 65);
    // C options: the exact allowance alone and the campaign's deadline have no field of the 0.5 Options. The work
    // budget, which may be 0 units (no work), and the deadline, a callback on the whole run's clock, are the options'
    // internal fields: no public call sets either since stage 3b-3 (0 selects the automatic budget, and the
    // callback went).
    const auto options = certpolysolve_test::raw_options();
    accept(cps_options_set_precision(options.get(), static_cast<uint64_t>(bits)), "precision refused");
    accept(cps_options_set_max_work_bits(options.get(), maximum), "precision cap refused");
    accept(cps_options_set_max_iterations(options.get(), 0), "iteration budget refused");
    options.get()->work_units = work_units;
    options.get()->work_units_default = false;
    accept(cps_options_set_live_cells(options.get(), 65536), "live cells refused");
    accept(certpolysolve_test::set_exact_allowance(options.get(), generator() % 2 ? 0 : 32768), "exact allowance refused");
    const bool cluster_goal = generator() % 3 == 0 && !explicit_method;
    if (cluster_goal)
        accept(cps_options_set_cluster_radius(options.get(), ("1@" + std::to_string(-cluster_bits)).c_str()),
               "cluster radius refused");
    if (explicit_method) {
        accept(cps_options_set_method(options.get(), configuration.method), "method refused");
        accept(cps_options_set_goal(options.get(), CPS_ROOT_APPROXIMATION), "goal refused");
    }
    options.get()->cancel = cancel_callback;
    options.get()->cancel_data = &deadline;
    ComplexPolynomial polynomial;
    std::vector<Complex> known;
    Complex multiplier{1};
    std::vector<ComplexPolynomial> family_samples;
    certpolysolve::Polynomial input_polynomial;
    if (kind <= 6) {
        const int degree = 1 + static_cast<int>(generator() % configuration.max_degree);
        const bool real = generator() % 3 == 0;
        for (int j = 0; j < degree; ++j) {
            if (j && generator() % 4 == 0) known.push_back(known[generator() % known.size()]);
            else {
                const int scale = index % 17 == 0 ? signed_small(generator, 128) : 0;
                const Rational unit = power2(scale - static_cast<int>(generator() % 5));
                known.emplace_back(signed_small(generator) * unit, real ? Rational(0) : signed_small(generator) * unit);
            }
        }
        polynomial = from_roots(known);
        multiplier = {power2(signed_small(generator, 128)), 0};
        if (!real && generator() % 2) multiplier.imag = multiplier.real;
        for (auto& coefficient : polynomial) coefficient = coefficient * multiplier;
        input_polynomial = input(polynomial);
    } else if (kind <= 8) {
        const int degree = 2 + static_cast<int>(generator() % 7);
        polynomial.resize(static_cast<size_t>(degree) + 1);
        for (auto& coefficient : polynomial)
            coefficient = {Rational(signed_small(generator)) / (1 + generator() % 5),
                           Rational(signed_small(generator)) / (1 + generator() % 5)};
        polynomial.back() = {1};
        input_polynomial = input(polynomial);
    } else {
        // Independent coefficient balls around a cubic with a double root.
        // Four exact corners are checked separately; production certificates
        // must be uniform over the ENTIRE family, not just these samples.
        const Complex center{Rational(signed_small(generator)) / 4, Rational(signed_small(generator)) / 4};
        polynomial = from_roots({center, center, {center.real + 4, center.imag - 1}});
        const Rational uncertainty = power2(-16 - static_cast<int>(generator() % 65));
        input_polynomial = input(polynomial);
        input_polynomial.set(0, certpolysolve::ball(certpolysolve::Number::complex(certpolysolve::Number(polynomial[0].real.str()),
                                                                           certpolysolve::Number(polynomial[0].imag.str())),
                                                certpolysolve::Number(uncertainty.str()), certpolysolve::Number(uncertainty.str())));
        for (int real_sign : {-1, 1}) for (int imag_sign : {-1, 1}) {
            auto corner = polynomial;
            corner[0] = corner[0] + Complex(real_sign * uncertainty, imag_sign * uncertainty);
            family_samples.push_back(std::move(corner));
        }
    }
    try {
        require(input_polynomial.ok(), "the planted input recorded a bad value");
        const auto result = [&] {
            if (!portfolio) return certpolysolve_test::solve(input_polynomial, options.get());
#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
            certpolysolve::detail::EnginePolicy policy;policy.candidates=certpolysolve::detail::CandidateEngine::Portfolio;
            if (portfolio==2) {
                using M=certpolysolve::detail::IterationMethod;
                constexpr std::array methods{M::Aberth,M::Secular,M::DurandKerner,M::StructuredQR,M::StructuredQZ,M::JenkinsTraub,M::Laguerre,M::DenseEigen};
                policy.iteration=methods[(index/10)%methods.size()];
            }
            certpolysolve::detail::SolveState state(certpolysolve::detail::dense(*input_polynomial.get()),policy);
            return certpolysolve_test::published(state.advance(options.get()));
#else
            throw std::invalid_argument("--portfolio needs the research build");
#endif
        }();
        require(result.stop() != certpolysolve::Stop::error, "the solve reported an error");
        if (cps_result_termination(result.get()) == CPS_CANCELLED) return false;
        statistics.largest_outer_bits = std::max(statistics.largest_outer_bits, result.work_bits());
        try {
            if (family_samples.empty()) verify(result, polynomial, known.empty() ? nullptr : &known, bits, cluster_bits,
                                               statistics, explicit_method);
            else for (const auto& sample : family_samples) verify(result, sample, nullptr, bits, cluster_bits, statistics);
        } catch (...) {
            // Keep exact exported endpoints in a failed replay, including
            // endpoints that exceed the independent checker's storage budget.
            for (const auto& root : result.roots()) {
                std::cerr << "enclosure count=" << root.count << " certificate=" << root.certificate;
                for (const auto *x : {&root.box.real.lo, &root.box.real.hi, &root.box.imag.lo, &root.box.imag.hi})
                    std::cerr << ' ' << x->mantissa() << '@' << x->exponent();
                std::cerr << '\n';
            }
            throw;
        }
        statistics.goals_met += result.goal_met() ? 1u : 0u;
        // Exact transforms of a planted-root case: each variant's roots are an
        // exact image of the planted roots, and its count checks are exact.
        if (configuration.transforms && !known.empty()) {
            // equivariant: the transform preserves every root's relative
            // conditioning (a shift can move a root towards 0, where the same
            // relative accuracy costs more precision, so it is checked for
            // soundness only).
            struct Variant { const char* name; std::vector<Complex> roots; Complex scalar; bool equivariant; };
            std::vector<Variant> variants{{"rotation by (3+4i)/5", {}, {}, true}, {"shift by sqrt2", {}, {}, false},
                                          {"conjugation", {}, {}, true}, {"scalar pi", {}, {}, true},
                                          {"reversal", {}, {}, true}};
            const bool reversible = std::none_of(known.begin(), known.end(), [](const Complex& z) {
                return z.real == 0 && z.imag == 0;
            });
            for (const auto& z : known) {
                variants[0].roots.push_back(z * conjugate(unit)); // roots of p(unit x)
                variants[1].roots.push_back(z - Complex{sqrt2_64}); // roots of p(x + sqrt2)
                variants[2].roots.push_back(conjugate(z));
                variants[3].roots.push_back(z);
                // Roots of x^n p(1/x): 1/z = conj(z) / |z|^2.
                if (reversible) {
                    const Rational norm = z.real * z.real + z.imag * z.imag;
                    variants[4].roots.push_back({z.real / norm, -z.imag / norm});
                }
            }
            if (!reversible) variants.pop_back();
            for (auto& variant : variants) variant.scalar = multiplier;
            variants[3].scalar = multiplier * Complex{pi64};
            for (const auto& variant : variants) {
                auto transformed = from_roots(variant.roots);
                for (auto& coefficient : transformed) coefficient = coefficient * variant.scalar;
                const auto transformed_input = input(transformed);
                require(transformed_input.ok(), "the transformed input recorded a bad value");
                const auto again = certpolysolve_test::solve(transformed_input, options.get());
                require(again.stop() != certpolysolve::Stop::error, "the transformed solve reported an error");
                if (cps_result_termination(again.get()) == CPS_CANCELLED) return false;
                verify(again, transformed, &variant.roots, bits, cluster_bits, statistics, explicit_method);
                // Equivariance (principled audit, fix 19): an exact transform
                // changes nothing mathematical, so under the same budgets it
                // must meet a goal the original met, given headroom above the
                // goal. Without it, exactness decides: dyadic coefficients
                // reach 128 bits at a 128-bit budget, their rotation by
                // (3+4i)/5 (denominators 5^n) cannot. 64 guard bits, as the
                // solver's own arithmetic formulas use. An exact-size budget
                // is representation-dependent by nature (the rotation's
                // denominators grow the exact preparation), so a variant whose
                // preparation that budget limited is exempt.
                if (variant.equivariant && result.goal_met() && !again.goal_met() &&
                    maximum >= static_cast<uint64_t>(bits) + 64 && !again.preparation_limited())
                    throw std::runtime_error(std::string("metamorphic: ") + variant.name +
                                             " lost the goal the original met: " + again.message());
                statistics.goals_met += again.goal_met() ? 1u : 0u;
                ++statistics.transformed;
            }
        }
        ++statistics.cases;
    } catch (const std::exception& error) {
        std::cerr << "REPLAY: certpolysolve_stress --seed " << seed << " --start " << index
                  << " --cases 1 --portfolio " << portfolio << " --work-units " << work_units
                  << " --method " << configuration.method << " --transforms " << configuration.transforms
                  << " --max-degree " << configuration.max_degree << '\n'
                  << "kind=" << kind << " bits=" << bits << " max_outer=" << maximum
                  << " cluster_goal=" << cluster_goal << " cluster_bits=" << cluster_bits << '\n';
        for (size_t j = 0; j < polynomial.size(); ++j)
            std::cerr << "coefficient[" << j << "]=" << polynomial[j].real << "," << polynomial[j].imag << '\n';
        throw;
    }
    return true;
}
void report(const Configuration& configuration, const Statistics& statistics, Clock::time_point started, const char *label) {
    std::cout << label << " corpus=1 seed=" << configuration.seed << " start=" << configuration.start
              << " portfolio=" << configuration.portfolio << " work_units=" << configuration.work_units
              << " cases=" << statistics.cases << " next=" << configuration.start + statistics.cases
              << " regions_checked=" << statistics.regions << " known_counts=" << statistics.known_counts
              << " winding_counts=" << statistics.winding_counts << " boundary_inconclusive=" << statistics.boundary_checks_skipped
              << " goals_met=" << statistics.goals_met << " pellet_regions=" << statistics.local_regions
              << " max_outer_bits=" << statistics.largest_outer_bits << " method=" << configuration.method
              << " transformed=" << statistics.transformed << " unfinished=" << statistics.unfinished << " elapsed_s="
              << std::chrono::duration<double>(Clock::now() - started).count() << '\n' << std::flush;
}
}

int main(int argc, char **argv) {
    certpolysolve::ThreadScope thread_lifetime;
    try {
        Configuration configuration;
        for (int j = 1; j < argc; ++j) {
            const std::string option(argv[j]);
            require(j + 1 < argc, "usage: certpolysolve_stress [--seed N] [--start N] [--cases N] [--seconds N] "
                                  "[--work-units N] [--portfolio 0|1|2] [--method N] [--transforms 0|1] "
                                  "[--max-degree N]");
            const uint64_t value = integer(argv[++j]);
            if (option == "--seed") configuration.seed = value;
            else if (option == "--start") configuration.start = value;
            else if (option == "--cases") configuration.cases = value;
            else if (option == "--seconds") configuration.seconds = value;
            else if (option == "--work-units") configuration.work_units = value;
            else if (option == "--portfolio" && value<=2) configuration.portfolio = value;
            // cps_method values: 0 auto, 1 classical, 2 secular, 3 piecewise.
            else if (option == "--method" && value <= CPS_METHOD_PIECEWISE) configuration.method = static_cast<cps_method>(value);
            else if (option == "--transforms" && value <= 1) configuration.transforms = value == 1;
            else if (option == "--max-degree" && value >= 1 && value <= 512) configuration.max_degree = value;
            else throw std::invalid_argument("unknown stress option");
        }
        require(configuration.cases >= 1 && configuration.cases <= 10000000 && configuration.seconds <= 43200
                && configuration.start <= UINT64_MAX - configuration.cases, "stress limits out of range");
        require(configuration.method == CPS_METHOD_AUTO || configuration.portfolio == 0,
                "portfolio policies apply to Auto only");
        const auto started = Clock::now();
        auto deadline = configuration.seconds ? started + std::chrono::seconds(configuration.seconds) : Clock::time_point::max();
        auto next_report = started + std::chrono::seconds(30);
        Statistics statistics;
        report(configuration, statistics, started, "start");
        for (uint64_t offset = 0; offset < configuration.cases && Clock::now() < deadline; ++offset) {
            if (!run_case(configuration, configuration.start + offset, deadline, statistics)) break;
            if (Clock::now() >= next_report) {
                report(configuration, statistics, started, "progress");
                next_report = Clock::now() + std::chrono::seconds(30);
            }
        }
        report(configuration, statistics, started, "passed");
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "FAILED: " << error.what() << '\n';
        return 1;
    }
}
