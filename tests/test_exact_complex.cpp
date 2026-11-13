#include "exact_winding_api.hpp"
#include "api_test_support.hpp"
#include "numeric/exact/exact_complex.hpp"
#include "numeric/exact/modular_squarefree.hpp"
#include <flint/ulong_extras.h>
#include <future>
#include <iostream>
#include <map>
#include <random>
#include <string>
#include <vector>

using namespace reference;
namespace detail = certpolysolve::detail;

static std::vector<detail::Coefficient> coefficients(const ComplexPolynomial &p) {
    std::vector<detail::Coefficient> result(p.size());
    for (size_t i = 0; i < p.size(); ++i) {
        result[i].real = detail::Number(p[i].real.str());
        result[i].imag = detail::Number(p[i].imag.str());
    }
    return result;
}
static Rational reference_value(const detail::Rational &q) {
    return Rational(detail::integer_string(fmpq_numref(q.value()))) /
           Rational(detail::integer_string(fmpq_denref(q.value())));
}
static ComplexPolynomial reference_polynomial(const detail::GaussianPolynomial &p) {
    ComplexPolynomial result;
    for (const auto &z : p)
        result.emplace_back(reference_value(z.real), reference_value(z.imag));
    return result;
}
static void verify_reconstruction(const ComplexPolynomial &original,
                                  const detail::GaussianDecomposition &decomposition) {
    require(decomposition.available(), "Q(i) preprocessing unexpectedly unavailable");
    ComplexPolynomial product{{1}};
    for (const auto &factor : decomposition.factors) {
        require(factor.multiplicity > 0, "zero factor multiplicity");
        const auto polynomial = reference_polynomial(factor.polynomial);
        for (uint64_t i = 0; i < factor.multiplicity; ++i)
            product = multiply(product, polynomial);
    }
    require(product.size() == original.size(), "factor reconstruction degree mismatch");
    for (size_t i = 0; i < product.size(); ++i)
        require(zero(product[i] * original.back() - original[i] * product.back()),
                "independent exact factor reconstruction failed");
}
static void public_roots(const std::vector<Complex> &known, int bits = 80) {
    const auto polynomial = from_roots(known);
    certpolysolve::Options options;
    options.bits = static_cast<uint64_t>(bits);
    options.max_bits = std::max<uint64_t>(4096, static_cast<uint64_t>(bits) + 128);
    auto result = certpolysolve::solve(input(polynomial), options);
    if (!result.goal_met())
        throw std::runtime_error("exact complex goal unmet: " + result.message());
    const auto& regions = result.roots();
    check_disjoint(regions);
    size_t total = 0;
    for (const auto &r : regions) {
        size_t count = 0;
        std::optional<Complex> representative;
        for (const auto &z : known)
            if (endpoint(r.box.real.lo) <= z.real && z.real <= endpoint(r.box.real.hi) &&
                endpoint(r.box.imag.lo) <= z.imag && z.imag <= endpoint(r.box.imag.hi)) {
                ++count;
                if (representative)
                    require(zero(z - *representative), "distinct exact roots mislabeled as a multiple root");
                else
                    representative = z;
            }
        require(count == r.count && count == r.multiplicity.value_or(0), "independent exact multiplicity mismatch");
        check_accuracy(r, bits);
        total += count;
    }
    require(total == known.size(), "lost exact complex roots");
}
static void algebra_and_bad_primes() {
    // Force both supported residue fields to be inconclusive when testing the
    // characteristic-zero fallback. Failure at only one prime needs no exact GCD.
    const Rational prime = Rational("2147483647") * Rational("998244353");
    const std::vector<std::vector<Complex>> constructions{
        {{0, 1}, {0, 1}},
        {{0, 1}, {0, 1}, {0, -1}, {0, -1}, {0, -1}},
        {{0}, {0}, {0}, {Rational(1), Rational(2)}},
    };
    cps_options options;
    for (const auto &roots : constructions) {
        auto p = from_roots(roots);
        for (auto &z : p)
            z = z * Complex(Rational(-17) / 3, Rational(5) / 7);
        const auto decomposition = detail::exact_complex_factors(coefficients(p), options);
        verify_reconstruction(p, decomposition);
    }
    // Distinct roots coincide modulo p. A modular failure cannot establish
    // multiplicity; the exact fallback must retain a squarefree degree-2 factor.
    auto p = from_roots({{0, 1}, {prime, 1}});
    auto decomposition = detail::exact_complex_factors(coefficients(p), options);
    verify_reconstruction(p, decomposition);
    require(decomposition.factors.size() == 1 && decomposition.factors[0].multiplicity == 1,
            "modular collision falsely inferred multiplicity");
    require(decomposition.status.find("reconstruction") != std::string::npos,
            "bad prime did not use exact fallback");
    p = from_roots({{Rational(1) / prime, 1}, {2, 3}});
    decomposition = detail::exact_complex_factors(coefficients(p), options);
    verify_reconstruction(p, decomposition);
    require(decomposition.status.find("reconstruction") != std::string::npos,
            "noninvertible denominator was accepted modulo p");
    for (auto &z : p)
        z = z * Complex(prime);
    verify_reconstruction(p, detail::exact_complex_factors(coefficients(p), options));
    p = from_roots({{Rational(1) / 3, Rational(2) / 7}, {2, 3}});
    decomposition = detail::exact_complex_factors(coefficients(p), options);
    verify_reconstruction(p, decomposition);
    require(decomposition.status.find("modulo") != std::string::npos,
            "simple input did not use modular certificate");
    // Checker negative control: reconstruction must actually reject corruption.
    fmpq_add_ui(decomposition.factors[0].polynomial[0].real.value(),
                decomposition.factors[0].polynomial[0].real.value(), 1);
    bool rejected = false;
    try {
        verify_reconstruction(p, decomposition);
    } catch (const std::runtime_error &) {
        rejected = true;
    }
    require(rejected, "independent reconstruction accepted corrupted factors");

    // A collision at the cheap split prime is recovered by the independent
    // extension field, including its existing denominator checks.
    p = from_roots({{0, 1}, {Rational("998244353"), 1}});
    decomposition = detail::exact_complex_factors(coefficients(p), options);
    verify_reconstruction(p, decomposition);
    require(decomposition.status.find("2147483647") != std::string::npos,
            "split-prime collision did not reach the extension-field proof");
    p = from_roots({{Rational(1)/Rational("998244353"), 1}, {2, 3}});
    decomposition = detail::exact_complex_factors(coefficients(p), options);
    verify_reconstruction(p, decomposition);
    require(decomposition.status.find("2147483647") != std::string::npos,
            "split-prime denominator was not refused");
    // Unlike the piecewise wrapper, the shared proof must retain repeated zeros.
    p = from_roots({{0}, {0}, {0,1}});
    require(detail::squarefree_modular(coefficients(p), options) == detail::ModularSquarefreeResult::Inconclusive,
            "shared modular proof silently removed repeated zero roots");
    // A nonzero Gaussian leading coefficient can vanish under i -> sqrt(-1).
    uint64_t root = 1, base = 3, exponent = (UINT64_C(998244353)-1)/4;
    for (; exponent; exponent >>= 1, base = base*base % 998244353)
        if (exponent & 1) root = root*base % 998244353;
    require(root*root % 998244353 == 998244352, "invalid test field embedding");
    p = {{1}, {0,1}, {-Rational(root), 1}};
    require(detail::squarefree_modular(coefficients(p), options) == detail::ModularSquarefreeResult::Inconclusive,
            "split-field degree loss was accepted as a square-free proof");
    verify_reconstruction(p, detail::exact_complex_factors(coefficients(p), options));
    int polls = 0;
    options.cancel_data = &polls;
    options.cancel = [](void *data) -> int { return ++*static_cast<int *>(data) == 3; };
    const auto stopped = detail::exact_complex_factors(coefficients(p), options);
    require(stopped.cancelled && !stopped.available(), "one-shot modular cancellation was lost");
}
static void limits() {
    auto input_coefficients = coefficients(from_roots({{0, 1}, {0, 1}}));
    for (uint64_t budget : {0u, 1u, 8u}) {
        cps_options options;
        options.exact_budget = budget;
        require(!detail::exact_complex_factors(input_coefficients, options).available(),
                "ignored exact preprocessing budget");
    }
    cps_options options;
    input_coefficients[0].real_radius = detail::Number("1/1000000");
    require(!detail::exact_complex_factors(input_coefficients, options).available(),
            "uncertain input got exact factors");
    input_coefficients[0].real_radius = detail::Number("0");
    input_coefficients[0].real = detail::Number("1@100000000000000000000");
    require(!detail::exact_complex_factors(input_coefficients, options).available(),
            "compressed exponent expanded outside budget");
    input_coefficients = coefficients(from_roots({{0, 1}, {0, 1}}));
    options.cancel = [](void *) { return 1; };
    require(detail::exact_complex_factors(input_coefficients, options).status.find("cancelled") !=
                std::string::npos,
            "ignored preprocessing cancellation");
    options.cancel = nullptr;
    input_coefficients.resize(258);
    input_coefficients.back().real = detail::Number("1");
    require(detail::exact_complex_factors(input_coefficients, options).available(),
            "former degree boundary rejected budgeted Q(i) input");
    options.exact_budget = 8;
    require(!detail::exact_complex_factors(input_coefficients, options).available(),
            "large Q(i) input ignored arithmetic budget");
    auto disabled = certpolysolve_test::raw_options();
    require(certpolysolve_test::set_exact_allowance(disabled.get(), 0) == CPS_OK, "setting the exact allowance failed");
    require(cps_options_set_precision(disabled.get(), 53) == CPS_OK && cps_options_set_max_work_bits(disabled.get(), 256) == CPS_OK,
            "setting the precision failed");
    auto partial = certpolysolve_test::solve(input(from_roots({{0, 1}, {0, 1}})), disabled.get());
    require(partial.complete() && !partial.goal_met(),
            "disabled preprocessing invented exact complex multiplicity");
    auto one_shot = certpolysolve_test::raw_options();
    require(cps_testing_options_stop_at_poll(one_shot.get(), 3) == CPS_OK, "setting the cancellation poll failed");
    // A quintic: Auto answers degree 3 and below by its closed forms in fewer polls than three, without exact
    // preparation.
    auto cancelled =
        certpolysolve_test::solve(input(from_roots({{0, 1}, {0, 1}, {1, 0}, {2, 0}, {-1, 0}})), one_shot.get());
    require(cps_result_termination(cancelled.get()) == CPS_CANCELLED, "one-shot Q(i) cancellation was lost");
}
static void large_modular_repeats() {
    // Exercise the backend's nontrivial polynomial GCD with a repeated
    // high-degree factor. Reconstruction alone would not detect a false
    // square-free claim, so also check the known factor exponents/degrees.
    ComplexPolynomial factor(130);
    factor[0] = Complex(1);
    factor[1] = Complex(0,1);
    factor[129] = Complex(1);
    const ComplexPolynomial other{{3,1}, {0}, {1}};
    const auto p = multiply(multiply(factor, factor), other);
    cps_options options;
    const auto decomposition = detail::exact_complex_factors(coefficients(p), options);
    verify_reconstruction(p, decomposition);
    require(decomposition.factors.size() == 2, "large modular gcd lost repeated factor");
    bool simple = false, repeated = false;
    for (const auto &f : decomposition.factors) {
        simple = simple || (f.multiplicity == 1 && f.polynomial.size() == 3);
        repeated = repeated || (f.multiplicity == 2 && f.polynomial.size() == 130);
    }
    require(simple && repeated, "large modular gcd fabricated square-freeness");
}
static void seeded_properties() {
    std::mt19937 generator(2026090802);
    cps_options options;
    for (int sample = 0; sample < 64; ++sample) {
        std::vector<Complex> known;
        const unsigned distinct = 1 + generator() % 4;
        for (unsigned i = 0; i < distinct; ++i) {
            const Complex z(Rational(static_cast<int>(generator() % 17) - 8) / (1 + generator() % 7),
                            Rational(1 + generator() % 11) / (1 + generator() % 5));
            const unsigned copies = 1 + generator() % 4;
            for (unsigned j = 0; j < copies; ++j)
                known.push_back(z);
        }
        auto p = from_roots(known);
        const auto decomposition = detail::exact_complex_factors(coefficients(p), options);
        verify_reconstruction(p, decomposition);
        public_roots(known);
    }
}
static void algebraic_repeats() {
    // (x^2-i)^3 has two non-rational complex roots, each of multiplicity three.
    ComplexPolynomial factor{{0, -1}, {0}, {1}};
    auto p = multiply(multiply(factor, factor), factor);
    auto result = certpolysolve::solve(input(p));
    require(result.goal_met() && result.roots().size() == 2, "non-rational repeated roots not isolated");
    for (const auto &r : result.roots()) {
        require(r.multiplicity.value_or(0) == 3 && r.count == 3, "non-rational multiplicity not proved");
        const auto count = rectangle_count(p, endpoint(r.box.real.lo), endpoint(r.box.real.hi),
                                           endpoint(r.box.imag.lo), endpoint(r.box.imag.hi));
        require(count && *count == 3, "independent winding count disagrees on algebraic repeats");
    }
}
// The decomposition of prod (x - z)^m over distinct nonzero z: one monic factor
// for each multiplicity m, equal to the product of (x - z) over the roots of
// that multiplicity. Reconstruction alone cannot see a repeated root left
// inside a factor labelled square-free.
static void verify_structure(const std::vector<std::pair<Complex, unsigned>> &roots,
                             const detail::GaussianDecomposition &decomposition) {
    std::map<unsigned, std::vector<Complex>> by_multiplicity;
    for (const auto &[z, multiplicity] : roots)
        by_multiplicity[multiplicity].push_back(z);
    require(decomposition.factors.size() == by_multiplicity.size(), "wrong number of square-free factors");
    for (const auto &factor : decomposition.factors) {
        const auto expected = by_multiplicity.find(static_cast<unsigned>(factor.multiplicity));
        require(expected != by_multiplicity.end(), "a factor has a multiplicity no root has");
        const auto wanted = from_roots(expected->second);
        const auto found = reference_polynomial(factor.polynomial);
        require(found.size() == wanted.size(), "a square-free factor has the wrong degree");
        for (size_t i = 0; i < found.size(); ++i)
            require(zero(found[i] * wanted.back() - wanted[i] * found.back()),
                    "a square-free factor is not the product of the roots of its multiplicity");
    }
}
static std::vector<Complex> expand(const std::vector<std::pair<Complex, unsigned>> &roots) {
    std::vector<Complex> all;
    for (const auto &[z, multiplicity] : roots)
        all.insert(all.end(), multiplicity, z);
    return all;
}
// The first primes the modular gcd uses: p = 1 (mod 4) above 2^61, in order.
static std::vector<Rational> gcd_primes(int count) {
    std::vector<Rational> primes;
    for (ulong p = UINT64_C(1) << 61; static_cast<int>(primes.size()) < count;) {
        p = n_nextprime(p, 1);
        if (p % 4 == 1) primes.emplace_back(std::to_string(p));
    }
    return primes;
}
// The gcd over Q(i) behind the decomposition is modular (split primes, the
// Chinese remainder theorem, rational reconstruction, exact division as the
// proof). Its answers are checked structurally, on inputs built to defeat
// the first primes it tries.
static void modular_gcd() {
    cps_options options;
    // Seeded: distinct nonzero Gaussian-rational roots, multiplicities 1 to 4.
    std::mt19937 generator(2026093001);
    for (int sample = 0; sample < 48; ++sample) {
        std::vector<std::pair<Complex, unsigned>> roots;
        const unsigned distinct = 2 + generator() % 9;
        while (roots.size() < distinct) {
            const Complex z(Rational(static_cast<int>(generator() % 2001) - 1000) / (1 + generator() % 97),
                            Rational(static_cast<int>(generator() % 2001) - 1000) / (1 + generator() % 89));
            bool fresh = !zero(z);
            for (const auto &[w, m] : roots) fresh = fresh && !zero(z - w);
            if (fresh) roots.emplace_back(z, 1 + generator() % 4);
        }
        // A sample with no repeated root would take the modular certificate.
        roots.front().second = std::max(roots.front().second, 2u);
        const auto p = from_roots(expand(roots));
        const auto decomposition = detail::exact_complex_factors(coefficients(p), options);
        verify_reconstruction(p, decomposition);
        verify_structure(roots, decomposition);
    }
    const auto primes = gcd_primes(3);
    const Rational p1 = primes[0], p2 = primes[1];
    // An unlucky prime: z and z + p1 coincide modulo p1, so gcd(f, f') has
    // degree 2 there and degree 1 in truth. The larger image must be dropped.
    {
        const Complex z(Rational(3) / 5, Rational(-7) / 11);
        const std::vector<std::pair<Complex, unsigned>> roots{{z, 2}, {z + Complex(p1), 1}, {Complex(2, 5), 1}};
        const auto p = from_roots(expand(roots));
        const auto decomposition = detail::exact_complex_factors(coefficients(p), options);
        verify_reconstruction(p, decomposition);
        verify_structure(roots, decomposition);
    }
    // Unlucky at the first two primes in a row.
    {
        const Complex z(Rational(-4) / 9, Rational(1) / 2);
        const std::vector<std::pair<Complex, unsigned>> roots{{z, 3}, {z + Complex(p1 * p2), 2}, {Complex(-1, 3), 1}};
        const auto p = from_roots(expand(roots));
        const auto decomposition = detail::exact_complex_factors(coefficients(p), options);
        verify_reconstruction(p, decomposition);
        verify_structure(roots, decomposition);
    }
    // A prime dividing a denominator, and one dividing the leading
    // coefficient's norm under one of its two maps, are skipped.
    {
        const std::vector<std::pair<Complex, unsigned>> roots{{Complex(Rational(1) / p1, Rational(2)), 2},
                                                              {Complex(Rational(5) / 3, Rational(1) / p2), 1}};
        auto p = from_roots(expand(roots));
        auto decomposition = detail::exact_complex_factors(coefficients(p), options);
        verify_reconstruction(p, decomposition);
        verify_structure(roots, decomposition);
        // Leading coefficient r + i with r^2 = -1 (mod p1): its image under
        // i -> -r is zero.
        ulong prime = UINT64_C(1) << 61;
        do prime = n_nextprime(prime, 1); while (prime % 4 != 1);
        const ulong unit = n_sqrtmod(prime - 1, prime);
        require(unit && n_mulmod2(unit, unit, prime) == prime - 1, "no square root of -1 for the test prime");
        const Complex leading(Rational(std::to_string(unit)), Rational(1));
        const std::vector<std::pair<Complex, unsigned>> simple{{Complex(1, 2), 2}, {Complex(-3, 1), 1}};
        p = from_roots(expand(simple));
        for (auto &c : p) c = c * leading;
        decomposition = detail::exact_complex_factors(coefficients(p), options);
        verify_reconstruction(p, decomposition);
        verify_structure(simple, decomposition);
    }
}
// Both modular square-free proofs are admitted by one rule, a unit of
// exact_operations an imported element and their residue polynomials within
// temporary_bytes, and a refused proof is a resource refusal, not evidence
// about the input: below the coefficient count no stage runs (the exact
// decomposition could not fit either), the extension field is refused on its
// own larger storage estimate, and a refusal reports no factor and no
// multiplicity (14 October 2026). The stages entered are
// read from the trace.
static std::vector<std::string> stages_entered;
static detail::GaussianDecomposition traced(const ComplexPolynomial &p, const cps_options &options) {
    detail::EnginePolicy policy;
    policy.trace = [](void *, const detail::TraceEvent &event) {
        // The four stages the admission decides between; the exact decomposition's own gcd stages
        // (gaussian_modular_gcd_images) are not at issue here.
        static const char *const stages[] = {"gaussian_expand", "gaussian_split_modular_gcd", "gaussian_modular_gcd",
                                             "gaussian_exact_decomposition"};
        if (std::string(event.phase) != "begin") return;
        for (const char *stage : stages)
            if (std::string(event.stage) == stage) stages_entered.emplace_back(stage);
    };
    stages_entered.clear();
    return detail::exact_complex_factors(coefficients(p), options, &policy);
}
static void modular_proof_admission() {
    using Stages = std::vector<std::string>;
    const Stages expand{"gaussian_expand"}, split{"gaussian_expand", "gaussian_split_modular_gcd"},
        both{"gaussian_expand", "gaussian_split_modular_gcd", "gaussian_modular_gcd"},
        all{"gaussian_expand", "gaussian_split_modular_gcd", "gaussian_modular_gcd", "gaussian_exact_decomposition"},
        split_then_exact{"gaussian_expand", "gaussian_split_modular_gcd", "gaussian_exact_decomposition"};
    const auto refused = [](const detail::GaussianDecomposition &d) {
        return !d.available() && d.factors.empty() && d.status.find("budget") != std::string::npos;
    };
    // x^2 + i: three coefficients, square-free at both primes. With a denominator the split prime divides it
    // is inconclusive there and proved in the extension field. (x + i)^2 is inconclusive at both.
    const ComplexPolynomial quadratic{{0, 1}, {0}, {1}};
    const ComplexPolynomial bad_split{{Rational(0), Rational(1) / Rational("998244353")}, {0}, {1}};
    const auto repeated = from_roots({{0, -1}, {0, -1}});
    cps_options options;
    options.exact_operations = 0;
    require(refused(traced(quadratic, options)) && stages_entered.empty(), "a zero allowance entered a stage");
    for (uint64_t allowance : {1, 2}) {
        options.exact_operations = allowance;
        for (const auto *p : {&quadratic, &bad_split, &repeated}) {
            const auto decomposition = traced(*p, options);
            require(refused(decomposition) && decomposition.status.find("refused") != std::string::npos,
                    "an allowance below the coefficient count was not refused as such");
            require(stages_entered == expand, "an allowance below the coefficient count entered a proof");
        }
    }
    options.exact_operations = 3;
    auto decomposition = traced(quadratic, options);
    require(decomposition.available() && decomposition.status.find("998244353") != std::string::npos &&
                stages_entered == split,
            "at the coefficient count the split prime did not prove x^2 + i alone");
    decomposition = traced(bad_split, options);
    require(decomposition.available() && decomposition.status.find("2147483647") != std::string::npos &&
                stages_entered == both,
            "a split prime inconclusive at a denominator was not followed by the extension field");
    decomposition = traced(repeated, options);
    require(refused(decomposition) && stages_entered == all,
            "a repeated root at the coefficient count was not an honest budget refusal after every stage");
    options.exact_operations = UINT64_MAX;
    decomposition = traced(repeated, options);
    verify_reconstruction(repeated, decomposition);
    verify_structure({{Complex(0, -1), 2}}, decomposition);
    require(stages_entered == all, "the unbounded decomposition skipped a stage");
    // Storage. The extension field's residue polynomials take more than the split prime's (whose estimate is
    // the entry guard's, three Gaussian rationals, or less): between the two the split prime is tried and is
    // inconclusive at its prime, the extension field is refused, and the exact decomposition decides, or names
    // the refusal when the allowance stops it.
    const uint64_t between = std::max<uint64_t>(3 * detail::split_proof_words * sizeof(ulong),
                                                3 * sizeof(detail::GaussianRational)),
                   extension_bytes = 3 * detail::extension_proof_words * sizeof(ulong);
    require(between < extension_bytes, "the extension field's estimate is not the larger");
    cps_options storage;
    storage.temporary_bytes = between;
    decomposition = traced(bad_split, storage);
    require(decomposition.available() && decomposition.status.find("reconstruction") != std::string::npos &&
                stages_entered == split_then_exact,
            "the extension field ran beyond its storage allowance, or its refusal was taken for evidence");
    storage.temporary_bytes = extension_bytes;
    decomposition = traced(bad_split, storage);
    require(decomposition.available() && decomposition.status.find("2147483647") != std::string::npos &&
                stages_entered == both,
            "the extension field was refused at its own estimate");
    storage.temporary_bytes = between;
    storage.exact_operations = 3;
    decomposition = traced(repeated, storage);
    require(refused(decomposition) && decomposition.status.find("temporary_bytes") != std::string::npos &&
                stages_entered == split_then_exact,
            "a modular proof refused for storage was not named in the budget refusal");
    // The shared proof refuses on its own below the coefficient count.
    options.exact_operations = 2;
    require(detail::squarefree_modular(coefficients(quadratic), options) == detail::ModularSquarefreeResult::Refused,
            "the shared proof ran below the coefficient count");
    options.exact_operations = 3;
    require(detail::squarefree_modular(coefficients(quadratic), options) == detail::ModularSquarefreeResult::Proved,
            "the shared proof was refused at the coefficient count");
}
// The exact budget counts word operations, so it bounds time as well as the
// number of rational operations (principled audit, fuzz finding F4).
static void exact_budget_bounds_the_work() {
    // 46 distinct Gaussian-rational roots with dyadic coordinates of up to 24
    // bits, one of them doubled: the shape of F4's degree-47 case, which took
    // 10.8 s in a Euclidean remainder sequence over Q(i). The modular gcd
    // needs about 2^18 word operations for it; the remainder sequence, charged
    // by the same unit, does not fit in eight times that.
    std::mt19937 generator(20260930);
    const auto coordinate = [&] {
        const unsigned bits = 1 + generator() % 24;
        const long range = 1L << bits;
        return Rational(static_cast<long>(generator() % (2 * range + 1)) - range) / Rational(range);
    };
    std::vector<std::pair<Complex, unsigned>> roots;
    while (roots.size() < 46) {
        const Complex z(coordinate(), coordinate());
        bool fresh = !zero(z);
        for (const auto &[w, m] : roots) fresh = fresh && !zero(z - w);
        if (fresh) roots.emplace_back(z, 1);
    }
    roots.front().second = 2;
    const auto p = from_roots(expand(roots));
    cps_options options;
    options.exact_operations = UINT64_C(1) << 21;
    const auto decomposition = detail::exact_complex_factors(coefficients(p), options);
    require(decomposition.available(), "a degree-47 double root did not decompose within 2^21 word operations");
    verify_reconstruction(p, decomposition);
    verify_structure(roots, decomposition);
    // The charge grows with the operands. The same three-root structure with
    // coordinates of about 40,000 bits in place of 4 must not fit a budget
    // that the small one fits with room to spare.
    const auto triple = [&](const Rational &scale) {
        const std::vector<std::pair<Complex, unsigned>> small{{Complex(scale * 3, scale * 5), 2},
                                                              {Complex(scale * -7, scale * 2), 1}};
        return from_roots(expand(small));
    };
    Rational large(1);
    for (int i = 0; i < 40000; ++i) large *= 2;
    options.exact_operations = UINT64_C(1) << 12;
    require(detail::exact_complex_factors(coefficients(triple(Rational(1))), options).available(),
            "a small double root did not fit 2^12 word operations");
    const auto refused = detail::exact_complex_factors(coefficients(triple(large + 1)), options);
    require(!refused.available() && refused.status.find("budget") != std::string::npos,
            "40,000-bit operands were charged like 4-bit ones");
}
// Exact zero roots are split off before the modular certificates. Otherwise
// x^m (m > 1) makes gcd(p, p') nontrivial, both certificates fail by
// construction and the exact Euclidean fallback runs at full degree: complex
// copies of Runnels11 (x^512 times degree 853) spent minutes there.
static void zero_roots_before_modular_proof() {
    std::vector<Complex> roots;
    for (int k = 1; k <= 24; ++k) roots.push_back({Rational(k, 7), Rational(3 - k % 5, 4)});
    cps_options options;
    for (int m : {1, 9}) {
        auto p = from_roots(roots);
        p.insert(p.begin(), static_cast<size_t>(m), Complex{});
        const auto decomposition = detail::exact_complex_factors(coefficients(p), options);
        verify_reconstruction(p, decomposition);
        require(decomposition.status.find("proved modulo") != std::string::npos,
                "exact zero roots defeated the modular squarefree certificate");
        bool zero_factor = false;
        for (const auto &f : decomposition.factors)
            zero_factor = zero_factor || (f.multiplicity == static_cast<uint64_t>(m) && f.polynomial.size() == 2 &&
                                          f.polynomial[0].zero() && reference_value(f.polynomial[1].real) == 1 &&
                                          reference_value(f.polynomial[1].imag) == 0);
        require(zero_factor, "exact zero-root factor x^m missing");
    }
    ComplexPolynomial monomial(6);
    monomial.push_back({Rational(3), Rational(4)});
    const auto decomposition = detail::exact_complex_factors(coefficients(monomial), options);
    require(decomposition.factors.size() == 1 && decomposition.factors[0].multiplicity == 6,
            "a complex monomial should have exactly the zero root");
}

int main() {
    certpolysolve::ThreadScope scope;
    try {
        zero_roots_before_modular_proof();
        algebra_and_bad_primes();
        modular_gcd();
        modular_proof_admission();
        exact_budget_bounds_the_work();
        large_modular_repeats();
        limits();
        seeded_properties();
        algebraic_repeats();
        for (int bits : {53, 128, 512, 4096})
            public_roots({{0, 1}, {0, 1}}, bits);
        std::vector<std::future<void>> tasks;
        for (int i = 0; i < 8; ++i)
            tasks.push_back(std::async(std::launch::async, [] {
                certpolysolve::ThreadScope worker_scope;
                public_roots(
                    {{Rational(1) / 3, Rational(2) / 7}, {Rational(1) / 3, Rational(2) / 7}, {2, -1}});
            }));
        for (auto &task : tasks)
            task.get();
        std::cout << "Q(i) reconstruction, bad-prime, budget, multiplicity, winding and concurrency checks "
                     "passed\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
