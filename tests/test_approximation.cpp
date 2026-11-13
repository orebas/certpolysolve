#include "solve/policies.hpp"
#include "numeric/arithmetic/polynomial_arithmetic.hpp"
#include "exact_winding_api.hpp"
#include "api_test_support.hpp"
#include "internal.hpp"
#include <iostream>

using namespace reference;
namespace detail = certpolysolve::detail;

// An internal solver's answer, published as the header's own Result (the C entry points do the same).
static certpolysolve::Result published(cps_result answer) {
    return certpolysolve_test::published(std::move(answer));
}

static bool contains(const certpolysolve::Root &r, const Complex &z) {
    return endpoint(r.box.real.lo) <= z.real && z.real <= endpoint(r.box.real.hi) &&
           endpoint(r.box.imag.lo) <= z.imag && z.imag <= endpoint(r.box.imag.hi);
}

static void verify_geometry(const certpolysolve::Root &r, int bits,
                            const Rational &absolute_tolerance) {
    const Rational lo = endpoint(r.box.real.lo), hi = endpoint(r.box.real.hi);
    const Rational ilo = endpoint(r.box.imag.lo), ihi = endpoint(r.box.imag.hi);
    const Rational real_lower = lo > 0 ? lo : hi < 0 ? -hi : Rational(0);
    const Rational imag_lower = ilo > 0 ? ilo : ihi < 0 ? -ihi : Rational(0);
    const Rational error_squared = ((hi - lo) * (hi - lo) + (ihi - ilo) * (ihi - ilo)) / 4;
    require(error_squared <= absolute_tolerance * absolute_tolerance ||
                error_squared * power2(2 * bits) <= real_lower * real_lower + imag_lower * imag_lower,
            "independent squared-distance approximation check failed");
}

// Exact known roots give independent counts, including multiplicity, without
// calling the production certificate or another numerical root solver.
static void verify(const certpolysolve::Result &result, const std::vector<Complex> &known,
                   int bits, const Rational &absolute_tolerance = 0) {
    require(result.complete() && result.goal_met(), "approximation goal did not complete");
    const auto regions = result.roots();
    check_disjoint(regions);
    size_t total = 0;
    for (const auto &r : regions) {
        const auto count = std::count_if(known.begin(), known.end(), [&](const auto &z) { return contains(r, z); });
        require(static_cast<uint64_t>(count) == r.count, "independent approximation count disagrees");
        total += r.count;
        verify_geometry(r, bits, absolute_tolerance);
        if (r.multiplicity) {
            const auto first = std::find_if(known.begin(), known.end(), [&](const auto &z) { return contains(r, z); });
            require(first != known.end(), "claimed root absent");
            require(static_cast<uint64_t>(std::count_if(known.begin(), known.end(), [&](const auto &z) {
                        return z.real == first->real && z.imag == first->imag;
                    })) == *r.multiplicity, "cluster count was relabeled as multiplicity");
        }
    }
    require(total == known.size(), "approximation lost a root");
}

static std::vector<Complex> close_roots() {
    return {{1 - power2(-200)}, {1 + power2(-200)}, {4}};
}

static void newton_point_fast_path() {
    std::vector<Complex> roots{{0}, {0}};
    for (int j = 0; j < 20; ++j)
        roots.emplace_back(Rational(3*j-29)/8, Rational((j*j)%17+1)/8);
    auto polynomial = input(from_roots(roots));
    auto policy = detail::production_policy();
    policy.candidates = detail::CandidateEngine::Aberth;
    policy.study = detail::StudyMethod::Aberth;
    policy.refinement = detail::RefinementMethod::Aberth;
    policy.newton_point_refinement = true;
    policy.candidate_allow_underflow = true;
    detail::SolveState state(certpolysolve::detail::dense(*polynomial.get()), policy);
    cps_options options;
    options.goal = CPS_ROOT_APPROXIMATION;
    options.bits = 53;
    options.max_work_bits = 2048;
    auto answer = state.advance(&options);
    verify(published(answer), roots, 53);
    require(std::any_of(answer.roots.begin(), answer.roots.end(), [](const auto &r) {
        return r.certificate.find("Newton linear model") != std::string::npos;
    }), "new original-source Newton certificate was not exercised");
    // 96 bits, each a fresh solve (one advance a solve state since 18 October 2026): cancelled, it claims nothing;
    // otherwise it meets the goal.
    options.bits = 96;
    options.cancel = [](void *) -> int { return 1; };
    detail::SolveState cancelled(certpolysolve::detail::dense(*polynomial.get()), policy);
    require(!cancelled.advance(&options).goal_met, "a cancelled 96-bit Newton query claimed its goal");
    options.cancel = nullptr;
    detail::SolveState tighter(certpolysolve::detail::dense(*polynomial.get()), policy);
    verify(published(tighter.advance(&options)), roots, 96);
    roots = {{0}, {0}, {-3}, {-3}, {1-power2(-100)}, {1+power2(-100)}, {5}};
    polynomial = input(from_roots(roots));
    detail::SolveState clustered(certpolysolve::detail::dense(*polynomial.get()), policy);
    options.bits = 53;
    verify(published(clustered.advance(&options)), roots, 53);
}

static void classical_cancellation(bool deferred = false) {
    std::vector<Complex> roots{{0}, {0}};
    for (int j = 0; j < 24; ++j)
        roots.emplace_back(1 + (j + 1) * power2(-18), 1 + (j * j + 1) * power2(-20));
    auto polynomial = input(from_roots(roots));
    auto policy = detail::production_policy();
    // Exercise the ordinary packet owner; small Gaussian inputs otherwise
    // enter the existing local-component route after their machine packets.
    policy.candidates = detail::CandidateEngine::Aberth;
    policy.study = detail::StudyMethod::Aberth;
    policy.refinement = detail::RefinementMethod::Aberth;
    policy.classical_error_stopping = true;
    policy.defer_classical_validation = deferred;
    policy.study_statistics = std::make_shared<detail::StudyStatistics>();
    uint64_t pauses = 0;
    policy.trace_data = &pauses;
    policy.trace = [](void *data, const detail::TraceEvent &event) {
        if (std::string_view(event.stage) == "classical_arithmetic_stopping")
            *static_cast<uint64_t *>(data) += event.model_error_stops;
    };
    detail::SolveState state(certpolysolve::detail::dense(*polynomial.get()), policy);
    cps_options options;
    options.bits = 53;
    options.goal = CPS_ROOT_APPROXIMATION;
    options.max_work_bits = 2048;
    const auto result = state.advance(&options);
    require(pauses != 0, "expanded complex regression did not pause at an arithmetic floor");
    verify(published(result), roots, 53);
    if (deferred) require(policy.study_statistics->deferred_validations != 0,
                          "validation scheduling fixture did not exercise a deferred packet");
}

#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
// Through the research build's Persistent test engine (the default suite's restart coverage: test_refinement_state).
static void bounded_restart_proposals() {
    // Two unresolved tight clusters in a larger polynomial exercise local
    // restart proposals. The planted roots independently check every count
    // and the true-relative accuracy of the final original-input cover.
    std::vector<Complex> roots;
    for (int sign : {-1, 1})
        for (int j = 0; j < 4; ++j)
            roots.emplace_back(sign * power2(-16) + j * power2(-100));
    for (int re = -7; re <= 7; re += 2)
        for (int im = -7; im <= 7; im += 2) roots.emplace_back(re, im);
    auto polynomial = input(from_roots(roots));
    auto policy = detail::production_policy();
    policy.study = detail::StudyMethod::Aberth;
    policy.refinement = detail::RefinementMethod::Aberth;
    // The structured (component) route with the coefficient evaluator, where
    // the restarts happen. Until 7 October 2026 an explicit evaluation=
    // Coefficients override opened it; that private control is gone, and the
    // Persistent engine opens the same route with the same source.
    policy.candidates = detail::CandidateEngine::Persistent;
    policy.study_statistics = std::make_shared<detail::StudyStatistics>();
    detail::SolveState state(certpolysolve::detail::dense(*polynomial.get()), policy);
    cps_options options;
    options.bits = 53;
    options.goal = CPS_ROOT_APPROXIMATION;
    options.max_work_bits = 4096;
    const auto result = state.advance(&options);
    require(policy.study_statistics->candidate_restarts != 0 &&
                policy.study_statistics->candidates_replaced != 0,
            "cluster regression did not use local restart proposals");
    verify(published(result), roots, 53);
}
#endif

static void exact_clusters() {
    // Degree 5 (two far roots beside each cluster): Auto's closed forms answer cubics, separating these pairs at
    // 256 bits, and the general route's counted cluster is what this exercises.
    const auto with_far = [](std::vector<Complex> roots) {
        roots.emplace_back(Rational(-7));
        roots.emplace_back(Rational(5), Rational(3));
        return roots;
    };
    for (const auto &known : {with_far(close_roots()),
            with_far(std::vector<Complex>{{2 - power2(-220), 1}, {2 + power2(-220), 1}, {-3, 2}}),
            with_far(std::vector<Complex>{{power2(-300) - power2(-550)},
                                 {power2(-300) + power2(-550)}, {4}})}) {
        certpolysolve::Options options;
        options.goal = certpolysolve::Goal::approximation();
        options.bits = 53;
        options.max_bits = 256;
        auto result = certpolysolve::solve(input(from_roots(known)), options);
        if (!result.goal_met()) throw std::runtime_error(result.message());
        require(cps_result_goal(result.get()) == CPS_ROOT_APPROXIMATION, "wrong approximation goal tag");
        verify(result, known, 53);
        const auto regions = result.roots();
        require(std::any_of(regions.begin(), regions.end(),
                            [](const auto &r) { return r.count > r.multiplicity.value_or(0); }),
                "close roots did not exercise an unresolved counted cluster");
        // The individual approximation of a cluster is refused (the C call returns an error for it).
        bool refused = false;
        for (size_t i = 0; i < regions.size(); ++i) {
            double re = 0, im = 0;
            if (cps_result_get_approximation(result.get(), i, &re, &im) != CPS_OK) refused = true;
        }
        require(refused, "individual approximation conversion accepted a cluster");
    }
    const auto known = close_roots();
    certpolysolve::Options strict;
    strict.bits = 53;
    strict.max_bits = 128;
    const auto result = certpolysolve::solve(input(from_roots(known)), strict);
    require(cps_result_goal(result.get()) == CPS_ROOT_ACCURACY && !result.goal_met(),
            "default isolation accepted a cluster below its separation precision");
}

static void uncertain_clusters() {
    // (x-1)^2+t, |t|<=2^-160. The whole family has two nearby roots but
    // neither separation nor an exact double root holds uniformly.
    certpolysolve::Polynomial p;
    p.set(0, certpolysolve::ball(certpolysolve::Number("1"), certpolysolve::Number("1@-160")))
        .set(1, certpolysolve::Number("-2"))
        .set(2, certpolysolve::Number("1"));
    certpolysolve::Options options;
    options.goal = certpolysolve::Goal::approximation();
    options.bits = 53;
    options.max_bits = 256;
    auto result = certpolysolve::solve(p, options);
    for (const auto &known : {std::vector<Complex>{{1 - power2(-80)}, {1 + power2(-80)}},
                             std::vector<Complex>{{1}, {1}},
                             std::vector<Complex>{{1, -power2(-80)}, {1, power2(-80)}}})
        verify(result, known, 53);
    require(result.roots().size() == 1 && result.roots()[0].multiplicity.value_or(0) == 0,
            "uncertain collision invented a fixed multiplicity");
    options.goal = certpolysolve::Goal::accuracy();
    require(!certpolysolve::solve(p, options).goal_met(), "uncertainty yielded uniform individual isolation");

    // A ball around zero cannot imply any fixed true-relative accuracy.
    p.set(0, certpolysolve::ball(certpolysolve::Number("0"), certpolysolve::Number("1@-160"))).set(1, certpolysolve::Number("0"));
    options.goal = certpolysolve::Goal::approximation();
    options.work_budget = 2000;
    options.cell_budget = CPS_DEFAULT_LIVE_CELLS;
    require(!certpolysolve::solve(p, options).goal_met(), "approximation silently used an absolute floor");
    options.absolute_tolerance = certpolysolve::Number("1@-70");
    result = certpolysolve::solve(p, options);
    for (const auto &known : {std::vector<Complex>{{-power2(-80)}, {power2(-80)}},
                             std::vector<Complex>{{0}, {0}},
                             std::vector<Complex>{{0, -power2(-80)}, {0, power2(-80)}}})
        verify(result, known, 53, power2(-70));
}

static void several_tight_clusters() {
    std::vector<Complex> known{{0}, {-4}};
    for (int center = 1; center <= 8; ++center) {
        known.push_back({center, power2(-200)});
        known.push_back({center, -power2(-200)});
    }
    certpolysolve::Options options;
    options.goal = certpolysolve::Goal::approximation();
    options.bits = 53;
    options.max_bits = 192;
    const auto result = certpolysolve::solve(input(from_roots(known)), options);
    verify(result, known, 53);
    const auto regions = result.roots();
    require(std::any_of(regions.begin(), regions.end(), [](const auto &region) {
        return region.count > region.multiplicity.value_or(0);
    }), "tight conjugate pairs did not exercise counted approximation regions");
}

static void zeros_multiplicity_and_real() {
    std::vector<Complex> known{{0}, {0}, {1}, {1}, {4}};
    certpolysolve::Options options;
    options.goal = certpolysolve::Goal::approximation();
    options.bits = 80;
    options.max_bits = 512;
    auto result = certpolysolve::solve(input(from_roots(known)), options);
    verify(result, known, 80);
    const auto regions = result.roots();
    require(std::any_of(regions.begin(), regions.end(), [](const auto &r) {
        return r.multiplicity.value_or(0) == 2 && endpoint(r.box.real.lo) == 0 && endpoint(r.box.real.hi) == 0;
    }), "exact zero multiplicity lost");
    result = certpolysolve::real_roots(input(from_roots(known)), options);
    verify(result, known, 80);
    for (const auto &r : result.roots())
        require(r.realness == certpolysolve::Realness::real && r.count == r.multiplicity.value_or(0),
                "real query lost its individual proof");
}

// The approximation goal on close roots: a one-unit solve yields; 53, 128 and 300 bits are met; a cancelled solve
// does not claim success. (Until 18 October 2026 one session was tightened through these goals; the 0.5 API has
// one-shot solves only.)
static void budgets_and_cancellation() {
    const auto known = close_roots();
    const auto p = input(from_roots(known));
    certpolysolve::Options options;
    options.goal = certpolysolve::Goal::approximation();
    options.bits = 53;
    options.max_bits = 256;
    options.work_budget = 1;
    options.cell_budget = CPS_DEFAULT_LIVE_CELLS;
    // The final algebraic attempt can retain a counted cluster and meet
    // approximation within one unit, without resolving individual roots.
    const auto one_unit = certpolysolve::solve(p, options);
    verify(one_unit, known, 53);
    require(one_unit.statistics().work_units == 1 && one_unit.engine_counters().isolation_passes == 0,
            "one-unit algebraic approximation spent fallback work");
    const auto clustered = one_unit.roots();
    require(std::any_of(clustered.begin(), clustered.end(), [](const auto &root) {
                return root.count == 2 && !root.multiplicity;
            }), "one-unit approximation lost its counted cluster");
    // A degree-five input still exercises an actual general-route work stop.
    auto larger = known;
    larger.emplace_back(-7);
    larger.emplace_back(5, 3);
    const auto limited = certpolysolve::solve(input(from_roots(larger)), options);
    require(cps_result_termination(limited.get()) == CPS_WORK_LIMIT && !limited.goal_met() &&
                limited.statistics().work_units == 1, "general-route one-unit approximation did not yield");
    options.work_budget = 100000;
    verify(certpolysolve::solve(p, options), known, 53);
    options.bits = 128;
    options.max_bits = 512;
    verify(certpolysolve::solve(p, options), known, 128);
    options.bits = 300;
    options.max_bits = 2048;
    const auto result = certpolysolve::solve(p, options);
    verify(result, known, 300);
    require(cps_result_goal(result.get()) == CPS_ROOT_APPROXIMATION, "the approximation solve's goal");
    // Cancelled before the solve begins (the old always-true callback stopped it at its first poll).
    std::stop_source source;
    options.stop = source.get_token();
    source.request_stop();
    auto cancelled = certpolysolve::solve(p, options);
    require(cps_result_termination(cancelled.get()) == CPS_CANCELLED && !cancelled.goal_met(),
            "cancelled approximation claimed success");
}

static void native_and_degenerate_inputs() {
    const auto known = close_roots();
    Polynomial polynomial{1};
    for (const auto &root : known) polynomial = multiply(polynomial, {-root.real, 1});
    const std::vector<Rational> poles{0, 2, 6};
    std::vector<certpolysolve::SecularEquation::Term> terms;
    for (size_t i = 0; i < poles.size(); ++i) {
        Rational derivative = 1;
        for (size_t j = 0; j < poles.size(); ++j)
            if (i != j) derivative *= poles[i] - poles[j];
        const Rational residue = -evaluate(polynomial, poles[i]) / derivative;
        terms.push_back({.residue = certpolysolve::Number(residue.str()), .pole = certpolysolve::Number(poles[i].str())});
    }
    certpolysolve::Options options;
    options.goal = certpolysolve::Goal::approximation();
    options.bits = 53;
    options.max_bits = 512;
    const certpolysolve::SecularEquation secular(terms);
    verify(certpolysolve::solve(secular, options), known, 53);
    auto constant = certpolysolve::solve(certpolysolve::Polynomial{7}, options);
    require(cps_result_goal(constant.get()) == CPS_ROOT_APPROXIMATION && constant.goal_met() && constant.complete() &&
                constant.roots().size() == 0, "constant approximation contract failed");
    auto zero = certpolysolve::solve(certpolysolve::Polynomial{0}, options);
    require(cps_result_termination(zero.get()) == CPS_ZERO_POLYNOMIAL && !zero.goal_met(), "zero polynomial approximation succeeded");
    certpolysolve::Polynomial indeterminate;
    indeterminate.set(0, certpolysolve::Number("1")).set(1, certpolysolve::ball(certpolysolve::Number("0"), certpolysolve::Number("1")));
    auto degree = certpolysolve::solve(indeterminate, options);
    require(cps_result_termination(degree.get()) == CPS_INDETERMINATE_DEGREE && !degree.goal_met(),
            "uncertain degree approximation invented a finite root count");
}

static void euclidean_accuracy_boundary() {
    detail::Root root;
    root.count = 2;
    root.multiplicity = 0;
    arb_one(acb_realref(root.box.value()));
    arb_one(acb_imagref(root.box.value()));
    arb_add_error_2exp_si(acb_realref(root.box.value()), -3);
    arb_add_error_2exp_si(acb_imagref(root.box.value()), -3);
    root.publish(128);
    const detail::Number no_absolute;
    require(detail::approximation_met(root, 2, no_absolute, 128),
            "Euclidean diagonal enclosure was unnecessarily rejected");
    require(!detail::approximation_met(root, 4, no_absolute, 128),
            "Euclidean enclosure invented additional relative accuracy");
    require(detail::approximation_met(root, 4, detail::Number("3@-4"), 128) &&
                !detail::approximation_met(root, 4, detail::Number("1@-3"), 128),
            "Euclidean absolute tolerance boundary is incorrect");
    root.multiplicity = root.count;
    require(!detail::accuracy_met(root, 2, no_absolute, 128),
            "individual-root stopping rule changed with approximation geometry");
    acb_indeterminate(root.box.value());
    require(!detail::approximation_met(root, 2, detail::Number("1"), 128),
            "nonfinite enclosure passed the Euclidean check");
}

int main() {
    certpolysolve::ThreadScope lifetime;
    try {
        exact_clusters();
#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
        bounded_restart_proposals();
#endif
        classical_cancellation();
        classical_cancellation(true);
        newton_point_fast_path();
        uncertain_clusters();
        several_tight_clusters();
        zeros_multiplicity_and_real();
        budgets_and_cancellation();
        native_and_degenerate_inputs();
        euclidean_accuracy_boundary();
        std::cout << "Certified relative approximation: independent counts, accuracy and uncertainty passed\n";
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
