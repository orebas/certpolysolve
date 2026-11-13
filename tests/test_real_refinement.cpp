#include "solve/real_query.hpp"
#include "solve/prepared_piecewise.hpp"
#include "solve/policies.hpp"
#include "numeric/arithmetic/polynomial_arithmetic.hpp"
#include "exact_reference_api.hpp"
#include "api_test_support.hpp"
#include "core/heuristics.hpp"
#include "internal.hpp"
#include "methods/piecewise/piecewise_roots.hpp"
#include "solve/real_only_route.hpp"
#include "numeric/evaluate/sparse_real_evaluator.hpp"
#include <flint/fmpq.h>
#include <iostream>
#include <sstream>

using namespace reference;
static void require(bool condition, const std::string &message) { reference::require(condition, message.c_str()); }
static void verify(const Polynomial &p, const certpolysolve::Result &r, int count, int bits,
                   const Rational &absolute_tolerance = 0) {
    if (!r.complete() || !r.goal_met())
        throw std::runtime_error("real refinement incomplete: " + r.message());
    check_disjoint(r.roots());
    int total = 0;
    for (const auto &root : r.roots()) {
        require(root.realness == certpolysolve::Realness::real, "realness missing");
        require(roots_closed(p, endpoint(root.box.real.lo), endpoint(root.box.real.hi)) ==
                    static_cast<int>(root.count),
                "independent real Sturm disagreement");
        total += static_cast<int>(root.count);
        if ((endpoint(root.box.real.hi) - endpoint(root.box.real.lo)) / 2 > absolute_tolerance)
            check_accuracy(root, bits);
    }
    require(total == count, "wrong total real count");
}
// An absolute tolerance on real roots: met at 53 bits; a one-unit solve stops honestly at its budget; a cancelled
// solve says so; and the tolerance cleared, 160 bits are met. (Until 18 October 2026 one session was resumed a
// unit at a time and then tightened; the 0.5 API has one-shot solves only.)
static void absolute_separation() {
    const std::vector<std::pair<Polynomial, std::string>> cases{
        {{-2, 0, 1}, "2"}, {{-2, 0, 1}, "1/1000"}, {{-1, 0, Rational("100000000000000")}, "1e-6"}};
    for (const auto &[p, tolerance] : cases) {
        certpolysolve::Options options;
        options.bits = 53;
        options.max_bits = 4096;
        options.absolute_tolerance = certpolysolve::Number(tolerance);
        options.work_budget = 100000;
        const Rational tolerance_value = tolerance == "1e-6" ? Rational(1) / 1000000 : Rational(tolerance);
        verify(p, certpolysolve::real_roots(input(p), options), 2, 53, tolerance_value);
        options.work_budget = 1;
        const auto unit = certpolysolve::real_roots(input(p), options);
        if (!unit.goal_met()) {
            require(cps_result_termination(unit.get()) == CPS_WORK_LIMIT, "absolute-tolerance separation stopped early");
            require(unit.statistics().work_units <= 1, "real separation exceeded its quota");
        }
        std::stop_source cancel;
        cancel.request_stop();
        options.stop = cancel.get_token();
        options.work_budget = 100000;
        require(cps_result_termination(certpolysolve::real_roots(input(p), options).get()) == CPS_CANCELLED,
                "an absolute-tolerance real solve ignored cancellation");
        options.stop = {};
        options.absolute_tolerance = certpolysolve::Number("0");
        options.bits = 160;
        options.max_bits = 4096;
        verify(p, certpolysolve::real_roots(input(p), options), 2, 160);
    }
}
// Real roots at 53 bits, a one-unit solve stopping honestly at its budget, and 256 bits under a larger cap.
static void real_solves(const Polynomial &p, int count) {
    certpolysolve::Options o;
    o.work_budget = 1;
    o.bits = 53;
    o.max_bits = 4096;
    const auto unit = certpolysolve::real_roots(input(p), o);
    if (!unit.goal_met()) {
        require(cps_result_termination(unit.get()) == CPS_WORK_LIMIT, "a one-unit real solve stopped for another reason");
        require(unit.statistics().work_units <= 1, "real refinement exceeded controller quota");
    }
    o.work_budget = 100000;
    verify(p, certpolysolve::real_roots(input(p), o), count, 53);
    o.bits = 256;
    o.max_bits = 8192;
    verify(p, certpolysolve::real_roots(input(p), o), count, 256);
}
// (3x - 1)(5x - 7)(7x + 17) prod_{i=1..pairs} (x^2 + 9^i): the real roots 1/3, 7/5 and -17/7 and the pairs +-3^i i
// over many scales. At 126 pairs (degree 255) the selector sends its complex query to the piecewise method, so a
// real-only query tries the axis charts first (session.cpp).
static Polynomial graded_pairs(int pairs) {
    auto p = multiply(multiply({-1, 3}, {-7, 5}), {17, 7});
    Rational scale = 1;
    for (int i = 1; i <= pairs; ++i) {
        scale *= 9;
        p = multiply(p, {scale, 0, 1});
    }
    return p;
}
// An answer holding exactly the known real roots (with their multiplicities), one in each region; by construction
// no other root is real. (Sturm on the degree-255 input with its 25,000-bit coefficients would take too long.)
static void verify_known(const certpolysolve::Result &r, const std::vector<std::pair<Rational, unsigned>> &roots,
                         int bits) {
    if (!r.complete() || !r.goal_met())
        throw std::runtime_error("real answer incomplete: " + r.message());
    check_disjoint(r.roots());
    require(r.roots().size() == roots.size(), "a known real root missing or an extra region");
    for (const auto &root : r.roots()) {
        require(root.realness == certpolysolve::Realness::real, "realness missing");
        const Rational a = endpoint(root.box.real.lo), b = endpoint(root.box.real.hi);
        int inside = 0;
        unsigned count = 0;
        for (const auto &[x, multiplicity] : roots)
            if (a <= x && x <= b) {
                ++inside;
                count = multiplicity;
            }
        require(inside == 1 && root.count == count, "a region without its known root");
        check_accuracy(root, bits);
    }
}
static const std::vector<std::pair<Rational, unsigned>> graded_real{
    {Rational(1) / 3, 1}, {Rational(7) / 5, 1}, {Rational(-17) / 7, 1}};

// Drive the private axis state directly: small public budgets deliberately skip
// this route. Retain an exact zero factor across cold and mid-round cancellation,
// then complete one work unit at a time and check the independently planted roots.
static void axis_interruptions(const Polynomial &nonzero) {
    namespace d = certpolysolve::detail;
    const auto polynomial = input(multiply(nonzero, {0, 0, 1}));
    auto config = d::maintained_piecewise_config();
    config.axis_only = true;
    config.axis_rounds = d::heuristics::axis_attempt_rounds;
    d::PiecewiseRootState state(d::dense(*polynomial.get()), config);
    cps_options query;
    // The real-only adapter supplies a complex query to the private axis state.
    query.goal = CPS_ROOT_APPROXIMATION;
    query.bits = 53;
    query.max_work_bits = 4096;
    query.work_units = 1;
    const uint64_t budget = d::heuristics::default_work_units_per_root * (nonzero.size() + 1);
    uint64_t charged = 0;
    const auto retained_zero = [&](const cps_result &r) {
        require(!r.complete && !r.goal_met && r.roots.size() == 1 && state.certified_regions() == 1,
                "an interrupted axis pass published tentative roots or lost its exact zero");
        const auto &zero = r.roots.front();
        require(zero.count == 2 && zero.multiplicity == 2 && zero.realness == CPS_REAL &&
                    acb_is_zero(zero.box.value()),
                "an interrupted axis pass changed its exact zero certificate");
    };
    const auto cancel = [&] {
        const auto before = state.statistics();
        query.cancel = [](void *) { return 1; };
        const auto stopped = state.advance(query);
        query.cancel = nullptr;
        retained_zero(stopped);
        require(stopped.termination == CPS_CANCELLED && stopped.progress.work_units_this_call == 0 &&
                    stopped.statistics_v2.work_units == charged,
                "an axis cancellation performed work or reported the wrong stop");
        const auto &after = state.statistics();
        require(after.work == before.work && after.rounds == before.rounds && after.models == before.models &&
                    after.fft_transforms == before.fft_transforms && after.coefficients == before.coefficients,
                "an axis cancellation advanced the retained model stream");
    };
    const auto advance = [&] {
        auto r = state.advance(query);
        require(r.progress.work_units_this_call <= 1, "an axis call exceeded its one-unit quota");
        charged += r.progress.work_units_this_call;
        require(r.statistics_v2.work_units == charged, "axis resume reset or duplicated charged work");
        if (!r.goal_met) {
            retained_zero(r);
            require(r.termination == CPS_WORK_LIMIT && r.progress.work_units_this_call == 1,
                    "the axis interruption fixture stopped before completing its query");
        }
        return r;
    };
    cancel();
    cps_result result;
    // Two charts per ring: reaching a second pair exercises retained charts
    // and a neighbour boundary, rather than cancelling only before construction.
    while (state.statistics().models < 4 && !result.goal_met && charged < budget) result = advance();
    require(state.statistics().models >= 4 && !result.goal_met,
            "the axis interruption fixture missed the active model stream");
    cancel();
    while (!result.goal_met && charged < budget) result = advance();
    require(result.goal_met, "the axis state did not finish within the default work budget");
    auto known = graded_real;
    known.emplace_back(0, 2);
    verify_known(certpolysolve_test::published(std::move(result)), known, 53);
    // A tighter query uses a fresh axis state, as the public adapter does.
    ++query.bits;
    bool rejected = false;
    try { (void)state.advance(query); }
    catch (const std::logic_error &) { rejected = true; }
    require(rejected, "an axis state accepted a different query");
}

// The axis charts' attempt (real_only_route.cpp) under the default work budget: a real-only query it answers, a
// 200-bit one, and the same query under a budget below the default's units a root,
// which goes to the complex route (an attempt could spend a budget the route needs: the fuzz campaign's tight
// budgets lost 15 completions so).
static void axis_attempt() {
    const auto p = graded_pairs(126);
    certpolysolve::Options options;
    options.bits = 53;
    options.max_bits = 4096;
    auto result = certpolysolve::real_roots(input(p), options);
    verify_known(result, graded_real, 53);
    require(result.message().find("real roots by the axis charts") != std::string::npos,
            "the axis charts did not answer a real-only query under the default budget");
    require(result.message().find("sparse recursion") == std::string::npos,
            "the sparse recursion ran on a dense input");
    // A 200-bit query is answered too, by whichever route the selector gives it. (Until 18 October 2026 a later
    // call refined the axis charts' answer in place; a solve state is advanced once since, and the real-only
    // route's refinements are fresh solves.)
    {
        certpolysolve::Options tighter;
        tighter.bits = 200;
        tighter.max_bits = 4096;
        verify_known(certpolysolve::real_roots(input(p), tighter), graded_real, 200);
    }
    // One unit below the attempt's threshold, the default's units a root.
    certpolysolve::Options budgeted;
    budgeted.bits = 53;
    budgeted.max_bits = 4096;
    budgeted.work_budget = certpolysolve::detail::heuristics::default_work_units_per_root * (p.size() - 1) - 1;
    const auto routed = certpolysolve::real_roots(input(p), budgeted);
    verify_known(routed, graded_real, 53);
    require(routed.message().find("axis charts") == std::string::npos &&
                routed.message().find("real roots by the complex route") != std::string::npos,
            "a real-only query under a tight budget did not go straight to the complex route");
}
// Where the selector chooses the global recipes for the complex query there is no attempt: in the first full check
// it cost 10% to 15% there whether it answered or not.
static void axis_not_attempted() {
    Polynomial far(21, 0); // x^20 + 2
    far[0] = 2;
    far[20] = 1;
    const auto p = multiply(multiply(multiply({-1, 1}, {-2, 1}), {3, 1}), far);
    certpolysolve::Options options;
    options.bits = 53;
    options.max_bits = 4096;
    const auto result = certpolysolve::real_roots(input(p), options);
    verify(p, result, 3, 53);
    require(result.message().find("axis charts") == std::string::npos &&
                result.message().find("real roots by the complex route") != std::string::npos,
            "an attempt ran where the selector chooses the global recipes");
}
// A factored input (a repeated or a power factor) ends the axis charts' attempt at once. The factor wrapper counts a
// factor's answer complete when it covers all of the factor's roots, which an axis answer (the real ones) never
// does, so the attempt raised the goal to the precision cap before it stopped (prepared_piecewise.cpp). Now its work
// is its preparation's unit: the same query just below the attempt's budget threshold, which goes straight to the
// complex route, spends one unit less.
static void axis_factored() {
    const auto p = multiply(graded_pairs(126), {-1, 3}); // 1/3 twice
    const std::vector<std::pair<Rational, unsigned>> real{
        {Rational(1) / 3, 2}, {Rational(7) / 5, 1}, {Rational(-17) / 7, 1}};
    const uint64_t threshold = certpolysolve::detail::heuristics::default_work_units_per_root * (p.size() - 1);
    const auto solve = [&](uint64_t budget) {
        certpolysolve::Options options;
        options.bits = 53;
        options.max_bits = 4096;
        options.work_budget = budget;
        auto result = certpolysolve::real_roots(input(p), options);
        verify_known(result, real, 53);
        return result;
    };
    const auto attempted = solve(threshold), direct = solve(threshold - 1);
    require(attempted.message().find("the axis charts' attempt stopped") != std::string::npos &&
                direct.message().find("axis charts") == std::string::npos,
            "the attempt's threshold moved");
    require(attempted.statistics().work_units <= direct.statistics().work_units + 1,
            "the attempt on a factored input spent more than its preparation");
}
// B6: the sparse recursion's attempt (session.cpp) on an exact real input with few terms. The answer's roots are
// checked by the input itself: each enclosure holds a sign change of the exact input (a certified evaluation at
// its dyadic ends), is real with count 1, and meets the accuracy; the count is known by construction.
static certpolysolve::detail::RealBall arb_of(const Rational &x) {
    certpolysolve::detail::RealBall out;
    fmpq_t q;
    fmpq_init(q);
    fmpq_set_str(q, x.str().c_str(), 10);
    arb_set_fmpq(out.value(), q, 4096);
    fmpq_clear(q);
    require(arb_is_exact(out.value()), "an endpoint is not a short dyadic");
    return out;
}
static void check_sparse_answer(const certpolysolve::Result &r, const std::vector<std::pair<slong, const char *>> &terms,
                                size_t count, int bits, const std::string &label) {
    using namespace certpolysolve::detail;
    require(r.complete() && r.goal_met(), label + ": " + r.message());
    require(r.message().find("real roots by the sparse recursion") != std::string::npos,
            label + ": the sparse recursion did not answer: " + r.message());
    check_disjoint(r.roots());
    require(r.roots().size() == count, label + ": " + std::to_string(r.roots().size()) + " roots");
    std::vector<SparseRealTerm> sparse;
    for (const auto &[e, c] : terms) sparse.push_back({e, Number(c)});
    const SparseRealEvaluator q(sparse);
    for (const auto &root : r.roots()) {
        require(root.realness == certpolysolve::Realness::real && root.count == 1, label + ": realness or count");
        check_accuracy(root, bits);
        RealBall value, slope;
        const auto lo = arb_of(endpoint(root.box.real.lo)), hi = arb_of(endpoint(root.box.real.hi));
        q.evaluate2(value.value(), slope.value(), lo.value(), 512);
        const int at_lo = arb_is_positive(value.value()) ? 1 : arb_is_negative(value.value()) ? -1 : 0;
        q.evaluate2(value.value(), slope.value(), hi.value(), 512);
        const int at_hi = arb_is_positive(value.value()) ? 1 : arb_is_negative(value.value()) ? -1 : 0;
        require(at_lo * at_hi < 0, label + ": an enclosure without a sign change of the input at its ends");
    }
}
static void sparse_attempt() {
    Polynomial p(4097, 0); // x^4096 - 2: the roots +-2^(1/4096)
    p[0] = -2;
    p[4096] = 1;
    const std::vector<std::pair<slong, const char *>> terms{{0, "-2"}, {4096, "1"}};
    certpolysolve::Options options;
    options.bits = 53;
    options.max_bits = 4096;
    auto result = certpolysolve::real_roots(input(p), options);
    check_sparse_answer(result, terms, 2, 53, "x^4096 - 2");
    options.bits = 200;
    options.max_bits = 4096;
    result = certpolysolve::real_roots(input(p), options);
    check_sparse_answer(result, terms, 2, 200, "x^4096 - 2 at 200 bits");
    require(result.work_bits() >= 200, "the 200-bit answer did not report its working precision");
    // One unit below the attempt's budget threshold: no attempt of either kind.
    certpolysolve::Options budgeted;
    budgeted.bits = 53;
    budgeted.max_bits = 4096;
    budgeted.work_budget = certpolysolve::detail::heuristics::default_work_units_per_root * 4096 - 1;
    const auto routed = certpolysolve::real_roots(input(p), budgeted);
    require(routed.message().find("sparse recursion") == std::string::npos &&
                routed.message().find("real roots by the complex route") != std::string::npos,
            "a tight budget did not go straight to the complex route");
}
// 1025 x^2049 - 12294 x^1025 + 18902025 x + 2100225 = 2100225 (x^2049/2049 - 6 x^1025/1025 + 9 x + 1): its
// derivative is 2100225 (x^1024 - 3)^2, whose double root is the next level's root, so the recursion stops
// (a possible common root of consecutive derivatives) and the complex route answers: p is increasing, with one
// real root (negative).
static void sparse_stop() {
    Polynomial p(2050, 0);
    p[0] = 2100225;
    p[1] = 18902025;
    p[1025] = -12294;
    p[2049] = 1025;
    certpolysolve::Options options;
    options.bits = 53;
    options.max_bits = 4096;
    const auto result = certpolysolve::real_roots(input(p), options);
    require(result.complete() && result.goal_met() && result.roots().size() == 1 &&
                endpoint(result.roots()[0].box.real.hi) < 0,
            "the planted input's one negative root: " + result.message());
    require(result.message().find("the sparse recursion stopped") != std::string::npos &&
                result.message().find("real roots by the complex route") != std::string::npos,
            "the stop did not hand over: " + result.message());
}
// The derivative of this four-term polynomial is n(m+1)(x^m - 3)^2.
// It is strictly increasing with one real root, but the repeated derivative
// root stops sparse isolation. At this degree the axis charts then succeed.
static void sparse_axis_accounting() {
    namespace d = certpolysolve::detail;
    constexpr uint64_t m = 4096, n = 2 * m + 1;
    Polynomial p(n + 1, 0);
    p[0] = n * (m + 1); p[1] = 9 * n * (m + 1);
    p[m + 1] = -6 * Rational(n); p[n] = m + 1;
    const auto source = input(p);
    certpolysolve::Options options;
    options.goal = certpolysolve::Goal::accuracy();
    options.bits = 53;
    options.max_bits = 256;
    options.work_budget = d::heuristics::default_work_units_per_root * n;
    certpolysolve_test::COptions raw(options, true);
    auto query = d::real_query::complex_options(*raw.get());
    d::SparseRealAttempt sparse(*source.get());
    const auto first = sparse.advance(query);
    require(!first.goal_met && first.progress.work_units_this_call > 0 &&
                first.progress.work_units_this_call < query.work_units,
            "accounting fixture did not stop a nonempty sparse attempt");
    query.work_units -= first.progress.work_units_this_call;
    auto config = d::maintained_piecewise_config();
    config.axis_only = true;
    config.axis_rounds = d::heuristics::axis_attempt_rounds;
    d::PreparedPiecewiseState axis(d::dense(*source.get()), config);
    const auto second = axis.advance(query);
    require(second.goal_met, "accounting fixture did not finish through the axis charts");
    const auto result = certpolysolve::real_roots(source, options);
    require(result.goal_met() && result.complete() &&
                result.message().find("the sparse recursion stopped") != std::string::npos &&
                result.message().find("real roots by the axis charts") != std::string::npos,
            "public accounting fixture did not take both real attempts");
    require(result.statistics().work_units == first.statistics_v2.work_units + second.statistics_v2.work_units &&
                result.engine_counters().evaluations == first.statistics_v3.evaluations + second.statistics_v3.evaluations,
            "successful axis handoff discarded the sparse attempt's work");

    // The derivative identity proves uniqueness. Independent exact dyadic
    // endpoint signs prove that the published interval contains that root;
    // no production sign evaluator or degree-8193 Sturm chain is used.
    const auto power = [](Rational base, uint64_t exponent) {
        Rational value = 1;
        while (exponent) {
            if (exponent & 1) value *= base;
            exponent >>= 1;
            if (exponent) base *= base;
        }
        return value;
    };
    const auto sign_value = [&](const Rational &x) -> Rational {
        return (m + 1) * power(x, n) - 6 * n * power(x, m + 1) + 9 * n * (m + 1) * x + n * (m + 1);
    };
    const auto roots = result.roots();
    require(roots.size() == 1 && roots[0].count == 1 && roots[0].multiplicity == 1 &&
                roots[0].realness == certpolysolve::Realness::real, "monotone fixture lost its sole simple real root");
    const auto lo = endpoint(roots[0].box.real.lo), hi = endpoint(roots[0].box.real.hi);
    require(-1 < lo && hi < 0 && sign_value(lo) < 0 && sign_value(hi) > 0,
            "monotone fixture failed independent endpoint signs");
    check_accuracy(roots[0], 53);
}

int main(int argc, char **argv) {
    certpolysolve::ThreadScope lifetime;
    try {
        if (argc == 2 && std::string(argv[1]) == "sparse-axis-accounting") {
            sparse_axis_accounting();
            return 0;
        }
        if (argc == 2 && std::string(argv[1]) == "axis-interruption") {
            axis_interruptions(graded_pairs(126));
            return 0;
        }
        require(argc == 1, "unknown real-refinement test variant");
        sparse_attempt();
        sparse_stop();
        sparse_axis_accounting();
        axis_attempt();
        axis_not_attempted();
        axis_factored();
        absolute_separation();
        real_solves({-2, 0, 1}, 2);
        real_solves(multiply(multiply({0, 1}, {-1, 1}), {-2, 0, 1}), 4);
        real_solves(multiply({-Rational(1) / 3, 1}, {-1 - power2(-150), 1}), 2);
        const auto close_dyadics = multiply({-1, 1}, {-1 - power2(-200), 1});
        real_solves(close_dyadics, 2);
        certpolysolve::Options close_options;
        close_options.bits = 128;
        close_options.max_bits = 4096;
        verify(close_dyadics, certpolysolve::real_roots(input(close_dyadics), close_options), 2, 128);
        // Fixed FLINT 3.6.0 arb_fmpz_poly/test/t-real_roots.c regressions.
        // Coefficients are data; counts below come from the independent Sturm code.
        const char *regressions[] = {
            "-1 -1 0 -1 -1 1 -1 0 0 0 0 0 1 0 0 0 1 -1 0 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 1 -1 1 -1 -1 0 0 0 0 "
            "1 0 1 0 0 0 1 0 1 0 1 0 -1 -1 0 -1 0 0 -1 -1 0 1 0 0 0 0 0 0 -1 1",
            "3 0 -1 0 -2", "0 -5 -1"};
        for (const auto text : regressions) {
            std::istringstream stream(text);
            Polynomial p;
            int coefficient;
            while (stream >> coefficient)
                p.emplace_back(coefficient);
            real_solves(p, roots_closed(p, -8, 8));
        }
        // The two frozen application families that stalled at a 53-bit request.
        const Rational a = 1 + power2(-100);
        real_solves({a * a - 1, 0, 1 - 2 * a, 0, 1}, 4);
        Polynomial previous{0}, current{1};
        const Rational b = Rational(1) / 512;
        for (int j = 0; j < 32; ++j) {
            auto next = multiply(current, {-Rational(j + 1) / 32, 1});
            for (size_t k = 0; k < previous.size(); ++k)
                next[k] -= b * b * previous[k];
            previous = current;
            current = std::move(next);
        }
        real_solves(current, 32);
        // Chebyshev T_64: 64 real roots, coefficients up to 2^63. Near a root
        // its value cancels about 64 bits, and the interval Newton step's
        // slope over the whole interval holds zero until the interval is far
        // narrower than the goal: bisection took 2,731 refinement steps (one
        // bit each, about 43 a root). Refinement by cells (refine_by_cells)
        // gains bits quadratically: five confirmed cells (4, 16, 256, 2^16,
        // 2^32 of them) narrow an interval by 62 bits. It took 527 steps.
        {
            Polynomial before{1}, chebyshev{0, 1};
            for (int k = 1; k < 64; ++k) {
                auto next = multiply({0, 2}, chebyshev);
                for (size_t j = 0; j < before.size(); ++j)
                    next[j] -= before[j];
                before = std::move(chebyshev);
                chebyshev = std::move(next);
            }
            certpolysolve::Options options;
            const auto result = certpolysolve::real_roots(input(chebyshev), options);
            verify(chebyshev, result, 64, 53);
            require(result.engine_counters().real_refinement_steps <= 12 * 64,
                    "real refinement of Chebyshev T_64 took more than 12 steps a root");
        }
        // A Mignotte cubic, x^3 - 2 (2^19 x - 1)^2: two roots about 2^-48
        // apart near 2^-19, and one near 2^39. At a 64-bit cap the large
        // root's interval reaches the cap first; the close pair, which still
        // shared an isolation endpoint, must go on being refined, so that
        // the precision stop returns a complete disjoint cover (it returned
        // none when one component at the cap stopped the others). Auto's
        // real-only query takes the complex route and, since 15 October 2026,
        // stops there at the cap (test_real_handoff); no public query on coefficient
        // input reaches the real engine any more (Auto's real-only queries on
        // native secular input do), so this case drives the D5 recipe
        // (production_policy()), whose real-only queries go to the engine,
        // directly (until W8's S2 through a session that skipped the selector).
        {
            const Polynomial mignotte{-2, 2097152, Rational("-549755813888"), 1};
            const auto source = input(mignotte);
            certpolysolve::detail::SolveState engine(certpolysolve::detail::dense(*source.get()),
                                                 certpolysolve::detail::maintained_method_policy(CPS_METHOD_AUTO));
            cps_options options;
            options.real_only = true;
            options.bits = 100;
            options.max_work_bits = 64;
            const auto result = certpolysolve_test::published(engine.advance(&options));
            require(result.message().find("real roots by the complex route") == std::string::npos,
                    "the D5 recipe did not run the real engine");
            require(cps_result_termination(result.get()) == CPS_PRECISION_LIMIT && result.complete() && result.roots().size() == 3,
                    "a component at the cap stopped the refinement of the others");
            check_disjoint(result.roots());
            for (const auto &root : result.roots())
                require(roots_closed(mignotte, endpoint(root.box.real.lo), endpoint(root.box.real.hi)) == 1,
                        "a region of the Mignotte cubic does not hold exactly one root");
        }
        certpolysolve::Options low;
        low.bits = 256;
        low.max_bits = 32;
        const auto r = certpolysolve::real_roots(input({-2, 0, 1}), low);
        require(!r.goal_met() && cps_result_termination(r.get()) == CPS_PRECISION_LIMIT, "precision cap ignored");
        low.bits = 256;
        low.max_bits = 4096;
        verify({-2, 0, 1}, certpolysolve::real_roots(input({-2, 0, 1}), low), 2, 256);
        certpolysolve::Options o;
        std::stop_source stopped;
        stopped.request_stop();
        o.stop = stopped.get_token();
        require(cps_result_termination(certpolysolve::real_roots(input({-2, 0, 1}), o).get()) == CPS_CANCELLED,
                "real cancellation ignored");
        o.stop = {};
        verify({-2, 0, 1}, certpolysolve::real_roots(input({-2, 0, 1}), o), 2, 53);
        std::cout << "real Sturm checks passed\n";
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
