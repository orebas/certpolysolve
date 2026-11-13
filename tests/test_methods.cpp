#include "api_test_support.hpp"
#include "exact_winding_api.hpp"
#include <iostream>
#include <stdexcept>

using namespace reference;
// The solve's termination as the C result reports it. An error result (no C result behind it) fails here with its
// message: the 0.5 solve returns an error where the old one threw.
static cps_termination termination(const certpolysolve::Result &result) {
    if (!result.get()) throw std::runtime_error("the solve returned an error: " + result.message());
    return cps_result_termination(result.get());
}
static void verify(const certpolysolve::Result &result, const std::vector<Complex> &known, int bits) {
    auto regions = result.roots();
    check_disjoint(regions);
    uint64_t total = 0;
    for (const auto &root : regions) {
        uint64_t count = 0;
        for (const auto &z : known)
            count += endpoint(root.box.real.lo) <= z.real && z.real <= endpoint(root.box.real.hi) &&
                     endpoint(root.box.imag.lo) <= z.imag && z.imag <= endpoint(root.box.imag.hi);
        require(count == root.count, "partial/final method cover has an incorrect original-input count");
        if (root.multiplicity) require(*root.multiplicity == count, "incorrect restored multiplicity");
        if (result.goal_met()) check_accuracy(root, bits);
        total += count;
    }
    if (result.complete()) require(total == known.size(), "complete method cover lost roots");
}
// Each explicit method on simple and repeated roots: a one-unit solve returns a correct partial cover, a full
// solve meets the 53-bit goal, and one at 160 bits meets that. (Until 18 October 2026 this was one session
// advanced a unit at a time, then tightened; the 0.5 API has one-shot solves only.)
static void method_covers(certpolysolve::Method method, const certpolysolve::Goal &goal, bool repeated) {
    std::vector<Complex> known{{-2, 1}, {Rational(1)/3, Rational(1)/5}, {2, -1}, {3, 2}};
    if (repeated) {
        auto copy = known;
        known.insert(known.end(), copy.begin(), copy.end());
        known.insert(known.end(), 3, Complex(0));
        known.insert(known.end(), 3, Complex(5));
    }
    const auto p = input(from_roots(known));
    certpolysolve::Options options;
    options.method = method;
    options.goal = goal;
    options.bits = 53;
    options.max_bits = 1024;
    options.work_budget = 1;
    // The 0.5 Options cannot hold an invalid method, so the refusal is checked on C options.
    const auto raw = certpolysolve_test::raw_options();
    require(cps_options_set_method(raw.get(), static_cast<cps_method>(99)) == CPS_INVALID_ARGUMENT,
            "invalid method accepted");
    const auto first = certpolysolve::solve(p, options);
    require(first.statistics().work_units <= 1, "a solve exceeded its work allowance");
    verify(first, known, 53);
    options.work_budget = 100000;
    const auto result = certpolysolve::solve(p, options);
    if (!result.goal_met()) throw std::runtime_error("method failed: " + result.message());
    verify(result, known, 53);
    options.bits = 160;
    options.max_bits = 2048;
    const auto refined = certpolysolve::solve(p, options);
    require(refined.goal_met(), "method failed at 160 bits");
    verify(refined, known, 160);
}
static void cancelled_and_unsupported() {
    const auto p = input(from_roots(std::vector<Complex>{{1, 1}, {2, -1}, {3, 2}}));
    certpolysolve::Options o;
    o.method = certpolysolve::Method::piecewise;
    // A real-only query with an explicit method is a usage error (C7): the solve's error, no result.
    const auto unsupported = certpolysolve::real_roots(p, o);
    require(unsupported.stop() == certpolysolve::Stop::error && !unsupported.get() && !unsupported.message().empty(),
            "a real-only query with an explicit method was not a usage error");
    // A real-only query on nonreal input can never be answered: it is
    // unsupported at every work budget, with the input's stop code on every
    // route (fuzz F8: a one-unit budget reported a work limit, as if more
    // work could help, and Auto left the stop code unset).
    certpolysolve::Polynomial nonreal;
    nonreal.set(0, certpolysolve::Number::complex(certpolysolve::Number("0"), certpolysolve::Number("-1")))
        .set(1, certpolysolve::Number("1")); // x - i
    for (const uint64_t work : {uint64_t(1), uint64_t(2), uint64_t(100000)}) {
        certpolysolve::Options real_query;
        real_query.work_budget = work;
        const auto refused = certpolysolve::real_roots(nonreal, real_query);
        require(termination(refused) == CPS_UNSUPPORTED_QUERY &&
                    cps_result_stop_reason(refused.get()) == CPS_STOP_INPUT && !refused.complete() &&
                    refused.roots().empty(),
                "a real-only query on nonreal input was not refused as unsupported input");
    }
    // A cancellation at the first poll (certpolysolve_testing.h). Until stage 3b-3 removed the callback, a one-shot
    // callback answered here, so that the solve had to latch it; test_solve_state keeps that latch's test, on the
    // internal callback.
    const auto raw = certpolysolve_test::raw_options();
    require(cps_options_set_method(raw.get(), CPS_METHOD_PIECEWISE) == CPS_OK, "piecewise method refused");
    require(cps_testing_options_stop_at_poll(raw.get(), 1) == CPS_OK, "the poll counter refused");
    const auto stopped = certpolysolve_test::solve(p, raw.get());
    require(termination(stopped) == CPS_CANCELLED && !stopped.complete(), "a cancellation at the first poll was lost");
    o.method = certpolysolve::Method::classical;
    certpolysolve::Polynomial uncertain;
    uncertain.set(0, certpolysolve::ball(certpolysolve::Number("-1"), certpolysolve::Number("1@-100"))).set(1, certpolysolve::Number("1"));
    require(termination(certpolysolve::solve(uncertain, o)) == CPS_UNSUPPORTED_QUERY,
            "explicit method silently accepted uncertain coefficients");
    certpolysolve::Polynomial zero;
    o.method = certpolysolve::Method::piecewise;
    require(termination(certpolysolve::solve(zero, o)) == CPS_ZERO_POLYNOMIAL, "zero polynomial lost its status");
    zero.set(0, certpolysolve::Number("7"));
    require(certpolysolve::solve(zero, o).goal_met(), "nonzero constant failed");
    const auto number = [](const char *re, const char *im) {
        return certpolysolve::Number::complex(certpolysolve::Number(re), certpolysolve::Number(im));
    };
    certpolysolve::SecularEquation native{{.residue = number("1", "0"), .pole = number("0", "0")},
                                      {.residue = number("1", "0"), .pole = number("2", "0")}};
    require(native.ok(), "the native secular input was rejected");
    require(certpolysolve::solve(native, o).stop() == certpolysolve::Stop::error,
            "piecewise on native secular input was not a usage error (it would expand the input)");
    o.method = certpolysolve::Method::secular;
    require(certpolysolve::solve(native, o).goal_met(), "explicit native secular query failed");
}
// (The zero-work check of review S01 moved to test_method_handoff in stage 3b-3: zero work is an internal setting
// since then, as no public call sets it.)
// Explicit piecewise hands every stop that is not the caller's to secular:
// its precision budget and a storage refusal as well as the regime exit
// (principled audit, fix 3). Below piecewise's minimum precision it used to
// return no regions at all; Wilkinson-20 at 160 bits stopped with 18 of 20.
static void piecewise_hands_off_budget_stops() {
    std::vector<Complex> known;
    for (int k = 1; k <= 20; ++k) known.emplace_back(k);
    const auto p = input(from_roots(known));
    // At 128 bits the cap is below piecewise's first working precision (149
    // for 53 bits): it refuses before any work and secular finishes. At 160
    // bits the ladder's last rung is cut to the cap (queue item 3.2) and
    // piecewise finishes there itself.
    for (const uint64_t cap : {128u, 160u}) {
        certpolysolve::Options o;
        o.method = certpolysolve::Method::piecewise;
        o.goal = certpolysolve::Goal::approximation();
        o.bits = 53;
        o.max_bits = cap;
        const auto result = certpolysolve::solve(p, o);
        const bool handed = result.message().find("piecewise precision budget to secular") != std::string::npos;
        require(result.goal_met() && handed == (cap == 128),
                cap == 128 ? "piecewise did not hand its precision budget to secular"
                           : "piecewise did not finish at a cap its ladder can reach");
        verify(result, known, 53);
        // The call's work budget is charged once across the handoff.
        o.work_budget = 40;
        const auto limited = certpolysolve::solve(p, o);
        require(limited.statistics().work_units <= 40,
                "the handoff gave secular a second work budget in the same call");
    }
    // Coefficient exponents beyond piecewise's envelope range (|height| >
    // 2^50) are a refusal of the input, not an internal error that poisons
    // the solve (principled audit, fix 14): x^3 + 2^(2^50 + 6) x - 1.
    certpolysolve::Polynomial huge;
    huge.set(0, certpolysolve::Number("-1")).set(1, certpolysolve::Number("1@1125899906842630")).set(3, certpolysolve::Number("1"));
    certpolysolve::Options o;
    o.method = certpolysolve::Method::piecewise;
    o.goal = certpolysolve::Goal::approximation();
    o.work_budget = 200;
    const auto refused = certpolysolve::solve(huge, o);
    require(refused.message().find("piecewise range refusal to secular") != std::string::npos &&
                refused.complete(),
            "piecewise's exponent range limit was not a clean refusal");
    uint64_t total = 0;
    for (const auto &region : refused.roots()) total += region.count;
    require(total == 3, "the refused input's cover lost roots");
    // The starting points of such an input stand at its roots' own scale. Their
    // radius used to be clamped to 2^(+-1e12), a factor 2^(5.6e14) inside this
    // input's two large roots: secular crept outward and never converged (no
    // budget was enough), and classical needed 553 units through Subdivision.
    for (const auto method : {certpolysolve::Method::automatic, certpolysolve::Method::classical, certpolysolve::Method::secular,
                              certpolysolve::Method::piecewise}) {
        o.method = method;
        const auto solved = certpolysolve::solve(huge, o);
        require(solved.goal_met() && solved.roots().size() == 3,
                "an input with exponents near 2^50 did not meet its goal within 200 units");
    }
    // Past 2^53 the exponents are no longer exact as doubles and no starting
    // scale is derived from them: the solve goes through Subdivision instead
    // (375 units measured), and still meets its goal.
    certpolysolve::Polynomial larger;
    larger.set(0, certpolysolve::Number("-1")).set(1, certpolysolve::Number("1@1152921504606846982")).set(3, certpolysolve::Number("1")); // x^3 + 2^(2^60 + 6) x - 1
    for (const auto method : {certpolysolve::Method::classical, certpolysolve::Method::secular}) {
        o.method = method;
        o.work_budget = 4000;
        const auto solved = certpolysolve::solve(larger, o);
        require(solved.goal_met() && solved.roots().size() == 3,
                "an input with exponents near 2^60 did not meet its goal");
    }
}
// The square-free factorisation cannot be cancelled once started, so it is
// charged to exact_operations like the import before it (principled audit,
// fix 14). (x - 1)^2 (x^9 + B) with B = 3^6300, about 10,000 bits: the import is
// 12^2 operations, the factorisation about 12^2 x 157. Its square-free twin
// (x - 1)(x^9 + B) needs no factorisation since S8's B2: the modular proof, one
// operation an imported coefficient, shows it square-free within the same
// budget, and its one factor is the import's primitive part.
static void squarefree_is_charged() {
    // Plain integer text: a mantissa@exponent coefficient would take the
    // compact-exponent route, which is never expanded without an allowance.
    Rational power(1);
    for (int i = 0; i < 6300; ++i) power *= 3; // about 2^9985
    const std::string big = power.str(), twice = (power * Rational(2)).str();
    certpolysolve::Polynomial repeated, square_free;
    // (x^2 - 2x + 1)(x^9 + B) = x^11 - 2 x^10 + x^9 + B x^2 - 2B x + B
    repeated.set(0, certpolysolve::Number(big)).set(1, certpolysolve::Number("-" + twice)).set(2, certpolysolve::Number(big))
        .set(9, certpolysolve::Number("1")).set(10, certpolysolve::Number("-2")).set(11, certpolysolve::Number("1"));
    square_free.set(0, certpolysolve::Number("-" + big)).set(1, certpolysolve::Number(big)).set(9, certpolysolve::Number("-1"))
        .set(10, certpolysolve::Number("1"));
    certpolysolve::Options o;
    o.method = certpolysolve::Method::classical;
    o.goal = certpolysolve::Goal::approximation();
    const auto prepared = certpolysolve::solve(repeated, o);
    require(prepared.goal_met() &&
                prepared.message().find("budgeted exact real squarefree preprocessing") != std::string::npos,
            "an affordable square-free preparation did not run");
    certpolysolve::ResourceLimits limits;
    limits.exact_operations = 1000;
    o.resources = limits;
    const auto limited = certpolysolve::solve(repeated, o);
    require(limited.message().find("budgeted exact real squarefree preprocessing") == std::string::npos,
            "the square-free factorisation ran beyond the exact operation budget");
    uint64_t total = 0;
    for (const auto &region : limited.roots()) total += region.count;
    require(limited.complete() && total == 11, "the unprepared solve lost its cover");
    const auto proved = certpolysolve::solve(square_free, o);
    require(proved.goal_met() &&
                proved.message().find("budgeted exact real squarefree preprocessing") != std::string::npos,
            "the modular proof did not prepare a square-free input within the budget");
}
// Every piecewise stop short of the goal returns a complete cover (principled
// audit fix 3, queue item 3.2): the regions it has certified, and the roots it
// has not in one disk about the origin whose count is proved on the input; or,
// when no such disk separates them, Cauchy's global disk. The counts sum to
// the degree, every region's count is right (the oracle checks it against the
// known roots), the goal is not met, and the stop keeps its reason. The input
// is fuzz finding F6's: 22 real roots and a four-root cluster 2^-38 wide.
static void piecewise_stops_with_a_complete_cover() {
    std::vector<Complex> known;
    const auto at = [&](long numerator, long denominator) { known.emplace_back(Rational(numerator) / denominator); };
    at(199, 256); at(-37, 772); at(-105, 233); at(751, 1024); at(247, 256); at(-20, 17); at(214, 109);
    at(-53875, 65536); at(25, 32); at(1189, 2048); at(27489, 131072); at(27839, 32768); at(-373, 512);
    at(51, 64); at(-1, 1); at(-42467, 65536); at(28, 45); at(-42175, 65536); at(0, 1); at(-27, 128);
    at(83, 512); at(-1415, 2048); at(3, 4);
    for (const int offset : {3, 1, -1}) known.emplace_back(Rational(3) / 4 + offset * power2(-38));
    const auto p = input(from_roots(known));
    bool coarse_seen = false, informative_seen = false;
    for (const uint64_t budget : {1u, 20u, 120u, 300u}) {
        certpolysolve::Options o;
        o.method = certpolysolve::Method::piecewise;
        o.goal = certpolysolve::Goal::approximation();
        o.bits = 100;
        o.max_bits = 256;
        o.work_budget = budget;
        const auto stopped = certpolysolve::solve(p, o);
        if (stopped.goal_met()) continue;
        require(stopped.complete() && termination(stopped) == CPS_WORK_LIMIT,
                "a piecewise work-budget stop did not return a complete cover with its reason");
        verify(stopped, known, 100);
        const auto regions = stopped.roots();
        if (regions.size() == 1) coarse_seen = true;
        if (regions.size() > 1) informative_seen = true;
        for (const auto &region : regions)
            require(region.count < known.size() || regions.size() == 1, "a coarse region beside others holds every root");
    }
    require(coarse_seen || informative_seen, "no work budget stopped the solve short of its goal");
    certpolysolve::Options o;
    o.method = certpolysolve::Method::piecewise;
    o.goal = certpolysolve::Goal::approximation();
    o.bits = 100;
    o.max_bits = 256;
    o.work_budget = 100000;
    const auto finished = certpolysolve::solve(p, o);
    require(finished.goal_met(), "the input was not solved within a large budget");
    verify(finished, known, 100);
}
int main() {
    certpolysolve::ThreadScope lifetime;
    try {
        for (const auto method : {certpolysolve::Method::classical, certpolysolve::Method::secular, certpolysolve::Method::piecewise})
            for (const auto &goal : {certpolysolve::Goal::approximation(), certpolysolve::Goal::accuracy()})
                for (bool repeated : {false, true}) method_covers(method, goal, repeated);
        cancelled_and_unsupported();
        piecewise_hands_off_budget_stops();
        piecewise_stops_with_a_complete_cover();
        squarefree_is_charged();
        std::cout << "explicit method ownership, counts, work budgets and cancellation pass\n";
        return 0;
    } catch (const std::exception &error) { std::cerr << error.what() << '\n'; return 1; }
}
