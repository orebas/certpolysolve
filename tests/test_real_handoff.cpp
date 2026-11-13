// A real-only query on exact real input is answered by an attempt of the axis
// charts (piecewise's axis mode, two rounds) or, where it stops, by the complex
// route (session.cpp): its individual roots proved real by reflection. Checked
// against known roots: every real root in exactly one region, with its
// multiplicity; no region without a real root. Where the complex route stops
// at its precision cap the query stops there (since 15 October 2026; until
// then the old real engine was a second attempt), and the stop names a cap at
// which the query finishes. Input that is not exact and real stays
// unsupported; the solve's work and iteration budgets hold.
#include "exact_winding_api.hpp"
#include <certpolysolve/certpolysolve.hpp>
#include <iostream>

using namespace reference;

namespace {
const char *const handed_off = "real roots by the complex route";
const char *const by_axis = "real roots by the axis charts";

// Each expected real root (with its multiplicity) lies in exactly one region,
// that region is a real interval, and no region is left over.
void check_real(const certpolysolve::Result &result, const std::vector<std::pair<Rational, uint64_t>> &expected,
                const char *what) {
    require(result.goal_met() && cps_result_termination(result.get()) == CPS_SUCCESS && result.complete(), what);
    const auto &regions = result.roots();
    require(regions.size() == expected.size(), "the real cover does not have one region for each real root");
    for (const auto &region : regions) {
        require(region.realness == certpolysolve::Realness::real && endpoint(region.box.imag.lo) == 0 &&
                    endpoint(region.box.imag.hi) == 0,
                "a real-only region is not an interval of the axis");
        uint64_t inside = 0, multiplicity = 0;
        for (const auto &[root, m] : expected)
            if (endpoint(region.box.real.lo) <= root && root <= endpoint(region.box.real.hi)) {
                ++inside;
                multiplicity = m;
            }
        require(inside == 1 && region.count == multiplicity, "a real-only region does not hold exactly its real root");
    }
}
bool took_handoff(const certpolysolve::Result &result) {
    return result.message().find(handed_off) != std::string::npos;
}
// Answered by one of the two real-only routes (the axis charts' attempt, or the complex route after it).
bool real_route(const certpolysolve::Result &result) {
    return took_handoff(result) || result.message().find(by_axis) != std::string::npos;
}
// The cap a precision stop names ("so a cap of N bits suffices"), or 0.
uint64_t named_cap(const certpolysolve::Result &result) {
    const std::string said = result.message(), marker = "so a cap of ";
    const auto at = said.find(marker);
    return at == std::string::npos ? 0 : std::stoull(said.substr(at + marker.size()));
}
// The polynomial times 2^-100000, written in compact dyadic form: expanding
// it to integers takes far more storage than the input owns, which the
// default exponent policy refuses (closure-fresh-flat-1543-0 and the large
// historical inputs are refused so).
certpolysolve::Polynomial compact(const ComplexPolynomial &p) {
    certpolysolve::Polynomial result;
    for (size_t j = 0; j < p.size(); ++j) {
        require(p[j].imag == 0 && denominator(p[j].real) == 1, "the fixture needs integer coefficients");
        result.set(j, certpolysolve::Number(p[j].real.str() + "@-100000"));
    }
    return result;
}
} // namespace

int main() {
    certpolysolve::ThreadScope lifetime;
    try {
        // Three real roots and a complex pair.
        const auto p = from_roots({Complex(Rational(1)), Complex(Rational(2)), Complex(Rational(3)),
                                   Complex(Rational(0), Rational(1)), Complex(Rational(0), Rational(-1))});
        const std::vector<std::pair<Rational, uint64_t>> real{{Rational(1), 1}, {Rational(2), 1}, {Rational(3), 1}};
        for (const auto goal : {certpolysolve::Goal::accuracy(), certpolysolve::Goal::approximation()}) {
            certpolysolve::Options options;
            options.goal = goal;
            options.bits = 53;
            options.max_bits = 4096;
            // Compact input (which the real engine refuses), integer input
            // and compact input with the expansion allowed: all three are
            // answered by a real-only route.
            const auto refused = certpolysolve::real_roots(compact(p), options);
            check_real(refused, real, "no real-only route answered a real-only query on compact input");
            require(real_route(refused),
                    "the real-only query on compact input did not name its route");
            const auto integer = certpolysolve::real_roots(reference::input(p), options);
            check_real(integer, real, "no real-only route answered integer input");
            require(real_route(integer),
                    "the real-only query on integer input did not name its route");
            certpolysolve::ResourceLimits limits;
            limits.flags |= CPS_RESOURCE_ALLOW_EXPONENT_EXPANSION;
            certpolysolve::Options expanding;
            expanding.goal = goal;
            expanding.bits = 53;
            expanding.max_bits = 4096;
            expanding.resources = limits;
            const auto allowed = certpolysolve::real_roots(compact(p), expanding);
            check_real(allowed, real, "no real-only route answered input whose expansion was allowed");
            require(real_route(allowed),
                    "the real-only query with the expansion allowed did not name its route");
        }
        // A Mignotte cubic, x^3 - 2 (2^19 x - 1)^2, has two real roots about 2^-48 apart: at a 64-bit cap
        // the complex route cannot separate them. The query stops at the cap (the real engine's second attempt
        // isolated them exactly until 15 October 2026) and names a cap (Mahler's separation bound of the
        // square-free part, solver.cpp); at that cap the same query meets its goal, for every root at the
        // accuracy goal and for the real roots alike.
        {
            certpolysolve::Polynomial mignotte;
            mignotte.set(0, certpolysolve::Number("-2"))
                .set(1, certpolysolve::Number("2097152"))
                .set(2, certpolysolve::Number("-549755813888"))
                .set(3, certpolysolve::Number("1"));
            for (const bool real_only : {false, true}) {
                certpolysolve::Options options;
                options.goal = certpolysolve::Goal::accuracy();
                options.bits = 100;
                options.max_bits = 64;
                const auto query = [&](const certpolysolve::Options &o) {
                    return real_only ? certpolysolve::real_roots(mignotte, o) : certpolysolve::solve(mignotte, o);
                };
                const auto stopped = query(options);
                const uint64_t cap = named_cap(stopped);
                require(cps_result_termination(stopped.get()) == CPS_PRECISION_LIMIT && !stopped.goal_met() &&
                            cps_result_stop_reason(stopped.get()) == CPS_STOP_PRECISION &&
                            stopped.message().find("stopped at the precision cap of 64 bits") != std::string::npos &&
                            cap > 64,
                        "a precision stop did not name the cap that suffices");
                options.max_bits = cap;
                const auto finished = query(options);
                require(finished.goal_met() && cps_result_termination(finished.get()) == CPS_SUCCESS &&
                            finished.roots().size() == 3 &&
                            finished.work_bits() <= cap && named_cap(finished) == 0,
                        "the cap a precision stop named did not suffice");
                if (real_only) // all three roots are real: x^3 - 2^39 x^2 + 2^21 x - 2 (the real engine found three)
                    for (const auto &region : finished.roots())
                        require(region.realness == certpolysolve::Realness::real && region.count == 1,
                                "a real-only answer at the named cap published a region that is not one real root");
            }
        }
        // A triple real root in compact form: without exact preparation the
        // complex route cannot establish its multiplicity, and stops at the
        // precision cap: the simple real root, incomplete.
        {
            const auto triple = from_roots({Complex(Rational(1)), Complex(Rational(1)), Complex(Rational(1)),
                                            Complex(Rational(-2)), Complex(Rational(0), Rational(1)),
                                            Complex(Rational(0), Rational(-1))});
            certpolysolve::Options options;
            options.bits = 53;
            options.max_bits = 256;
            const auto result = certpolysolve::real_roots(compact(triple), options);
            require(!result.complete() && !result.goal_met() && cps_result_termination(result.get()) == CPS_PRECISION_LIMIT &&
                        cps_result_stop_reason(result.get()) == CPS_STOP_PRECISION && took_handoff(result),
                    "an unresolved real cluster did not stop honestly at the precision cap");
            for (const auto &region : result.roots())
                require(region.realness == certpolysolve::Realness::real && region.count == 1 && endpoint(region.box.real.lo) <= -2 &&
                            -2 <= endpoint(region.box.real.hi),
                        "an incomplete real-only answer published a region that is not a proved real root");
        }
        // The same triple root under small work allowances: a solve is never an unsupported query, and every
        // region it publishes is the proved simple root (written for the real engine's hand-over, review of 14
        // October 2026, S02, as a session advanced call by call; one-shot solves since 18 October 2026).
        {
            const auto triple = from_roots({Complex(Rational(1)), Complex(Rational(1)), Complex(Rational(1)),
                                            Complex(Rational(-2)), Complex(Rational(0), Rational(1)),
                                            Complex(Rational(0), Rational(-1))});
            for (uint64_t allowance : {1, 2, 3, 5, 8, 13}) {
                certpolysolve::Options options;
                options.bits = 53;
                options.max_bits = 256;
                options.work_budget = allowance;
                options.cell_budget = CPS_DEFAULT_LIVE_CELLS;
                const auto result = certpolysolve::real_roots(compact(triple), options);
                require(cps_result_termination(result.get()) != CPS_UNSUPPORTED_QUERY,
                        ("work " + std::to_string(allowance) + ": the solve lost the complex route's answer").c_str());
                for (const auto &region : result.roots())
                    require(region.realness == certpolysolve::Realness::real && region.count == 1 && endpoint(region.box.real.lo) <= -2 &&
                                -2 <= endpoint(region.box.real.hi),
                            ("work " + std::to_string(allowance) + ": a published region is not the proved simple root")
                                .c_str());
            }
        }
        // The same close pair under iteration budgets up to the sweeps the solve needs: each solve stops for its
        // budget or at the cap, and the one with the sweeps it needs ends at the cap with the stop that names one.
        {
            certpolysolve::Polynomial mignotte;
            mignotte.set(0, certpolysolve::Number("-2"))
                .set(1, certpolysolve::Number("2097152"))
                .set(2, certpolysolve::Number("-549755813888"))
                .set(3, certpolysolve::Number("1"));
            certpolysolve::Options unlimited;
            unlimited.bits = 100;
            unlimited.max_bits = 64;
            const uint64_t needed = certpolysolve::real_roots(mignotte, unlimited).statistics().sweeps;
            for (uint64_t budget = 1; budget <= needed; ++budget) {
                certpolysolve::Options options;
                options.bits = 100;
                options.max_bits = 64;
                options.iteration_budget = budget;
                const auto result = certpolysolve::real_roots(mignotte, options);
                const auto stop = cps_result_stop_reason(result.get());
                require(stop == CPS_STOP_ITERATION_BUDGET || (stop == CPS_STOP_PRECISION && named_cap(result) > 64),
                        "a solve under an iteration budget ended with another stop, or without a named cap");
                if (budget == needed)
                    require(stop == CPS_STOP_PRECISION, "a solve with the sweeps it needs did not end at the cap");
            }
        }
        // Input that is not real stays unsupported.
        {
            const auto complex = from_roots({Complex(Rational(1)), Complex(Rational(0), Rational(1))});
            certpolysolve::Options options;
            const auto result = certpolysolve::real_roots(reference::input(complex), options);
            require(cps_result_termination(result.get()) == CPS_UNSUPPORTED_QUERY && result.roots().size() == 0 && !took_handoff(result),
                    "a real-only query on complex input was answered");
        }
        // One work unit: the solve spends no more, and stops honestly or answers.
        {
            certpolysolve::Options options;
            options.bits = 53;
            options.max_bits = 4096;
            options.work_budget = 1;
            options.cell_budget = CPS_DEFAULT_LIVE_CELLS;
            const auto result = certpolysolve::real_roots(compact(p), options);
            require(result.statistics().work_units <= 1, "a solve spent more than its work budget");
            if (result.goal_met())
                check_real(result, real, "the handoff under a one-unit work budget answered wrongly");
            else
                require(cps_result_termination(result.get()) == CPS_WORK_LIMIT, "an unfinished solve reported another stop");
        }
        // The complex route is under the solve's iteration budget: under each budget up to the sweeps it needs,
        // no solve makes more, and the last one answers.
        {
            certpolysolve::Options unlimited;
            unlimited.bits = 53;
            unlimited.max_bits = 4096;
            const uint64_t needed = certpolysolve::real_roots(compact(p), unlimited).statistics().sweeps;
            for (uint64_t budget = 1; budget <= needed; ++budget) {
                certpolysolve::Options options;
                options.bits = 53;
                options.max_bits = 4096;
                options.iteration_budget = budget;
                const auto result = certpolysolve::real_roots(compact(p), options);
                require(result.statistics().sweeps <= budget,
                        "a solve made more sweeps than its iteration budget");
                if (result.goal_met())
                    check_real(result, real, "the handoff under an iteration budget answered wrongly");
                if (budget == needed) require(result.goal_met(), "the handoff did not answer with the sweeps it needs");
            }
        }
        std::cout << "real-only handoff tests passed\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "FAILED: " << error.what() << '\n';
        return 1;
    }
}
