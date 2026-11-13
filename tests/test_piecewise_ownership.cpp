#include "solve/prepared_piecewise.hpp"
#include "numeric/arithmetic/polynomial_arithmetic.hpp"
#include "exact_winding_api.hpp"
#include "methods/piecewise/piecewise_ownership.hpp"
#include "methods/piecewise/piecewise_roots.hpp"
#include <iostream>

namespace d = certpolysolve::detail;
using namespace reference;
static Rational endpoint(const d::Endpoint &e) {
    return reference::endpoint(certpolysolve::Dyadic{e.mantissa_text(), e.exponent_text()});
}
static d::Root rectangle(const ComplexPolynomial &p, double left, double right,
                         double bottom, double top) {
    d::Root r;
    d::Float lo, hi;
    arf_set_d(lo.value(), left); arf_set_d(hi.value(), right);
    arb_set_interval_arf(acb_realref(r.box.value()), lo.value(), hi.value(), 256);
    arf_set_d(lo.value(), bottom); arf_set_d(hi.value(), top);
    arb_set_interval_arf(acb_imagref(r.box.value()), lo.value(), hi.value(), 256);
    r.publish(256);
    const auto count = rectangle_count(p, endpoint(r.endpoints[0]), endpoint(r.endpoints[1]),
                                      endpoint(r.endpoints[2]), endpoint(r.endpoints[3]));
    require(count && *count > 0, "test rectangle has no independent count proof");
    r.count = static_cast<uint64_t>(*count);
    r.certificate = "independent rational winding in ownership test";
    return r;
}
static void checked_cover(const cps_result &result, const std::vector<Complex> &known,
                          const char *stage = "solve") {
    if (!result.goal_met || !d::disjoint(result.roots)) throw std::runtime_error(std::string(stage) +
        ": ownership cover incomplete; degree=" + std::to_string(result.degree) +
        " regions=" + std::to_string(result.roots.size()) + " precision=" + std::to_string(result.work_bits) +
        " termination=" + std::to_string(result.termination));
    uint64_t count = 0;
    for (const auto &r : result.roots) {
        uint64_t matched = 0;
        for (const auto &z : known)
            if (endpoint(r.endpoints[0]) <= z.real && z.real <= endpoint(r.endpoints[1]) &&
                endpoint(r.endpoints[2]) <= z.imag && z.imag <= endpoint(r.endpoints[3])) ++matched;
        require(matched == r.count, "ownership output disagrees with exact planted roots");
        count += matched;
    }
    require(count == known.size(), "ownership cover lost multiplicity");
}
static d::PiecewiseRootConfig configuration() {
    d::PiecewiseRootConfig c;
    c.reuse_guards = c.reuse_polish = c.blocked_polish = true;
    c.conjugate_symmetry = c.recover_discovery = c.durable_ownership = true;
    c.adaptive_taylor = c.interior_proposals = true;
    return c;
}
static void close_root_progress() {
    for (const auto &separation : std::array<Rational, 3>{power2(-40), power2(-55), 7 * power2(-60)}) {
        const std::vector<Complex> known{{1}, {1 + separation, separation}, {-2}, {0, 3}, {-1, -2}, {2, 2}};
        auto input = reference::input(from_roots(known));
        cps_options query;
        query.goal = CPS_ROOT_APPROXIMATION;
        query.max_work_bits = 2048;
        d::PiecewiseRootStatistics stats;
        checked_cover(d::solve_piecewise_roots(certpolysolve::detail::dense(*input.get()), query, stats, configuration()), known,
                      "close pair among six roots");
        require(stats.seed_guards>0,"unpublished candidate points never entered singleton seed guards");
    }
    const auto separation = 7 * power2(-60);
    const std::vector<Complex> known{{1}, {1 + separation, separation}};
    auto input = reference::input(from_roots(known));
    for (const auto goal : {CPS_ROOT_APPROXIMATION, CPS_ROOT_ACCURACY}) {
        cps_options query;
        query.goal = goal;
        query.max_work_bits = 2048;
        d::PiecewiseRootStatistics stats;
        checked_cover(d::solve_piecewise_roots(certpolysolve::detail::dense(*input.get()), query, stats, configuration()), known,
                      "overlapping distinct-root certificates");
    }
}
static void interrupted_pending_work() {
    const std::vector<Complex> known{{1, 2}, {-3, 1}, {2, -1}, {5}, {0, 3}, {-2, -3}};
    auto input = reference::input(from_roots(known));
    // Cancel inside different chart/proposal/certification stages, then resume
    // the same state. A fresh solve cannot hide lost pending work.
    for (const unsigned stop : {40U, 100U, 200U, 400U, 800U}) {
        struct Cancellation { unsigned calls = 0, stop = 0; } cancellation{0, stop};
        cps_options query;
        query.goal = CPS_ROOT_APPROXIMATION;
        query.cancel_data = &cancellation;
        query.cancel = [](void *data) {
            auto &c = *static_cast<Cancellation *>(data);
            return static_cast<int>(++c.calls >= c.stop);
        };
        d::PiecewiseRootState state(certpolysolve::detail::dense(*input.get()), configuration());
        const auto partial = state.advance(query);
        require(partial.complete || partial.termination == CPS_CANCELLED,
                "pending-work interruption did not return cancellation");
        query.cancel = nullptr;
        const auto resumed = state.advance(query);
        checked_cover(resumed, known, "resumed pending certificate");
    }
}
int main() {
    certpolysolve::ThreadScope lifetime;
    try {
        // Publication rounds exact endpoints outward, then encloses them in
        // a ball. Near a binade boundary that final ball can be strictly larger
        // than the certified rectangle. Re-admitting the same proof must still
        // identify the same owner, using the published region on both sides.
        {
            const auto p=from_roots({Complex(1)});
            d::Root r;
            acb_one(r.box.value());
            arb_add_error_2exp_si(acb_realref(r.box.value()),-300);
            arb_add_error_2exp_si(acb_imagref(r.box.value()),-300);
            r.publish(128);
            r.count=1;
            d::Float actual_lower;
            arb_get_lbound_arf(actual_lower.value(),acb_realref(r.box.value()),256);
            d::RealBall published_lower;
            d::Number(r.endpoints[0].mantissa_text()+"@"+r.endpoints[0].exponent_text()).enclose(published_lower.value(),256);
            require(arf_cmp(actual_lower.value(),arb_midref(published_lower.value()))<0,
                    "publication regression no longer exercises strict ball enlargement");
            require(rectangle_count(p,endpoint(r.endpoints[0]),endpoint(r.endpoints[1]),
                                    endpoint(r.endpoints[2]),endpoint(r.endpoints[3]))==1,
                    "publication reproduction lacks an independent count");
            d::PiecewiseOwnership publication_ledger(1);
            const auto first=publication_ledger.admit(r);
            require(first && publication_ledger.identifies(*first,r) &&
                    publication_ledger.admit(r)==first && publication_ledger.count()==1,
                    "a published proof does not identify itself after outward ball rounding");
        }
        close_root_progress();
        interrupted_pending_work();
        const auto polynomial = from_roots({Complex(1), Complex(3), Complex(-2)});
        d::PiecewiseOwnership ledger(3);
        const auto first = rectangle(polynomial, .5, 2.5, -.25, .25);
        const auto a = ledger.admit(first);
        require(a && ledger.count() == 1, "first ownership missing");
        const auto overlap = rectangle(polynomial, 2, 3.5, -.25, .25);
        require(!ledger.admit(overlap) && ledger.count() == 1 && ledger.active(*a),
                "partial overlap silently identified different roots or mutated count");
        require(!ledger.identifies(*a, overlap), "intersection substituted for ownership identity");
        const auto refined = rectangle(polynomial, .75, 1.25, -.125, .125);
        {
            d::PiecewiseOwnership nesting(3);
            const auto narrow=nesting.admit(refined);
            d::ComplexBall midpoint;
            acb_get_mid(midpoint.value(),first.box.value());
            require(narrow && !nesting.containing(midpoint.value()) && nesting.identifying(first)==narrow,
                    "reverse proof nesting relied on the wider proof midpoint");
        }
        require(ledger.admit(refined) == a && ledger.count() == 1,
                "contained refinement duplicated an existing owner");
        const auto second = rectangle(polynomial, 2.75, 3.25, -.125, .125);
        const auto b = ledger.admit(second);
        require(b && *b != *a && ledger.count() == 2 && ledger.singletons(), "disjoint ownership count failed");
        const auto cluster = rectangle(polynomial, .25, 3.75, -.5, .5);
        const auto c = ledger.admit(cluster);
        require(c && *c != *a && *c != *b && !ledger.active(*a) && !ledger.active(*b) &&
                ledger.active(*c) && ledger.count() == 2 && !ledger.singletons(),
                "cluster replacement lost count conservation or stable IDs");
        d::RootModel chart;
        acb_set_si(chart.center.value(), 2); arb_set_si(chart.radius.value(), 2);
        require(ledger.contained_in(chart, 256) == 2, "contained cluster count not available to chart");
        arb_set_d(chart.radius.value(), .5);
        require(ledger.contained_in(chart, 256) == 0, "midpoint containment retired a partial region");
        // Exact endpoint storage must survive widths far below 64-bit rounding
        // at a nonzero center; no midpoint-only or approximate ownership check.
        d::Root narrow;
        acb_set_si(narrow.box.value(), -2);
        arb_add_error_2exp_si(acb_realref(narrow.box.value()), -180);
        arb_add_error_2exp_si(acb_imagref(narrow.box.value()), -180);
        narrow.publish(256);
        narrow.count = 1;
        require(rectangle_count(polynomial, endpoint(narrow.endpoints[0]), endpoint(narrow.endpoints[1]),
                               endpoint(narrow.endpoints[2]), endpoint(narrow.endpoints[3])) == 1,
                "narrow test proof failed winding");
        const auto n = ledger.admit(narrow);
        d::ComplexBall point;
        acb_set_si(point.value(), -2);
        require(n && ledger.contains(*n, point.value()) && ledger.containing(point.value()) == n &&
                ledger.count() == 3, "exact narrow ownership was lost");

        for (const auto &known : {std::vector<Complex>{{0}, {0}, {-2}, {1}, {0, 1}, {0, -1}, {2, 1}, {2, -1}},
                                 std::vector<Complex>{{1, 2}, {-3, 1}, {2, -1}, {5}}}) {
            auto input = reference::input(from_roots(known));
            d::PiecewiseRootState state(certpolysolve::detail::dense(*input.get()), configuration());
            cps_options query;
            query.goal = CPS_ROOT_APPROXIMATION;
            query.work_units = 1;
            cps_result partial;
            do {
                partial = state.advance(query);
                require(partial.complete || partial.progress.work_units_this_call > 0,
                        "ownership continuation exhausted precision");
            } while (!partial.complete && state.statistics().discovery_count <= 2);
            const auto before = state.statistics().discovery_count;
            query.cancel = [](void *) { return 1; };
            state.advance(query);
            require(state.statistics().discovery_count == before, "cancellation discarded proven ownership");
            query.cancel = nullptr;
            query.work_units = 100000;
            const auto first_answer = state.advance(query);
            checked_cover(first_answer, known, "initial session");
            require(state.statistics().discovery_count == known.size(), "discovery coverage did not account for degree");
            query.bits = 192;
            query.goal = CPS_ROOT_ACCURACY;
            const auto finer = state.advance(query);
            checked_cover(finer, known, "tighter session");
            require(state.statistics().discovery_count == known.size(), "goal/precision change discarded ownership");
            for (const auto &r : finer.roots)
                require(d::meets_goal(r, query, static_cast<slong>(finer.work_bits)), "changed goal kept inaccurate roots");
        }
        {
            const auto epsilon = power2(-100);
            const std::vector<Complex> known{{Rational(1)+epsilon}, {Rational(1)-epsilon}};
            auto input = reference::input(from_roots(known));
            d::PiecewiseRootState state(certpolysolve::detail::dense(*input.get()), configuration());
            cps_options query;
            query.goal = CPS_ROOT_APPROXIMATION;
            checked_cover(state.advance(query), known, "initial cluster");
            query.bits = 80; // A tighter counted cluster, still wider than the pair separation.
            checked_cover(state.advance(query), known, "tighter cluster");
            require(state.statistics().discovery_count == known.size(), "tightening a cluster lost its count proof");
        }
        {
            std::vector<Complex> known(7, Complex(0));
            known.insert(known.end(), 5, Complex(1, 2));
            known.insert(known.end(), 3, Complex(-2, 1));
            known.insert(known.end(), 3, Complex(3, -2));
            auto input = reference::input(from_roots(known));
            cps_options query;
            query.goal = CPS_ROOT_APPROXIMATION;
            d::PiecewiseRootStatistics stats;
            checked_cover(d::solve_prepared_piecewise_roots(certpolysolve::detail::dense(*input.get()), query, stats, configuration()), known);
            require(stats.discovery_count == known.size() && stats.final_accurate_count == known.size(),
                    "prepared ownership telemetry lost restored multiplicities");
        }
        std::cout << "ownership identity, cluster ledger, persistence, exact endpoints and prepared multiplicity passed\n";
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
