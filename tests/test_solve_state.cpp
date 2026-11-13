#include "numeric/arithmetic/polynomial_arithmetic.hpp"
#include "exact_winding_api.hpp"
#include "api_test_support.hpp"
#include "internal.hpp"
#include "methods/subdivision/subdivision.hpp"
#include <chrono>
#include <future>
#include <iostream>
#include <random>

using namespace reference;
namespace detail = certpolysolve::detail;
static uint64_t checked = 0;
static void fair_component_cells() {
    // The first component must split; the other can raise precision and refine
    // without allocating children. A cell stop must wait for both components.
    std::vector<detail::Root> cover(2);
    acb_zero(cover[0].box.value());
    arb_add_error_2exp_si(acb_realref(cover[0].box.value()), 1);
    arb_add_error_2exp_si(acb_imagref(cover[0].box.value()), 1);
    cover[0].count = 2;
    acb_set_si(cover[1].box.value(), 100);
    arb_add_error_2exp_si(acb_realref(cover[1].box.value()), 0);
    arb_add_error_2exp_si(acb_imagref(cover[1].box.value()), 0);
    cover[1].count = cover[1].multiplicity = 1;
    for (auto &r : cover) r.publish(64);
    certpolysolve::detail::SearchStatistics statistics{};
    detail::Subdivision subdivision(cover, 1, 64, statistics);
    cps_options options;
    options.bits = 300;
    options.max_work_bits = 512;
    const auto enclose = [](acb_poly_t p, slong) {
        acb_poly_zero(p);
        acb_poly_set_coeff_si(p, 0, 100);
        acb_poly_set_coeff_si(p, 1, -1);
        acb_poly_set_coeff_si(p, 2, -100);
        acb_poly_set_coeff_si(p, 3, 1);
    };
    bool progressed_beside_limit = false, cell_limit = false;
    for (int i = 0; i < 32; ++i) {
        const auto reason = subdivision.step(options, enclose, 0);
        std::vector<cps_component_progress> progress;
        subdivision.progress(progress, 0);
        const bool waiting = std::any_of(progress.begin(), progress.end(), [](const auto &c) {
            return c.state == CPS_COMPONENT_CELL_LIMIT;
        });
        progressed_beside_limit = progressed_beside_limit || (waiting &&
            std::any_of(progress.begin(), progress.end(), [](const auto &c) {
                return c.state == CPS_COMPONENT_FINISHED || c.work_bits > 128;
            }));
        if (reason == CPS_STOP_CELL_BUDGET) { cell_limit = true; break; }
        require(reason == CPS_STOP_NONE, "unexpected component stop");
    }
    require(cell_limit && progressed_beside_limit, "cell-limited component starved another component");
    for (int i = 0; i < 4000 && subdivision.pending(); ++i)
        require(subdivision.step(options, enclose, 65536) == CPS_STOP_NONE,
                "increased cell allowance failed to resume");
    std::vector<detail::Root> roots;
    require(!subdivision.pending() && subdivision.complete_cover(roots), "cell-limited search lost coverage");
    cps_result result;
    result.roots = std::move(roots); result.complete = result.goal_met = true;
    const auto published = certpolysolve_test::published(std::move(result));
    for (const auto &r : published.roots()) {
        unsigned count = 0;
        for (int root : {-1, 1, 100})
            count += endpoint(r.box.real.lo) <= root && root <= endpoint(r.box.real.hi) &&
                     endpoint(r.box.imag.lo) <= 0 && 0 <= endpoint(r.box.imag.hi);
        require(count == r.count, "resumed component returned an incorrect root count");
        check_accuracy(r, 300);
    }
}
static certpolysolve::Result snapshot(const cps_result &r) {
    return certpolysolve_test::published(r);
}
// Subdivision is the complete method, and production starts it at the
// precision ceiling. A failed Newton contraction must not stop it. With roots
// 16^-j, each component holding the smaller roots is separated from the next
// larger one, so it is certified. Its k-fold Newton step then fails for
// several levels per root (the roots span its disc), and failures reset only
// on a Newton success. Bisection still separates every root. The search must
// isolate all of them, each certified count checked against the known roots.
static void subdivision_bisects_past_newton_failures() {
    std::vector<Complex> roots;
    for (int j = 0; j < 8; ++j) roots.emplace_back(power2(-4 * j));
    const auto exact = from_roots(roots);
    std::vector<detail::Root> cover(1);
    acb_zero(cover[0].box.value());
    arb_add_error_2exp_si(acb_realref(cover[0].box.value()), 1);
    arb_add_error_2exp_si(acb_imagref(cover[0].box.value()), 1);
    cover[0].count = static_cast<uint64_t>(roots.size());
    cover[0].publish(64);
    cps_options options;
    options.bits = 53;
    options.max_work_bits = 128;
    certpolysolve::detail::SearchStatistics statistics{};
    detail::Subdivision subdivision(cover, 1, options.max_work_bits, statistics);
    const auto enclose = [&exact](acb_poly_t p, slong bits) {
        acb_poly_zero(p);
        detail::ComplexBall c;
        for (size_t k = 0; k < exact.size(); ++k) {
            detail::Number(exact[k].real.str()).enclose(acb_realref(c.value()), bits);
            detail::Number(exact[k].imag.str()).enclose(acb_imagref(c.value()), bits);
            acb_poly_set_coeff_acb(p, static_cast<slong>(k), c.value());
        }
    };
    cps_stop_reason reason = CPS_STOP_NONE;
    for (int i = 0; i < 20000 && subdivision.pending() && reason == CPS_STOP_NONE; ++i)
        reason = subdivision.step(options, enclose, 65536);
    require(reason == CPS_STOP_NONE && !subdivision.pending(), "subdivision stopped before isolating the roots");
    std::vector<detail::Root> found;
    require(subdivision.complete_cover(found), "subdivision lost its complete cover");
    require(found.size() == roots.size(), "subdivision did not isolate every root");
    for (const auto &r : found) {
        cps_result single;
        single.roots = {r};
        single.complete = single.goal_met = true;
        const auto published = certpolysolve_test::published(std::move(single));
        const auto box = published.roots().front();
        uint64_t inside = 0;
        for (const auto &z : roots)
            inside += endpoint(box.box.real.lo) <= z.real && z.real <= endpoint(box.box.real.hi) &&
                      endpoint(box.box.imag.lo) <= z.imag && z.imag <= endpoint(box.box.imag.hi);
        require(inside == box.count, "subdivision certified a wrong count");
    }
}

static void verify_known(const certpolysolve::Result &result, const std::vector<Complex> &known,
                         bool require_goal, int bits) {
    require(result.complete(), "lost complete coverage during a solve");
    if (require_goal && !result.goal_met())
        throw std::runtime_error("solve goal unmet: " + result.message());
    const auto regions = result.roots();
    check_disjoint(regions);
    size_t total = 0;
    for (const auto &r : regions) {
        size_t count = 0;
        for (const auto &z : known)
            if (endpoint(r.box.real.lo) <= z.real && z.real <= endpoint(r.box.real.hi) &&
                endpoint(r.box.imag.lo) <= z.imag && z.imag <= endpoint(r.box.imag.hi))
                ++count;
        require(count == r.count, "independent known-root count disagrees with the solve");
        if (result.goal_met() && cps_result_goal(result.get()) == CPS_ROOT_ACCURACY)
            check_accuracy(r, bits);
        total += count;
        ++checked;
    }
    require(total == known.size(), "solve lost roots");
}
// A one-unit work budget returns the Cauchy cover as a complete answer, and a result owns its data: it outlives
// the input it was solved from. (Until 18 October 2026 this was a public session's lifetime and resume test; the
// 0.5 API has one-shot solves only: docs/session-state-design-note-20261017.md, the decision at its top.)
static void one_unit_cover_and_owned_results() {
    // Degree 5: Auto's closed forms answer degree 3 and below within the unit.
    const std::vector<Complex> known{{-2}, {Rational(1) / 3}, {2, 1}, {Rational(5) / 2}, {-1, -1}};
    certpolysolve::Options o;
    o.work_budget = 1;
    o.cell_budget = 65536;
    o.bits = 53;
    o.max_bits = 128;
    const auto first = certpolysolve::solve(input(from_roots(known)), o);
    require(cps_result_termination(first.get()) == CPS_WORK_LIMIT && first.complete(), "a one-unit solve did not return Cauchy coverage");
    auto survivor = [&] {
        const auto temporary = input(from_roots(known));
        return certpolysolve::solve(temporary, certpolysolve::Options{});
    }();
    verify_known(survivor, known, true, 53);
}
// A tighter goal must not hand refined nodes back to binary64. Four real
// roots 1000 + k 2^-66 / 3 (k = 1..4), closer together than binary64 resolves
// at 1000 (its spacing there is 2^-43), and one at 2^-50 / 3: separating them
// for the 53-bit goal takes the secular iteration's precision past binary64's
// initial stage. Raised to 1,024 bits, that stage's bound (the goal's bits
// plus 32) rose past the refined nodes, and a binary64 packet rounded the
// cluster back onto 1000: about 150 work units to separate it again, where a
// fresh 1,024-bit solve needs 13 and the raise now needs about 20.
// A session tightened twice whose roots the certificate attempt keeps for the new goal (17 October 2026): on
// (x - 1)((x - 1)^2 + 9) the calls at 16 and 32 bits answer; at 64 bits the attempt certified all three roots for the
// new goal and left no node active, the secular model, unchanged, was taken for the same model and never validated
// again, and the call ran empty packets to its work budget (49,999 of them for 100,000 units), where a fresh 64-bit
// solve takes three units. Since 18 October 2026 a solve state is advanced once, so no later call tightens a
// model: every goal of the ladder 16, 32, ..., 4096 bits is a fresh solve that must answer within a small
// budget. The fix (secular_refinement.cpp, model_bits_.clear() at a goal change) stays, since a solve still
// changes its own goal (the staged goal, separation rounds), but no test reaches it any more: reverted on
// 18 October, only the earlier version of this test failed.
static void tightened_goal_validates_kept_certificates() {
    certpolysolve::Polynomial cubic;
    cubic.set(0, certpolysolve::Number("-10"));
    cubic.set(1, certpolysolve::Number("12"));
    cubic.set(2, certpolysolve::Number("-3"));
    cubic.set(3, certpolysolve::Number("1"));
    for (uint64_t bits = 16; bits <= 4096; bits *= 2) {
        certpolysolve::Options options;
        options.goal = certpolysolve::Goal::accuracy();
        options.bits = bits;
        options.max_bits = 4096;
        options.work_budget = 1000;
        options.cell_budget = 65536;
        const auto result = certpolysolve::solve(cubic, options);
        require(result.goal_met() && cps_result_termination(result.get()) == CPS_SUCCESS,
                ("a " + std::to_string(bits) + "-bit solve did not answer within 1,000 work units").c_str());
        require(result.roots().size() == 3, "a solve of the cubic lost a root");
    }
}
// Since 18 October 2026 (one advance a solve state) the same transition happens inside one solve: the secular
// route stages the default goal first, then the caller's. A fresh 1,024-bit solve must not spend its budget
// re-separating the cluster after the staged 53 bits.
static void tighter_goal_keeps_refined_nodes() {
    std::vector<Complex> known{{Rational(1) / 3 * power2(-50)}};
    for (int k = 1; k <= 4; ++k) known.emplace_back(1000 + Rational(k) / 3 * power2(-66));
    const auto p = input(from_roots(known));
    certpolysolve::Options o;
    o.method = certpolysolve::Method::secular;
    o.goal = certpolysolve::Goal::accuracy();
    o.bits = 53;
    o.max_bits = 4096;
    verify_known(certpolysolve::solve(p, o), known, true, 53);
    o.bits = 1024;
    o.max_bits = 4096;
    o.work_budget = 64;
    o.cell_budget = 65536;
    const auto refined = certpolysolve::solve(p, o);
    require(refined.goal_met(), "a 1,024-bit solve spent its work budget re-separating refined nodes");
    verify_known(refined, known, true, 1024);
}
// The forced subdivision path: a small budget stops honestly with a correct partial cover, and the solve with
// the 100,000 units the small calls of a session used to share meets the goal (one advance a solve state since
// 18 October 2026; until then one state was resumed call by call).
static void forced(const std::vector<Complex> &known, uint64_t quota = 100000) {
    auto p = input(from_roots(known));
    cps_options o;
    o.bits = 53;
    o.max_work_bits = 1024;
    if (quota < 100000) {
        o.work_units = quota;
        detail::SolveState limited(certpolysolve::detail::dense(*p.get()), {.backend = false});
        const auto stopped = limited.advance(&o);
        if (!stopped.goal_met)
            require(stopped.termination == CPS_WORK_LIMIT, "forced subdivision under a small budget stopped otherwise");
        verify_known(snapshot(stopped), known, false, 53);
    }
    o.work_units = 100000;
    detail::SolveState state(certpolysolve::detail::dense(*p.get()), {.backend = false});
    auto r = state.advance(&o);
    if (!r.goal_met) {
        std::cerr << "forced roots:";
        for (const auto &z : known)
            std::cerr << " (" << z.real << "," << z.imag << ")";
        std::cerr << " splits=" << r.statistics_v2.boxes_split
                  << " excluded=" << r.statistics_v2.boxes_excluded
                  << " newton=" << r.statistics_v2.newton_successes << "/" << r.statistics_v2.newton_attempts
                  << '\n';
    }
    verify_known(snapshot(r), known, true, 53);
    require(r.statistics_v2.backend_calls == 0, "forced subdivision used backend numerical roots");
    if (r.statistics.exact_factors < r.roots.size())
        require(r.statistics_v2.disk_tests > 0, "forced path did not certify disk counts");
}
static void cancellation_and_cells() {
    // A cell limit, a cancellation at the seventh poll, and the whole solve: three fresh solves (one advance a solve
    // state since 18 October 2026).
    const std::vector<Complex> known{{-2}, {-1}, {1}, {2}};
    auto p = input(from_roots(known));
    cps_options o;
    o.live_cells = 1;
    detail::SolveState cell_limited(certpolysolve::detail::dense(*p.get()), {.backend = false});
    auto limited = cell_limited.advance(&o);
    require(limited.termination == CPS_WORK_LIMIT && limited.progress.stop_reason == CPS_STOP_CELL_BUDGET,
            "cell limit not explicit");
    verify_known(snapshot(limited), known, false, 53);
    int calls = 0;
    o.live_cells = 65536;
    o.cancel = [](void *p) { return ++*static_cast<int *>(p) == 7 ? 1 : 0; };
    o.cancel_data = &calls;
    detail::SolveState cancelled(certpolysolve::detail::dense(*p.get()), {.backend = false});
    auto stopped = cancelled.advance(&o);
    require(stopped.termination == CPS_CANCELLED, "one-shot cancellation was not latched");
    verify_known(snapshot(stopped), known, false, 53);
    require(calls == 7, "the cancellation callback was polled after it cancelled");
    detail::SolveState whole(certpolysolve::detail::dense(*p.get()), {.backend = false});
    auto done = whole.advance(nullptr);
    verify_known(snapshot(done), known, true, 53);
    require(!done.components.empty(), "missing structured component progress");
}
static void algebraic_and_uncertain() {
    certpolysolve::Polynomial p{-2, 0, 1};
    detail::SolveState state(certpolysolve::detail::dense(*p.get()), {.backend = false});
    cps_options o;
    auto r = state.advance(&o);
    auto result = snapshot(r);
    require(result.goal_met(), "unhinted sqrt2 did not converge");
    check_disjoint(result.roots());
    for (const auto &root : result.roots()) {
        require(roots_closed({-2, 0, 1}, endpoint(root.box.real.lo), endpoint(root.box.real.hi)) == 1,
                "independent Sturm count failed");
        check_accuracy(root, 53);
        ++checked;
    }
    certpolysolve::Polynomial quartic{1, 1, 0, 0, 1};
    detail::SolveState complex_algebraic(certpolysolve::detail::dense(*quartic.get()), {.backend = false});
    auto complex_result = snapshot(complex_algebraic.advance(&o));
    require(complex_result.goal_met(), "unhinted algebraic complex roots did not converge");
    check_disjoint(complex_result.roots());
    uint64_t algebraic_count = 0;
    for (const auto &root : complex_result.roots()) {
        const auto count =
            rectangle_count({{1}, {1}, {0}, {0}, {1}}, endpoint(root.box.real.lo), endpoint(root.box.real.hi),
                            endpoint(root.box.imag.lo), endpoint(root.box.imag.hi));
        require(count && *count == static_cast<int>(root.count), "algebraic complex winding count failed");
        check_accuracy(root, 53);
        algebraic_count += root.count;
        ++checked;
    }
    require(algebraic_count == 4, "algebraic complex coverage count failed");
    certpolysolve::Polynomial family;
    family.set(0, certpolysolve::ball(certpolysolve::Number("-1"), certpolysolve::Number("1/100")));
    family.set(2, certpolysolve::Number("1"));
    detail::SolveState uncertain(certpolysolve::detail::dense(*family.get()), {.backend = false});
    o.max_work_bits = 128;
    o.work_units = 2000;
    auto u = snapshot(uncertain.advance(&o));
    require(u.complete() && !u.goal_met(), "fixed uncertainty was discarded");
    for (const Rational constant : {Rational(-101) / 100, Rational(-1), Rational(-99) / 100}) {
        for (const auto &root : u.roots()) {
            if (root.realness == certpolysolve::Realness::real)
                require(roots_closed({constant, 0, 1}, endpoint(root.box.real.lo),
                                     endpoint(root.box.real.hi)) == static_cast<int>(root.count),
                        "uncertainty sample Sturm count failed");
            else {
                auto count = rectangle_count({{constant}, {0}, {1}}, endpoint(root.box.real.lo),
                                             endpoint(root.box.real.hi), endpoint(root.box.imag.lo),
                                             endpoint(root.box.imag.hi));
                require(count && *count == static_cast<int>(root.count),
                        "uncertainty sample winding count failed");
            }
            ++checked;
        }
    }
}
// The default work budget grows with the degree (call_work_units,
// since 7 October 2026): 64 units a root, never fewer than the fixed 100,000; a
// budget the caller sets is kept exactly, whatever the degree.
static void default_work_budget_scales() {
    cps_options defaults;
    require(detail::call_work_units(defaults, 1000) == 100000, "a small input lost the fixed default budget");
    require(detail::call_work_units(defaults, 8192) == 64 * 8192, "the default budget did not grow with the degree");
    cps_options set;
    require(cps_options_set_work_budget(&set, 100000) == CPS_OK, "work budget rejected");
    require(detail::call_work_units(set, 8192) == 100000, "an explicit budget was changed");
    require(cps_options_set_work_budget(&set, 5) == CPS_OK, "work budget rejected");
    require(detail::call_work_units(set, 1 << 20) == 5, "an explicit small budget was changed");
}
int main(int argc, char **argv) {
    certpolysolve::ThreadScope scope;
    try {
        const uint64_t count = argc == 3 && std::string(argv[1]) == "--cases" ? std::stoull(argv[2]) : 32;
        default_work_budget_scales();
        one_unit_cover_and_owned_results();
        tighter_goal_keeps_refined_nodes();
        tightened_goal_validates_kept_certificates();
        cancellation_and_cells();
        fair_component_cells();
        subdivision_bisects_past_newton_failures();
        algebraic_and_uncertain();
        forced({{-1}, {0, 1}, {1}, {0, -1}}, 1);
        forced({{0}, {1}, {2}, {-1}}, 1);
        forced({{1}, {1 + power2(-200)}, {2, 1}});
        forced({{power2(-5000)}, {-power2(-5000)}});
        forced({{power2(5000)}, {-power2(5000)}});
        forced({{0, 1}, {0, 1}, {0, -1}, {0, -1}, {0, -1}});
        const Rational first = Rational(1) / 3, second = first + power2(-200);
        forced({{first}, {first}, {second}, {second}, {second}});
        std::vector<std::future<void>> workers;
        for (int i = 0; i < 8; ++i)
            workers.push_back(std::async(std::launch::async, [] {
                certpolysolve::ThreadScope worker_scope;
                const certpolysolve::Polynomial p{-2, 0, 1};
                const auto low = certpolysolve::solve(p, certpolysolve::Options{});
                certpolysolve::Options high;
                high.bits = 256;
                const auto refined = certpolysolve::solve(p, high);
                require(low.goal_met() && refined.goal_met(), "concurrent solves failed");
            }));
        for (auto &worker : workers)
            worker.get();
        std::mt19937_64 generator(2026090901);
        for (uint64_t i = 0; i < count; ++i) {
            std::vector<Complex> roots;
            const size_t degree = 2 + generator() % 5;
            for (size_t j = 0; j < degree; ++j)
                roots.emplace_back(Rational(static_cast<int>(generator() % 21) - 10) / 8,
                                   Rational(static_cast<int>(generator() % 21) - 10) / 8);
            try {
                forced(roots, i % 3 == 0 ? 1 : 100000);
            } catch (const std::exception &e) {
                std::cerr << "solve_state seed=2026090901 index=" << i << '\n';
                throw;
            }
            if (i % 1000 == 999)
                std::cout << "cases=" << i + 1 << " checked_regions=" << checked << std::endl;
        }
        std::cout << "solve_state/subdivision cases=" << count << " independent_regions=" << checked
                  << " passed\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
