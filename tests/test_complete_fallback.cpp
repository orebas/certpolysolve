#include "exact_winding_api.hpp"
#include "internal.hpp"
#include "numeric/arithmetic/polynomial_arithmetic.hpp"
#include "solve/local_model_solver.hpp"
#include "support/fuzz/fuzz_case.hpp"
#include <algorithm>
#include <iostream>
#include <set>
#include <string_view>

namespace d = certpolysolve::detail;
using reference::Complex;
using reference::Rational;
using reference::require;

namespace {
struct Transition {
    uint64_t ticks = 0, ceiling = 0, sweeps = 0, regenerations = 0;
    std::string reason;
};
struct Observed {
    std::vector<Transition> transitions;
    unsigned cancel_polls = 0, cancelling_polls = 0;
    bool cancel_after_transition = false;
    static void trace(void *data, const d::TraceEvent &event) {
        if (std::string_view(event.stage) != "complete_fallback" ||
            std::string_view(event.phase) != "transition") return;
        auto &self = *static_cast<Observed *>(data);
        self.transitions.push_back({event.algorithm_cutoff, event.exact_budget,
            event.candidates ? event.candidates->candidate_sweeps : 0,
            event.candidates ? event.candidates->secular_regenerations : 0, event.reason});
    }
    static int cancel(void *data) {
        auto &self = *static_cast<Observed *>(data);
        ++self.cancel_polls;
        if (self.cancel_after_transition && !self.transitions.empty())
            return ++self.cancelling_polls == 1; // deliberately one-shot
        return 0;
    }
};

d::EnginePolicy policy(Observed &observed, uint64_t packets) {
    auto result = d::complete_global_policy(CPS_METHOD_SECULAR);
    result.heuristic_packets = packets;
    result.heuristic_precision_multiplier = 2;
    result.trace = Observed::trace;
    result.trace_data = &observed;
    return result;
}

cps_options options(cps_goal goal = CPS_ROOT_APPROXIMATION) {
    cps_options result;
    result.method = CPS_METHOD_SECULAR;
    result.goal = goal;
    result.bits = 53;
    result.max_work_bits = 4096;
    result.work_units = 100000;
    result.work_units_default = false;
    result.live_cells = 4096;
    result.temporary_bytes = uint64_t(64) << 20;
    return result;
}

reference::Rectangle rectangle(const d::Root &root) {
    reference::Rectangle result;
    Rational *ends[] = {&result.real_lower, &result.real_upper, &result.imag_lower, &result.imag_upper};
    for (size_t i = 0; i < 4; ++i) {
        const auto text = root.endpoints[i].text();
        *ends[i] = reference::dyadic(text.mantissa, text.exponent);
    }
    result.count = root.count;
    return result;
}

void check_published(const cps_result &result, const std::vector<Complex> &roots, bool success,
                     cps_goal expected_goal = CPS_ROOT_APPROXIMATION, unsigned bits = 53) {
    require(result.degree == static_cast<int64_t>(roots.size()), "published degree changed");
    require(result.goal == expected_goal, "returned contract differs from the requested goal");
    if (success)
        require(result.complete && result.goal_met && result.termination == CPS_SUCCESS,
                "complete fallback did not attain the requested goal");
    std::vector<unsigned> ownership(roots.size(), 0);
    std::vector<reference::Rectangle> boxes;
    for (const auto &region : result.roots) {
        const auto box = rectangle(region);
        uint64_t count = 0;
        std::set<std::pair<Rational, Rational>> distinct;
        for (size_t i = 0; i < roots.size(); ++i)
            if (box.real_lower <= roots[i].real && roots[i].real <= box.real_upper &&
                box.imag_lower <= roots[i].imag && roots[i].imag <= box.imag_upper) {
                ++ownership[i];
                ++count;
                distinct.emplace(roots[i].real, roots[i].imag);
                if (result.goal_met) {
                    const Rational dr = (box.real_lower + box.real_upper) / 2 - roots[i].real;
                    const Rational di = (box.imag_lower + box.imag_upper) / 2 - roots[i].imag;
                    const Rational norm = roots[i].real * roots[i].real + roots[i].imag * roots[i].imag;
                    require(dr * dr + di * di <= reference::power2(-2 * static_cast<int>(bits)) * norm,
                            "published midpoint misses an independently known root");
                }
            }
        require(count == region.count && count > 0, "published count disagrees with independent exact roots");
        if (result.goal_met) {
            if (result.goal == CPS_ROOT_ACCURACY) {
                reference::check_accuracy(box, bits);
                require(distinct.size() == 1 && region.multiplicity == count,
                        "individual accuracy was replaced by a counted cluster");
            }
            // Counted approximation also permits the Euclidean predicate;
            // the exact per-root inequality above checks its actual promise.
        }
        boxes.push_back(box);
    }
    reference::check_disjoint(boxes);
    for (const auto count : ownership)
        require(count <= 1 && (!result.complete || count == 1), "published cover lost unique root ownership");
}

std::vector<Complex> complex_roots() {
    return {{-2, 1}, {-1, -1}, {1, 2}, {2, -2}, {3, 1}, {4, -1}};
}

struct Run {
    cps_result result;
    std::optional<d::PreparedStructure> prepared;
};
Run solve(const std::vector<Complex> &roots, cps_options opts, d::EnginePolicy recipe) {
    const auto polynomial = reference::input(reference::from_roots(roots));
    d::SolveState state(d::dense(*polynomial.get()), std::move(recipe));
    auto result = state.advance(&opts);
    return {std::move(result), state.prepared_structure()};
}

void immediate_entry() {
    const Rational near = reference::power2(-72);
    const std::vector<std::vector<Complex>> fixtures{
        complex_roots(),
        {{-2}, {-2}, {-2}, {1}, {1}, {3}}, // distinct square-free factors/multiplicities
        {{1}, {-1}, {0, 1}, {0, -1}, {2}, {-2}, {0, 2}, {0, -2}}, // q(x^4)
        {{0}, {0}, {near}, {near}, {2}, {-3}} // extracted zero near a different factor
    };
    for (size_t fixture = 0; fixture < fixtures.size(); ++fixture)
        for (const auto goal : {CPS_ROOT_APPROXIMATION, CPS_ROOT_ACCURACY}) {
            Observed observed;
            auto opts = options(goal);
            opts.cancel = Observed::cancel;
            opts.cancel_data = &observed;
            const auto done = solve(fixtures[fixture], opts, policy(observed, 0));
            check_published(done.result, fixtures[fixture], true, goal);
            require(!observed.transitions.empty() && observed.cancel_polls > 0,
                    "K=0 did not force complete entry or exercise caller polling");
            for (const auto &transition : observed.transitions)
                require(transition.ticks == 0 && transition.sweeps == 0 && transition.regenerations == 0 &&
                            transition.reason == "heuristic-packets",
                        "K=0 performed heuristic factor work before complete entry");
            require(done.result.statistics_v3.candidate_sweeps == 0 &&
                        done.result.statistics_v3.secular_regenerations == 0,
                    "mandatory service continued candidate iteration");
            require(observed.cancelling_polls == 0 && done.result.progress.stop_reason == CPS_STOP_NONE,
                    "internal budget transition latched caller cancellation");
            require(done.prepared && done.prepared->accounted_degree() == fixtures[fixture].size(),
                    "prepared factors lost the original degree");
            if (fixture == 1)
                require(done.prepared->factors.size() >= 3 && observed.transitions.size() >= 3,
                        "multiplicity fixture did not exercise multiple complete workers");
            if (fixture == 2)
                require(std::any_of(done.prepared->factors.begin(), done.prepared->factors.end(),
                                    [](const auto &factor) { return factor.power > 1; }),
                        "power-map fixture was not reduced and lifted");
            if (fixture == 3)
                require(done.prepared->zero_order == 2, "zero-near-factor fixture did not extract its zero roots");
        }
}

void bounded_prefix_without_event() {
    std::vector<Complex> roots;
    for (int i = 0; i < 6; ++i)
        roots.emplace_back(Rational(1) / 3 + i * reference::power2(-96), Rational(i % 2 + 1) / 7);
    for (const uint64_t packets : {1u, 2u}) {
        Observed observed;
        auto recipe = policy(observed, packets);
        recipe.research_dovetail_start = d::EnginePolicy::DovetailStart::Never;
        auto opts = options(CPS_ROOT_ACCURACY);
        const auto done = solve(roots, opts, std::move(recipe));
        check_published(done.result, roots, true, CPS_ROOT_ACCURACY);
        require(!observed.transitions.empty() && done.result.statistics_v3.candidate_sweeps > 0,
                "small positive K did not exercise a numerical prefix followed by complete service");
        require(observed.transitions.front().reason == "heuristic-packets" &&
                    observed.transitions.front().ticks >= packets * roots.size(),
                "mandatory service depended on a practical precision trigger or renewed its allowance");
    }
}

void caller_limits() {
    const auto roots = complex_roots();
    for (unsigned limit = 0; limit < 4; ++limit) {
        Observed observed;
        auto opts = options();
        cps_stop_reason expected = CPS_STOP_WORK_BUDGET;
        cps_termination termination = CPS_WORK_LIMIT;
        if (limit == 0) opts.work_units = 0;
        if (limit == 1) opts.work_units = 3; // degree proof, preparation, mandatory transition
        if (limit == 2) { opts.live_cells = 0; expected = CPS_STOP_CELL_BUDGET; }
        if (limit == 3) { opts.max_work_bits = 32; expected = CPS_STOP_PRECISION; termination = CPS_PRECISION_LIMIT; }
        const auto stopped = solve(roots, opts, policy(observed, 0));
        require(!stopped.result.goal_met && stopped.result.termination == termination &&
                    stopped.result.progress.stop_reason == expected,
                "complete fallback converted a caller resource limit into success or a different stop");
        require(stopped.result.progress.work_units_this_call <= opts.work_units &&
                    stopped.result.progress.max_outer_bits <= opts.max_work_bits,
                "mandatory transition bypassed a caller work or precision cap");
        check_published(stopped.result, roots, false);
        if (limit == 1) require(!observed.transitions.empty(), "work cap fixture did not reach mandatory transition");
    }
    Observed observed;
    observed.cancel_after_transition = true;
    auto opts = options();
    opts.cancel = Observed::cancel;
    opts.cancel_data = &observed;
    const auto cancelled = solve(roots, opts, policy(observed, 0));
    require(!cancelled.result.goal_met && cancelled.result.termination == CPS_CANCELLED &&
                cancelled.result.progress.stop_reason == CPS_STOP_CANCELLED &&
                observed.transitions.size() == 1 && observed.cancelling_polls == 1,
            "one-shot external cancellation after transition was lost or polled after latching");
    check_published(cancelled.result, roots, false);
}

void borrowed_budget_does_not_renew() {
    const auto roots = complex_roots();
    d::HeuristicBudget precision(roots.size(), 32, 53, 16, 1);
    auto precision_options = options();
    precision_options.heuristic_budget = &precision;
    precision.observe_precision(precision.ceiling());
    require(!precision.exhausted() && precision.peak() == precision.ceiling(),
            "touching Q expired the parent's unspent packet allowance");
    precision.finish_at_precision_limit();
    require(precision.exhausted() && precision.stop() == d::HeuristicBudget::Stop::Precision &&
                precision.ticks() == 0 && !precision_options.cancelled(),
            "terminal Q stop was confused with spent packets or caller cancellation");
    d::HeuristicBudget outer(roots.size(), 32, 53, 1, 2);
    auto parent = options();
    parent.heuristic_budget = &outer;
    auto child = parent;
    child.heuristic_budget->linear();
    require(!parent.heuristic_budget->exhausted(), "one linear action spent a full-size packet");
    child.heuristic_budget->quadratic(); // finish the admitted packet, including its charged overrun
    require(outer.ticks() == roots.size() + 1 && outer.exhausted() && !parent.cancelled(),
            "shared soft overrun was lost or became cancellation");
    const auto ticks = outer.ticks();
    const auto polynomial = reference::input(reference::from_roots(roots));
    child.exact_budget = 0; // maintained local factory's ordinary preparation contract
    const auto local = d::solve_local_model_priced(d::dense(*polynomial.get()), child);
    require(!local.goal_met && local.termination == CPS_WORK_LIMIT &&
                local.progress.stop_reason == CPS_STOP_WORK_BUDGET &&
                local.statistics_v3.candidate_sweeps == 0 && outer.ticks() == ticks,
            "a local solver renewed or ignored its borrowed exhausted heuristic allowance");
    check_published(local, roots, false);
    Observed observed;
    const auto sibling = solve(roots, parent, policy(observed, 1024));
    require(!sibling.result.goal_met && sibling.result.progress.stop_reason == CPS_STOP_WORK_BUDGET &&
                observed.transitions.empty() && outer.ticks() == ticks,
            "a borrowing state created its own renewed allowance or complete worker");
    check_published(sibling.result, roots, false);
}

void capped_child_preserves_parent_allowance() {
    const auto roots = complex_roots();
    d::HeuristicBudget outer(roots.size(), 32, 53, 1024, 1);
    auto child = options();
    child.bits = outer.ceiling() + 32;
    child.max_work_bits = outer.ceiling();
    child.heuristic_budget = &outer;
    child.exact_budget = 0;
    const auto polynomial = reference::input(reference::from_roots(roots));
    const auto local = d::solve_local_model_priced(d::dense(*polynomial.get()), child);
    require(outer.ticks() > 0 && outer.peak() == outer.ceiling() && !outer.exhausted(),
            "a capped optional local solve expired the outer query's remaining allowance");
    require(local.progress.max_outer_bits <= child.max_work_bits,
            "a child escaped its admitted precision ceiling");
    check_published(local, roots, false, child.goal, static_cast<unsigned>(child.bits));
}

void precision_transition() {
    std::vector<Complex> roots;
    for (int i = 0; i < 6; ++i)
        roots.emplace_back(Rational(1) / 3 + i * reference::power2(-96), Rational(i % 2 + 1) / 7);
    Observed observed;
    auto recipe = policy(observed, 1024);
    // Test override: the minimum 64-bit heuristic ceiling. The caller still
    // allows 4096 bits, which must remain available to complete service.
    recipe.heuristic_precision_multiplier = 0;
    recipe.research_dovetail_start = d::EnginePolicy::DovetailStart::Never;
    const auto done = solve(roots, options(CPS_ROOT_ACCURACY), std::move(recipe));
    check_published(done.result, roots, true, CPS_ROOT_ACCURACY);
    require(!observed.transitions.empty() && observed.transitions.front().reason == "heuristic-precision" &&
                observed.transitions.front().ceiling == 64 && done.result.progress.max_outer_bits > 64,
            "heuristic Q became the complete worker's precision cap or waited for the practical trigger");
}

void terminal_handoff_at_caller_cap() {
    // Four exact geometries from the standing axis fuzzer: an improved final
    // numerical cover, retained-precision refinement, a zero-metered warm
    // state, and a warm worker reaching its guard before numerical iteration.
    // None needs the mandatory K/Q transition or a larger caller allowance.
    const auto p2 = reference::power2;
    const Rational a = p2(31), c = Rational(-92987, 131072), h(1, 2), t(64969, 65536);
    struct Fixture {
        std::vector<Complex> roots;
        unsigned bits;
        uint64_t cap, work, cells;
        cps_goal goal;
    };
    const std::vector<Fixture> fixtures{
        {{{a, p2(-480)}, {a, -p2(-480)}, {a}, {7823360}, {a + p2(-192)},
          {a + p2(-222)}, {a - p2(-160)}, {a}, {a}, {a}},
         100, 512, 100, 256, CPS_ROOT_APPROXIMATION},
        {{{-p2(-20), c}, {p2(-20), c}, {-p2(-38)}, {p2(-38)}, {0, c},
          {0, Rational(502701, 1048576)}, {0, c}, {0, c}},
         100, 96, 100, 16, CPS_ROOT_ACCURACY},
        {{{h, p2(-66)}, {h, -p2(-66)}, {Rational(6055, 8192)}, {h + p2(-48)},
          {h}, {h}, {h}, {h - p2(-59)},
          {Rational(68271, 131072), Rational(-966371, 1048576)},
          {Rational(68271, 131072), Rational(966371, 1048576)},
          {Rational(43, 25), Rational(3, 8)}, {Rational(43, 25), Rational(-3, 8)}},
         16, 80, 100000, 256, CPS_ROOT_APPROXIMATION},
        {{{Rational(15, 16), p2(-42)}, {Rational(15, 16), -p2(-42)},
          {t, p2(-50)}, {t, -p2(-50)}, {t}, {Rational(141, 256)}, {t},
          {Rational(259877, 262144)}, {t}},
         16, 80, 100000, 65536, CPS_ROOT_APPROXIMATION}
    };
    for (const auto &fixture : fixtures) {
        auto opts = options(fixture.goal);
        opts.method = CPS_METHOD_AUTO;
        opts.bits = fixture.bits;
        opts.max_work_bits = fixture.cap;
        opts.work_units = fixture.work;
        opts.live_cells = fixture.cells;
        Observed observed;
        auto recipe = d::complete_global_policy(CPS_METHOD_SECULAR);
        recipe.trace = Observed::trace;
        recipe.trace_data = &observed;
        const auto done = solve(fixture.roots, opts, std::move(recipe));
        check_published(done.result, fixture.roots, true, fixture.goal, fixture.bits);
        require(observed.transitions.empty() && done.result.progress.max_outer_bits <= fixture.cap,
                "terminal handoff ignored the caller cap or relied on mandatory complete service");
    }
    // This warm search reaches the trial guard only after the final numerical
    // handoff. Reuse the versioned exact generator rather than duplicating its
    // seventy planted roots and coefficient expansion.
    fuzz::high_degree_family = true;
    const auto late = fuzz::generate(UINT64_C(3695115102503025219));
    fuzz::high_degree_family = false;
    require(late.roots.size() == 70 && late.bits == 16 && late.max_work_bits == 512,
            "the fixed late-guard regression fixture changed");
    auto opts = options(late.goal);
    opts.method = CPS_METHOD_AUTO;
    opts.bits = late.bits;
    opts.max_work_bits = late.max_work_bits;
    opts.work_units = late.work_units;
    opts.live_cells = late.live_cells;
    Observed observed;
    auto recipe = d::complete_global_policy(CPS_METHOD_SECULAR);
    recipe.trace = Observed::trace;
    recipe.trace_data = &observed;
    const auto done = solve(late.roots, opts, std::move(recipe));
    check_published(done.result, late.roots, true, late.goal, late.bits);
    require(observed.transitions.empty() && done.result.progress.max_outer_bits <= late.max_work_bits,
            "late trial guard used a larger cap or mandatory complete service");
}

#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
void bounded_state_cannot_escape_into_independent_fork() {
    Observed observed;
    const auto input = reference::input(reference::from_roots(complex_roots()));
    d::SolveState state(d::dense(*input.get()), policy(observed, 1024));
    state.research_start(options());
    for (unsigned step = 0; step < 4 && !state.research_checkpoint().prepared; ++step)
        state.research_step();
    require(state.research_checkpoint().prepared, "fork-lifetime fixture did not prepare its owned budget");
    bool refused = false;
    try { (void)state.research_warm_fork(4096); }
    catch (const std::logic_error &) { refused = true; }
    require(refused, "an independently owned fork retained a pointer to its parent's heuristic budget");
}
#endif
} // namespace

int main() {
    certpolysolve::ThreadScope lifetime;
    try {
        immediate_entry();
        bounded_prefix_without_event();
        caller_limits();
        borrowed_budget_does_not_renew();
        capped_child_preserves_parent_allowance();
        precision_transition();
        terminal_handoff_at_caller_cap();
#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
        bounded_state_cannot_escape_into_independent_fork();
#endif
        std::cout << "complete fallback: soft boundary, exact contracts, caller limits and shared budget passed\n";
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
