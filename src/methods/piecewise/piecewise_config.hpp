#pragma once
// The piecewise recipe, its precision ladder entry and observer contracts. These
// describe a solve without depending on its controller, models or statistics layout.
#include "core/cover.hpp"
#include "core/resource_limits.hpp"
#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>

namespace certpolysolve::detail {
struct PiecewiseRootStatistics;
struct PiecewiseRootSnapshot;
struct PiecewisePopulationSnapshot;
struct RootModel;
struct PiecewiseTransformCapture;
struct PiecewiseTransformEvent {
    uint64_t length = 0, model_bits = 0, arithmetic_bits = 0, row = 0, execution_ns = 0;
};
struct PiecewiseCandidateEvent {
    uint64_t candidate = UINT64_MAX, model_bits = 0, arithmetic_bits = 0;
    bool carried = false, model_stopped = false, newly_certified = false, improved = false;
    int failure = 0;
    double error_log2 = 0, margin_log2 = 0, correction_log2 = 0;
    bool component_terms_available = false;
    std::array<double, 5> component_log2{}; // band, model tail, jet tail, arithmetic, constant
    uint64_t model = 0, local_degree = 0, owner = UINT64_MAX;
    uint64_t seed_ns = 0, final_ns = 0;
    // The failure as the certificate reported it, -1 when proved (the same as `failure` since W8's S2 deleted the one
    // reclassification, aligned_escalation; kept for the observers' records).
    int raw_failure = 0;
    uint64_t source = 0;
};
struct PiecewiseSelectionEvent {
    uint64_t source = 0, candidate = 0, model = 0, model_bits = 0, arithmetic_bits = 0;
    bool moved = false, indexed_in_disk = false, current_in_disk = false;
    double radius_log2 = 0, indexed_norm = 0, current_norm = 0;
};
struct PiecewisePolishEvent {
    uint64_t candidate = UINT64_MAX, model = 0, local_degree = 0;
    uint64_t model_bits = 0, arithmetic_bits = 0, evaluations = 0, nanoseconds = 0;
    // Read-only state at entry; zero means no active carried-owner proof.
    uint64_t owner_count = 0, source = 0;
    uint64_t model_stops = 0, rounding_stops = 0, step_stops = 0, reused = 0;
    uint64_t nonimproving_steps = 0, iteration_limit = 0, domain_stops = 0, derivative_stops = 0;
};
struct PiecewiseFinisherEvent {
    uint64_t source = 0, arithmetic_bits = 0, covered = 0;
    uint64_t owners = 0, position = 0;
    // Counted work (counted_work's unit): the cohort's budget, what it has
    // spent, and the parts of that.
    double budget = 0, spent = 0, pilot = 0, preparation = 0, step = 0;
    double admission = 0, marginal_pilot = 0;
    bool certificate = false;
};
struct PiecewiseRootConfig {
    // Explicit construction scratch allowance, independent of degree. The
    // caller may restrict it; unspecified means the addressable range.
    size_t model_storage_bytes = std::numeric_limits<size_t>::max();
    // Value initialization retains the explicit baseline for private ablations.
    // The prepared standalone entry uses the maintained configuration below.
    bool reuse_guards = false;
    bool recover_discovery = false;
    bool reuse_polish = false;
    bool conjugate_symmetry = false;
    // C1 prototype (docs/c1-axis-charts-note-20261015.md): build only the charts whose centre lies on the
    // real axis (sector 0 and S/2 of each ring); their rows' outputs as a sum and an alternating sum.
    bool axis_only = false;
    // Rounds an axis-mode solve may take before it stops short (0: no limit): Auto's attempt before the complex
    // route takes heuristics::axis_attempt_rounds, since on the corpus the axis mode answered faster than Auto
    // only in its first rounds. History: docs/code-history-20261021.md, piecewise_roots.hpp 1.
    uint64_t axis_rounds = 0;
    bool blocked_polish = false;
    bool detailed_accounting = false;
    bool error_components = false;
    bool local_band_error = false;
    bool aligned_polish = false;
    // Qualified original-coefficient block/disk refinement and certificate.
    // Historical representation/evaluator ablations live outside the core.
    bool source_finisher = false;
    // Exact p(x) = x^s product q_j(x^k_j)^m_j preparation and outward lifting;
    // on in the maintained recipe (maintained_piecewise_config).
    bool power_composition = false;
    bool durable_ownership = false;
    // A model's transform runs in fixed-point integers (fixed_transform.hpp) where a word of at most
    // FixedTransform::maximum_limbs limbs reaches the working precision; FLINT's ball transform otherwise. It stands in for the ball transform at the same precision: its error is at most one
    // unit in the last place of the largest input at that precision, by a bound known before it runs.
    bool fixed_fft = false;
    // Skip a local error calculation only when it cannot change polishing's
    // stopping predicate. Final certificates still calculate their own bound.
    bool reuse_polish_bound = false;
    // Degree and local remainder use the same proved Taylor-order provenance.
    // The analytic remainder always covers the actual omitted Taylor terms.
    bool adaptive_taylor = false;
    // Private experiment: once every root is owned, a model-limited
    // escalation skips the doubled rung unless enough of this round's
    // model-remainder refusals would pass there (round_shortfall_log2).
    // Building the doubled rung costs c(2m) for every refused root and then
    // c(4m) for those it cannot finish; skipping costs c(4m) for all. So the
    // doubled rung pays only if the fraction it finishes exceeds c(2m)/c(4m),
    // taken as the ratio of the two rungs' arithmetic precisions
    // (2 model bits + goal bits + 64 guard bits, as the ladder sets them).
    // At most one rung is skipped per round.
    bool predicted_refinement = false;
    // Leave piecewise for a global recipe once its counted work exceeds the
    // predicted cost of a global solve: regime_exit_units complex coefficient
    // operations (set by the prepared state from the whole input:
    // initial_sweep_allowance sweeps of degree x (degree + support), the
    // global recipes' minimum budget, as in auto_selection.cpp), at one of
    // two prices. No clock is read: the rule gives the same answer on a
    // loaded machine and an idle one. Checked before every work unit; the
    // owner of the state performs the handoff.
    //
    // regime_exit_estimate > 0 (Auto: the selector chose this route because
    // its estimate of piecewise's work, that many operations, is below the
    // budget): the state leaves once its counted work exceeds what the
    // estimate predicts (estimate_operation_overhead binary64 operations
    // each) plus the whole budget at its own price, a binary64 operation
    // each (the global sweeps begin in binary64). The overhead factor is a
    // measured calibration, not a universal upper bound. A solve within
    // that predicted allowance does not leave; an overrun is compared with
    // the global minimum budget. Work between checks can overshoot the
    // threshold, so this rule alone gives no whole-call overhead guarantee.
    //
    // The caller asked for piecewise (caller_method): the session asks the
    // selector as Auto would, and where it would choose piecewise the state
    // takes its estimate and leaves exactly as Auto's route does (above).
    // Where it would not, the state leaves once its counted work exceeds the
    // global route's minimum budget at a binary64 operation each: the
    // estimate's term is zero, the budget's price the same. Priced instead at
    // the state's own counted work per model-construction unit (as below), far
    // more, the budget would let a forced solve leave only far outside its
    // regime. History: docs/code-history-20261021.md, piecewise_roots.hpp 3.
    //
    // Otherwise (neither an estimate nor the caller's method: a direct caller
    // of the state, such as the piecewise worker or a test; the real-only
    // route's axis attempt has no regime exit, configure()): each operation at
    // this state's own counted work per model-construction unit.
    bool regime_exit = false;
    bool caller_method = false;
    // Both counts of operations, exact (saturating at 2^64 - 1).
    uint64_t regime_exit_units = 0, regime_exit_estimate = 0;
    // Seed the roots a round left undiscovered when it owned no new root or
    // when every owned root already meets the goal, the owned roots being the
    // majority: with every owned root held, (degree - owned) free points iterate on the
    // source by Aberth's correction, which is Newton's step on the source
    // divided by the held roots, so a free point cannot settle on a root
    // already owned. The converged points join the carried pool and the next
    // round certifies them. Without it, discovery rests on the first round's
    // binary64 proposals: in a cluster below their residual floor several
    // proposals stop at one spot, Newton polishing may send two to one root,
    // and a chart that cannot count its roots never proposes again
    // (clustered Mahler polynomials).
    // History: docs/code-history-20261021.md, piecewise_roots.hpp 4.
    bool deflated_discovery = false;
    // When a round's deflated discovery has seeded missing roots and a discovery recovery (a round that starts
    // every chart's discovery again) is due as well, the recovery waits one round: the seeds are polished and
    // certified first, and the recovery follows only if roots are still missing after it.
    bool seeds_before_recovery = false;
    // Begin with the certificate attempt the global recipes begin with
    // (machine_newton_cover, its binary64 packet only), on the squarefree
    // nonzero part, before any model is built. The prepared state makes it
    // once; refused, piecewise runs as without it. Set by the session when
    // the selector's cost comparison says so (auto_selection.cpp); a forced
    // piecewise solve does not make it.
    bool certificate_attempt = false;
    // What piecewise is predicted to cost on this input, in the selector's operations (its admission estimate,
    // auto_selection.cpp): the attempt before it may spend no more (machine_newton_cover's fallback_operations).
    // Zero: the global recipes' minimum budget.
    uint64_t certificate_attempt_fallback = 0;
    // Qualified singleton guard optimization; false retains the private control.
    // Refusals use value/slope first; successful proofs keep the complete jet.
    bool lazy_guard_jet = false;
    // Search covered interior disks with a geometric buffer. Inadequate models
    // are strengthened before a deferred full-chart discovery retry.
    bool interior_proposals = false;
    // A chart whose disk holds no root of its polynomial (Pellet's test after
    // root-squaring steps, chart_excluded) runs no proposal iteration.
    bool chart_exclusion = false;
    // Private experiment: counted singleton coverage makes unowned discovery
    // seeds unnecessary. Preserve them for partial/cluster ownership states.
    bool owned_only_after_discovery = false;
    // Owned private diagnostic callback. Must not throw, re-enter or retain the
    // borrowed statistics reference. It does not participate in scheduling.
    std::function<void(const PiecewiseRootStatistics &, const char *)> observe;
    std::function<void(const PiecewiseCandidateEvent &)> observe_candidate;
    std::function<void(const PiecewiseSelectionEvent &)> observe_selection;
    std::function<void(PiecewisePopulationSnapshot)> observe_population;
    std::function<void(const PiecewisePolishEvent &)> observe_polish;
    // Read-only borrowed batches; an observer must own any retained copies.
    std::function<void(const RootModel &, const RootVector &, slong, const cps_options &, slong)> observe_polish_batch;
    std::function<void(const RootModel &, const Root &, const ComplexBall &, const cps_options &, slong)> observe_owned_refinement;
    // Read-only cost diagnostics; observation never supplies a decision.
    std::function<void(const PiecewiseFinisherEvent &, const char *)> observe_finisher;
    // Read-only borrows valid only during the callback. An observer retaining
    // a model/point must copy it; the model's error envelope owns its source.
    std::function<void(const RootModel &, const ComplexBall &, const cps_options &, slong, bool)> observe_guard;
    std::function<void(const RootModel &, const ComplexBall &, const cps_options &, slong, bool, int, bool)> observe_certificate;
    bool snapshot_discovery = false, snapshot_rounds = false;
    // Private diagnostics own every retained value. Callbacks must not re-enter
    // or mutate this solve. No snapshot state participates in numerical policy.
    uint64_t diagnostic_source_id = 0, diagnostic_multiplicity = 1;
    std::function<bool()> snapshot_due;
    std::function<void(PiecewiseRootSnapshot)> observe_snapshot;
    std::function<void(const PiecewiseTransformEvent &)> observe_transform;
    std::function<bool(const PiecewiseTransformEvent &)> want_transform;
    std::function<void(PiecewiseTransformCapture)> capture_transform;
};
// Qualified standalone recipe. Numerical proposals still require original-source
// certificates. Source finishing is admitted by measured construction cost;
// extra escalation, external FFTs and detailed timers stay opt-in. This does not
// change public automatic algorithm admission.
// The model bits a state starts with, and its first working precision for a
// goal: the ladder's base rung, max(128, 2 x model bits + goal bits + 64 guard
// bits). A precision cap below it is refused before any work; the router reads
// this so it does not send such a query to piecewise (queue item 3.2).
inline constexpr slong piecewise_initial_model_bits = 16;
inline slong piecewise_start_precision(uint64_t goal_bits, slong model_bits = piecewise_initial_model_bits) {
    return std::max<slong>(128, resources::precision_sum(resources::precision_sum(
        resources::precision_product(model_bits, 2), goal_bits), 64));
}
inline PiecewiseRootConfig maintained_piecewise_config() {
    PiecewiseRootConfig config;
    config.reuse_guards = config.reuse_polish = true;
    config.conjugate_symmetry = config.recover_discovery = config.blocked_polish = true;
    config.durable_ownership = true;
    config.local_band_error = config.aligned_polish = true;
    config.source_finisher = true;
    config.fixed_fft = true;
    config.owned_only_after_discovery = true;
    config.adaptive_taylor = config.interior_proposals = true;
    config.chart_exclusion = true;
    config.lazy_guard_jet = config.reuse_polish_bound = true;
    config.power_composition = true;
    config.predicted_refinement = true;
    config.regime_exit = true;
    config.deflated_discovery = true;
    config.seeds_before_recovery = true;
    return config;
}
} // namespace certpolysolve::detail
