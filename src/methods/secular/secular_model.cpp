// Refreshing the secular interpolation model: compare proposals with the
// retained model, regenerate when needed, and certify before scheduling growth.
// Packet arithmetic lives in secular_packets.cpp; low-level residue
// construction and inclusion proofs keep their existing implementations.
#include "methods/secular/secular_refinement.hpp"
#include "core/telemetry.hpp"
#include <algorithm>
#include <cmath>
#include <limits>

namespace certpolysolve::detail {

SecularRefinement::ModelMatch SecularRefinement::model_match(const cps_options &options) const {
    // The model there already answers this request when every node, every node's iteration bits and the
    // precision budget are the ones it was built for: the same residues, the same cover, the same
    // residues the budget cut short. It is kept and nothing is computed. That is the case after a packet
    // that moved no node, which is what a pair of roots closer than the working precision gives on every
    // rung of the precision's climb, where a rebuild would repeat every finished node's residues and the
    // cover on each rung. History: docs/code-history-20261021.md, secular_refinement.cpp 50.
    bool same_model = nodes_.size() == proposals_.size() && model_bits_ == iteration_bits_ &&
                      model_budget_ == options.max_work_bits;
    for (size_t i = 0; same_model && i < proposals_.size(); ++i)
        same_model = acb_equal(proposals_[i].value(), nodes_[i].value());
    if (!same_model) return ModelMatch::Different;
    return exhausted_before_refresh_ && policy_.model_driven_refresh
        ? ModelMatch::Exhausted : ModelMatch::Unchanged;
}

void SecularRefinement::validate_model(const cps_options &options, SearchStatistics &search,
                                       CandidateStatistics &statistics) {
    TraceScope trace(&policy_, "secular_validation", maximum_precision_, active(), &statistics,
                     &search);
    StageTimer timer(statistics.validation_ns);
    ++statistics.validation_calls;
    std::vector<Root> cover;
    std::vector<uint8_t> inactive;
    if (secular_inclusion_cover(cover, inactive, nodes_, residues_, source_->real(), multiplicity_,
                                static_cast<slong>(maximum_precision_), options, coordinates_, kept_roots_,
                                kept_))
        certified_.consider(std::move(cover), std::move(inactive),
                            static_cast<slong>(maximum_precision_));
    // A node the new cover returns to the active set takes the growths it missed while certified
    // (missed_growths_). Left at its old precision, a lone node within binary64's stage among nodes grown
    // past it would make every later round a mixed one (EnginePolicy::mixed_packets) of single sweeps,
    // where the others' packets would rebuild the model. History:
    // docs/code-history-20261021.md, secular_refinement.cpp 51.
    size_t returned = 0;
    for (size_t i = 0; i < missed_growths_.size(); ++i) {
        returned += !certified_.inactive[i] && missed_growths_[i];
        for (; !certified_.inactive[i] && missed_growths_[i]; --missed_growths_[i])
            iteration_bits_[i] = std::min<slong>(static_cast<slong>(options.max_work_bits),
                                                 resources::precision_product(iteration_bits_[i], 2));
    }
    if (returned && policy_.trace) {
        TraceEvent event{"secular_growth_catch_up", "validation", maximum_precision_, returned, &statistics,
                         &search};
        event.engine = "secular";
        event.reason = "missed-growths";
        policy_.trace(policy_.trace_data, event);
    }
    trace.progress(maximum_precision_, active());
}

SecularRefinement::ModelRefresh SecularRefinement::refresh_model(
    const cps_options &options, SearchStatistics &search, CandidateStatistics &statistics,
    RegenerationWork &regeneration) {
    if (policy_.regeneration.remember_precision && source_precision_hints_.size() != proposals_.size())
        source_precision_hints_.assign(proposals_.size(), 0);
    snapshot("before_regeneration");
    std::vector<ComplexBall> next;
    size_t budget_limited = 0;
    // The verdicts belong to the model they were taken with: they replace
    // the kept ones only when the regeneration delivers its model.
    std::vector<uint8_t> input_limited;
    // The goals are relative to the original origin: a point nearer it than its own magnitude says
    // (RootCoordinates::origin_excess_log2) keeps the goal's bits at that nearer scale. When its
    // coordinates no longer hold them (the excess has used up their guard), they get the bits of a
    // goal that much finer, guard included. A crawl that only passes near the original origin keeps
    // the binary64 stage.
    if (policy_.centred_origin_basis)
        for (size_t i = 0; i < proposals_.size(); ++i)
            if (!certified_.inactive[i]) {
                // Only a node already at the goal's stage: below it the precision growth is still raising
                // the node toward the goal, and the excess is added once it gets there. (Raising every node
                // to the goal at once would send low-precision nodes straight to the goal's packets.)
                // History: docs/code-history-20261021.md, secular_refinement.cpp 49.
                if (iteration_bits_[i] < SecularMachineWorkspace::stage_bits(options.bits))
                    continue;
                const double excess = coordinates_.origin_excess_log2(proposals_[i].value());
                if (!std::isfinite(excess) ||
                    static_cast<double>(options.bits) + excess <= static_cast<double>(iteration_bits_[i]))
                    continue;
                iteration_bits_[i] = std::min<slong>(
                    static_cast<slong>(options.max_work_bits),
                    SecularMachineWorkspace::stage_bits(options.bits + static_cast<uint64_t>(std::ceil(excess))));
            }
    // Separation can put coincident proposals back at their previous
    // interpolation nodes. An exhausted packet followed by this identical
    // model would repeat the same arithmetic, even if exhaustion took more
    // than one evaluation per root. Compare here, after separation, then use
    // the existing validation/growth/reseed path. Consume the observation so
    // a subsequent arithmetic-tier change still gets its own packet.
    const auto match = model_match(options);
    exhausted_before_refresh_ = false;
    if (match == ModelMatch::Different) {
        if (!regenerate_residues(next, *source_, proposals_, nodes_, residues_, iteration_bits_, options,
                                 statistics, maximum_precision_, &policy_, true, source_precision_hints_,
                                 &regeneration, &budget_limited, &input_limited,
                                 policy_.centred_origin_basis ? coordinates_.exact_offset() : nullptr,
                                 restart_estimated_)) {
            if (options.cancelled())
                return ModelRefresh::Cancelled;
            growth_stage_ = GrowthStage::RaisePrecision;
            pending_grow_ = false;
            return ModelRefresh::StepComplete;
        }
        nodes_ = proposals_;
        residues_ = std::move(next);
        input_limited_ = std::move(input_limited);
        model_bits_ = iteration_bits_;
        model_budget_ = options.max_work_bits;
        budget_limited_ = budget_limited;
    }
    // The restarts' estimates are spent: from here the hints are precisions that served.
    std::fill(restart_estimated_.begin(), restart_estimated_.end(), 0);
    snapshot("regenerated_model");
    if (match == ModelMatch::Different)
        validate_model(options, search, statistics);
    need_regeneration_ = false;
    packets_ = 0;
    ++rounds_;
    previous_correction_ = std::numeric_limits<double>::infinity();
    if (!active())
        return ModelRefresh::StepComplete;
    if (policy_.geometric_local_restart && geometric_restart(options, search, statistics)) {
        snapshot("after_geometric_restart");
        separate_interpolation_nodes(proposals_, certified_.inactive, iteration_bits_);
        need_regeneration_ = true;
        return options.cancelled() ? ModelRefresh::Cancelled : ModelRefresh::StepComplete;
    }
    // Certify the newly converged proposals BEFORE changing precision or
    // restarting. Restarting from the old cover here discards valid roots
    // and can make an otherwise convergent secular solve cycle forever.
    if (pending_grow_ || match == ModelMatch::Exhausted) {
        pending_grow_ = false;
        growth_stage_ = GrowthStage::RaisePrecision;
        return ModelRefresh::StepComplete;
    }
    // Regeneration delivers every residue to its demanded accuracy unless
    // the precision budget stops it, so only a model with a residue the
    // budget cut short can leave a step undetermined. Nodes whose step it
    // does not determine take none. When that is every active node, the
    // fresh model is exhausted at once, which is the signal that this
    // precision has nothing more to give: growth below the budget, the
    // precision stop at it. Without this the nodes of a cluster the budget
    // cannot resolve move by noise every round and the solve never stops.
    // History: docs/code-history-20261021.md, secular_refinement.cpp 52.
    const bool any_input_limited =
        std::any_of(input_limited_.begin(), input_limited_.end(), [](uint8_t v) { return v != 0; });
    auto steering = budget_limited_ || any_input_limited
                        ? unsteered_nodes(nodes_, residues_, certified_.inactive)
                        : std::vector<uint8_t>{};
    model_stopping_high_.steering(rounds_, steering);
    model_stopping_.steering(rounds_, std::move(steering));
    return ModelRefresh::PacketReady;
}

} // namespace certpolysolve::detail
