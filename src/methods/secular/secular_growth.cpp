// Precision growth and resumable candidate restarts for the secular controller.
// The caller's meter spans this phase. A cancelled restart retains its stage so
// a later step resumes the restart without charging another precision increase.
#include "methods/secular/secular_refinement.hpp"
#include "numeric/candidates/candidates.hpp"
#include "core/heuristics.hpp"
#include "core/telemetry.hpp"
#include <algorithm>
#include <limits>
#include <optional>

namespace certpolysolve::detail {

std::optional<cps_stop_reason> SecularRefinement::advance_growth(
    const cps_options &options, SearchStatistics &search, CandidateStatistics &statistics) {
    snapshot(growth_stage_ == GrowthStage::RaisePrecision ? "before_precision_growth" : "before_restart_resume");
    TraceScope growth_trace(&policy_, "secular_precision_growth", precision(), active(), &statistics, &search);
    if (growth_stage_ == GrowthStage::RaisePrecision) {
        // A node at the input's floor: its residue is as accurate as the
        // input allows (input-limited), and the model does not move it:
        // the residue has no correct bit (the node is a root of some
        // polynomial the input admits), or the node's step has none
        // (unsteered: its neighbours' residues are noise). No precision
        // of the iteration changes either. An input-limited residue that
        // still steers its node has not reached the floor: the model was
        // exhausted because the step is below this precision's
        // resolution, and the precision is raised for it as for any other.
        // The steering verdict is taken here, on the model this growth
        // decides for: the one handed to the packets belongs to the model
        // before the regeneration that led here.
        const bool any_input_limited =
            std::any_of(input_limited_.begin(), input_limited_.end(), [](uint8_t v) { return v != 0; });
        const auto unsteered = any_input_limited && residues_.size() == nodes_.size()
                                   ? unsteered_nodes(nodes_, residues_, certified_.inactive)
                                   : std::vector<uint8_t>{};
        const auto at_input_floor = [&](size_t i) {
            return i < input_limited_.size() && input_limited_[i] && i < residues_.size() &&
                   (residue_unresolved(residues_[i].value()) || (i < unsteered.size() && unsteered[i]));
        };
        bool raised = false, every_input_limited = active() > 0;
        for (size_t i = 0; i < iteration_bits_.size(); ++i)
            if (!certified_.inactive[i]) {
                if (at_input_floor(i))
                    continue;
                every_input_limited = false;
                const auto next =
                    std::min<slong>(static_cast<slong>(options.max_work_bits), resources::precision_product(iteration_bits_[i], 2));
                raised = raised || next > iteration_bits_[i];
                iteration_bits_[i] = next;
            }
        if (raised)
            for (size_t i = 0; i < missed_growths_.size(); ++i)
                if (certified_.inactive[i] && missed_growths_[i] < std::numeric_limits<uint8_t>::max())
                    ++missed_growths_[i];
        if (!raised && !every_input_limited && last_packet_machine_ && !machine_retired_) {
            // The binary64 machine is exhausted and the budget allows no
            // more precision than every active node already has: the
            // multiprecision engine has not iterated at that precision.
            // (Not when every node is at the input's floor: the model,
            // not binary64, stopped them, and no engine moves them.)
            // x^8 - 2 (64 x - 1)^2 under an 80-bit budget with a 53-bit
            // goal: the machine's floor leaves the pair near 1/64 at
            // 2^-60 and the goal needs 2^-59; 80-bit sweeps reach it.
            machine_retired_ = true;
            if (policy_.trace) {
                TraceEvent event{"secular_machine_retired", "growth", maximum_precision_, active(), &statistics,
                                 &search};
                event.engine = "secular";
                event.reason = "precision-budget";
                policy_.trace(policy_.trace_data, event);
            }
            growth_stage_ = GrowthStage::None;
            need_regeneration_ = true;
            return CPS_STOP_NONE;
        }
        if (!raised) {
            // At the precision limit the model steers nothing it settled,
            // but a cluster's certificate still depends on its geometry.
            // Try the certificate-optimal placement once; then the stop
            // is honest: no geometry certifies more at this precision.
            bool moved = false;
            if (!limit_reseeded_) {
                limit_reseeded_ = true;
                for (const auto &component : certified_.cover) {
                    const auto members = restart_members(component, nodes_);
                    if (unresolved_members(members))
                        moved = certificate_reseed(members, static_cast<slong>(maximum_precision_)) || moved;
                }
            }
            // Every active node's residue the input's: the input, not
            // the precision budget, bounds this factor's result.
            if (!moved && staged_ && !every_input_limited) {
                // The default goal staged first is not met within the
                // budget's precision: the caller's goal continues from
                // here as it would have from the start.
                staged_ = false;
                growth_stage_ = GrowthStage::None;
                need_regeneration_ = true;
                return CPS_STOP_NONE;
            }
            if (!moved)
                return every_input_limited ? CPS_STOP_INPUT : CPS_STOP_PRECISION;
            if (policy_.trace) {
                TraceEvent event{"secular_certificate_reseed", "limit", maximum_precision_, active(),
                                 &statistics, &search};
                event.engine = "secular";
                policy_.trace(policy_.trace_data, event);
            }
            growth_stage_ = GrowthStage::None;
            need_regeneration_ = true;
            return CPS_STOP_NONE;
        }
        ++search.precision_increases;
        // Restart work can be cancelled. Resume it at this precision,
        // including at the ceiling, without repeating the increase.
        growth_stage_ = GrowthStage::RestartCandidates;
    }
    // A packet can stall because several approximations have collapsed
    // into the same numerical neighbourhood. Restart those approximations
    // on their CERTIFIED component scale, rather than perturbing them by
    // an arbitrarily tiny working-precision ulp and repeatedly doubling.
    std::fill(source_precision_hints_.begin(),source_precision_hints_.end(),0);
    std::fill(restart_estimated_.begin(), restart_estimated_.end(), 0);
    for (const auto &component : certified_.cover) {
        const auto members = restart_members(component, nodes_);
        if (members.empty())
            continue;
        // Residues with no correct bits: the precision limit, where
        // local seeds land at the numerical roots, bunched inside the
        // pseudo-zero set (the worst certificate). Place the members
        // for the certificate instead.
        if (unresolved_members(members) &&
            certificate_reseed(members, static_cast<slong>(maximum_precision_))) {
            if (policy_.trace) {
                TraceEvent event{"secular_certificate_reseed", "growth", maximum_precision_, members.size(),
                                 &statistics, &search};
                event.engine = "secular";
                policy_.trace(policy_.trace_data, event);
            }
            continue;
        }
        slong prec = 64;
        for (const auto i : members)
            prec = std::max(prec, iteration_bits_[i]);
        if (policy_.study_statistics) {
            ++policy_.study_statistics->candidate_restarts;
            policy_.study_statistics->candidates_replaced += members.size();
        }
        TraceScope restart_trace(&policy_, "secular_candidate_restart", static_cast<uint64_t>(prec),
                                 members.size(), &statistics, &search);
        snapshot("before_candidate_restart");
        switch (local_component_restart(component, members, prec, options, statistics, search)) {
        case LocalRestart::Cancelled:
            return CPS_STOP_CANCELLED;
        case LocalRestart::Initialized:
            continue;
        case LocalRestart::None:
            break;
        }
        // A component of more members than a local model takes, whose nodes the model itself shows in
        // groups about their own roots, is restarted by those groups (regroup_restart).
        if (members.size() > static_cast<size_t>(heuristics::local_component_degree)) {
            if (regroup_restart(members, prec, options, statistics, search))
                continue;
            if (options.cancelled())
                return CPS_STOP_CANCELLED;
        }
        RealBall radius;
        Magnitude error;
        mag_add(error.value(), arb_radref(acb_realref(component.box.value())),
                arb_radref(acb_imagref(component.box.value())));
        arf_set_mag(arb_midref(radius.value()), error.value());
        if (arb_is_zero(radius.value()))
            continue;
        arb_mul_2exp_si(radius.value(), radius.value(), -1);
        if (policy_.trace) {
            TraceEvent event{"secular_component_reseed", "circle", static_cast<uint64_t>(prec), members.size(),
                             &statistics, &search};
            event.engine = "secular";
            policy_.trace(policy_.trace_data, event);
        }
        ComplexBall center;
        acb_get_mid(center.value(), component.box.value());
        // On the circle at the placement angles (candidates.hpp), not
        // the unit roots: those are conjugation-invariant about a real
        // centre, and a real source's sweep never breaks that symmetry.
        RootVector placed(static_cast<slong>(members.size()));
        circle_candidates(placed.data(), static_cast<slong>(members.size()), center.value(), radius.value(), prec);
        for (size_t j = 0; j < members.size(); ++j)
            acb_set(proposals_[members[j]].value(), placed.data() + static_cast<slong>(j));
    }
    growth_stage_ = GrowthStage::None;
    need_regeneration_ = true;
    return std::nullopt;
}

} // namespace certpolysolve::detail
