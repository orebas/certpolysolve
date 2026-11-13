// Dispatch one secular candidate packet after model validation, then adopt its
// proposals and schedule the next refresh. This phase owns packet arithmetic and
// exhaustion bookkeeping; the caller retains the meter for the whole step.
#include "methods/secular/secular_refinement.hpp"
#include "methods/secular/secular_internal.hpp"
#include "core/heuristics.hpp"
#include "core/service_meter.hpp"
#include "core/telemetry.hpp"
#include <algorithm>
#include <limits>

namespace certpolysolve::detail {

cps_stop_reason SecularRefinement::advance_packet(const cps_options &options, SearchStatistics &search,
                                                CandidateStatistics &statistics, bool interpolated_model,
                                                PacketWork &work) {
    // The model is regenerated and validated; what follows needs sweeps. With
    // the solve's iteration budget spent, the step returns the budget's stop.
    const unsigned wanted = heuristics::secular_packet_sweeps;
    const unsigned packet = options.packet_sweeps(wanted, statistics);
    if (!packet)
        return CPS_STOP_ITERATION_BUDGET;
    // The packet's sweeps: fewer than wanted when the budget cuts the packet short. A mixed round's sweep is a
    // pair of counted sweeps and may end after its first half.
    unsigned swept = packet;
    TraceScope trace(&policy_, "secular_packet", maximum_precision_, active(), &statistics, &search);
    RootVector initial(source_->degree()), output(source_->degree());
    for (size_t i = 0; i < proposals_.size(); ++i)
        acb_set(initial.data() + i, proposals_[i].value());
    auto *polynomial_source = polynomial_;
    auto *direct = polynomial_source && (polynomial_source->local() || polynomial_source->direct())
                       ? source_.get()
                       : nullptr;
    // The existing sparse coefficient evaluator can be cheaper than a
    // degree-sized secular sum. Retain prepared local models
    // and use regenerated updates for the expensive exact-Taylor path.
    if (prepared_source_policy(&policy_) && polynomial_source &&
        polynomial_source->sparse_coefficients() && !polynomial_source->exact_local_models() &&
        !polynomial_source->local())
        direct = source_.get();
    if (interpolated_model)
        direct = nullptr;
    if (policy_.refinement == RefinementMethod::Aberth)
        direct = source_.get();
    if (policy_.refinement == RefinementMethod::Secular)
        direct = nullptr;
    const bool native = dynamic_cast<SecularSource *>(source_.get()) != nullptr;
    CandidateProgress progress;
    work.terms = statistics.evaluation_terms;
    work.interactions = statistics.root_interactions;
    const bool use_tree =
        source_->degree() >= heuristics::interaction_tree_degree || tree_engine(policy_.candidates) ||
        (polynomial_source && polynomial_source->direct());
    const bool jacobi = !interpolated_model &&
        (!native || tree_engine(policy_.candidates) || persistent_engine(policy_.candidates));
    // A mixed round (EnginePolicy::mixed_packets): the active nodes within the initial stage sweep in binary64
    // while the others are held, then the others sweep in the multiprecision engine while the first group is held
    // where it moved, one sweep each in turn, so that each group always meets the other's latest points. A
    // binary64 refusal falls back to one multiprecision packet for every node, as without the policy.
    // Each half is one counted sweep (candidate_sweeps) under the call's iteration budget: a mixed round needs
    // two sweeps left to start, and a budget that ends inside a pair ends it after its binary64 half (the held
    // nodes keep their points), so a call stopped for the budget has spent it. History:
    // docs/code-history-20261021.md, secular_refinement.cpp 53.
    bool mixed = false;
    std::vector<uint8_t> low_inactive, high_inactive;
    if (policy_.mixed_packets && policy_.study_machine && !direct && !machine_retired_) {
        const slong stage = SecularMachineWorkspace::stage_bits(options.bits);
        low_inactive = high_inactive = certified_.inactive;
        size_t low = 0, high = 0;
        for (size_t i = 0; i < iteration_bits_.size(); ++i) {
            if (certified_.inactive[i]) continue;
            if (iteration_bits_[i] > stage) { low_inactive[i] = 1; ++high; }
            else { high_inactive[i] = 1; ++low; }
        }
        mixed = low && high && options.sweeps_left(statistics) >= 2;
    }
    bool machine = false;
    // The binary64 packet's sub-floor steps at a fresh model's floor (EnginePolicy::machine_node_offsets), except in
    // a higher goal's staged default solve: there the multiprecision packet the floor leads to is not waste, its
    // points start the higher goal's certificate attempt, and floor-stepped points (about 106 bits against that
    // packet's 170) cost the attempt more than the stage saves. History:
    // docs/code-history-20261021.md, secular_refinement.cpp 54.
    const bool floor_steps = policy_.machine_node_offsets && !staged_;
    work.bits = std::numeric_limits<double>::digits; // a binary64 packet
    if (mixed) {
        for (size_t i = 0; i < iteration_bits_.size(); ++i)
            if (!high_inactive[i]) work.bits = std::max<int64_t>(work.bits, iteration_bits_[i]);
        RootVector moved(source_->degree());
        _acb_vec_set(moved.data(), initial.data(), source_->degree());
        CandidateProgress low, high;
        machine = true;
        swept = 0;
        for (unsigned sweep = 0; sweep < packet && machine && options.sweeps_left(statistics); ++sweep) {
            ++swept;
            const auto before = statistics.evaluation_terms + statistics.root_interactions;
            machine = machine_.advance(low, output.data(), moved.data(), nodes_, residues_, low_inactive,
                                       iteration_bits_, rounds_, 1, options, statistics, policy_, &model_stopping_,
                                       policy_.machine_keeps_unmoved, floor_steps);
            work.binary64_units += statistics.evaluation_terms + statistics.root_interactions - before;
            work.service_scaled_binary64 += low.service_scaled_binary64;
            if (!machine || low.cancelled) break;
            if (!options.sweeps_left(statistics)) {
                _acb_vec_set(moved.data(), output.data(), source_->degree());
                break;
            }
            high = secular_candidate_packet(moved.data(), output.data(), nodes_, residues_, high_inactive,
                                            iteration_bits_, 1, options, statistics, nullptr, use_tree, jacobi,
                                            nullptr, rounds_, &policy_, &model_stopping_high_);
            work.service_evaluation += high.service_evaluation;
            work.service_interaction += high.service_interaction;
            if (high.cancelled || !high.finite) break;
            if (low.converged && high.converged) break;
        }
        if (machine) {
            if (policy_.study_statistics) ++policy_.study_statistics->adaptive_packets;
            _acb_vec_set(output.data(), moved.data(), source_->degree());
            progress.cancelled = low.cancelled || high.cancelled;
            progress.finite = low.finite && high.finite;
            progress.correction_log2 = std::max(low.correction_log2, high.correction_log2);
            progress.converged = low.converged && high.converged;
            progress.model_exhausted = low.model_exhausted && high.model_exhausted;
            progress.model_immediate = low.model_immediate && high.model_immediate;
        } else {
            // The binary64 sweeps refused: the whole packet again from the round's entry points, every node at
            // its own precision, as without the policy.
            // Failed low sweeps still consumed binary64 work. Legacy replay
            // retains its historical reset; the service meter pays both attempts.
            if (!policy_.research_service_meter)
                work.binary64_units = 0;
            mixed = false;
        }
    } else if (policy_.study_machine && !direct && !machine_retired_) {
        const auto machine_before = statistics.evaluation_terms + statistics.root_interactions;
        machine = machine_.advance(progress, output.data(), initial.data(), nodes_, residues_, certified_.inactive,
                                   iteration_bits_, rounds_, packet, options, statistics, policy_, &model_stopping_,
                                   policy_.machine_keeps_unmoved && options.bits <= cps_options{}.bits, floor_steps);
        if (policy_.research_service_meter)
            work.binary64_units += statistics.evaluation_terms + statistics.root_interactions - machine_before;
        work.service_scaled_binary64 += progress.service_scaled_binary64;
    }
    last_packet_machine_ = machine;
    if(policy_.study_machine && policy_.study==StudyMethod::Aberth && direct && polynomial_source &&
       !polynomial_source->local()) {
        if(!direct_machine_)direct_machine_=std::make_unique<StudyCandidates>(policy_);
        const slong precision=*std::max_element(iteration_bits_.begin(),iteration_bits_.end());
        const auto coefficients=polynomial_source->coefficient_enclosure(precision);
        progress=direct_machine_->advance(output.data(),coefficients->polynomial.value(),initial.data(),
            precision,packet,options,statistics,certified_.inactive,true);
        for(size_t i=0;i<certified_.inactive.size();++i)if(certified_.inactive[i])acb_set(output.data()+static_cast<slong>(i),initial.data()+static_cast<slong>(i));
        machine=true; // This fixed method owns its arithmetic continuation.
    }
    if (!machine) {
    // A refused binary64 packet may have made counted sweeps already: the multiprecision packet that replaces it
    // makes what the call's iteration budget leaves (the whole packet without a budget).
    swept = options.packet_sweeps(packet, statistics);
    if (!swept)
        return CPS_STOP_ITERATION_BUDGET;
    for (size_t i = 0; i < iteration_bits_.size(); ++i)
        if (!certified_.inactive[i]) work.bits = std::max<int64_t>(work.bits, iteration_bits_[i]);
    if (policy_.study_statistics) ++policy_.study_statistics->adaptive_packets;
    progress = secular_candidate_packet(
        output.data(), initial.data(), nodes_, residues_, certified_.inactive, iteration_bits_, swept, options, statistics,
        direct, use_tree, jacobi, nullptr, rounds_, &policy_, &model_stopping_);
    work.service_evaluation += progress.service_evaluation;
    work.service_interaction += progress.service_interaction;
    }
    if (progress.cancelled)
        return CPS_STOP_CANCELLED;
    const bool adopt_proposals = progress.finite && !(policy_.test_frozen_secular && certified_.available);
    // A repeated model only describes this packet if its proposals were kept.
    // The frozen-iteration control deliberately discards them to exercise the
    // fallback scheduler independently of the iteration's convergence.
    exhausted_before_refresh_ = adopt_proposals && !direct && progress.model_exhausted;
    if (adopt_proposals)
        for (size_t i = 0; i < proposals_.size(); ++i)
            if (!certified_.inactive[i])
                acb_get_mid(proposals_[i].value(), output.data() + i);
    snapshot("after_secular_packet");
    // The solve's iteration budget cut the packet short (the budget is spent),
    // and the packet did not end by itself: the step returns the budget's stop,
    // and the cut packet is not one for the stagnation and refresh rules below,
    // which compare whole packets.
    if (swept < wanted && !options.sweeps_left(statistics) && progress.finite && !progress.converged &&
        !progress.model_exhausted)
        return CPS_STOP_ITERATION_BUDGET;
    ++packets_;
    const bool stagnated = packets_ >= 2 && progress.correction_log2 >= previous_correction_ - 0.5;
    need_regeneration_ = progress.converged || !progress.finite || stagnated || packets_ >= 4;
    pending_grow_ =
        need_regeneration_ && (progress.converged || !progress.finite || (stagnated && rounds_ % 4 == 0));
    if (!direct && policy_.model_stopping == SecularModelPolicy::Continuation && progress.model_exhausted) {
        // Model convergence can request a fresh interpolation at the same
        // coordinate precision. Immediate exhaustion of a fresh model signals
        // the arithmetic floor; original-input certification still runs first.
        pending_grow_ = progress.model_immediate;
    }
    if (!direct && policy_.model_driven_refresh)
        pending_grow_ = need_regeneration_ &&
            (!progress.finite || (progress.model_exhausted && progress.model_immediate));
    if (pending_grow_)
        ++statistics.candidate_stagnations;
    previous_correction_ = progress.correction_log2;
    if (progress.finite)
        last_correction_log2_ = progress.correction_log2;
    return CPS_STOP_NONE;
}
} // namespace certpolysolve::detail
