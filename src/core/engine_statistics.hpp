#pragma once
// The engines' counters as they accumulate during a solve. The C API publishes them through cps_statistics and
// cps_engine_counters (c_api.cpp). These are the layouts of its former public blocks (cps_statistics_v1, _v2, _v3
// and cps_progress_v1), value-initialised with {} as those were.
// History: docs/code-history-20261021.md, engine_statistics.hpp 1.
#include <certpolysolve/certpolysolve.h>
#include <cstdint>
#include <iterator>

namespace certpolysolve::detail {
// Stage times and counts of the engine that answers. Durations are measured wall-clock nanoseconds in disjoint
// stages, not CPU time. Isolation includes backend certification, local recovery and endpoint export. A hand-off
// carries an earlier engine's apart (CarriedStatistics, session_statistics.hpp).
struct StageStatistics {
    uint64_t preprocessing_ns;
    uint64_t input_enclosure_ns;
    uint64_t isolation_ns;
    uint64_t coverage_validation_ns;
    uint64_t precision_attempts;
    uint64_t isolation_passes;
    uint64_t exact_factors;
};
// The native search's cumulative stages and counts; backend_ns stays opaque, and counting excludes taylor_shift_ns.
struct SearchStatistics {
    uint64_t backend_ns, taylor_shift_ns, counting_ns, refinement_ns, search_ns;
    uint64_t disk_tests, graeffe_steps, taylor_shifts, taylor_cache_hits;
    uint64_t boxes_split, boxes_excluded, newton_attempts, newton_successes;
    uint64_t work_units, backend_calls, precision_increases;
};
// The candidate engines' cumulative observations (not certificates). Candidate time includes evaluation,
// interactions and regeneration; validation is separate. double_seed_calls counts the research build's DoubleSeed
// starts: zero in the default build, kept so that no public layout changes.
struct CandidateStatistics {
    uint64_t candidate_ns, evaluation_ns, interaction_ns, validation_ns, regeneration_ns;
    uint64_t candidate_sweeps, evaluations, evaluation_terms, root_interactions;
    uint64_t validation_calls, secular_regenerations, double_seed_calls, candidate_stagnations;
    uint64_t real_refinement_steps, real_bisections;
};
// Descriptions shared by nested/attempt aggregation, engine handoff and publication.
// exact_factors describes the answering preparation, so it is never carried.
namespace statistics {
inline constexpr uint64_t StageStatistics::*carried_stage_fields[] = {
    &StageStatistics::preprocessing_ns, &StageStatistics::input_enclosure_ns, &StageStatistics::isolation_ns,
    &StageStatistics::coverage_validation_ns, &StageStatistics::precision_attempts,
    &StageStatistics::isolation_passes};
inline constexpr uint64_t SearchStatistics::*search_fields[] = {
    &SearchStatistics::backend_ns, &SearchStatistics::taylor_shift_ns, &SearchStatistics::counting_ns,
    &SearchStatistics::refinement_ns, &SearchStatistics::search_ns, &SearchStatistics::disk_tests,
    &SearchStatistics::graeffe_steps, &SearchStatistics::taylor_shifts, &SearchStatistics::taylor_cache_hits,
    &SearchStatistics::boxes_split, &SearchStatistics::boxes_excluded, &SearchStatistics::newton_attempts,
    &SearchStatistics::newton_successes, &SearchStatistics::work_units, &SearchStatistics::backend_calls,
    &SearchStatistics::precision_increases};
inline constexpr uint64_t CandidateStatistics::*candidate_fields[] = {
    &CandidateStatistics::candidate_ns, &CandidateStatistics::evaluation_ns, &CandidateStatistics::interaction_ns,
    &CandidateStatistics::validation_ns, &CandidateStatistics::regeneration_ns, &CandidateStatistics::candidate_sweeps,
    &CandidateStatistics::evaluations, &CandidateStatistics::evaluation_terms, &CandidateStatistics::root_interactions,
    &CandidateStatistics::validation_calls, &CandidateStatistics::secular_regenerations,
    &CandidateStatistics::double_seed_calls, &CandidateStatistics::candidate_stagnations,
    &CandidateStatistics::real_refinement_steps, &CandidateStatistics::real_bisections};
// A field added to a struct must be added to its list (v1: every field but exact_factors).
static_assert(sizeof(StageStatistics) == (std::size(carried_stage_fields) + 1) * sizeof(uint64_t));
static_assert(sizeof(SearchStatistics) == std::size(search_fields) * sizeof(uint64_t));
static_assert(sizeof(CandidateStatistics) == std::size(candidate_fields) * sizeof(uint64_t));

} // namespace statistics

// Where the search stood when the call ended.
struct SearchProgress {
    uint64_t work_units_this_call, pending_components, live_cells;
    uint64_t min_pending_bits, max_pending_bits, max_outer_bits;
    uint64_t preprocessing_limited, input_uncertainty;
    cps_stop_reason stop_reason;
};
} // namespace certpolysolve::detail
