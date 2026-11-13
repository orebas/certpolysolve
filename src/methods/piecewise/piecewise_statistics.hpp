#pragma once
// Observations and counted work of a piecewise solve, independent of its recipe
// and controller. Factor aggregation is defined in piecewise_statistics.cpp.
#include "certify/chart_statistics.hpp"
#include "methods/piecewise/piecewise_work.hpp"
#include <array>
#include <cstdint>
#include <iterator>

namespace certpolysolve::detail {
struct PiecewiseRootStatistics {
    PiecewiseWork work;
    // The part of `work` done while a finishing cohort was being formed or
    // finished (its certificates charge the certificate kinds), so the
    // finisher's cost can be compared with its own timer.
    PiecewiseWork finishing_work;
    ChartCertificateStatistics certificate;
    uint64_t setup_ns = 0, construction_ns = 0, local_solve_ns = 0;
    uint64_t certification_ns = 0, deduplication_ns = 0;
    uint64_t rings = 0, models = 0, coefficients = 0, maximum_degree = 0;
    uint64_t proposals = 0, certificates = 0, duplicates = 0, rounds = 0;
    uint64_t overlap_attempts = 0, overlap_identified = 0, overlap_resolved = 0;
    uint64_t unproved_seeds_retained = 0;
    uint64_t model_bits = 0, arithmetic_bits = 0, fft_transforms = 0;
    uint64_t guards = 0, guard_attempts = 0, retired_points = 0;
    uint64_t model_limited = 0, rounding_limited = 0, candidate_limited = 0;
    uint64_t seed_guards = 0, seed_retired = 0, seed_rejected = 0;
    uint64_t skipped_models = 0, storage_refusals = 0, axis_gap_failures = 0, axis_undecided = 0, axis_unrefined = 0, axis_sign_conflicts = 0;
    // Coefficient exponents beyond the envelope's working range (|height| >
    // 2^50, so index x height products stay exact): a refusal of the input.
    uint64_t range_refusals = 0;
    uint64_t accepted_moves = 0, reseeded_models = 0, local_counts = 0;
    uint64_t batches = 0;
    uint64_t machine_ns = 0, polish_ns = 0, seed_certificate_ns = 0, final_certificate_ns = 0;
    uint64_t one_root_models = 0, other_counted_models = 0;
    uint64_t machine_sweeps = 0, machine_points = 0, reused_guards = 0;
    // Points the charts' binary64 iterations started (machine_points counts the ones they kept).
    uint64_t machine_seeds = 0;
    uint64_t excluded_charts = 0;
    uint64_t local_count_ns = 0;
    uint64_t unsettled_points = 0;
    uint64_t polish_points = 0, polish_evaluations = 0, polish_model_stops = 0;
    uint64_t polish_rounding_stops = 0, polish_step_stops = 0, polish_reused = 0;
    uint64_t polish_local_bounds = 0, model_stopped_candidates = 0;
    uint64_t polish_bounds_skipped = 0;
    uint64_t polish_evaluation_ns = 0, polish_error_bound_ns = 0, diagnostic_ns = 0;
    uint64_t construction_units = 0;
    uint64_t finisher_ns = 0, finisher_admission_ns = 0, finisher_evaluations = 0, finisher_terms = 0;
    uint64_t finisher_preparation_ns = 0;
    uint64_t finisher_block_evaluations = 0;
    uint64_t finisher_attempts = 0, finisher_certificates = 0, finisher_refusals = 0;
    uint64_t reflected_roots = 0, symmetry_skipped_models = 0;
    uint64_t recovery_rounds = 0, repeated_recoveries_avoided = 0;
    uint64_t stationary_rounds = 0;
    uint64_t predicted_rungs_skipped = 0;
    uint64_t regime_exits = 0;
    // Deflated discovery (PiecewiseRootConfig::deflated_discovery): stalled
    // rounds that ran it, the free points' sweeps, the seeds it added.
    uint64_t deflation_rounds = 0, deflation_sweeps = 0, deflation_seeds = 0, deflation_ns = 0;
    // The certificate attempt before piecewise (PiecewiseRootConfig): made,
    // and finished the solve.
    uint64_t certificate_attempts = 0, certificate_attempt_finishes = 0;
    uint64_t candidate_improvements = 0;
    uint64_t preprocessing_ns = 0, prepared_factors = 0, preprocessing_limited = 0;
    uint64_t coefficient_preparation_ns = 0, taylor_rows_ns = 0, fft_setup_ns = 0;
    uint64_t fft_execution_ns = 0, model_assembly_ns = 0;
    // Transforms the fixed-point transform delivered and declined (a declined one is a ball transform), and its
    // time: conversion, butterflies and enclosure.
    uint64_t fft_fixed_accepted = 0, fft_fixed_declined = 0, fft_fixed_ns = 0;
    uint64_t ownership_ns = 0, ownership_reuses = 0, ownership_conflicts = 0;
    uint64_t ownership_escapes = 0, owned_proposals_retired = 0, counted_discovery_retired = 0;
    // Fresh proposals dropped before polishing as copies of a held root (counted in retired_points or
    // owned_proposals_retired as well).
    uint64_t unpolished_copies = 0;
    uint64_t discovery_count = 0, final_accurate_count = 0;
    uint64_t core_proposal_models = 0, core_proposal_setup_ns = 0;
    PiecewiseRootStatistics &operator+=(const PiecewiseRootStatistics &);
};

// Factor aggregation adds ordinary counters and takes maxima. These descriptions
// live beside their type so a new scalar field requires an explicit classification.
namespace statistics {
inline constexpr uint64_t PiecewiseRootStatistics::*piecewise_sum_fields[] = {
    &PiecewiseRootStatistics::core_proposal_models,
    &PiecewiseRootStatistics::core_proposal_setup_ns,
    &PiecewiseRootStatistics::polish_evaluation_ns,
    &PiecewiseRootStatistics::polish_error_bound_ns,
    &PiecewiseRootStatistics::diagnostic_ns,
    &PiecewiseRootStatistics::polish_bounds_skipped,
    &PiecewiseRootStatistics::fft_fixed_accepted,
    &PiecewiseRootStatistics::fft_fixed_declined,
    &PiecewiseRootStatistics::fft_fixed_ns,
    &PiecewiseRootStatistics::setup_ns,
    &PiecewiseRootStatistics::construction_ns,
    &PiecewiseRootStatistics::local_solve_ns,
    &PiecewiseRootStatistics::certification_ns,
    &PiecewiseRootStatistics::deduplication_ns,
    &PiecewiseRootStatistics::rings,
    &PiecewiseRootStatistics::models,
    &PiecewiseRootStatistics::coefficients,
    &PiecewiseRootStatistics::proposals,
    &PiecewiseRootStatistics::certificates,
    &PiecewiseRootStatistics::duplicates,
    &PiecewiseRootStatistics::rounds,
    &PiecewiseRootStatistics::fft_transforms,
    &PiecewiseRootStatistics::guards,
    &PiecewiseRootStatistics::guard_attempts,
    &PiecewiseRootStatistics::retired_points,
    &PiecewiseRootStatistics::model_limited,
    &PiecewiseRootStatistics::rounding_limited,
    &PiecewiseRootStatistics::candidate_limited,
    &PiecewiseRootStatistics::seed_guards,
    &PiecewiseRootStatistics::seed_retired,
    &PiecewiseRootStatistics::seed_rejected,
    &PiecewiseRootStatistics::skipped_models,
    &PiecewiseRootStatistics::storage_refusals,
    &PiecewiseRootStatistics::range_refusals,
    &PiecewiseRootStatistics::accepted_moves,
    &PiecewiseRootStatistics::reseeded_models,
    &PiecewiseRootStatistics::local_counts,
    &PiecewiseRootStatistics::batches,
    &PiecewiseRootStatistics::machine_ns,
    &PiecewiseRootStatistics::polish_ns,
    &PiecewiseRootStatistics::seed_certificate_ns,
    &PiecewiseRootStatistics::final_certificate_ns,
    &PiecewiseRootStatistics::one_root_models,
    &PiecewiseRootStatistics::other_counted_models,
    &PiecewiseRootStatistics::machine_sweeps,
    &PiecewiseRootStatistics::machine_seeds,
    &PiecewiseRootStatistics::machine_points,
    &PiecewiseRootStatistics::excluded_charts,
    &PiecewiseRootStatistics::reused_guards,
    &PiecewiseRootStatistics::local_count_ns,
    &PiecewiseRootStatistics::unsettled_points,
    &PiecewiseRootStatistics::polish_points,
    &PiecewiseRootStatistics::polish_evaluations,
    &PiecewiseRootStatistics::polish_model_stops,
    &PiecewiseRootStatistics::polish_rounding_stops,
    &PiecewiseRootStatistics::polish_step_stops,
    &PiecewiseRootStatistics::polish_reused,
    &PiecewiseRootStatistics::polish_local_bounds,
    &PiecewiseRootStatistics::model_stopped_candidates,
    &PiecewiseRootStatistics::construction_units,
    &PiecewiseRootStatistics::finisher_ns,
    &PiecewiseRootStatistics::finisher_admission_ns,
    &PiecewiseRootStatistics::finisher_preparation_ns,
    &PiecewiseRootStatistics::finisher_evaluations,
    &PiecewiseRootStatistics::finisher_block_evaluations,
    &PiecewiseRootStatistics::finisher_terms,
    &PiecewiseRootStatistics::finisher_attempts,
    &PiecewiseRootStatistics::finisher_certificates,
    &PiecewiseRootStatistics::finisher_refusals,
    &PiecewiseRootStatistics::reflected_roots,
    &PiecewiseRootStatistics::symmetry_skipped_models,
    &PiecewiseRootStatistics::recovery_rounds,
    &PiecewiseRootStatistics::repeated_recoveries_avoided,
    &PiecewiseRootStatistics::stationary_rounds,
    &PiecewiseRootStatistics::overlap_attempts,
    &PiecewiseRootStatistics::overlap_identified,
    &PiecewiseRootStatistics::overlap_resolved,
    &PiecewiseRootStatistics::unproved_seeds_retained,
    &PiecewiseRootStatistics::candidate_improvements,
    &PiecewiseRootStatistics::preprocessing_ns,
    &PiecewiseRootStatistics::prepared_factors,
    &PiecewiseRootStatistics::preprocessing_limited,
    &PiecewiseRootStatistics::coefficient_preparation_ns,
    &PiecewiseRootStatistics::taylor_rows_ns,
    &PiecewiseRootStatistics::fft_setup_ns,
    &PiecewiseRootStatistics::fft_execution_ns,
    &PiecewiseRootStatistics::model_assembly_ns,
    &PiecewiseRootStatistics::ownership_ns,
    &PiecewiseRootStatistics::ownership_reuses,
    &PiecewiseRootStatistics::ownership_conflicts,
    &PiecewiseRootStatistics::ownership_escapes,
    &PiecewiseRootStatistics::owned_proposals_retired,
    &PiecewiseRootStatistics::unpolished_copies,
    &PiecewiseRootStatistics::counted_discovery_retired,
    &PiecewiseRootStatistics::discovery_count,
    &PiecewiseRootStatistics::final_accurate_count,
    &PiecewiseRootStatistics::regime_exits,
    &PiecewiseRootStatistics::deflation_rounds,
    &PiecewiseRootStatistics::deflation_sweeps,
    &PiecewiseRootStatistics::deflation_seeds,
    &PiecewiseRootStatistics::deflation_ns,
    &PiecewiseRootStatistics::certificate_attempts,
    &PiecewiseRootStatistics::certificate_attempt_finishes,
};
inline constexpr uint64_t PiecewiseRootStatistics::*piecewise_max_fields[] = {
    &PiecewiseRootStatistics::maximum_degree,
    &PiecewiseRootStatistics::model_bits,
    &PiecewiseRootStatistics::arithmetic_bits,
};
// These diagnostics have not participated in factor aggregation. Extending their
// reporting is a separate behavior change; the current round's histogram is local too.
inline constexpr uint64_t PiecewiseRootStatistics::*piecewise_local_fields[] = {
    &PiecewiseRootStatistics::axis_gap_failures,
    &PiecewiseRootStatistics::axis_undecided,
    &PiecewiseRootStatistics::axis_unrefined,
    &PiecewiseRootStatistics::axis_sign_conflicts,
    &PiecewiseRootStatistics::predicted_rungs_skipped,
};
static_assert(sizeof(PiecewiseRootStatistics) == 2 * sizeof(PiecewiseWork) + sizeof(ChartCertificateStatistics) +
    (std::size(piecewise_sum_fields) + std::size(piecewise_max_fields) + std::size(piecewise_local_fields)) * sizeof(uint64_t));
} // namespace statistics

} // namespace certpolysolve::detail
