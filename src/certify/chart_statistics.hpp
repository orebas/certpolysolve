#pragma once
// Chart-certificate observations and precision-weighted work. The caller owns these
// values; certification borrows them only for the call and never prices work or schedules a method.
#include "numeric/models/model_error_statistics.hpp"
#include <array>
#include <cstdint>
#include <iterator>

namespace certpolysolve::detail {
struct ChartCertificateWork {
    double bounds = 0;  // local model error bounds
    double jets = 0;    // jet terms, tails included
    double shrinks = 0; // rectangles tried
    bool operator==(const ChartCertificateWork &) const = default;
};
struct ChartCertificateStatistics {
    ModelErrorStatistics model_error;
    uint64_t error_samples = 0, omitted_band_limited = 0, model_tail_limited = 0;
    uint64_t jet_tail_limited = 0, arithmetic_limited = 0, omitted_band_decisive = 0;
    uint64_t jets_reused = 0, radius_shrinks = 0, width_limited = 0;
    // Largest refusal shortfall in model bits: ceil(log2(error / allowance)) + 1.
    // The caller resets this gauge and histogram at each round boundary.
    uint64_t round_model_shortfall_bits = 0;
    // Bucket k holds shortfalls in (2^(k-1), 2^k]; current-round observations are not aggregated.
    std::array<uint64_t, 16> round_shortfall_log2{};
};
namespace statistics {
inline constexpr uint64_t ChartCertificateStatistics::*chart_sum_fields[] = {
    &ChartCertificateStatistics::error_samples,
    &ChartCertificateStatistics::omitted_band_limited,
    &ChartCertificateStatistics::model_tail_limited,
    &ChartCertificateStatistics::jet_tail_limited,
    &ChartCertificateStatistics::arithmetic_limited,
    &ChartCertificateStatistics::omitted_band_decisive,
    &ChartCertificateStatistics::jets_reused,
    &ChartCertificateStatistics::radius_shrinks,
    &ChartCertificateStatistics::width_limited,
};
inline constexpr uint64_t ChartCertificateStatistics::*chart_local_fields[] = {
    &ChartCertificateStatistics::round_model_shortfall_bits,
};
static_assert(sizeof(ChartCertificateStatistics) == sizeof(ModelErrorStatistics) +
    (std::size(chart_sum_fields) + std::size(chart_local_fields)) * sizeof(uint64_t) +
    sizeof(ChartCertificateStatistics::round_shortfall_log2));
} // namespace statistics
} // namespace certpolysolve::detail
