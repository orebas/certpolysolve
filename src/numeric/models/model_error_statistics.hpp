#pragma once
// Diagnostics of local model-error bounds, independent of the method requesting them.
#include <cstdint>
#include <iterator>

namespace certpolysolve::detail {
struct ModelErrorStatistics {
    uint64_t local_band_attempts = 0, local_band_improvements = 0;
    uint64_t local_taylor_attempts = 0, local_taylor_improvements = 0;
};
namespace statistics {
inline constexpr uint64_t ModelErrorStatistics::*model_error_sum_fields[] = {
    &ModelErrorStatistics::local_band_attempts,
    &ModelErrorStatistics::local_band_improvements,
    &ModelErrorStatistics::local_taylor_attempts,
    &ModelErrorStatistics::local_taylor_improvements,
};
static_assert(sizeof(ModelErrorStatistics) == std::size(model_error_sum_fields) * sizeof(uint64_t));
} // namespace statistics
} // namespace certpolysolve::detail
