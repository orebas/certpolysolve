#pragma once
#include <cstdint>

namespace certpolysolve::detail {
// Optional observations of the existing budget checks, not a new allowance.
// A rejected operation reports a lower bound, not an invented exact work cost.
struct PreparationDiagnostic {
    const char *reason = "not-attempted";
    uint64_t required_bits_lower_bound = 0;
    uint64_t input_bits = 0, output_bits = 0;
    bool accepted = false;
    bool reject(const char *why, uint64_t budget) {
        reason = why;
        required_bits_lower_bound = budget == UINT64_MAX ? budget : budget + 1;
        return false;
    }
};
} // namespace certpolysolve::detail
