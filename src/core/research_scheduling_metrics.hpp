#pragma once
// Private scheduling observations. The owner outlives every synchronous action
// that writes these fields; a nested local solve owns a separate value. Work is
// accumulated at the executed action's precision, before a search can be moved
// or discarded. These prices are empirical service units, not a bit-work bound.
#include <cstdint>

namespace certpolysolve::detail {
struct ResearchSchedulingMetrics {
    double numerical_work = 0, subdivision_work = 0;
    uint64_t numerical_ns = 0, subdivision_ns = 0;
    uint64_t subdivision_turns = 0, engagements = 0, fallback_goals = 0;

    // Nested local solves are part of the enclosing numerical action. These
    // fields describe that inclusive cost; never add it to numerical_work/ns
    // a second time. accepted counts only a final local solution retained as
    // restart hints, after the parent checks its geometry and conditioning.
    uint64_t local_calls = 0, local_goal_calls = 0, local_accepted_calls = 0;
    uint64_t local_ns = 0, local_subdivision_ns = 0, local_subdivision_turns = 0;
    double local_subdivision_work = 0;
};
} // namespace certpolysolve::detail
