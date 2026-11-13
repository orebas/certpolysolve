#pragma once
// What the parts of the secular refinement share and nothing else reads: secular_refinement.cpp,
// secular_regeneration.cpp, secular_validation.cpp, secular_restarts.cpp (split in W8's S5,
// docs/w8-plan-20261021.md). Each function is defined in the part named beside it.
#include "core/engine_policy.hpp"
#include "core/number.hpp"
#include "numeric/evaluate/structured_source.hpp"
#include <cstdint>
#include <span>

namespace certpolysolve::detail {
// The bits an evaluation carries beyond the accuracy it is asked for, for what it cancels: a node's first
// evaluation in a regeneration beyond its residue's bits (its iteration bits, or after a restart the precision
// its component's model showed), and a restart's centre steps beyond their aim. A measured shortfall raises it.
inline constexpr slong regeneration_guard_bits = 32;
bool exact_nodes(std::span<const ComplexBall> nodes);                             // secular_validation.cpp
bool prepared_source_policy(const EnginePolicy *policy);                          // secular_refinement.cpp
bool conditioned_secular(const MonicSource &source, const EnginePolicy *policy);  // secular_refinement.cpp
slong adequate_accuracy_bits(slong bits, uint64_t goal, bool iteration_model);    // secular_regeneration.cpp
} // namespace certpolysolve::detail
