#pragma once
#include "core/cover.hpp"
#include <span>

namespace certpolysolve::detail {
// Each input rectangle contains its stated root set. On a rectangle avoiding
// zero and the chosen branch cut, the k root branches are injective and their
// images exhaust its preimages. Output rectangles can overlap; the ordinary
// complete-cover union ledger must merge those before publishing a result.
bool lift_power_cover(std::vector<Root> &out, std::span<const Root> reduced, ulong power, slong precision);
} // namespace certpolysolve::detail
