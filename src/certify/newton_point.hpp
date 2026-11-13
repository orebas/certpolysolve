#pragma once
// Certify one proposed point against its original monic source. This predicate
// does not run candidate packets or retry at another precision.
#include "core/cover.hpp"
#include "core/query.hpp"

namespace certpolysolve::detail {
class MonicSource;
// A linear local Taylor model with a derivative-variation remainder. The same
// published-rectangle Rouche predicate certifies the original source.
bool certify_newton_point(Root &, MonicSource &, const acb_t, const cps_options &, slong precision);
} // namespace certpolysolve::detail
