#pragma once
#include "core/number.hpp"
#include <vector>

namespace certpolysolve::detail {
struct SecularTerm {
    GaussianRational residue, pole;
};
} // namespace certpolysolve::detail
struct cps_secular {
    std::vector<certpolysolve::detail::SecularTerm> terms;
};
