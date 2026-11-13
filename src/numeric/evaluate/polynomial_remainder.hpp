#pragma once
#include "core/raii.hpp"
#include <algorithm>
#include <vector>

namespace certpolysolve::detail {
// A coefficient majorant for the linear Taylor remainder. Construction borrows
// the coefficient view only for the call; the value retains just owned
// magnitudes of the second-derivative coefficients, never the polynomial.
class PolynomialLinearRemainder {
    std::vector<Magnitude> second_;

  public:
    explicit PolynomialLinearRemainder(ConstComplexBallView values) {
        second_.resize(static_cast<size_t>(std::max<slong>(0, values.size() - 2)));
        for (slong j = 2; j < values.size(); ++j) {
            auto &m = second_[static_cast<size_t>(j - 2)];
            acb_get_mag(m.value(), values[j]);
            mag_mul_ui(m.value(), m.value(), static_cast<ulong>(j));
            mag_mul_ui(m.value(), m.value(), static_cast<ulong>(j - 1));
        }
    }
    // For |z-c| <= r, Taylor's integral remainder is at most
    // (r²/2) sum_{j>=2} j(j-1)|a_j| (|c|+r)^(j-2).
    void linear_remainder(arb_t out, const acb_t center, const arb_t radius) const {
        Magnitude norm, r, bound;
        acb_get_mag(norm.value(), center);
        arb_get_mag(r.value(), radius);
        mag_add(norm.value(), norm.value(), r.value());
        for (size_t j = second_.size(); j-- > 0;) {
            mag_mul(bound.value(), bound.value(), norm.value());
            mag_add(bound.value(), bound.value(), second_[j].value());
        }
        mag_mul(bound.value(), bound.value(), r.value());
        mag_mul(bound.value(), bound.value(), r.value());
        mag_mul_2exp_si(bound.value(), bound.value(), -1);
        arb_zero(out);
        arf_set_mag(arb_midref(out), bound.value());
    }
};
} // namespace certpolysolve::detail
