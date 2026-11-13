#pragma once
// Immutable concave coefficient-height envelope. Geometry may estimate with doubles;
// bounds use directed ball arithmetic and orientation uses exact integers.
#include "core/number.hpp"
#include "core/wide_int.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace certpolysolve::detail {
// A concave, piecewise linear upper bound for log2 |a_j|. Integer heights
// round OUTWARD; hull orientation uses exact integer arithmetic. Geometry may
// use double estimates, while every tail bound uses directed ball arithmetic.
class CoefficientEnvelope {
    struct Vertex { slong index, height; };
    std::vector<Vertex> hull_;
    slong constant_lower_ = 0, leading_lower_ = 0;

    // The lower bounds of the constant and the leading coefficient.
    void ends(const acb_t constant, const acb_t leading) {
        Float bound;
        acb_get_abs_lbound_arf(bound.value(), constant, 64);
        constant_lower_ = arf_abs_bound_lt_2exp_si(bound.value()) - 1;
        acb_get_abs_lbound_arf(bound.value(), leading, 64);
        leading_lower_ = arf_abs_bound_lt_2exp_si(bound.value()) - 1;
    }
    // Coefficient j, in increasing order: one that is zero (its upper bound 0) adds nothing.
    void add(slong j, const acb_t c) {
        Float bound;
        acb_get_abs_ubound_arf(bound.value(), c, 64);
        if (arf_is_zero(bound.value()))
            return;
        if (!arf_is_finite(bound.value()))
            throw std::invalid_argument("nonfinite coefficient envelope");
        const slong height = arf_abs_bound_lt_2exp_si(bound.value());
        // Keep index/height products inside the exact orientation adapter.
        if (height < -(slong(1) << 50) || height > (slong(1) << 50))
            throw std::length_error("piecewise exponent envelope exceeds private working range");
        while (hull_.size() > 1) {
            const auto a = hull_[hull_.size() - 2], b = hull_.back();
            const int128 left = static_cast<int128>(b.height - a.height) * (j - b.index);
            const int128 right = static_cast<int128>(height - b.height) * (b.index - a.index);
            if (left > right)
                break;
            hull_.pop_back();
        }
        hull_.push_back({j, height});
    }

  public:
    explicit CoefficientEnvelope(const acb_poly_t p) {
        ends(p->coeffs, p->coeffs + acb_poly_degree(p));
        for (slong j = 0; j < acb_poly_length(p); ++j)
            add(j, p->coeffs + j);
    }
    // The same from the nonzero terms (exponent, ball) in increasing order, the first the constant and the last the
    // leading coefficient: a polynomial's zero coefficients add nothing to the hull, so it is the envelope of the
    // polynomial they make (C40's second step: the selector on the input object's terms).
    explicit CoefficientEnvelope(std::span<const std::pair<slong, ComplexBall>> terms) {
        ends(terms.front().second.value(), terms.back().second.value());
        for (const auto &[j, c] : terms)
            add(j, c.value());
    }
    slong degree() const { return hull_.back().index; }
    std::vector<double> negative_heights(slong n) const {
        std::vector<double> result(static_cast<size_t>(n + 1));
        for (size_t k = 1; k < hull_.size(); ++k) {
            const auto a = hull_[k - 1], b = hull_[k];
            for (slong j = a.index; j <= b.index; ++j)
                result[static_cast<size_t>(j)] =
                    -static_cast<double>(a.height) - static_cast<double>(b.height - a.height) *
                                                         static_cast<double>(j - a.index) /
                                                         static_cast<double>(b.index - a.index);
        }
        return result;
    }
    std::vector<slong> indices() const {
        std::vector<slong> result;
        for (auto v : hull_)
            result.push_back(v.index);
        return result;
    }
    // The hull's largest excess over the chord of its ends, exact: the excess times the ends' index span, and the
    // span. The bend, in bits, is their quotient.
    struct Bend {
        int128 excess = 0;
        slong span = 1;
    };
    Bend bend() const {
        if (hull_.size() < 2)
            return {};
        const auto a = hull_.front(), b = hull_.back();
        int128 maximum = 0;
        for (const auto v : hull_) {
            const int128 excess = static_cast<int128>(v.height - a.height) * (b.index - a.index) -
                                  static_cast<int128>(b.height - a.height) * (v.index - a.index);
            maximum = std::max(maximum, excess);
        }
        return {maximum, b.index - a.index};
    }
    std::pair<double, double> root_window(const std::vector<double> &negative_heights) const {
        const slong n = hull_.back().index;
        double inner = std::numeric_limits<double>::infinity(), outer = -inner;
        // Include the truncated-domain endpoints j=1 and j=n-1 as well as
        // hull vertices. A single envelope scan remains linear in degree.
        for (slong j = 0; j <= n; ++j) {
            const double height = -negative_heights[static_cast<size_t>(j)];
            if (j)
                inner =
                    std::min(inner, (static_cast<double>(constant_lower_) - height) / static_cast<double>(j));
            if (j < n)
                outer = std::max(outer,
                                 (height - static_cast<double>(leading_lower_)) / static_cast<double>(n - j));
        }
        return {std::floor(inner) - 2, std::ceil(outer) + 2};
    }
    // Bound sum_{j=lo}^hi |a_j| r^(j-shift), without scanning coefficients.
    // A concave envelope has a unimodal support function. Exact comparisons
    // of dyadic log bounds locate its maximum in O(log degree) operations.
    void bound(arb_t out, slong lo, slong hi, const arb_t r, slong shift, slong scale) const {
        arb_zero(out);
        if (lo > hi)
            return;
        if (!arb_is_positive(r)) {
            arb_indeterminate(out);
            return;
        }
        RealBall log, log_two, height, exponent;
        Float x, difference, maximum;
        arb_log(log.value(), r, 96);
        arb_const_log2(log_two.value(), 96);
        arb_div(log.value(), log.value(), log_two.value(), 96);
        if (hi <= shift)
            arb_get_lbound_arf(x.value(), log.value(), 96);
        else if (lo >= shift)
            arb_get_ubound_arf(x.value(), log.value(), 96);
        else {
            RealBall left, right;
            bound(left.value(), lo, shift, r, shift, scale);
            bound(right.value(), shift + 1, hi, r, shift, scale);
            arb_add(out, left.value(), right.value(), 96);
            return;
        }
        size_t a = 0, b = hull_.size() - 1;
        while (a < b) {
            const size_t mid = a + (b - a) / 2;
            const auto u = hull_[mid], v = hull_[mid + 1];
            arf_mul_si(difference.value(), x.value(), v.index - u.index, ARF_PREC_EXACT, ARF_RND_NEAR);
            arf_add_si(difference.value(), difference.value(), v.height - u.height, ARF_PREC_EXACT, ARF_RND_NEAR);
            if (arf_sgn(difference.value()) > 0)
                a = mid + 1;
            else
                b = mid;
        }
        const slong j = std::clamp(hull_[a].index, lo, hi);
        auto it =
            std::lower_bound(hull_.begin(), hull_.end(), j, [](Vertex v, slong i) { return v.index < i; });
        if (it->index == j)
            arb_set_si(height.value(), it->height);
        else {
            const auto u = *(it - 1), v = *it;
            arb_set_si(height.value(), v.height - u.height);
            arb_mul_si(height.value(), height.value(), j - u.index, 96);
            arb_div_si(height.value(), height.value(), v.index - u.index, 96);
            arb_add_si(height.value(), height.value(), u.height, 96);
        }
        arb_set_arf(exponent.value(), x.value());
        arb_mul_si(exponent.value(), exponent.value(), j - shift, 96);
        arb_add(exponent.value(), exponent.value(), height.value(), 96);
        arb_sub_si(exponent.value(), exponent.value(), scale, 96);
        arb_get_ubound_arf(maximum.value(), exponent.value(), 96);
        Integer e;
        arf_get_fmpz(e.value(), maximum.value(), ARF_RND_CEIL);
        arb_one(out);
        arb_mul_2exp_fmpz(out, out, e.value());
        arb_mul_si(out, out, hi - lo + 1, 96);
    }
};

} // namespace certpolysolve::detail
