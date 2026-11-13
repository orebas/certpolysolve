#pragma once

#include "platform/parallel.hpp"
#include "numeric/evaluate/evaluator.hpp"
#include <concepts>
#include <algorithm>
#include <cmath>

namespace certpolysolve::detail {
// FLINT block dot products, with disk errors between blocks.
// This object owns both polynomials and its reusable scratch. Calls on the
// same object must be exclusive; there is no mutable shared/global workspace.
// The fixed-precision snapshot encloses every coefficient/point combination.
// Both storage choices own their source. The snapshot variant shares immutable
// coefficients with the remainder evaluator; the ordinary variant keeps its
// existing unique polynomial storage without adding a shared allocation.
template<class Polynomial>
    requires (std::same_as<Polynomial, ComplexPolynomial> || std::same_as<Polynomial, PolynomialSnapshot>)
class OwnedDiskBlockedPolynomialEvaluator {
    Polynomial polynomial_;
    ComplexPolynomial derivative_;
    slong precision_, block_;
    RootVector powers_;
    // The blocks' dot products, one slot a block, where W9's split computes them before the fold (below).
    RootVector terms_;
    ComplexBall term_, factor_midpoint_;
    Magnitude factor_norm_, factor_error_;
    bool finite_ = true;

    static slong block_size(const acb_poly_struct *polynomial, slong precision) {
        if (precision < 2) throw std::invalid_argument("evaluation precision must be at least two bits");
        return std::max<slong>(1, static_cast<slong>(std::sqrt(
            static_cast<double>(acb_poly_length(polynomial)))) + 1);
    }

    static slong block_count(const acb_poly_struct *polynomial, slong block) {
        const slong length = acb_poly_length(polynomial);
        return length ? (length - 1) / block + 1 : 1;
    }
    // W9's split (parallel.hpp): from split_min_length coefficients, on a solve with workers, the dot products run
    // first, in fixed groups of split_group blocks, each into its block's slot; the fold below then reads them in
    // today's order. Each dot product is the same FLINT call on the same terms, so the bits are today's. Both
    // constants fix the groups by the length alone; at 4,096 coefficients a block holds 65 terms.
    static constexpr slong split_min_length = 4096, split_group = 8;
    void evaluate_polynomial(acb_t out, const acb_poly_struct *polynomial) {
        const slong length = acb_poly_length(polynomial);
        acb_zero(out);
        if (!length) return;
        const slong blocks = (length - 1) / block_ + 1;
        const bool split = length >= split_min_length && parallel_width() > 1 && blocks <= terms_.size();
        if (split)
            parallel_for(SplitSite::BlockedEvaluation, static_cast<size_t>((blocks + split_group - 1) / split_group),
                         [&](size_t g) {
                             const slong first = static_cast<slong>(g) * split_group;
                             for (slong b = first; b < std::min(blocks, first + split_group); ++b)
                                 acb_dot(terms_[b], nullptr, 0, polynomial->coeffs + b * block_, 1,
                                         powers_.data(), 1, std::min(block_, length - b * block_), precision_);
                         });
        Magnitude error, midpoint_norm, local_error;
        for (slong offset = ((length - 1) / block_) * block_; offset >= 0; offset -= block_) {
            if (split)
                acb_set(term_.value(), terms_[offset / block_]);
            else
                acb_dot(term_.value(), nullptr, 0, polynomial->coeffs + offset, 1,
                        powers_.data(), 1, std::min(block_, length - offset), precision_);

            // If v is in D(vhat,E) and T is in D(that,E_T), then
            // |v T - vhat that| <= E |T| + |vhat| E_T. The power table
            // includes argument uncertainty; no exact-point assumption is made.
            acb_get_mag(midpoint_norm.value(), out);
            mag_mul(error.value(), error.value(), factor_norm_.value());
            mag_addmul(error.value(), midpoint_norm.value(), factor_error_.value());
            acb_mul(out, out, factor_midpoint_.value(), precision_);
            acb_add(out, out, term_.value(), precision_);

            // The l1 radius dominates the Euclidean disk radius. Include
            // coefficient/power uncertainty and arithmetic rounding, then
            // recenter to avoid repeated rectangle rotation in the next block.
            mag_add(local_error.value(), arb_radref(acb_realref(out)), arb_radref(acb_imagref(out)));
            mag_add(error.value(), error.value(), local_error.value());
            acb_get_mid(out, out);
        }
        acb_add_error_mag(out, error.value());
    }

  public:
    OwnedDiskBlockedPolynomialEvaluator(Polynomial polynomial, slong precision)
        : polynomial_(std::move(polynomial)), precision_(precision), block_(block_size(polynomial_.value(), precision)),
          powers_(block_ + 1), terms_(block_count(polynomial_.value(), block_)) {
        acb_poly_derivative(derivative_.value(), polynomial_.value(), precision_);
        for (slong i = 0; i < acb_poly_length(polynomial_.value()); ++i)
            finite_ = finite_ && acb_is_finite(polynomial_.value()->coeffs + i);
    }
    OwnedDiskBlockedPolynomialEvaluator(const OwnedDiskBlockedPolynomialEvaluator &) = delete;
    OwnedDiskBlockedPolynomialEvaluator &operator=(const OwnedDiskBlockedPolynomialEvaluator &) = delete;
    OwnedDiskBlockedPolynomialEvaluator(OwnedDiskBlockedPolynomialEvaluator &&) noexcept = default;
    OwnedDiskBlockedPolynomialEvaluator &operator=(OwnedDiskBlockedPolynomialEvaluator &&) = delete;

    bool evaluate(acb_t value, acb_t derivative, const acb_t point) {
        // Invalid aliasing is refused before touching any caller-owned value.
        if (value == derivative || value == point || derivative == point) return false;
        if (!finite_ || !acb_is_finite(point) || powers_.size() != block_ + 1) return false;
        _acb_vec_set_powers(powers_.data(), point, powers_.size(), precision_);
        acb_get_mid(factor_midpoint_.value(), powers_[block_]);
        acb_get_mag(factor_norm_.value(), powers_[block_]);
        mag_add(factor_error_.value(), arb_radref(acb_realref(powers_[block_])),
                arb_radref(acb_imagref(powers_[block_])));
        evaluate_polynomial(value, polynomial_.value());
        evaluate_polynomial(derivative, derivative_.value());
        return acb_is_finite(value) && acb_is_finite(derivative);
    }
    // The value alone.
    bool evaluate_value(acb_t value, const acb_t point) {
        if (value == point) return false;
        if (!finite_ || !acb_is_finite(point) || powers_.size() != block_ + 1) return false;
        _acb_vec_set_powers(powers_.data(), point, powers_.size(), precision_);
        acb_get_mid(factor_midpoint_.value(), powers_[block_]);
        acb_get_mag(factor_norm_.value(), powers_[block_]);
        mag_add(factor_error_.value(), arb_radref(acb_realref(powers_[block_])),
                arb_radref(acb_imagref(powers_[block_])));
        evaluate_polynomial(value, polynomial_.value());
        return acb_is_finite(value);
    }
};
using DiskBlockedPolynomialEvaluator = OwnedDiskBlockedPolynomialEvaluator<ComplexPolynomial>;
using SharedDiskBlockedPolynomialEvaluator = OwnedDiskBlockedPolynomialEvaluator<PolynomialSnapshot>;
} // namespace certpolysolve::detail
