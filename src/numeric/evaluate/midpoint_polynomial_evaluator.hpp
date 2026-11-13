#pragma once
#include "core/raii.hpp"
#include <flint/acf.h>
#include <algorithm>
#include <cmath>

namespace certpolysolve::detail {
// Candidate arithmetic only. A fixed-precision midpoint snapshot and reusable
// scratch; no returned value is an enclosure. Original-input certificates must
// check any proposals obtained from it. Sparse gap-Horner remains a separate
// representation, so this dense evaluator never expands sparse powers.
class MidpointPolynomialEvaluator {
    class Buffer {
        struct Clear {
            slong size;
            void operator()(acf_ptr p) const noexcept { if (p) _acf_vec_clear(p, size); }
        };
        slong size_;
        std::unique_ptr<acf_struct[], Clear> values_;
        static slong checked(slong size) {
            checked_vector_size(size, sizeof(acf_struct));
            if (size > std::numeric_limits<slong>::max() / 2)
                throw std::length_error("complex midpoint vector exceeds backend indexing");
            return size;
        }
      public:
        explicit Buffer(slong size) : size_(checked(size)), values_(_acf_vec_init(size_), Clear{size_}) {}
        slong size() const noexcept { return size_; }
        acf_ptr operator[](slong i) & {
            if (i < 0 || i >= size_) throw std::out_of_range("complex midpoint vector index");
            return values_.get() + i;
        }
        acf_srcptr operator[](slong i) const & {
            if (i < 0 || i >= size_) throw std::out_of_range("complex midpoint vector index");
            return values_.get() + i;
        }
    };
    static constexpr auto rounding = ARF_RND_NEAR;
    slong precision_, block_;
    Buffer coefficients_, derivative_, powers_;

    void evaluate_polynomial(acf_t out, const Buffer &coefficients) {
        arf_zero(acf_realref(out));
        arf_zero(acf_imagref(out));
        const slong length = coefficients.size();
        if (!length) return;
        // p(z) = P₀(z) + zᵇ P₁(z) + ... . A dot product accumulates each
        // block before rounding, and Horner combines the blocks using zᵇ.
        slong offset = ((length - 1) / block_) * block_;
        acf_approx_dot(out, nullptr, 0, coefficients[offset], 1, powers_[0], 1,
                       length - offset, precision_, rounding);
        while (offset >= block_) {
            offset -= block_;
            acf_mul(out, out, powers_[block_], precision_, rounding);
            acf_approx_dot(out, out, 0, coefficients[offset], 1, powers_[0], 1,
                           block_, precision_, rounding);
        }
    }
  public:
    MidpointPolynomialEvaluator(const acb_poly_t polynomial, slong precision)
        : precision_(precision),
          block_(std::max<slong>(1, static_cast<slong>(std::sqrt(
              static_cast<double>(acb_poly_length(polynomial)))))),
          coefficients_(acb_poly_length(polynomial)),
          derivative_(std::max<slong>(0, acb_poly_degree(polynomial))), powers_(block_ + 1) {
        if (precision < 2) throw std::invalid_argument("midpoint evaluation needs positive precision");
        for (slong j = 0; j < coefficients_.size(); ++j) {
            arf_set(acf_realref(coefficients_[j]), arb_midref(acb_realref(polynomial->coeffs + j)));
            arf_set(acf_imagref(coefficients_[j]), arb_midref(acb_imagref(polynomial->coeffs + j)));
            if (j) {
                arf_mul_ui(acf_realref(derivative_[j - 1]), acf_realref(coefficients_[j]),
                           static_cast<ulong>(j), precision_, rounding);
                arf_mul_ui(acf_imagref(derivative_[j - 1]), acf_imagref(coefficients_[j]),
                           static_cast<ulong>(j), precision_, rounding);
            }
        }
        arf_one(acf_realref(powers_[0]));
    }
    // Outputs must be distinct, and cannot alias the input point.
    void evaluate(acf_t value, acf_t derivative, const acf_t point) {
        acf_set(powers_[1], point);
        for (slong j = 2; j <= block_; ++j)
            acf_mul(powers_[j], powers_[j / 2], powers_[j - j / 2], precision_, rounding);
        evaluate_polynomial(value, coefficients_);
        evaluate_polynomial(derivative, derivative_);
    }
};
} // namespace certpolysolve::detail
