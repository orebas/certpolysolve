#pragma once
#include "core/number.hpp"
#include <algorithm>
#include <numeric>

namespace certpolysolve::detail {
// Value-only description of an owned input, independent of any solver. A
// snapshot must be rebuilt when its coefficient sequence changes. Compressed
// exponents are never expanded and no arithmetic objects are allocated here.
struct InputStructure {
    size_t start = 0, length = 0;
    uint64_t support = 0, support_stride = 0;
    bool exact = true, real = true;

    // Degree after removing an exact zero factor, before square-free algebra.
    // For uncertain leading coefficients this is an upper bound, not a proof
    // that every polynomial in the coefficient box has that degree.
    uint64_t degree() const { return length > start ? length - start - 1 : 0; }
    uint64_t nonzero_part_power() const {
        return exact ? std::max<uint64_t>(1, support_stride) : 1;
    }
    uint64_t reduced_degree() const { return degree() / nonzero_part_power(); }
    // The existing power-lifting path requires a nonzero constant. Zero roots
    // must be treated separately before using the more general x^s q(x^k).
    uint64_t unshifted_power() const {
        return start == 0 && length > 2 ? nonzero_part_power() : 1;
    }

    // Visit in increasing coefficient-index order. Skipped indices are exact
    // zeros; this also allows an importer to gather metadata as it copies.
    void observe(size_t index, const Coefficient &coefficient) {
        const bool coefficient_exact = coefficient.exact();
        const bool real_zero = coefficient.real.zero();
        const bool imag_zero = coefficient.imag.zero();
        observe(index, coefficient_exact, imag_zero && coefficient.imag_radius.zero(),
                coefficient_exact && real_zero && imag_zero);
    }
    // The exact integer/Gaussian adapters already know these properties and
    // need not allocate a Coefficient just to describe an owned factor.
    void observe(size_t index, bool coefficient_exact, bool coefficient_real, bool coefficient_zero) {
        exact = exact && coefficient_exact;
        real = real && coefficient_real;
        if (coefficient_exact && coefficient_zero) return;
        if (!support) start = index;
        ++support;
        length = index + 1;
        if (support_stride != 1)
            support_stride = std::gcd(support_stride, static_cast<uint64_t>(index - start));
    }
};

inline InputStructure input_structure(std::span<const Coefficient> coefficients) {
    InputStructure result;
    for (size_t i = 0; i < coefficients.size(); ++i) result.observe(i, coefficients[i]);
    return result;
}
// The same from a polynomial's nonzero terms (exponent, coefficient) in increasing exponent order, the input object's
// storage (C40): the exponents it skips are exact zeros, which observe() leaves out, so the structure is the one of
// the expanded coefficient list.
template <class Terms> InputStructure input_structure_of_terms(const Terms &terms) {
    InputStructure result;
    for (const auto &[exponent, coefficient] : terms) result.observe(static_cast<size_t>(exponent), coefficient);
    return result;
}
} // namespace certpolysolve::detail
