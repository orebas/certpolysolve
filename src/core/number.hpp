#pragma once
#include "core/raii.hpp"
#include "core/resource_limits.hpp"
#include <cstdint>
#include <string_view>
#include <span>

namespace certpolysolve::detail {
struct GaussianRational {
    Rational real, imag;
    bool zero() const { return fmpq_is_zero(real.value()) && fmpq_is_zero(imag.value()); }
};
// Exact value = significand * radix^exponent. Exponents stay compressed.
// Only an explicitly budgeted exact-algebra path materializes powers.
struct Number {
    Rational significand;
    Integer exponent;
    unsigned radix = 2;
    Number() = default;
    explicit Number(std::string_view source);
    static Number from_double(double x);
    bool zero() const { return fmpq_is_zero(significand.value()); }
    int sign() const { return fmpq_sgn(significand.value()); }
    void enclose(arb_t result, slong precision) const;
    uint64_t owned_storage_bits() const;
    uint64_t expanded_storage_bits() const;
    bool rational_with_budget(fmpq_t result, uint64_t& remaining_bits,
                              bool allow_exponent_expansion = true) const;
};
struct Coefficient {
    Number real, imag, real_radius, imag_radius;
    bool exact() const { return real_radius.zero() && imag_radius.zero(); }
    bool zero() const { return exact() && real.zero() && imag.zero(); }
    bool is_real() const { return imag.zero() && imag_radius.zero(); }
    void enclose(acb_t result, slong precision) const;
};
// The compact-exponent rule (docs/resource-policy.md): numbers written out as rationals may exceed the storage the
// input owns by at most `allowance` bits for each nonzero number written out. Exact preparation passes one precision
// cap (exact_expansion_allowance in query.hpp); no degree or absolute byte threshold.
bool expansion_fits(uint64_t owned, uint64_t expanded, uint64_t numbers, uint64_t allowance);
bool compact_expansion_fits(std::span<const Coefficient>, uint64_t allowance);
} // namespace certpolysolve::detail
