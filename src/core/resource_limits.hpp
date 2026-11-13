#pragma once
#include <flint/flint.h>
#include <algorithm>
#include <cstdint>
#include <concepts>
#include <limits>
#include <stdexcept>

namespace certpolysolve::detail::resources {
inline constexpr uint64_t unlimited = UINT64_MAX;
// WORD_MAX is FLINT's exact-precision sentinel, not a working precision.
inline constexpr uint64_t precision_max = static_cast<uint64_t>(WORD_MAX) - 1;
inline constexpr uint64_t add(uint64_t a, uint64_t b, uint64_t limit = unlimited) {
    return a > limit || b > limit - a ? limit : a + b;
}
inline constexpr uint64_t multiply(uint64_t a, uint64_t b, uint64_t limit = unlimited) {
    return b && a > limit / b ? limit : a * b;
}
template <std::integral T> uint64_t nonnegative(T value) {
    if constexpr (std::signed_integral<T>)
        if (value < 0) throw std::invalid_argument("negative size or working precision");
    return static_cast<uint64_t>(value);
}
template <std::integral A, std::integral B> slong precision_sum(A a, B b) {
    return static_cast<slong>(add(nonnegative(a), nonnegative(b), precision_max));
}
template <std::integral A, std::integral B> slong precision_product(A a, B b) {
    return static_cast<slong>(multiply(nonnegative(a), nonnegative(b), precision_max));
}
// A precision rounded up to the limbs it occupies: arithmetic at that many bits costs no more.
template <std::integral A> slong whole_limbs(A bits) {
    return precision_product(add(nonnegative(bits), FLINT_BITS - 1, precision_max) / FLINT_BITS, FLINT_BITS);
}
inline bool charge(uint64_t &remaining, uint64_t amount) {
    if (amount > remaining) return false;
    remaining -= amount;
    return true;
}
// Saturated estimates must still be refused at a representability boundary.
inline bool storage_fits(uint64_t estimate, uint64_t allowance) {
    return estimate != unlimited && estimate <= allowance &&
           estimate <= static_cast<uint64_t>(std::numeric_limits<ptrdiff_t>::max());
}
} // namespace certpolysolve::detail::resources
