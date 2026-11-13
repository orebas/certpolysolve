#pragma once
// 128-bit integers and checked 64-bit arithmetic on every compiler the library builds with.
//
// GCC and Clang have unsigned __int128 and __int128 and the built-ins __builtin_mul_overflow and
// __builtin_sub_overflow; MSVC has neither. The sites (fixed_transform.cpp, limb_horner_evaluator.cpp,
// node_products.cpp and .hpp, piecewise_models.cpp; candidates.cpp, root_index.hpp, machine_adapter.hpp) name the
// aliases and call the helpers of this header. Where the compiler has the native types and built-ins they are
// those, so the code the library compiles on Linux is the one before this header, byte for byte; elsewhere they
// are the portable classes below, two 64-bit limbs with exactly the two's complement semantics the sites rely on:
// construction from any integer with sign extension, truncating conversions to the integer types, arithmetic
// modulo 2^128, logical shifts of the unsigned type and arithmetic right shifts of the signed one, comparisons,
// a signed reading of an unsigned value and the reverse by explicit casts.
//
// CERTPOLYSOLVE_PORTABLE_INT128 (the CMake option of the same name) forces the portable classes and helpers on GCC
// and Clang, where tests/test_wide_int.cpp compares them with the native types and built-ins and the suite runs
// the library on them. A shift by an amount outside 0..127 is undefined for the native types; the portable
// classes move every bit out (zero, or the sign of the signed type).
#include <bit>
#include <cmath>
#include <compare>
#include <concepts>
#include <cstdint>
#include <limits>
#include <type_traits>
#if defined(_MSC_VER) && !defined(__clang__) && (defined(_M_X64) || defined(_M_ARM64))
#include <intrin.h>
#define CERTPOLYSOLVE_WIDE_INT_MSVC_INTRINSICS 1
#endif
#if defined(__GNUC__) || defined(__clang__)
#define CERTPOLYSOLVE_WIDE_INT_ALWAYS_INLINE [[gnu::always_inline]]
#else
#define CERTPOLYSOLVE_WIDE_INT_ALWAYS_INLINE
#endif

namespace certpolysolve::detail {
// A limb is 64 bits, and a wide integer is two limbs.
constexpr int wide_limb_bits = std::numeric_limits<std::uint64_t>::digits;
constexpr int uint128_bits = 2 * wide_limb_bits;

namespace portable {
class Int128;
// An unsigned 128-bit integer as unsigned __int128 behaves: arithmetic modulo 2^128, division, logical shifts, unsigned
// comparisons. Construction from an integer sign-extends it (a negative value reads as its two's complement);
// a conversion to an integer type keeps the low bits; a conversion to a floating type rounds once, to nearest.
class UInt128 {
    std::uint64_t low_ = 0, high_ = 0;
    static constexpr int half_bits = wide_limb_bits / 2;
    static constexpr std::uint64_t half_mask = (std::uint64_t{1} << half_bits) - 1;

    // The quotient of a by b, and the remainder where `remainder` is given (operator/ and operator%).
    static constexpr UInt128 divide(UInt128 a, UInt128 b, UInt128 *remainder) {
        if (a.high_ == 0 && b.high_ == 0) {
            if (remainder) *remainder = UInt128(a.low_ % b.low_);
            return UInt128(a.low_ / b.low_);
        }
        UInt128 q, r;
        for (int i = uint128_bits - 1; i >= 0; --i) {
            const bool carry = (r.high_ >> (wide_limb_bits - 1)) != 0;
            r = (r << 1) | ((a >> i) & UInt128(1));
            if (carry || r >= b) {
                r -= b;
                q |= UInt128(1) << i;
            }
        }
        if (remainder) *remainder = r;
        return q;
    }

    // The exact product of two limbs.
    static constexpr UInt128 limb_product(std::uint64_t a, std::uint64_t b) {
#ifdef CERTPOLYSOLVE_WIDE_INT_MSVC_INTRINSICS
        if (!std::is_constant_evaluated()) {
#if defined(_M_X64)
            std::uint64_t high = 0;
            const std::uint64_t low = _umul128(a, b, &high);
            return from_limbs(low, high);
#else
            return from_limbs(a * b, __umulh(a, b));
#endif
        }
#endif
        // Four half-limb products; the middle sum is at most 3 (2^32 - 1), no carry out of it.
        const std::uint64_t a0 = a & half_mask, a1 = a >> half_bits, b0 = b & half_mask, b1 = b >> half_bits;
        const std::uint64_t p00 = a0 * b0, p01 = a0 * b1, p10 = a1 * b0, p11 = a1 * b1;
        const std::uint64_t middle = (p00 >> half_bits) + (p01 & half_mask) + (p10 & half_mask);
        return from_limbs((middle << half_bits) | (p00 & half_mask),
                          p11 + (p01 >> half_bits) + (p10 >> half_bits) + (middle >> half_bits));
    }
    // A shift amount as the native types take it (0..127), every other amount moving every bit out.
    template <std::integral S>
    static constexpr int amount(S n) {
        if constexpr (std::is_signed_v<S>)
            if (n < 0) return uint128_bits;
        const auto count = static_cast<std::uint64_t>(n);
        return count >= static_cast<std::uint64_t>(uint128_bits) ? uint128_bits : static_cast<int>(count);
    }

  public:
    constexpr UInt128() = default;
    template <std::integral T>
    constexpr UInt128(T value) : low_(static_cast<std::uint64_t>(value)) { // NOLINT: implicit, as the native type
        if constexpr (std::is_signed_v<T>)
            if (value < 0) high_ = ~std::uint64_t{0};
    }
    constexpr explicit UInt128(Int128 value);
    static constexpr UInt128 from_limbs(std::uint64_t low, std::uint64_t high) {
        UInt128 r;
        r.low_ = low;
        r.high_ = high;
        return r;
    }
    constexpr std::uint64_t low() const { return low_; }
    constexpr std::uint64_t high() const { return high_; }

    template <std::integral T>
        requires(!std::same_as<T, bool>)
    constexpr explicit operator T() const {
        return static_cast<T>(low_);
    }
    constexpr explicit operator bool() const { return (low_ | high_) != 0; }
    // Correctly rounded (to nearest, ties to even), as the native conversion.
    template <std::floating_point T>
    explicit operator T() const {
        if (high_ == 0) return static_cast<T>(low_);
        if constexpr (std::numeric_limits<T>::digits >= wide_limb_bits) {
            // Each limb is exact in T, and their sum rounds once.
            constexpr T limb_unit = static_cast<T>(std::uint64_t{1} << (wide_limb_bits - 1)) * 2;
            return static_cast<T>(high_) * limb_unit + static_cast<T>(low_);
        } else {
            // The leading 64 bits with a sticky bit for the rest: T has fewer than 63 digits, so rounding that
            // limb to T rounds the whole value, and the power of two is exact.
            const int excess = static_cast<int>(std::bit_width(high_)); // bits above the leading 64: 1..64
            std::uint64_t leading = high_;
            bool sticky = low_ != 0;
            if (excess < wide_limb_bits) {
                leading = (high_ << (wide_limb_bits - excess)) | (low_ >> excess);
                sticky = (low_ & ((std::uint64_t{1} << excess) - 1)) != 0;
            }
            return std::ldexp(static_cast<T>(leading | static_cast<std::uint64_t>(sticky)), excess);
        }
    }

    friend constexpr bool operator==(const UInt128 &a, const UInt128 &b) = default;
    friend constexpr std::strong_ordering operator<=>(const UInt128 &a, const UInt128 &b) {
        if (a.high_ != b.high_) return a.high_ <=> b.high_;
        return a.low_ <=> b.low_;
    }
    friend constexpr UInt128 operator+(UInt128 a, UInt128 b) {
        const std::uint64_t low = a.low_ + b.low_;
        return from_limbs(low, a.high_ + b.high_ + static_cast<std::uint64_t>(low < a.low_));
    }
    friend constexpr UInt128 operator-(UInt128 a, UInt128 b) {
        return from_limbs(a.low_ - b.low_, a.high_ - b.high_ - static_cast<std::uint64_t>(a.low_ < b.low_));
    }
    friend constexpr UInt128 operator*(UInt128 a, UInt128 b) {
        const UInt128 low = limb_product(a.low_, b.low_);
        return from_limbs(low.low_, low.high_ + a.low_ * b.high_ + a.high_ * b.low_);
    }
    friend constexpr UInt128 operator&(UInt128 a, UInt128 b) { return from_limbs(a.low_ & b.low_, a.high_ & b.high_); }
    friend constexpr UInt128 operator|(UInt128 a, UInt128 b) { return from_limbs(a.low_ | b.low_, a.high_ | b.high_); }
    friend constexpr UInt128 operator^(UInt128 a, UInt128 b) { return from_limbs(a.low_ ^ b.low_, a.high_ ^ b.high_); }
    constexpr UInt128 operator-() const { return UInt128{} - *this; }
    constexpr UInt128 operator~() const { return from_limbs(~low_, ~high_); }
    template <std::integral S>
    constexpr UInt128 operator<<(S n) const {
        const int k = amount(n);
        if (k == 0) return *this;
        if (k >= uint128_bits) return {};
        if (k >= wide_limb_bits) return from_limbs(0, low_ << (k - wide_limb_bits));
        return from_limbs(low_ << k, (high_ << k) | (low_ >> (wide_limb_bits - k)));
    }
    template <std::integral S>
    constexpr UInt128 operator>>(S n) const {
        const int k = amount(n);
        if (k == 0) return *this;
        if (k >= uint128_bits) return {};
        if (k >= wide_limb_bits) return from_limbs(high_ >> (k - wide_limb_bits), 0);
        return from_limbs((low_ >> k) | (high_ << (wide_limb_bits - k)), high_ >> k);
    }
    constexpr UInt128 &operator+=(UInt128 b) { return *this = *this + b; }
    constexpr UInt128 &operator-=(UInt128 b) { return *this = *this - b; }
    constexpr UInt128 &operator*=(UInt128 b) { return *this = *this * b; }
    constexpr UInt128 &operator&=(UInt128 b) { return *this = *this & b; }
    constexpr UInt128 &operator|=(UInt128 b) { return *this = *this | b; }
    constexpr UInt128 &operator^=(UInt128 b) { return *this = *this ^ b; }
    // Division and remainder, as the native type's: the divisor must not be zero, as there. Where both operands fit a
    // limb, the limb's own; otherwise restoring long division, one quotient bit a step from the top. The remainder
    // stays below the divisor, so after its shift it holds at most 129 bits: a bit shifted out of it means it exceeds
    // every divisor, and the subtraction, taken modulo 2^128, is due.
    friend constexpr UInt128 operator/(UInt128 a, UInt128 b) { return divide(a, b, nullptr); }
    friend constexpr UInt128 operator%(UInt128 a, UInt128 b) {
        UInt128 r;
        divide(a, b, &r);
        return r;
    }
    constexpr UInt128 &operator/=(UInt128 b) { return *this = *this / b; }
    constexpr UInt128 &operator%=(UInt128 b) { return *this = *this % b; }
    template <std::integral S>
    constexpr UInt128 &operator<<=(S n) {
        return *this = *this << n;
    }
    template <std::integral S>
    constexpr UInt128 &operator>>=(S n) {
        return *this = *this >> n;
    }
    // The arithmetic right shift of the bits read as a signed integer (Int128's operator>>).
    template <std::integral S>
    constexpr UInt128 shifted_right_arithmetic(S n) const {
        const int k = amount(n);
        const std::uint64_t fill = (high_ >> (wide_limb_bits - 1)) ? ~std::uint64_t{0} : 0;
        if (k == 0) return *this;
        if (k >= uint128_bits) return from_limbs(fill, fill);
        const auto signed_high = static_cast<std::int64_t>(high_);
        if (k >= wide_limb_bits)
            return from_limbs(static_cast<std::uint64_t>(signed_high >> (k - wide_limb_bits)), fill);
        return from_limbs((low_ >> k) | (high_ << (wide_limb_bits - k)), static_cast<std::uint64_t>(signed_high >> k));
    }
};

// A signed 128-bit integer in two's complement as __int128 behaves: the bits of a UInt128, read with bit 127 as
// the sign; arithmetic modulo 2^128 (an overflow wraps, where the native type's is undefined), arithmetic right
// shifts, signed comparisons.
class Int128 {
    UInt128 bits_;

  public:
    constexpr Int128() = default;
    template <std::integral T>
    constexpr Int128(T value) : bits_(value) {} // NOLINT: implicit, as the native type
    constexpr explicit Int128(UInt128 bits) : bits_(bits) {}
    static constexpr Int128 from_limbs(std::uint64_t low, std::uint64_t high) { return Int128(UInt128::from_limbs(low, high)); }
    constexpr UInt128 bits() const { return bits_; }
    constexpr std::uint64_t low() const { return bits_.low(); }
    constexpr std::uint64_t high() const { return bits_.high(); }
    constexpr bool negative() const { return (bits_.high() >> (wide_limb_bits - 1)) != 0; }

    template <std::integral T>
        requires(!std::same_as<T, bool>)
    constexpr explicit operator T() const {
        return static_cast<T>(bits_.low());
    }
    constexpr explicit operator bool() const { return static_cast<bool>(bits_); }
    template <std::floating_point T>
    explicit operator T() const {
        return negative() ? -static_cast<T>(-bits_) : static_cast<T>(bits_);
    }

    friend constexpr bool operator==(const Int128 &a, const Int128 &b) = default;
    friend constexpr std::strong_ordering operator<=>(const Int128 &a, const Int128 &b) {
        const auto a_high = static_cast<std::int64_t>(a.high()), b_high = static_cast<std::int64_t>(b.high());
        if (a_high != b_high) return a_high <=> b_high;
        return a.low() <=> b.low();
    }
    friend constexpr Int128 operator+(Int128 a, Int128 b) { return Int128(a.bits_ + b.bits_); }
    friend constexpr Int128 operator-(Int128 a, Int128 b) { return Int128(a.bits_ - b.bits_); }
    friend constexpr Int128 operator*(Int128 a, Int128 b) { return Int128(a.bits_ * b.bits_); }
    friend constexpr Int128 operator&(Int128 a, Int128 b) { return Int128(a.bits_ & b.bits_); }
    friend constexpr Int128 operator|(Int128 a, Int128 b) { return Int128(a.bits_ | b.bits_); }
    friend constexpr Int128 operator^(Int128 a, Int128 b) { return Int128(a.bits_ ^ b.bits_); }
    constexpr Int128 operator-() const { return Int128(-bits_); }
    constexpr Int128 operator~() const { return Int128(~bits_); }
    template <std::integral S>
    constexpr Int128 operator<<(S n) const {
        return Int128(bits_ << n);
    }
    template <std::integral S>
    constexpr Int128 operator>>(S n) const {
        return Int128(bits_.shifted_right_arithmetic(n));
    }
    constexpr Int128 &operator+=(Int128 b) { return *this = *this + b; }
    constexpr Int128 &operator-=(Int128 b) { return *this = *this - b; }
    constexpr Int128 &operator*=(Int128 b) { return *this = *this * b; }
    constexpr Int128 &operator&=(Int128 b) { return *this = *this & b; }
    constexpr Int128 &operator|=(Int128 b) { return *this = *this | b; }
    constexpr Int128 &operator^=(Int128 b) { return *this = *this ^ b; }
    template <std::integral S>
    constexpr Int128 &operator<<=(S n) {
        return *this = *this << n;
    }
    template <std::integral S>
    constexpr Int128 &operator>>=(S n) {
        return *this = *this >> n;
    }
};
constexpr UInt128::UInt128(Int128 value) : UInt128(value.bits()) {}

// The checked arithmetic of the GCC built-ins on the signed types of at most 64 bits: *out is the result
// modulo the type's range, and the return is whether the exact result does not fit it.
template <std::signed_integral T>
constexpr bool multiply_overflows(T a, T b, T *out) {
    static_assert(std::numeric_limits<T>::digits < wide_limb_bits, "the exact product of two values fits 128 bits");
    const Int128 product = Int128(a) * Int128(b);
    *out = static_cast<T>(product);
    return product != Int128(*out);
}
template <std::signed_integral T>
constexpr bool subtract_overflows(T a, T b, T *out) {
    using Unsigned = std::make_unsigned_t<T>;
    const T difference = static_cast<T>(static_cast<Unsigned>(a) - static_cast<Unsigned>(b));
    *out = difference;
    // The exact difference leaves the range only when the operands have different signs and the result has b's.
    return ((a < 0) != (b < 0)) && ((difference < 0) == (b < 0));
}
} // namespace portable

#if defined(__SIZEOF_INT128__) && !defined(CERTPOLYSOLVE_PORTABLE_INT128)
__extension__ using uint128 = unsigned __int128;
__extension__ using int128 = __int128;
#else
using uint128 = portable::UInt128;
using int128 = portable::Int128;
#endif

#if !defined(CERTPOLYSOLVE_PORTABLE_INT128) && defined(__has_builtin)
#if __has_builtin(__builtin_mul_overflow) && __has_builtin(__builtin_sub_overflow)
#define CERTPOLYSOLVE_WIDE_INT_OVERFLOW_BUILTINS 1
#endif
#endif
#ifdef CERTPOLYSOLVE_WIDE_INT_OVERFLOW_BUILTINS
template <std::signed_integral T>
CERTPOLYSOLVE_WIDE_INT_ALWAYS_INLINE inline bool multiply_overflows(T a, T b, T *out) {
    return __builtin_mul_overflow(a, b, out);
}
template <std::signed_integral T>
CERTPOLYSOLVE_WIDE_INT_ALWAYS_INLINE inline bool subtract_overflows(T a, T b, T *out) {
    return __builtin_sub_overflow(a, b, out);
}
#else
using portable::multiply_overflows;
using portable::subtract_overflows;
#endif
} // namespace certpolysolve::detail
