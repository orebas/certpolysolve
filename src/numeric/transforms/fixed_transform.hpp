#pragma once
#include "numeric/transforms/fft_backend.hpp"
#include "core/raii.hpp"
#include <cstdint>
#include <memory>
#include <optional>
#include <span>

namespace certpolysolve::detail {
// A radix-2 discrete Fourier transform in fixed-point integers of 64 K bits, with one error bound for the whole
// transform that is known before it runs. It computes what acb_dft_rad2_precomp computes,
//
//     X_k = sum_j x_j w^(j k),   w = exp(-2 pi i / N),   N = 2^L,
//
// and returns balls that contain every X_k for every choice of the x_j in their input balls.
//
// Why: FLINT's ball transform spends most of each butterfly on rounding, radii and normalisation whatever its
// precision, and the transform is a large part of a large piecewise solve. Integer butterflies need no rounding
// analysis of their own: every operation below is exact except three named truncations.
// History: docs/code-history-20261021.md, fixed_transform.hpp 1.
//
// Representation. Let 2^e exceed every midpoint component of the input (e from the components' exponents) and
// W = 64 K. A component v is held as the integer trunc(v 2^(F - e)), F = W - L - headroom_bits, in two's
// complement; one unit is u = 2^(e - F). A root of unity is held with G = W - 2 fraction bits, each component
// within twiddle_units of its exact value times 2^G (a truncated ball midpoint; Impl::table checks the ball when it
// builds a word size's table, at that size's first transform).
//
// Algorithm (decimation in frequency, in place; the result is in bit-reversed order): at level l the blocks have
// length M = N / 2^l, and for a pair (A, B) half a block apart and the root r = exp(-2 pi i q / M),
//
//     A' = A + B,   B' = floor((A - B) r 2^-G)
//
// with the sum, the difference and the complex product exact in integers (the product in 2 W bits) and the floor
// taken a component (an arithmetic shift).
//
// No overflow. Let R_l bound the modulus of every integer pair before level l; R_0 < sqrt(2) 2^F. A sum or
// difference has modulus at most 2 R_l, and |B'| <= 2 R_l (1 + sqrt(2) twiddle_units 2^-G) + sqrt(2), so
// R_L < 2^(F + L + 1) = 2^(W - 2): every component, sum and difference fits W bits, and every product of a
// difference (below 2^(W - 2)) with a root component (at most 2^G + twiddle_units) fits 2 W bits with its partner.
//
// Error. Let E_l bound the modulus of (held pair) - (exact pair of the same algorithm run on the midpoints), in
// units. E_0 <= sqrt(2): two truncations below one unit. A' is exact, so its error is at most 2 E_l. For B', with
// D = A - B as held, D* exact and r~ the held root,
//
//     |B' - D* r| <= |D| |r~ 2^-G - r| + |D - D*| |r| + sqrt(2) <= sqrt(2) twiddle_units + 2 E_l + sqrt(2),
//
// since |D| <= 2^G, |r~ - r 2^G| <= sqrt(2) twiddle_units, |r| = 1, and the two floors move a component by less
// than one unit each. So E_(l+1) <= 2 E_l + sqrt(2) (1 + twiddle_units) and
//
//     E_L < 2^L sqrt(2) (2 + twiddle_units) <= error_units_per_point N units.
//
// The input radii: the transform is linear and |w| = 1, so moving every x_j inside its ball moves an output by at
// most rho = sum_j (radius of x_j's real part + radius of its imaginary part). Each output is therefore within
//
//     bound = error_units_per_point N 2^(e - F) + rho
//
// of the held value times u, and a rectangle of that half-width in each coordinate contains it.
//
// Word size. apply() takes the smallest K whose arithmetic part of the bound is at most 2^(e - precision): one
// unit in the last place of the largest input component at the caller's working precision, which a ball transform
// at that precision cannot improve on at that scale. So the transform stands in for the ball transform at the
// same working precision; it does not trade accuracy for speed. Beyond maximum_limbs it declines.
struct FixedTransformStatistics {
    uint64_t accepted = 0, declined = 0;
    uint64_t limbs = 0; // of the accepted transform
    uint64_t conversion_ns = 0, execution_ns = 0, enclosure_ns = 0;
};
class FixedTransform {
    struct Impl;
    std::unique_ptr<Impl> impl_;

  public:
    // Spare leading bits of a held component: one for the sign, one for the modulus of a pair against its
    // components (sqrt(2)), one for the last level's sum (R_L < 2^(F + L + 1)).
    static constexpr int headroom_bits = 3;
    // A held root component is a truncated midpoint (below one unit) of a ball whose radius Impl::table checks
    // to be far below one unit (at most 2^-16 units; a table that fails declines its transforms): within two
    // units in all.
    static constexpr int twiddle_units = 2;
    // ceil(sqrt(2) (2 + twiddle_units)) = ceil(5.66).
    static constexpr int error_units_per_point = 6;
    static constexpr int minimum_limbs = 2, maximum_limbs = 8;
    // The smallest K in [minimum_limbs, maximum_limbs] with error_units_per_point 2^L 2^-F <= 2^-precision, that
    // is 64 K >= precision + 2 L + headroom_bits + ceil(log2(error_units_per_point)); zero if there is none (a
    // precision beyond 64 maximum_limbs is declined before any arithmetic on it).
    static int limbs_for(slong precision, int log2_length);
    // error_units_per_point 2^L 2^(e - F) as a power of two times a small integer: the arithmetic part of the
    // bound for inputs below 2^e, in the caller's unit. For tests and for the derivation's reader.
    static void arithmetic_bound(mag_t bound, slong input_exponent, int log2_length, int limbs);

    // Throws std::invalid_argument for an exponent outside [1, 31], before anything is computed from it.
    explicit FixedTransform(int log2_length);
    ~FixedTransform();
    FixedTransform(const FixedTransform &) = delete;
    FixedTransform &operator=(const FixedTransform &) = delete;
    // The transform of `input` (length 2^log2_length) as enclosures rounded to `precision` bits, or nothing when
    // an input is not finite or no supported word size reaches the precision.
    std::optional<EnclosedTransform> apply(const RootVector &input, slong precision, FixedTransformStatistics &stats);
    // The same with K words, whatever the precision asks for: for tests of each word size.
    std::optional<EnclosedTransform> apply_with(const RootVector &input, slong precision, int limbs,
                                                FixedTransformStatistics &stats);
    // apply() into the caller's vector, publishing only the outputs listed (distinct indices; every output when the
    // list is empty) and leaving the others as they were: a ring's models take one output a sector of interest, and
    // turning a fixed-point value into a ball costs as much as its share of the transform. False when declined.
    bool apply_into(RootVector &output, const RootVector &input, slong precision, std::span<const slong> wanted,
                    FixedTransformStatistics &stats);

  private:
    bool transform(RootVector &output, const RootVector &input, slong precision, int limbs,
                   std::span<const slong> wanted, FixedTransformStatistics &stats);
};
} // namespace certpolysolve::detail
