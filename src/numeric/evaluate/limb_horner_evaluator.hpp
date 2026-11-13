#pragma once
#include "core/raii.hpp"
#include <cstdint>
#include <limits>
#include <vector>

namespace certpolysolve::detail {
// A dense polynomial's value at a point whose two parts are short numbers (at most two limbs each: a binary64
// number, or the 115-bit nodes of a last regeneration), by Horner's rule on integer limbs, with its error as
// one disk radius.
//
// Why: a secular regeneration evaluates the polynomial at every node, at 117 to 1,872 bits (the precision
// follows the cancellation at the node), and in every regeneration but the last the nodes are binary64 numbers.
// The other two evaluators do not use that. The block dot products (DiskBlockedPolynomialEvaluator) multiply
// the coefficients by powers of the node, which are long numbers however short the node is. Above the degree
// the source's own evaluation is FLINT's ball Horner, whose work a term is mostly not limb arithmetic. With the
// node products in integers these evaluations are a large part of a hard secular solve. Horner's rule with a
// multiplier of one limb needs four passes over a long mantissa a term, and no ball. Its error is a disk, so it
// also loses nothing to the rectangles of ball Horner (half a bit a term).
// History: docs/code-history-20261021.md, limb_horner_evaluator.hpp 1.
//
// The point. x = (X_r + i X_i) 2^q with X_r, X_i integers: the two parts' mantissas (odd, at most 128 bits) at
// their common exponent, each held as at most three limbs at a limb offset (a part far above the other, as at
// a node beside the real axis, is a short multiplier placed high, not a long one). Where the real part is
// negative the evaluation runs at -x on the coefficients (-1)^k a_k, the same value, so that the real
// multiplier is never negative and one sign decides between adding and subtracting the cross products.
//
// The window. The running value t is (T_r + i T_i) 2^e, T_r and T_i two's complement integers of L limbs with
// |T| <= 2^capacity, capacity = 64 L - 2, in both parts. A step t <- t x + a_k forms the products exactly in
// L + m limbs, m the limbs the larger multiplier reaches (|X_r T_r| + |X_i T_i| < 2^(capacity + 1 + 64 m), which
// is what L + m limbs hold), and moves the window by whole limbs to where the larger of the product and the
// coefficient ends within one bit of the capacity (their sum is then within it): a floor of each part of the
// product when the window moves up, and a truncation toward zero of each part of the coefficient's midpoint at
// the window's unit. No bit is shifted in the product, and a step that drops nothing is exact.
//
// The precision. The window is placed by the larger of the bounds 2^s on the product and the coefficient,
// and s ends above capacity - 65; what that bound is on is at least half of it, so the larger operand of a
// step is at least 2^(capacity - 66) = 2^(64 (L - 1) - 4) units. L = limbs_for(precision) makes a unit at
// most 2^-(precision + guard_bits) of it: a step's error, below four units, is then below the 2^-precision of
// its operands that a ball operation at that precision may lose.
//
// The bound. A floor moves each part by less than a unit and so the value by less than sqrt(2) units, which
// is counted as truncation_units = 1.5; the coefficient's truncation the same; the coefficient's ball adds
// its two radii. With E the bound on |t - exact| before a step, after it E |x| + (dropped) + r_k, kept as a
// binary64 mantissa and a separate exponent with every operation rounded up. The result is the window's
// value with that radius in both coordinates. The radius is zero where nothing was dropped and the
// coefficients are exact.
class LimbHornerEvaluator {
    // One part of a coefficient's midpoint: mantissa 2^(top - 64 count), the mantissa `count` limbs with the
    // top bit of the last one set (FLINT's own storage, borrowed from polynomial_); count zero for a zero part.
    struct Term {
        const ulong *limbs = nullptr;
        int count = 0;
        bool negative = false;
        slong top = 0;
    };
    struct Coefficient {
        Term real, imaginary;
        // The ball's two radii added: radius 2^radius_exponent, zero for an exact coefficient.
        double radius = 0;
        slong radius_exponent = 0;
    };
    ComplexPolynomial polynomial_;
    std::vector<Coefficient> coefficients_;
    std::vector<ulong> real_[2], imaginary_[2], scratch_;
    int limbs_ = 0;
    bool usable_ = false, real_coefficients_ = true;

  public:
    // Bits a window keeps beyond the precision: a step loses less than four units (two truncations).
    static constexpr int guard_bits = 2;
    // Bits of a window's limbs that a value never reaches: the sign, and one so that the sum of the two cross
    // products of a step stays within the limbs a product has.
    static constexpr int headroom_bits = 2;
    // Bits between the window's capacity and the least size its larger operand has after a placement, beyond the
    // limb a placement moves by: one because the larger operand is placed a bit below the capacity (its sum with
    // the other then fits), one because an operand may be half the bound it is placed by.
    static constexpr int placement_bits = 2;
    // The largest window. Its limbs are counted in int, with a multiplier's reach (a few limbs) and the buffers'
    // slack added to the count, and its bits must stay far inside a long beside coefficient_exponent_range: half
    // of int (2^36 bits) leaves room for both.
    // History: docs/code-history-20261021.md, limb_horner_evaluator.hpp 2.
    static constexpr int maximum_limbs = std::numeric_limits<int>::max() / 2;
    // The largest precision limbs_for serves: the one whose window is maximum_limbs.
    static constexpr slong maximum_precision =
        slong{FLINT_BITS} * (maximum_limbs - 1) - (guard_bits + headroom_bits + placement_bits);
    // The modulus of a part's error after a floor of both coordinates, sqrt(2) units, rounded up.
    static constexpr double truncation_units = 1.5;
    // Limbs a part of the point may have.
    static constexpr int part_limbs = 2;
    // The largest exponent of two a part of the point may have in size (and the least, negated): the point's
    // modulus is then a binary64 number.
    static constexpr int exponent_range = 1000;
    // The largest exponent of two a coefficient's part or radius may have in size (and the least, negated),
    // so that every exponent of an evaluation stays far inside a long (a step moves the unit by at most the
    // point's exponent and the window's length).
    static constexpr slong coefficient_exponent_range = slong{1} << 40;

    // The window's limbs for a precision: the least L with 64 (L - 1) - 4 >= precision + guard_bits; zero for a
    // precision below one bit or above maximum_precision, declined before anything is added to it.
    static int limbs_for(slong precision);
    // Whether a point is one this evaluator serves: an exact ball, both parts zero or of at most part_limbs
    // limbs and within exponent_range.
    static bool serves(const acb_t point);

    // The polynomial's coefficients are balls of at most `precision` bits (longer midpoints are refused).
    LimbHornerEvaluator(ComplexPolynomial polynomial, slong precision);
    LimbHornerEvaluator(const LimbHornerEvaluator &) = delete;
    LimbHornerEvaluator &operator=(const LimbHornerEvaluator &) = delete;
    // Whether the coefficients are ones the evaluation takes: finite, at most the window's length less one
    // limb, exponents within coefficient_exponent_range; and the precision one limbs_for serves (nothing is
    // sized for one it declines).
    bool usable() const { return usable_; }
    // An enclosure of the polynomial's values at the point, over the coefficients' balls. False, with nothing
    // written, for a point or coefficients that are not served. Calls on one object must be exclusive.
    bool evaluate_value(acb_t value, const acb_t point);
};
} // namespace certpolysolve::detail
