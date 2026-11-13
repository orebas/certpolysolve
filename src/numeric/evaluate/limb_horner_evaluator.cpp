// LimbHornerEvaluator (limb_horner_evaluator.hpp derives its window and bound): a dense complex polynomial's value
// at an exact point with short parts (serves), by Horner's rule on integer limbs with a disk error (evaluate_value).
//
// The constructor borrows each coefficient's midpoint limbs from FLINT (Term) and sums its two radii into a Bound.
// evaluate_value reads the point (multiplier, from short_part and factor) as x = (X_r + i X_i) 2^q, X_r and X_i
// integers of at most three limbs at a limb offset; where the real part is negative it runs at -x on the coefficients
// (-1)^k a_k. A Horner step t <- t x + a_k holds each part of t as a two's complement integer of L limbs
// (limbs_for(precision)) times 2^unit: it forms t x exactly in L + m limbs, m the multiplier's reach (multiply and its
// variants, or accumulate), moves the window by whole limbs to fit the larger of the product and the coefficient
// (size_bits; place, a floor), and adds the coefficient's midpoint truncated toward zero at the window's unit
// (add_part). The error bound is a Bound (a binary64 mantissa with a separate exponent, rounded up by rounding_slack):
// a step scales it by a bound on |x| and adds the coefficient's radius and truncation_units units for each floor or
// truncation that dropped bits. The integers are exact, so the midpoint is the same on every CPU.
//
// State: the owned polynomial_, whose limbs coefficients_ borrows, and the window buffers real_, imaginary_ and
// scratch_, which every evaluate_value call reuses (so calls on one object must be exclusive).
// Entry points: serves, the constructor, usable and evaluate_value, from DiskValues in secular_refinement.cpp
// (choose_limb_horner, value): regenerate_residues' values of a dense source at short nodes, under
// RegenerationPolicy::limb_horner (session.cpp's maintained_method_policy sets it for the secular method).
// Not here: the other disk evaluators (disk_blocked_polynomial_evaluator.hpp, disk_sparse_polynomial_evaluator.hpp),
// the choice among them and the division by the leading coefficient (DiskValues).
#include "numeric/evaluate/limb_horner_evaluator.hpp"
#include "core/wide_int.hpp"
#include <algorithm>
#include <bit>
#include <cmath>
#include <limits>

namespace certpolysolve::detail {
namespace {
using Wide = uint128;
constexpr int limb_bits = 64;
static_assert(FLINT_BITS == limb_bits, "the evaluation is written for 64-bit limbs");
constexpr int headroom_bits = LimbHornerEvaluator::headroom_bits, placement_bits = LimbHornerEvaluator::placement_bits;
static_assert((LimbHornerEvaluator::maximum_precision + LimbHornerEvaluator::guard_bits + headroom_bits + placement_bits +
               limb_bits - 1) / limb_bits + 1 == LimbHornerEvaluator::maximum_limbs,
              "maximum_precision is the precision of the largest window");

slong
ceil_div(slong a, slong b) { // b > 0
    return a >= 0 ? (a + b - 1) / b : -((-a) / b);
}

// 1 + 2^-50: more than the three roundings to nearest (each at most 1 + 2^-53) an operation below makes.
constexpr double rounding_slack = 1 + 0x1p-50;
// A nonnegative bound mantissa 2^exponent, the mantissa zero or in [1, 2^stride), every operation rounded up.
struct Bound {
    static constexpr int stride = 256;
    // Exponents farther apart than this are not added: the smaller term is below 2^-(reach - stride) of the
    // larger, which the rounding slack covers.
    static constexpr slong reach = 400;
    double mantissa = 0;
    slong exponent = 0;
    void normalize() {
        while (mantissa >= 0x1p256) {
            mantissa *= 0x1p-256;
            exponent += stride;
        }
    }
    // Times factor 2^shift, factor at least one.
    void scale(double factor, slong shift) {
        mantissa = mantissa * factor * rounding_slack;
        exponent += shift;
        normalize();
    }
    // Plus amount 2^amount_exponent, amount zero or in [1, 2^stride).
    void add(double amount, slong amount_exponent) {
        if (amount == 0) return;
        const slong gap = amount_exponent - exponent;
        if (mantissa == 0 || gap > reach) {
            mantissa = mantissa == 0 ? amount : amount * rounding_slack;
            exponent = amount_exponent;
        } else if (gap < -reach)
            mantissa *= rounding_slack;
        else
            mantissa = (mantissa + amount * power_of_two(gap)) * rounding_slack;
        normalize();
    }
    // 2^e for an exponent within reach: the binary64 number with that biased exponent and a zero fraction.
    static double power_of_two(slong e) {
        constexpr int fraction_bits = std::numeric_limits<double>::digits - 1;
        constexpr slong bias = std::numeric_limits<double>::max_exponent - 1;
        return std::bit_cast<double>(static_cast<uint64_t>(e + bias) << fraction_bits);
    }
};
static_assert(Bound::reach + Bound::stride < std::numeric_limits<double>::max_exponent - Bound::stride,
              "a sum of two bounds within reach of each other is a binary64 number");

// A part of the point as magnitude 2^exponent, the magnitude odd and at most part_limbs limbs (zero for zero).
struct Part {
    Wide magnitude = 0;
    bool negative = false;
    slong exponent = 0;
    int bits() const {
        const auto high = static_cast<ulong>(magnitude >> limb_bits), low = static_cast<ulong>(magnitude);
        return high ? limb_bits + std::bit_width(high) : std::bit_width(low);
    }
};
static_assert(LimbHornerEvaluator::part_limbs * limb_bits == uint128_bits,
              "a part's magnitude is held in one double-limb integer");
bool
short_part(Part &out, const arf_t x) {
    out = {};
    if (arf_is_zero(x)) return true;
    if (!arf_is_finite(x) || arf_bits(x) > LimbHornerEvaluator::part_limbs * limb_bits) return false;
    const slong size = arf_abs_bound_lt_2exp_si(x);
    if (size > LimbHornerEvaluator::exponent_range || size < -LimbHornerEvaluator::exponent_range) return false;
    // One or two limbs, the top bit of the last set and the first nonzero: the value is those limbs times
    // 2^(size - 64 count).
    mp_srcptr limbs;
    mp_size_t count;
    ARF_GET_MPN_READONLY(limbs, count, x);
    const Wide mantissa = count == 2 ? (Wide{limbs[1]} << limb_bits) | limbs[0] : Wide{limbs[0]};
    const int trailing = std::countr_zero(limbs[0]);
    out.magnitude = mantissa >> trailing;
    out.exponent = size - limb_bits * static_cast<slong>(count) + trailing;
    out.negative = ARF_SGNBIT(x);
    return true;
}
// One part of the multiplier: sum limb[j] 2^(64 (offset + j)) over `count` limbs, in units of the point's
// common exponent.
struct Factor {
    static constexpr int most = LimbHornerEvaluator::part_limbs + 1; // a magnitude shifted by less than a limb
    ulong limb[most] = {0, 0, 0};
    int count = 0, offset = 0;
    int reach() const { return offset + count; }
};
Factor
factor(const Part &part, slong exponent) {
    Factor out;
    if (!part.magnitude) return out;
    const slong shift = part.exponent - exponent; // not negative
    out.offset = static_cast<int>(shift / limb_bits);
    const int bits = static_cast<int>(shift % limb_bits);
    const auto low = static_cast<ulong>(part.magnitude), high = static_cast<ulong>(part.magnitude >> limb_bits);
    out.limb[0] = low << bits;
    out.limb[1] = bits ? (high << bits) | (low >> (limb_bits - bits)) : high;
    out.limb[2] = bits ? high >> (limb_bits - bits) : 0;
    out.count = out.limb[2] ? 3 : out.limb[1] ? 2 : 1;
    return out;
}
// The point as (real + i imaginary) 2^exponent, or its negative where `negated` (the header's "the point"):
// the real factor is the real part's size, and the imaginary one carries its sign apart.
struct Multiplier {
    Factor real, imaginary;
    bool imaginary_negative = false, negated = false, zero = false;
    slong exponent = 0;
    // |x| <= modulus 2^modulus_exponent, modulus in [1, 2] and a little.
    double modulus = 0;
    slong modulus_exponent = 0;
    // Limbs a product has beyond the window's.
    int reach() const { return std::max(real.reach(), imaginary.reach()); }
    // Both factors one limb at the common exponent: a product is two passes.
    bool single() const { return real.reach() <= 1 && imaginary.reach() <= 1; }
};
bool
multiplier(Multiplier &out, const acb_t point) {
    out = {};
    if (!acb_is_exact(point)) return false;
    Part real, imaginary;
    if (!short_part(real, arb_midref(acb_realref(point))) || !short_part(imaginary, arb_midref(acb_imagref(point))))
        return false;
    if (!real.magnitude && !imaginary.magnitude) {
        out.zero = true;
        return true;
    }
    out.exponent = !real.magnitude        ? imaginary.exponent
                   : !imaginary.magnitude ? real.exponent
                                          : std::min(real.exponent, imaginary.exponent);
    if (real.magnitude && imaginary.magnitude) {
        // Where the higher part would straddle two limbs, a common exponent a little lower puts it on a limb
        // boundary, if the lower part still fits its limb.
        const Part &higher = real.exponent > imaginary.exponent ? real : imaginary;
        const Part &lower = real.exponent > imaginary.exponent ? imaginary : real;
        const int within = static_cast<int>((higher.exponent - out.exponent) % limb_bits);
        if (within && within + higher.bits() > limb_bits && higher.bits() <= limb_bits &&
            lower.bits() + (limb_bits - within) <= limb_bits)
            out.exponent -= limb_bits - within;
    }
    out.real = factor(real, out.exponent);
    out.imaginary = factor(imaginary, out.exponent);
    out.negated = real.magnitude && real.negative;
    out.imaginary_negative = imaginary.magnitude && imaginary.negative != out.negated;
    // The modulus from the parts rounded to binary64: two roundings and the hypotenuse's, within the slack.
    const double size = std::hypot(arf_get_d(arb_midref(acb_realref(point)), ARF_RND_NEAR),
                                   arf_get_d(arb_midref(acb_imagref(point)), ARF_RND_NEAR));
    int exponent = 0;
    out.modulus = 2 * std::frexp(size, &exponent) * rounding_slack;
    out.modulus_exponent = exponent - 1;
    return true;
}

// Products of a two's complement integer of `limbs` limbs with a limb, as two's complement integers of one
// limb more: GMP's unsigned product reads a negative t as t + 2^(64 limbs), so u leaves the top limb once.
inline void
multiply(ulong *out, const ulong *t, int limbs, ulong u) {
    const ulong high = mpn_mul_1(out, t, limbs, u);
    out[limbs] = high - ((t[limbs - 1] >> (limb_bits - 1)) ? u : 0);
}
inline void
multiply_add(ulong *out, const ulong *t, int limbs, ulong u) {
    const ulong high = mpn_addmul_1(out, t, limbs, u);
    out[limbs] += high - ((t[limbs - 1] >> (limb_bits - 1)) ? u : 0);
}
inline void
multiply_subtract(ulong *out, const ulong *t, int limbs, ulong u) {
    const ulong high = mpn_submul_1(out, t, limbs, u);
    out[limbs] -= high - ((t[limbs - 1] >> (limb_bits - 1)) ? u : 0);
}
// The same into a longer two's complement integer of `total` limbs at a limb position: out += t u 2^(64
// position), or -= where `subtract`. The carry and the correction for a negative t run up through the limbs
// above (position + limbs < total).
inline void
accumulate(ulong *out, int total, const ulong *t, int limbs, ulong u, int position, bool subtract) {
    ulong *at = out + position, *above = at + limbs;
    const int rest = total - position - limbs;
    const bool negative = t[limbs - 1] >> (limb_bits - 1);
    if (subtract) {
        const ulong borrow = mpn_submul_1(at, t, limbs, u);
        mpn_sub_1(above, above, rest, borrow);
        if (negative) mpn_add_1(above, above, rest, u);
    } else {
        const ulong carry = mpn_addmul_1(at, t, limbs, u);
        mpn_add_1(above, above, rest, carry);
        if (negative) mpn_sub_1(above, above, rest, u);
    }
}
inline void
accumulate(ulong *out, int total, const ulong *t, int limbs, const Factor &factor, bool subtract) {
    for (int j = 0; j < factor.count; ++j)
        if (factor.limb[j]) accumulate(out, total, t, limbs, factor.limb[j], factor.offset + j, subtract);
}
// An exponent b with |x| <= 2^b for a two's complement integer of `count` limbs (zero for zero): its limbs
// less the leading bits that repeat the sign.
slong
size_bits(const ulong *x, int count) {
    const ulong sign = x[count - 1] >> (limb_bits - 1);
    const ulong extension = sign ? ~ulong{0} : 0;
    int top = count - 1;
    while (top > 0 && x[top] == extension && (x[top - 1] >> (limb_bits - 1)) == sign) --top;
    const ulong pattern = x[top] ^ extension;
    return slong{limb_bits} * (top + 1) - (pattern ? std::countl_zero(pattern) : limb_bits);
}
inline bool
is_zero(const ulong *x, int limbs) {
    return x[limbs - 1] == 0 && mpn_zero_p(x, limbs);
}
// floor(product / 2^(64 move)) as `limbs` limbs, for a product of `total` limbs that fits them after the
// move. A move within the product is a view of it, and so is a move of one limb down: the limb below a
// product is slack, which that move zeroes.
ulong *
place(ulong *product, int limbs, int total, slong move) {
    if (move >= 0 && move <= total - limbs) return product + move;
    if (move == -1) {
        product[-1] = 0;
        return product - 1;
    }
    const ulong extension = (product[total - 1] >> (limb_bits - 1)) ? ~ulong{0} : 0;
    if (move > 0) {
        for (slong j = 0; j < limbs; ++j) product[j] = j + move < total ? product[j + move] : extension;
    } else {
        for (slong j = limbs - 1; j >= 0; --j) product[j] = j + move >= 0 ? product[j + move] : 0;
    }
    return product;
}
// Adds a part of a coefficient's midpoint to the window (subtracts it where `subtract`), truncated toward
// zero at the window's unit 2^unit. The window is placed so that the part ends within it. Whether bits of the
// part were below the unit (a mantissa's lowest limb is not zero, so whole limbs below it are).
template <class Term>
bool
add_part(ulong *window, int limbs, const Term &term, bool subtract, slong unit, ulong *scratch) {
    if (!term.count) return false;
    const slong shift = term.top - slong{limb_bits} * term.count - unit; // of the mantissa's unit against the window's
    const ulong *source = term.limbs;
    slong count = term.count, offset = 0;
    bool below = false;
    if (shift >= 0) {
        offset = shift / limb_bits;
        const int bits = static_cast<int>(shift % limb_bits);
        if (bits) {
            scratch[count] = mpn_lshift(scratch, source, count, bits);
            source = scratch;
            ++count;
        }
    } else {
        const slong down = -shift, skip = down / limb_bits;
        const int bits = static_cast<int>(down % limb_bits);
        if (skip >= count) return true; // the whole part is below one unit
        below = skip > 0 || (source[0] & ((ulong{1} << bits) - 1)) != 0;
        source += skip;
        count -= skip;
        if (bits) {
            mpn_rshift(scratch, source, count, bits);
            source = scratch;
        }
    }
    // Limbs of the shifted part beyond the window are zero by the placement.
    const slong used = std::min(count, limbs - offset);
    if (used <= 0) return below;
    const ulong carry = subtract ? mpn_sub_n(window + offset, window + offset, source, used)
                                 : mpn_add_n(window + offset, window + offset, source, used);
    const slong above = limbs - offset - used;
    if (above > 0) {
        if (subtract)
            mpn_sub_1(window + offset + used, window + offset + used, above, carry);
        else
            mpn_add_1(window + offset + used, window + offset + used, above, carry);
    }
    return below;
}
} // namespace

int
LimbHornerEvaluator::limbs_for(slong precision) {
    // Zero outside [1, maximum_precision]: the sum below then stays within slong and the count within
    // int. History: docs/code-history-20261021.md, limb_horner_evaluator.cpp 1.
    if (precision < 1 || precision > maximum_precision) return 0;
    const slong bits = precision + guard_bits + headroom_bits + placement_bits;
    return static_cast<int>((bits + limb_bits - 1) / limb_bits) + 1;
}
bool
LimbHornerEvaluator::serves(const acb_t point) {
    Multiplier unused;
    return multiplier(unused, point);
}
LimbHornerEvaluator::LimbHornerEvaluator(ComplexPolynomial polynomial, slong precision)
    : polynomial_(std::move(polynomial)), limbs_(limbs_for(precision)) {
    const acb_poly_struct *p = polynomial_.value();
    const slong length = acb_poly_length(p);
    // A precision limbs_for declines sizes nothing.
    if (!limbs_ || length < 1) return;
    coefficients_.resize(static_cast<size_t>(length));
    const auto term = [&](Term &out, const arf_t x) {
        if (arf_is_zero(x)) return true;
        if (!arf_is_finite(x) || COEFF_IS_MPZ(*ARF_EXPREF(x))) return false;
        const slong top = ARF_EXP(x);
        if (top > coefficient_exponent_range || top < -coefficient_exponent_range) return false;
        mp_srcptr limbs;
        mp_size_t count;
        ARF_GET_MPN_READONLY(limbs, count, x);
        if (count > limbs_ - 1) return false;
        out.limbs = limbs;
        out.count = static_cast<int>(count);
        out.negative = ARF_SGNBIT(x);
        out.top = top;
        return true;
    };
    const auto radius = [&](Bound &sum, const mag_t r) {
        if (mag_is_zero(r)) return true;
        if (mag_is_inf(r) || COEFF_IS_MPZ(*MAG_EXPREF(r))) return false;
        const slong exponent = MAG_EXP(r);
        if (exponent > coefficient_exponent_range || exponent < -coefficient_exponent_range) return false;
        // A magnitude is its 30-bit mantissa times 2^(exponent - 30).
        sum.add(static_cast<double>(MAG_MAN(r)), exponent - MAG_BITS);
        return true;
    };
    for (slong k = 0; k < length; ++k) {
        const acb_struct *a = p->coeffs + k;
        Coefficient &c = coefficients_[static_cast<size_t>(k)];
        Bound sum;
        if (!term(c.real, arb_midref(acb_realref(a))) || !term(c.imaginary, arb_midref(acb_imagref(a))) ||
            !radius(sum, arb_radref(acb_realref(a))) || !radius(sum, arb_radref(acb_imagref(a))))
            return;
        c.radius = sum.mantissa;
        c.radius_exponent = sum.exponent;
        if (c.imaginary.count || !mag_is_zero(arb_radref(acb_imagref(a)))) real_coefficients_ = false;
    }
    scratch_.assign(static_cast<size_t>(limbs_) + 2, 0);
    usable_ = true;
}

bool
LimbHornerEvaluator::evaluate_value(acb_t value, const acb_t point) {
    Multiplier m;
    if (!usable_ || !multiplier(m, point)) return false;
    if (m.zero) { // the constant term
        acb_set(value, polynomial_.value()->coeffs);
        return true;
    }
    const int limbs = limbs_;
    // A product has the window's limbs and the multiplier's reach; it is written one limb into its buffer
    // (the slack a view one limb down takes).
    const int total = limbs + m.reach();
    for (auto *buffer : {&real_[0], &real_[1], &imaginary_[0], &imaginary_[1]})
        if (buffer->size() < static_cast<size_t>(total) + 1) buffer->assign(static_cast<size_t>(total) + 1, 0);
    const slong capacity = slong{limb_bits} * limbs - headroom_bits;
    // A real polynomial at a real point: the imaginary part is zero throughout.
    const bool real_only = real_coefficients_ && !m.imaginary.count;
    const bool single = m.single();
    constexpr slong nothing = std::numeric_limits<slong>::min();
    int side = 0;
    ulong *tr = real_[0].data() + 1, *ti = imaginary_[0].data() + 1;
    bool zero = true; // the running value is exactly zero (its window is then not read)
    slong unit = 0;
    Bound error;
    for (slong k = static_cast<slong>(coefficients_.size()) - 1; k >= 0; --k) {
        const Coefficient &c = coefficients_[static_cast<size_t>(k)];
        ulong *pr = real_[1 - side].data() + 1, *pi = imaginary_[1 - side].data() + 1;
        const slong product_unit = unit + m.exponent;
        slong size = nothing;
        error.scale(m.modulus, m.modulus_exponent);
        if (!zero) {
            size = 0;
            if (single) {
                multiply(pr, tr, limbs, m.real.limb[0]);
                if (!real_only) {
                    multiply(pi, ti, limbs, m.real.limb[0]);
                    if (m.imaginary.count) {
                        if (m.imaginary_negative) {
                            multiply_add(pr, ti, limbs, m.imaginary.limb[0]);
                            multiply_subtract(pi, tr, limbs, m.imaginary.limb[0]);
                        } else {
                            multiply_subtract(pr, ti, limbs, m.imaginary.limb[0]);
                            multiply_add(pi, tr, limbs, m.imaginary.limb[0]);
                        }
                    }
                }
            } else {
                mpn_zero(pr, total);
                accumulate(pr, total, tr, limbs, m.real, false);
                if (!real_only) {
                    mpn_zero(pi, total);
                    accumulate(pr, total, ti, limbs, m.imaginary, !m.imaginary_negative);
                    accumulate(pi, total, ti, limbs, m.real, false);
                    accumulate(pi, total, tr, limbs, m.imaginary, m.imaginary_negative);
                }
            }
            if (!real_only) size = size_bits(pi, total);
            size = std::max(size, size_bits(pr, total));
        }
        if (c.real.count) size = std::max(size, c.real.top - product_unit);
        if (c.imaginary.count) size = std::max(size, c.imaginary.top - product_unit);
        if (size == nothing) { // a zero midpoint on a zero value
            if (c.radius != 0) error.add(c.radius, c.radius_exponent);
            unit = product_unit;
            continue;
        }
        // The window moves so that the larger bound ends within a bit of the capacity.
        const slong move = ceil_div(size - (capacity - 1), limb_bits);
        const slong placed = product_unit + limb_bits * move;
        bool dropped = false;
        if (zero) {
            mpn_zero(pr, limbs);
            if (!real_only) mpn_zero(pi, limbs);
            tr = pr;
            ti = pi;
        } else {
            // A floor drops something only where the limbs below the window are not all zero.
            const int low = static_cast<int>(std::min<slong>(move, total));
            dropped = move > 0 && (!mpn_zero_p(pr, low) || (!real_only && !mpn_zero_p(pi, low)));
            tr = place(pr, limbs, total, move);
            ti = real_only ? pi : place(pi, limbs, total, move);
        }
        side = 1 - side;
        const bool flip = m.negated && (k & 1);
        bool truncated = add_part(tr, limbs, c.real, c.real.negative != flip, placed, scratch_.data());
        if (!real_only)
            truncated = add_part(ti, limbs, c.imaginary, c.imaginary.negative != flip, placed, scratch_.data()) || truncated;
        if (dropped || truncated)
            error.add(truncation_units * ((dropped ? 1 : 0) + (truncated ? 1 : 0)), placed);
        if (c.radius != 0) error.add(c.radius, c.radius_exponent);
        unit = placed;
        zero = is_zero(tr, limbs) && (real_only || is_zero(ti, limbs));
    }
    const auto set = [&](arf_t out, const ulong *t) {
        const bool negative = t[limbs - 1] >> (limb_bits - 1);
        const ulong *magnitude = t;
        if (negative) {
            mpn_neg(scratch_.data(), t, limbs);
            magnitude = scratch_.data();
        }
        slong count = limbs;
        while (count > 0 && !magnitude[count - 1]) --count;
        if (!count) {
            arf_zero(out);
            return;
        }
        arf_set_mpn(out, magnitude, count, negative);
        arf_mul_2exp_si(out, out, unit);
    };
    Magnitude radius;
    mag_set_d(radius.value(), error.mantissa);
    mag_mul_2exp_si(radius.value(), radius.value(), error.exponent);
    if (zero)
        arf_zero(arb_midref(acb_realref(value)));
    else
        set(arb_midref(acb_realref(value)), tr);
    mag_set(arb_radref(acb_realref(value)), radius.value());
    if (real_only || zero)
        arf_zero(arb_midref(acb_imagref(value)));
    else
        set(arb_midref(acb_imagref(value)), ti);
    if (real_only)
        mag_zero(arb_radref(acb_imagref(value)));
    else
        mag_set(arb_radref(acb_imagref(value)), radius.value());
    return acb_is_finite(value);
}
} // namespace certpolysolve::detail
