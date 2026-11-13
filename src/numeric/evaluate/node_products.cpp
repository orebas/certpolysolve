// ShortNodeProducts (node_products.hpp derives its bound): a secular model's node products prod_{j != i} (b_i - b_j)
// (product), and products over a subset of this or another node set (product_over), for exact nodes with short parts
// (binary64 or at most two limbs), in integer mantissas with a relative error bound known in advance (accuracy_bits).
//
// The constructor reads each node's parts as LongPart and, where all have at most 53 bits, as Part (decompose). Each
// product runs a Chain from one. Chain::multiply forms the factor b_i - b_j part by part (difference: a signed 128-bit
// integer from binary64 parts, a signed 256-bit Product from long ones) and floors both to one unit, the larger with
// factor_bits bits (scaled). Chain::step multiplies the held product, two signed 128-bit integers with one exponent,
// by the factor exactly in 256 bits (wide_product, the + and - of Product) and floors it back to product_bits bits
// (top_bit_length, shifted): a relative error below 2^step_error_log2 a factor. Chain::enclose gives the held product
// the radius m 2^(step_error_log2 + 1) times its modulus in both coordinates after m factors, refusing a ball that
// holds zero. Chain::least_modulus bounds the least |b_i - b_j| from below by a binary64 estimate of correctly rounded
// products, a sum and a square root; all else is exact or an integer floor, so the results are the same on every CPU.
//
// State: the nodes' parts, read once (long_real_, long_imaginary_; real_, imaginary_ where binary64_); product and
// product_over are const, their running product a local Chain.
// Entry points: the constructor, usable, accuracy_bits, product and product_over, from regenerate_residues
// (secular_refinement.cpp) under RegenerationPolicy::integer_node_products and reuse_products (session.cpp's
// maintained_method_policy sets both for the secular method): product for a residue's denominator, product_over for
// the update of an unmoved node's residue when other nodes move.
// Not here: the ball products used otherwise (product_derivative, secular_refinement.cpp).
#include "numeric/evaluate/node_products.hpp"
#include <algorithm>
#include <bit>
#include <cmath>
#include <limits>

namespace certpolysolve::detail {
namespace {
using Limb = ulong;
using Wide = uint128;
using Signed = int128;
constexpr int limb_bits = 64;
static_assert(FLINT_BITS == limb_bits, "the node products are written for 64-bit limbs");

using Part = ShortNodeProducts::Part;
Part
decompose(double x) {
    if (x == 0) return {};
    int exponent = 0;
    const double fraction = std::frexp(x, &exponent); // exact: |fraction| in [1/2, 1)
    return {static_cast<int64_t>(std::ldexp(fraction, std::numeric_limits<double>::digits)),
            exponent - std::numeric_limits<double>::digits};
}
int
bit_length(Wide magnitude) {
    const auto high = static_cast<Limb>(magnitude >> limb_bits), low = static_cast<Limb>(magnitude);
    return high ? 2 * limb_bits - std::countl_zero(high) : limb_bits - std::countl_zero(low);
}
Wide
magnitude(Signed x) {
    return x < 0 ? Wide{0} - static_cast<Wide>(x) : static_cast<Wide>(x);
}
// floor(x 2^shift), for a shift up that keeps x within 127 bits or any shift down: within one unit of x 2^shift.
// (A shift of a negative two's complement integer is a multiplication or a floor in C++20.)
inline Signed
scaled(Signed x, slong shift) {
    constexpr slong all = 2 * limb_bits - 1;
    return shift >= 0 ? x << shift : x >> std::min(-shift, all);
}
// a - b for two parts, as D 2^q0 (the header's "a factor"): exact, or within one unit 2^q0 of the exact
// difference when the exponents are more than difference_shift apart.
struct Difference {
    Signed value = 0;
    slong exponent = 0;
};
inline Difference
difference(Part a, Part b) {
    if (!a.mantissa) return {-static_cast<Signed>(b.mantissa), b.exponent};
    if (!b.mantissa) return {a.mantissa, a.exponent};
    const slong low = std::min(a.exponent, b.exponent), high = std::max(a.exponent, b.exponent);
    const slong base = std::max(low, high - ShortNodeProducts::difference_shift);
    return {scaled(a.mantissa, a.exponent - base) - scaled(b.mantissa, b.exponent - base), base};
}

// A signed 256-bit integer in two's complement as two halves, low first.
struct Product {
    Wide low = 0, high = 0;
};
// The exact product of two signed 128-bit integers: the unsigned product of the two's complement readings, then
// the correction of its high half (a negative operand is its unsigned reading minus 2^128, so the other operand
// leaves the high half once for each negative one, modulo 2^256). Held in registers throughout.
inline Product
wide_product(Signed a, Signed b) {
    const Wide x = static_cast<Wide>(a), y = static_cast<Wide>(b);
    const Limb x0 = static_cast<Limb>(x), x1 = static_cast<Limb>(x >> limb_bits);
    const Limb y0 = static_cast<Limb>(y), y1 = static_cast<Limb>(y >> limb_bits);
    const Wide p00 = static_cast<Wide>(x0) * y0, p01 = static_cast<Wide>(x0) * y1;
    const Wide p10 = static_cast<Wide>(x1) * y0, p11 = static_cast<Wide>(x1) * y1;
    // At most 3 (2^64 - 1): no carry out of the middle sum.
    const Wide middle = (p00 >> limb_bits) + static_cast<Limb>(p01) + static_cast<Limb>(p10);
    Product r;
    r.low = (middle << limb_bits) | static_cast<Limb>(p00);
    r.high = p11 + (p01 >> limb_bits) + (p10 >> limb_bits) + (middle >> limb_bits);
    r.high -= (a < 0 ? y : Wide{0}) + (b < 0 ? x : Wide{0});
    return r;
}
inline Product
operator+(Product a, Product b) {
    Product r;
    r.low = a.low + b.low;
    r.high = a.high + b.high + (r.low < a.low ? 1 : 0);
    return r;
}
inline Product
operator-(Product a, Product b) {
    Product r;
    r.low = a.low - b.low;
    r.high = a.high - b.high - (a.low < b.low ? 1 : 0);
    return r;
}
// The bit length of the magnitude of a product whose magnitude reaches its high half's upper limb (zero when it
// does not): that limb of the magnitude is the complement of the limb, and one more when all below it is zero.
inline int
top_bit_length(Product x) {
    const Limb top = static_cast<Limb>(x.high >> limb_bits);
    const Limb negative = Limb{0} - (top >> (limb_bits - 1));
    const bool below_zero = x.low == 0 && static_cast<Limb>(x.high) == 0;
    const Limb magnitude_top = (top ^ negative) + (negative & static_cast<Limb>(below_zero));
    return magnitude_top ? 4 * limb_bits - std::countl_zero(magnitude_top) : 0;
}
// A long part as a signed 256-bit integer, and the helpers of a difference of two.
using LongPart = ShortNodeProducts::LongPart;
inline Product
negate(Product x) {
    return {Wide{0} - x.low, ~x.high + (x.low == 0 ? Wide{1} : Wide{0})};
}
inline bool
negative(Product x) {
    return static_cast<Signed>(x.high) < 0;
}
inline Product
signed_wide(LongPart part, bool negated = false) {
    const Product size{part.magnitude, 0};
    return part.negative != negated ? negate(size) : size;
}
inline Product
shift_left(Product x, slong shift) { // 0 <= shift < 128, the result within 255 bits
    if (!shift) return x;
    return {x.low << shift, (x.high << shift) | (x.low >> (2 * limb_bits - shift))};
}
inline int
bit_length(Product size) { // of a magnitude
    return size.high ? 2 * limb_bits + bit_length(size.high) : bit_length(size.low);
}
// floor(x 2^shift) as a signed 128-bit integer, for a shift up that keeps x within 127 bits or any shift down.
inline Signed
scaled(Product x, slong shift) {
    constexpr slong half = 2 * limb_bits;
    if (shift >= 0) return static_cast<Signed>(x.low << shift);
    const slong down = -shift;
    if (down >= 2 * half) return negative(x) ? Signed{-1} : Signed{0};
    if (down >= half) return static_cast<Signed>(x.high) >> (down - half);
    return static_cast<Signed>((x.low >> down) | (x.high << (half - down)));
}
// a - b for two long parts, as D 2^q0 (the header's "a factor between longer parts"): exact, or the larger
// part alone where the other is below a quarter of the held factor's unit.
struct LongDifference {
    Product value;
    slong exponent = 0;
};
inline LongDifference
difference(LongPart a, LongPart b) {
    if (!a.magnitude) return {signed_wide(b, true), b.exponent};
    if (!b.magnitude) return {signed_wide(a), a.exponent};
    const slong gap = a.exponent - b.exponent;
    if (gap > ShortNodeProducts::long_difference_shift) return {signed_wide(a), a.exponent};
    if (gap < -ShortNodeProducts::long_difference_shift) return {signed_wide(b, true), b.exponent};
    if (gap >= 0) return {shift_left(signed_wide(a), gap) - signed_wide(b), b.exponent};
    return {signed_wide(a) - shift_left(signed_wide(b), -gap), a.exponent};
}
// floor(x / 2^shift) for a shift in (64, 128): the bits of x from `shift` up, 128 of them, taken limb by limb
// (the amounts below are then in (0, 64), and a shift by a whole limb never arises).
inline Signed
shifted(Product x, int shift) {
    const int down = shift - limb_bits, up = 2 * limb_bits - shift;
    const Limb l1 = static_cast<Limb>(x.low >> limb_bits), l2 = static_cast<Limb>(x.high),
               l3 = static_cast<Limb>(x.high >> limb_bits);
    const Limb low = (l1 >> down) | (l2 << up), high = (l2 >> down) | (l3 << up);
    return static_cast<Signed>((static_cast<Wide>(high) << limb_bits) | low);
}
} // namespace

ShortNodeProducts::ShortNodeProducts(std::span<const ComplexBall> nodes) {
    if (nodes.size() < 2) return;
    long_real_.reserve(nodes.size());
    long_imaginary_.reserve(nodes.size());
    bool doubles = true;
    // A part of at most two limbs: its limbs (the last one's leading bit set) times 2^(size - 64 count).
    const auto part = [&](const arf_t x, LongPart &out) {
        out = {};
        if (arf_is_zero(x)) return true;
        constexpr slong two_limbs = 2 * limb_bits;
        if (!arf_is_finite(x) || arf_bits(x) > two_limbs) return false;
        const slong size = arf_abs_bound_lt_2exp_si(x);
        if (size > exponent_range || size < -exponent_range) return false;
        mp_srcptr limbs;
        mp_size_t count;
        ARF_GET_MPN_READONLY(limbs, count, x);
        out.magnitude = count == 2 ? (Wide{limbs[1]} << limb_bits) | limbs[0] : Wide{limbs[0]} << limb_bits;
        out.exponent = size - two_limbs;
        out.negative = ARF_SGNBIT(x);
        doubles = doubles && arf_bits(x) <= std::numeric_limits<double>::digits;
        return true;
    };
    for (const auto &node : nodes) {
        LongPart re, im;
        if (!acb_is_exact(node.value()) || !part(arb_midref(acb_realref(node.value())), re) ||
            !part(arb_midref(acb_imagref(node.value())), im))
            return;
        long_real_.push_back(re);
        long_imaginary_.push_back(im);
    }
    if (doubles) {
        real_.reserve(nodes.size());
        imaginary_.reserve(nodes.size());
        for (const auto &node : nodes) {
            // Exact: at most 53 bits, within the normal range.
            real_.push_back(decompose(arf_get_d(arb_midref(acb_realref(node.value())), ARF_RND_DOWN)));
            imaginary_.push_back(decompose(arf_get_d(arb_midref(acb_imagref(node.value())), ARF_RND_DOWN)));
        }
    }
    binary64_ = doubles;
    usable_ = true;
}
slong
ShortNodeProducts::accuracy_bits(size_t nodes) {
    if (nodes < 2) return -(step_error_log2 + 1);
    return -(step_error_log2 + 1) - static_cast<slong>(std::bit_width(nodes - 2));
}
namespace {
// The held product of a chain of factors (the header's "the product"): (zr + i zi) 2^exponent, the larger part
// in [2^(product_bits - 1), 2^product_bits); and the least modulus among its factors.
class Chain {
    static constexpr int factor_bits = ShortNodeProducts::factor_bits, product_bits = ShortNodeProducts::product_bits;
    Signed zr_ = Signed{1} << (product_bits - 1), zi_ = 0;
    slong exponent_ = -(product_bits - 1);
    size_t factors_ = 0;
    // The least |factor| so far, as nearest_square 2^(2 nearest_exponent) in units of 2^128 squared: an
    // estimate within 2^-50 of the factor's modulus squared.
    double nearest_square_ = 0;
    slong nearest_exponent_ = 0;

    // Times the held factor (dr + i di) 2^unit.
    bool step(Signed dr, Signed di, slong unit) {
        {
            // Its modulus squared from the parts' high limbs: each within 2^64 of its part where the larger part
            // is at least 2^124.
            const double high_real = static_cast<double>(static_cast<int64_t>(dr >> limb_bits)),
                         high_imaginary = static_cast<double>(static_cast<int64_t>(di >> limb_bits));
            const double square = high_real * high_real + high_imaginary * high_imaginary;
            // Squares of held factors differ by less than 2^3 at one unit, so two units of exponent decide.
            constexpr double powers_of_four[] = {0x1p-4, 0x1p-2, 1, 0x1p2, 0x1p4};
            const slong gap = unit - nearest_exponent_;
            if (!factors_ || gap < -2 || (gap <= 2 && square * powers_of_four[gap + 2] < nearest_square_)) {
                nearest_square_ = square;
                nearest_exponent_ = unit;
            }
        }
        const Product next_real = wide_product(zr_, dr) - wide_product(zi_, di);
        const Product next_imaginary = wide_product(zr_, di) + wide_product(zi_, dr);
        const int length = std::max(top_bit_length(next_real), top_bit_length(next_imaginary));
        const int shift = length - product_bits;
        // |held product| >= 2^(product_bits - 1) and |held factor| >= 2^(factor_bits - 1) put the larger part of
        // their product between these lengths; anything else would be a defect of this file.
        if (length < product_bits + factor_bits - 2 || length > product_bits + factor_bits + 2) return false;
        zr_ = shifted(next_real, shift);
        zi_ = shifted(next_imaginary, shift);
        exponent_ += unit + shift;
        ++factors_;
        return true;
    }

  public:
    // Times (own - other), both parts binary64. False for a zero factor (the two points coincide).
    bool multiply(Part own_real, Part own_imaginary, Part other_real, Part other_imaginary) {
        const Difference real = difference(own_real, other_real);
        const Difference imaginary = difference(own_imaginary, other_imaginary);
        const int real_length = bit_length(magnitude(real.value)),
                  imaginary_length = bit_length(magnitude(imaginary.value));
        if (!real_length && !imaginary_length) return false;
        const slong real_top = real_length ? real.exponent + real_length : std::numeric_limits<slong>::min();
        const slong imaginary_top =
            imaginary_length ? imaginary.exponent + imaginary_length : std::numeric_limits<slong>::min();
        const slong unit = std::max(real_top, imaginary_top) - factor_bits;
        return step(real_length ? scaled(real.value, real.exponent - unit) : 0,
                    imaginary_length ? scaled(imaginary.value, imaginary.exponent - unit) : 0, unit);
    }
    // The same between parts of at most two limbs.
    bool multiply(LongPart own_real, LongPart own_imaginary, LongPart other_real, LongPart other_imaginary) {
        const LongDifference real = difference(own_real, other_real);
        const LongDifference imaginary = difference(own_imaginary, other_imaginary);
        const int real_length = bit_length(negative(real.value) ? negate(real.value) : real.value),
                  imaginary_length = bit_length(negative(imaginary.value) ? negate(imaginary.value) : imaginary.value);
        if (!real_length && !imaginary_length) return false;
        const slong real_top = real_length ? real.exponent + real_length : std::numeric_limits<slong>::min();
        const slong imaginary_top =
            imaginary_length ? imaginary.exponent + imaginary_length : std::numeric_limits<slong>::min();
        const slong unit = std::max(real_top, imaginary_top) - factor_bits;
        return step(real_length ? scaled(real.value, real.exponent - unit) : 0,
                    imaginary_length ? scaled(imaginary.value, imaginary.exponent - unit) : 0, unit);
    }
    // The product as a ball: the held value with the header's radius for the factors taken, in both
    // coordinates. False if the ball holds zero.
    bool enclose(acb_t product) const {
        Integer integer;
        const auto high = [](Signed x) { return static_cast<Limb>(static_cast<Wide>(x) >> limb_bits); };
        fmpz_set_signed_uiui(integer.value(), high(zr_), static_cast<Limb>(static_cast<Wide>(zr_)));
        arb_set_fmpz(acb_realref(product), integer.value());
        fmpz_set_signed_uiui(integer.value(), high(zi_), static_cast<Limb>(static_cast<Wide>(zi_)));
        arb_set_fmpz(acb_imagref(product), integer.value());
        acb_mul_2exp_si(product, product, exponent_);
        // The exact product is within m 2^(step_error_log2 + 1) of the held one for m factors, in units of its
        // modulus.
        Magnitude radius;
        acb_get_mag(radius.value(), product);
        mag_mul_ui(radius.value(), radius.value(), static_cast<ulong>(factors_));
        mag_mul_2exp_si(radius.value(), radius.value(), ShortNodeProducts::step_error_log2 + 1);
        arb_add_error_mag(acb_realref(product), radius.value());
        arb_add_error_mag(acb_imagref(product), radius.value());
        return !acb_contains_zero(product);
    }
    // A lower bound of the least modulus among the factors: the estimate less 2^-40 of itself, far more than
    // its own error and the factor's.
    void least_modulus(mag_t separation) const {
        constexpr double separation_slack = 0x1p-40;
        mag_set_d_lower(separation, std::sqrt(nearest_square_) * (1 - separation_slack));
        mag_mul_2exp_si(separation, separation, nearest_exponent_ + limb_bits);
    }
};
} // namespace
bool
ShortNodeProducts::product(acb_t product, mag_t separation, size_t i) const {
    if (!usable_ || i >= long_real_.size()) return false;
    Chain chain;
    if (binary64_) {
        for (size_t j = 0; j < real_.size(); ++j)
            if (j != i && !chain.multiply(real_[i], imaginary_[i], real_[j], imaginary_[j])) return false; // two nodes coincide
    } else
        for (size_t j = 0; j < long_real_.size(); ++j)
            if (j != i && !chain.multiply(long_real_[i], long_imaginary_[i], long_real_[j], long_imaginary_[j])) return false;
    chain.least_modulus(separation);
    return chain.enclose(product);
}
bool
ShortNodeProducts::product_over(acb_t product, size_t i, std::span<const size_t> subset,
                                 const ShortNodeProducts &others) const {
    if (!usable_ || !others.usable_ || i >= long_real_.size()) return false;
    Chain chain;
    const bool doubles = binary64_ && others.binary64_;
    for (const size_t j : subset) {
        if (j >= others.long_real_.size()) return false;
        if (doubles ? !chain.multiply(real_[i], imaginary_[i], others.real_[j], others.imaginary_[j])
                    : !chain.multiply(long_real_[i], long_imaginary_[i], others.long_real_[j], others.long_imaginary_[j]))
            return false;
    }
    return chain.enclose(product);
}
} // namespace certpolysolve::detail
