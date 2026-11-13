// The shared arithmetic implementations used to evaluate exact program input. Dual balls
// enclose a value and derivative; top coefficients determine degree; modular images prove
// square-freeness; majorants bound coefficients. Each arithmetic owns its working values.
// Called by program preparation and ProgramSource through the generic interpreter.
#include "numeric/evaluate/program_arithmetic.hpp"
#include "numeric/evaluate/program_interpreter.hpp"
#include <flint/nmod.h>
#include <flint/ulong_extras.h>
#include <algorithm>
#include <limits>
#include <stdexcept>

namespace certpolysolve::detail {

// ------------------------------------------------------------------------------------------- dual balls

DualBallArithmetic::DualBallArithmetic(const acb_t point, slong precision) : precision_(precision) {
    acb_set(point_.value(), point);
}
DualBallArithmetic::Value DualBallArithmetic::variable() {
    Value v;
    acb_set(v.value.value(), point_.value());
    acb_one(v.derivative.value());
    return v;
}
DualBallArithmetic::Value DualBallArithmetic::constant(const Coefficient &c) {
    Value v;
    c.enclose(v.value.value(), precision_);
    return v;
}
DualBallArithmetic::Value DualBallArithmetic::add(const Value &a, const Value &b) {
    Value v;
    acb_add(v.value.value(), a.value.value(), b.value.value(), precision_);
    acb_add(v.derivative.value(), a.derivative.value(), b.derivative.value(), precision_);
    return v;
}
DualBallArithmetic::Value DualBallArithmetic::subtract(const Value &a, const Value &b) {
    Value v;
    acb_sub(v.value.value(), a.value.value(), b.value.value(), precision_);
    acb_sub(v.derivative.value(), a.derivative.value(), b.derivative.value(), precision_);
    return v;
}
DualBallArithmetic::Value DualBallArithmetic::multiply(const Value &a, const Value &b) {
    Value v;
    ComplexBall term;
    acb_mul(v.value.value(), a.value.value(), b.value.value(), precision_);
    acb_mul(v.derivative.value(), a.derivative.value(), b.value.value(), precision_);
    acb_mul(term.value(), a.value.value(), b.derivative.value(), precision_);
    acb_add(v.derivative.value(), v.derivative.value(), term.value(), precision_);
    return v;
}
DualBallArithmetic::Value DualBallArithmetic::negate(const Value &a) {
    Value v;
    acb_neg(v.value.value(), a.value.value());
    acb_neg(v.derivative.value(), a.derivative.value());
    return v;
}
DualBallArithmetic::Value DualBallArithmetic::power(const Value &a, uint64_t exponent) {
    // (u^e)' = e u^(e-1) u'.
    Value v;
    ComplexBall lower;
    acb_pow_ui(lower.value(), a.value.value(), exponent - 1, precision_);
    acb_mul(v.value.value(), lower.value(), a.value.value(), precision_);
    acb_mul(v.derivative.value(), lower.value(), a.derivative.value(), precision_);
    acb_mul_ui(v.derivative.value(), v.derivative.value(), exponent, precision_);
    return v;
}

// ------------------------------------------------------------------------------- top coefficients

namespace {
bool gaussian_zero(const GaussianRational &x) { return x.zero(); }
GaussianRational gaussian_add(const GaussianRational &a, const GaussianRational &b) {
    GaussianRational r;
    fmpq_add(r.real.value(), a.real.value(), b.real.value());
    fmpq_add(r.imag.value(), a.imag.value(), b.imag.value());
    return r;
}
GaussianRational gaussian_negate(const GaussianRational &a) {
    GaussianRational r;
    fmpq_neg(r.real.value(), a.real.value());
    fmpq_neg(r.imag.value(), a.imag.value());
    return r;
}
// r += a b
void gaussian_addmul(GaussianRational &r, const GaussianRational &a, const GaussianRational &b) {
    Rational t;
    fmpq_mul(t.value(), a.real.value(), b.real.value());
    fmpq_add(r.real.value(), r.real.value(), t.value());
    fmpq_mul(t.value(), a.imag.value(), b.imag.value());
    fmpq_sub(r.real.value(), r.real.value(), t.value());
    fmpq_mul(t.value(), a.real.value(), b.imag.value());
    fmpq_add(r.imag.value(), r.imag.value(), t.value());
    fmpq_mul(t.value(), a.imag.value(), b.real.value());
    fmpq_add(r.imag.value(), r.imag.value(), t.value());
}
GaussianRational exact_constant(const Coefficient &c, uint64_t &budget) {
    if (!c.exact()) throw std::invalid_argument("program constants must be exact");
    GaussianRational r;
    if (!c.real.rational_with_budget(r.real.value(), budget) || !c.imag.rational_with_budget(r.imag.value(), budget))
        throw std::out_of_range("a program constant exceeds the exact storage budget");
    return r;
}
constexpr uint64_t unbounded = std::numeric_limits<uint64_t>::max();
} // namespace

TopArithmetic::TopArithmetic(size_t width, uint64_t exact_bits) : width_(std::max<size_t>(width, 1)), exact_bits_(exact_bits) {}

TopArithmetic::Value TopArithmetic::variable() {
    Value v;
    v.degree = 1;
    v.top.resize(std::min<size_t>(2, width_));
    fmpq_one(v.top[0].real.value());
    v.complete = width_ >= 2;
    return v;
}
TopArithmetic::Value TopArithmetic::constant(const Coefficient &c) {
    Value v;
    auto budget = exact_bits_;
    auto value = exact_constant(c, budget);
    if (value.zero()) {
        v.zero = v.complete = true;
        return v;
    }
    v.top.push_back(std::move(value));
    v.complete = true;
    return v;
}
TopArithmetic::Value TopArithmetic::add(const Value &a, const Value &b) {
    if (a.undetermined || b.undetermined) return Value{false, true};
    if (a.zero) return b;
    if (b.zero) return a;
    const uint64_t degree = std::max(a.degree, b.degree);
    const uint64_t sa = degree - a.degree, sb = degree - b.degree;
    // How far down from x^degree each operand's coefficients are known: all of them when complete.
    const uint64_t ka = a.complete ? unbounded : sa + a.top.size();
    const uint64_t kb = b.complete ? unbounded : sb + b.top.size();
    const uint64_t full = degree == unbounded ? unbounded : degree + 1;
    // At most the width: the higher operand (shift 0) knows at most `width` coefficients, or is complete
    // and then has at most `width` in all.
    const uint64_t known = std::min({ka, kb, full});
    std::vector<GaussianRational> sum(static_cast<size_t>(known));
    for (uint64_t k = 0; k < known; ++k) {
        if (k >= sa && k - sa < a.top.size()) sum[k] = gaussian_add(sum[k], a.top[k - sa]);
        if (k >= sb && k - sb < b.top.size()) sum[k] = gaussian_add(sum[k], b.top[k - sb]);
    }
    size_t lead = 0;
    while (lead < sum.size() && gaussian_zero(sum[lead])) ++lead;
    Value v;
    if (lead == sum.size()) {
        if (a.complete && b.complete && known == full) {
            v.zero = v.complete = true; // every coefficient cancelled: the zero polynomial
            return v;
        }
        v.undetermined = true;
        return v;
    }
    v.degree = degree - lead;
    v.complete = a.complete && b.complete && known == full;
    for (size_t k = lead; k < sum.size() && v.top.size() < width_; ++k) v.top.push_back(std::move(sum[k]));
    if (v.complete && v.top.size() < sum.size() - lead) v.complete = false; // truncated to the width
    return v;
}
TopArithmetic::Value TopArithmetic::negate(const Value &a) {
    Value v = a;
    for (auto &c : v.top) c = gaussian_negate(c);
    return v;
}
TopArithmetic::Value TopArithmetic::subtract(const Value &a, const Value &b) { return add(a, negate(b)); }
TopArithmetic::Value TopArithmetic::multiply(const Value &a, const Value &b) {
    if (a.undetermined || b.undetermined) return Value{false, true};
    if (a.zero || b.zero) {
        Value v;
        v.zero = v.complete = true;
        return v;
    }
    Value v;
    if (a.degree > unbounded - b.degree) throw std::out_of_range("program degree overflows");
    v.degree = a.degree + b.degree;
    // Coefficient k from the top needs a's and b's first k + 1: known while both are, or past a complete
    // operand's end, where its coefficients are zero.
    const uint64_t ka = a.complete ? unbounded : a.top.size();
    const uint64_t kb = b.complete ? unbounded : b.top.size();
    const uint64_t full = v.degree + 1;
    const uint64_t known = std::min({ka, kb, full, static_cast<uint64_t>(width_)});
    v.top.resize(static_cast<size_t>(known));
    for (uint64_t k = 0; k < known; ++k)
        for (uint64_t i = 0; i <= k; ++i)
            if (i < a.top.size() && k - i < b.top.size()) gaussian_addmul(v.top[k], a.top[i], b.top[k - i]);
    v.complete = a.complete && b.complete && known == full;
    return v;
}
TopArithmetic::Value TopArithmetic::power(const Value &a, uint64_t exponent) {
    Value result, base = a;
    bool first = true;
    for (; exponent; exponent >>= 1) {
        if (exponent & 1) {
            result = first ? base : multiply(result, base);
            first = false;
        }
        if (exponent > 1) base = multiply(base, base);
    }
    return result;
}

ProgramTop top_of_program(const Program &program, size_t max_width, uint64_t exact_bits) {
    ProgramTop out;
    for (size_t width = 1; width <= std::max<size_t>(max_width, 1); width *= 2) {
        TopArithmetic arithmetic(width, exact_bits);
        const auto value = run_program(program, arithmetic);
        out.width = width;
        if (value.zero) {
            out.verdict = ProgramTop::Verdict::ZeroPolynomial;
            return out;
        }
        if (!value.undetermined && !value.top.empty()) {
            out.verdict = ProgramTop::Verdict::Determined;
            out.degree = value.degree;
            out.leading = value.top.front();
            return out;
        }
    }
    out.verdict = ProgramTop::Verdict::Undetermined;
    return out;
}

// ---------------------------------------------------------------------------------------- modular

GaussianPrime gaussian_prime(unsigned index) {
    ulong p = UWORD(1) << 62;
    for (unsigned j = 0; j <= index; ++j) {
        do p = n_nextprime(p, 1);
        while (p % 4 != 1);
    }
    GaussianPrime out;
    out.prime = p;
    out.imaginary_unit = n_sqrtmod(p - 1, p);
    return out;
}

std::optional<ulong> reduce_modulo(const Number &number, ulong prime) {
    nmod_t mod;
    nmod_init(&mod, prime);
    const ulong denominator = fmpz_fdiv_ui(fmpq_denref(number.significand.value()), prime);
    if (!denominator) return std::nullopt;
    const ulong numerator = fmpz_fdiv_ui(fmpq_numref(number.significand.value()), prime);
    const ulong radix = number.radix % prime;
    if (!radix) return std::nullopt;
    // radix^exponent with the exponent reduced modulo p - 1 (Fermat), nonnegative after fdiv.
    const ulong exponent = fmpz_fdiv_ui(number.exponent.value(), prime - 1);
    ulong out = nmod_mul(numerator, nmod_inv(denominator, mod), mod);
    return nmod_mul(out, nmod_pow_ui(radix, exponent, mod), mod);
}

ModularArithmetic::Value ModularArithmetic::variable() {
    Value v(prime_.prime);
    nmod_poly_set_coeff_ui(v.value(), 1, 1);
    return v;
}
ModularArithmetic::Value ModularArithmetic::constant(const Coefficient &c) {
    Value v(prime_.prime);
    const auto real = reduce_modulo(c.real, prime_.prime), imag = reduce_modulo(c.imag, prime_.prime);
    if (!c.exact() || !real || !imag) {
        failed_ = true;
        return v;
    }
    nmod_t mod;
    nmod_init(&mod, prime_.prime);
    nmod_poly_set_coeff_ui(v.value(), 0, nmod_add(*real, nmod_mul(*imag, prime_.imaginary_unit, mod), mod));
    return v;
}
ModularArithmetic::Value ModularArithmetic::add(const Value &a, const Value &b) {
    Value v(prime_.prime);
    nmod_poly_add(v.value(), a.value(), b.value());
    return v;
}
ModularArithmetic::Value ModularArithmetic::subtract(const Value &a, const Value &b) {
    Value v(prime_.prime);
    nmod_poly_sub(v.value(), a.value(), b.value());
    return v;
}
ModularArithmetic::Value ModularArithmetic::multiply(const Value &a, const Value &b) {
    Value v(prime_.prime);
    nmod_poly_mul(v.value(), a.value(), b.value());
    return v;
}
ModularArithmetic::Value ModularArithmetic::negate(const Value &a) {
    Value v(prime_.prime);
    nmod_poly_neg(v.value(), a.value());
    return v;
}
ModularArithmetic::Value ModularArithmetic::power(const Value &a, uint64_t exponent) {
    Value v(prime_.prime);
    nmod_poly_pow(v.value(), a.value(), exponent);
    return v;
}

// --------------------------------------------------------------------------------------- majorant

MajorantArithmetic::MajorantArithmetic(const mag_t radius) { mag_set(radius_.value(), radius); }
MajorantArithmetic::Value MajorantArithmetic::variable() { return radius_; }
MajorantArithmetic::Value MajorantArithmetic::constant(const Coefficient &c) {
    // An upper bound of the modulus from a ball enclosure; 64 bits are enough for a bound.
    constexpr slong bound_bits = 64;
    ComplexBall enclosure;
    c.enclose(enclosure.value(), bound_bits);
    Value v;
    acb_get_mag(v.value(), enclosure.value());
    return v;
}
MajorantArithmetic::Value MajorantArithmetic::add(const Value &a, const Value &b) {
    Value v;
    mag_add(v.value(), a.value(), b.value());
    return v;
}
MajorantArithmetic::Value MajorantArithmetic::subtract(const Value &a, const Value &b) { return add(a, b); }
MajorantArithmetic::Value MajorantArithmetic::multiply(const Value &a, const Value &b) {
    Value v;
    mag_mul(v.value(), a.value(), b.value());
    return v;
}
MajorantArithmetic::Value MajorantArithmetic::negate(const Value &a) { return a; }
MajorantArithmetic::Value MajorantArithmetic::power(const Value &a, uint64_t exponent) {
    Value v;
    mag_pow_ui(v.value(), a.value(), exponent);
    return v;
}

} // namespace certpolysolve::detail
