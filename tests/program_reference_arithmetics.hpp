#pragma once
// Two arithmetics of the program interpreter (run_program, run_prepared: program.hpp) that only the tests and a
// benchmark use, test support since W8's S2 (23 October 2026; until then in src/program.{hpp,cpp}, which no library
// code called): ComplexDoubleArithmetic, value and derivative in binary64, unchecked (tests/test_program.cpp and
// benchmarks/program/eval_overhead.cpp, which measures the interpreter's overhead against the batched evaluator of
// program_machine.hpp); ExpansionArithmetic, the exact expansion over Q(i) (tests/test_program.cpp's exact oracle: the
// expansion the program input exists to avoid). Moved unchanged but for inline; ExpansionArithmetic::exact_constant is
// a copy of program.cpp's (an anonymous helper the maintained top-coefficient arithmetic keeps).
#include "input/program/program.hpp"
#include <flint/fmpq_poly.h>
#include <algorithm>
#include <complex>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace certpolysolve::detail {

// Value and derivative in binary64 through the generic interpreter, unchecked: for measuring the
// interpreter's overhead against the batched evaluator (program_machine.hpp), and for tests. The program's
// constants are rounded once, at construction, and taken by index (constant_at).
class ComplexDoubleArithmetic {
    std::complex<double> point_;
    std::vector<std::complex<double>> constants_;

  public:
    struct Value {
        std::complex<double> value, derivative;
    };
    ComplexDoubleArithmetic(const Program &program, std::complex<double> point);
    void set_point(std::complex<double> point) { point_ = point; }
    Value variable() const { return {point_, 1.0}; }
    Value constant_at(uint64_t index) const { return {constants_[static_cast<size_t>(index)], 0.0}; }
    Value add(const Value &a, const Value &b) const { return {a.value + b.value, a.derivative + b.derivative}; }
    Value subtract(const Value &a, const Value &b) const { return {a.value - b.value, a.derivative - b.derivative}; }
    Value multiply(const Value &a, const Value &b) const {
        return {a.value * b.value, a.derivative * b.value + a.value * b.derivative};
    }
    Value negate(const Value &a) const { return {-a.value, -a.derivative}; }
    Value power(const Value &a, uint64_t exponent) const {
        // a^(e-1) by binary powering, no squaring past the last exponent bit (x^3 at |x| = 2^400 would
        // otherwise form x^4 and overflow though the result fits).
        std::complex<double> lower = 1.0, base = a.value;
        for (uint64_t e = exponent - 1; e;) {
            if (e & 1) lower *= base;
            e >>= 1;
            if (e) base *= base;
        }
        return {lower * a.value, static_cast<double>(exponent) * lower * a.derivative};
    }
};
inline ComplexDoubleArithmetic::ComplexDoubleArithmetic(const Program &program, std::complex<double> point)
    : point_(point) {
    constexpr slong rounding_bits = 64;
    ComplexBall c;
    for (const auto &constant : program.constants) {
        constant.enclose(c.value(), rounding_bits);
        constants_.emplace_back(arf_get_d(arb_midref(acb_realref(c.value())), ARF_RND_NEAR),
                                arf_get_d(arb_midref(acb_imagref(c.value())), ARF_RND_NEAR));
    }
}

// The exact expansion over Q(i), as real and imaginary parts (tests only: it is the expansion the program
// input exists to avoid).
class GaussianExpansion {
    fmpq_poly_t real_, imag_;

  public:
    GaussianExpansion() { fmpq_poly_init(real_); fmpq_poly_init(imag_); }
    ~GaussianExpansion() { fmpq_poly_clear(real_); fmpq_poly_clear(imag_); }
    GaussianExpansion(const GaussianExpansion &x) : GaussianExpansion() {
        fmpq_poly_set(real_, x.real_);
        fmpq_poly_set(imag_, x.imag_);
    }
    GaussianExpansion(GaussianExpansion &&x) noexcept : GaussianExpansion() {
        fmpq_poly_swap(real_, x.real_);
        fmpq_poly_swap(imag_, x.imag_);
    }
    GaussianExpansion &operator=(GaussianExpansion x) noexcept {
        fmpq_poly_swap(real_, x.real_);
        fmpq_poly_swap(imag_, x.imag_);
        return *this;
    }
    fmpq_poly_struct *real() & noexcept { return real_; }
    const fmpq_poly_struct *real() const & noexcept { return real_; }
    fmpq_poly_struct *imag() & noexcept { return imag_; }
    const fmpq_poly_struct *imag() const & noexcept { return imag_; }
    slong degree() const { return std::max(fmpq_poly_degree(real_), fmpq_poly_degree(imag_)); }
};
class ExpansionArithmetic {
    uint64_t exact_bits_;
    static GaussianRational exact_constant(const Coefficient &c, uint64_t &budget) {
        if (!c.exact()) throw std::invalid_argument("program constants must be exact");
        GaussianRational r;
        if (!c.real.rational_with_budget(r.real.value(), budget) ||
            !c.imag.rational_with_budget(r.imag.value(), budget))
            throw std::out_of_range("a program constant exceeds the exact storage budget");
        return r;
    }

  public:
    using Value = GaussianExpansion;
    explicit ExpansionArithmetic(uint64_t exact_bits) : exact_bits_(exact_bits) {}
    Value variable();
    Value constant(const Coefficient &c);
    Value add(const Value &a, const Value &b);
    Value subtract(const Value &a, const Value &b);
    Value multiply(const Value &a, const Value &b);
    Value negate(const Value &a);
    Value power(const Value &a, uint64_t exponent);
};
inline ExpansionArithmetic::Value ExpansionArithmetic::variable() {
    Value v;
    fmpq_poly_set_coeff_si(v.real(), 1, 1);
    return v;
}
inline ExpansionArithmetic::Value ExpansionArithmetic::constant(const Coefficient &c) {
    auto budget = exact_bits_;
    const auto value = exact_constant(c, budget);
    Value v;
    fmpq_poly_set_coeff_fmpq(v.real(), 0, value.real.value());
    fmpq_poly_set_coeff_fmpq(v.imag(), 0, value.imag.value());
    return v;
}
inline ExpansionArithmetic::Value ExpansionArithmetic::add(const Value &a, const Value &b) {
    Value v;
    fmpq_poly_add(v.real(), a.real(), b.real());
    fmpq_poly_add(v.imag(), a.imag(), b.imag());
    return v;
}
inline ExpansionArithmetic::Value ExpansionArithmetic::subtract(const Value &a, const Value &b) {
    Value v;
    fmpq_poly_sub(v.real(), a.real(), b.real());
    fmpq_poly_sub(v.imag(), a.imag(), b.imag());
    return v;
}
inline ExpansionArithmetic::Value ExpansionArithmetic::multiply(const Value &a, const Value &b) {
    // (ar + i ai)(br + i bi) = (ar br - ai bi) + i (ar bi + ai br)
    Value v;
    fmpq_poly_t t;
    fmpq_poly_init(t);
    fmpq_poly_mul(v.real(), a.real(), b.real());
    fmpq_poly_mul(t, a.imag(), b.imag());
    fmpq_poly_sub(v.real(), v.real(), t);
    fmpq_poly_mul(v.imag(), a.real(), b.imag());
    fmpq_poly_mul(t, a.imag(), b.real());
    fmpq_poly_add(v.imag(), v.imag(), t);
    fmpq_poly_clear(t);
    return v;
}
inline ExpansionArithmetic::Value ExpansionArithmetic::negate(const Value &a) {
    Value v;
    fmpq_poly_neg(v.real(), a.real());
    fmpq_poly_neg(v.imag(), a.imag());
    return v;
}
inline ExpansionArithmetic::Value ExpansionArithmetic::power(const Value &a, uint64_t exponent) {
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

} // namespace certpolysolve::detail
