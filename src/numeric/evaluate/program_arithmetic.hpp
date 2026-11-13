#pragma once
// Numerical interpretations of exact program input: dual balls, exact leading coefficients,
// modular polynomial images and coefficient majorants. No solver state or dispatch.
#include "input/program/program.hpp"
#include "core/raii.hpp"
#include <flint/nmod_poly.h>

namespace certpolysolve::detail {
// Value and derivative in complex balls at a point (forward-mode differentiation).
class DualBallArithmetic {
    ComplexBall point_;
    slong precision_;

  public:
    struct Value {
        ComplexBall value, derivative;
    };
    DualBallArithmetic(const acb_t point, slong precision);
    Value variable();
    Value constant(const Coefficient &c);
    Value add(const Value &a, const Value &b);
    Value subtract(const Value &a, const Value &b);
    Value multiply(const Value &a, const Value &b);
    Value negate(const Value &a);
    Value power(const Value &a, uint64_t exponent);
};

// The top coefficients of each value, exactly, over Q(i): the exact degree and leading coefficient of the
// program's polynomial without its expansion. A value keeps its degree and its highest `known`
// coefficients (at most `width`), and whether those are all of them (`complete`). A sum whose leading
// coefficients cancel loses that many known coefficients; when none is left the value is undetermined,
// and a wider run decides it (top_of_program).
class TopArithmetic {
    size_t width_;
    uint64_t exact_bits_;

  public:
    struct Value {
        bool zero = false, undetermined = false, complete = false;
        uint64_t degree = 0;
        std::vector<GaussianRational> top = {}; // top[0] is the coefficient of x^degree
    };
    // `exact_bits`: the storage a constant may take when made rational (Number::rational_with_budget).
    TopArithmetic(size_t width, uint64_t exact_bits);
    Value variable();
    Value constant(const Coefficient &c);
    Value add(const Value &a, const Value &b);
    Value subtract(const Value &a, const Value &b);
    Value multiply(const Value &a, const Value &b);
    Value negate(const Value &a);
    Value power(const Value &a, uint64_t exponent);
};

struct ProgramTop {
    enum class Verdict { Determined, ZeroPolynomial, Undetermined };
    Verdict verdict = Verdict::Undetermined;
    uint64_t degree = 0;
    GaussianRational leading;
    size_t width = 0; // the width that decided it
};
// Runs TopArithmetic with width 1, 2, 4, ... up to `max_width` until the output's degree is decided.
ProgramTop top_of_program(const Program &program, size_t max_width, uint64_t exact_bits);

// An owning nmod_poly.
class ModularPolynomial {
    nmod_poly_t value_;

  public:
    explicit ModularPolynomial(ulong modulus) { nmod_poly_init(value_, modulus); }
    ~ModularPolynomial() { nmod_poly_clear(value_); }
    ModularPolynomial(const ModularPolynomial &x) : ModularPolynomial(x.value_->mod.n) { nmod_poly_set(value_, x.value_); }
    ModularPolynomial(ModularPolynomial &&x) noexcept : ModularPolynomial(x.value_->mod.n) { nmod_poly_swap(value_, x.value_); }
    ModularPolynomial &operator=(ModularPolynomial x) noexcept {
        nmod_poly_swap(value_, x.value_);
        return *this;
    }
    nmod_poly_struct *value() & noexcept { return value_; }
    const nmod_poly_struct *value() const & noexcept { return value_; }
    nmod_poly_struct *value() && = delete;
    const nmod_poly_struct *value() const && = delete;
};

// A prime p = 1 (mod 4) just above 2^62, so of 63 bits, and a square root of -1 modulo p: the `index`-th such
// prime from 2^62 upward (a fixed sequence, so every run reduces the same way).
struct GaussianPrime {
    ulong prime = 0, imaginary_unit = 0;
};
GaussianPrime gaussian_prime(unsigned index);
// The residue of an exact number modulo p; nullopt when p divides its denominator.
std::optional<ulong> reduce_modulo(const Number &number, ulong prime);

// The image modulo p, with i mapped to imaginary_unit. A constant whose denominator p divides makes the
// arithmetic fail (failed()); the caller then takes the next prime.
class ModularArithmetic {
    GaussianPrime prime_;
    bool failed_ = false;

  public:
    using Value = ModularPolynomial;
    explicit ModularArithmetic(GaussianPrime prime) : prime_(prime) {}
    bool failed() const { return failed_; }
    Value variable();
    Value constant(const Coefficient &c);
    Value add(const Value &a, const Value &b);
    Value subtract(const Value &a, const Value &b);
    Value multiply(const Value &a, const Value &b);
    Value negate(const Value &a);
    Value power(const Value &a, uint64_t exponent);
};

// A majorant: with every constant replaced by an upper bound of its modulus and subtraction by addition,
// evaluated at a real r >= 0, the program bounds the sum of |a_i| r^i of every node's coefficients.
class MajorantArithmetic {
    Magnitude radius_;

  public:
    using Value = Magnitude;
    explicit MajorantArithmetic(const mag_t radius);
    Value variable();
    Value constant(const Coefficient &c);
    Value add(const Value &a, const Value &b);
    Value subtract(const Value &a, const Value &b);
    Value multiply(const Value &a, const Value &b);
    Value negate(const Value &a);
    Value power(const Value &a, uint64_t exponent);
};

} // namespace certpolysolve::detail
