#pragma once

// Owning numerical values. Backend handles are explicit, const-correct borrows.
// A temporary cannot lend storage; no global precision or thread mutation.
#ifdef CERTPOLYSOLVE_EXPERIMENTAL_ARITHMETIC
#include <mpfr.h> // FLINT exposes MPFR conversions only when this precedes arf.h.
#endif
#include <flint/acb_poly.h>
#include <flint/arb_poly.h>
#include <flint/fmpq.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_poly_factor.h>
#include <memory>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

#if !FLINT_USES_TLS
#error "CertPolySolve currently requires a FLINT build with thread-local storage enabled"
#endif
// Objects may be destroyed on another thread than the one that made them, also after its cleanup
// (docs/api-spec-0.5.md section 8, T3): FLINT's memory manager supports that only when built with thread support.
#if !FLINT_USES_PTHREAD
#error "CertPolySolve requires a FLINT build with thread support (FLINT_USES_PTHREAD)"
#endif

namespace certpolysolve::detail {

class Integer {
    fmpz_t value_;

  public:
    Integer() { fmpz_init(value_); }
    ~Integer() { fmpz_clear(value_); }
    Integer(const Integer &x) : Integer() { fmpz_set(value_, x.value_); }
    Integer(Integer &&x) noexcept : Integer() { fmpz_swap(value_, x.value_); }
    Integer &operator=(Integer x) noexcept {
        fmpz_swap(value_, x.value_);
        return *this;
    }
    fmpz * value() & noexcept { return value_; }
    const fmpz * value() const & noexcept { return value_; }
    fmpz * value() && = delete;
    const fmpz * value() const && = delete;
};
class Rational {
    fmpq_t value_;

  public:
    Rational() { fmpq_init(value_); }
    ~Rational() { fmpq_clear(value_); }
    Rational(const Rational &x) : Rational() { fmpq_set(value_, x.value_); }
    Rational(Rational &&x) noexcept : Rational() { fmpq_swap(value_, x.value_); }
    Rational &operator=(Rational x) noexcept {
        fmpq_swap(value_, x.value_);
        return *this;
    }
    fmpq * value() & noexcept { return value_; }
    const fmpq * value() const & noexcept { return value_; }
    fmpq * value() && = delete;
    const fmpq * value() const && = delete;
};
class Magnitude {
    mag_t value_;

  public:
    Magnitude() { mag_init(value_); }
    ~Magnitude() { mag_clear(value_); }
    Magnitude(const Magnitude &x) : Magnitude() { mag_set(value_, x.value_); }
    Magnitude(Magnitude &&x) noexcept : Magnitude() { mag_swap(value_, x.value_); }
    Magnitude &operator=(Magnitude x) noexcept {
        mag_swap(value_, x.value_);
        return *this;
    }
    mag_ptr value() & noexcept { return value_; }
    mag_srcptr value() const & noexcept { return value_; }
    mag_ptr value() && = delete;
    mag_srcptr value() const && = delete;
};
class RealBall {
    arb_t value_;

  public:
    RealBall() { arb_init(value_); }
    ~RealBall() { arb_clear(value_); }
    RealBall(const RealBall &x) : RealBall() { arb_set(value_, x.value_); }
    RealBall(RealBall &&x) noexcept : RealBall() { arb_swap(value_, x.value_); }
    RealBall &operator=(RealBall x) noexcept {
        arb_swap(value_, x.value_);
        return *this;
    }
    arb_ptr value() & noexcept { return value_; }
    arb_srcptr value() const & noexcept { return value_; }
    arb_ptr value() && = delete;
    arb_srcptr value() const && = delete;
};
class Float {
    arf_t value_;

  public:
    Float() { arf_init(value_); }
    ~Float() { arf_clear(value_); }
    Float(const Float &) = delete;
    Float &operator=(const Float &) = delete;
    arf_ptr value() & noexcept { return value_; }
    arf_srcptr value() const & noexcept { return value_; }
    arf_ptr value() && = delete;
    arf_srcptr value() const && = delete;
};
class ComplexBall {
    acb_t value_;

  public:
    ComplexBall() { acb_init(value_); }
    ~ComplexBall() { acb_clear(value_); }
    ComplexBall(const ComplexBall &x) : ComplexBall() { acb_set(value_, x.value_); }
    ComplexBall(ComplexBall &&x) noexcept : ComplexBall() { acb_swap(value_, x.value_); }
    ComplexBall &operator=(ComplexBall x) noexcept {
        acb_swap(value_, x.value_);
        return *this;
    }
    acb_ptr value() & noexcept { return value_; }
    acb_srcptr value() const & noexcept { return value_; }
    acb_ptr value() && = delete;
    acb_srcptr value() const && = delete;
};
// Move-only: polynomial storage is never copied implicitly. Move assignment swaps
// contents, leaving the source with the destination's old polynomial. A value()
// handle refers to this owner; reacquire it from the destination after a move.
class RealPolynomial {
    arb_poly_t value_;

  public:
    RealPolynomial() { arb_poly_init(value_); }
    ~RealPolynomial() { arb_poly_clear(value_); }
    RealPolynomial(const RealPolynomial &) = delete;
    RealPolynomial &operator=(const RealPolynomial &) = delete;
    RealPolynomial(RealPolynomial &&x) noexcept : RealPolynomial() { arb_poly_swap(value_, x.value_); }
    RealPolynomial &operator=(RealPolynomial &&x) noexcept {
        arb_poly_swap(value_, x.value_);
        return *this;
    }
    arb_poly_struct *value() & noexcept { return value_; }
    const arb_poly_struct *value() const & noexcept { return value_; }
    arb_poly_struct *value() && = delete;
    const arb_poly_struct *value() const && = delete;
};
// The checked read-only ball range also serves coefficient storage. Root
// callers keep their existing spelling; both names share one implementation.
class ConstRootView;
using ConstComplexBallView = ConstRootView;

class ComplexPolynomial {
    acb_poly_t value_;

  public:
    ComplexPolynomial() { acb_poly_init(value_); }
    ~ComplexPolynomial() { acb_poly_clear(value_); }
    ComplexPolynomial(const ComplexPolynomial &x) : ComplexPolynomial() { acb_poly_set(value_, x.value_); }
    ComplexPolynomial(ComplexPolynomial &&x) noexcept : ComplexPolynomial() { acb_poly_swap(value_, x.value_); }
    ComplexPolynomial &operator=(ComplexPolynomial x) noexcept {
        acb_poly_swap(value_, x.value_);
        return *this;
    }
    acb_poly_struct * value() & noexcept { return value_; }
    const acb_poly_struct * value() const & noexcept { return value_; }
    acb_poly_struct * value() && = delete;
    const acb_poly_struct * value() const && = delete;
    // Borrows the current coefficient storage and extent, without copying or
    // retaining it. Resizing/reassignment that replaces storage invalidates the
    // view; ownership transfers preserve it while the destination retains it.
    ConstComplexBallView coefficients() const &;
    ConstComplexBallView coefficients() const && = delete;
};
class IntegerPolynomial {
    fmpz_poly_t value_;

  public:
    IntegerPolynomial() { fmpz_poly_init(value_); }
    ~IntegerPolynomial() { fmpz_poly_clear(value_); }
    IntegerPolynomial(const IntegerPolynomial &x) : IntegerPolynomial() { fmpz_poly_set(value_, x.value_); }
    IntegerPolynomial(IntegerPolynomial &&x) noexcept : IntegerPolynomial() { fmpz_poly_swap(value_, x.value_); }
    IntegerPolynomial &operator=(IntegerPolynomial x) noexcept {
        fmpz_poly_swap(value_, x.value_);
        return *this;
    }
    fmpz_poly_struct * value() & noexcept { return value_; }
    const fmpz_poly_struct * value() const & noexcept { return value_; }
    fmpz_poly_struct * value() && = delete;
    const fmpz_poly_struct * value() const && = delete;
};
class Factorization {
    fmpz_poly_factor_t value_;

  public:
    Factorization() { fmpz_poly_factor_init(value_); }
    ~Factorization() { fmpz_poly_factor_clear(value_); }
    Factorization(const Factorization &) = delete;
    Factorization &operator=(const Factorization &) = delete;
    fmpz_poly_factor_struct * value() & noexcept { return value_; }
    const fmpz_poly_factor_struct * value() const & noexcept { return value_; }
    fmpz_poly_factor_struct * value() && = delete;
    const fmpz_poly_factor_struct * value() const && = delete;
};
// Backend vectors are contiguous FLINT values, not arrays of C++ owning values.
// Bounds and allocation arithmetic live here so mathematical callers need not
// manipulate ownership or reinterpret object layouts.
inline slong checked_vector_size(slong n, size_t element_size) {
    if (n < 0 || static_cast<uint64_t>(n) > std::numeric_limits<size_t>::max() / element_size)
        throw std::length_error("numerical vector size exceeds addressable storage");
    return n;
}
class ConstRootView {
    acb_srcptr data_;
    slong size_;

  public:
    ConstRootView(acb_srcptr data, slong size)
        : data_(data), size_(checked_vector_size(size, sizeof(acb_struct))) {
        if (size && !data)
            throw std::invalid_argument("null root view");
    }
    slong size() const { return size_; }
    acb_srcptr operator[](slong i) const {
        if (i < 0 || i >= size_)
            throw std::out_of_range("root view index");
        return data_ + i;
    }
};
inline ConstComplexBallView ComplexPolynomial::coefficients() const & {
    return {value_->coeffs, acb_poly_length(value_)};
}
// These buffers own their allocation. Handles and views borrow it; create them
// only from a named owner and do not retain them past its lifetime or reassignment.
class RootVector {
    acb_ptr data_;
    slong size_;

  public:
    explicit RootVector(slong n) : data_(_acb_vec_init(checked_vector_size(n, sizeof(acb_struct)))), size_(n) {}
    ~RootVector() { _acb_vec_clear(data_, size_); }
    RootVector(const RootVector &) = delete;
    RootVector &operator=(const RootVector &) = delete;
    RootVector(RootVector &&x) noexcept
        : data_(std::exchange(x.data_, nullptr)), size_(std::exchange(x.size_, 0)) {}
    slong size() const noexcept { return size_; }
    acb_ptr data() & noexcept { return data_; }
    acb_srcptr data() const & noexcept { return data_; }
    acb_ptr data() && = delete;
    acb_srcptr data() const && = delete;
    acb_ptr operator[](slong i) & {
        if (i < 0 || i >= size_)
            throw std::out_of_range("root vector index");
        return data_ + i;
    }
    acb_srcptr operator[](slong i) const & { return view()[i]; }
    acb_ptr operator[](slong) && = delete;
    acb_srcptr operator[](slong) const && = delete;
    ConstRootView view() const & { return {data_, size_}; }
    ConstRootView view() const && = delete;
};
class RealVector {
    arb_ptr data_;
    slong size_;

  public:
    explicit RealVector(slong n) : data_(_arb_vec_init(checked_vector_size(n, sizeof(arb_struct)))), size_(n) {}
    ~RealVector() { _arb_vec_clear(data_, size_); }
    RealVector(const RealVector &) = delete;
    RealVector &operator=(const RealVector &) = delete;
    slong size() const noexcept { return size_; }
    arb_ptr data() & noexcept { return data_; }
    arb_srcptr data() const & noexcept { return data_; }
    arb_ptr data() && = delete;
    arb_srcptr data() const && = delete;
    arb_ptr operator[](slong i) & {
        if (i < 0 || i >= size_)
            throw std::out_of_range("real vector index");
        return data_ + i;
    }
    arb_srcptr operator[](slong i) const & {
        if (i < 0 || i >= size_)
            throw std::out_of_range("real vector index");
        return data_ + i;
    }
    arb_ptr operator[](slong) && = delete;
    arb_srcptr operator[](slong) const && = delete;
};

inline std::string integer_string(const fmpz_t value) {
    std::unique_ptr<char, decltype(&flint_free)> s(fmpz_get_str(nullptr, 10, value), &flint_free);
    return s.get();
}
} // namespace certpolysolve::detail
