#pragma once
#include "numeric/arithmetic/machine_arithmetic.hpp"
#include "core/raii.hpp"
#include <array>
#include <flint/gr.h>
#include <flint/nfloat.h>
#include <mpfr.h>

namespace certpolysolve::detail::machine {
// Fixed 128-bit candidate experiments. These types do not change any FLINT or
// MPFR default precision, exponent bounds, rounding mode or certificate path.
class MPReal {
    mpfr_t value_;

  public:
    MPReal() {
        mpfr_init2(value_, 128);
        mpfr_set_zero(value_, 1);
    }
    MPReal(double x) : MPReal() { mpfr_set_d(value_, x, MPFR_RNDN); }
    ~MPReal() { mpfr_clear(value_); }
    MPReal(const MPReal &x) : MPReal() { mpfr_set(value_, x.value_, MPFR_RNDN); }
    MPReal(MPReal &&x) noexcept : MPReal() { mpfr_swap(value_, x.value_); }
    MPReal &operator=(MPReal x) noexcept {
        mpfr_swap(value_, x.value_);
        return *this;
    }
    mpfr_ptr data() { return value_; }
    mpfr_srcptr data() const { return value_; }
    explicit operator double() const { return mpfr_get_d(value_, MPFR_RNDN); }
    friend MPReal operator-(const MPReal &a) {
        MPReal r;
        mpfr_neg(r.value_, a.value_, MPFR_RNDN);
        return r;
    }
    friend MPReal operator+(const MPReal &a, const MPReal &b) {
        MPReal r;
        mpfr_add(r.value_, a.value_, b.value_, MPFR_RNDN);
        return r;
    }
    friend MPReal operator-(const MPReal &a, const MPReal &b) {
        MPReal r;
        mpfr_sub(r.value_, a.value_, b.value_, MPFR_RNDN);
        return r;
    }
    friend MPReal operator*(const MPReal &a, const MPReal &b) {
        MPReal r;
        mpfr_mul(r.value_, a.value_, b.value_, MPFR_RNDN);
        return r;
    }
    friend MPReal operator/(const MPReal &a, const MPReal &b) {
        MPReal r;
        mpfr_div(r.value_, a.value_, b.value_, MPFR_RNDN);
        return r;
    }
    friend bool operator==(const MPReal &a, const MPReal &b) { return mpfr_equal_p(a.value_, b.value_); }
    friend bool operator<(const MPReal &a, const MPReal &b) { return mpfr_less_p(a.value_, b.value_); }
};
inline MPReal sqrt(const MPReal &a) {
    MPReal r;
    mpfr_sqrt(r.data(), a.data(), MPFR_RNDN);
    return r;
}
inline MPReal abs(const MPReal &a) {
    MPReal r;
    mpfr_abs(r.data(), a.data(), MPFR_RNDN);
    return r;
}
inline bool finite(const MPReal &a) { return mpfr_number_p(a.data()); }
inline double logabs(const MPReal &a) {
    if (mpfr_zero_p(a.data()))
        return -std::numeric_limits<double>::infinity();
    mpfr_exp_t e;
    const double m = mpfr_get_d_2exp(&e, a.data(), MPFR_RNDN);
    return std::log2(std::abs(m)) + static_cast<double>(e);
}
class NFloat {
    std::array<ulong, NFLOAT_HEADER_LIMBS + 128 / FLINT_BITS> value_{};
    struct Context {
        gr_ctx_t value;
        Context() {
            if (nfloat_ctx_init(value, 128, NFLOAT_ALLOW_INF | NFLOAT_ALLOW_NAN) != GR_SUCCESS)
                throw std::runtime_error("nfloat context");
        }
        ~Context() { gr_ctx_clear(value); }
    };
    static void check(NFloat &r, int status) {
        if (status != GR_SUCCESS)
            nfloat_nan(r.data(), context());
    }

  public:
    // Context initialization is synchronized by C++; every later operation
    // reads the same immutable precision/flags. Values own their packed limbs.
    static gr_ctx_struct *context() {
        static Context c;
        return c.value;
    }
    NFloat() { nfloat_init(data(), context()); }
    NFloat(double x) : NFloat() { check(*this, nfloat_set_d(data(), x, context())); }
    nfloat_ptr data() { return value_.data(); }
    nfloat_srcptr data() const { return value_.data(); }
    explicit operator double() const {
        double m;
        slong e;
        if (nfloat_get_d_2exp_si(&m, &e, data(), context()) != GR_SUCCESS)
            return NAN;
        return static_cast<double>(ScaledDouble(m, e));
    }
    friend NFloat operator-(const NFloat &a) {
        NFloat r;
        check(r, nfloat_neg(r.data(), a.data(), context()));
        return r;
    }
    friend NFloat operator+(const NFloat &a, const NFloat &b) {
        NFloat r;
        check(r, nfloat_add(r.data(), a.data(), b.data(), context()));
        return r;
    }
    friend NFloat operator-(const NFloat &a, const NFloat &b) {
        NFloat r;
        check(r, nfloat_sub(r.data(), a.data(), b.data(), context()));
        return r;
    }
    friend NFloat operator*(const NFloat &a, const NFloat &b) {
        NFloat r;
        check(r, nfloat_mul(r.data(), a.data(), b.data(), context()));
        return r;
    }
    friend NFloat operator/(const NFloat &a, const NFloat &b) {
        NFloat r;
        check(r, nfloat_div(r.data(), a.data(), b.data(), context()));
        return r;
    }
    friend bool operator==(const NFloat &a, const NFloat &b) {
        return nfloat_equal(a.data(), b.data(), context()) == T_TRUE;
    }
    friend bool operator<(const NFloat &a, const NFloat &b) {
        int c = 0;
        return nfloat_cmp(&c, a.data(), b.data(), context()) == GR_SUCCESS && c < 0;
    }
    friend NFloat sqrt(const NFloat &a) {
        NFloat r;
        check(r, nfloat_sqrt(r.data(), a.data(), context()));
        return r;
    }
};
inline NFloat abs(const NFloat &a) { return a < NFloat(0) ? -a : a; }
inline bool finite(const NFloat &a) { return !NFLOAT_IS_INF(a.data()) && !NFLOAT_IS_NAN(a.data()); }
inline double logabs(const NFloat &a) {
    if (NFLOAT_IS_ZERO(a.data()))
        return -std::numeric_limits<double>::infinity();
    double m;
    slong e;
    if (nfloat_get_d_2exp_si(&m, &e, a.data(), NFloat::context()) != GR_SUCCESS)
        return INFINITY;
    return std::log2(std::abs(m)) + static_cast<double>(e);
}

#if defined(CERTPOLYSOLVE_HAVE_BINARY128)
extern "C" __float128 sqrtf128(__float128);
extern "C" __float128 frexpf128(__float128, int *);
struct Binary128 {
    __float128 value = 0;
    Binary128() = default;
    Binary128(double x) : value(x) {}
    static Binary128 raw(__float128 x) {
        Binary128 r;
        r.value = x;
        return r;
    }
    explicit operator double() const { return static_cast<double>(value); }
    friend Binary128 operator-(Binary128 a) { return raw(-a.value); }
    friend Binary128 operator+(Binary128 a, Binary128 b) { return raw(a.value + b.value); }
    friend Binary128 operator-(Binary128 a, Binary128 b) { return raw(a.value - b.value); }
    friend Binary128 operator*(Binary128 a, Binary128 b) { return raw(a.value * b.value); }
    friend Binary128 operator/(Binary128 a, Binary128 b) { return raw(a.value / b.value); }
    friend bool operator==(Binary128 a, Binary128 b) { return a.value == b.value; }
    friend bool operator<(Binary128 a, Binary128 b) { return a.value < b.value; }
};
inline Binary128 sqrt(Binary128 a) { return Binary128::raw(sqrtf128(a.value)); }
inline Binary128 abs(Binary128 a) { return a < Binary128(0) ? -a : a; }
inline bool finite(Binary128 a) { return __builtin_isfinite(a.value); }
inline double logabs(Binary128 a) {
    if (a.value == 0)
        return -INFINITY;
    int e;
    const auto m = frexpf128(a.value, &e);
    return std::log2(std::abs(static_cast<double>(m))) + e;
}
#endif
} // namespace certpolysolve::detail::machine
namespace certpolysolve::detail {
inline bool from_arf(machine::MPReal &out, const arf_t in) {
    arf_get_mpfr(out.data(), in, MPFR_RNDN);
    return machine::finite(out) && (out != machine::MPReal(0) || arf_is_zero(in));
}
inline void to_arf(arf_t out, const machine::MPReal &in) { arf_set_mpfr(out, in.data()); }
inline bool from_arf(machine::NFloat &out, const arf_t in) {
    return nfloat_set_arf(out.data(), in, machine::NFloat::context()) == GR_SUCCESS &&
           machine::finite(out) && (out != machine::NFloat(0) || arf_is_zero(in));
}
inline void to_arf(arf_t out, const machine::NFloat &in) {
    if (nfloat_get_arf(out, in.data(), machine::NFloat::context()) != GR_SUCCESS)
        arf_nan(out);
}
#if defined(CERTPOLYSOLVE_HAVE_BINARY128)
inline bool from_arf(machine::Binary128 &out, const arf_t in) {
    machine::MPReal converted;
    mpfr_prec_round(converted.data(), 113, MPFR_RNDN);
    arf_get_mpfr(converted.data(), in, MPFR_RNDN);
    out.value = mpfr_get_float128(converted.data(), MPFR_RNDN);
    return machine::finite(out) && (out.value != 0 || arf_is_zero(in));
}
inline void to_arf(arf_t out, machine::Binary128 in) {
    machine::MPReal converted;
    mpfr_set_float128(converted.data(), in.value, MPFR_RNDN);
    arf_set_mpfr(out, converted.data());
}
#endif
} // namespace certpolysolve::detail
