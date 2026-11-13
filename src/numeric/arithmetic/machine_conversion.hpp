#pragma once
// FLINT ball midpoints and machine values. A successful conversion preserves
// nonzero coordinates in the target arithmetic; it does not certify the ball.
// Precision and rounding remain explicit in the conversion implementations.
#include "numeric/arithmetic/arf_binary64.hpp"
#include "numeric/arithmetic/machine_arithmetic.hpp"
#include "core/raii.hpp"
#ifdef CERTPOLYSOLVE_EXPERIMENTAL_ARITHMETIC
#include "numeric/arithmetic/experimental_arithmetic.hpp"
#endif

namespace certpolysolve::detail {
inline bool from_arf(double &out, const arf_t in) {
    out = arf_to_binary64(in, ARF_RND_NEAR);
    return std::isfinite(out) && (out != 0 || arf_is_zero(in));
}
inline bool from_arf(machine::DoubleDouble &out, const arf_t in) {
    if (!from_arf(out.hi, in))
        return false;
    Float high, remainder;
    arf_set_d(high.value(), out.hi);
    arf_sub(remainder.value(), in, high.value(), 128, ARF_RND_NEAR);
    return from_arf(out.lo, remainder.value());
}
inline bool from_arf(machine::ScaledDouble &out, const arf_t in) {
    if (!arf_is_finite(in))
        return false;
    if (arf_is_zero(in)) {
        out = {};
        return true;
    }
    const slong e = arf_abs_bound_lt_2exp_si(in);
    if (e < -(int64_t(1) << 59) || e > (int64_t(1) << 59))
        return false;
    Float scaled;
    arf_mul_2exp_si(scaled.value(), in, -e);
    out = {arf_to_binary64(scaled.value(), ARF_RND_NEAR), e};
    return machine::finite(out) && out.mantissa != 0;
}
inline void to_arf(arf_t out, double x) { arf_set_d(out, x); }
inline void to_arf(arf_t out, machine::DoubleDouble x) {
    Float low;
    arf_set_d(out, x.hi);
    arf_set_d(low.value(), x.lo);
    arf_add(out, out, low.value(), 128, ARF_RND_NEAR);
}
inline void to_arf(arf_t out, machine::ScaledDouble x) {
    arf_set_d(out, x.mantissa);
    arf_mul_2exp_si(out, out, static_cast<slong>(x.exponent));
}
template <class R> bool from_ball(machine::Complex<R> &out, const acb_t in, slong scale = 0) {
    Float coordinate;
    arf_mul_2exp_si(coordinate.value(), arb_midref(acb_realref(in)), scale);
    if (!from_arf(out.re, coordinate.value()))
        return false;
    arf_mul_2exp_si(coordinate.value(), arb_midref(acb_imagref(in)), scale);
    return from_arf(out.im, coordinate.value());
}
template <class R> void to_ball(acb_t out, machine::Complex<R> z, slong scale = 0) {
    Float coordinate;
    to_arf(coordinate.value(), z.re);
    arf_mul_2exp_si(coordinate.value(), coordinate.value(), scale);
    arb_set_arf(acb_realref(out), coordinate.value());
    to_arf(coordinate.value(), z.im);
    arf_mul_2exp_si(coordinate.value(), coordinate.value(), scale);
    arb_set_arf(acb_imagref(out), coordinate.value());
}
} // namespace certpolysolve::detail
