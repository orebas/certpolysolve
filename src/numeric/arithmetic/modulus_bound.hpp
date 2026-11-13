#pragma once
// The upper bound of a complex ball's modulus that the certificates compare, at modulus_bound_bits bits (upper()):
// the proposals, the piecewise state, the chart certificate and the certificate attempt use it.
#include "core/number.hpp"

namespace certpolysolve::detail {
// Bits of the upper bounds on moduli that the certificates compare (upper()).
constexpr slong modulus_bound_bits = 64;
inline void upper(arb_t out, const acb_t in) {
    Float f;
    acb_get_abs_ubound_arf(f.value(), in, modulus_bound_bits);
    arb_set_arf(out, f.value());
}
} // namespace certpolysolve::detail
