#pragma once
// Exact values as text, for the C API's additions of stage 3a (docs/api-spec-0.5.md; the C++ layer's numbers,
// root ends and polynomials print, compare and convert through them, so a caller needs no big-number library):
// a root end m 2^e parsed from its decimal parts, formatted, written as a fraction or as "m@e"; a coefficient
// written exactly in the coefficient grammar, so that the text sets the same coefficient again.
#include <certpolysolve/certpolysolve.h>
#include "core/number.hpp"
#include <stdexcept>
#include <string>

namespace certpolysolve::detail {
// m 2^e from a root end's decimal mantissa and exponent; std::invalid_argument for a malformed or null part.
void parse_dyadic(arf_t out, cps_dyadic_view end);
// Decimal with `digits` significant digits, correctly rounded to nearest with ties to even (the exact value set in
// MPFR and printed with "%.*RNg").
std::string dyadic_format(const arf_t x, unsigned digits);
// The exact value as "p/q" or "p" in lowest terms.
std::string dyadic_rational(const arf_t x);
// The exact value as "m@e" (m odd), or "0".
std::string dyadic_text(const arf_t x);
// A coefficient part's exact value in the coefficient grammar, a text that parses back to the same value: "p/q" or "p"
// in lowest terms, or "m@e" (m odd) for a binary exponent beyond 64 bits. The form follows how the value is stored, so
// one value entered in two forms (1@100 and its decimal digits) can print as two texts.
std::string number_text(const Number &n);
} // namespace certpolysolve::detail
