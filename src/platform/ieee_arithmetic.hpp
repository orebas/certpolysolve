#pragma once
// The library's floating-point arithmetic is IEEE 754 binary64 as written: each operation rounded once, to nearest,
// in the order of the source; no product and sum fused into one rounding unless written (std::fma); infinities and
// NaNs kept. The certificates' error-free transformations (a sum's error as (a + b) - a and the like) and the
// candidate stage's finiteness tests rest on it. -ffast-math, -ffinite-math-only and /fp:fast let the compiler assume
// otherwise: it may simplify an error term to zero or remove a finiteness test, and a certificate could be false. Such
// a build is refused here. Contraction is turned off by the build's -ffp-contract=off. Clang also takes that from the
// pragma below in its default mode (not under an explicit -ffp-contract=fast); GCC has no such pragma and contracts
// C++ by default, so a GCC build must keep the option.
#if defined(__FAST_MATH__) || (defined(__FINITE_MATH_ONLY__) && __FINITE_MATH_ONLY__) || defined(_M_FP_FAST)
#error "CertPolySolve needs IEEE 754 arithmetic as written: build it without -ffast-math, -ffinite-math-only or /fp:fast"
#endif
#if defined(__clang__)
#pragma STDC FP_CONTRACT OFF
#endif
