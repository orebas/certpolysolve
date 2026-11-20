#pragma once
// Which build of a floating-point kernel a call runs: the
// reference, ordinary IEEE operations in the order the source writes them, or a variant for instructions this CPU
// has. A variant makes the reference's operations in the reference's order at every point, so it returns the
// reference's bits and the choice changes the time of a method and nothing else in it; tests/test_kernel_variants.cpp
// runs every variant the machine supports against the reference. One decision reads the variant: Auto's order of
// the certificate attempt and the piecewise method, through the variant's prices (price_variant, below).
//
// The variant today is Avx2Fma (x86-64, GCC and Clang builds): four points at a time in AVX2 instructions, std::fma
// as the FMA instruction. Its kernels are the compensated certificate's evaluation
// (certify/attempts/compensated_kernel.hpp), the machine iteration's dense Horner evaluation
// (numeric/candidates/horner_kernel.hpp) and the terms of its direct interaction sum
// (numeric/candidates/interaction_kernel.hpp); the measured gains are in dev/docs/simd-kernels-20261117.md.
// A kernel's file takes the instructions by a function's target attribute, never by a compiler flag for the file:
// a flag would also build the inline functions of every header the file includes for those instructions, and the
// linker may keep such a copy for the whole library.
#include "platform/ieee_arithmetic.hpp"
#include <cstdint>

// Whether this build has the Avx2Fma kernels at all.
#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
#define CERTPOLYSOLVE_HAVE_AVX2_FMA_KERNELS 1
#endif

namespace certpolysolve::detail {
enum class KernelVariant : uint8_t { Reference, Avx2Fma };
// Whether this build has the variant and this CPU, under this operating system, runs it.
bool kernel_variant_available(KernelVariant variant) noexcept;
// The variant the kernels run: the best available one, unless set_kernel_variant chose another.
KernelVariant kernel_variant() noexcept;
// For the tests: the kernels run `variant` from now on, on every thread; an unavailable variant leaves the choice as
// it was. Returns the variant in force before. A solve under way on another thread goes on with the new choice,
// which changes none of its results.
KernelVariant set_kernel_variant(KernelVariant variant) noexcept;
// The variant whose prices Auto's selector reads (core/heuristics.hpp, AttemptPrices): the variant the kernels run,
// or the reference where the prices are pinned. For the tests: set_reference_prices pins or frees them, on every
// thread, and returns the pin before. Nothing else pins them (the library reads no environment): a caller who
// needs the routes of every machine forces a method.
KernelVariant price_variant() noexcept;
bool set_reference_prices(bool pinned) noexcept;
} // namespace certpolysolve::detail
