#pragma once
// The two kinds of a discrete Fourier transform's output, kept apart by type: ApproximateTransform (rounded
// estimates, the established FFT libraries of the giant-FFT screen, benchmarks/giant-fft/) and EnclosedTransform
// (balls that contain the transform: FixedTransform's, fixed_transform.hpp). The bounded binary64 transform BEFFT's
// adapter was here until W8's S2 removed BEFFT support (23 October 2026; research/rejected/befft/README.md).
#include "core/raii.hpp"

namespace certpolysolve::detail {
// Ordinary transforms cannot be passed off as certified model coefficients.
struct ApproximateTransform {
    // Exact storage of rounded floating-point estimates, NOT error enclosures.
    RootVector values;
    explicit ApproximateTransform(slong length)
      : values(length) {}
};
struct EnclosedTransform {
    RootVector values;
    explicit EnclosedTransform(slong length)
      : values(length) {}
};
} // namespace certpolysolve::detail
