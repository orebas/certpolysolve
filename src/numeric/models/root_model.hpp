#pragma once
// A local model of p on a chart (RootModel: a Taylor polynomial in the chart's coordinate with a
// uniform error bound on its disk), the error on a rectangle (root_model_error, local_omitted_band_bound) and the
// restriction to a smaller disk (restricted_root_model). Model bounds and coordinate changes are independent
// of the method constructing the model and of the certificate that uses it.
#include "core/number.hpp"
#include <memory>
#include <optional>

namespace certpolysolve::detail {
struct ModelErrorStatistics;
// p(z)/(z^lower_index * 2^scale) = polynomial(w) + error(w),
// z = center + radius*w, |w| <= 1. The error bound is uniform on that disk.
// Centers enclose exact roots-of-unity chart centers; all subsequent coordinate
// conversions retain their uncertainty. An omitted monomial requires z != 0.
class CoefficientEnvelope;
struct RootModel {
    ComplexBall center;
    RealBall radius, error;
    // Decomposed error supports diagnostics and the private local-band experiment.
    // The global error remains the fallback. A tighter certificate bound may also
    // change failure classification and subsequent precision selection.
    RealBall omitted_band_error, taylor_model_error;
    // With positive order k, E_T(w) = w^k H(w), |H| <= taylor_model_error
    // on the unit disk. Coefficient/FFT uncertainty remains in polynomial.
    // Zero means this provenance is unavailable, so retain the uniform bound.
    slong taylor_error_order = 0;
    bool error_components = false;
    // Shared immutable snapshot, independent of cursor/input lifetime. Present
    // when the local omitted-band bound is enabled.
    std::shared_ptr<const CoefficientEnvelope> error_envelope;
    ComplexPolynomial polynomial;
    slong lower_index = 0, upper_index = 0, scale = 0;
    // Optional normalized proposal disk covering the annulus's cells. Zero
    // means unavailable. It never shrinks this model's certificate domain.
    RealBall proposal_radius;
};
// A uniform bound on the entire original-coordinate rectangle. No value is
// returned when a required negative power could meet the origin.
std::optional<RealBall> local_omitted_band_bound(const RootModel &, const acb_t rectangle);
// Same outward error computation for certification and proposal stopping.
// Only the certificate's rectangle can establish a proof.
// For proof use, the caller first checks that its chart coordinates lie in
// |w| <= 1, where the stored Taylor remainder applies.
RealBall root_model_error(const RootModel &, const acb_t rectangle, slong precision,
                         RealBall *band = nullptr, ModelErrorStatistics * = nullptr,
                         RealBall *taylor = nullptr);
// Owned representation of the same function under w = fraction*v, 0 < fraction
// <= 1. Its error stays uniform on the new unit disk; proof provenance follows
// the coordinate change, independently of the proposal/certificate caller.
RootModel restricted_root_model(const RootModel &, const RealBall &fraction, slong precision);
} // namespace certpolysolve::detail
