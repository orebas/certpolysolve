// Local omitted-coefficient and Taylor-remainder bounds, and exact coordinate
// restriction of an existing model. These operations neither propose nor publish roots.
// Called by chart certification, piecewise axis work and proposal stopping.
#include "numeric/models/root_model.hpp"
#include "numeric/models/coefficient_envelope.hpp"
#include "numeric/models/model_error_statistics.hpp"
#include <stdexcept>

namespace certpolysolve::detail {
std::optional<RealBall> local_omitted_band_bound(const RootModel &model, const acb_t rectangle) {
    if (!model.error_envelope || !acb_is_finite(rectangle)) return std::nullopt;
    const auto &envelope = *model.error_envelope;
    const slong degree = envelope.degree();
    if (model.lower_index < 0 || model.upper_index < model.lower_index || model.upper_index > degree)
        return std::nullopt;
    Float lower, upper;
    acb_get_abs_lbound_arf(lower.value(), rectangle, 96);
    acb_get_abs_ubound_arf(upper.value(), rectangle, 96);
    if (model.lower_index && arf_sgn(lower.value()) <= 0) return std::nullopt;
    RealBall lower_radius, upper_radius, bound, part;
    arb_set_arf(lower_radius.value(), lower.value());
    arb_set_arf(upper_radius.value(), upper.value());
    envelope.bound(bound.value(), 0, model.lower_index - 1, lower_radius.value(), model.lower_index, model.scale);
    envelope.bound(part.value(), model.upper_index + 1, degree, upper_radius.value(), model.lower_index, model.scale);
    arb_add(bound.value(), bound.value(), part.value(), 96);
    if (!arb_is_finite(bound.value())) return std::nullopt;
    return bound;
}

RealBall root_model_error(const RootModel &model, const acb_t rectangle, slong precision,
                         RealBall *band, ModelErrorStatistics *statistics, RealBall *taylor) {
    RealBall error = model.error;
    RealBall applied_band = model.omitted_band_error;
    RealBall applied_taylor = model.taylor_model_error;
    if (model.error_envelope) {
        if (statistics) ++statistics->local_band_attempts;
        if (auto local = local_omitted_band_bound(model, rectangle)) {
            RealBall total;
            arb_add(total.value(), local->value(), model.taylor_model_error.value(), precision);
            if (arb_lt(total.value(), error.value())) {
                error = std::move(total);
                applied_band = *local;
                if (statistics) ++statistics->local_band_improvements;
            }
        }
    }
    if (model.taylor_error_order > 0 && arb_is_positive(model.radius.value()) &&
        acb_is_finite(rectangle)) {
        if (statistics) ++statistics->local_taylor_attempts;
        ComplexBall local;
        acb_sub(local.value(), rectangle, model.center.value(), precision);
        acb_div_arb(local.value(), local.value(), model.radius.value(), precision);
        Float outward;
        acb_get_abs_ubound_arf(outward.value(), local.value(), precision);
        if (arf_is_finite(outward.value()) && arf_cmp_ui(outward.value(), 1) <= 0) {
            // The whole rectangle is inside |w| <= r <= 1. The omitted
            // Taylor terms start at k, so their absolute sum is <= T*r^k.
            // For r > 1 this inequality is invalid; retain the old bound.
            RealBall radial, tail, total;
            arb_set_arf(radial.value(), outward.value());
            arb_pow_ui(tail.value(), radial.value(), static_cast<ulong>(model.taylor_error_order), precision);
            arb_mul(tail.value(), tail.value(), model.taylor_model_error.value(), precision);
            arb_add(total.value(), applied_band.value(), tail.value(), precision);
            if (arb_lt(total.value(), error.value())) {
                error = std::move(total);
                applied_taylor = std::move(tail);
                if (statistics) ++statistics->local_taylor_improvements;
            }
        }
    }
    if (band) *band = std::move(applied_band);
    if (taylor) *taylor = std::move(applied_taylor);
    return error;
}

RootModel restricted_root_model(const RootModel &model, const RealBall &fraction, slong precision) {
    RealBall one;
    arb_one(one.value());
    if (!arb_is_positive(fraction.value()) || !arb_le(fraction.value(), one.value()))
        throw std::invalid_argument("model restriction requires 0 < fraction <= 1");
    RootModel inner = model;
    arb_mul(inner.radius.value(), inner.radius.value(), fraction.value(), precision);
    RealBall power;
    arb_one(power.value());
    for (slong j = 0; j < acb_poly_length(inner.polynomial.value()); ++j) {
        acb_mul_arb(inner.polynomial.value()->coeffs + j,
                    inner.polynomial.value()->coeffs + j, power.value(), precision);
        arb_mul(power.value(), power.value(), fraction.value(), precision);
    }
    if (inner.taylor_error_order > 0) {
        arb_pow_ui(power.value(), fraction.value(), static_cast<ulong>(inner.taylor_error_order), precision);
        arb_mul(inner.taylor_model_error.value(), inner.taylor_model_error.value(), power.value(), precision);
        // This remains a uniform bound on |v| <= 1, not a pointwise bound
        // reused by unit-circle counts or derivative Cauchy estimates.
        arb_add(inner.error.value(), inner.omitted_band_error.value(), inner.taylor_model_error.value(), precision);
    }
    return inner;
}

} // namespace certpolysolve::detail
