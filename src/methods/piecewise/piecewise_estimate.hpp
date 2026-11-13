#pragma once
// The selector's piecewise estimate (PiecewiseModelEstimate, estimate_piecewise_models on a coefficient list or on
// the input object's terms): piecewise's rings walked without building a model, for routing (auto_selection.cpp,
// piecewise_admission.cpp). Defined in piecewise_models.cpp. Part of piecewise_roots.hpp.
#include "core/number.hpp"
#include <cstdint>
#include <functional>
#include <optional>
#include <span>
#include <utility>

namespace certpolysolve::detail {
// Scheduling estimate only. This reuses the solver's ring/band geometry but
// builds no local polynomials and makes no mathematical root-count claim.
struct PiecewiseModelEstimate {
    bool available = false;
    // The walk stopped at the caller's bound: the sums cover the rings so far.
    bool truncated = false;
    uint64_t rings = 0, models = 0, coefficients = 0, maximum_degree = 0;
    // Counts of operations, exact; a sum saturates at 2^64 - 1, beyond any bound a caller sets.
    uint64_t construction_work = 0, local_work = 0;
    // The envelope's bend in bits (its hull's largest excess over the chord of its ends), a statistic; the cap is
    // decided exactly, in integers (bend_exceeded).
    double envelope_bend = 0;
    bool bend_exceeded = false;
};
PiecewiseModelEstimate estimate_piecewise_models(
    const acb_poly_t polynomial, slong model_bits,
    // The bend above which the estimate stops (bend_exceeded); none: no cap.
    std::optional<uint64_t> maximum_envelope_bend = std::nullopt,
    // Called with construction_work + local_work after each ring; true stops
    // the walk (the estimate is then `truncated`).
    const std::function<bool(uint64_t)> &beyond = {});
// The same from a polynomial's nonzero terms (exponent, ball) in increasing order, the first at exponent 0: the
// estimate of the polynomial they make, without its zero coefficients (C40's second step).
PiecewiseModelEstimate estimate_piecewise_models(
    std::span<const std::pair<slong, ComplexBall>> terms, slong model_bits,
    std::optional<uint64_t> maximum_envelope_bend = std::nullopt, const std::function<bool(uint64_t)> &beyond = {});
} // namespace certpolysolve::detail
