// The global driver's work on one factor (FactorState, factor_state.hpp): exact centring (prepare_secular_center,
// refine_center), tracing exact preparation (observe_preparation), the precision stop's bound (separation_bound), and
// the enclosers of its sources (owned_encloser, owned_exact_taylor, shared_factor_polynomial, attach_origin).
//
// prepare_secular_center moves the origin of a factor of degree at least centering_min_degree, real or (under
// complex_centering) Gaussian, to its root centroid rounded to initial_center_bits (round_center) by an exact Taylor
// shift (center_integer, center_gaussian), at most once (center_rejected). Before the shift it refuses, in this order,
// a centre whose translate cannot gain enough height (center_height_allows_gain, or the scale-invariant admission) or
// cannot be stored, one that would cost roots near the old origin (unseparated_origin_roots, under its flag) or (under scale_invariant_centering) moves too little to gain enough, and one whose translate's
// constant term already bounds its height too high (translated_height_lower_bound). It keeps the translate when the
// monic coefficient height falls by at least the bit length of the degree, at scale zero or, under
// scale_invariant_centering, at the root scale without rising as much at any scale (height_at, worst_height_excess).
// refine_center (under centering_refinement) then recentres at the residual centroid, rounded alike, while that at
// least halves the geometric mean of the root distances: their product is the modulus of the constant coefficient over
// the leading one, so the exact test is |q(t)| 2^n <= |q(0)|. Only candidate coordinates change; certificates remain
// exact. separation_bound gives Mahler's bound on root separation and Cauchy's on root size, from bit lengths, for an
// integer polynomial of degree 2 or more.
//
// State: none here. The centring writes the FactorState it is passed (centre, centred polynomials, candidate centres),
// one per factor, owned by SolveState::Impl; the enclosers keep their own or shared copies of its polynomial.
// Entry points, all called from solver.cpp only: observe_preparation (preprocess, factor_step); prepare_secular_center
// (factor_step, before a dense factor's secular stage is built); refine_center (factor_step's classical centring, and
// prepare_secular_center); the enclosers and attach_origin (factor_step's sources); separation_bound (advance).
// Not here: shifts and height estimates (centering.cpp); the flags' settings (maintained_method_policy, session.cpp).
// Design notes: docs/equivalence-screen-findings-20260928.md (section B, the refinement),
// docs/selector-campaign-20260928.md (scale-invariant centering), docs/api.md (sufficient_work_bits).
#include "solve/factor_state.hpp"
#include "numeric/exact/centering.hpp"
#include "core/telemetry.hpp"

namespace certpolysolve::detail {
void observe_preparation(const EnginePolicy &policy, const char *stage,
                         const PreparationDiagnostic &diagnostic, uint64_t budget,
                         std::chrono::steady_clock::time_point started, uint64_t active) {
    if (!policy.trace) return;
    TraceEvent event{stage, "end", 0, active, nullptr, nullptr};
    event.reason = diagnostic.reason;
    event.exact_budget = budget;
    event.required_bits_lower_bound = diagnostic.required_bits_lower_bound;
    event.exact_input_bits = diagnostic.input_bits;
    event.exact_output_bits = diagnostic.output_bits;
    event.preparation_accepted = diagnostic.accepted;
    event.preparation_ns = static_cast<uint64_t>(std::chrono::duration_cast<std::chrono::nanoseconds>(
        std::chrono::steady_clock::now() - started).count());
    policy.trace(policy.trace_data, event);
}
std::optional<SeparationBound> separation_bound(const fmpz_poly_t p) {
    const slong n = fmpz_poly_degree(p);
    if (n < 2)
        return std::nullopt;
    uint64_t height = 0;
    for (slong i = 0; i <= n; ++i)
        height = std::max<uint64_t>(height, fmpz_bits(p->coeffs + i));
    const double log_norm = static_cast<double>(height) + 0.5 * std::log2(static_cast<double>(n + 1));
    const double separation = static_cast<double>(n - 1) * log_norm +
                              0.5 * static_cast<double>(n + 2) * std::log2(static_cast<double>(n));
    // |a_i / a_n| < 2^(height - bits(a_n) + 1), so a root is below 2^(height - bits(a_n) + 2).
    const uint64_t lead = fmpz_bits(p->coeffs + n);
    return SeparationBound{static_cast<uint64_t>(std::ceil(separation)), height + 2 - std::min(height + 1, lead)};
}
namespace {
// Exact |q(t)|^2 for a Gaussian-rational polynomial q, by Horner's rule.
void gaussian_norm_at(Rational &out, const GaussianPolynomial &q, const GaussianRational &t) {
    Rational re, im, a, b;
    for (size_t k = q.size(); k-- > 0;) {
        fmpq_mul(a.value(), re.value(), t.real.value());
        fmpq_submul(a.value(), im.value(), t.imag.value());
        fmpq_mul(b.value(), re.value(), t.imag.value());
        fmpq_addmul(b.value(), im.value(), t.real.value());
        fmpq_add(re.value(), a.value(), q[k].real.value());
        fmpq_add(im.value(), b.value(), q[k].imag.value());
    }
    fmpq_mul(out.value(), re.value(), re.value());
    fmpq_addmul(out.value(), im.value(), im.value());
}
} // namespace
// An accepted centre is the root centroid rounded to initial_center_bits, so
// that the exact shift stays cheap. When the centroid is not such a short
// dyadic, the rounding leaves an offset, and a tight root cloud can sit many
// of its own radii away from the new origin: the monomial basis then needs
// about n log2(offset/spread) extra bits. The coefficient-height rule cannot
// see this, because every coefficient of a polynomial whose roots lie well
// inside the unit disk is below its leading one.
//
// |q_0/q_n| is the product of the root distances from the origin, and moving
// the origin to t makes it |q(t)/q_n| (q_n is unchanged). So recentre at the
// residual centroid, again rounded to initial_center_bits, while that at
// least halves the geometric mean of the root distances: |q(t)| 2^n <= |q(0)|.
// For an off-centre cloud the product falls from about offset^n towards
// spread^n; a step towards a few distant roots, or towards a single root
// near the centroid, gains far less than one bit per root and is refused.
// Every comparison is exact, so decisions are equivariant under x -> 2^s x.
// Only candidate coordinates change; certificates remain exact.
void refine_center(FactorState &factor, const cps_options &options, certpolysolve::detail::SearchStatistics &statistics,
                   const EnginePolicy &policy) {
    const bool real = factor.kind == FactorState::Real;
    for (;;) {
        const slong n = real ? fmpz_poly_degree(factor.centered.value())
                             : static_cast<slong>(factor.centered_gaussian.size()) - 1;
        if (n < 1 || !options.exact_effort_fits(resources::multiply(resources::nonnegative(n + 1),
                                                                     resources::nonnegative(n + 1))))
            return;
        // Residual centroid -q_{n-1}/(n q_n), rounded like the first centre.
        GaussianRational step;
        if (real) {
            fmpz_neg(fmpq_numref(step.real.value()), factor.centered.value()->coeffs + n - 1);
            fmpz_mul_si(fmpq_denref(step.real.value()), factor.centered.value()->coeffs + n, n);
            fmpq_canonicalise(step.real.value());
        } else {
            const auto &lead = factor.centered_gaussian[static_cast<size_t>(n)];
            const auto &next = factor.centered_gaussian[static_cast<size_t>(n - 1)];
            Rational norm, temp;
            fmpq_mul(norm.value(), lead.real.value(), lead.real.value());
            fmpq_addmul(norm.value(), lead.imag.value(), lead.imag.value());
            fmpq_mul_si(norm.value(), norm.value(), -n);
            // -(a+bi)/(n(c+di)) = -((ac+bd) + (bc-ad)i) / (n(c^2+d^2)).
            fmpq_mul(temp.value(), next.real.value(), lead.real.value());
            fmpq_addmul(temp.value(), next.imag.value(), lead.imag.value());
            fmpq_div(step.real.value(), temp.value(), norm.value());
            fmpq_mul(temp.value(), next.imag.value(), lead.real.value());
            fmpq_submul(temp.value(), next.real.value(), lead.imag.value());
            fmpq_div(step.imag.value(), temp.value(), norm.value());
        }
        ComplexBall target;
        for (auto [part, coordinate] : {std::pair{step.real.value(), acb_realref(target.value())},
                                        std::pair{step.imag.value(), acb_imagref(target.value())}}) {
            arf_set_fmpq(arb_midref(coordinate), part, heuristics::initial_center_bits, ARF_RND_NEAR);
            arf_get_fmpq(part, arb_midref(coordinate));
        }
        if (step.zero()) return;
        const auto owned = real ? center_owned_bits(factor.integer.value()) : center_owned_bits(factor.gaussian);
        if (!center_materialization_fits(target.value(), options.allow_exponent_expansion ?
                options.exact_storage_limit() : std::min(options.exact_storage_limit(), owned)))
            return;
        const auto started = std::chrono::steady_clock::now();
        PreparationDiagnostic diagnostic;
        bool gain;
        if (real) {
            Rational at_step, at_origin;
            fmpz_poly_evaluate_fmpq(at_step.value(), factor.centered.value(), step.real.value());
            fmpq_abs(at_step.value(), at_step.value());
            fmpq_set_fmpz(at_origin.value(), factor.centered.value()->coeffs);
            fmpq_abs(at_origin.value(), at_origin.value());
            fmpq_mul_2exp(at_step.value(), at_step.value(), static_cast<ulong>(n));
            gain = !fmpq_is_zero(at_step.value()) && fmpq_cmp(at_step.value(), at_origin.value()) <= 0;
        } else {
            Rational at_step, at_origin;
            gaussian_norm_at(at_step, factor.centered_gaussian, step);
            gaussian_norm_at(at_origin, factor.centered_gaussian, GaussianRational{});
            fmpq_mul_2exp(at_step.value(), at_step.value(), 2 * static_cast<ulong>(n));
            gain = !fmpq_is_zero(at_step.value()) && fmpq_cmp(at_step.value(), at_origin.value()) <= 0;
        }
        if (!gain) {
            diagnostic.reason = "centering-refinement-declined";
            observe_preparation(policy, "exact-centering", diagnostic, options.exact_budget, started,
                                resources::nonnegative(n));
            return;
        }
        IntegerPolynomial shifted;
        GaussianPolynomial shifted_gaussian;
        {
            StageTimer timer(statistics.taylor_shift_ns);
            TraceScope construction(&policy, "exact-centering-construction", factor.precision,
                                    resources::nonnegative(n));
            const bool converted = real
                ? center_integer(shifted.value(), factor.centered.value(), step.real.value(),
                                 options.exact_storage_limit(), policy.trace ? &diagnostic : nullptr)
                : center_gaussian(shifted_gaussian, factor.centered_gaussian, step, options.exact_storage_limit(),
                                  policy.trace ? &diagnostic : nullptr);
            if (!converted) {
                if (!real) diagnostic.reason = "gaussian-centering-unavailable";
                observe_preparation(policy, "exact-centering", diagnostic, options.exact_budget, started,
                                    resources::nonnegative(n));
                return;
            }
        }
        diagnostic.reason = "centering-refined";
        diagnostic.accepted = true;
        observe_preparation(policy, "exact-centering", diagnostic, options.exact_budget, started,
                            resources::nonnegative(n));
        fmpq_add(factor.center.real.value(), factor.center.real.value(), step.real.value());
        fmpq_add(factor.center.imag.value(), factor.center.imag.value(), step.imag.value());
        if (real) factor.centered = std::move(shifted);
        else factor.centered_gaussian = std::move(shifted_gaussian);
        if (factor.candidates.centers)
            for (slong i = 0; i < factor.candidates.centers->size(); ++i) {
                acb_sub(factor.candidates.centers->data() + i, factor.candidates.centers->data() + i,
                        target.value(), static_cast<slong>(factor.precision));
                acb_get_mid(factor.candidates.centers->data() + i, factor.candidates.centers->data() + i);
            }
        if (factor.candidates.study) factor.candidates.study->coordinates_changed();
        ++statistics.taylor_shifts;
    }
}
void prepare_secular_center(FactorState &factor, const acb_poly_t polynomial,
                            const cps_options &options, certpolysolve::detail::SearchStatistics &statistics,
                            const EnginePolicy &policy) {
    const slong degree = acb_poly_degree(polynomial);
    if ((factor.kind != FactorState::Real && !(policy.complex_centering && factor.kind == FactorState::Gaussian)) ||
        factor.have_center || factor.center_rejected ||
        degree < heuristics::centering_min_degree || acb_is_zero(polynomial->coeffs + degree - 1) ||
        !options.exact_effort_fits(resources::multiply(resources::nonnegative(degree + 1), resources::nonnegative(degree + 1))))
        return;
    factor.center_rejected = true;
    const slong required_gain = static_cast<slong>(FLINT_BIT_COUNT(static_cast<ulong>(degree)));
    const bool real = factor.kind == FactorState::Real;
    const bool plain_admits = real ? center_height_allows_gain(factor.integer.value()) :
                                     center_height_allows_gain(factor.gaussian);
    // With the flag, a scale-invariant route may also admit and accept the
    // centering. It measures heights at the original polynomial's root scale.
    HeightProfile before;
    double scale = 0, before_height = 0;
    if (policy.scale_invariant_centering || policy.centering_keeps_origin_scales)
        before = real ? height_profile(factor.integer.value()) : height_profile(factor.gaussian);
    if (policy.scale_invariant_centering) {
        scale = root_scale(before);
        before_height = height_at(before, scale);
    }
    const bool scaled_admits = policy.scale_invariant_centering && before_height >= static_cast<double>(required_gain);
    if (!plain_admits && !scaled_admits) {
        PreparationDiagnostic diagnostic;
        diagnostic.reason="centering-height-gain-impossible";
        observe_preparation(policy,"exact-centering",diagnostic,options.exact_budget,
                            std::chrono::steady_clock::now(),resources::nonnegative(degree));
        return;
    }
    ComplexBall target;
    acb_div(target.value(), polynomial->coeffs + degree - 1, polynomial->coeffs + degree, 64);
    acb_div_ui(target.value(), target.value(), static_cast<ulong>(degree), 64);
    acb_neg(target.value(), target.value());
    round_center(target.value(), heuristics::initial_center_bits);
    const auto owned = factor.kind == FactorState::Real ? center_owned_bits(factor.integer.value()) : center_owned_bits(factor.gaussian);
    const auto center_limit = options.allow_exponent_expansion ? options.exact_storage_limit() :
        std::min(options.exact_storage_limit(), owned);
    if (!center_materialization_fits(target.value(), center_limit)) return;
    // With the flag, a translate that carries roots out of the binary64 stage's reach of the new origin with
    // no certified group to hold them is refused before the exact shift is built
    // (unseparated_origin_roots). The reach is the stage's guard: a node's coordinates are raised past the
    // stage once it lies nearer the original origin than 2^-(stage_bits - bits) of the centre
    // (SecularRefinement's coordinate rule).
    if (policy.centering_keeps_origin_scales) {
        Magnitude size;
        acb_get_mag(size.value(), target.value());
        const double reach = static_cast<double>(SecularMachineWorkspace::stage_bits(options.bits)) -
                             static_cast<double>(options.bits);
        if (!mag_is_zero(size.value()) &&
            unseparated_origin_roots(before, mag_get_d_log2_approx(size.value()), reach) > 0) {
            PreparationDiagnostic diagnostic;
            diagnostic.reason = "centering-origin-scales";
            observe_preparation(policy, "exact-centering", diagnostic, options.exact_budget,
                                std::chrono::steady_clock::now(), resources::nonnegative(degree));
            return;
        }
    }
    // With the flag, skip a shift that no route can accept. Translating by c
    // changes the coefficient of x^m, at scale 2^t, by at most
    // sum_{j>m} |c_j| C(j,m) |c|^(j-m) 2^(tm) <= T ((1-e)^-(n+1) - 1), where
    // T = max_k |c_k| 2^(tk) is the largest term and e = |c| 2^-t. When
    // 2(n+1)|c| <= 2^t the factor is below exp(1/2 / (1 - 1/130)) - 1 < 0.66
    // (n >= centering_min_degree), so the largest term keeps over a third of
    // its size: the true gain at scale t is under 1.6 bits, and the height
    // estimates (within 2.5 bits each) cannot show required_gain >=
    // bitcount(centering_min_degree) = 7. The plain rule measures at t = 0,
    // the scale-invariant route at its root scale.
    bool plain_futile = false, scaled_futile = false;
    if (policy.scale_invariant_centering) {
        slong center_bits = WORD_MIN / 2; // |c| < 2^(center_bits + 1), from the larger part
        for (const arf_struct *part : {arb_midref(acb_realref(target.value())), arb_midref(acb_imagref(target.value()))})
            if (!arf_is_zero(part)) center_bits = std::max(center_bits, arf_abs_bound_lt_2exp_si(part));
        const slong limit = center_bits + 2 + static_cast<slong>(FLINT_BIT_COUNT(static_cast<ulong>(degree + 1)));
        plain_futile = limit <= 0;                              // 2(n+1)|c| < 2^limit
        scaled_futile = static_cast<double>(limit) <= scale;
        if ((!plain_admits || plain_futile) && (!scaled_admits || scaled_futile)) {
            PreparationDiagnostic diagnostic;
            diagnostic.reason = "centering-shift-too-small";
            observe_preparation(policy, "exact-centering", diagnostic, options.exact_budget,
                                std::chrono::steady_clock::now(), resources::nonnegative(degree));
            return;
        }
    }
    GaussianRational center;
    arf_get_fmpq(center.real.value(), arb_midref(acb_realref(target.value())));
    arf_get_fmpq(center.imag.value(), arb_midref(acb_imagref(target.value())));
    const auto height = [](const fmpz_poly_t q) {
        return FLINT_ABS(fmpz_poly_max_bits(q)) - static_cast<slong>(fmpz_bits(fmpz_poly_lead(q)));
    };
    // A translate that neither rule below can select is not built. Its constant term is p(center)
    // (up to the translate's common factor), so one ball evaluation bounds its height from below
    // (translated_height_lower_bound): the plain rule needs height(p) - height(translate) >=
    // required_gain, and the scale-invariant one the same at the root scale, where the constant
    // term's share of the translate's height is that bound less scale * degree. A centroid dragged
    // far from most roots by a few distant ones is the case. Every translate refused here would be
    // refused by the same rules after its construction, so only the time changes. The profile's bit
    // lengths are exact for an integer polynomial; one bit of slack covers the rounding of
    // scale * degree. History: docs/code-history-20261021.md, factor_state.cpp 1.
    if (real) {
        if (const auto floor_height = translated_height_lower_bound(factor.integer.value(), center.real.value())) {
            const bool plain_cannot = height(factor.integer.value()) - *floor_height < required_gain;
            const double scaled_rounding = 1;
            const bool scaled_cannot = !policy.scale_invariant_centering ||
                before_height - (static_cast<double>(*floor_height) - scale * static_cast<double>(degree)) <
                    static_cast<double>(required_gain) - scaled_rounding;
            if ((!plain_admits || plain_futile || plain_cannot) && (!scaled_admits || scaled_futile || scaled_cannot)) {
                PreparationDiagnostic diagnostic;
                diagnostic.reason = "centering-constant-term-bound";
                observe_preparation(policy, "exact-centering", diagnostic, options.exact_budget,
                                    std::chrono::steady_clock::now(), resources::nonnegative(degree));
                return;
            }
        }
    }
    IntegerPolynomial shifted;
    GaussianPolynomial shifted_gaussian;
    StageTimer timer(statistics.taylor_shift_ns);
    PreparationDiagnostic diagnostic;
    const auto started = std::chrono::steady_clock::now();
    TraceScope construction(&policy,"exact-centering-construction",factor.precision,resources::nonnegative(degree));
    const bool converted = factor.kind == FactorState::Real
        ? center_integer(shifted.value(), factor.integer.value(), center.real.value(),
                         options.exact_storage_limit(), policy.trace ? &diagnostic : nullptr)
        : center_gaussian(shifted_gaussian, factor.gaussian, center, options.exact_storage_limit(),
                          policy.trace ? &diagnostic : nullptr);
    const bool plain_gains = converted && plain_admits && !plain_futile &&
        (real ? height(factor.integer.value()) - height(shifted.value()) >= required_gain
              : gaussian_center_improves_height(factor.gaussian, shifted_gaussian));
    std::optional<HeightProfile> after;
    const auto translated_profile = [&]() -> const HeightProfile & {
        if (!after) after = real ? height_profile(shifted.value()) : height_profile(shifted_gaussian);
        return *after;
    };
    const bool plain_selected = plain_gains;
    // The scale-invariant route measures gains and losses in the plain rule's
    // unit, required_gain: the translate must gain at least that much height
    // at the original's root scale, and lose less than that at every scale.
    // An offset cloud gains everywhere. Centering Laguerre at its mean root
    // gains at the root scale but loses more at small scales, where
    // |p(c)| >> |p(0)|, so it keeps the plain refusal.
    // History: docs/code-history-20261021.md, factor_state.cpp 2.
    bool scaled_selected = false;
    if (converted && scaled_admits && !scaled_futile && !plain_selected) {
        const double gain = before_height - height_at(translated_profile(), scale);
        scaled_selected = gain >= static_cast<double>(required_gain) &&
                          worst_height_excess(translated_profile(), before) < static_cast<double>(required_gain);
    }
    const bool selected = converted && (plain_selected || scaled_selected);
    if (converted)
        diagnostic.reason = selected ? "centering-selected" : "centering-no-height-improvement";
    diagnostic.accepted = selected;
    observe_preparation(policy, "exact-centering", diagnostic, options.exact_budget, started,
                        static_cast<uint64_t>(degree));
    if (!selected)
        return;
    factor.center = std::move(center);
    factor.centered = std::move(shifted);
    factor.centered_gaussian = std::move(shifted_gaussian);
    factor.have_center = true;
    if (factor.candidates.centers)
        for (slong i = 0; i < degree; ++i) {
            acb_sub(factor.candidates.centers->data() + i, factor.candidates.centers->data() + i,
                    target.value(), static_cast<slong>(factor.precision));
            acb_get_mid(factor.candidates.centers->data() + i, factor.candidates.centers->data() + i);
        }
    ++statistics.taylor_shifts;
    if (policy.centering_refinement) refine_center(factor, options, statistics, policy);
}
PolynomialEncloser owned_encloser(const FactorState &factor, const DensePolynomial &input) {
    if (factor.kind == FactorState::Real)
        return [polynomial = factor.integer](acb_poly_t out, slong bits) {
            acb_poly_set_fmpz_poly(out, polynomial.value(), bits);
        };
    if (factor.kind == FactorState::Gaussian)
        return [polynomial = factor.gaussian](acb_poly_t out, slong bits) {
            enclose_gaussian(out, polynomial, bits);
        };
    return [coefficients = input.coefficients](acb_poly_t out, slong bits) {
        ComplexBall coefficient;
        acb_poly_zero(out);
        for (size_t i = 0; i < coefficients.size(); ++i) {
            coefficients[i].enclose(coefficient.value(), bits);
            acb_poly_set_coeff_acb(out, static_cast<slong>(i), coefficient.value());
        }
    };
}
ExactTaylorEncloser owned_exact_taylor(const FactorState &factor, uint64_t budget, bool allow_expansion, uint64_t operations
                                     , std::function<double()> *work_reader
                                     ) {
    if (work_reader) *work_reader = {};
    std::shared_ptr<ExactTaylorCache> model;
    if (factor.kind == FactorState::Real)
        model = std::make_shared<ExactTaylorCache>(factor.integer, budget, allow_expansion, operations);
    else if (factor.kind == FactorState::Gaussian)
        model = std::make_shared<ExactTaylorCache>(factor.gaussian, budget, allow_expansion, operations);
    else
        return {};
    if (work_reader) {
        model->research_service_meter(true);
        *work_reader = [model] { return model->research_work(); };
    }
    return [model = std::move(model)](acb_poly_t out, const acb_t center, slong precision) {
        return model->enclose(out, center, precision);
    };
}
std::optional<SharedFactorPolynomial> shared_factor_polynomial(const FactorState &factor, bool centered) {
    if (factor.kind == FactorState::Real)
        return SharedFactorPolynomial{std::make_shared<const IntegerPolynomial>(centered ? factor.centered : factor.integer), {}};
    if (factor.kind == FactorState::Gaussian)
        return SharedFactorPolynomial{{}, std::make_shared<const GaussianPolynomial>(centered ? factor.centered_gaussian : factor.gaussian)};
    return std::nullopt;
}
// The uncentred factor and its centre for a centred source (PolynomialSourceOptions::origin), with both
// polynomials' log2 coefficient sizes. Nothing when the centre is not exactly representable (centres are
// low-denominator dyadics, round_center).
void attach_origin(PolynomialSourceOptions &options, const FactorState &factor) {
    if (factor.kind != FactorState::Real && factor.kind != FactorState::Gaussian)
        return;
    ComplexBall offset;
    for (bool imaginary : {false, true}) {
        const fmpq *part = imaginary ? factor.center.imag.value() : factor.center.real.value();
        const auto bits = static_cast<slong>(fmpz_bits(fmpq_numref(part)) + fmpz_bits(fmpq_denref(part))) + 64;
        arb_set_fmpq(imaginary ? acb_imagref(offset.value()) : acb_realref(offset.value()), part, bits);
    }
    if (!acb_is_exact(offset.value()))
        return;
    const auto log2_sizes = [](const HeightProfile &profile) {
        std::vector<double> sizes(profile.bits.size(), -HUGE_VAL);
        for (size_t k = 0; k < sizes.size(); ++k)
            if (profile.present[k])
                sizes[k] = static_cast<double>(profile.bits[k]);
        return sizes;
    };
    const bool real = factor.kind == FactorState::Real;
    options.centred_log2 = log2_sizes(real ? height_profile(factor.centered.value()) : height_profile(factor.centered_gaussian));
    options.origin_log2 = log2_sizes(real ? height_profile(factor.integer.value()) : height_profile(factor.gaussian));
    options.origin = shared_factor_polynomial(factor, false)->encloser();
    acb_set(options.origin_offset.value(), offset.value());
}
} // namespace certpolysolve::detail
