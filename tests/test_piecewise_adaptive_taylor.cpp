#include "numeric/arithmetic/polynomial_arithmetic.hpp"
#include "exact_winding_api.hpp"
#include "methods/piecewise/piecewise_roots.hpp"
#include <iostream>

namespace d = certpolysolve::detail;
using reference::require;

// Compare with a full, independently computed Taylor shift, including all the
// terms the streaming construction omits. Input balls exercise coefficient
// uncertainty as well as cancellation and complex coefficients.
static void full_shift_oracle(bool uncertain, bool fixed_fft = false) {
    d::ComplexPolynomial input;
    d::ComplexBall coefficient;
    constexpr slong degree = 257, precision = 256;
    for (slong j=0; j<=degree; ++j) {
        acb_set_si_si(coefficient.value(), j%7-3, j%5-2);
        if (j==0 || j==degree) acb_one(coefficient.value());
        if (uncertain && j && j<degree) {
            arb_add_error_2exp_si(acb_realref(coefficient.value()), -220);
            arb_add_error_2exp_si(acb_imagref(coefficient.value()), -219);
        }
        acb_poly_set_coeff_acb(input.value(), j, coefficient.value());
    }
    auto config = d::maintained_piecewise_config();
    config.adaptive_taylor = true;
    config.error_components = true;
    config.source_finisher = false;
    config.fixed_fft = fixed_fft;
    d::PiecewiseRootStatistics stats;
    cps_options options;
    size_t checked = 0, truncated = 0;
    const auto inspect = [&](const d::RootModel &model) {
        d::ComplexPolynomial band, shifted;
        for (slong j=model.lower_index; j<=model.upper_index; ++j) {
            acb_poly_get_coeff_acb(coefficient.value(), input.value(), j);
            acb_mul_2exp_si(coefficient.value(), coefficient.value(), -model.scale);
            acb_poly_set_coeff_acb(band.value(), j-model.lower_index, coefficient.value());
        }
        acb_poly_taylor_shift(shifted.value(), band.value(), model.center.value(), precision+128);
        d::RealBall power, norm, tail;
        d::Float upper;
        arb_one(power.value());
        const slong local_degree = acb_poly_degree(model.polynomial.value());
        for (slong j=0; j<=acb_poly_degree(shifted.value()); ++j) {
            acb_poly_get_coeff_acb(coefficient.value(), shifted.value(), j);
            acb_mul_arb(coefficient.value(), coefficient.value(), power.value(), precision+128);
            if (j<=local_degree) {
                d::ComplexBall streamed;
                acb_poly_get_coeff_acb(streamed.value(), model.polynomial.value(), j);
                require(acb_overlaps(coefficient.value(), streamed.value()),
                        "adaptive coefficient disagrees with full Taylor shift");
            } else {
                acb_get_abs_ubound_arf(upper.value(), coefficient.value(), precision+128);
                arb_set_arf(norm.value(), upper.value());
                arb_add(tail.value(), tail.value(), norm.value(), precision+128);
            }
            arb_mul(power.value(), power.value(), model.radius.value(), precision+128);
        }
        require(arb_le(tail.value(), model.taylor_model_error.value()),
                "stored adaptive tail excludes a full-shift coefficient tail");
        require(model.taylor_error_order > local_degree,
                "Taylor provenance does not cover the constructed polynomial");
        // Check the complete original polynomial, including both omitted
        // coefficient bands. Coefficient-ball uncertainty is distinct from
        // the analytic remainder and is charged explicitly below.
        for (slong radial_shift : {0, -1, -2}) for (unsigned direction=0;direction<4;++direction) {
            d::ComplexBall w,z,source,value,midpoint,monomial,difference;
            acb_set_si_si(w.value(), direction==0?1:direction==2?-1:0,
                                      direction==1?1:direction==3?-1:0);
            acb_mul_2exp_si(w.value(),w.value(),radial_shift);
            acb_mul_arb(z.value(),w.value(),model.radius.value(),precision+128);
            acb_add(z.value(),z.value(),model.center.value(),precision+128);
            acb_poly_evaluate(source.value(),input.value(),z.value(),precision+128);
            acb_pow_ui(monomial.value(),z.value(),static_cast<ulong>(model.lower_index),precision+128);
            acb_div(source.value(),source.value(),monomial.value(),precision+128);
            acb_mul_2exp_si(source.value(),source.value(),-model.scale);
            acb_poly_evaluate(value.value(),model.polynomial.value(),w.value(),precision+128);
            acb_get_mid(midpoint.value(),value.value());
            acb_sub(difference.value(),source.value(),midpoint.value(),precision+128);
            d::RealBall applied_tail;
            d::RealBall bound=d::root_model_error(model,z.value(),precision+128,nullptr,nullptr,&applied_tail);
            d::RealBall real_radius,imag_radius,residual;
            arb_get_rad_arb(real_radius.value(),acb_realref(value.value()));
            arb_get_rad_arb(imag_radius.value(),acb_imagref(value.value()));
            arb_add(bound.value(),bound.value(),real_radius.value(),precision+128);
            arb_add(bound.value(),bound.value(),imag_radius.value(),precision+128);
            acb_get_abs_ubound_arf(upper.value(),difference.value(),precision+128);
            arb_set_arf(residual.value(),upper.value());
            require(arb_le(residual.value(),bound.value()),
                    "full original polynomial lies outside analytic plus arithmetic model error");
            if (radial_shift<0) {
                // Directly sum the full-shift omitted coefficients at this
                // radius. This oracle does not use the geometric row bound.
                d::RealBall actual_tail,radial_power,norm;
                for (slong j=local_degree+1;j<=acb_poly_degree(shifted.value());++j) {
                    acb_poly_get_coeff_acb(coefficient.value(),shifted.value(),j);
                    acb_get_abs_ubound_arf(upper.value(),coefficient.value(),precision+128);
                    arb_set_arf(norm.value(),upper.value());
                    arb_pow_ui(radial_power.value(),model.radius.value(),static_cast<ulong>(j),precision+128);
                    arb_mul_2exp_si(radial_power.value(),radial_power.value(),j*radial_shift);
                    arb_mul(norm.value(),norm.value(),radial_power.value(),precision+128);
                    arb_add(actual_tail.value(),actual_tail.value(),norm.value(),precision+128);
                }
                require(arb_le(actual_tail.value(),applied_tail.value()),
                        "pointwise Taylor remainder excludes exact shifted tail");
            }
        }
        // Coordinate invariance is tested on the same original-coordinate
        // rectangle, including uncertain geometry, through both owned models.
        d::RealBall fraction;arb_set_d(fraction.value(),0.5);
        auto uncertain_model=model;
        arb_add_error_2exp_si(acb_realref(uncertain_model.center.value()),-220);
        arb_add_error_2exp_si(uncertain_model.radius.value(),-221);
        const auto restricted=d::restricted_root_model(uncertain_model,fraction,precision);
        d::ComplexBall rectangle=uncertain_model.center;
        d::RealBall offset;
        arb_mul_2exp_si(offset.value(),uncertain_model.radius.value(),-3);
        arb_add(acb_realref(rectangle.value()),acb_realref(rectangle.value()),offset.value(),precision);
        arb_add_error_2exp_si(acb_imagref(rectangle.value()),-150);
        const auto outer_bound=d::root_model_error(uncertain_model,rectangle.value(),precision);
        const auto inner_bound=d::root_model_error(restricted,rectangle.value(),precision);
        require(arb_overlaps(outer_bound.value(),inner_bound.value()),
                "Taylor bound changed under exact dyadic chart restriction");
        if (local_degree < std::min<slong>(model.upper_index-model.lower_index, 128)) ++truncated;
        ++checked;
        return true;
    };
    // Real and imaginary axis charts from every radial band exercise the FFT
    // sign/index convention as well as dynamic range, without thousands of
    // dense shifts. A positive-real chart alone cannot detect a phase reversal.
    require(d::piecewise_root_models(input.value(),16,precision,options,stats,inspect,
        [](const auto &m) { return (arb_is_zero(acb_imagref(m.center.value())) &&
                                   arb_is_positive(acb_realref(m.center.value()))) ||
                                  (arb_is_zero(acb_realref(m.center.value())) &&
                                   arb_is_positive(acb_imagref(m.center.value()))); }, config),
        "adaptive model construction failed");
    require(checked>0 && truncated>0, "adaptive oracle never exercised early truncation");
    // The fixed-point transform takes every transform (256 bits at these lengths are within its largest word), so
    // no ball transform remains; without it none is fixed-point.
    if (fixed_fft) require(stats.fft_fixed_accepted>0 && stats.fft_fixed_declined==0 &&
                           stats.fft_fixed_accepted==stats.fft_transforms,
                           "full-source oracle left a ball transform beside the fixed-point one");
    else require(stats.fft_fixed_accepted==0, "a fixed-point transform ran with its switch off");
}

static void exact_roots_and_continuation() {
    // Exact, independently known roots include a close pair and complex affine
    // transforms. Resumption and a stronger query must retain the same roots.
    using namespace reference;
    const Rational epsilon=power2(-18);
    std::vector<Complex> roots{{1}, {1+epsilon}, {-2,1}, {3,-2}, {0,2}, {-3,-3}};
    for (auto &z:roots) z=z*Complex(Rational(3)/8,Rational(1)/4)+Complex(2,-1);
    const auto input=reference::input(from_roots(roots));
    {
        auto config=d::maintained_piecewise_config();
        config.adaptive_taylor = config.interior_proposals = true;
        d::PiecewiseRootState state(certpolysolve::detail::dense(*input.get()),config);
        cps_options options; options.goal=CPS_ROOT_APPROXIMATION; options.work_units=7;
        (void)state.advance(options);
        options.work_units=UINT64_MAX;
        for (auto bits:{53u,96u}) {
            options.bits=bits;
            const auto result=state.advance(options);
            require(result.goal_met && result.roots.size()==roots.size(), "adaptive continuation lost roots");
            for (const auto &r:result.roots) {
                size_t count=0;
                for (const auto &z:roots) {
                    d::Rational real,imag;
                    fmpq_set_str(real.value(),z.real.str().c_str(),10);
                    fmpq_set_str(imag.value(),z.imag.str().c_str(),10);
                    count+=arb_contains_fmpq(acb_realref(r.box.value()),real.value()) &&
                           arb_contains_fmpq(acb_imagref(r.box.value()),imag.value());
                }
                require(count==r.count,"adaptive certificate excludes its exact planted root");
            }
        }
    }
}

int main() {
    certpolysolve::ThreadScope lifetime;
    try {
        for (bool uncertain:{false,true}) for (bool fixed:{false,true}) full_shift_oracle(uncertain,fixed);
        exact_roots_and_continuation();
        std::cout << "adaptive Taylor full-shift, uncertain-input and continuation checks passed\n";
    } catch (const std::exception &e) { std::cerr << e.what() << '\n'; return 1; }
}
