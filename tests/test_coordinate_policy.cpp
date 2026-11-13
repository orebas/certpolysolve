#include "certify/attempts/newton_cover.hpp"
#include "core/engine_policy.hpp"
#include "core/engine_statistics.hpp"
#include "solve/policies.hpp"
#include "internal.hpp"
#include <bit>
#include <certpolysolve/certpolysolve.hpp>
#include "numeric/exact/centering.hpp"
#include "core/heuristics.hpp"
#include "methods/piecewise/piecewise_roots.hpp"
#include "numeric/evaluate/structured_source.hpp"
#include "numeric/candidates/study_candidates.hpp"
#include "exact_winding_api.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <stdexcept>

namespace d = certpolysolve::detail;
static void require(bool ok, const char *why) { if (!ok) throw std::runtime_error(why); }

static void height() {
    d::GaussianPolynomial p(3), q(3);
    fmpq_set_si(p[0].real.value(), -65, 1);
    fmpq_set_si(p[1].imag.value(), -16, 1);
    fmpq_one(p[2].real.value());
    fmpq_set_si(q[0].real.value(), -1, 1);
    fmpq_one(q[2].real.value());
    require(d::gaussian_center_improves_height(p,q), "useful complex translation refused");
    require(d::center_height_allows_gain(p), "precheck refused useful translation");
    require(!d::center_height_allows_gain(q), "zero monic height can improve");
    require(!d::gaussian_center_improves_height(q,p), "height growth accepted");
    require(!d::gaussian_center_improves_height(q,q), "unchanged coordinates accepted");
    // Multiplication by 2i cannot change a monic-height decision.
    for (auto &z : q) {
        fmpq_swap(z.real.value(),z.imag.value());
        fmpq_mul_si(z.real.value(),z.real.value(),-2);
        fmpq_mul_si(z.imag.value(),z.imag.value(),2);
    }
    require(d::gaussian_center_improves_height(p,q), "scalar-dependent height decision");
}

static void exact_shift_identities() {
    using reference::Complex;
    using reference::Rational;
    const std::vector<Complex> roots{{Rational(1)/3,Rational(2)/7},{-4,2},{2,-3},{0,1},{5,0}};
    auto p=reference::from_roots(roots);
    d::GaussianPolynomial source(p.size());
    for (size_t j=0;j<p.size();++j) {
        p[j]=p[j]*Complex(3,-2);
        require(!fmpq_set_str(source[j].real.value(),p[j].real.str().c_str(),10),"test rational import");
        require(!fmpq_set_str(source[j].imag.value(),p[j].imag.str().c_str(),10),"test rational import");
    }
    const auto value=[](const d::GaussianRational &z) {
        const auto rational=[](const d::Rational &q) {
            return Rational(d::integer_string(fmpq_numref(q.value())))/Rational(d::integer_string(fmpq_denref(q.value())));
        };
        return Complex(rational(z.real),rational(z.imag));
    };
    for (const auto c : std::vector<Complex>{{0,0},{Rational(5)/3,0},{0,Rational(-7)/5},{Rational(5)/3,Rational(-7)/5}}) {
        d::GaussianRational center;
        fmpq_set_str(center.real.value(),c.real.str().c_str(),10);
        fmpq_set_str(center.imag.value(),c.imag.str().c_str(),10);
        d::GaussianPolynomial shifted;
        require(d::center_gaussian(shifted,source,center,UINT64_MAX),"rational complex shift failed");
        auto shifted_roots=roots;
        for (auto &root:shifted_roots) root=root-c;
        const auto expected=reference::from_roots(shifted_roots);
        for (size_t j=0;j<expected.size();++j)
            require(reference::zero(value(shifted[j])*expected.back()-expected[j]*value(shifted.back())),
                    "FLINT complex shift differs from independent planted-root expansion");
        const auto retained=value(shifted[0]);
        require(!d::center_gaussian(shifted,source,center,1),"tiny exact budget ignored");
        require(reference::zero(value(shifted[0])-retained),"refused shift mutated retained output");
    }
}

static void translated_roots(cps_method method) {
    constexpr slong n = 96;
    // x^3 ((x-3i/8)^96-1)^2. Independent analytic roots, exact repeated
    // factors and a zero factor exercise both prepared-factor coordinates
    // and restoration. No fixture roots enter the candidate initializer.
    d::GaussianPolynomial base(static_cast<size_t>(n+1)), shifted;
    fmpq_set_si(base[0].real.value(),-1,1);
    fmpq_one(base[n].real.value());
    d::GaussianRational center;
    fmpq_set_si(center.imag.value(),-3,8);
    require(d::center_gaussian(shifted,base,center,UINT64_MAX), "exact test translation failed");
    d::ComplexPolynomial exact, squared;
    d::enclose_gaussian(exact.value(),shifted,4096);
    acb_poly_mul(squared.value(),exact.value(),exact.value(),4096);
    for (slong j=0;j<=2*n;++j)
        require(acb_is_exact(squared.value()->coeffs+j), "test coefficient arithmetic not exact");
    certpolysolve::detail::DensePolynomial input;
    input.coefficients.resize(static_cast<size_t>(2*n+4));
    for (slong j=0;j<=2*n;++j) {
        auto &c=input.coefficients[static_cast<size_t>(j+3)];
        arf_get_fmpq(c.real.significand.value(),arb_midref(acb_realref(squared.value()->coeffs+j)));
        arf_get_fmpq(c.imag.significand.value(),arb_midref(acb_imagref(squared.value()->coeffs+j)));
    }
    auto policy=d::maintained_method_policy(method);
    require(policy.complex_centering,"maintained method omitted qualified complex centering");
    d::SolveState state(std::move(input),policy);
    cps_options options;
    options.method=method;options.goal=CPS_ROOT_APPROXIMATION;
    // One solve with the units the one-unit calls of a session used to share (one advance a solve state since 18 October 2026).
    options.bits=53;options.max_work_bits=4096;options.work_units=1000;
    cps_result result=state.advance(&options);
    require(result.goal_met,"translated prepared factors failed to finish");
    require(d::disjoint(result.roots),"translated published boxes overlap");
    d::RootVector known(n);
    _acb_vec_unit_roots(known.data(),n,n,512);
    d::ComplexBall offset;
    arb_set_si(acb_imagref(offset.value()),3);
    acb_mul_2exp_si(offset.value(),offset.value(),-3);
    for (slong i=0;i<n;++i) acb_add(known[i],known[i],offset.value(),512);
    uint64_t total=0;
    for (const auto &r:result.roots) {
        uint64_t count=acb_contains_zero(r.box.value()) ? 3 : 0;
        for (slong i=0;i<n;++i) if (acb_contains(r.box.value(),known[i])) count+=2;
        require(count==r.count,"analytic translated-root count differs");
        require(d::meets_goal(r,options,static_cast<slong>(result.work_bits)),"original-coordinate goal lost");
        total+=count;
    }
    require(total==2*n+3,"factor restoration lost roots");
}

static void refusal_reason() {
    certpolysolve::detail::DensePolynomial p;p.coefficients.resize(3);
    // Past the compact-exponent rule's allowance (one precision cap a number), so the import is refused.
    p.coefficients[0].real=d::Number("1@1000000");
    p.coefficients[2].real=d::Number("1");
    cps_options o;o.work_units=2;o.goal=CPS_ROOT_APPROXIMATION;
    d::SolveState state(p);
    auto result=state.advance(&o);
    require(result.diagnostics.find("compact-exponent-expansion-needs-explicit-allowance")!=std::string::npos,
            "real import refusal hidden by inapplicable Gaussian preparation");
}

static void failed_cover_retains_only_proposals() {
    std::unique_ptr<d::RootVector> hints;
    std::vector<d::Root> cover(1);
    cover.front().count=37;
    {
        d::ComplexPolynomial p;
        acb_poly_set_coeff_si(p.value(),0,-1);
        acb_poly_set_coeff_si(p.value(),2,1);
        arb_add_error_2exp_si(acb_realref(p.value()->coeffs),-20);
        d::PolynomialSource source([&](acb_poly_t out,slong) {acb_poly_set(out,p.value());},2,true);
        cps_options o;o.goal=CPS_ROOT_APPROXIMATION;o.bits=53;o.max_work_bits=85;
        certpolysolve::detail::CandidateStatistics statistics{};
        require(!d::machine_newton_cover(cover,source,nullptr,o,85,true,statistics,true,nullptr,&hints),
                "machine guesses certified irreducibly uncertain coefficients");
    }
    require(cover.size()==1 && cover.front().count==37,"failed shortcut changed proof cover");
    require(hints && hints->size()==2,"failed proof discarded useful machine proposals");
    require(acb_is_finite((*hints)[0]) && acb_is_finite((*hints)[1]),"retained proposals lost ownership");
}

// Maintained secular centering at the original root scale. A flat polynomial
// offset by about its own root radius has no plain monic height to remove
// (every coefficient is below the leading one), so the plain rule refuses it;
// its translate is flat at the original root scale and must be centered. A
// partial centering that leaves most of that height (Laguerre, whose roots
// spread over orders of magnitude) keeps the plain rule's refusal.
static std::vector<std::pair<std::string,bool>> secular_centering(const fmpz_poly_t p,bool scale_invariant,
                                                                 bool origin_scales=true) {
    certpolysolve::detail::DensePolynomial input;
    input.coefficients.resize(static_cast<size_t>(fmpz_poly_length(p)));
    for (slong k=0;k<fmpz_poly_length(p);++k)
        fmpq_set_fmpz(input.coefficients[static_cast<size_t>(k)].real.significand.value(),p->coeffs+k);
    auto policy=d::maintained_method_policy(CPS_METHOD_SECULAR);
    require(policy.scale_invariant_centering,"maintained secular omitted scale-invariant centering");
    policy.scale_invariant_centering=scale_invariant;
    require(policy.centering_keeps_origin_scales,"maintained secular omitted the origin-scales rule");
    policy.centering_keeps_origin_scales=origin_scales;
    std::vector<std::pair<std::string,bool>> seen;
    policy.trace_data=&seen;
    policy.trace=[](void *data,const d::TraceEvent &e) {
        if (std::string_view(e.stage)=="exact-centering")
            static_cast<std::vector<std::pair<std::string,bool>> *>(data)->emplace_back(e.reason,e.preparation_accepted);
    };
    d::SolveState state(std::move(input),policy);
    cps_options options;
    options.method=CPS_METHOD_SECULAR;options.goal=CPS_ROOT_APPROXIMATION;
    // One solve with the units the one-unit calls of a session used to share (one advance a solve state since 18 October 2026).
    options.bits=53;options.max_work_bits=4096;options.work_units=1000;
    cps_result result=state.advance(&options);
    require(result.goal_met,"centering fixture failed to finish");
    require(d::disjoint(result.roots),"centering fixture boxes overlap");
    uint64_t total=0;
    for (const auto &r:result.roots) total+=r.count;
    require(total==static_cast<uint64_t>(fmpz_poly_degree(p)),"centering fixture lost roots");
    return seen;
}

// A centre from ball arithmetic (centering.hpp, round_center): each part to
// `bits` significant bits, and a part that may be zero is zero.
static void center_rounding() {
    const slong bits=d::heuristics::initial_center_bits;
    const auto midpoint_is=[](const arb_t part,const d::Rational &expected) {
        d::Rational value;
        arf_get_fmpq(value.value(),arb_midref(part));
        return fmpq_equal(value.value(),expected.value())!=0;
    };
    // Parts known to be nonzero keep their own rounding, as before the rule:
    // one part alone, either way round (ties included: 511/512 and 513/1024
    // at 8 bits), and two parts of very different size (the centroid of
    // near_multiple-193 scaled and moved by 3i/8: the -3/8 must survive
    // beside 1,485.7, because the translation by it makes the polynomial real).
    for (const char *text:{"1/3","-5/7","255/256","511/512","513/1024","-1023/4096","3/1048576","7","1000001/7","-257/2"}) {
        d::Rational exact;
        require(fmpq_set_str(exact.value(),text,10)==0,"bad fixture");
        d::Float expected;
        for (const bool imaginary:{false,true}) {
            d::ComplexBall ball;
            const auto part=imaginary ? acb_imagref(ball.value()) : acb_realref(ball.value());
            const auto other=imaginary ? acb_realref(ball.value()) : acb_imagref(ball.value());
            arb_set_fmpq(part,exact.value(),64);
            arf_set_round(expected.value(),arb_midref(part),bits,ARF_RND_NEAR);
            d::round_center(ball.value(),bits);
            require(arf_equal(arb_midref(part),expected.value()) && arb_is_zero(other),
                    "a centre with one part is not rounded as before");
        }
    }
    {
        d::ComplexBall ball;
        arb_set_d(acb_realref(ball.value()),1485.7130829015543);
        arb_set_si(acb_imagref(ball.value()),-3);
        arb_mul_2exp_si(acb_imagref(ball.value()),acb_imagref(ball.value()),-3);
        d::round_center(ball.value(),bits);
        d::Rational real,imaginary;
        fmpq_set_si(real.value(),1488,1);
        fmpq_set_si(imaginary.value(),-3,8);
        require(midpoint_is(acb_realref(ball.value()),real) &&
                    midpoint_is(acb_imagref(ball.value()),imaginary),
                "a small part known to be nonzero did not keep its own rounding");
    }
    // A part whose enclosure contains zero is zero: the centre of resultant_50
    // times 3+4i (2,157/2^18, and an imaginary part of -1.7e-22 that is
    // rounding noise). 2,157/2^18 at 8 bits is 135/2^14.
    {
        d::ComplexBall ball;
        arb_set_si(acb_realref(ball.value()),2157);
        arb_mul_2exp_si(acb_realref(ball.value()),acb_realref(ball.value()),-18);
        arb_set_si(acb_imagref(ball.value()),-202);
        arb_mul_2exp_si(acb_imagref(ball.value()),acb_imagref(ball.value()),-80);
        d::ComplexBall known;
        acb_set(known.value(),ball.value());
        mag_set_ui_2exp_si(arb_radref(acb_imagref(ball.value())),1,-72);
        d::round_center(ball.value(),bits);
        d::Rational expected,tiny;
        fmpq_set_si(expected.value(),135,1<<14);
        require(midpoint_is(acb_realref(ball.value()),expected) &&
                    arb_is_zero(acb_imagref(ball.value())),
                "a part that may be zero survived the rounding");
        // The same small part with no uncertainty stays: what the enclosure
        // says decides, not the part's size.
        d::round_center(known.value(),bits);
        fmpq_set_si(tiny.value(),-202,1);
        fmpq_div_2exp(tiny.value(),tiny.value(),80);
        require(midpoint_is(acb_realref(known.value()),expected) &&
                    midpoint_is(acb_imagref(known.value()),tiny),
                "a small part known exactly was dropped");
    }
    // Zero stays zero; a non-finite centre is left alone.
    {
        d::ComplexBall ball;
        d::round_center(ball.value(),bits);
        require(acb_is_zero(ball.value()),"a zero centre moved");
        arb_pos_inf(acb_realref(ball.value()));
        d::round_center(ball.value(),bits);
        require(!acb_is_finite(ball.value()),"a non-finite centre was changed");
    }
}

// The exact centering of a complex multiple of a real polynomial costs what
// the real polynomial's does. The centroid of c p is p's, and real; computed
// from ball coefficients its imaginary part is rounding noise, which a
// rounding of each part to its own 8 bits kept as a centre with a 70-bit
// denominator: the shifted coefficients then carried degree x 70 bits more.
struct CenteringEvent { std::string reason; bool accepted; uint64_t output_bits; };
static std::vector<CenteringEvent> centering_events(certpolysolve::detail::DensePolynomial input) {
    const auto degree=input.coefficients.size()-1;
    auto policy=d::maintained_method_policy(CPS_METHOD_SECULAR);
    std::vector<CenteringEvent> seen;
    policy.trace_data=&seen;
    policy.trace=[](void *data,const d::TraceEvent &e) {
        if (std::string_view(e.stage)=="exact-centering")
            static_cast<std::vector<CenteringEvent> *>(data)->push_back({e.reason,e.preparation_accepted,e.exact_output_bits});
    };
    d::SolveState state(std::move(input),policy);
    cps_options options;
    options.method=CPS_METHOD_SECULAR;options.goal=CPS_ROOT_APPROXIMATION;
    // One solve with the units the one-unit calls of a session used to share (one advance a solve state since 18 October 2026).
    options.bits=53;options.max_work_bits=4096;options.work_units=1000;
    cps_result result=state.advance(&options);
    require(result.goal_met,"scaled centering fixture failed to finish");
    require(d::disjoint(result.roots),"scaled centering fixture boxes overlap");
    uint64_t total=0;
    for (const auto &r:result.roots) total+=r.count;
    require(total==degree,"scaled centering fixture lost roots");
    return seen;
}
static void complex_multiple_centering() {
    // Roots on two circles about 2^-b (the fixture of scale_invariant_centering),
    // times an integer of 96 bits so that the leading coefficients do not fit
    // the 64 bits the centroid is computed with.
    const slong m=32,b=8;
    d::IntegerPolynomial shifted_linear,power,inner,outer,circles;
    fmpz_poly_set_coeff_si(shifted_linear.value(),0,-2);
    fmpz_poly_set_coeff_si(shifted_linear.value(),1,WORD(1)<<(b+1));
    fmpz_poly_pow(power.value(),shifted_linear.value(),static_cast<ulong>(m));
    d::Integer four_power,scale;
    fmpz_ui_pow_ui(four_power.value(),4,static_cast<ulong>(m));
    fmpz_poly_set(outer.value(),power.value());
    fmpz_sub(outer.value()->coeffs,outer.value()->coeffs,four_power.value());
    fmpz_poly_set(inner.value(),power.value());
    fmpz_sub_ui(inner.value()->coeffs,inner.value()->coeffs,1);
    fmpz_poly_mul(circles.value(),outer.value(),inner.value());
    fmpz_ui_pow_ui(scale.value(),3,60);
    fmpz_poly_scalar_mul_fmpz(circles.value(),circles.value(),scale.value());
    const auto length=static_cast<size_t>(fmpz_poly_length(circles.value()));
    certpolysolve::detail::DensePolynomial real_input,complex_input;
    real_input.coefficients.resize(length);
    complex_input.coefficients.resize(length);
    d::Integer part;
    for (size_t k=0;k<length;++k) {
        const auto coefficient=circles.value()->coeffs+static_cast<slong>(k);
        fmpq_set_fmpz(real_input.coefficients[k].real.significand.value(),coefficient);
        fmpz_mul_ui(part.value(),coefficient,3);
        fmpq_set_fmpz(complex_input.coefficients[k].real.significand.value(),part.value());
        fmpz_mul_ui(part.value(),coefficient,4);
        fmpq_set_fmpz(complex_input.coefficients[k].imag.significand.value(),part.value());
    }
    const auto selected=[](const std::vector<CenteringEvent> &seen) {
        const auto found=std::find_if(seen.begin(),seen.end(),[](const CenteringEvent &e) {
            return e.reason=="centering-selected" && e.accepted;
        });
        require(found!=seen.end(),"the offset two-circle cloud was not centered");
        return found->output_bits;
    };
    const uint64_t real_bits=selected(centering_events(std::move(real_input)));
    const uint64_t complex_bits=selected(centering_events(std::move(complex_input)));
    require(real_bits>0 && complex_bits>0,"the centering reported no coefficient sizes");
    // Each shifted coefficient s becomes 3s and 4s: two parts, two and three bits more.
    require(complex_bits<=2*real_bits+5*length,
            "the exact centering of a complex multiple grew beyond the real polynomial's");
}

static bool reason_seen(const std::vector<std::pair<std::string,bool>> &seen,std::string_view reason,bool accepted) {
    return std::any_of(seen.begin(),seen.end(),[&](const auto &e) {return e.first==reason && e.second==accepted;});
}

static void scale_invariant_centering() {
    const slong n=65;
    const slong required_gain=static_cast<slong>(FLINT_BIT_COUNT(static_cast<ulong>(n)));
    // Flat q with q[n-1]=0 (so the mean root, the center, is exactly 2^-a), and
    // p(x)=q(2^a x-1). Coefficient n-j of p over the leading one is at most
    // C(n+1,j) 2^-aj <= ((n+1)/2^a)^j < 1 once 2^a > n+1: no plain height.
    const slong a=static_cast<slong>(FLINT_BIT_COUNT(static_cast<ulong>(n+1)));
    d::IntegerPolynomial q,linear,cloud,translate;
    for (slong k=0;k<=n;++k) {
        // Rudin-Shapiro signs: +-1 by the parity of adjacent 1-bit pairs in k.
        const int pairs=std::popcount(static_cast<ulong>(k)&(static_cast<ulong>(k)>>1));
        if (k!=n-1) fmpz_poly_set_coeff_si(q.value(),k,pairs%2 ? -1 : 1);
    }
    fmpz_poly_set_coeff_si(linear.value(),0,-1);
    fmpz_poly_set_coeff_si(linear.value(),1,WORD(1)<<a);
    fmpz_poly_compose(cloud.value(),q.value(),linear.value());
    require(!fmpz_is_zero(cloud.value()->coeffs),"cloud fixture has a root at zero");
    require(!d::center_height_allows_gain(cloud.value()),"cloud fixture has plain monic height");
    // The exact translate p(x+2^-a)=q(2^a x).
    fmpz_poly_set(translate.value(),q.value());
    for (slong k=0;k<=n;++k) fmpz_mul_2exp(translate.value()->coeffs+k,translate.value()->coeffs+k,static_cast<ulong>(a*k));
    const auto cloud_profile=d::height_profile(cloud.value()),translate_profile=d::height_profile(translate.value());
    const double scale=d::root_scale(cloud_profile);
    require(d::height_at(translate_profile,scale)<required_gain &&
                d::height_at(cloud_profile,scale)-d::height_at(translate_profile,scale)>=required_gain,
            "cloud translate is not flat at the original root scale");
    require(reason_seen(secular_centering(cloud.value(),true),"centering-selected",true),
            "offset flat cloud was not centered at the original root scale");
    require(reason_seen(secular_centering(cloud.value(),false),"centering-height-gain-impossible",false),
            "plain rule unexpectedly admitted the offset flat cloud");

    // Roots on two circles about 2^-b, of radii 2^(1-b) and 2^(-1-b):
    // with z=2^(b+1) x-2, p=(z^m-4^m)(z^m-1). At the original root scale the
    // translate is not flat (the circles straddle it: about m/2 bits remain),
    // yet it gains far more there and loses at no scale. A flatness test
    // would refuse it.
    const slong m=32,b=8;
    d::IntegerPolynomial shifted_linear,power,inner,outer,circles;
    fmpz_poly_set_coeff_si(shifted_linear.value(),0,-2);
    fmpz_poly_set_coeff_si(shifted_linear.value(),1,WORD(1)<<(b+1));
    fmpz_poly_pow(power.value(),shifted_linear.value(),static_cast<ulong>(m));
    d::Integer four_power;
    fmpz_ui_pow_ui(four_power.value(),4,static_cast<ulong>(m));
    fmpz_poly_set(outer.value(),power.value());
    fmpz_sub(outer.value()->coeffs,outer.value()->coeffs,four_power.value());
    fmpz_poly_set(inner.value(),power.value());
    fmpz_sub_ui(inner.value()->coeffs,inner.value()->coeffs,1);
    fmpz_poly_mul(circles.value(),outer.value(),inner.value());
    const slong circles_gain=static_cast<slong>(FLINT_BIT_COUNT(static_cast<ulong>(2*m)));
    require(!d::center_height_allows_gain(circles.value()),"two-circle fixture has plain monic height");
    {
        // The exact translate: z^m-4^m and z^m-1 with z=2^(b+1) y.
        d::IntegerPolynomial zpower,outer_shift,inner_shift,shifted;
        fmpz_poly_set_coeff_si(zpower.value(),m,1);
        fmpz_poly_set(outer_shift.value(),zpower.value());
        fmpz_sub(outer_shift.value()->coeffs,outer_shift.value()->coeffs,four_power.value());
        fmpz_poly_set(inner_shift.value(),zpower.value());
        fmpz_sub_ui(inner_shift.value()->coeffs,inner_shift.value()->coeffs,1);
        fmpz_poly_mul(shifted.value(),outer_shift.value(),inner_shift.value());
        for (slong k=0;k<=2*m;++k)
            fmpz_mul_2exp(shifted.value()->coeffs+k,shifted.value()->coeffs+k,static_cast<ulong>((b+1)*k));
        const auto before=d::height_profile(circles.value()),after=d::height_profile(shifted.value());
        const double circles_scale=d::root_scale(before);
        require(d::height_at(after,circles_scale)>=circles_gain,"two-circle translate is flat");
        require(d::height_at(before,circles_scale)-d::height_at(after,circles_scale)>=circles_gain &&
                    d::worst_height_excess(after,before)<circles_gain,
                "two-circle translate does not dominate");
    }
    require(reason_seen(secular_centering(circles.value(),true),"centering-selected",true),
            "offset two-circle cloud was not centered");

    // n! L_n(x) for Laguerre degree 64: centering at its mean root removes
    // height at the original root scale, but adds more at small scales, where
    // |L(64)| >> |L(0)|; the plain rule refuses (the plain height grows), and
    // so must the scale-invariant route.
    const slong degree=64;
    d::IntegerPolynomial laguerre;
    d::Integer coefficient;
    for (slong k=0;k<=degree;++k) {
        fmpz_bin_uiui(coefficient.value(),static_cast<ulong>(degree),static_cast<ulong>(k));
        for (slong j=k+1;j<=degree;++j) fmpz_mul_ui(coefficient.value(),coefficient.value(),static_cast<ulong>(j));
        if (k%2) fmpz_neg(coefficient.value(),coefficient.value());
        fmpz_poly_set_coeff_fmpz(laguerre.value(),k,coefficient.value());
    }
    const auto seen=secular_centering(laguerre.value(),true);
    require(!reason_seen(seen,"centering-selected",true),"partial Laguerre centering was accepted");
    require(reason_seen(seen,"centering-no-height-improvement",false),"Laguerre centering was not considered");

    // Laguerre centered at its mean root (so its x^(n-1) coefficient vanishes),
    // scaled by 2^shift and given x^(n-1) coefficient 1: the new center is
    // below 2^-shift, far too small for any route to gain, so the exact shift
    // is skipped. Without the flag the plain rule declines it too, since 12 October 2026 before the
    // exact shift is built: the translate's constant term shows that its height cannot fall.
    const slong shift=64;
    d::IntegerPolynomial tiny;
    d::Integer mean;
    fmpz_set_si(mean.value(),degree);
    fmpz_poly_taylor_shift(tiny.value(),laguerre.value(),mean.value());
    require(fmpz_is_zero(tiny.value()->coeffs+degree-1),"Laguerre mean root is not its center");
    fmpz_poly_scalar_mul_2exp(tiny.value(),tiny.value(),static_cast<ulong>(shift));
    fmpz_poly_set_coeff_si(tiny.value(),degree-1,1);
    require(reason_seen(secular_centering(tiny.value(),true),"centering-shift-too-small",false),
            "a translation too small to gain was built");
    require(reason_seen(secular_centering(tiny.value(),false),"centering-constant-term-bound",false),
            "the plain rule did not refuse the small translation before building it");
}

// The centroid of a polynomial with one distant root lies far from every root. The translate there
// has enormous coefficients and no rule can accept it; it must be refused before it is built (found
// on a monic polynomial of degree 512 with 65,536-bit coefficients, 12 October 2026). The bound that
// predicts the refusal must never exceed the height of the translate actually built.
static void distant_root_centering() {
    const slong degree=96, far_bits=4096;
    flint_rand_t state;
    flint_rand_init(state);
    d::IntegerPolynomial p;
    d::Integer coefficient;
    for (slong k=0;k<degree;++k) {
        fmpz_randbits(coefficient.value(),state,far_bits);
        if (fmpz_is_zero(coefficient.value())) fmpz_one(coefficient.value());
        fmpz_poly_set_coeff_fmpz(p.value(),k,coefficient.value());
    }
    fmpz_poly_set_coeff_si(p.value(),degree,1);
    for (const bool scale_invariant : {false,true}) {
        const auto seen=secular_centering(p.value(),scale_invariant);
        require(reason_seen(seen,"centering-constant-term-bound",false),"the translate to a distant centroid was built");
        require(!reason_seen(seen,"centering-selected",true),"the translate to a distant centroid was taken");
    }
    const auto height=[](const fmpz_poly_t q) {
        return FLINT_ABS(fmpz_poly_max_bits(q))-static_cast<slong>(fmpz_bits(fmpz_poly_lead(q)));
    };
    d::IntegerPolynomial q,shifted;
    d::Rational center;
    slong bounded=0;
    for (int trial=0;trial<400;++trial) {
        fmpz_poly_randtest_not_zero(q.value(),state,2+static_cast<slong>(n_randint(state,40)),
                                    1+static_cast<flint_bitcnt_t>(n_randint(state,300)));
        if (fmpz_poly_degree(q.value())<1) continue;
        fmpq_randtest(center.value(),state,1+static_cast<flint_bitcnt_t>(n_randint(state,80)));
        const auto bound=d::translated_height_lower_bound(q.value(),center.value());
        if (!bound || !d::center_integer(shifted.value(),q.value(),center.value(),UINT64_MAX,nullptr)) continue;
        ++bounded;
        require(height(shifted.value())>=*bound,"the predicted height exceeds the built translate's");
    }
    require(bounded>=200,"too few translates were bounded to test the prediction");
    flint_rand_clear(state);
}

// A translate that carries roots out of the binary64 stage's reach with no gap in scale between them and
// the others (EnginePolicy::centering_keeps_origin_scales, unseparated_origin_roots). Forty-eight roots 2^-k
// (k = 0..47) and 144 on the unit circle about 2: the mean root is about 3/2 and the height rules take the
// translate, which would carry the sixteen roots below 2^-32 to one binary64 point of the new coordinates,
// with nothing for a cluster restart to hold: refused. A group of sixteen roots on one circle of radius 2^-40,
// far below the rest, is carried too, but as one certified group: it stays centred, and so does a cloud with
// no root near the origin.
static void origin_scales_centering() {
    const slong small=48,cloud=144,reach=d::SecularMachineWorkspace::stage_bits(53)-53;
    d::IntegerPolynomial p,factor,linear,power,offset;
    fmpz_poly_one(p.value());
    for (slong k=0;k<small;++k) {
        d::Integer lead;
        fmpz_one(lead.value());
        fmpz_mul_2exp(lead.value(),lead.value(),static_cast<ulong>(k));
        fmpz_poly_zero(factor.value());
        fmpz_poly_set_coeff_si(factor.value(),0,-1);
        fmpz_poly_set_coeff_fmpz(factor.value(),1,lead.value());
        fmpz_poly_mul(p.value(),p.value(),factor.value());
    }
    fmpz_poly_set_coeff_si(linear.value(),0,-2);
    fmpz_poly_set_coeff_si(linear.value(),1,1);
    fmpz_poly_pow(power.value(),linear.value(),static_cast<ulong>(cloud));
    fmpz_poly_set(offset.value(),power.value());
    fmpz_sub_ui(offset.value()->coeffs,offset.value()->coeffs,1);
    fmpz_poly_mul(p.value(),p.value(),offset.value());
    const slong carried=d::unseparated_origin_roots(d::height_profile(p.value()),std::log2(1.5),static_cast<double>(reach));
    require(carried>=small-reach-4 && carried<=small-reach+4,"the roots carried out of binary64's reach were not counted");
    require(reason_seen(secular_centering(p.value(),true,false),"centering-selected",true),
            "the fixture's translate is not one the height rules take");
    const auto refused=secular_centering(p.value(),true);
    require(reason_seen(refused,"centering-origin-scales",false),"the translate over the origin's scales was not refused");
    require(!reason_seen(refused,"centering-selected",true),"the translate over the origin's scales was taken");
    // One certified group: sixteen roots of modulus 2^-40 beside 48 about 2.
    d::IntegerPolynomial group,grouped;
    fmpz_poly_set_coeff_si(group.value(),0,-1);
    d::Integer lead;
    fmpz_one(lead.value());
    fmpz_mul_2exp(lead.value(),lead.value(),16*40);
    fmpz_poly_set_coeff_fmpz(group.value(),16,lead.value());
    fmpz_poly_pow(power.value(),linear.value(),48);
    fmpz_poly_set(offset.value(),power.value());
    fmpz_sub_ui(offset.value()->coeffs,offset.value()->coeffs,1);
    fmpz_poly_mul(grouped.value(),group.value(),offset.value());
    require(d::unseparated_origin_roots(d::height_profile(grouped.value()),std::log2(1.5),static_cast<double>(reach))==0,
            "a certified origin group counted as unseparated");
    require(d::unseparated_origin_roots(d::height_profile(offset.value()),std::log2(2.0),static_cast<double>(reach))==0,
            "an offset cloud has roots at the origin");
}

// Option 2 for a root group at the original origin (5 October 2026): the centred source evaluates the
// uncentred polynomial where its terms are smaller, so the group keeps the accuracy its goals (relative to
// the original origin) ask at the precision of its own spacing. The fixture: 16 roots of modulus 2^-40
// beside 48 about 2, centred at about 3/2 by the height rules.
static std::pair<uint64_t,bool> origin_cluster_solve(bool origin_basis) {
    const slong group=16,cloud=48,exponent=40;
    d::IntegerPolynomial small,linear,power,p;
    fmpz_poly_set_coeff_si(small.value(),0,-1);
    d::Integer lead;
    fmpz_one(lead.value());
    fmpz_mul_2exp(lead.value(),lead.value(),static_cast<ulong>(group*exponent));
    fmpz_poly_set_coeff_fmpz(small.value(),group,lead.value());
    fmpz_poly_set_coeff_si(linear.value(),0,-2);
    fmpz_poly_set_coeff_si(linear.value(),1,1);
    fmpz_poly_pow(power.value(),linear.value(),static_cast<ulong>(cloud));
    fmpz_sub_ui(power.value()->coeffs,power.value()->coeffs,1);
    fmpz_poly_mul(p.value(),small.value(),power.value());
    certpolysolve::detail::DensePolynomial input;
    input.coefficients.resize(static_cast<size_t>(fmpz_poly_length(p.value())));
    for (slong k=0;k<fmpz_poly_length(p.value());++k)
        fmpq_set_fmpz(input.coefficients[static_cast<size_t>(k)].real.significand.value(),p.value()->coeffs+k);
    auto policy=d::maintained_method_policy(CPS_METHOD_SECULAR);
    require(policy.centred_origin_basis,"maintained secular omitted the origin basis");
    policy.centred_origin_basis=origin_basis;
    bool centred=false;
    policy.trace_data=&centred;
    policy.trace=[](void *data,const d::TraceEvent &e) {
        if (std::string_view(e.stage)=="exact-centering" && std::string_view(e.reason)=="centering-selected" &&
            e.preparation_accepted)
            *static_cast<bool *>(data)=true;
    };
    d::SolveState state(std::move(input),policy);
    cps_options options;
    options.method=CPS_METHOD_SECULAR;options.goal=CPS_ROOT_APPROXIMATION;
    // One solve with the units the one-unit calls of a session used to share (one advance a solve state since 18 October 2026).
    options.bits=53;options.max_work_bits=65536;options.work_units=5000;
    cps_result result=state.advance(&options);
    require(result.goal_met,"origin cluster fixture failed to finish");
    uint64_t total=0;
    for (const auto &r:result.roots) total+=r.count;
    require(total==static_cast<uint64_t>(fmpz_poly_degree(p.value())),"origin cluster fixture lost roots");
    return {result.work_bits,centred};
}
// A cluster at the original origin narrower than binary64 resolves at the centre's distance (fresh
// clustered-719's 180 roots 2^-94 apart about 0, centred at about -1/2): 16 real roots k 2^-90 beside 48
// about 1, centred at 3/4. The restart's model about the cluster's mean has that mean only to about 2^-53
// of its circle; about the original origin itself it is exact.
static std::pair<bool,uint64_t> origin_cluster_restart(bool origin_basis,slong group=16,slong cloud=48) {
    const slong spacing=90,radius=4;
    d::IntegerPolynomial cluster,factor,power,ring,p;
    fmpz_poly_one(cluster.value());
    for (slong k=1;k<=group;++k) {
        fmpz_poly_zero(factor.value());
        d::Integer scale;
        fmpz_one(scale.value());
        fmpz_mul_2exp(scale.value(),scale.value(),static_cast<ulong>(spacing));
        fmpz_poly_set_coeff_fmpz(factor.value(),1,scale.value());
        fmpz_poly_set_coeff_si(factor.value(),0,-k);
        fmpz_poly_mul(cluster.value(),cluster.value(),factor.value());
    }
    fmpz_poly_zero(factor.value());
    fmpz_poly_set_coeff_si(factor.value(),0,-1);
    fmpz_poly_set_coeff_si(factor.value(),1,1);
    fmpz_poly_pow(power.value(),factor.value(),static_cast<ulong>(cloud));
    fmpz_poly_scalar_mul_2exp(ring.value(),power.value(),static_cast<ulong>(cloud*radius));
    fmpz_sub_ui(ring.value()->coeffs,ring.value()->coeffs,1);
    fmpz_poly_mul(p.value(),cluster.value(),ring.value());
    certpolysolve::detail::DensePolynomial input;
    input.coefficients.resize(static_cast<size_t>(fmpz_poly_length(p.value())));
    for (slong k=0;k<fmpz_poly_length(p.value());++k)
        fmpq_set_fmpz(input.coefficients[static_cast<size_t>(k)].real.significand.value(),p.value()->coeffs+k);
    auto policy=d::maintained_method_policy(CPS_METHOD_SECULAR);
    policy.centred_origin_basis=origin_basis;
    d::SolveState state(std::move(input),policy);
    cps_options options;
    options.method=CPS_METHOD_SECULAR;options.goal=CPS_ROOT_APPROXIMATION;
    // Both policies get enough work to finish, including their nested models.
    options.bits=53;options.max_work_bits=65536;options.work_units=1000;
    options.work_units_default=false;
    cps_result result=state.advance(&options);
    const auto calls=result.progress.work_units_this_call;
    if (result.goal_met) {
        uint64_t total=0;
        for (const auto &r:result.roots) total+=r.count;
        require(total==static_cast<uint64_t>(fmpz_poly_degree(p.value())),"origin cluster restart fixture lost roots");
    }
    return {result.goal_met,calls};
}
// The coordinate rule at a high goal: x^383 - 2 (2^128 x^2 - 1)^6 (the matrix's mignotte_generalized_383_64_6),
// centred, at 3,322 bits. Its nodes iterate at 170 bits while the precision growth raises them toward the goal;
// the rule compared the goal with their bits and raised every node, near the original origin or not, to the
// goal's 3,354 bits at once: 4 s became over 300 s. Only a node already at the goal's stage takes the excess.
static void origin_rule_at_a_high_goal() {
    d::IntegerPolynomial pair,power,p;
    d::Integer square;
    fmpz_one(square.value());
    fmpz_mul_2exp(square.value(),square.value(),128);
    fmpz_poly_set_coeff_fmpz(pair.value(),2,square.value());
    fmpz_poly_set_coeff_si(pair.value(),0,-1);
    fmpz_poly_pow(power.value(),pair.value(),6);
    fmpz_poly_scalar_mul_si(p.value(),power.value(),-2);
    d::Integer one;
    fmpz_one(one.value());
    d::IntegerPolynomial monomial;
    fmpz_poly_set_coeff_fmpz(monomial.value(),383,one.value());
    fmpz_poly_add(p.value(),p.value(),monomial.value());
    certpolysolve::detail::DensePolynomial input;
    input.coefficients.resize(static_cast<size_t>(fmpz_poly_length(p.value())));
    for (slong k=0;k<fmpz_poly_length(p.value());++k)
        fmpq_set_fmpz(input.coefficients[static_cast<size_t>(k)].real.significand.value(),p.value()->coeffs+k);
    d::SolveState state(std::move(input),d::maintained_method_policy(CPS_METHOD_SECULAR));
    cps_options options;
    options.method=CPS_METHOD_SECULAR;options.goal=CPS_ROOT_APPROXIMATION;
    // One solve with the units the one-unit calls of a session used to share (one advance a solve state since 18 October 2026).
    options.bits=3322;options.max_work_bits=65536;options.work_units=60;
    cps_result result=state.advance(&options);
    require(result.goal_met,"the coordinate rule raised every node to a high goal at once");
}
static void origin_basis() {
    // Evaluated through the translate, the group's 16 roots of modulus 2^-40 cost the cancellation of an
    // expansion about 3/2: about 16 x 40 bits. Through the uncentred polynomial they cost their own spacing.
    const auto [with,centred]=origin_cluster_solve(true);
    const auto [without,also_centred]=origin_cluster_solve(false);
    require(centred && also_centred,"the origin cluster fixture was not centred");
    const uint64_t translate_cost=16*40;
    require(without>=translate_cost,"the fixture no longer pays for the translate");
    require(with<translate_cost,"the origin group paid for the translate's cancellation");
    // 48 roots k 2^-90 at the origin beside 80 about 1. Compare the maintained
    // origin basis with the disabled-policy control, charging nested work in
    // both. The old cutoff of 30 units had excluded that work and could stop
    // the working policy before its goal; it did not directly test the benefit.
    const auto [met,calls]=origin_cluster_restart(true,48,80);
    const auto [control_met,control_calls]=origin_cluster_restart(false,48,80);
    require(met && control_met,"an origin-basis comparison did not finish");
    require(calls<control_calls,"the origin basis no longer reduces work on the origin cluster");
}

int main() {
    certpolysolve::ThreadScope lifetime;
    try {
        height();exact_shift_identities();refusal_reason();failed_cover_retains_only_proposals();
        translated_roots(CPS_METHOD_CLASSICAL);
        translated_roots(CPS_METHOD_SECULAR);
        scale_invariant_centering();
        distant_root_centering();
        origin_scales_centering();
        origin_basis();
        origin_rule_at_a_high_goal();
        center_rounding();complex_multiple_centering();
        std::cout<<"coordinate policy analytic/count/resource checks passed\n";
    } catch (const std::exception &e) {std::cerr<<e.what()<<'\n';return 1;}
}
