#pragma once
#include "methods/piecewise/piecewise_roots.hpp"

namespace piecewise_performance {
namespace d=certpolysolve::detail;
inline void upper(arb_t out,const acb_t z,slong bits) {
    d::Float value;acb_get_abs_ubound_arf(value.value(),z,bits);arb_set_arf(out,value.value());
}
inline void lower(arb_t out,const acb_t z,slong bits) {
    d::Float value;acb_get_abs_lbound_arf(value.value(),z,bits);arb_set_arf(out,value.value());
}
struct Rectangle {
    d::Root root;
    d::ComplexBall center,local,local_box;
    d::RealBall inner,outer;
};
inline std::optional<Rectangle> rectangle(const d::RootModel &m,const d::ComplexBall &point,
    const cps_options &o,slong bits,const d::RealBall *proposed_radius=nullptr) {
    if(!acb_is_finite(point.value()) || !arb_is_positive(m.radius.value()))return {};
    d::RealBall point_norm,unit;arb_one(unit.value());upper(point_norm.value(),point.value(),bits);
    if(!arb_lt(point_norm.value(),unit.value()))return {};
    Rectangle r;
    acb_mul_arb(r.center.value(),point.value(),m.radius.value(),bits);
    acb_add(r.center.value(),r.center.value(),m.center.value(),bits);
    acb_get_mid(r.center.value(),r.center.value());
    if(!acb_is_finite(r.center.value()) || acb_is_zero(r.center.value()))return {};
    acb_sub(r.local.value(),r.center.value(),m.center.value(),bits);
    acb_div_arb(r.local.value(),r.local.value(),m.radius.value(),bits);
    d::RealBall radius,one,distance,endpoint;
    upper(radius.value(),r.center.value(),bits);
    arb_mul_2exp_si(radius.value(),radius.value(),-static_cast<slong>(o.bits)-5);
    if(!o.absolute_tolerance.zero()) {
        d::RealBall tolerance;o.absolute_tolerance.enclose(tolerance.value(),bits);
        arb_mul_2exp_si(tolerance.value(),tolerance.value(),-3);
        if(arb_gt(tolerance.value(),radius.value()))radius=tolerance;
    }
    if(proposed_radius)radius=*proposed_radius;
    if(!arb_is_positive(radius.value()))return {};
    r.root.box=r.center;
    arb_add_error(acb_realref(r.root.box.value()),radius.value());
    arb_add_error(acb_imagref(r.root.box.value()),radius.value());
    r.root.count=r.root.multiplicity=1;r.root.realness=CPS_REAL_UNKNOWN;
    r.root.publish(bits);
    if(!d::meets_goal(r.root,o,bits))return {};
    acb_sub(r.local_box.value(),r.root.box.value(),m.center.value(),bits);
    acb_div_arb(r.local_box.value(),r.local_box.value(),m.radius.value(),bits);
    arb_one(one.value());upper(distance.value(),r.local_box.value(),bits);
    if(!arb_lt(distance.value(),one.value()) || (m.lower_index && acb_contains_zero(r.root.box.value())))return {};
    for(size_t k=0;k<4;++k) {
        const auto &e=r.root.endpoints[k];
        d::Number(e.mantissa_text()+"@"+e.exponent_text()).enclose(endpoint.value(),bits);
        arb_sub(distance.value(),endpoint.value(),k<2?acb_realref(r.center.value()):acb_imagref(r.center.value()),bits);
        if(!(k%2))arb_neg(distance.value(),distance.value());
        if(!arb_is_positive(distance.value()))return {};
        if(!k)r.inner=distance;else arb_min(r.inner.value(),r.inner.value(),distance.value(),bits);
    }
    d::ComplexBall displacement;
    acb_sub(displacement.value(),r.root.box.value(),r.center.value(),bits);
    upper(r.outer.value(),displacement.value(),bits);
    arb_div(r.inner.value(),r.inner.value(),m.radius.value(),bits);
    arb_div(r.outer.value(),r.outer.value(),m.radius.value(),bits);
    return r;
}

// q'' on the whole convex rectangle. Interval Horner retains coefficient
// cancellation; the previous integral majorant instead sums coefficient norms.
inline void second_derivative(acb_t out,const d::ComplexPolynomial &q,const acb_t box,slong bits) {
    d::ComplexBall term;acb_zero(out);
    for(slong j=acb_poly_degree(q.value());j>=2;--j) {
        acb_mul(out,out,box,bits);
        acb_mul_si(term.value(),q.value()->coeffs+j,j,bits);
        acb_mul_si(term.value(),term.value(),j-1,bits);
        acb_add(out,out,term.value(),bits);
    }
}
inline bool derivative_rectangle(d::Root &root,const d::RootModel &model,const d::ComplexBall &point,
    const cps_options &o,slong bits) {
    auto box=rectangle(model,point,o,bits);if(!box)return false;
    d::ComplexBall value,derivative,second;
    acb_poly_evaluate2(value.value(),derivative.value(),model.polynomial.value(),box->local.value(),bits);
    d::RealBall slope,constant,rounding,imaginary,halfwidth;
    lower(slope.value(),derivative.value(),bits);upper(constant.value(),value.value(),bits);
    if(!arb_is_positive(slope.value()))return false;
    arb_get_rad_arb(rounding.value(),acb_realref(value.value()));
    arb_get_rad_arb(imaginary.value(),acb_imagref(value.value()));arb_add(rounding.value(),rounding.value(),imaginary.value(),bits);
    for(slong shrink=0;shrink<=bits;++shrink) {
        second_derivative(second.value(),model.polynomial,box->local_box.value(),bits);
        d::RealBall lhs,rhs=constant,term,power;
        arb_mul(lhs.value(),slope.value(),box->inner.value(),bits);
        upper(term.value(),second.value(),bits);
        arb_sqr(power.value(),box->outer.value(),bits);
        arb_mul(term.value(),term.value(),power.value(),bits);arb_mul_2exp_si(term.value(),term.value(),-1);
        arb_add(rhs.value(),rhs.value(),term.value(),bits);
        const auto error=d::root_model_error(model,box->root.box.value(),bits);
        arb_add(rhs.value(),rhs.value(),error.value(),bits);
        if(arb_gt(lhs.value(),rhs.value())) {
            box->root.certificate="piecewise linear Rouche with interval second derivative on published rectangle";
            root=std::move(box->root);return true;
        }
        // Match the maintained proposal stopping policy. These checks do NOT
        // assert impossibility; only the strict Rouche inequality proves roots.
        arb_mul_2exp_si(lhs.value(),lhs.value(),-2);
        if((!arb_is_zero(error.value()) && !arb_lt(error.value(),lhs.value())) ||
           (!arb_is_zero(rounding.value()) && !arb_lt(rounding.value(),lhs.value())) ||
           !arb_lt(constant.value(),lhs.value()))return false;
        arb_get_rad_arb(halfwidth.value(),acb_realref(box->root.box.value()));
        arb_mul_2exp_si(halfwidth.value(),halfwidth.value(),-1);
        box=rectangle(model,point,o,bits,&halfwidth);if(!box)return false;
    }
    return false;
}

// A sufficient rejection of every translated-monomial count predicate, not
// the old post-failure linear/4 classification. The unchanged Taylor floor
// persists under shrinking; R*qhat'(|w|+R) bounds every nonconstant jet term.
inline bool immutable_tail_rejects(const d::RootModel &model,const d::ComplexBall &point,
    const cps_options &o,slong bits) {
    if(model.taylor_error_order > 0 || !arb_is_positive(model.taylor_model_error.value()) ||
       !arb_ge(model.error.value(),model.taylor_model_error.value()))return false;
    auto box=rectangle(model,point,o,bits);if(!box)return false;
    d::RealBall radial;
    upper(radial.value(),box->local.value(),bits);arb_add(radial.value(),radial.value(),box->outer.value(),bits);
    d::Magnitude r,coefficient,derivative,term,outer;
    arb_get_mag(r.value(),radial.value());arb_get_mag(outer.value(),box->outer.value());
    for(slong j=acb_poly_degree(model.polynomial.value());j>=1;--j) {
        acb_get_mag(coefficient.value(),model.polynomial.value()->coeffs+j);
        mag_mul_ui(term.value(),coefficient.value(),static_cast<ulong>(j));
        mag_mul(derivative.value(),derivative.value(),r.value());mag_add(derivative.value(),derivative.value(),term.value());
    }
    mag_mul(derivative.value(),derivative.value(),outer.value());
    d::Float bound;arf_set_mag(bound.value(),derivative.value());
    d::RealBall variation;arb_set_arf(variation.value(),bound.value());
    return arb_ge(model.taylor_model_error.value(),variation.value());
}

// A natural isolating disc, with an explicit proof of the exported rectangle.
// For F=q+E, Cauchy gives |E'|<=epsilon/delta and
// |E''|<=2*epsilon/delta^2. Rouche proves one root in both r and 3r;
// the published rectangle is checked to lie between those two discs.
inline bool derivative_disc(d::Root &root,const d::RootModel &model,const d::ComplexBall &point,
    const cps_options &o,slong bits) {
    if(!acb_is_finite(point.value()) || !arb_is_positive(model.radius.value()))return false;
    d::RealBall point_norm,unit;arb_one(unit.value());upper(point_norm.value(),point.value(),bits);
    if(!arb_lt(point_norm.value(),unit.value()))return false;
    d::RealBall center_lower;
    lower(center_lower.value(),model.center.value(),bits);
    if(model.lower_index && !arb_gt(center_lower.value(),model.radius.value()))return false;
    d::ComplexBall center,local,value,derivative;
    acb_mul_arb(center.value(),point.value(),model.radius.value(),bits);
    acb_add(center.value(),center.value(),model.center.value(),bits);acb_get_mid(center.value(),center.value());
    if(!acb_is_finite(center.value()) || acb_is_zero(center.value()))return false;
    acb_sub(local.value(),center.value(),model.center.value(),bits);
    acb_div_arb(local.value(),local.value(),model.radius.value(),bits);
    d::RealBall norm,delta,one,residual,slope,term,radius;
    upper(norm.value(),local.value(),bits);arb_one(one.value());
    arb_sub(delta.value(),one.value(),norm.value(),bits);
    if(!arb_is_positive(delta.value()))return false;
    acb_poly_evaluate2(value.value(),derivative.value(),model.polynomial.value(),local.value(),bits);
    upper(residual.value(),value.value(),bits);arb_add(residual.value(),residual.value(),model.error.value(),bits);
    lower(slope.value(),derivative.value(),bits);
    arb_div(term.value(),model.error.value(),delta.value(),bits);arb_sub(slope.value(),slope.value(),term.value(),bits);
    if(!arb_is_positive(slope.value()))return false;
    if(arb_is_zero(residual.value())) {
        auto requested=rectangle(model,point,o,bits);if(!requested)return false;
        radius=requested->inner;arb_mul_2exp_si(radius.value(),radius.value(),-1);
    } else {
        arb_div(radius.value(),residual.value(),slope.value(),bits);
        arb_mul_2exp_si(radius.value(),radius.value(),1);
    }
    d::Float upper_radius;arb_get_ubound_arf(upper_radius.value(),radius.value(),bits);
    arb_set_arf(radius.value(),upper_radius.value());
    if(!arb_is_positive(radius.value()) || !arb_is_finite(radius.value()))return false;
    d::RealBall outer,second_bound,power;
    arb_mul_ui(outer.value(),radius.value(),3,bits);
    arb_sub(delta.value(),one.value(),norm.value(),bits);arb_sub(delta.value(),delta.value(),outer.value(),bits);
    if(!arb_is_positive(delta.value()))return false;
    d::ComplexBall box=local,second;
    arb_add_error(acb_realref(box.value()),outer.value());arb_add_error(acb_imagref(box.value()),outer.value());
    second_derivative(second.value(),model.polynomial,box.value(),bits);upper(second_bound.value(),second.value(),bits);
    arb_sqr(power.value(),delta.value(),bits);arb_div(term.value(),model.error.value(),power.value(),bits);
    arb_mul_2exp_si(term.value(),term.value(),1);arb_add(second_bound.value(),second_bound.value(),term.value(),bits);
    for(const auto *r:{&radius,&outer}) {
        d::RealBall lhs,rhs;
        arb_mul(lhs.value(),slope.value(),r->value(),bits);
        arb_sqr(power.value(),r->value(),bits);arb_mul(rhs.value(),second_bound.value(),power.value(),bits);
        arb_mul_2exp_si(rhs.value(),rhs.value(),-1);arb_add(rhs.value(),rhs.value(),residual.value(),bits);
        if(!arb_gt(lhs.value(),rhs.value()))return false;
    }
    d::RealBall halfwidth;
    arb_mul(halfwidth.value(),radius.value(),model.radius.value(),bits);arb_mul_2exp_si(halfwidth.value(),halfwidth.value(),1);
    auto published=rectangle(model,point,o,bits,&halfwidth);if(!published)return false;
    if(!arb_ge(published->inner.value(),radius.value()) || !arb_le(published->outer.value(),outer.value()))return false;
    published->root.certificate="piecewise natural isolating discs enclosing the published rectangle";
    root=std::move(published->root);return true;
}
} // namespace piecewise_performance
