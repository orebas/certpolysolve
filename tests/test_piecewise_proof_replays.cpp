#include "exact_winding_api.hpp"
#include "support/piecewise-performance/proof_experiments.hpp"
#include <iostream>

namespace d=certpolysolve::detail;
namespace experiment=piecewise_performance;
using namespace reference;

static d::RootModel model(const std::vector<Complex> &roots) {
    d::RootModel result;
    acb_set_si_si(result.center.value(),1,1);arb_set_d(result.radius.value(),0.5);
    const auto q=from_roots(roots);
    result.scale=-static_cast<slong>(roots.size());
    for(size_t j=0;j<q.size();++j) {
        d::ComplexBall c;
        d::Number(q[j].real.str()).enclose(acb_realref(c.value()),256);
        d::Number(q[j].imag.str()).enclose(acb_imagref(c.value()),256);
        acb_poly_set_coeff_acb(result.polynomial.value(),static_cast<slong>(j),c.value());
    }
    return result;
}
static d::ComplexBall ball(const Complex &z) {
    d::ComplexBall result;d::Number(z.real.str()).enclose(acb_realref(result.value()),256);
    d::Number(z.imag.str()).enclose(acb_imagref(result.value()),256);return result;
}
static Rational coordinate(const d::Endpoint &e) {return endpoint(certpolysolve::Dyadic{e.mantissa_text(),e.exponent_text()});}
static void check_count(const d::Root &r,const reference::ComplexPolynomial &p) {
    const auto count=rectangle_count(p,coordinate(r.endpoints[0]),coordinate(r.endpoints[1]),
                                      coordinate(r.endpoints[2]),coordinate(r.endpoints[3]));
    require(count && *count==static_cast<int>(r.count),"derivative certificate disagrees with exact winding");
}
static void check_deferred_guard(const d::RootModel &m,const d::ComplexBall &point,const cps_options &o) {
    d::Root full,deferred;
    d::ChartFailure full_failure=d::ChartFailure::Domain, deferred_failure=full_failure;
    const bool a=d::certify_root_model(full,m,point.value(),o,256,&full_failure);
    const bool b=d::certify_root_model(deferred,m,point.value(),o,256,&deferred_failure,
                                     nullptr,nullptr,nullptr,0,true);
    require(a==b && (a || full_failure==deferred_failure),
            "deferred guard changed the full-jet outcome or refusal reason");
    if(a) {
        require(full.count==deferred.count,"deferred jet changed root count");
        for(size_t k=0;k<4;++k)
            require(full.endpoints[k].mantissa_text()==deferred.endpoints[k].mantissa_text() &&
                    full.endpoints[k].exponent_text()==deferred.endpoints[k].exponent_text(),
                    "deferred guard changed the published rectangle");
    }
}
static void deferred_guard_boundaries() {
    cps_options options;options.goal=CPS_ROOT_ACCURACY;
    options.absolute_tolerance=d::Number("1/16");
    auto m=model({{0},{-1}});
    const auto refusal=[&](const d::RootModel &input,const d::ComplexBall &point,
                           const cps_options &settings,d::ChartFailure expected) {
        check_deferred_guard(input,point,settings);
        d::Root root;d::ChartFailure why=d::ChartFailure::Domain;
        require(!d::certify_root_model(root,input,point.value(),settings,256,&why) && why==expected,
                "deferred-jet control did not reach its intended refusal");
    };
    // Both floors exceed the linear allowance before higher derivatives help.
    refusal(m,ball({Rational(1)/8}),options,d::ChartFailure::Candidate);
    d::ComplexBall uncertain;
    arb_add_error_2exp_si(acb_realref(uncertain.value()),-3);
    acb_poly_set_coeff_acb(m.polynomial.value(),0,uncertain.value());
    refusal(m,ball({0}),options,d::ChartFailure::Rounding);

    // A valid but very broad requested tolerance leaves every attempted box
    // outside this chart. Deferring the jet must preserve the Width refusal.
    auto broad=options;broad.absolute_tolerance=d::Number(power2(300).str());
    refusal(model({{0},{-1}}),ball({0}),broad,d::ChartFailure::Width);

    // Higher derivatives matter here: the initial box contains both nearby
    // roots; the ordinary full jet eventually proves one after shrinking.
    const auto pair=model({{0},{Rational(-1)/64}});
    options.absolute_tolerance=d::Number("1/4");
    const auto point=ball({0});
    d::Root root;d::PiecewiseRootStatistics statistics;
    require(d::certify_root_model(root,pair,point.value(),options,256,nullptr,&statistics.certificate) &&
            statistics.certificate.radius_shrinks>0 && root.count==1,
            "deferred-jet control did not exercise full-jet shrink then proof");
    check_deferred_guard(pair,point,options);
    check_count(root,from_roots({{1,1},{Rational(1)-Rational(1)/128,1}}));
}
static void local_tail_witness() {
    // True polynomial w - 1/8 + 2^-16 w^17. Its first omitted term saturates
    // the local bound T*r^17 on |w|=r. The uniform T cannot prove 53 bits.
    d::RootModel m;
    arb_one(m.radius.value());
    d::ComplexBall coefficient;
    acb_set_d(coefficient.value(),-0.125);
    acb_poly_set_coeff_acb(m.polynomial.value(),0,coefficient.value());
    acb_one(coefficient.value());
    acb_poly_set_coeff_acb(m.polynomial.value(),1,coefficient.value());
    arb_one(m.taylor_model_error.value());
    arb_mul_2exp_si(m.taylor_model_error.value(),m.taylor_model_error.value(),-16);
    m.error=m.taylor_model_error;
    const auto point=ball({Rational(1)/8});
    cps_options options;options.goal=CPS_ROOT_APPROXIMATION;
    d::Root root;
    require(!d::certify_root_model(root,m,point.value(),options,256),
            "uniform-tail control unexpectedly proved 53 bits");
    m.taylor_error_order=17;
    require(d::certify_root_model(root,m,point.value(),options,256),
            "local Taylor remainder did not recover its exact witness");
    reference::ComplexPolynomial exact(18);
    exact[0]={Rational(-1)/8};exact[1]={1};exact[17]={power2(-16)};
    check_count(root,exact);
    require(!experiment::immutable_tail_rejects(m,point,options,256),
            "immutable floor ignored pointwise Taylor provenance");
    // Saturation catches a too-large exponent, which would be unsound.
    for (const Rational &r : {Rational(1)/8,Rational(1)/2,Rational(3)/4}) {
        const auto z=ball({r});
        const auto error=d::root_model_error(m,z.value(),256);
        Rational expected=power2(-16);
        for (unsigned j=0;j<17;++j) expected*=r;
        d::Rational value;fmpq_set_str(value.value(),expected.str().c_str(),10);
        require(arb_contains_fmpq(error.value(),value.value()),"Taylor exponent lost its saturating term");
    }
    const auto outside=ball({Rational(9)/8});
    const auto outside_error=d::root_model_error(m,outside.value(),256);
    require(arb_equal(outside_error.value(),m.error.value()),
            "pointwise Taylor bound escaped its analytic disk");
    m.taylor_error_order=0;
    const auto fallback=d::root_model_error(m,point.value(),256);
    require(arb_equal(fallback.value(),m.error.value()),
            "missing Taylor provenance did not retain uniform fallback");
}
int main() {
    certpolysolve::ThreadScope lifetime;
    try {
        local_tail_witness();
        deferred_guard_boundaries();
        std::vector<Complex> roots{{Rational(1)/8,Rational(1)/16},{Rational(-1)/4,Rational(1)/8},
                                   {Rational(1)/4,Rational(-1)/4},{2,-1}};
        const auto m=model(roots);
        std::vector<Complex> original;
        for(const auto &z:roots)original.push_back(z*Complex(Rational(1)/2)+Complex(1,1));
        const auto polynomial=from_roots(original);
        for(bool coarse:{true,false}) {
            cps_options options;options.goal=coarse?CPS_ROOT_ACCURACY:CPS_ROOT_APPROXIMATION;
            if(coarse)options.absolute_tolerance=d::Number("1/16");
            for(size_t i=0;i<3;++i) {
                auto point=ball(roots[i]+Complex(power2(coarse?-30:-90),power2(coarse?-31:-91)));
                check_deferred_guard(m,point,options);
                check_deferred_guard(m,ball(roots[i]+Complex(Rational(1)/32)),options);
                auto inaccurate=m;
                arb_set_d(inaccurate.error.value(),0.25);
                check_deferred_guard(inaccurate,point,options);
                d::Root rectangle,disc;
                require(experiment::derivative_rectangle(rectangle,m,point,options,256),"derivative rectangle failed separated exact root");
                require(experiment::derivative_disc(disc,m,point,options,256),"natural disc failed separated exact root");
                check_count(rectangle,polynomial);check_count(disc,polynomial);
                require(!experiment::immutable_tail_rejects(m,point,options,256),"zero tail rejected a valid proof");
            }
        }
        // The early rejection must allow a counted double root, even though
        // its linear derivative vanishes and a singleton disc cannot prove it.
        auto multiple=model({{0},{0}});
        arb_one(multiple.taylor_model_error.value());
        arb_mul_2exp_si(multiple.taylor_model_error.value(),multiple.taylor_model_error.value(),-180);
        multiple.error=multiple.taylor_model_error;
        cps_options options;options.goal=CPS_ROOT_APPROXIMATION;
        d::ComplexBall point;d::Root counted;
        require(d::certify_root_model(counted,multiple,point.value(),options,256) && counted.count==2,
                "control failed to exercise higher-order count");
        check_deferred_guard(multiple,point,options);
        require(!experiment::immutable_tail_rejects(multiple,point,options,256),"early floor excluded a valid higher-order count");
        arb_mul_2exp_si(multiple.taylor_model_error.value(),multiple.taylor_model_error.value(),148);
        multiple.error=multiple.taylor_model_error;
        require(experiment::immutable_tail_rejects(multiple,point,options,256),"immutable remainder floor did not reject impossible count");
        require(!d::certify_root_model(counted,multiple,point.value(),options,256),"control proved a floor-limited rectangle");
        // Nonzero analytic remainder: compare with the actual perturbed root,
        // not the zero of the truncated model.
        const Rational epsilon=power2(-100);
        auto perturbed=model({{Rational(1)/8}});
        d::Number(epsilon.str()).enclose(perturbed.error.value(),256);
        const auto estimate=ball({Rational(1)/8});
        d::Root disc;
        require(experiment::derivative_disc(disc,perturbed,estimate,options,256),"disc omitted nonzero model error");
        const auto exact=from_roots({Complex(1,1)+Complex(Rational(1)/8-epsilon)*Complex(Rational(1)/2)});
        check_count(disc,exact);
        auto invalid=estimate;acb_set_d(invalid.value(),1.1);
        require(!experiment::derivative_disc(disc,perturbed,invalid,options,256),"disc left the analytic model domain");
        std::cout<<"derivative rectangles/discs and higher-order-safe early rejection passed exact winding\n";
    }catch(const std::exception &e){std::cerr<<e.what()<<'\n';return 1;}
}
