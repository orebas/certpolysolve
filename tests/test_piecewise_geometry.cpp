#include "numeric/arithmetic/polynomial_arithmetic.hpp"
#include "address_space_cap.hpp"
#include "exact_winding_api.hpp"
#include "methods/piecewise/piecewise_roots.hpp"
#include <iostream>

namespace d=certpolysolve::detail;
using reference::require;

int main() {
    certpolysolve::ThreadScope lifetime;
    try {
        d::RealBall inner,outer,gamma,radius,uncertainty;
        arb_one(inner.value());arb_set_si(outer.value(),2);
        arb_set_d(gamma.value(),1.5);arb_set_d(radius.value(),0.75);
        require(!d::piecewise_annulus_covered(inner,outer,gamma,radius,8,uncertainty,128),
                "eight charts leave uncovered corners of this annulus");
        require(d::piecewise_annulus_covered(inner,outer,gamma,radius,16,uncertainty,128),
                "sixteen charts should cover the entire annulus");
        arb_set_d(uncertainty.value(),0.5);
        require(!d::piecewise_annulus_covered(inner,outer,gamma,radius,16,uncertainty,128),
                "coverage ignored uncertain chart centers");
        // Strict coverage still succeeds with small radial and center balls.
        arb_set_d(uncertainty.value(),0.001);
        arb_add_error_2exp_si(inner.value(),-10);arb_add_error_2exp_si(outer.value(),-10);
        require(d::piecewise_annulus_covered(inner,outer,gamma,radius,16,uncertainty,128),
                "outward coverage rejected a well-separated valid cover");

        d::ComplexPolynomial polynomial;
        for(slong j=0;j<=257;++j)acb_poly_set_coeff_si(polynomial.value(),j,j%5-2);
        acb_poly_set_coeff_si(polynomial.value(),0,1);acb_poly_set_coeff_si(polynomial.value(),257,1);
        d::PiecewiseRootStatistics stats;
        auto config=d::maintained_piecewise_config();
        config.interior_proposals = config.adaptive_taylor = true;
        cps_options options;
        require(d::piecewise_root_models(polynomial.value(),16,192,options,stats,
            [](const auto &m){
                d::RealBall one;arb_one(one.value());
                require(arb_is_positive(m.proposal_radius.value()) && arb_lt(m.proposal_radius.value(),one.value()),
                        "covered proposal radius was not proved inside the chart");
                return acb_poly_degree(m.polynomial.value())>=0;
            },{},config),
            "proved cover failed to construct models");
        // Use independent planted roots, not the local model, to check counts.
        using namespace reference;
        const std::vector<Complex> roots{{1},{-2,1},{3,-2},{1,1},{-1,-1},{4,2},{Rational(1)/8,0}};
        const auto input=reference::input(from_roots(roots));
        options.goal=CPS_ROOT_APPROXIMATION;options.work_units=UINT64_MAX;
        const auto result=d::solve_piecewise_roots(certpolysolve::detail::dense(*input.get()),options,stats,config);
        require(result.goal_met && result.roots.size()==roots.size(),"covered proposal disks lost roots");
        require(stats.core_proposal_models>0,"covered proposal representation was not exercised");
        for(const auto &r:result.roots) {
            size_t count=0;
            for(const auto &z:roots) {
                d::Rational re,im;fmpq_set_str(re.value(),z.real.str().c_str(),10);fmpq_set_str(im.value(),z.imag.str().c_str(),10);
                count+=arb_contains_fmpq(acb_realref(r.box.value()),re.value()) && arb_contains_fmpq(acb_imagref(r.box.value()),im.value());
            }
            require(count==r.count,"angular proposal bypassed original-input count");
        }
        d::PiecewiseRootState state(certpolysolve::detail::dense(*input.get()), config);
        options.work_units = 7;
        (void)state.advance(options);
        options.work_units = UINT64_MAX;
        options.bits = 96;
        const auto admitted = state.advance(options);
        require(admitted.goal_met && admitted.roots.size()==roots.size(),
                "resumed interior proposals lost exact roots under a tighter query");
        for(const auto &r:admitted.roots) {
            size_t count=0;
            for(const auto &z:roots) {
                d::Rational re,im;fmpq_set_str(re.value(),z.real.str().c_str(),10);fmpq_set_str(im.value(),z.imag.str().c_str(),10);
                count+=arb_contains_fmpq(acb_realref(r.box.value()),re.value()) && arb_contains_fmpq(acb_imagref(r.box.value()),im.value());
            }
            require(count==r.count,"resumed interior proposals bypassed original-input count");
        }
        // Rings are made where terms meet, not across stretches where a single
        // term dominates (there the polynomial has no root). x^3 + 2^e x - 1,
        // e = 2^30 + 6: roots near 2^-e and 2^(e/2), and 1.6e9 octaves between
        // them. With a ring for every octave the list could not be allocated
        // (explicit piecewise ended in an allocation error); now there are
        // about a hundred rings.
        {
            const AddressSpaceCap cap(4ull << 30);
            d::ComplexPolynomial wide;
            d::ComplexBall big;
            acb_one(big.value());
            acb_mul_2exp_si(big.value(), big.value(), (slong(1) << 30) + 6);
            acb_poly_set_coeff_si(wide.value(), 0, -1);
            acb_poly_set_coeff_acb(wide.value(), 1, big.value());
            acb_poly_set_coeff_si(wide.value(), 3, 1);
            const auto estimate = d::estimate_piecewise_models(wide.value(), 16);
            require(estimate.available && estimate.rings > 0 && estimate.rings < 1000,
                    "the rings of a wide root window still follow its width");
            uint64_t models = 0;
            d::PiecewiseRootStatistics wide_stats;
            require(d::piecewise_root_models(wide.value(), 16, 192, cps_options{}, wide_stats,
                        [&](const auto &) {
                            ++models;
                            return true;
                        }, {}, d::maintained_piecewise_config()),
                    "the models of a wide root window were not constructed");
            require(models == estimate.models, "the estimate and the construction disagree on the rings");
            // The public route: no error, a budget stop or an answer.
            certpolysolve::Polynomial input;
            input.set(0, certpolysolve::Number("-1")).set(1, certpolysolve::Number("1@1073741830")).set(3, certpolysolve::Number("1"));
            certpolysolve::Options query;
            query.method = certpolysolve::Method::piecewise;
            query.goal = certpolysolve::Goal::approximation();
            query.work_budget = 2000;
            query.cell_budget = 65536;
            const auto answer = certpolysolve::solve(input, query);
            require(answer.stop() != certpolysolve::Stop::error && answer.get() != nullptr,
                    "explicit piecewise returned an error on a wide root window");
            const auto termination = cps_result_termination(answer.get());
            require(termination == CPS_SUCCESS || termination == CPS_WORK_LIMIT ||
                        termination == CPS_PRECISION_LIMIT,
                    "explicit piecewise did not end a wide root window on an answer or a budget");
        }
        // Without such a stretch nothing changes: the dense input above has
        // the rings it had (15 at 16 model bits, 1,728 models).
        {
            const auto estimate = d::estimate_piecewise_models(polynomial.value(), 16);
            require(estimate.rings == 15 && estimate.models == 1728, "the rings of a dense input changed");
        }
        std::cout<<"annulus coverage, center uncertainty and original roots passed\n";
    } catch(const std::exception &e){std::cerr<<e.what()<<'\n';return 1;}
}
