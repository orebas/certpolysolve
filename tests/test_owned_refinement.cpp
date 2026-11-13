#include "exact_winding_api.hpp"
#include "support/owned-refinement/experiments.hpp"
#include <iostream>
namespace d = certpolysolve::detail;
using namespace reference;

static d::RootModel model(const std::vector<Complex> &roots) {
    auto input = reference::input(from_roots(roots));
    d::RootModel result;
    arb_one(result.radius.value());
    d::ComplexBall c;
    for (const auto &[i, coefficient] : input.get()->terms) {
        coefficient.enclose(c.value(), 256);
        acb_poly_set_coeff_acb(result.polynomial.value(), static_cast<slong>(i), c.value());
    }
    return result;
}
static void verify(const d::Root &proof, const std::vector<Complex> &roots) {
    auto bound = [&](size_t k) { return reference::endpoint(certpolysolve::Dyadic{
        proof.endpoints[k].mantissa_text(), proof.endpoints[k].exponent_text()}); };
    size_t count = 0;
    for (const auto &r : roots)
        count += bound(0) < r.real && r.real < bound(1) && bound(2) < r.imag && r.imag < bound(3);
    require(count == 1 && proof.count == 1, "residual-radius proof has wrong exact root count");
}
int main() {
    certpolysolve::ThreadScope lifetime;
    try {
        const std::vector<Complex> roots{{Rational(1)/4,Rational(1)/8},
                                       {-Rational(1)/3,Rational(1)/5}, {Rational(3)/4}};
        auto m = model(roots);
        d::ComplexBall point;
        acb_set_d_d(point.value(), .25 + 0x1p-30, .125 - 0x1p-31);
        cps_options options;
        d::Root proof;
        require(owned_refinement::residual_guard(proof,m,point,options,256), "residual certificate refused simple root");
        verify(proof,roots);
        // Exact affine model coordinates; verify against independently mapped roots.
        acb_set_si_si(m.center.value(),1,-2);
        arb_set_si(m.radius.value(),3);
        std::vector<Complex> moved;
        for (const auto &r:roots) moved.push_back(r*Complex(3)+Complex(1,-2));
        require(owned_refinement::residual_guard(proof,m,point,options,256), "affine certificate failed");
        verify(proof,moved);
        auto repeated = model({roots[0],roots[0]});
        acb_set_d_d(point.value(),.25,.125);
        require(!owned_refinement::residual_guard(proof,repeated,point,options,256), "multiple root passed singleton proof");
        auto perturbed = model({roots[0]});
        arb_set_d(perturbed.error.value(),0x1p-20);
        require(owned_refinement::residual_guard(proof,perturbed,point,options,256), "bounded model error failed");
        for (int re : {-1,1}) for (int im : {-1,1}) {
            // |constant perturbation| <= 2^-20. Check four distinct admissible sources.
            verify(proof,{roots[0]+Complex(Rational(re)/4194304,Rational(im)/4194304)});
        }
        arb_one(perturbed.error.value());
        require(!owned_refinement::residual_guard(proof,perturbed,point,options,256), "uncertainty escaped chart proof");
        acb_indeterminate(point.value());
        require(!owned_refinement::residual_guard(proof,m,point,options,256), "nonfinite point admitted");
        std::cout << "Residual-radius linear Rouche checks passed\n";
    } catch (const std::exception &e) { std::cerr << e.what() << '\n'; return 1; }
}
