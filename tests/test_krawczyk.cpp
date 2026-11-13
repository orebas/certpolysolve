#include "numeric/arithmetic/polynomial_arithmetic.hpp"
#include "exact_winding_api.hpp"
#include "support/giant-fft/krawczyk.hpp"
#include <iostream>
namespace d = certpolysolve::detail;
using reference::require;
static void
check(const std::vector<reference::Complex> &roots, const reference::Complex &saved, bool expected) {
    auto input = reference::input(reference::from_roots(roots));
    giant_fft::OriginalCoefficientVerifier verifier(certpolysolve::detail::dense(*input.get()));
    d::ComplexBall point;
    d::Number(saved.real.str()).enclose(acb_realref(point.value()), 512);
    d::Number(saved.imag.str()).enclose(acb_imagref(point.value()), 512);
    d::ComplexBall original = point;
    const auto result = verifier.verify(point.value(), 53, 512);
    require(result.accurate == expected, "independent verifier disagrees with planted roots");
    require(acb_equal(point.value(), original.value()), "saved point was modified");
    if (result.root) {
        size_t contained = 0;
        for (const auto &z : roots) {
            const auto ep = [&](size_t i) {
                const auto &v = result.root->endpoints[i];
                return reference::endpoint(certpolysolve::Dyadic(v.mantissa_text(), v.exponent_text()));
            };
            if (ep(0) <= z.real && z.real <= ep(1) && ep(2) <= z.imag && z.imag <= ep(3)) ++contained;
        }
        require(contained == 1, "Krawczyk singleton contains the wrong planted root count");
    }
}
int
main() {
    certpolysolve::ThreadScope lifetime;
    try {
        const auto tiny = reference::power2(-55);
        check({ { 1 }, { 2, 3 }, { -4, 1 } }, { 1 + tiny, tiny }, true);
        check({ { 1 }, { 2, 3 }, { -4, 1 } }, { 1 + reference::power2(-50), 0 }, false);
        check({ { 1 }, { 1 } }, { 1 }, false);
        check({ { 1 }, { 1 + reference::power2(-80) } }, { 1 }, false);
        check({ { 0 }, { 2, 3 } }, { 0 }, true);
        check({ { 0 }, { 0 }, { 2, 3 } }, { 0 }, false);
        check({ { 0, 1 }, { -3, 2 }, { -1, -4 } }, { tiny, 1 + tiny }, true);
        // A complex point that repeatedly rotates interval rectangles. The
        // exact root 1+i of z^1024-2^512 exposes artificial Horner wrapping.
        certpolysolve::Polynomial high;
        high.set(0, certpolysolve::Number("-1@512"));
        high.set(1024, certpolysolve::Number("1"));
        giant_fft::OriginalCoefficientVerifier verifier(certpolysolve::detail::dense(*high.get()));
        d::ComplexBall point;
        acb_set_d_d(point.value(), 1, 1);
        auto result = verifier.verify(point.value(), 53, 512);
        require(result.accurate && result.precision <= 256, "complex high-degree evaluation wrapped excessively");
        require(!verifier.verify(point.value(), 53, 512, [] { return true; }).accurate,
                "cancelled verification published a proof");
        std::cout << "immutable-point Krawczyk and exact corner checks passed\n";
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
