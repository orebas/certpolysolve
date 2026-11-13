// box_contains (certification.cpp) decides exactly as acb_contains: random boxes about the boundary, inside and
// outside, with midpoints of up to 128 bits (as the compensated certificate's published rectangles have) at scales
// across binary64's range and beyond it; then the edges: a box in itself, an exact point on the boundary, a box a
// rounding outside, zero and infinite radii, NaN, and midpoints 2^10000 apart.
#include "core/cover.hpp"
#include "core/raii.hpp"
#include <flint/acb.h>
#include <cmath>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>

namespace d = certpolysolve::detail;

namespace {
void require(bool ok, const std::string &why) {
    if (!ok) throw std::runtime_error(why);
}
// A part's midpoint: c 2^scale plus a word 2^-60 below it and a third word below that (a midpoint of up to about
// 128 bits).
void midpoint(arf_t out, double c, double low, double lower, slong scale) {
    d::Float word;
    arf_set_d(out, c);
    arf_set_d(word.value(), low);
    arf_mul_2exp_si(word.value(), word.value(), -60);
    arf_add(out, out, word.value(), ARF_PREC_EXACT, ARF_RND_DOWN);
    arf_set_d(word.value(), lower);
    arf_mul_2exp_si(word.value(), word.value(), -120);
    arf_add(out, out, word.value(), ARF_PREC_EXACT, ARF_RND_DOWN);
    arf_mul_2exp_si(out, out, scale);
}
uint64_t compare(const acb_t outer, const acb_t inner, const std::string &what) {
    require(d::box_contains(outer, inner) == (acb_contains(outer, inner) != 0), what + ": differs from acb_contains");
    return 1;
}
} // namespace

int main() {
    try {
        std::mt19937_64 rng(20261026);
        std::uniform_real_distribution<double> unit(-1, 1), share(0, 1.25), shift(-1.5, 1.5);
        uint64_t checked = 0, contained = 0;
        for (int trial = 0; trial < 100000; ++trial) {
            const slong scale = static_cast<slong>(rng() % 2400) - 1200;
            d::ComplexBall outer, inner;
            for (int part = 0; part < 2; ++part) {
                arb_ptr o = part ? acb_imagref(outer.value()) : acb_realref(outer.value());
                arb_ptr i = part ? acb_imagref(inner.value()) : acb_realref(inner.value());
                const double c = unit(rng);
                midpoint(arb_midref(o), c, unit(rng), unit(rng), scale);
                // The outer radius about 2^-58 of the midpoint, the inner box shifted by up to 1.5 of it, its radius
                // up to 1.25 of it: about half of the boxes lie inside, many within a rounding of the boundary.
                const double radius = std::ldexp(1.0 + share(rng), -58);
                d::Float r;
                arf_set_d(r.value(), radius);
                arf_mul_2exp_si(r.value(), r.value(), scale);
                arf_get_mag(arb_radref(o), r.value());
                d::Float offset;
                arf_set_d(offset.value(), shift(rng) * radius * (trial % 3 ? 0.5 : 1));
                arf_mul_2exp_si(offset.value(), offset.value(), scale);
                arf_add(arb_midref(i), arb_midref(o), offset.value(), ARF_PREC_EXACT, ARF_RND_DOWN);
                arf_set_d(r.value(), radius * share(rng) * (trial % 3 ? 0.5 : 1));
                arf_mul_2exp_si(r.value(), r.value(), scale);
                arf_get_mag(arb_radref(i), r.value());
            }
            checked += compare(outer.value(), inner.value(), "random boxes");
            contained += acb_contains(outer.value(), inner.value()) != 0;
        }
        require(contained > checked / 10 && contained < checked - checked / 10,
                "the random boxes lie both inside and outside");
        // The edges.
        d::ComplexBall a, b;
        acb_set_d_d(a.value(), 1.5, -0.25);
        mag_set_d(arb_radref(acb_realref(a.value())), 0x1p-58);
        mag_set_d(arb_radref(acb_imagref(a.value())), 0x1p-58);
        checked += compare(a.value(), a.value(), "a box in itself");
        // An exact point on the boundary, then one a rounding beyond it.
        acb_set(b.value(), a.value());
        mag_zero(arb_radref(acb_realref(b.value())));
        mag_zero(arb_radref(acb_imagref(b.value())));
        d::Float end;
        arf_set_mag(end.value(), arb_radref(acb_realref(a.value())));
        arf_add(arb_midref(acb_realref(b.value())), arb_midref(acb_realref(a.value())), end.value(), ARF_PREC_EXACT,
                ARF_RND_DOWN);
        checked += compare(a.value(), b.value(), "an exact point on the boundary");
        require(d::box_contains(a.value(), b.value()), "the boundary point lies in the box");
        arf_mul_2exp_si(end.value(), end.value(), -200);
        arf_add(arb_midref(acb_realref(b.value())), arb_midref(acb_realref(b.value())), end.value(), ARF_PREC_EXACT,
                ARF_RND_DOWN);
        checked += compare(a.value(), b.value(), "a point 2^-258 beyond the boundary");
        require(!d::box_contains(a.value(), b.value()), "the point beyond the boundary lies outside");
        // Zero radii: an exact point in an exact box, and another point.
        acb_set_d_d(a.value(), 2.0, 3.0);
        acb_set(b.value(), a.value());
        checked += compare(a.value(), b.value(), "an exact point in itself");
        acb_set_d_d(b.value(), 2.0, 3.5);
        checked += compare(a.value(), b.value(), "an exact point in another");
        // Infinite radius and NaN.
        mag_inf(arb_radref(acb_realref(a.value())));
        checked += compare(a.value(), b.value(), "an infinite radius outside");
        checked += compare(b.value(), a.value(), "an infinite radius inside");
        acb_indeterminate(b.value());
        checked += compare(a.value(), b.value(), "NaN inside");
        checked += compare(b.value(), a.value(), "NaN outside");
        // Midpoints 2^10000 apart, and a huge box about a tiny one.
        acb_set_d_d(a.value(), 1.0, 1.0);
        arb_mul_2exp_si(acb_realref(a.value()), acb_realref(a.value()), 5000);
        mag_set_d(arb_radref(acb_realref(a.value())), 1.0);
        mag_set_d(arb_radref(acb_imagref(a.value())), 1.0);
        acb_set_d_d(b.value(), 1.0, 1.0);
        arb_mul_2exp_si(acb_realref(b.value()), acb_realref(b.value()), -5000);
        checked += compare(a.value(), b.value(), "midpoints 2^10000 apart");
        mag_set_d(arb_radref(acb_realref(a.value())), 1.0);
        mag_mul_2exp_si(arb_radref(acb_realref(a.value())), arb_radref(acb_realref(a.value())), 5001);
        checked += compare(a.value(), b.value(), "a huge box about a tiny one");
        std::cout << "box contains: " << checked << " boxes decided as acb_contains decides them (" << contained
                  << " random boxes inside)\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << "test_box_contains: " << e.what() << '\n';
        return 1;
    }
}
