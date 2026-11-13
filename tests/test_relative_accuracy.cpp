// The goal tests' relative accuracy in its two forms (certification.cpp): the binary64 form, tried first, may hold only
// where the ball form holds, so every goal test, and every answer, is the same. Random boxes about the threshold, at
// the goals and working precisions the solver uses, real and complex, with one part near zero or both large; then the
// edges: exact points, boxes holding zero, parts beyond the binary64 form's range. Boxes well inside the bound must
// be decided in binary64 (the saving), and none it decides may be one the ball form refuses.
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
void set_part(arb_t part, double mid, double radius) {
    arb_set_d(part, mid);
    d::Magnitude r;
    mag_set_d(r.value(), radius);
    mag_set(arb_radref(part), r.value());
}
} // namespace

int main() {
    try {
        std::mt19937_64 rng(20261026);
        std::uniform_real_distribution<double> unit(-1, 1), offset(-3, 3), scale(-60, 60);
        uint64_t checked = 0, decided = 0, wide_margin = 0, wide_decided = 0;
        for (const uint64_t bits : {uint64_t(10), uint64_t(53), uint64_t(80), uint64_t(105)})
            for (const slong precision : {slong(53), slong(85), slong(128), slong(320)})
                for (int trial = 0; trial < 4000; ++trial) {
                    const double size = std::exp2(scale(rng));
                    const double re = unit(rng) * size, im = (trial % 3 ? unit(rng) : unit(rng) * 0x1p-40) * size;
                    const double gap = offset(rng); // log2 of the box's error against the bound
                    const double bound = std::ldexp(std::max(std::abs(re), std::abs(im)), -static_cast<int>(bits));
                    const double total = bound * std::exp2(gap);
                    const double share = std::uniform_real_distribution<double>(0, 1)(rng);
                    d::ComplexBall box;
                    set_part(acb_realref(box.value()), re, total * share);
                    set_part(acb_imagref(box.value()), im, total * (1 - share));
                    const bool fast = d::relative_accuracy_binary64(box.value(), bits, precision);
                    const bool ball = d::relative_accuracy_ball(box.value(), bits, precision);
                    require(!fast || ball, "the binary64 form holds where the ball form does not (bits " +
                                               std::to_string(bits) + ", precision " + std::to_string(precision) + ")");
                    ++checked;
                    decided += fast;
                    if (gap < -1) {
                        ++wide_margin;
                        wide_decided += fast;
                    }
                }
        require(wide_decided == wide_margin, "a box well inside the bound was left to the ball form");
        // Edges: an exact point (both forms hold), a box holding zero (neither may hold), parts beyond 2^1000 (left
        // to the ball form), a goal past the form's range (left to the ball form).
        {
            d::ComplexBall point;
            acb_set_d_d(point.value(), 1.5, -0.25);
            require(d::relative_accuracy_binary64(point.value(), 53, 85) && d::relative_accuracy_ball(point.value(), 53, 85),
                    "an exact point");
            d::ComplexBall zero;
            set_part(acb_realref(zero.value()), 1e-20, 1e-10);
            set_part(acb_imagref(zero.value()), -1e-20, 1e-10);
            require(!d::relative_accuracy_binary64(zero.value(), 53, 85), "a box holding zero");
            d::ComplexBall huge;
            arb_set_d(acb_realref(huge.value()), 1.25);
            arb_mul_2exp_si(acb_realref(huge.value()), acb_realref(huge.value()), 4000);
            arb_zero(acb_imagref(huge.value()));
            require(!d::relative_accuracy_binary64(huge.value(), 53, 85) && d::relative_accuracy_ball(huge.value(), 53, 85),
                    "a part beyond 2^1000 is the ball form's");
            require(!d::relative_accuracy_binary64(point.value(), 5000, 6000), "a goal past the binary64 form's range");
        }
        std::cout << "relative accuracy: " << checked << " boxes, " << decided << " decided in binary64 (" << wide_decided
                  << " of " << wide_margin << " more than twice inside the bound), none against the ball form\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << "test_relative_accuracy: " << e.what() << '\n';
        return 1;
    }
}
