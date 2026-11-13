// Lift a proved cover through the exact map x -> x^k. Branches avoid zero and
// the branch cut; their images exhaust the preimages with counts preserved.
// Callers merge any overlapping images before publishing the complete cover.
#include "certify/power_cover.hpp"

namespace certpolysolve::detail {
bool lift_power_cover(std::vector<Root> &out, std::span<const Root> reduced, ulong power, slong precision) {
    if (power < 2) {
        out.assign(reduced.begin(), reduced.end());
        return true;
    }
    if (reduced.size() > out.max_size() / power)
        return false;
    RootVector unity(static_cast<slong>(power));
    _acb_vec_unit_roots(unity.data(), static_cast<slong>(power), static_cast<slong>(power), precision);
    std::vector<Root> proposed;
    proposed.reserve(reduced.size() * power);
    for (const auto &root : reduced) {
        if (!acb_is_finite(root.box.value()) || acb_contains_zero(root.box.value()))
            return false;
        ComplexBall principal, rotation;
        acb_set(principal.value(), root.box.value());
        const bool rotate = arb_is_negative(acb_realref(principal.value()));
        if (rotate)
            acb_neg(principal.value(), principal.value());
        // If a rectangle avoids zero but crosses the negative real axis, its
        // real part is strictly negative, so the rotation above removes it.
        acb_root_ui(principal.value(), principal.value(), power, precision);
        if (rotate) {
            acb_unit_root(rotation.value(), 2 * power, precision);
            acb_mul(principal.value(), principal.value(), rotation.value(), precision);
        }
        for (ulong j = 0; j < power; ++j) {
            Root lifted;
            acb_mul(lifted.box.value(), principal.value(), unity[static_cast<slong>(j)], precision);
            if (!acb_is_finite(lifted.box.value()))
                return false;
            // Enclose the analytic image in a balanced square with an interior
            // margin. Axis-aligned branches can otherwise produce zero-width
            // rectangles, which are awkward for independent contour checks
            // and subsequent two-dimensional component refinement.
            Magnitude radius;
            mag_max(radius.value(), arb_radref(acb_realref(lifted.box.value())),
                    arb_radref(acb_imagref(lifted.box.value())));
            mag_mul_2exp_si(radius.value(), radius.value(), 2);
            mag_set(arb_radref(acb_realref(lifted.box.value())), radius.value());
            mag_set(arb_radref(acb_imagref(lifted.box.value())), radius.value());
            lifted.count = root.count;
            lifted.multiplicity = root.multiplicity;
            lifted.realness =
                !arb_contains_zero(acb_imagref(lifted.box.value())) ? CPS_NONREAL : CPS_REAL_UNKNOWN;
            lifted.certificate = "exact-power-composition;" + root.certificate;
            lifted.publish(precision);
            proposed.push_back(std::move(lifted));
        }
    }
    out = std::move(proposed);
    return true;
}
} // namespace certpolysolve::detail
