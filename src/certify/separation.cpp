// The separation of factor covers that meet (separation.hpp): the separated goal, the rule that raises a factor's
// separation bits, and the choice of the factor to tighten.
//
// State: none. Entry points: separation.hpp's three functions. Not here: the callers' loops (solver.cpp's advance,
// prepared_piecewise.cpp's factor loop). Design notes: docs/mathematics.md (the global driver's factors; piecewise's
// factor loop).
#include "certify/separation.hpp"
#include "core/raii.hpp"

namespace certpolysolve::detail {
cps_options separated_options(const cps_options &base, uint64_t separation) {
    auto o = base;
    if (separation) {
        o.bits = std::min<uint64_t>(resources::precision_max, o.bits + separation);
        o.absolute_tolerance = Number();
        if (o.goal == CPS_CLUSTER_RADIUS)
            fmpq_div_2exp(o.cluster_radius.significand.value(), o.cluster_radius.significand.value(),
                          static_cast<ulong>(separation));
    }
    return o;
}

uint64_t next_separation(uint64_t bits, uint64_t separation) {
    return std::min<uint64_t>(resources::precision_max - bits,
                              std::max(resources::add(separation, 16),
                                       resources::add(bits, resources::multiply(separation, 2))));
}

std::optional<size_t> factor_to_separate(std::span<const std::vector<Root> *const> covers) {
    for (size_t i = 0; i < covers.size(); ++i)
        for (size_t j = 0; j < i; ++j)
            for (const auto &a : *covers[i])
                for (const auto &b : *covers[j])
                    if (acb_overlaps(a.box.value(), b.box.value())) {
                        RealBall ar, br;
                        mag_add(arb_radref(ar.value()), arb_radref(acb_realref(a.box.value())),
                                arb_radref(acb_imagref(a.box.value())));
                        mag_add(arb_radref(br.value()), arb_radref(acb_realref(b.box.value())),
                                arb_radref(acb_imagref(b.box.value())));
                        return mag_cmp(arb_radref(ar.value()), arb_radref(br.value())) >= 0 ? i : j;
                    }
    return std::nullopt;
}
} // namespace certpolysolve::detail
