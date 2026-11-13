#pragma once
// Argument-principle count from original-source enclosures. No proposed points
// or candidate iteration state enter this proof.
#include "core/engine_statistics.hpp"
#include "core/query.hpp"
#include <cstdint>
#include <optional>

namespace certpolysolve::detail {
class MonicSource;
// The number of roots, counted with multiplicity, of the source's polynomial in the open disk |z - centre| < radius,
// by the argument principle from evaluations over balls. The circle is cut into arcs with dyadic ends. On an arc,
// the image lies in the mean-value enclosure p(z_m) + p'(B) (B - z_m) and in the direct enclosure p(B), where B is
// a ball about the arc's midpoint z_m that holds the arc; one that excludes 0 is a convex set missing the origin, so
// the argument changes along the arc by less than a half-turn, by the principal value Arg(p(z_b) conj(p(z_a))). An
// arc whose enclosures both meet 0 is bisected. The changes, summed in balls and divided by 2 pi, hold a unique
// integer: the count. No count (`roots` empty) when an evaluation fails, when a value at an arc's end has no
// argument to a few bits at this precision (more arcs cannot help), or when the arcs run out (the circle passes too
// close to a root). Cancellation leaves `roots` empty and says so.
struct WindingCount {
    std::optional<uint64_t> roots;
    uint64_t arcs = 0, evaluations = 0;
    bool cancelled = false;
    // No count because a value at an arc's end was not readable at this precision (more precision may help, more
    // arcs cannot); otherwise a refused count ran out of arcs or met a failed evaluation.
    bool imprecise = false;
};
// `centre` is taken at its midpoint and `radius` exactly; at most `max_arcs` arcs are evaluated.
WindingCount certified_root_count(MonicSource &source, const acb_t centre, const mag_t radius, slong precision,
                                  uint64_t max_arcs, const cps_options &options, certpolysolve::detail::CandidateStatistics &statistics);

} // namespace certpolysolve::detail
