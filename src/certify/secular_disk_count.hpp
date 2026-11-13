#pragma once
#include "core/raii.hpp"
#include <cstdint>
#include <optional>
#include <span>
#include <vector>

namespace certpolysolve::detail {
// secular_disk_count: the exact number of roots of p in the open disk
// |z - c| < rho, from the model alone. For monic p with exact distinct nodes
// b_j and residues w_j = -p(b_j) / Q'(b_j), Q = prod (z - b_j),
//     p(z) = Q(z) S(z),   S(z) = 1 - sum_j w_j / (z - b_j).
// Let the m member nodes lie strictly inside the circle and every other node
// strictly outside. S is meromorphic in the disk with poles only at the
// member nodes, so by the argument principle its winding number along the
// circle is (roots of p in the disk) - m. (A root that coincides with a
// member node cancels that node's pole and changes neither side.)
//
// The winding number is computed with enclosures, without evaluating p:
// - The part of S from the outside nodes beyond twice the radius, E(z) = 1 -
//   sum_{j beyond} w_j / (z - b_j), is analytic in the disk. It is expanded
//   about c, each node's series cut where its geometric remainder on the
//   circle is negligible (one bit a term at least; far nodes need few terms),
//   and the remainders are bounded together.
// - The circle is cut into arcs. On each arc S is enclosed in a rectangle
//   (the expansion, the terms of the m members and of the outside nodes
//   within twice the radius exactly, the remainder), and the arc
//   is bisected until its rectangle lies in one open half-plane Re > 0,
//   Im > 0, Re < 0 or Im < 0. Neighbouring arcs share a value of S, so their
//   half-planes differ by a quarter turn at most, and the quarter turns
//   around the circle sum to four times the winding number.
// The count is exact when it is returned. nullopt means undecided (a root
// too near the circle for the arc budget, or enclosures too wide), and
// callers must then do nothing.
// History: docs/code-history-20261021.md, secular_refinement.hpp 2.
struct DiskCountScratch {
    std::vector<ComplexBall> inverse, moments;
    std::vector<RealBall> weight, ratio;
    // One-limb operations the counts made (series terms, then a term an arc
    // for the expansion and each direct node), accumulated for the meter.
    double operations = 0;
};
std::optional<slong> secular_disk_count(std::span<const ComplexBall> nodes, std::span<const ComplexBall> residues,
                                        std::span<const uint8_t> member, const acb_t centre, const arf_t radius,
                                        DiskCountScratch &scratch);

} // namespace certpolysolve::detail
