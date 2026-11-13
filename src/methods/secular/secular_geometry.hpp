#pragma once
// The geometry of a secular restart about a centre: functions of the model's nodes and residues alone. The
// widest gap between the nodes' distances from the centre and the ladder of circles through it; an isolating
// plateau of the model's exact counts (secular_disk_count.hpp does the counting); the mean
// of the roots inside it; the span a restart's hints need; and the roots about a tight group of nodes (the
// geometric restart's and the split's clusters).
// History: docs/code-history-20261021.md, secular_geometry.hpp 1.
#include "core/raii.hpp"
#include "certify/secular_disk_count.hpp"
#include <cstdint>
#include <optional>
#include <span>
#include <vector>

namespace certpolysolve::detail {
// A circle about a centre between the radii 2^low_log2 and 2^high_log2: the
// point of [low_log2, high_log2] farthest, in log2, from every node's
// distance to the centre (distances_log2, sorted), so that no node lies near
// the circle (a node near it lengthens an exact count's exterior expansion
// and its arcs, and costs the moments' arithmetic its accuracy).
double widest_gap_log2(std::span<const double> distances_log2, double low, double high);
// A ladder's circle near the radius 2^target_log2, within half a rung.
double ladder_radius_log2(std::span<const double> distances_log2, double target_log2);
// The radius 2^radius_log2 as an arf.
void set_radius_log2(arf_t radius, double radius_log2);
// An isolating plateau of the model's exact counts about a centre: walking
// down circles from about 2^top_log2 by halves (ladder_radius_log2), the
// outermost two consecutive circles with the same count K >= 2, extended
// inward while the count stays. Their annulus r_in <= |z - c| < r_out holds
// no root: the disk of radius r_in holds exactly K roots and every other root
// lies at r_out or beyond. The walk ends at a count of at most one, at an
// undecided circle after a decided one, or at two undecided circles in a row
// (inside a cluster whose nodes have not arrived, the model's terms cancel
// beyond the count's one limb, and deeper circles are undecided too), and
// after a limb's worth of halvings.
struct Plateau {
    slong count = 0;
    double inner_log2 = 0, outer_log2 = 0;
};
// log2 |b_j - c| for every node (-inf at the centre), in node order.
std::vector<double> node_distances_log2(std::span<const ComplexBall> nodes, const acb_t centre);
std::optional<Plateau> isolating_plateau(std::span<const ComplexBall> nodes, std::span<const ComplexBall> residues,
                                         const acb_t centre, std::span<const double> distance,
                                         std::span<const double> sorted, double top_log2, DiskCountScratch &scratch);
// The mean of the K roots inside a plateau's inner circle about c: the first
// moment sum (z_i - c) over the roots in the disk, (1/2 pi i) \oint (z - c)
// p'/p dz, divided by K. p'/p = sum_j 1/(z - b_j) + S'/S from the model has
// no pole in the plateau's annulus, so the trapezoid rule on N points of a
// circle of radius rho there converges as (r_in / rho)^N + (rho / r_out)^N
// at least; the circle is the widest gap between the nodes' distances (the
// poles of the two parts, whose sum the arithmetic cancels) in the annulus'
// middle half.
// The points double from 16 until successive estimates agree to binary64's
// resolution of the circle (the nested solve resolves no finer) or to the
// arithmetic's own rounding, or until the certified rate guarantees that
// resolution; the zeroth moment, which is K, checks the rule. A heuristic (no
// enclosure): the count is proved, the mean only centres the local model.
// Returns false when the check fails or the mean leaves the inner circle.
bool root_mean(acb_t mean, std::span<const ComplexBall> nodes, std::span<const ComplexBall> residues,
               const acb_t centre, const Plateau &plateau, std::span<const double> distances_log2,
               double &operations, uint64_t maximum_points = UINT64_MAX);
// The bits between a restart's hints' largest magnitude and their closest
// pair: the coordinates must carry that many to keep the hints apart. The
// hints are exact; their offsets from the centre are taken at `precision`
// (the model's, which formed them) and scaled into binary64's range before
// the pairs are compared, all of them: quadratic in the group's size, a
// fraction of the model that produced the hints. Zero for fewer than two
// hints.
double hint_span_log2(std::span<const ComplexBall> hints, const acb_t centre, slong precision);
// The roots about a tight group of nodes. First as split_by_exact_counts counts them: secular_disk_count on the
// circles of radius r and 2 r about the group's mean, r twice the members' largest distance from it, with every
// other node outside the second circle; equal counts put that many roots within r of the mean and none between
// the circles. Where those two circles do not agree (a root of the group's cluster that none of its nodes
// serves lies between them, or on one), the isolating plateau of the counts from just inside the nearest other
// node inward, as the geometric restart finds its clusters: its count is the cluster's, within the plateau's
// inner circle. `radius` is then that circle's, or r where r is the larger. `member` is all zero on entry and on
// return.
struct CountedGroup {
    slong count = 0;
    ComplexBall centre;
    RealBall radius; // exact
};
std::optional<CountedGroup> counted_group(const std::vector<size_t> &group, std::span<const ComplexBall> nodes,
                                          std::span<const ComplexBall> residues, std::vector<uint8_t> &member,
                                          slong precision, DiskCountScratch &scratch);
} // namespace certpolysolve::detail
