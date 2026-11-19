#pragma once
// Which points of a set mirror each other in the real axis: the guess the certificate attempt makes on a polynomial
// with real coefficients before its first pass (newton_cover.cpp, "Conjugate pairs"; docs/mathematics.md, "The
// certificate attempt"). Nothing is proved by it: a pair it misses is certified point by point, and a wrong pair
// leaves two rectangles on one root, which the cover's disjointness test refuses.
//
// Two points are partners where
// - one lies strictly above the real axis and one strictly below;
// - they mirror each other: their real parts, and the sizes of their imaginary parts, each differ by at most
//   2^-bits of the larger point's size (the larger of the four coordinate sizes), the same test from either side;
// - each is the only point of the other side that mirrors the other: where a point has two candidates, as in a
//   cluster tighter than the tolerance, neither it nor they are paired.
// The answer is a function of the set of points alone, whatever their order.
//
// The candidates of a point are looked for among the points of the other side whose real part lies within twice
// the tolerance of its own, or among those whose imaginary size does, whichever are fewer (two sorted orders, a
// binary search each). A point with more than `scan_limit` of both is left without a partner: that bounds the work
// at scan_limit comparisons a point, and such a point sits in a lattice or a cluster of that many points a side.
#include <cstddef>
#include <span>
#include <vector>

namespace certpolysolve::detail {
struct PlanePoint {
    double re = 0, im = 0;
};
// By point: the index of its partner, or -1. A point with a part that is not finite, or with imaginary part zero,
// has none.
std::vector<std::ptrdiff_t> conjugate_partners(std::span<const PlanePoint> points, int bits, std::size_t scan_limit);
} // namespace certpolysolve::detail
