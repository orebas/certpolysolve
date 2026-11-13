// The interaction tree (InteractionTreeOf, declared in interaction_tree.hpp): a Barnes-Hut-style treecode, in binary64,
// for the sums S_t = sum_{j != t} 1/(z_t - z_j) in Aberth's correction N/(1 - N S_t). Distant points enter as truncated
// multipole series, one tree walk per target and no local expansions (not the fast multipole method); candidate-only.
//
// build splits the points at the median along the longer side of their bounding box, down to leaves of at most 12
// points (or of coincident ones); each node keeps its centre c, radius r, its points' largest rounding error and the
// moments m_k = sum ((z_j - c)/r)^k, k < Order (build_two_part: the same for points given as two doubles). evaluate
// (visit) takes a node without the target t as a far field when r <= FarField::acceptance |t - c| and the rounding
// errors are at most interaction_truncation x (1 - acceptance) |t - c|, adding sum_k m_k (r/(t - c))^k / (t - c) to
// `far` and its truncation bound to `tail`; the other leaves it reaches list their points in `near`, for the caller to
// sum directly. move (binary64 points only) updates in place the moments, and if need be the radii, on a moved point's
// path from the root.
//
// State: per object only: the points' copies (points_, lows_, errors_), the permutation (indices_, positions_), nodes_.
// Entry points (the three instantiations at the end of this file): InteractionTree (Order 20) in candidates.cpp
// (prepare_interactions, for secular_candidate_packet and candidate_packet); MachineInteractionTree in
// machine_candidates.cpp (MachineEngine::advance, also move) and program_certificate.cpp (SeparatedRootsPass::advance);
// and the program route's tree (heuristics::program_tree_order, heuristics::ProgramFarField) in program_machine.cpp
// (ProgramMachineStage::sweep).
// Design notes: docs/fast-methods-assessment-20261008.md (the tree among fast methods; why its sums are proposals),
// docs/slp-structural-seeding-note-20261018.md (the program stage's coarser rule, and why the rule is a type).
#include "numeric/candidates/interaction_tree.hpp"
#include "core/heuristics.hpp"
#include <algorithm>
#include <cmath>
#include <numeric>

namespace certpolysolve::detail {
template <size_t Order, class FarField>
InteractionTreeOf<Order, FarField>::InteractionTreeOf(std::span<const std::complex<double>> points,
                                                       std::span<const double> errors,
                                                       std::span<const std::complex<double>> lows)
    : points_(points.begin(), points.end()), lows_(lows.begin(), lows.end()), errors_(errors.begin(), errors.end()),
      indices_(points.size()), positions_(points.size()) {
    std::iota(indices_.begin(), indices_.end(), size_t(0));
    nodes_.reserve(2 * points.size());
    if (!points.empty())
        build(0, points.size());
    for (size_t j = 0; j < indices_.size(); ++j)
        positions_[indices_[j]] = j;
}
template <size_t Order, class FarField>
size_t InteractionTreeOf<Order, FarField>::build(size_t begin, size_t end) {
    const size_t index = nodes_.size();
    nodes_.emplace_back();
    auto& node = nodes_[index]; // reserve above prevents invalidation
    node.begin = begin;
    node.end = end;
    if (!lows_.empty())
        return build_two_part(index, begin, end);
    double xmin = points_[indices_[begin]].real(), xmax = xmin;
    double ymin = points_[indices_[begin]].imag(), ymax = ymin;
    for (size_t j = begin; j < end; ++j) {
        const auto z = points_[indices_[j]];
        xmin = std::min(xmin, z.real());
        xmax = std::max(xmax, z.real());
        ymin = std::min(ymin, z.imag());
        ymax = std::max(ymax, z.imag());
    }
    node.center = {xmin / 2 + xmax / 2, ymin / 2 + ymax / 2};
    for (size_t j = begin; j < end; ++j)
        node.radius = std::max(node.radius, std::abs(points_[indices_[j]] - node.center));
    if (!errors_.empty())
        for (size_t j = begin; j < end; ++j)
            node.error = std::max(node.error, errors_[indices_[j]]);
    node.moments[0] = static_cast<double>(end - begin);
    build_terms_ += (end - begin) * order;
    if (node.radius > 0)
        for (size_t j = begin; j < end; ++j) {
            const auto displacement = (points_[indices_[j]] - node.center) / node.radius;
            std::complex<double> power = 1;
            for (size_t k = 1; k < order; ++k) {
                power *= displacement;
                node.moments[k] += power;
            }
        }
    if (end - begin > 12 && node.radius > 0) {
        const size_t middle = begin + (end - begin) / 2;
        const bool real_axis = xmax - xmin >= ymax - ymin;
        std::nth_element(indices_.begin() + static_cast<ptrdiff_t>(begin),
                         indices_.begin() + static_cast<ptrdiff_t>(middle),
                         indices_.begin() + static_cast<ptrdiff_t>(end), [&](size_t a, size_t b) {
                             return real_axis ? points_[a].real() < points_[b].real()
                                              : points_[a].imag() < points_[b].imag();
                         });
        node.left = build(begin, middle);
        node.right = build(middle, end);
    }
    return index;
}
// The same with points of two doubles: the box, the radius, the moments and
// the split come from displacements from the range's first point, which keep
// the points' second parts where their first parts coincide. The centre is
// that point plus the box's middle.
template <size_t Order, class FarField>
size_t InteractionTreeOf<Order, FarField>::build_two_part(size_t index, size_t begin, size_t end) {
    auto& node = nodes_[index];
    const size_t first = indices_[begin];
    const auto reference = points_[first], reference_low = lows_[first];
    double xmin = 0, xmax = 0, ymin = 0, ymax = 0;
    for (size_t j = begin; j < end; ++j) {
        const auto d = displacement(indices_[j], reference, reference_low);
        xmin = std::min(xmin, d.real());
        xmax = std::max(xmax, d.real());
        ymin = std::min(ymin, d.imag());
        ymax = std::max(ymax, d.imag());
    }
    node.center = reference;
    node.center_low = reference_low + std::complex<double>(xmin / 2 + xmax / 2, ymin / 2 + ymax / 2);
    for (size_t j = begin; j < end; ++j)
        node.radius = std::max(node.radius, std::abs(displacement(indices_[j], node.center, node.center_low)));
    for (size_t j = begin; j < end; ++j)
        node.error = std::max(node.error, errors_[indices_[j]]);
    node.moments[0] = static_cast<double>(end - begin);
    build_terms_ += (end - begin) * order;
    if (node.radius > 0)
        for (size_t j = begin; j < end; ++j) {
            const auto scaled = displacement(indices_[j], node.center, node.center_low) / node.radius;
            std::complex<double> power = 1;
            for (size_t k = 1; k < order; ++k) {
                power *= scaled;
                node.moments[k] += power;
            }
        }
    if (end - begin > 12 && node.radius > 0) {
        const size_t middle = begin + (end - begin) / 2;
        const bool real_axis = xmax - xmin >= ymax - ymin;
        std::nth_element(indices_.begin() + static_cast<ptrdiff_t>(begin),
                         indices_.begin() + static_cast<ptrdiff_t>(middle),
                         indices_.begin() + static_cast<ptrdiff_t>(end), [&](size_t a, size_t b) {
                             const auto da = displacement(a, reference, reference_low);
                             const auto db = displacement(b, reference, reference_low);
                             return real_axis ? da.real() < db.real() : da.imag() < db.imag();
                         });
        node.left = build(begin, middle);
        node.right = build(middle, end);
    }
    return index;
}
template <size_t Order, class FarField>
bool InteractionTreeOf<Order, FarField>::visit(size_t index, size_t target, std::complex<double>& far,
                                               std::vector<size_t>& near, double& tail, size_t& far_fields) const {
    const auto& node = nodes_[index];
    const bool contains_self = positions_[target] >= node.begin && positions_[target] < node.end;
    const auto delta = lows_.empty() ? points_[target] - node.center
                                     : displacement(target, node.center, node.center_low);
    const double distance = std::abs(delta);
    // A far field sums the node's points where their binary64 copies are. The
    // copies' rounding (the target's e_t, the node's e_N) changes each term by
    // at most (e_t + e_N) / ((1 - q) distance) of itself, q the acceptance
    // ratio; it may not add more than the series' truncation bound (both
    // relative to the node's sum of magnitudes). Points that round to one
    // position while they stand apart are therefore neighbours: a far field
    // would take them as one point mass and, at a target near them, have no
    // correct digit.
    // History: docs/code-history-20261021.md, interaction_tree.cpp 1.
    const double rounding = node.error + (errors_.empty() ? 0.0 : errors_[target]);
    if (!contains_self && distance > 0 && node.radius <= FarField::acceptance * distance &&
        rounding <= interaction_truncation(order, FarField::acceptance) * (1 - FarField::acceptance) * distance) {
        const auto ratio = node.radius / delta;
        auto sum = node.moments.back();
        for (size_t k = order - 1; k-- > 0;)
            sum = sum * ratio + node.moments[k];
        far += sum / delta;
        ++far_fields;
        const double q = node.radius / distance;
        tail += static_cast<double>(node.end - node.begin) / distance * std::pow(q, order) / (1 - q);
        return std::isfinite(far.real()) && std::isfinite(far.imag()) && std::isfinite(tail);
    }
    if (node.left)
        return visit(node.left, target, far, near, tail, far_fields) &&
               visit(node.right, target, far, near, tail, far_fields);
    for (size_t j = node.begin; j < node.end; ++j)
        if (indices_[j] != target)
            near.push_back(indices_[j]);
    return true;
}
template <size_t Order, class FarField>
bool InteractionTreeOf<Order, FarField>::evaluate(size_t target, std::complex<double>& far,
                                                  std::vector<size_t>& near, double& tail) const {
    size_t far_fields = 0;
    return evaluate(target, far, near, tail, far_fields);
}
template <size_t Order, class FarField>
bool InteractionTreeOf<Order, FarField>::evaluate(size_t target, std::complex<double>& far,
                                                  std::vector<size_t>& near, double& tail,
                                                  size_t& far_fields) const {
    far = 0;
    near.clear();
    tail = 0;
    far_fields = 0;
    return target < points_.size() && visit(0, target, far, near, tail, far_fields);
}
template <size_t Order, class FarField>
size_t InteractionTreeOf<Order, FarField>::move(size_t point, std::complex<double> to) {
    if (point >= points_.size() || nodes_.empty() || to == points_[point])
        return 0;
    const auto from = points_[point];
    const size_t position = positions_[point];
    size_t terms = 0;
    // From the root down to the leaf that holds the point.
    for (size_t index = 0;;) {
        auto& node = nodes_[index];
        const double reach = std::abs(to - node.center);
        if (reach > node.radius) {
            // sum ((z - c) / r')^k = (r / r')^k sum ((z - c) / r)^k. A node of
            // radius zero holds its points at its centre: its moments beyond
            // the count are zero, and stay so.
            const double scale = node.radius / reach;
            double factor = 1;
            for (size_t k = 1; k < order; ++k) {
                factor *= scale;
                node.moments[k] *= factor;
            }
            node.radius = reach;
        }
        // A node still of radius zero keeps the point at its centre.
        if (node.radius > 0) {
            const auto before = (from - node.center) / node.radius, after = (to - node.center) / node.radius;
            std::complex<double> old_power = 1, new_power = 1;
            for (size_t k = 1; k < order; ++k) {
                old_power *= before;
                new_power *= after;
                node.moments[k] += new_power - old_power;
            }
        }
        terms += order;
        if (!node.left)
            break;
        index = position < nodes_[node.left].end ? node.left : node.right;
    }
    points_[point] = to;
    return terms;
}
template class InteractionTreeOf<InteractionTree::order>;
template class InteractionTreeOf<MachineInteractionTree::order>;
// The program route's binary64 stage (program_machine.cpp).
template class InteractionTreeOf<heuristics::program_tree_order, heuristics::ProgramFarField>;
} // namespace certpolysolve::detail
