// The secular refinement's validation: the coordinates of a centred factor (RootCoordinates: a region in the original
// coordinates, its goal there), the exact disk counts of a cluster of nodes (secular_disk_count), the split of a
// component by exact counts (split_by_exact_counts), and the inclusion cover of the secular model at its nodes
// (secular_inclusion_cover), with the realness of its regions (classify_cover_realness). Part of the secular
// refinement (secular_refinement.cpp), split from it verbatim in W8's S5 (docs/w8-plan-20261021.md); state: none.
// Entry points: secular_inclusion_cover (SecularRefinement::step); RootCoordinates (the refinement and factor_state).
#include "methods/secular/secular_refinement.hpp"
#include "core/query.hpp"
#include "methods/secular/secular_internal.hpp"
#include "methods/secular/secular_geometry.hpp"
#include "numeric/evaluate/disk_blocked_polynomial_evaluator.hpp"
#include "numeric/evaluate/disk_sparse_polynomial_evaluator.hpp"
#include "core/dovetail.hpp"
#include "numeric/evaluate/limb_horner_evaluator.hpp"
#include "numeric/evaluate/node_products.hpp"
#include "core/telemetry.hpp"
#include <algorithm>
#include <cmath>
#include <complex>
#include <functional>
#include <limits>
#include <map>
#include <numeric>
#include <optional>

namespace certpolysolve::detail {
namespace {
// A proof sees the complete, immutable geometry. Publication may narrow or
// move regions only after this function returns all classifications.
std::vector<cps_realness> classify_cover_realness(
    std::span<const Root> cover, std::span<const std::vector<size_t>> members, bool real) {
    if (cover.size() != members.size())
        throw std::logic_error("cover membership size mismatch");
    std::vector<cps_realness> result(cover.size(), CPS_REAL_UNKNOWN);
    for (size_t i = 0; i < cover.size(); ++i) {
        if (members[i].empty()) continue;
        const auto &root = cover[i];
        result[i] = arb_contains_zero(acb_imagref(root.box.value())) ? CPS_REAL_UNKNOWN : CPS_NONREAL;
        if (real && members[i].size() == 1 && result[i] == CPS_REAL_UNKNOWN) {
            ComplexBall conjugate;
            acb_conj(conjugate.value(), root.box.value());
            bool other = false;
            for (size_t j = 0; j < cover.size(); ++j)
                if (i != j && !members[j].empty() && acb_overlaps(conjugate.value(), cover[j].box.value()))
                    other = true;
            // The conjugate belongs to this unique-root region: the complete
            // disjoint cover contains every root, and no other region meets it.
            if (!other) result[i] = CPS_REAL;
        }
    }
    return result;
}
} // namespace
bool exact_nodes(std::span<const ComplexBall> nodes) {
    for (const auto &node : nodes)
        if (!acb_is_finite(node.value()) || !acb_is_exact(node.value()))
            return false;
    return true;
}
RootCoordinates::RootCoordinates(GaussianRational offset) : offset_(std::move(offset)) {
    for (bool imaginary : {false, true}) {
        const fmpq *part = imaginary ? offset_.imag.value() : offset_.real.value();
        const auto bits = static_cast<slong>(fmpz_bits(fmpq_numref(part)) + fmpz_bits(fmpq_denref(part))) + 64;
        arb_set_fmpq(imaginary ? acb_imagref(exact_.value()) : acb_realref(exact_.value()), part, bits);
    }
    exact_available_ = acb_is_exact(exact_.value());
}
double RootCoordinates::origin_excess_log2(const acb_t z) const {
    const acb_struct *offset = exact_offset();
    if (!offset || !acb_is_finite(z))
        return 0;
    ComplexBall original;
    acb_add(original.value(), z, offset, ARF_PREC_EXACT);
    Magnitude local, from_origin;
    acb_get_mag(local.value(), z);
    acb_get_mag(from_origin.value(), original.value());
    if (mag_is_zero(local.value()))
        return 0;
    if (mag_is_zero(from_origin.value()))
        return HUGE_VAL;
    return std::max(0.0, mag_get_d_log2_approx(local.value()) - mag_get_d_log2_approx(from_origin.value()));
}
Root RootCoordinates::original(const Root &local, slong precision) const {
    Root root = local;
    if (identity())
        return root;
    ComplexBall offset;
    arb_set_fmpq(acb_realref(offset.value()), offset_.real.value(), precision);
    arb_set_fmpq(acb_imagref(offset.value()), offset_.imag.value(), precision);
    acb_add(root.box.value(), root.box.value(), offset.value(), precision);
    root.certificate += ";exact-rational-variable-shift";
    root.publish(precision);
    if (!fmpq_is_zero(offset_.imag.value()))
        root.realness = arb_is_zero(acb_imagref(root.box.value())) ? CPS_REAL
            : arb_contains_zero(acb_imagref(root.box.value())) ? CPS_REAL_UNKNOWN : CPS_NONREAL;
    return root;
}
bool RootCoordinates::meets(const Root &root, const cps_options &options, slong precision) const {
    if (identity())
        return meets_goal(root, options, precision);
    return meets_goal(original(root, precision), options, precision);
}
std::optional<slong> secular_disk_count(std::span<const ComplexBall> nodes, std::span<const ComplexBall> residues,
                                        std::span<const uint8_t> member, const acb_t centre, const arf_t radius,
                                        DiskCountScratch &scratch) {
    const size_t n = nodes.size();
    if (residues.size() != n || member.size() != n || !exact_nodes(nodes) || !acb_is_exact(centre) ||
        !arf_is_finite(radius) || arf_sgn(radius) <= 0)
        return std::nullopt;
    // Enclosures here need a few correct digits of values that are compared
    // with zero, so one limb suffices; too little precision can only make
    // the count inconclusive. (Subtracting exact numbers rounds the exact
    // difference, so offsets b_j - c keep full relative accuracy.)
    constexpr slong proof_bits = FLINT_BITS;
    // The truncation error of the exterior field, relative to its constant
    // term 1: the resolution of a magnitude bound.
    constexpr slong tail_bits = MAG_BITS;
    auto &inverse = scratch.inverse;
    auto &weight = scratch.weight;
    auto &ratio = scratch.ratio;
    if (inverse.size() < n) {
        inverse.resize(n);
        weight.resize(n);
        ratio.resize(n);
    }
    // The operations are added to the scratch's tally on every return.
    struct Tally {
        double &out;
        double operations = 0;
        ~Tally() { out += operations; }
    } tally{scratch.operations};
    RealBall rho, distance, exterior, doubled;
    arb_set_arf(rho.value(), radius);
    arb_mul_2exp_si(doubled.value(), rho.value(), 1);
    // Terms evaluated directly on each arc: the inside nodes, and the outside
    // nodes within twice the radius. (Expanded, a node at distance rho (1 +
    // eta) needs about 1 / eta terms.) The expansion holds the nodes beyond,
    // whose series gain at least a bit a term. History:
    // docs/code-history-20261021.md, secular_refinement.cpp 14.
    std::vector<ComplexBall> direct_offset, direct_residue;
    std::vector<uint8_t> expanded(n, 0);
    ComplexBall offset;
    // First pass: every node's side of the circle, and for an expanded node
    // 1 / (b_j - c), a_j = |w_j| / |b_j - c| and q_j = rho / |b_j - c| < 1/2.
    for (size_t j = 0; j < n; ++j) {
        if (!acb_is_finite(residues[j].value()))
            return std::nullopt;
        acb_sub(offset.value(), nodes[j].value(), centre, proof_bits);
        acb_abs(distance.value(), offset.value(), proof_bits);
        if (member[j] || !arb_gt(distance.value(), doubled.value())) {
            if (member[j] ? !arb_lt(distance.value(), rho.value()) : !arb_gt(distance.value(), rho.value()))
                return std::nullopt;
            direct_offset.emplace_back();
            acb_swap(direct_offset.back().value(), offset.value());
            direct_residue.emplace_back();
            acb_set_round(direct_residue.back().value(), residues[j].value(), proof_bits);
        } else {
            expanded[j] = 1;
            acb_inv(inverse[j].value(), offset.value(), proof_bits);
            acb_abs(weight[j].value(), residues[j].value(), proof_bits);
            arb_div(weight[j].value(), weight[j].value(), distance.value(), proof_bits);
            arb_add(exterior.value(), exterior.value(), weight[j].value(), proof_bits);
            arb_div(ratio[j].value(), rho.value(), distance.value(), proof_bits);
        }
    }
    const slong m = static_cast<slong>(std::count_if(member.begin(), member.end(), [](uint8_t k) { return k != 0; }));
    const size_t direct = direct_offset.size();
    // The exterior field E(z) = 1 - sum_{j expanded} w_j / (z - b_j) about the
    // centre: E = sum_k e_k (z - c)^k, e_k = [k = 0] + sum_j w_j / (b_j - c)^(k+1).
    // Node j's series converges at the rate q_j, so it is cut after K_j
    // terms, where its remainder on the circle, a_j q_j^K_j / (1 - q_j),
    // is at most a_j 2^-tail_bits / X with X = sum_j a_j: the remainders then
    // sum to 2^-tail_bits at most. K_j is chosen in binary64 (far nodes need
    // few terms); the remainder actually left is summed with enclosures.
    Float bound;
    arb_get_ubound_arf(bound.value(), exterior.value(), proof_bits);
    const double x = arf_get_d(bound.value(), ARF_RND_UP);
    if (!std::isfinite(x))
        return std::nullopt;
    const double target = x > 0 ? -static_cast<double>(tail_bits) - std::log2(x) : 0;
    auto &moments = scratch.moments;
    if (moments.empty())
        moments.emplace_back();
    size_t used = 1, series_terms = 0;
    acb_one(moments[0].value());
    RealBall tail, remainder, complement;
    ComplexBall power;
    for (size_t j = 0; j < n; ++j) {
        if (!expanded[j] || arb_is_zero(weight[j].value()))
            continue;
        arb_get_ubound_arf(bound.value(), ratio[j].value(), proof_bits);
        const double q = arf_get_d(bound.value(), ARF_RND_UP);
        if (!(q > 0 && q < 1))
            return std::nullopt;
        // q^K / (1 - q) <= 2^target: K >= (target + log2(1 - q)) / log2(q).
        const auto terms = static_cast<size_t>(
            std::max(1.0, std::ceil((target + std::log2(1 - q)) / std::log2(q))));
        series_terms += terms;
        tally.operations += static_cast<double>(terms);
        while (used < terms) {
            if (moments.size() == used)
                moments.emplace_back();
            acb_zero(moments[used++].value());
        }
        acb_mul(power.value(), residues[j].value(), inverse[j].value(), proof_bits);
        for (size_t k = 0; k < terms; ++k) {
            acb_add(moments[k].value(), moments[k].value(), power.value(), proof_bits);
            if (k + 1 < terms)
                acb_mul(power.value(), power.value(), inverse[j].value(), proof_bits);
        }
        arb_pow_ui(remainder.value(), ratio[j].value(), terms, proof_bits);
        arb_mul(remainder.value(), remainder.value(), weight[j].value(), proof_bits);
        arb_sub_ui(complement.value(), ratio[j].value(), 1, proof_bits);
        arb_neg(complement.value(), complement.value());
        if (!arb_is_positive(complement.value()))
            return std::nullopt;
        arb_div(remainder.value(), remainder.value(), complement.value(), proof_bits);
        arb_add(tail.value(), tail.value(), remainder.value(), proof_bits);
    }
    Magnitude tail_bound;
    arb_get_mag(tail_bound.value(), tail.value());
    if (!mag_is_finite(tail_bound.value()))
        return std::nullopt;
    // S = E - sum_{j direct} w_j / (z - b_j) on the arc of the circle between
    // the angles 2 pi t / 2^level and 2 pi (t + 1) / 2^level, as a rectangle.
    // The arc's points lie within pi / 2^level of its midpoint's direction.
    // Returns the open half-plane that holds the rectangle (0: Re > 0,
    // 1: Im > 0, 2: Re < 0, 3: Im < 0), or -1 when it may contain zero.
    // With `midpoint`, S at the arc's midpoint alone.
    ComplexBall direction, point, value, term, angle;
    Magnitude half_width;
    const auto half_plane = [&](ulong t, slong level, bool midpoint) {
        arb_set_ui(acb_realref(angle.value()), 2 * t + 1);
        arb_mul_2exp_si(acb_realref(angle.value()), acb_realref(angle.value()), -level);
        arb_zero(acb_imagref(angle.value()));
        acb_exp_pi_i(direction.value(), angle.value(), proof_bits);
        if (!midpoint) {
            mag_const_pi(half_width.value());
            mag_mul_2exp_si(half_width.value(), half_width.value(), -level);
            arb_add_error_mag(acb_realref(direction.value()), half_width.value());
            arb_add_error_mag(acb_imagref(direction.value()), half_width.value());
        }
        acb_mul_arb(point.value(), direction.value(), rho.value(), proof_bits);
        acb_set(value.value(), moments[used - 1].value());
        for (size_t k = used - 1; k-- > 0;) {
            acb_mul(value.value(), value.value(), point.value(), proof_bits);
            acb_add(value.value(), value.value(), moments[k].value(), proof_bits);
        }
        for (size_t j = 0; j < direct; ++j) {
            acb_sub(term.value(), point.value(), direct_offset[j].value(), proof_bits);
            acb_div(term.value(), direct_residue[j].value(), term.value(), proof_bits);
            acb_sub(value.value(), value.value(), term.value(), proof_bits);
        }
        arb_add_error_mag(acb_realref(value.value()), tail_bound.value());
        arb_add_error_mag(acb_imagref(value.value()), tail_bound.value());
        if (!acb_is_finite(value.value())) return -1;
        if (arb_is_positive(acb_realref(value.value()))) return 0;
        if (arb_is_positive(acb_imagref(value.value()))) return 1;
        if (arb_is_negative(acb_realref(value.value()))) return 2;
        if (arb_is_negative(acb_imagref(value.value()))) return 3;
        return -1;
    };
    // Arcs in increasing angle, bisected until each rectangle misses zero.
    // The budget is a cost cap, not a decision. An arc costs about as many
    // operations as the expansion has terms (used) plus one per direct
    // term, and all arcs together may cost what building the expansion did
    // (series_terms operations). The floor keeps small inputs decidable: a
    // group whose residues cancel to 1/60 of their size needs about a
    // thousand arcs, which cost little at any degree. History:
    // docs/code-history-20261021.md, secular_refinement.cpp 15.
    constexpr size_t minimum_arcs = 1024;
    const size_t arc_budget = std::max(minimum_arcs, series_terms / (used + direct));
    // Below this width an arc's endpoints are no longer distinct at proof_bits.
    const slong deepest = proof_bits - 2;
    std::vector<std::pair<ulong, slong>> pending{{3, 2}, {2, 2}, {1, 2}, {0, 2}};
    std::vector<int> planes;
    size_t evaluated = 0;
    while (!pending.empty()) {
        const auto [t, level] = pending.back();
        pending.pop_back();
        if (++evaluated > arc_budget)
            return std::nullopt;
        tally.operations += static_cast<double>(used + direct);
        const int plane = half_plane(t, level, false);
        if (plane >= 0) {
            planes.push_back(plane);
            continue;
        }
        // Both halves end at the arc's midpoint, so their rectangles hold S
        // there: when the arithmetic does not separate that value from zero
        // (a root on the circle, or terms that cancel beyond one limb, as on
        // a circle inside a cluster whose nodes have not arrived), no
        // bisection decides the count, and the arcs left would be spent for
        // nothing.
        tally.operations += static_cast<double>(used + direct);
        if (level >= deepest || half_plane(t, level, true) < 0)
            return std::nullopt;
        pending.emplace_back(2 * t + 1, level + 1);
        pending.emplace_back(2 * t, level + 1);
    }
    // Consecutive rectangles share a value of S, so their half-planes meet:
    // the step between them is a quarter turn at most. The quarter turns
    // around the circle add up to four times the winding number of S, which
    // is (roots of p in the disk) - (poles of S there, the m inside nodes).
    slong quarter_turns = 0;
    for (size_t i = 0; i < planes.size(); ++i) {
        const int step = (planes[(i + 1) % planes.size()] - planes[i] + 4) % 4;
        if (step == 2)
            return std::nullopt; // opposite half-planes cannot share a value
        quarter_turns += step == 3 ? -1 : step;
    }
    if (quarter_turns % 4 != 0)
        return std::nullopt;
    const slong count = m + quarter_turns / 4;
    return count >= 0 ? std::optional<slong>(count) : std::nullopt;
}
std::optional<slong> secular_disk_count(std::span<const ComplexBall> nodes, std::span<const ComplexBall> residues,
                                        std::span<const uint8_t> member, const acb_t centre, const arf_t radius) {
    DiskCountScratch scratch;
    return secular_disk_count(nodes, residues, member, centre, radius, scratch);
}
namespace {
// A component of the inclusion cover that fails the goal, certified piece by
// piece. Its count is proved (each component of the union of the disks
// D(b_j, n |w_j|) holds as many roots as disks), and an exact count of the
// roots in a square about each tight group of its nodes shows where they are:
// secular_disk_count on the circles of radius r and 2r about the group's
// centre, with the group's nodes inside the first and every other node
// outside the second; equal counts put that many roots in the square of
// half-width r between them. A square that meets no other component holds
// only this component's roots (every root lies in some component), so when
// the squares are pairwise disjoint and their counts add up to the
// component's, they hold all of its roots: the component is replaced by them.
//
// Why a component can fail where its groups meet the goal: a surplus node in
// a cluster serves no root of it, so its disk reaches a root the cluster's
// nodes do not serve, and the component spans clusters, while the clusters can
// meet a coarse goal long before the surplus moves.
//
// Linkage is scheduling only (the counts certify): nodes within half the
// goal's allowance at either of them, 2^-(bits+1) |b|, are grouped. Their
// distances come from the nodes rounded to the goal's bits plus 64 (at least
// two limbs): a cluster the goal accepts may be far narrower than binary64's
// spacing at its modulus, where nodes a few units apart in binary64's last
// place would never be linked. Binary64 midpoints only order the sweep. The
// counts cost about n one-limb operations each, so a component is split into
// at most max(split_group_floor, n / split_nodes_per_group) groups: about
// n^2 / 8 operations, below the n^2 of the round's regeneration. History:
// docs/code-history-20261021.md, secular_refinement.cpp 16.
inline constexpr size_t split_group_floor = 8;
inline constexpr size_t split_nodes_per_group = 16;
std::optional<std::vector<Root>> split_by_exact_counts(size_t component, const std::vector<Root> &cover,
                                                       const std::vector<std::vector<size_t>> &members,
                                                       std::span<const ComplexBall> nodes,
                                                       std::span<const ComplexBall> residues, uint64_t multiplicity,
                                                       slong precision, const cps_options &options,
                                                       const RootCoordinates &coordinates, DiskCountScratch &scratch) {
    if (options.goal != CPS_ROOT_APPROXIMATION || options.real_only || multiplicity == 0)
        return std::nullopt;
    const auto &component_nodes = members[component];
    const size_t n = nodes.size(), m = component_nodes.size();
    if (m < 2)
        return std::nullopt;
    const slong link_bits = std::max<slong>(2 * FLINT_BITS, resources::precision_sum(options.bits, 64));
    // The sweep's window in binary64: the allowance at the node plus sixteen
    // units in the last place of its modulus, for the rounding of both nodes.
    constexpr int window_ulps_log2 = std::numeric_limits<double>::digits - 1 - 4;
    std::vector<double> re(m), window(m);
    std::vector<ComplexBall> rounded(m);
    std::vector<RealBall> allowance(m);
    for (size_t k = 0; k < m; ++k) {
        const auto &z = nodes[component_nodes[k]].value();
        re[k] = arf_get_d(arb_midref(acb_realref(z)), ARF_RND_NEAR);
        const double im = arf_get_d(arb_midref(acb_imagref(z)), ARF_RND_NEAR);
        if (!std::isfinite(re[k]) || !std::isfinite(im))
            return std::nullopt;
        acb_set_round(rounded[k].value(), z, link_bits);
        acb_abs(allowance[k].value(), rounded[k].value(), link_bits);
        arb_mul_2exp_si(allowance[k].value(), allowance[k].value(), -static_cast<slong>(options.bits) - 1);
        const double modulus = std::hypot(re[k], im);
        window[k] = std::ldexp(modulus, -static_cast<int>(std::min<uint64_t>(options.bits, 2048)) - 1) +
                    std::ldexp(modulus, -window_ulps_log2);
    }
    std::vector<size_t> order(m), parent(m);
    std::iota(order.begin(), order.end(), size_t{0});
    std::iota(parent.begin(), parent.end(), size_t{0});
    std::sort(order.begin(), order.end(), [&](size_t a, size_t b) { return re[a] < re[b]; });
    const auto find = [&](size_t k) {
        while (parent[k] != k) k = parent[k] = parent[parent[k]];
        return k;
    };
    ComplexBall difference;
    RealBall gap, nearer;
    for (size_t a = 0; a < m; ++a)
        for (size_t b = a + 1; b < m && re[order[b]] - re[order[a]] <= window[order[a]]; ++b) {
            const size_t i = order[a], j = order[b];
            acb_sub(difference.value(), rounded[i].value(), rounded[j].value(), link_bits);
            acb_abs(gap.value(), difference.value(), link_bits);
            arb_min(nearer.value(), allowance[i].value(), allowance[j].value(), link_bits);
            if (arb_le(gap.value(), nearer.value()))
                parent[find(i)] = find(j);
        }
    std::vector<std::vector<size_t>> groups(m);
    for (size_t k = 0; k < m; ++k)
        groups[find(k)].push_back(component_nodes[k]);
    std::erase_if(groups, [](const auto &g) { return g.empty(); });
    if (groups.size() > std::max(split_group_floor, n / split_nodes_per_group) ||
        std::none_of(groups.begin(), groups.end(), [](const auto &g) { return g.size() > 1; }))
        return std::nullopt;
    std::vector<Root> boxes;
    std::vector<uint8_t> member(n, 0);
    uint64_t total = 0;
    ComplexBall centre, offset;
    RealBall distance, outer;
    Float radius, spread, upper;
    Magnitude bound;
    for (const auto &group : groups) {
        // The centre: the members' mean, made exact; the radius: twice the
        // members' largest distance from it (a single node: its own disk's).
        acb_zero(centre.value());
        for (const auto j : group)
            acb_add(centre.value(), centre.value(), nodes[j].value(), precision);
        acb_div_ui(centre.value(), centre.value(), static_cast<ulong>(group.size()), precision);
        acb_get_mid(centre.value(), centre.value());
        arf_zero(spread.value());
        for (const auto j : group) {
            acb_sub(offset.value(), nodes[j].value(), centre.value(), precision);
            acb_get_abs_ubound_arf(upper.value(), offset.value(), precision);
            arf_max(spread.value(), spread.value(), upper.value());
        }
        if (group.size() > 1)
            arf_mul_2exp_si(radius.value(), spread.value(), 1);
        else {
            acb_get_mag(bound.value(), residues[group.front()].value());
            mag_mul_ui(bound.value(), bound.value(), static_cast<ulong>(n));
            arf_set_mag(radius.value(), bound.value());
        }
        arf_set_round(radius.value(), radius.value(), MAG_BITS, ARF_RND_UP);
        if (!arf_is_finite(radius.value()) || arf_sgn(radius.value()) <= 0 ||
            arf_cmp(radius.value(), spread.value()) <= 0)
            return std::nullopt;
        arb_set_arf(outer.value(), radius.value());
        arb_mul_2exp_si(outer.value(), outer.value(), 1);
        for (const auto j : group) member[j] = 1;
        for (size_t j = 0; j < n; ++j) {
            if (member[j]) continue;
            acb_sub(offset.value(), nodes[j].value(), centre.value(), precision);
            acb_abs(distance.value(), offset.value(), precision);
            if (!arb_gt(distance.value(), outer.value())) {
                for (const auto k : group) member[k] = 0;
                return std::nullopt;
            }
        }
        Float doubled;
        arf_mul_2exp_si(doubled.value(), radius.value(), 1);
        const auto inner_count = secular_disk_count(nodes, residues, member, centre.value(), radius.value(), scratch);
        const auto outer_count = secular_disk_count(nodes, residues, member, centre.value(), doubled.value(), scratch);
        for (const auto k : group) member[k] = 0;
        if (!inner_count || !outer_count || *inner_count != *outer_count)
            return std::nullopt;
        if (*inner_count == 0)
            continue;
        const auto roots = static_cast<uint64_t>(*inner_count);
        if (roots > std::numeric_limits<uint64_t>::max() / multiplicity)
            return std::nullopt;
        Root box;
        acb_set(box.box.value(), centre.value());
        acb_add_error_arf(box.box.value(), radius.value());
        box.count = roots * multiplicity;
        box.multiplicity = roots == 1 ? multiplicity : 0;
        box.certificate = "secular-model-count";
        box.publish(precision);
        if (!acb_is_finite(box.box.value()) || !coordinates.meets(box, options, precision))
            return std::nullopt;
        // The published (outward rounded) square must still lie strictly in
        // the outer circle, as certify_cluster_box checks.
        acb_sub(offset.value(), box.box.value(), centre.value(), precision);
        acb_get_abs_ubound_arf(upper.value(), offset.value(), precision);
        if (!arf_is_finite(upper.value()) || arf_cmp(upper.value(), doubled.value()) >= 0)
            return std::nullopt;
        total += box.count;
        boxes.push_back(std::move(box));
    }
    if (total != cover[component].count)
        return std::nullopt;
    for (size_t a = 0; a < boxes.size(); ++a) {
        for (size_t b = a + 1; b < boxes.size(); ++b)
            if (acb_overlaps(boxes[a].box.value(), boxes[b].box.value()))
                return std::nullopt;
        for (size_t j = 0; j < cover.size(); ++j)
            if (j != component && !members[j].empty() && acb_overlaps(boxes[a].box.value(), cover[j].box.value()))
                return std::nullopt;
    }
    return boxes;
}
} // namespace
bool secular_inclusion_cover(std::vector<Root> &out, std::vector<uint8_t> &inactive,
                             std::span<const ComplexBall> nodes, std::span<const ComplexBall> residues,
                             bool real, uint64_t multiplicity, slong precision, const cps_options &options,
                             const RootCoordinates &coordinates, std::span<const Root> kept_roots,
                             std::span<const uint8_t> kept) {
    const size_t n = nodes.size();
    if (residues.size() != n || !exact_nodes(nodes))
        return false;
    const bool keeping = kept.size() == n && kept_roots.size() == n &&
                         std::any_of(kept.begin(), kept.end(), [](uint8_t k) { return k != 0; });
    std::vector<Root> cover(n);
    std::vector<std::vector<size_t>> members(n);
    Magnitude radius;
    for (size_t i = 0; i < n; ++i) {
        if (!acb_is_finite(residues[i].value()))
            return false;
        // det(xI-(diag(b)+a*1^T)) = Q(x)*(1-sum a/(x-b)).
        // Gershgorin row disks are contained in D(b_i,n*|a_i|). Every
        // component of these enlarged disks retains its exact eigenvalue count.
        acb_set(cover[i].box.value(), nodes[i].value());
        acb_get_mag(radius.value(), residues[i].value());
        mag_mul_ui(radius.value(), radius.value(), static_cast<ulong>(n));
        acb_add_error_mag(cover[i].box.value(), radius.value());
        cover[i].count = multiplicity;
        cover[i].multiplicity = multiplicity;
        cover[i].certificate = "secular-identity-gershgorin-cover";
        cover[i].publish(precision);
        members[i].push_back(i);
    }
    // Each node's own disk, before the merging: the rule for kept nodes in a mixed component (assemble) asks which
    // node's disk reached them.
    std::vector<ComplexBall> disks;
    if (keeping) {
        disks.resize(n);
        for (size_t i = 0; i < n; ++i)
            acb_set(disks[i].value(), cover[i].box.value());
    }
    // Merge EXPORTED rectangles to a fixed point, not just the original disks.
    std::vector<OutwardBounds> bounds;
    bounds.reserve(n);
    for (const auto &root : cover)
        bounds.emplace_back(root.box.value());
    bool changed = true;
    while (changed) {
        changed = false;
        for (size_t i = 0; i < cover.size(); ++i)
            if (!members[i].empty())
                for (size_t j = i + 1; j < cover.size(); ++j)
                    if (!members[j].empty() && !bounds[i].apart(bounds[j]) &&
                        acb_overlaps(cover[i].box.value(), cover[j].box.value())) {
                        acb_union(cover[i].box.value(), cover[i].box.value(), cover[j].box.value(), resources::precision_sum(precision, 16));
                        cover[i].count += cover[j].count;
                        cover[i].multiplicity = 0;
                        cover[i].publish(precision);
                        bounds[i] = OutwardBounds(cover[i].box.value());
                        members[i].insert(members[i].end(), members[j].begin(), members[j].end());
                        members[j].clear();
                        changed = true;
                    }
    }
    const auto realness = classify_cover_realness(cover, members, real);
    std::vector<uint8_t> frozen(n, 0);
    std::vector<Root> result;
    result.reserve(cover.size());
    std::vector<uint8_t> met(cover.size(), 0);
    for (size_t i = 0; i < cover.size(); ++i)
        if (!members[i].empty()) {
            auto &root = cover[i];
            root.realness = realness[i];
            if (root.realness == CPS_REAL) {
                arb_zero(acb_imagref(root.box.value()));
                root.publish(precision);
            }
            met[i] = coordinates.meets(root, options, precision);
        }
    // A component whose every member holds a kept certificate takes their
    // regions (the declaration says why that is sound).
    std::vector<uint8_t> replaced(cover.size(), 0);
    if (keeping)
        for (size_t i = 0; i < cover.size(); ++i)
            replaced[i] = !members[i].empty() &&
                          std::all_of(members[i].begin(), members[i].end(), [&](size_t m) { return kept[m] != 0; });
    // Components that fail the goal, split by exact counts on the whole cover
    // (split_by_exact_counts); a split whose squares meet an earlier split's
    // is not taken.
    std::vector<std::optional<std::vector<Root>>> splits(cover.size());
    std::vector<const Root *> taken;
    DiskCountScratch scratch;
    for (size_t i = 0; i < cover.size(); ++i)
        if (!members[i].empty() && !met[i] && !replaced[i])
            if (auto split = split_by_exact_counts(i, cover, members, nodes, residues, multiplicity, precision,
                                                   options, coordinates, scratch)) {
                const bool clear = std::none_of(split->begin(), split->end(), [&](const Root &box) {
                    return std::any_of(taken.begin(), taken.end(), [&](const Root *other) {
                        return acb_overlaps(box.box.value(), other->box.value());
                    });
                });
                if (!clear)
                    continue;
                splits[i] = std::move(split);
                for (const auto &box : *splits[i])
                    taken.push_back(&box);
            }
    // Whether an uncertified member's own disk holds a kept member of component i: a node far from its root, whose
    // wide disk merges kept certificates into its component.
    const auto reached = [&](size_t i) {
        for (const auto a : members[i])
            if (!kept[a])
                for (const auto k : members[i])
                    if (kept[k] && acb_contains(disks[a].value(), nodes[k].value()))
                        return true;
        return false;
    };
    // With kept regions, the cover is tried with them first; if a kept region
    // meets another region it is not used, and the cover is the components'.
    const auto assemble = [&](bool substitute) {
        result.clear();
        std::fill(frozen.begin(), frozen.end(), 0);
        for (size_t i = 0; i < cover.size(); ++i)
            if (!members[i].empty()) {
                if (substitute && replaced[i]) {
                    for (auto member : members[i]) {
                        result.push_back(kept_roots[member]);
                        frozen[member] = 1;
                    }
                    continue;
                }
                if (splits[i] || met[i])
                    for (auto member : members[i])
                        frozen[member] = 1;
                // A component that mixes kept certificates with nodes still uncertified, where an uncertified
                // node's own disk reaches a kept one (reached), keeps the kept nodes where they were certified:
                // only the others iterate, and the component is not restarted as a whole (restart_members takes
                // active nodes only). The cover is unchanged: the component is one region until it splits and its
                // kept members' singletons take their regions above. Before this, two nodes far from their roots
                // (Hermite 64's reduced polynomial: residues to 2^1.3, disks of radius about 80) merged 29 kept
                // certificates into their component, all 31 nodes were restarted at 234 bits and the solve climbed
                // to 335 (24 ms; degree 80 takes 3.6). Where the uncertified nodes' disks are small (a converged
                // close pair, its component formed through the kept node's side), the kept nodes iterate and the
                // component restarts whole, as before: holding them there left a cluster unrestarted and seven
                // fuzz solves stopped at their caps (28 October 2026).
                if (substitute && !replaced[i] && !met[i] && !splits[i] && reached(i))
                    for (auto member : members[i])
                        if (kept[member])
                            frozen[member] = 1;
                if (splits[i])
                    result.insert(result.end(), splits[i]->begin(), splits[i]->end());
                else
                    result.push_back(cover[i]);
            }
    };
    assemble(keeping);
    if (keeping && !disjoint(result))
        assemble(false);
    if (!disjoint(result))
        return false;
    if (!coordinates.identity()) {
        std::vector<Root> exported;
        exported.reserve(result.size());
        for (const auto &root : result)
            exported.push_back(coordinates.original(root, precision));
        // Translation is exact mathematically, but enclosing the rational
        // offset and publishing endpoints can merge close output rectangles.
        // Keep refining until the actual exported cover is disjoint as well.
        if (!disjoint(exported))
            std::fill(frozen.begin(), frozen.end(), 0);
    }
    out = std::move(result);
    inactive = std::move(frozen);
    return true;
}
} // namespace certpolysolve::detail
