// The geometry of a secular model about a centre c, a function of its nodes b_j and residues w_j, for the secular
// refinement's restarts: isolating plateaus of exact disk counts (isolating_plateau), the mean of the roots inside one
// (root_mean), the roots about a tight group of nodes (counted_group), the span of a restart's hints (hint_span_log2).
//
// For monic p the model is p = Q S, Q = prod (z - b_j), S = 1 - sum w_j / (z - b_j); every count is
// secular_disk_count's, the roots in an open disk by the argument principle on S, exact when it decides. Circles keep
// away from the nodes: widest_gap_log2 picks the radius in a log2 interval farthest from every log2 |b_j - c|
// (node_distances_log2). isolating_plateau walks such circles down by halves from 2^top_log2 (ladder_radius_log2, at
// most FLINT_BITS rungs) to the outermost two consecutive ones with one count K >= 2, extended inward while it holds:
// the annulus between holds no root, and the inner disk exactly K. root_mean estimates their mean by the trapezoid rule
// for the integral of (z - c) p'/p over a circle in the annulus' middle half, divided by 2 pi i, with p'/p = Q'/Q +
// S'/S, doubling the points from 16 to binary64's resolution and checking the zeroth moment against K: a centre for a
// local model, not an enclosure. counted_group takes equal counts on the circles of radius r and 2r about the group's
// centroid (r twice the members' largest distance from it, every other node beyond 2r), else the isolating plateau from
// half the distance to the nearest other node. hint_span_log2: the bits from the hints' largest magnitude to their
// closest pair.
//
// State: none; operation counts return to the caller's meter (DiskCountScratch::operations, root_mean's operations).
// Entry points, all from secular_refinement.cpp: node_distances_log2, isolating_plateau, root_mean, set_radius_log2
// from SecularRefinement::group_model, counted_group from SecularRefinement::regroup_restart, hint_span_log2 from
// SecularRefinement::geometric_restart. widest_gap_log2 and ladder_radius_log2 are used only here.
// Not here: the count itself (secular_disk_count, DiskCountScratch: secular_disk_count.hpp), split_by_exact_counts,
// locate_mean and the restarts (secular_refinement.cpp).
// Design notes: docs/clustered-roots-analysis-20261003.md (re-seeding a cluster; certifying one by an exact count).
#include "methods/secular/secular_geometry.hpp"
#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>

namespace certpolysolve::detail {
double widest_gap_log2(std::span<const double> distances_log2, double low, double high) {
    const double target_log2 = (low + high) / 2;
    const auto gap = [&](double x) {
        const auto above = std::lower_bound(distances_log2.begin(), distances_log2.end(), x);
        double g = std::numeric_limits<double>::infinity();
        if (above != distances_log2.end()) g = std::min(g, *above - x);
        if (above != distances_log2.begin()) g = std::min(g, x - *(above - 1));
        return g;
    };
    double best = target_log2, best_gap = gap(target_log2);
    const auto consider = [&](double x) {
        x = std::clamp(x, low, high);
        const double g = gap(x);
        if (g > best_gap) {
            best = x;
            best_gap = g;
        }
    };
    consider(low);
    consider(high);
    auto first = std::lower_bound(distances_log2.begin(), distances_log2.end(), low);
    if (first != distances_log2.begin()) --first;
    for (auto it = first; it != distances_log2.end() && *it <= high; ++it)
        if (it + 1 != distances_log2.end())
            consider((*it + *(it + 1)) / 2);
    return best;
}
double ladder_radius_log2(std::span<const double> distances_log2, double target_log2) {
    return widest_gap_log2(distances_log2, target_log2 - 0.5, target_log2 + 0.5);
}
void set_radius_log2(arf_t radius, double radius_log2) {
    const double whole = std::floor(radius_log2);
    arf_set_d(radius, std::exp2(radius_log2 - whole));
    arf_mul_2exp_si(radius, radius, static_cast<slong>(whole));
}
std::vector<double> node_distances_log2(std::span<const ComplexBall> nodes, const acb_t centre) {
    std::vector<double> distance(nodes.size());
    ComplexBall difference;
    Magnitude size;
    for (size_t j = 0; j < nodes.size(); ++j) {
        acb_sub(difference.value(), nodes[j].value(), centre, FLINT_BITS);
        acb_get_mag(size.value(), difference.value());
        distance[j] = mag_is_zero(size.value()) ? -std::numeric_limits<double>::infinity()
                                                : mag_get_d_log2_approx(size.value());
    }
    return distance;
}
std::optional<Plateau> isolating_plateau(std::span<const ComplexBall> nodes, std::span<const ComplexBall> residues,
                                         const acb_t centre, std::span<const double> distance,
                                         std::span<const double> sorted, double top_log2, DiskCountScratch &scratch) {
    const size_t n = nodes.size();
    std::vector<uint8_t> member(n);
    Float radius;
    std::optional<Plateau> plateau;
    std::optional<slong> previous;
    double previous_log2 = 0;
    bool decided = false;
    int undecided_run = 0;
    for (slong rung = 0; rung < FLINT_BITS; ++rung) {
        const double radius_log2 = ladder_radius_log2(sorted, top_log2 - static_cast<double>(rung));
        for (size_t j = 0; j < n; ++j) member[j] = distance[j] < radius_log2;
        set_radius_log2(radius.value(), radius_log2);
        const auto count = secular_disk_count(nodes, residues, member, centre, radius.value(), scratch);
        if (!count) {
            if (plateau || decided || ++undecided_run == 2)
                return plateau;
            continue;
        }
        decided = true;
        if (plateau) {
            if (*count != plateau->count)
                return plateau;
            plateau->inner_log2 = radius_log2;
            continue;
        }
        if (*count <= 1)
            return std::nullopt;
        if (previous && *previous == *count)
            plateau = Plateau{*count, radius_log2, previous_log2};
        previous = count;
        previous_log2 = radius_log2;
    }
    return plateau;
}
bool root_mean(acb_t mean, std::span<const ComplexBall> nodes, std::span<const ComplexBall> residues,
               const acb_t centre, const Plateau &plateau, std::span<const double> distances_log2,
               double &operations, uint64_t maximum_points) {
    constexpr slong precision = FLINT_BITS;
    constexpr int resolution_bits = std::numeric_limits<double>::digits;
    constexpr slong first_points = 16;
    // In the annulus' middle half, so that both of the rule's rates are a
    // quarter of its width at least.
    const double width_log2 = plateau.outer_log2 - plateau.inner_log2;
    const double circle_log2 = widest_gap_log2(distances_log2, plateau.inner_log2 + width_log2 / 4,
                                               plateau.outer_log2 - width_log2 / 4);
    const double rate_log2 = std::min(circle_log2 - plateau.inner_log2, plateau.outer_log2 - circle_log2);
    if (!(rate_log2 > 0))
        return false;
    slong most = first_points;
    while (static_cast<double>(most) * rate_log2 < resolution_bits) {
        // A nearly collapsed annulus can ask for arbitrarily many nodes.
        // This mean is only a restart hint: decline before exceeding its
        // admitted quadratic shell or overflowing the signed point count.
        if (static_cast<uint64_t>(most) > maximum_points / 2 || most > WORD_MAX / 2)
            return false;
        most *= 2;
    }
    if (static_cast<uint64_t>(most) > maximum_points) return false;
    Float radius;
    set_radius_log2(radius.value(), circle_log2);
    RealBall rho;
    arb_set_arf(rho.value(), radius.value());
    ComplexBall zeroth, first, angle, offset, z, inverse, s, slope, poles, term, value;
    const auto add_point = [&](slong numerator, slong points) {
        acb_set_si(angle.value(), 2 * numerator);
        acb_div_si(angle.value(), angle.value(), points, precision);
        acb_exp_pi_i(offset.value(), angle.value(), precision);
        acb_mul_arb(offset.value(), offset.value(), rho.value(), precision);
        acb_add(z.value(), centre, offset.value(), precision);
        acb_one(s.value());
        acb_zero(slope.value());
        acb_zero(poles.value());
        // Six one-limb operations a node: the inverse, three products, two sums.
        constexpr double operations_a_node = 6;
        operations += operations_a_node * static_cast<double>(nodes.size());
        for (size_t j = 0; j < nodes.size(); ++j) {
            acb_sub(inverse.value(), z.value(), nodes[j].value(), precision);
            acb_inv(inverse.value(), inverse.value(), precision);
            acb_add(poles.value(), poles.value(), inverse.value(), precision);
            acb_mul(term.value(), residues[j].value(), inverse.value(), precision);
            acb_sub(s.value(), s.value(), term.value(), precision);
            acb_mul(term.value(), term.value(), inverse.value(), precision);
            acb_add(slope.value(), slope.value(), term.value(), precision);
        }
        acb_div(value.value(), slope.value(), s.value(), precision);
        acb_add(value.value(), value.value(), poles.value(), precision);
        acb_mul(value.value(), value.value(), offset.value(), precision);
        acb_add(zeroth.value(), zeroth.value(), value.value(), precision);
        acb_mul(value.value(), value.value(), offset.value(), precision);
        acb_add(first.value(), first.value(), value.value(), precision);
    };
    slong points = first_points;
    for (slong l = 0; l < points; ++l)
        add_point(l, points);
    ComplexBall estimate, last, change;
    acb_div_si(last.value(), first.value(), points, precision);
    Magnitude tolerance, moved, floor;
    arb_get_mag(tolerance.value(), rho.value());
    mag_mul_ui(tolerance.value(), tolerance.value(), static_cast<ulong>(plateau.count));
    mag_mul_2exp_si(tolerance.value(), tolerance.value(), -resolution_bits);
    while (points < most) {
        for (slong l = 0; l < points; ++l)
            add_point(2 * l + 1, 2 * points);
        points *= 2;
        acb_div_si(estimate.value(), first.value(), points, precision);
        acb_sub(change.value(), estimate.value(), last.value(), precision);
        acb_get_mid(change.value(), change.value());
        acb_get_mag(moved.value(), change.value());
        // The arithmetic's floor: twice the estimate's own rounding radius.
        mag_max(floor.value(), arb_radref(acb_realref(estimate.value())), arb_radref(acb_imagref(estimate.value())));
        mag_mul_2exp_si(floor.value(), floor.value(), 1);
        mag_max(floor.value(), floor.value(), tolerance.value());
        acb_swap(last.value(), estimate.value());
        if (mag_cmp(moved.value(), floor.value()) <= 0)
            break;
    }
    // The zeroth moment is the count; a rule a quarter off has not converged.
    acb_div_si(zeroth.value(), zeroth.value(), points, precision);
    acb_get_mid(zeroth.value(), zeroth.value());
    acb_sub_si(zeroth.value(), zeroth.value(), plateau.count, precision);
    Magnitude miss, quarter;
    acb_get_mag(miss.value(), zeroth.value());
    mag_set_ui_2exp_si(quarter.value(), 1, -2);
    if (!mag_is_finite(miss.value()) || mag_cmp(miss.value(), quarter.value()) > 0)
        return false;
    acb_get_mid(last.value(), last.value());
    acb_div_si(last.value(), last.value(), plateau.count, precision);
    acb_get_mag(moved.value(), last.value());
    Float inner;
    set_radius_log2(inner.value(), plateau.inner_log2);
    Magnitude inner_bound;
    arf_get_mag_lower(inner_bound.value(), inner.value());
    if (!mag_is_finite(moved.value()) || mag_cmp(moved.value(), inner_bound.value()) >= 0)
        return false;
    acb_add(mean, centre, last.value(), precision);
    acb_get_mid(mean, mean);
    return true;
}

double hint_span_log2(std::span<const ComplexBall> hints, const acb_t centre, slong precision) {
    const size_t k = hints.size();
    if (k < 2)
        return 0;
    std::vector<ComplexBall> offsets(k);
    Magnitude spread, size;
    for (size_t j = 0; j < k; ++j) {
        acb_sub(offsets[j].value(), hints[j].value(), centre, precision);
        acb_get_mid(offsets[j].value(), offsets[j].value());
        Magnitude d;
        acb_get_mag(d.value(), offsets[j].value());
        mag_max(spread.value(), spread.value(), d.value());
    }
    if (mag_is_zero(spread.value()) || !mag_is_finite(spread.value()))
        return 0;
    const slong scale = static_cast<slong>(std::floor(mag_get_d_log2_approx(spread.value())));
    std::vector<std::complex<double>> u(k);
    for (size_t j = 0; j < k; ++j) {
        acb_mul_2exp_si(offsets[j].value(), offsets[j].value(), -scale);
        u[j] = {arf_get_d(arb_midref(acb_realref(offsets[j].value())), ARF_RND_NEAR),
                arf_get_d(arb_midref(acb_imagref(offsets[j].value())), ARF_RND_NEAR)};
    }
    double closest = std::numeric_limits<double>::infinity();
    for (size_t a = 0; a < k; ++a)
        for (size_t b = a + 1; b < k; ++b)
            closest = std::min(closest, std::abs(u[a] - u[b]));
    if (!(closest > 0) || !std::isfinite(closest))
        return 0;
    acb_get_mag(size.value(), centre);
    mag_add(size.value(), size.value(), spread.value());
    return std::max(0.0, mag_get_d_log2_approx(size.value()) - (static_cast<double>(scale) + std::log2(closest)));
}
std::optional<CountedGroup> counted_group(const std::vector<size_t> &group, std::span<const ComplexBall> nodes,
                                          std::span<const ComplexBall> residues, std::vector<uint8_t> &member,
                                          slong precision, DiskCountScratch &scratch) {
    CountedGroup counted;
    ComplexBall offset;
    RealBall distance, outer;
    Float spread, upper, doubled;
    const auto radius = arb_midref(counted.radius.value());
    for (const auto j : group)
        acb_add(counted.centre.value(), counted.centre.value(), nodes[j].value(), precision);
    acb_div_ui(counted.centre.value(), counted.centre.value(), static_cast<ulong>(group.size()), precision);
    acb_get_mid(counted.centre.value(), counted.centre.value());
    for (const auto j : group) {
        acb_sub(offset.value(), nodes[j].value(), counted.centre.value(), precision);
        acb_get_abs_ubound_arf(upper.value(), offset.value(), precision);
        arf_max(spread.value(), spread.value(), upper.value());
    }
    arf_mul_2exp_si(radius, spread.value(), 1);
    arf_set_round(radius, radius, MAG_BITS, ARF_RND_UP);
    if (!arf_is_finite(radius) || arf_sgn(radius) <= 0 || arf_cmp(radius, spread.value()) <= 0)
        return std::nullopt;
    arb_set_arf(outer.value(), radius);
    arb_mul_2exp_si(outer.value(), outer.value(), 1);
    for (const auto j : group)
        member[j] = 1;
    bool apart = true;
    for (size_t j = 0; j < nodes.size() && apart; ++j)
        if (!member[j]) {
            acb_sub(offset.value(), nodes[j].value(), counted.centre.value(), precision);
            acb_abs(distance.value(), offset.value(), precision);
            apart = arb_gt(distance.value(), outer.value());
        }
    std::optional<slong> inner_count, outer_count;
    if (apart) {
        arf_mul_2exp_si(doubled.value(), radius, 1);
        inner_count = secular_disk_count(nodes, residues, member, counted.centre.value(), radius, scratch);
        outer_count = secular_disk_count(nodes, residues, member, counted.centre.value(), doubled.value(), scratch);
    }
    if (inner_count && outer_count && *inner_count == *outer_count) {
        for (const auto j : group)
            member[j] = 0;
        counted.count = *inner_count;
        return counted;
    }
    // The plateau, from half the distance to the nearest node outside the group.
    const auto distances = node_distances_log2(nodes, counted.centre.value());
    double nearest = std::numeric_limits<double>::infinity();
    for (size_t j = 0; j < nodes.size(); ++j)
        if (!member[j])
            nearest = std::min(nearest, distances[j]);
    for (const auto j : group)
        member[j] = 0;
    if (!std::isfinite(nearest))
        return std::nullopt;
    auto sorted = distances;
    std::sort(sorted.begin(), sorted.end());
    const auto plateau = isolating_plateau(nodes, residues, counted.centre.value(), distances, sorted, nearest - 1, scratch);
    if (!plateau)
        return std::nullopt;
    // Every node of the group within the plateau's outer circle, and no other: the count is this group's.
    for (const auto j : group)
        if (!(distances[j] < plateau->outer_log2))
            return std::nullopt;
    Float inner;
    set_radius_log2(inner.value(), plateau->inner_log2);
    arf_max(radius, radius, inner.value());
    arf_set_round(radius, radius, MAG_BITS, ARF_RND_UP);
    counted.count = plateau->count;
    return counted;
}
} // namespace certpolysolve::detail
