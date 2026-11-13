// Certification of a program's roots from its own evaluations in balls (declared in program_certificate.hpp): the
// separated-roots pass (SeparatedRootsPass), which certifies a complete disjoint cover of Newton inclusion disks, with
// counted clusters where disks keep overlapping; the root count on a circle by the argument principle
// (certified_root_count) that the clusters rest on.
//
// certified_root_count cuts the circle into arcs with dyadic ends (ticks of the circle's parameter, circle_point) and
// reads each end's value once (value_at: its argument must be known to a few bits, or the count ends). On an arc, with
// B a ball about its midpoint z_m that holds it and h half its length, the image lies in the disk about p(z_m) of
// radius |p'(B)| h, in the mean-value box p(z_m) + p'(B) (B - z_m) and in p(B); one of them that misses 0 keeps the
// argument's change below a half-turn, so the principal value Arg(p(z_b) conj(p(z_a))) is the change, accepted when
// placed within 2^-3 radians. An arc not resolved is bisected, down to the smallest arc. The changes, summed in balls
// and divided by 2 pi, hold a unique integer: the count.
//
// SeparatedRootsPass::advance, one round. Each point whose disk is not current is evaluated once and its disk formed,
// n |p| / |p'| from an upper bound over a lower bound; a derivative without correct bits doubles the precision, to the
// cap. The overlap sweep (culprits) finds the overlapping pairs from binary64 squares, as guidance only. The pass
// certifies when the n regions meet the goal and are pairwise disjoint (disjoint): n disjoint disks of a polynomial of
// degree n, each holding a root, hold one each. Otherwise the cluster step (cluster_cover) runs once when only overlaps
// remain, again at the first round without progress, and wherever the pass would fail; then the steps: Aberth's
// N / (1 - N S) for a point whose disk overlaps (its sum in balls among few such points, in binary64 otherwise, by the
// interaction tree from machine_interaction_tree_degree points), Newton's for one that is only too wide. A lost point
// (lost_point) iterates alone toward a root no point holds, and a culprit without progress is nudged off its root and
// re-homed the same way. A value at its rounding asks for a raise of precision; while disks overlap, a raise waits for
// progress since the last one.
//
// cluster_cover groups the overlapping points into components (union-find), counts each component's circle about its
// centroid (certified_root_count, on the undivided polynomial where the circle lies far from zero), shrinks each circle
// toward the goal (the precision doubling while its circle is unreadable), and certifies the separated disks and the
// counted clusters when their counts add up to the degree and the cover is disjoint. Counted clusters answer the
// approximation goal over the complex roots and the cluster-radius goal only.
//
// State: the pass keeps its points, disks, precision and round counters across advance calls (a cancelled round is done
// again on the next call); the counts keep nothing.
// Entry points: SeparatedRootsPass (constructed from the points; advance, cover), from solver_program.cpp
// (program_advance: the Newton-disk pass between the preliminary phase and the refinement). certified_root_count is the
// cluster step's, and tests/test_program.cpp calls it.
// Not here: the evaluations (ProgramSource, program_source.cpp), the refinement the route falls back to
// (SecularRefinement), the publication of the cover (program_advance in solver_program.cpp).
// Design notes: docs/slp-second-version-note-20261017.md (section 2, the certificate; section 7, phase 2a as built;
// sections 8 and 9, phase 2b).
#include "certify/attempts/separated_roots.hpp"
#include "numeric/evaluate/monic_source.hpp"
#include "certify/winding_count.hpp"
#include "core/heuristics.hpp"
#include "numeric/candidates/interaction_tree.hpp"
#include "core/resource_limits.hpp"
#include "core/telemetry.hpp"
#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <memory>
#include <numeric>

namespace certpolysolve::detail {

namespace {
// Guidance only, never a certificate (the pass certifies by disjoint() over the published regions): the points
// whose disk overlaps a disk no larger than its own, the larger of each overlapping pair, from double
// approximations of the squares. A disk wider than a sixteenth of the points' span is marked outright and left out
// of the sweep, which keeps the sweep's active list short. Points not representable in binary64 are marked.
struct Overlaps {
    std::vector<uint8_t> culprit;  // the larger disk of an overlapping pair, or a disk marked outright
    std::vector<uint8_t> involved; // every point of an overlapping pair (each gets a step)
    double widest = -std::numeric_limits<double>::infinity(); // log2 of the widest culprit disk relative to its point
    std::vector<std::pair<size_t, size_t>> pairs; // the overlapping pairs the sweep found (the cluster step's components)
};
Overlaps culprits(const std::vector<ComplexBall> &points, const std::vector<Root> &regions) {
    const size_t n = points.size();
    std::vector<uint8_t> marked(n, 0), involved(n, 0);
    struct Square {
        double x, y, r;
        size_t i;
    };
    std::vector<Square> squares;
    squares.reserve(n);
    bool unrepresentable = false;
    double low_x = std::numeric_limits<double>::infinity(), high_x = -low_x, low_y = low_x, high_y = -low_x;
    for (size_t i = 0; i < n; ++i) {
        const double x = arf_get_d(arb_midref(acb_realref(points[i].value())), ARF_RND_NEAR);
        const double y = arf_get_d(arb_midref(acb_imagref(points[i].value())), ARF_RND_NEAR);
        // The published square's half-width: its real radius about its own centre, at least the disk's.
        const double r = mag_get_d(arb_radref(acb_realref(regions[i].box.value()))) +
                         std::abs(arf_get_d(arb_midref(acb_realref(regions[i].box.value())), ARF_RND_NEAR) - x);
        if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(r)) {
            marked[i] = 1;
            unrepresentable = true;
            continue;
        }
        squares.push_back({x, y, r, i});
        low_x = std::min(low_x, x), high_x = std::max(high_x, x);
        low_y = std::min(low_y, y), high_y = std::max(high_y, y);
    }
    const double span = std::max(high_x - low_x, high_y - low_y);
    std::vector<Square> sweep;
    sweep.reserve(squares.size());
    for (const auto &q : squares) {
        if (squares.size() > 16 && q.r > span / 16)
            marked[q.i] = 1;
        else
            sweep.push_back(q);
    }
    std::sort(sweep.begin(), sweep.end(), [](const Square &a, const Square &b) { return a.x - a.r < b.x - b.r; });
    std::vector<size_t> active;
    std::vector<std::pair<size_t, size_t>> pairs;
    for (size_t k = 0; k < sweep.size(); ++k) {
        const auto &q = sweep[k];
        active.erase(std::remove_if(active.begin(), active.end(),
                                    [&](size_t a) { return sweep[a].x + sweep[a].r < q.x - q.r; }),
                     active.end());
        for (const size_t a : active) {
            const auto &p = sweep[a];
            const double reach = (p.r + q.r) * (1 + 0x1p-40);
            if (std::abs(p.x - q.x) <= reach && std::abs(p.y - q.y) <= reach) {
                if (p.r >= q.r) marked[p.i] = 1;
                if (q.r >= p.r) marked[q.i] = 1;
                involved[p.i] = involved[q.i] = 1;
                pairs.emplace_back(std::min(p.i, q.i), std::max(p.i, q.i));
            }
        }
        active.push_back(k);
    }
    // Every disk omitted from squares was marked above. Remember that fact
    // instead of searching all squares again for every marked disk.
    double widest = unrepresentable ? std::numeric_limits<double>::infinity()
                                   : -std::numeric_limits<double>::infinity();
    for (const auto &q : squares)
        if (marked[q.i])
            widest = std::max(widest, std::log2(std::max(q.r, 0x1p-1074)) -
                                          std::log2(std::max(std::hypot(q.x, q.y), 0x1p-1074)));
    for (size_t i = 0; i < n; ++i)
        involved[i] = involved[i] || marked[i];
    return {std::move(marked), std::move(involved), widest, std::move(pairs)};
}

// A culprit far from converged: its disk wider than a 2^-20 of its point (a lost point the preliminary phase
// handed over), which no cluster holds.
bool lost_point(const ComplexBall &point, const Root &region) {
    Magnitude radius, size;
    mag_set(radius.value(), arb_radref(acb_realref(region.box.value())));
    acb_get_mag(size.value(), point.value());
    mag_mul_2exp_si(size.value(), size.value(), -20);
    return mag_cmp(radius.value(), size.value()) > 0;
}

} // namespace

SeparatedRootsPass::SeparatedRootsPass(acb_srcptr points, size_t n, slong precision)
    : points_(n), corrections_(n), regions_(n), current_(n, 0), limited_(n, 0), unreadable_(n, 0),
      precision_(precision) {
    for (size_t i = 0; i < n; ++i)
        acb_get_mid(points_[i].value(), points + static_cast<slong>(i));
}

SeparatedRootsPass::State SeparatedRootsPass::advance(MonicSource &source, const cps_options &options,
                                                      certpolysolve::detail::CandidateStatistics &statistics,
                                                      MonicSource *undivided) {
    if (state_ != State::Running)
        return state_;
    const size_t n = points_.size();
    const auto cap = static_cast<slong>(options.max_work_bits);
    if (n == 0 || static_cast<slong>(n) != source.degree()) {
        state_ = State::Failed;
        failure_ = "the points do not match the degree";
        return state_;
    }
    StageTimer timer(statistics.validation_ns);
    // The disks of the points whose disk is not current.
    ComplexBall value, derivative;
    Magnitude upper_value, lower_value, lower_derivative, radius;
    for (size_t i = 0; i < n; ++i) {
        if (current_[i])
            continue;
        if (options.cancelled())
            return state_; // the round again on the next call: the disks done so far stay current
        ++evaluations_;
        ++statistics.evaluations;
        statistics.evaluation_terms += source.terms();
        unreadable_[i] = !source.evaluate(value.value(), derivative.value(), points_[i].value(), precision_) ||
                         !acb_is_finite(value.value()) || acb_contains_zero(derivative.value());
        Root &region = regions_[i];
        region = Root{};
        if (unreadable_[i]) {
            current_[i] = 1;
            continue;
        }
        // n |p| / |p'|, rounded up: an upper bound of |p| over a lower bound of |p'|.
        acb_get_mag(upper_value.value(), value.value());
        acb_get_mag_lower(lower_derivative.value(), derivative.value());
        mag_div(radius.value(), upper_value.value(), lower_derivative.value());
        mag_mul_ui(radius.value(), radius.value(), static_cast<ulong>(n));
        acb_set(region.box.value(), points_[i].value());
        arb_add_error_mag(acb_realref(region.box.value()), radius.value());
        arb_add_error_mag(acb_imagref(region.box.value()), radius.value());
        region.count = region.multiplicity = 1;
        region.publish(precision_);
        region.realness = arb_contains_zero(acb_imagref(region.box.value())) ? CPS_REAL_UNKNOWN : CPS_NONREAL;
        // The Newton correction, for a step if the disk overlaps or is too wide; the value at its rounding (its
        // ball's lower bound under half its upper bound) makes that step noise.
        acb_get_mag_lower(lower_value.value(), value.value());
        mag_mul_2exp_si(lower_value.value(), lower_value.value(), 1);
        limited_[i] = mag_cmp(lower_value.value(), upper_value.value()) < 0;
        acb_div(corrections_[i].value(), value.value(), derivative.value(), precision_);
        current_[i] = 1;
    }
    ++rounds_;
    // A derivative without correct bits: the precision doubles, to the cap, and those points are evaluated again.
    const size_t unreadable = static_cast<size_t>(std::count(unreadable_.begin(), unreadable_.end(), uint8_t{1}));
    if (unreadable) {
        if (precision_ >= cap) {
            state_ = State::Failed;
            failure_ = "a derivative without correct bits at the precision cap";
            return state_;
        }
        precision_ = std::min(cap, resources::precision_product(precision_, 2));
        ++raises_;
        for (size_t i = 0; i < n; ++i)
            if (unreadable_[i])
                current_[i] = 0;
        return state_;
    }
    const auto overlaps = culprits(points_, regions_);
    const auto &conflict = overlaps.involved;
    overlapping_ = static_cast<size_t>(std::count(overlaps.culprit.begin(), overlaps.culprit.end(), uint8_t{1}));
    wide_ = 0;
    std::vector<uint8_t> step(n, 0);
    for (size_t i = 0; i < n; ++i) {
        const bool wide = !meets_goal(regions_[i], options, precision_);
        wide_ += wide;
        step[i] = conflict[i] || wide;
    }
    if (!overlapping_ && !wide_ && disjoint(regions_)) {
        for (auto &region : regions_)
            region.certificate = "original-source Newton inclusion disk in a complete disjoint cover";
        cover_ = regions_;
        clusters_ = 0;
        state_ = State::Certified;
        return state_;
    }
    // Every region meets the goal and only overlaps remain: points converged on a multiple root or on a cluster
    // tighter than the goal (or two points on one simple root, which the counts tell apart). The cluster step once.
    if (overlapping_ && !wide_ && !converged_tried_) {
        const auto early = cluster_cover(source, undivided, options, statistics, overlaps.pairs, overlaps.involved);
        if (early == ClusterOutcome::Certified)
            state_ = State::Certified;
        if (early != ClusterOutcome::Cancelled)
            converged_tried_ = true;
        if (early != ClusterOutcome::None)
            return state_; // certified, or cancelled (the round again on the next call)
    }
    // Where the pass would fail with overlapping disks, the cluster step first (phase 2b): counted clusters in a
    // complete disjoint cover, under the goals that admit them.
    const auto fail = [&](const char *reason) {
        switch (cluster_cover(source, undivided, options, statistics, overlaps.pairs, overlaps.involved)) {
        case ClusterOutcome::Certified:
            state_ = State::Certified;
            break;
        case ClusterOutcome::Cancelled:
            break; // the round again on the next call
        case ClusterOutcome::None:
            state_ = State::Failed;
            failure_ = reason;
            break;
        }
        return state_;
    };
    // The culprits must grow fewer: a cluster, or two points at one root that Aberth's repulsion does not part,
    // keeps them.
    // Progress: fewer culprits, or the widest culprit disk a bit narrower.
    bool stalled_raise = false;
    const bool narrower = overlaps.widest < best_widest_ - 1;
    if (narrower)
        best_widest_ = overlaps.widest;
    if (overlapping_ < best_overlapping_ || narrower) {
        best_overlapping_ = std::min(best_overlapping_, overlapping_);
        stalled_rounds_ = 0;
        raised_since_progress_ = false;
    } else if (overlapping_) {
        ++stalled_rounds_;
        if (stalled_rounds_ >= heuristics::program_disk_stalled_rounds)
            return fail("overlapping disks that do not separate");
        // At the first round without progress, before re-homing starts, the cluster step once: two points at a
        // multiple root stall there (Aberth's rate is linear), and re-homing cannot help them (no root is free),
        // while two points at one simple root leave the counts short of the degree, and re-homing then goes on.
        if (stalled_rounds_ == 1) {
            const auto early = cluster_cover(source, undivided, options, statistics, overlaps.pairs, overlaps.involved);
            if (early == ClusterOutcome::Certified)
                state_ = State::Certified;
            if (early != ClusterOutcome::None)
                return state_; // certified, or cancelled (the round again on the next call)
        }
        // Halfway, the precision doubles once: the disks of a close pair are n |p/p'| wide at least, which the
        // precision's own rounding bounds from below.
        stalled_raise = stalled_rounds_ == heuristics::program_disk_stalled_rounds / 2 && precision_ < cap;
    }
    if (rounds_ >= heuristics::program_disk_rounds)
        return fail("the rounds allowed");
    // The steps where they are not noise (the precision doubles where they are): for a point whose disk overlaps
    // another, Aberth's correction N / (1 - N S) with N = p / p' and S the sum of 1 / (z - z_j) over the other
    // points, whose repulsion keeps two points from one root; for a point that is only too wide, Newton's. S is a
    // proposal's direction, taken in binary64 from the points' midpoints (the interaction tree from
    // machine_interaction_tree_degree points on), and the step itself in balls at the pass's precision.
    std::vector<std::complex<double>> positions(n);
    bool representable = true;
    for (size_t i = 0; i < n; ++i) {
        positions[i] = {arf_get_d(arb_midref(acb_realref(points_[i].value())), ARF_RND_NEAR),
                        arf_get_d(arb_midref(acb_imagref(points_[i].value())), ARF_RND_NEAR)};
        representable = representable && std::isfinite(positions[i].real()) && std::isfinite(positions[i].imag());
    }
    std::unique_ptr<MachineInteractionTree> tree;
    if (representable && overlapping_ && n >= heuristics::machine_interaction_tree_degree)
        tree = std::make_unique<MachineInteractionTree>(positions);
    std::vector<size_t> near;
    const auto interaction = [&](size_t i) {
        std::complex<double> sum = 0;
        if (tree) {
            std::complex<double> far;
            double tail = 0;
            if (tree->evaluate(i, far, near, tail)) {
                sum = far;
                for (const size_t j : near)
                    if (positions[j] != positions[i]) sum += 1.0 / (positions[i] - positions[j]);
                return sum;
            }
        }
        for (size_t j = 0; j < n; ++j)
            if (j != i && positions[j] != positions[i]) sum += 1.0 / (positions[i] - positions[j]);
        return sum;
    };
    // With few points to step, their sums are taken in balls at the pass's precision over every point: binary64
    // midpoints cannot tell apart two points closer than its own precision, the close pairs these are.
    const size_t involved = static_cast<size_t>(std::count(conflict.begin(), conflict.end(), uint8_t{1}));
    const bool exact_sums = involved <= heuristics::program_disk_exact_sums;
    bool raise = stalled_raise;
    size_t moved = 0;
    const auto wide_point = [&](size_t i) { return lost_point(points_[i], regions_[i]); };
    ComplexBall delta, denominator, sum, difference, term;
    for (size_t i = 0; i < n; ++i) {
        if (!step[i])
            continue;
        if (limited_[i]) {
            // A value at its rounding: a converged point. Only a disk that must shrink (a culprit's, or one too
            // wide for the goal) asks for more precision; a smaller partner of an overlapping pair waits.
            if (overlaps.culprit[i] || !meets_goal(regions_[i], options, precision_))
                raise = true;
            continue;
        }
        acb_set(delta.value(), corrections_[i].value());
        bool have_sum = false;
        if (conflict[i] && exact_sums) {
            acb_zero(sum.value());
            for (size_t j = 0; j < n; ++j) {
                if (j == i)
                    continue;
                acb_sub(difference.value(), points_[i].value(), points_[j].value(), precision_);
                if (acb_contains_zero(difference.value()))
                    continue;
                acb_inv(term.value(), difference.value(), precision_);
                acb_add(sum.value(), sum.value(), term.value(), precision_);
            }
            have_sum = acb_is_finite(sum.value());
        } else if (conflict[i] && representable) {
            const auto s_i = interaction(i);
            if (std::isfinite(s_i.real()) && std::isfinite(s_i.imag())) {
                acb_set_d_d(sum.value(), s_i.real(), s_i.imag());
                have_sum = true;
            }
        }
        if (have_sum) {
            acb_mul(denominator.value(), corrections_[i].value(), sum.value(), precision_);
            acb_neg(denominator.value(), denominator.value());
            acb_add_ui(denominator.value(), denominator.value(), 1, precision_);
            if (!acb_contains_zero(denominator.value()))
                acb_div(delta.value(), corrections_[i].value(), denominator.value(), precision_);
        }
        // A culprit far from converged (its disk wider than a 2^-20 of its point): a lost point the preliminary
        // phase handed over. With the other points converged and fixed in its sum, Aberth's step is Newton's on
        // the polynomial with their roots divided out, which converges to a root no point holds; so it iterates
        // here, evaluating only itself, until its step is below the precision's resolution or the inner steps run
        // out. Its sum is in balls among few culprits, in binary64 otherwise (a lost point is far from the others,
        // where binary64 resolves them).
        // A narrow culprit whose overlap made no progress last round sits on a root another point holds (two points
        // at one root: the binary64 stage cannot keep two points apart closer than its own precision). It is
        // re-homed the same way, from a start nudged off that root by a 2^-8 of its modulus: with the occupant in its
        // sum the iteration heads for a root no point holds, and if its partner was a distinct close root instead,
        // the only free root is its own, so it comes back and the next raise of precision separates the two.
        // History: docs/slp-second-version-note-20261017.md, section 10.
        const bool lost = overlaps.culprit[i] && wide_point(i);
        const bool doubled = overlaps.culprit[i] && !lost && stalled_rounds_ >= 1;
        if (lost || doubled) {
            ComplexBall value, derivative, point;
            if (lost) {
                acb_sub(point.value(), points_[i].value(), delta.value(), precision_);
            } else {
                // The nudge: 2^-8 of the point's modulus, at an angle from the golden ratio and the point's index.
                const double angle = 2 * 3.141592653589793 * std::fmod(0.6180339887498949 * static_cast<double>(i + 1), 1.0);
                Magnitude size;
                acb_get_mag(size.value(), points_[i].value());
                const double modulus = std::max(mag_get_d(size.value()), 0x1p-1000);
                ComplexBall nudge;
                acb_set_d_d(nudge.value(), std::ldexp(modulus * std::cos(angle), -8), std::ldexp(modulus * std::sin(angle), -8));
                acb_add(point.value(), points_[i].value(), nudge.value(), precision_);
                ++rehomed_;
            }
            acb_get_mid(point.value(), point.value());
            for (unsigned inner = 0; inner < heuristics::program_disk_inner_steps; ++inner) {
                if (options.cancelled())
                    break;
                ++evaluations_;
                ++statistics.evaluations;
                statistics.evaluation_terms += source.terms();
                if (!source.evaluate(value.value(), derivative.value(), point.value(), precision_) ||
                    acb_contains_zero(derivative.value()))
                    break;
                acb_div(delta.value(), value.value(), derivative.value(), precision_);
                if (exact_sums) {
                    acb_zero(sum.value());
                    for (size_t j = 0; j < n; ++j) {
                        if (j == i)
                            continue;
                        acb_sub(difference.value(), point.value(), points_[j].value(), precision_);
                        if (acb_contains_zero(difference.value()))
                            continue;
                        acb_inv(term.value(), difference.value(), precision_);
                        acb_add(sum.value(), sum.value(), term.value(), precision_);
                    }
                } else {
                    // The far points in binary64; the near ones (within 2^-30 of the point, relative) in balls,
                    // exactly: a lost point that reaches a root another point holds must still be pushed off it, and
                    // binary64 would drop that term the moment the two positions coincide.
                    const std::complex<double> z(arf_get_d(arb_midref(acb_realref(point.value())), ARF_RND_NEAR),
                                                 arf_get_d(arb_midref(acb_imagref(point.value())), ARF_RND_NEAR));
                    const double near_distance = 0x1p-30 * std::max(std::abs(z), 0x1p-1000);
                    std::complex<double> total = 0;
                    acb_zero(sum.value());
                    for (size_t j = 0; j < n; ++j) {
                        if (j == i)
                            continue;
                        if (std::abs(z - positions[j]) > near_distance) {
                            total += 1.0 / (z - positions[j]);
                            continue;
                        }
                        acb_sub(difference.value(), point.value(), points_[j].value(), precision_);
                        if (acb_contains_zero(difference.value()))
                            continue;
                        acb_inv(term.value(), difference.value(), precision_);
                        acb_add(sum.value(), sum.value(), term.value(), precision_);
                    }
                    if (!std::isfinite(total.real()) || !std::isfinite(total.imag()))
                        break;
                    acb_set_d_d(term.value(), total.real(), total.imag());
                    acb_add(sum.value(), sum.value(), term.value(), precision_);
                }
                acb_mul(denominator.value(), delta.value(), sum.value(), precision_);
                acb_neg(denominator.value(), denominator.value());
                acb_add_ui(denominator.value(), denominator.value(), 1, precision_);
                if (!acb_contains_zero(denominator.value()))
                    acb_div(delta.value(), delta.value(), denominator.value(), precision_);
                acb_sub(point.value(), point.value(), delta.value(), precision_);
                acb_get_mid(point.value(), point.value());
                ++steps_;
                // Small enough: the next round's disk says the rest.
                Magnitude step_size, size;
                acb_get_mag(step_size.value(), delta.value());
                acb_get_mag(size.value(), point.value());
                mag_mul_2exp_si(size.value(), size.value(), -(precision_ / 2));
                if (mag_cmp(step_size.value(), size.value()) <= 0)
                    break;
            }
            acb_swap(points_[i].value(), point.value());
            // The binary64 position the later points of this round read, moved with the point.
            positions[i] = {arf_get_d(arb_midref(acb_realref(points_[i].value())), ARF_RND_NEAR),
                            arf_get_d(arb_midref(acb_imagref(points_[i].value())), ARF_RND_NEAR)};
            current_[i] = 0;
            ++moved;
            continue;
        }
        acb_sub(points_[i].value(), points_[i].value(), delta.value(), precision_);
        acb_get_mid(points_[i].value(), points_[i].value());
        current_[i] = 0;
        ++moved;
        ++steps_;
    }
    // A raise only after progress since the previous one: more precision separates close roots and narrows the
    // disks of converged points, but two points at one root stay together at any precision.
    if (raise && raised_since_progress_ && overlapping_)
        raise = false;
    if (raise) {
        if (precision_ >= cap)
            return moved ? state_ : fail("values at their rounding at the precision cap");
        precision_ = std::min(cap, resources::precision_product(precision_, 2));
        ++raises_;
        raised_since_progress_ = true;
        for (size_t i = 0; i < n; ++i)
            if (step[i])
                current_[i] = 0;
    }
    return state_;
}

SeparatedRootsPass::ClusterOutcome SeparatedRootsPass::cluster_cover(MonicSource &source, MonicSource *undivided,
                                                                     const cps_options &options,
                                                                     certpolysolve::detail::CandidateStatistics &statistics,
                                                                     const std::vector<std::pair<size_t, size_t>> &pairs,
                                                                     const std::vector<uint8_t> &involved) {
    const size_t n = points_.size();
    const auto cap = static_cast<slong>(options.max_work_bits);
    // Counted clusters answer the approximation goal over the complex roots and the cluster-radius goal. The
    // accuracy and isolation goals, and real-only queries, need every root alone, which no cluster region gives.
    const bool admitted =
        options.goal == CPS_CLUSTER_RADIUS || (options.goal == CPS_ROOT_APPROXIMATION && !options.real_only);
    if (!admitted || pairs.empty())
        return ClusterOutcome::None;
    // Before any evaluation: a lost point belongs to no cluster, and a separated disk too wide for the goal leaves
    // the cover short of it.
    size_t separated = 0;
    for (size_t i = 0; i < n; ++i) {
        if (involved[i]) {
            if (lost_point(points_[i], regions_[i]))
                return ClusterOutcome::None;
        } else {
            if (!meets_goal(regions_[i], options, precision_))
                return ClusterOutcome::None;
            ++separated;
        }
    }
    // The components of the overlapping pairs (union-find, each component named by its smallest index), with the
    // involved points no pair holds (disks marked outright by the sweep) as components of their own.
    std::vector<size_t> parent(n);
    std::iota(parent.begin(), parent.end(), size_t{0});
    const auto find = [&](size_t i) {
        while (parent[i] != i)
            i = parent[i] = parent[parent[i]];
        return i;
    };
    for (const auto &[i, j] : pairs) {
        const size_t a = find(i), b = find(j);
        if (a != b)
            parent[std::max(a, b)] = std::min(a, b);
    }
    std::vector<std::vector<size_t>> members;
    std::vector<size_t> component(n, SIZE_MAX);
    for (size_t i = 0; i < n; ++i) {
        if (!involved[i])
            continue;
        const size_t root = find(i);
        if (component[root] == SIZE_MAX) {
            component[root] = members.size();
            members.emplace_back();
        }
        members[component[root]].push_back(i);
    }
    // A count on a circle, charged to the pass; `expected` roots size its arcs. With an undivided polynomial z^v q
    // (advance()), a circle whose centre's midpoint (where certified_root_count puts the circle) lies at least
    // program_undivided_count_margin v + 1 radii from zero, by a lower bound of its modulus, is counted on it: its
    // disk keeps zero outside, so it holds q's roots there and its count is q's. Any other circle is counted on q.
    // One count only: far from zero q's enclosures are p's divided, so a count that fails on p fails on q too, and a
    // second count would only double the cost.
    // History: docs/slp-second-version-note-20261017.md, section 10.
    const auto count = [&](const ComplexBall &centre, const Magnitude &radius, uint64_t expected) {
        MonicSource *polynomial = &source;
        if (undivided) {
            const auto order = static_cast<ulong>(undivided->degree() - source.degree());
            ComplexBall middle;
            Magnitude distance, far;
            acb_get_mid(middle.value(), centre.value());
            acb_get_mag_lower(distance.value(), middle.value());
            mag_mul_ui(far.value(), radius.value(), heuristics::program_undivided_count_margin * order + 1);
            if (mag_cmp(distance.value(), far.value()) > 0)
                polynomial = undivided;
        }
        const auto counted = certified_root_count(
            *polynomial, centre.value(), radius.value(), precision_,
            heuristics::program_winding_arc_base + heuristics::program_winding_arc_scale * (expected + 2), options,
            statistics);
        ++counts_;
        count_arcs_ += counted.arcs;
        evaluations_ += counted.evaluations;
        return counted;
    };
    // Each component's first circle: about its points' centroid, program_cluster_circle times the reach of its
    // published squares (the farthest of their points from the centroid), so the circle holds every member disk
    // and so at least one root.
    struct Cluster {
        ComplexBall centre;
        Magnitude radius;
        uint64_t count = 0;
        bool meets = false;
    };
    std::vector<Cluster> clusters(members.size());
    uint64_t counted_roots = 0;
    bool raised = false; // the one raise allowed for the first counts
    ComplexBall sum, difference;
    Magnitude reach, distance;
    for (size_t k = 0; k < members.size(); ++k) {
        auto &cluster = clusters[k];
        acb_zero(sum.value());
        for (const size_t i : members[k])
            acb_add(sum.value(), sum.value(), points_[i].value(), precision_);
        acb_div_ui(sum.value(), sum.value(), static_cast<ulong>(members[k].size()), precision_);
        acb_get_mid(cluster.centre.value(), sum.value());
        mag_zero(reach.value());
        for (const size_t i : members[k]) {
            acb_sub(difference.value(), regions_[i].box.value(), cluster.centre.value(), precision_);
            acb_get_mag(distance.value(), difference.value());
            mag_max(reach.value(), reach.value(), distance.value());
        }
        mag_mul_ui(cluster.radius.value(), reach.value(), heuristics::program_cluster_circle);
        auto counted = count(cluster.centre, cluster.radius, members[k].size());
        // Values on the circle at their rounding (points converged to a multiple root at the pass's precision):
        // the precision doubles once, to the cap, and the count is taken again.
        if (!counted.roots && counted.imprecise && !counted.cancelled && !raised && precision_ < cap) {
            precision_ = std::min(cap, resources::precision_product(precision_, 2));
            ++raises_;
            raised = true;
            counted = count(cluster.centre, cluster.radius, members[k].size());
        }
        if (counted.cancelled)
            return ClusterOutcome::Cancelled;
        if (!counted.roots || *counted.roots == 0)
            return ClusterOutcome::None;
        cluster.count = *counted.roots;
        counted_roots += cluster.count;
    }
    // Complete only when the counts and the separated disks add up to the degree (two points at one root leave it
    // short; the refinement answers then).
    if (counted_roots + separated != n)
        return ClusterOutcome::None;
    // A cluster's region: the square of half-side `radius` about its centre, which holds the circle's closed disk,
    // rounded outward when published.
    const auto region = [&](const Cluster &cluster, mag_srcptr radius) {
        Root r;
        acb_set(r.box.value(), cluster.centre.value());
        arb_add_error_mag(acb_realref(r.box.value()), radius);
        arb_add_error_mag(acb_imagref(r.box.value()), radius);
        r.count = cluster.count;
        r.multiplicity = cluster.count == 1 ? 1 : 0;
        r.publish(precision_);
        r.realness = arb_contains_zero(acb_imagref(r.box.value())) ? CPS_REAL_UNKNOWN : CPS_NONREAL;
        r.certificate = "original-source argument-principle count in a complete disjoint cover";
        return r;
    };
    // Each circle shrunk toward the goal about its centre: the largest power of two whose square meets the goal. The
    // roots of a multiple root or a tight cluster stay where they are while two points converging on them gain
    // only about 1.6 bits a round (Aberth's linear rate at a double root), so the circle shrinks instead of waiting
    // for the points. A count at precision p resolves a circle of radius about 2^(-p/m) about an m-fold root, so
    // while a goal's circle is unreadable the precision doubles, to the cap.
    Magnitude trial;
    const auto power = [&](slong e) -> mag_srcptr {
        mag_one(trial.value());
        mag_mul_2exp_si(trial.value(), trial.value(), e);
        return trial.value();
    };
    for (;;) {
        bool all = true, imprecise = false;
        for (auto &cluster : clusters) {
            if (cluster.meets)
                continue;
            if (meets_goal(region(cluster, cluster.radius.value()), options, precision_)) {
                cluster.meets = true;
                continue;
            }
            const auto meets_at = [&](slong e) { return meets_goal(region(cluster, power(e)), options, precision_); };
            // The largest e at most `top` (2^top at most the radius) whose square meets the goal: down from `top` by
            // doubling steps, then by bisection. Never met (a cluster about 0 under a relative goal): no cover.
            const auto top = static_cast<slong>(std::floor(mag_get_d_log2_approx(cluster.radius.value())));
            slong hit = top, miss = top + 1;
            if (!meets_at(top)) {
                miss = top;
                bool found = false;
                for (slong step = 1; step <= (slong{1} << 24); step *= 2) {
                    hit = top - step;
                    if (meets_at(hit)) {
                        found = true;
                        break;
                    }
                    miss = hit;
                }
                if (!found)
                    return ClusterOutcome::None;
                while (miss - hit > 1) {
                    const slong middle = hit + (miss - hit) / 2;
                    if (meets_at(middle))
                        hit = middle;
                    else
                        miss = middle;
                }
            }
            Magnitude radius;
            mag_set(radius.value(), power(hit));
            const auto at_goal = count(cluster.centre, radius, cluster.count);
            if (at_goal.cancelled)
                return ClusterOutcome::Cancelled;
            if (at_goal.roots && *at_goal.roots == cluster.count) {
                mag_set(cluster.radius.value(), radius.value());
                cluster.meets = true;
                continue;
            }
            // Values at their rounding on the goal's circle: more precision may read them. A readable count that
            // differs, or a circle that cannot be resolved, says the roots are not one cluster at the goal (distinct
            // roots the pass can still separate, or a centre off the cluster): no cover now.
            if (!at_goal.imprecise)
                return ClusterOutcome::None;
            all = false;
            imprecise = true;
        }
        if (all)
            break;
        if (!imprecise || precision_ >= cap)
            return ClusterOutcome::None;
        precision_ = std::min(cap, resources::precision_product(precision_, 2));
        ++raises_;
    }
    // The cover: the separated disks and the clusters, pairwise disjoint (the certificate's own check, over the
    // published squares).
    std::vector<Root> cover;
    cover.reserve(separated + clusters.size());
    for (size_t i = 0; i < n; ++i)
        if (!involved[i]) {
            cover.push_back(regions_[i]);
            cover.back().certificate = "original-source Newton inclusion disk in a complete disjoint cover";
        }
    for (const auto &cluster : clusters)
        cover.push_back(region(cluster, cluster.radius.value()));
    if (!disjoint(cover))
        return ClusterOutcome::None;
    cover_ = std::move(cover);
    clusters_ = clusters.size();
    return ClusterOutcome::Certified;
}

} // namespace certpolysolve::detail
