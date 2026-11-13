// Certified subdivision and contraction (Subdivision, subdivision.hpp): a persistent quadtree search from the
// certified cover of a factor p to certified regions that meet the goal, one controller action a call of step.
// The code-to-paper progress and complexity obligations are recorded separately in
// dev/docs/m2a-fallback-design-20261105.md; sound certificates alone do not establish those bounds.
//
// The constructor makes each region of the cover a square of power-of-two half-width (or a point) that keeps the
// region as its certificate, merging overlapping squares into uncertified ones; step acts on the front of a
// round-robin queue of components, on p from the caller's Encloser. Its main cases, in order: a split in progress cuts
// each cell into four children and tests one a step, discarding a child only when Pellet's test with Graeffe
// iterations (pellet_count, counting.cpp) proves no root in the disc of 3/2 its half-width, which contains the closed
// square, and the retained children regroup into eight-neighbour connected components (connected); a certified one
// that meets the goal is finished; an uncertified one whose square of twice its radius misses every other component
// (separated: its certificate then overlaps none) is counted and certified by two discs (certify_cluster_box); a
// certified one tries the Newton step c - k p(c)/p'(c), k its certified count, to a square 2^-newton_bits as wide as
// its certificate's, kept only when certified with the same count (newton), newton_bits doubling on success and
// halving on failure; otherwise it starts a split, keeping its certificate in before_split for abandon_split to
// restore if the precision cap stops the split. Precision doubles per component (raise). Steps are metered (work, at
// operation_weight); a disk test stops before an operation that would pass the allowance, and the next step resumes it.
//
// State: Subdivision::Impl: active (the queue of Component), finished, total_count, multiplicity, the meter (work,
// suspended, idle, pending_work) and counting (CountContext: the allowance, a suspended disk test). A Component: cells
// (Cell: centre, half-width, grid indices), the split (untested, retained), precision, newton_bits, certificate.
// Entry points: Subdivision (step; complete_cover, the regions once all are certified and disjoint; reopen; counters),
// called from solver.cpp only (SolveState::Impl: candidate_handoff and factor_step build it from a factor's certified
// cover, as FactorState::search at a hand-over or a local refinement, or as FactorState::dovetail_search beside
// another method, as on the dovetail share) and from tests/test_sessions.cpp and tests/test_method_handoff.cpp.
// Not here: the disk tests (counting.cpp), the dovetail schedule (dovetail.hpp), meets_goal (certification.cpp).
// Design notes: docs/mathematics.md, "Persistent certified subdivision and contraction";
// docs/principled-audit-20260929.md, "Dovetailing design" and "Dovetailing as built and qualified".
#include "methods/subdivision/subdivision.hpp"
#include "core/telemetry.hpp"
#include "core/dovetail.hpp"
#include <algorithm>
#include <deque>
#include <map>
#include <numeric>
#include <optional>
#include <stdexcept>

namespace certpolysolve::detail {
namespace {
struct Cell {
    ComplexBall center;
    RealBall half; // exact positive power of two
    Integer x, y;  // centers on a common grid, adjacent indices differ by two
};
// A Newton step contracts a certificate of radius R to R 2^-newton_bits. A
// cell of half-width h that was split regroups with a certificate of radius
// 2h (the power-of-two bound of its diagonal), so two bits is the least
// contraction that is never coarser than the cells of a fresh bisection:
// one bit would return the cell that was split, and a search that took it
// after every bisection would cycle. It is the exponent a component starts
// with, and the least one a failed step is retried at (step's Newton
// branch). A one-bit step is still taken by a component that persists
// unsplit after a failure at two bits (its precision raised, or reopened at
// a larger budget): no bisection stands between it and its certificate.
inline constexpr uint64_t least_newton_bits = 2;
struct Component {
    std::vector<Cell> cells, retained;
    std::deque<Cell> untested;
    size_t split_index = 0;
    bool splitting = false, blocked = false;
    cps_component_state state = CPS_COMPONENT_ACTIVE;
    uint64_t precision = 64, newton_bits = least_newton_bits;
    std::optional<uint64_t> trial_precision;
    uint64_t working_precision() const { return trial_precision.value_or(precision); }
    ComplexBall bounds, center;
    RealBall radius;
    ComplexBall certificate_center;
    RealBall certificate_radius;
    std::optional<Root> certificate;
    // The certificate held before the current split. It is restored if the
    // split cannot finish within the precision budget.
    std::optional<Root> before_split;
};

void power_bound(arb_t out, const arf_t radius) {
    Integer exponent;
    arf_abs_bound_le_2exp_fmpz(exponent.value(), radius);
    arb_one(out);
    arb_mul_2exp_fmpz(out, out, exponent.value());
}
bool snap_center(acb_t center, const arb_t radius, slong precision) {
    if (arb_is_zero(radius))
        return true;
    Float grid, scaled;
    arf_mul_2exp_si(grid.value(), arb_midref(radius), -3);
    for (auto coordinate : {acb_realref(center), acb_imagref(center)}) {
        arf_div(scaled.value(), arb_midref(coordinate), grid.value(), resources::precision_sum(precision, 16), ARF_RND_NEAR);
        // Do not materialize an integer spanning an unbounded exponent gap.
        if (arf_cmpabs_2exp_si(scaled.value(), resources::precision_sum(precision, 16)) >= 0)
            return false;
        Integer index;
        arf_get_fmpz(index.value(), scaled.value(), ARF_RND_NEAR);
        arf_set_fmpz(arb_midref(coordinate), index.value());
        arf_mul(arb_midref(coordinate), arb_midref(coordinate), grid.value(), ARF_PREC_EXACT, ARF_RND_NEAR);
    }
    return true;
}
bool accurate_newton_proposal(const acb_t proposal, const arb_t radius) {
    // r/64 is the local-progress accuracy allowance. Together with the r/8
    // coordinate snap and exact geometric error, it leaves the roots in r/4.
    constexpr slong accuracy_bits = 6;
    if (!acb_is_finite(proposal))
        return false;
    Magnitude error;
    mag_add(error.value(), arb_radref(acb_realref(proposal)), arb_radref(acb_imagref(proposal)));
    Float upper, allowance;
    arf_set_mag(upper.value(), error.value());
    arf_mul_2exp_si(allowance.value(), arb_midref(radius), -accuracy_bits);
    return arf_cmp(upper.value(), allowance.value()) <= 0;
}
void geometry(Component &c, slong precision) {
    if (c.cells.empty())
        throw std::logic_error("empty subdivision component");
    // Every cell in a component has the same exact half-width. Find the
    // exact coordinate extrema of their centers, then round each hull
    // radius once. Repeated unions of rounded cell boxes can accumulate
    // radius inflation with the number of cells.
    ComplexBall lower(c.cells.front().center), upper(lower);
    for (size_t i = 1; i < c.cells.size(); ++i) {
        for (bool imaginary : {false, true}) {
            const auto coordinate = imaginary ? acb_imagref(c.cells[i].center.value())
                                              : acb_realref(c.cells[i].center.value());
            auto minimum = imaginary ? acb_imagref(lower.value()) : acb_realref(lower.value());
            auto maximum = imaginary ? acb_imagref(upper.value()) : acb_realref(upper.value());
            if (arf_cmp(arb_midref(coordinate), arb_midref(minimum)) < 0)
                arf_set(arb_midref(minimum), arb_midref(coordinate));
            if (arf_cmp(arb_midref(coordinate), arb_midref(maximum)) > 0)
                arf_set(arb_midref(maximum), arb_midref(coordinate));
        }
    }
    for (bool imaginary : {false, true}) {
        const auto minimum = imaginary ? acb_imagref(lower.value()) : acb_realref(lower.value());
        const auto maximum = imaginary ? acb_imagref(upper.value()) : acb_realref(upper.value());
        auto bound = imaginary ? acb_imagref(c.bounds.value()) : acb_realref(c.bounds.value());
        arf_add(arb_midref(bound), arb_midref(minimum), arb_midref(maximum), ARF_PREC_EXACT, ARF_RND_NEAR);
        arf_mul_2exp_si(arb_midref(bound), arb_midref(bound), -1);
        // Form the center span before adding the half-width. A single
        // cell can have a huge gap between its center and its radius;
        // materializing its exact endpoints would allocate across that gap.
        Float half;
        arf_sub(half.value(), arb_midref(maximum), arb_midref(minimum), ARF_PREC_EXACT, ARF_RND_NEAR);
        arf_mul_2exp_si(half.value(), half.value(), -1);
        arf_add(half.value(), half.value(), arb_midref(c.cells.front().half.value()), ARF_PREC_EXACT, ARF_RND_NEAR);
        arf_get_mag(arb_radref(bound), half.value());
    }
    acb_get_mid(c.center.value(), c.bounds.value());
    ComplexBall offset;
    acb_sub(offset.value(), c.bounds.value(), c.center.value(), resources::precision_sum(precision, 16));
    Float radius;
    acb_get_abs_ubound_arf(radius.value(), offset.value(), resources::precision_sum(precision, 16));
    power_bound(c.radius.value(), radius.value());
}
std::pair<std::string, std::string> key(const Integer &x, const Integer &y) {
    return {integer_string(x.value()), integer_string(y.value())};
}
std::vector<Component> connected(std::vector<Cell> cells, const Component &parent) {
    // All children share a grid. Eight-neighbor union includes corners, which
    // is essential for CLOSED squares and roots on subdivision lines.
    std::vector<size_t> representatives(cells.size());
    std::iota(representatives.begin(), representatives.end(), 0);
    auto find = [&](size_t i) {
        while (representatives[i] != i) {
            representatives[i] = representatives[representatives[i]];
            i = representatives[i];
        }
        return i;
    };
    std::map<std::pair<std::string, std::string>, size_t> index;
    for (size_t i = 0; i < cells.size(); ++i) {
        for (int dx : {-2, 0, 2})
            for (int dy : {-2, 0, 2}) {
                Integer x(cells[i].x), y(cells[i].y);
                fmpz_add_si(x.value(), x.value(), dx);
                fmpz_add_si(y.value(), y.value(), dy);
                auto it = index.find(key(x, y));
                if (it != index.end())
                    representatives[find(i)] = find(it->second);
            }
        index.emplace(key(cells[i].x, cells[i].y), i);
    }
    std::map<size_t, size_t> groups;
    std::vector<Component> result;
    for (size_t i = 0; i < cells.size(); ++i) {
        const auto root = find(i);
        auto [it, inserted] = groups.emplace(root, result.size());
        if (inserted) {
            Component c;
            c.precision = parent.precision;
            result.push_back(std::move(c));
        }
        result[it->second].cells.push_back(std::move(cells[i]));
    }
    for (auto &c : result)
        geometry(c, static_cast<slong>(c.precision));
    return result;
}
} // namespace
struct Subdivision::Impl {
    std::deque<Component> active;
    std::vector<Component> finished;
    // Binary64-equivalent operations spent in step() (dovetail.hpp).
    double work = 0;
    uint64_t multiplicity, total_count;
    certpolysolve::detail::SearchStatistics &statistics;
    CountContext counting;
    const SubdivisionPrecision precision_policy;
    bool complete_entry = false;
    // The last step's disk test stopped for its allowance (step), and the
    // weighted cost of the operation it waits for.
    bool suspended = false, idle = false;
    double pending_work = 0;
    Impl(const std::vector<Root>& cover, uint64_t m, uint64_t precision, certpolysolve::detail::SearchStatistics& s,
         SubdivisionPrecision policy)
        : multiplicity(m), total_count(0), statistics(s), counting(s), precision_policy(policy) {
        // Subdivision retries a test that arithmetic width obstructs
        // before it asks for precision: pellet_count retries on the rotated
        // real shift (class G). At a precision budget this is the difference
        // between a count and a precision stop.
        counting.real_shift = true;
        for (const auto &global : cover) {
            total_count += global.count;
            Cell cell;
            acb_get_mid(cell.center.value(), global.box.value());
            Float radius, extra;
            arf_set_mag(radius.value(), arb_radref(acb_realref(global.box.value())));
            arf_set_mag(extra.value(), arb_radref(acb_imagref(global.box.value())));
            if (arf_cmp(extra.value(), radius.value()) > 0)
                arf_set(radius.value(), extra.value());
            if (arf_is_zero(radius.value()))
                arb_zero(cell.half.value());
            else {
                power_bound(cell.half.value(), radius.value());
                ComplexBall snapped(cell.center);
                if (snap_center(snapped.value(), cell.half.value(), static_cast<slong>(precision))) {
                    cell.center = std::move(snapped);
                    // Contains the original rectangle despite rounding its center.
                    arb_mul_2exp_si(cell.half.value(), cell.half.value(), 1);
                }
            }
            Component c;
            c.precision = precision;
            c.cells.push_back(std::move(cell));
            geometry(c, static_cast<slong>(c.precision));
            c.certificate = global;
            c.certificate_center = c.cells.front().center;
            c.certificate_radius = c.cells.front().half;
            active.push_back(std::move(c));
        }
        // Initial certified rectangles may become overlapping squares. Merge
        // their supports before assigning root sets to independent searches.
        bool changed = true;
        while (changed) {
            changed = false;
            for (size_t i = 0; i < active.size(); ++i)
                for (size_t j = i + 1; j < active.size();) {
                    if (!acb_overlaps(active[i].bounds.value(), active[j].bounds.value())) {
                        ++j;
                        continue;
                    }
                    auto &c = active[i];
                    ComplexBall united;
                    acb_union(united.value(), c.bounds.value(), active[j].bounds.value(),
                              static_cast<slong>(resources::precision_sum(precision, 16)));
                    Cell cell;
                    acb_get_mid(cell.center.value(), united.value());
                    Float r, im;
                    arf_set_mag(r.value(), arb_radref(acb_realref(united.value())));
                    arf_set_mag(im.value(), arb_radref(acb_imagref(united.value())));
                    if (arf_cmp(im.value(), r.value()) > 0)
                        arf_set(r.value(), im.value());
                    power_bound(cell.half.value(), r.value());
                    c.cells.clear();
                    c.cells.push_back(std::move(cell));
                    c.certificate.reset();
                    geometry(c, static_cast<slong>(c.precision));
                    active.erase(active.begin() + static_cast<std::ptrdiff_t>(j));
                    changed = true;
                }
        }
    }
    // The outer square must miss every other component's cell support AND its
    // certificate box. A certificate issued inside this outer square is then
    // disjoint from every certificate issued before it, so the final cover can
    // always be assembled (complete_cover's disjointness never fails once all
    // components are certified).
    bool separated(const Component &c) const {
        ComplexBall outer(c.center);
        Float radius;
        arf_mul_2exp_si(radius.value(), arb_midref(c.radius.value()), 1);
        acb_add_error_arf(outer.value(), radius.value());
        const auto clear = [&](const Component &other) {
            const auto misses = [&](const std::optional<Root> &root) {
                return !root || !acb_overlaps(outer.value(), root->box.value());
            };
            if (acb_overlaps(outer.value(), other.bounds.value())) {
                // During a partial split, cells still holds the old support.
                // Empty support retains the conservative bounds check.
                if (other.cells.empty())
                    return false;
                for (const auto &cell : other.cells) {
                    ComplexBall box(cell.center);
                    acb_add_error_arf(box.value(), arb_midref(cell.half.value()));
                    if (acb_overlaps(outer.value(), box.value()))
                        return false;
                }
            }
            return misses(other.certificate) && misses(other.before_split);
        };
        for (const auto &other : active)
            if (&other != &c && !clear(other))
                return false;
        for (const auto &other : finished)
            if (!clear(other))
                return false;
        return true;
    }
    // A split that the precision budget cannot finish leaves the component at
    // its resolution before the split, with its certificate, stopped for the
    // precision budget. separated() kept that certificate's box clear.
    void abandon_split(Component &c) {
        c.untested.clear();
        c.retained.clear();
        c.split_index = 0;
        c.splitting = false;
        c.certificate = std::move(c.before_split);
        c.before_split.reset();
        c.blocked = true;
        c.state = CPS_COMPONENT_PRECISION_LIMIT;
    }
    bool raise(Component &c, const cps_options &o) {
        if (c.precision == o.max_work_bits)
            return false;
        c.precision = std::min(o.max_work_bits, c.precision * 2);
        ++statistics.precision_increases;
        return true;
    }
    bool raise_trial(Component &c, const cps_options &o) {
        auto &precision = *c.trial_precision;
        if (precision >= o.max_work_bits)
            return false;
        precision = std::min(o.max_work_bits, precision * 2);
        ++statistics.precision_increases;
        return true;
    }
    bool newton(Component &c, const acb_poly_t p, const cps_options &o, bool &needs_precision,
                NewtonAccuracy accuracy, slong prec) {
        ExclusiveTimer timer(statistics.refinement_ns, &statistics.counting_ns, &statistics.taylor_shift_ns);
        ++statistics.newton_attempts;
        const slong count = static_cast<slong>(c.certificate->count / multiplicity);
        ComplexBall old_box(c.certificate_center);
        acb_add_error_arf(old_box.value(), arb_midref(c.certificate_radius.value()));
        ComplexBall point(c.certificate_center), value, derivative, step, next;
        const auto distant_proposal = [&] {
            // k|F(x)| > (|x-c|+2R)|F'(x)| rules out the old square, whose
            // diagonal radius is sqrt(2)R. Use the actual rounded sample x.
            ComplexBall offset;
            RealBall reach, value_size, derivative_size;
            acb_sub(offset.value(), point.value(), c.certificate_center.value(), prec);
            acb_abs(reach.value(), offset.value(), prec);
            arb_addmul_ui(reach.value(), c.certificate_radius.value(), 2, prec);
            acb_abs(value_size.value(), value.value(), prec);
            arb_mul_si(value_size.value(), value_size.value(), count, prec);
            acb_abs(derivative_size.value(), derivative.value(), prec);
            arb_mul(derivative_size.value(), derivative_size.value(), reach.value(), prec);
            return arb_gt(value_size.value(), derivative_size.value()) != 0;
        };
        // Keep a useful center proposal. If its derivative is unresolved or
        // a multi-root proposal is outside (or a proposal is insufficiently accurate under
        // the exact-factor policy), try the annular point c + 3R/2.
        // A native equal-count (R, 2R) certificate makes this annulus root-free;
        // inherited covers still require the same strict candidate certificate.
        // Legacy singleton fallback keeps c + R. The prepared-factor refiner
        // uses the same empty annulus and proposal guard as isolation. Only
        // one candidate is counted per call.
        for (unsigned sample = 0; sample < 2; ++sample) {
            acb_set(point.value(), c.certificate_center.value());
            const bool annular = sample != 0 && (count > 1 || accuracy == NewtonAccuracy::PreparedExact);
            if (annular) {
                // Mid-annulus 3R/2, with R/4 clearance from either boundary.
                // Check the actual rounded displacement, not the intended one.
                RealBall offset(c.certificate_radius), actual, lower(c.certificate_radius), upper(c.certificate_radius);
                arb_mul_ui(offset.value(), offset.value(), 3, prec);
                arb_mul_2exp_si(offset.value(), offset.value(), -1);
                arb_add(acb_realref(point.value()), acb_realref(point.value()), offset.value(), prec);
                acb_get_mid(point.value(), point.value());
                arb_sub(actual.value(), acb_realref(point.value()), acb_realref(c.certificate_center.value()), prec);
                arb_mul_ui(lower.value(), lower.value(), 5, prec);
                arb_mul_2exp_si(lower.value(), lower.value(), -2);
                arb_mul_ui(upper.value(), upper.value(), 7, prec);
                arb_mul_2exp_si(upper.value(), upper.value(), -2);
                if (!acb_is_finite(point.value()) || !arb_ge(actual.value(), lower.value()) ||
                    !arb_le(actual.value(), upper.value())) {
                    needs_precision = true;
                    return false;
                }
            } else if (sample != 0 && count == 1) {
                arb_add(acb_realref(point.value()), acb_realref(point.value()), c.certificate_radius.value(), prec);
                acb_get_mid(point.value(), point.value());
            }
            acb_poly_evaluate2(value.value(), derivative.value(), p, point.value(), prec);
            if (acb_contains_zero(derivative.value())) {
                if (sample == 0)
                    continue;
                needs_precision = !acb_is_exact(derivative.value());
                // F is nonzero at an exact stationary point of a square-free
                // factor, so this becomes decisive without resolving F'.
                if (needs_precision && distant_proposal())
                    needs_precision = false;
                return false;
            }
            acb_div(step.value(), value.value(), derivative.value(), prec);
            acb_mul_si(step.value(), step.value(), count, prec);
            acb_sub(next.value(), point.value(), step.value(), prec);
            if (accuracy == NewtonAccuracy::PreparedExact) {
                RealBall trial_radius(c.certificate_radius);
                arb_mul_2exp_si(trial_radius.value(), trial_radius.value(), -static_cast<slong>(c.newton_bits));
                if (!accurate_newton_proposal(next.value(), trial_radius.value())) {
                    if (sample == 0)
                        continue;
                    needs_precision = !distant_proposal();
                    return false;
                }
            }
            acb_get_mid(next.value(), next.value());
            if (!acb_is_finite(next.value()))
                return false;
            // Snapping moves a coordinate by at most the new half-width / 8.
            // A contained new square needs at least that half-width of margin,
            // so an outside midpoint cannot become a usable proposal. Reject it
            // before snapping asks for precision to materialize a distant index.
            if (!acb_contains(old_box.value(), next.value())) {
                if (sample == 0 && (count > 1 || accuracy == NewtonAccuracy::PreparedExact))
                    continue;
                return false;
            }
            break;
        }
        RealBall radius(c.certificate_radius);
        arb_mul_2exp_si(radius.value(), radius.value(), -static_cast<slong>(c.newton_bits));
        if (!snap_center(next.value(), radius.value(), prec)) {
            needs_precision = true;
            return false;
        }
        ComplexBall new_box(next);
        acb_add_error_arf(new_box.value(), arb_midref(radius.value()));
        if (!acb_contains(old_box.value(), new_box.value()) || o.cancelled())
            return false;
        Root root;
        if (!certify_cluster_box(root, p, next.value(), arb_midref(radius.value()), count, multiplicity, prec,
                                 subdivision_graeffe_steps(acb_poly_degree(p)), &counting, &needs_precision)) {
            return false;
        }
        // Equal strict counts in nested regions preserve exactly the parent's
        // root set. Numerical Newton values alone never justify discarding cells.
        Cell cell;
        cell.center = std::move(next);
        cell.half = std::move(radius);
        c.cells.clear();
        c.cells.push_back(std::move(cell));
        geometry(c, prec);
        c.certificate = std::move(root);
        c.certificate_center = c.cells.front().center;
        c.certificate_radius = c.cells.front().half;
        c.newton_bits = std::min<uint64_t>(resources::precision_max, c.newton_bits * 2);
        ++statistics.newton_successes;
        return true;
    }
};
Subdivision::Subdivision(const std::vector<Root>& r, uint64_t m, uint64_t p, certpolysolve::detail::SearchStatistics& s,
                         SubdivisionPrecision policy)
    : impl_(std::make_unique<Impl>(r, m, p, s, policy)) {}
Subdivision::~Subdivision() = default;
namespace {
struct CompleteSquare {
    Root root;
    ComplexBall center;
    RealBall half;
    bool padded = false;
};

bool complete_square(CompleteSquare &out, const acb_t hull, const arf_t floor,
                     uint64_t count, uint64_t multiplicity, slong precision) {
    ComplexBall rounded;
    acb_set_round(rounded.value(), hull, precision);
    acb_get_mid(out.center.value(), rounded.value());
    Float radius, imaginary;
    arf_set_mag(radius.value(), arb_radref(acb_realref(rounded.value())));
    arf_set_mag(imaginary.value(), arb_radref(acb_imagref(rounded.value())));
    if (arf_cmp(imaginary.value(), radius.value()) > 0)
        arf_set(radius.value(), imaginary.value());
    // Four times the enclosing rectangle's half-width leaves its corners
    // strictly inside the R/2 disk, even after a center snap of at most R/16
    // in each coordinate. Thus both the R and 2R disks contain all roots.
    arf_mul_2exp_si(radius.value(), radius.value(), 2);
    out.padded = arf_cmp(radius.value(), floor) < 0;
    if (out.padded)
        arf_set(radius.value(), floor);
    power_bound(out.half.value(), radius.value());
    if (!snap_center(out.center.value(), out.half.value(), precision))
        return false;
    ComplexBall offset;
    Float reach, inner;
    acb_sub(offset.value(), rounded.value(), out.center.value(), precision);
    acb_get_abs_ubound_arf(reach.value(), offset.value(), precision);
    arf_mul_2exp_si(inner.value(), arb_midref(out.half.value()), -1);
    if (!arf_is_finite(reach.value()) || arf_cmp(reach.value(), inner.value()) >= 0)
        return false;
    acb_set(out.root.box.value(), out.center.value());
    for (auto coordinate : {acb_realref(out.root.box.value()), acb_imagref(out.root.box.value())})
        arf_get_mag(arb_radref(coordinate), arb_midref(out.half.value()));
    out.root.count = count;
    out.root.multiplicity = count == multiplicity ? multiplicity : 0;
    out.root.certificate = "complete-entry-all-root-count";
    out.root.publish(precision);
    // The published certificate also stays in the outer disk. This check
    // includes any outward endpoint rounding, as native local counts do.
    acb_sub(offset.value(), out.root.box.value(), out.center.value(), precision);
    acb_get_abs_ubound_arf(reach.value(), offset.value(), precision);
    arf_mul_2exp_si(inner.value(), arb_midref(out.half.value()), 1);
    return arf_is_finite(reach.value()) && arf_cmp(reach.value(), inner.value()) < 0;
}

bool complete_warm_cover(const Root &global, const std::vector<Root> &cover, uint64_t multiplicity) {
    if (cover.empty() || cover.size() > global.count / multiplicity)
        return false;
    uint64_t count = 0;
    for (const auto &root : cover) {
        if (!root.count || root.count % multiplicity || root.count > global.count - count ||
            (root.multiplicity && root.multiplicity != multiplicity) ||
            !acb_is_finite(root.box.value()) || !acb_contains(global.box.value(), root.box.value()))
            return false;
        count += root.count;
    }
    return count == global.count && disjoint(cover);
}
} // namespace

std::unique_ptr<Subdivision> Subdivision::complete_start(
    const Root &global, const std::vector<Root> &warm, uint64_t multiplicity,
    uint64_t seed_scale_bits, SearchStatistics &statistics, SubdivisionCompleteEntry *report) {
    if (report) *report = {};
    // The complete factor's Cauchy square has exact center zero. Do not accept
    // an arbitrary count-n warm rectangle in place of this independent bound.
    if (!multiplicity || !global.count || global.count % multiplicity ||
        !acb_is_finite(global.box.value()) ||
        !arf_is_zero(arb_midref(acb_realref(global.box.value()))) ||
        !arf_is_zero(arb_midref(acb_imagref(global.box.value()))) ||
        !mag_equal(arb_radref(acb_realref(global.box.value())), arb_radref(acb_imagref(global.box.value()))) ||
        mag_is_zero(arb_radref(acb_realref(global.box.value()))) ||
        seed_scale_bits > (resources::precision_max - 64) / 2)
        return nullptr;
    Float global_radius, source_radius;
    arf_set_mag(global_radius.value(), arb_radref(acb_realref(global.box.value())));
    arf_one(source_radius.value());
    arf_mul_2exp_si(source_radius.value(), source_radius.value(), static_cast<slong>(seed_scale_bits));
    if (arf_cmp(global_radius.value(), source_radius.value()) > 0)
        return nullptr;
    constexpr uint64_t initial_bits = 64;
    // With source-bounded coordinates and a 2^-s width floor, O(s) bits
    // retain useful translated geometry without importing EA's precision.
    // A large center/width ratio that still cannot be represented is refused
    // by snap_center; the origin-centered global alternative remains valid.
    const slong geometry_bits = resources::precision_sum(resources::precision_product(seed_scale_bits, 2), 64);
    SubdivisionCompleteEntry entry;
    entry.normalization_bits = static_cast<uint64_t>(geometry_bits);
    const bool warm_complete = complete_warm_cover(global, warm, multiplicity);
    if (warm_complete && warm.size() == 1 && warm.front().count == multiplicity &&
        warm.front().multiplicity == multiplicity && acb_is_exact(warm.front().box.value()) &&
        acb_bits(warm.front().box.value()) <= geometry_bits) {
        auto search = std::make_unique<Subdivision>(warm, multiplicity, initial_bits, statistics,
                                                    SubdivisionPrecision::TrialLocal);
        search->impl_->complete_entry = true;
        entry.used_warm_cover = entry.exact_point = true;
        if (report) *report = entry;
        return search;
    }
    Float floor;
    arf_one(floor.value());
    arf_mul_2exp_si(floor.value(), floor.value(), -static_cast<slong>(seed_scale_bits));
    CompleteSquare square;
    if (!complete_square(square, global.box.value(), floor.value(), global.count, multiplicity, geometry_bits))
        return nullptr;
    if (warm_complete) {
        ComplexBall hull;
        acb_set_round(hull.value(), warm.front().box.value(), geometry_bits);
        for (size_t i = 1; i < warm.size(); ++i)
            acb_union(hull.value(), hull.value(), warm[i].box.value(), geometry_bits);
        CompleteSquare candidate;
        if (complete_square(candidate, hull.value(), floor.value(), global.count, multiplicity, geometry_bits) &&
            arf_cmp(arb_midref(candidate.half.value()), arb_midref(square.half.value())) < 0) {
            square = std::move(candidate);
            entry.used_warm_cover = true;
        }
    }
    entry.padded_width = square.padded;
    // Bypass the generic warm constructor's support inflation and retained
    // certificate convention. This is one native cell with known all-root
    // containment in R/2, a published count-n square inside 2R, and fresh speed.
    auto search = std::make_unique<Subdivision>(std::vector<Root>{}, multiplicity, initial_bits, statistics,
                                                SubdivisionPrecision::TrialLocal);
    auto &state = *search->impl_;
    Component component;
    Cell cell;
    cell.center = square.center;
    cell.half = square.half;
    component.cells.push_back(std::move(cell));
    component.precision = initial_bits;
    geometry(component, static_cast<slong>(initial_bits));
    component.certificate = std::move(square.root);
    component.certificate_center = std::move(square.center);
    component.certificate_radius = std::move(square.half);
    state.total_count = global.count;
    state.complete_entry = true;
    state.active.push_back(std::move(component));
    if (report) *report = entry;
    return search;
}
bool Subdivision::complete_entry() const { return impl_->complete_entry; }
SubdivisionPrecision Subdivision::precision_policy() const { return impl_->precision_policy; }
uint64_t Subdivision::cells() const {
    uint64_t n = 0;
    for (const auto &c : impl_->active)
        n += c.cells.size() + c.retained.size() + c.untested.size();
    for (const auto &c : impl_->finished)
        n += c.cells.size();
    return n;
}
uint64_t Subdivision::pending() const {
    return impl_->active.size();
}
bool Subdivision::suspended() const { return impl_->suspended; }
bool Subdivision::idle() const { return impl_->idle; }
double Subdivision::pending_work() const { return impl_->pending_work; }
void Subdivision::research_service_meter(const EnginePolicy *policy) {
    impl_->counting.research_service_policy = policy && policy->research_service_meter ? policy : nullptr;
}
void Subdivision::research_identity_shift_meter(bool enabled) {
    impl_->counting.research_identity_shift_meter = enabled;
}
double Subdivision::work() const {
    return impl_->work;
}
uint64_t Subdivision::min_bits() const {
    uint64_t p = 0;
    for (const auto &c : impl_->active)
        if (!p || c.working_precision() < p)
            p = c.working_precision();
    return p;
}
uint64_t Subdivision::max_bits() const {
    uint64_t p = 0;
    for (const auto &c : impl_->active)
        p = std::max(p, c.working_precision());
    return p;
}
void Subdivision::reopen(const cps_options &o) {
    for (auto it = impl_->finished.begin(); it != impl_->finished.end();) {
        if (!meets_goal(*it->certificate, o, static_cast<slong>(it->precision))) {
            it->blocked = false;
            it->state = CPS_COMPONENT_ACTIVE;
            impl_->active.push_back(std::move(*it));
            it = impl_->finished.erase(it);
        } else
            ++it;
    }
    for (auto &c : impl_->active)
        if (c.blocked && c.working_precision() < o.max_work_bits) {
            c.blocked = false;
            c.state = CPS_COMPONENT_ACTIVE;
        }
}
void Subdivision::progress(std::vector<cps_component_progress> &out, uint64_t factor) const {
    for (const auto &c : impl_->active)
        out.push_back({factor, c.working_precision(), c.cells.size() + c.retained.size() + c.untested.size(),
                       c.certificate ? c.certificate->count : 0, c.state});
    for (const auto &c : impl_->finished)
        out.push_back({factor, c.precision, c.cells.size(), c.certificate ? c.certificate->count : 0,
                       CPS_COMPONENT_FINISHED});
}
void Subdivision::precision_state(std::vector<SubdivisionPrecisionState> &out) const {
    const auto add = [&](const Component &c, bool finished) {
        SubdivisionPrecisionState state;
        state.base_bits = c.precision;
        state.trial_bits = c.trial_precision.value_or(0);
        state.newton_bits = c.newton_bits;
        state.count = c.certificate ? c.certificate->count : (c.before_split ? c.before_split->count : 0);
        state.splitting = c.splitting;
        state.blocked = c.blocked;
        state.finished = finished;
        state.certificate_radius = c.certificate_radius;
        out.push_back(std::move(state));
    };
    for (const auto &c : impl_->active) add(c, false);
    for (const auto &c : impl_->finished) add(c, true);
}
bool Subdivision::complete_cover(std::vector<Root> &out) const {
    std::vector<Root> roots;
    uint64_t count = 0;
    uint64_t precision = 0;
    auto add = [&](const Component &c) {
        if (!c.certificate)
            return false;
        roots.push_back(*c.certificate);
        count += roots.back().count;
        precision = std::max(precision, c.precision);
        return true;
    };
    for (const auto &c : impl_->active)
        if (!add(c))
            return false;
    for (const auto &c : impl_->finished)
        if (!add(c))
            return false;
    if (count != impl_->total_count || !disjoint(roots, static_cast<slong>(precision)))
        return false;
    out = std::move(roots);
    return true;
}
cps_stop_reason Subdivision::step(const cps_options &o, const Encloser &enclose, uint64_t cells_available,
                                 double allowance, NewtonAccuracy accuracy) {
    auto &s = *impl_;
    // Meter this step: its Taylor shifts and Graeffe steps cost what FLINT's
    // block multiplication costs on their operands (counting.cpp charges the
    // context, dovetail.hpp) and a Newton step about n, weighted at the
    // precision of the front component as the step begins. A step that
    // rotates past blocked components may act on another one; the weight is
    // a scheduling quantum (dovetail.hpp), so that approximation changes no
    // result.
    const double metered = s.counting.metered;
    const uint64_t newton = s.statistics.newton_attempts;
    int64_t bits = s.active.empty() ? 0 : static_cast<int64_t>(s.active.front().working_precision());
    struct Meter {
        Impl &s;
        double metered;
        uint64_t newton;
        int64_t bits;
        ~Meter() {
            const double n = static_cast<double>(s.multiplicity ? s.total_count / s.multiplicity : s.total_count);
            const double operations =
                (s.counting.metered - metered) + static_cast<double>(s.statistics.newton_attempts - newton) * n;
            s.work += operations * operation_weight(bits);
        }
    } meter{s, metered, newton, bits};
    // The allowance, in the counting context's unweighted operations at the
    // weight this step is charged at; lifted again when the step returns.
    s.suspended = s.idle = false;
    s.counting.allowance = std::isfinite(allowance) && bits > 0
                               ? s.counting.metered + std::max(0.0, allowance) / operation_weight(bits)
                               : std::numeric_limits<double>::infinity();
    s.counting.cancellation = &o;
    struct Lift {
        CountContext &counting;
        ~Lift() {
            counting.allowance = std::numeric_limits<double>::infinity();
            counting.cancellation = nullptr;
        }
    } lift{s.counting};
    // A disk test that stopped for the allowance leaves the component as it
    // was, at the front, for the next step to resume.
    const auto stopped = [&] {
        if (!s.counting.suspended)
            return false;
        s.suspended = true;
        s.idle = s.counting.metered == metered && s.statistics.newton_attempts == newton;
        s.pending_work = s.counting.pending * operation_weight(bits);
        return true;
    };
    uint64_t enclosure_ns = 0;
    ExclusiveTimer search_timer(s.statistics.search_ns, &s.statistics.counting_ns,
                                &s.statistics.taylor_shift_ns, &s.statistics.refinement_ns, &enclosure_ns);
    if (s.active.empty())
        return CPS_STOP_NONE;
    if (o.cancelled())
        return CPS_STOP_CANCELLED;
    auto rotate = [&] {
        s.active.push_back(std::move(s.active.front()));
        s.active.pop_front();
    };
    // A component waiting for cells must not hide other useful work in this
    // factor. Only stop after a full rotation finds no eligible component.
    bool cell_limited = false;
    size_t skipped = 0;
    for (; skipped < s.active.size(); ++skipped) {
        auto &front = s.active.front();
        if (!front.blocked) {
            const bool needs_children = front.splitting && front.untested.empty() &&
                                        front.split_index < front.cells.size();
            if (!needs_children || cells_available >= 4) break;
            front.state = CPS_COMPONENT_CELL_LIMIT;
            cell_limited = true;
        }
        rotate();
    }
    if (skipped == s.active.size())
        return cell_limited ? CPS_STOP_CELL_BUDGET : CPS_STOP_PRECISION;
    auto &c = s.active.front();
    c.state = CPS_COMPONENT_ACTIVE;
    const bool local_trial = s.precision_policy == SubdivisionPrecision::TrialLocal;
    if (local_trial) {
        // A blocked front component may have rotated away. Price operations
        // at the selected component's actual temporary precision.
        bits = static_cast<int64_t>(c.working_precision());
        meter.bits = bits;
        s.counting.allowance = std::isfinite(allowance) && bits > 0
                                  ? metered + std::max(0.0, allowance) / operation_weight(bits)
                                  : std::numeric_limits<double>::infinity();
        accuracy = NewtonAccuracy::PreparedExact;
    }
    const auto retire = [&] {
        if (c.certificate->multiplicity == c.certificate->count && acb_is_exact(c.certificate->box.value())) {
            // An exact point needs no retained search cells, even when the
            // explicit isolation path retires it as soon as it is certified.
            c.cells.clear();
            c.bounds = c.certificate->box;
        }
        c.trial_precision.reset();
        s.finished.push_back(std::move(c));
        s.active.pop_front();
    };
    // The explicit isolation service retires a certified singleton without
    // requested-bits raises, refinement, or another coefficient delivery.
    if (local_trial && o.goal == CPS_ROOT_ISOLATION && c.certificate &&
        meets_goal(*c.certificate, o, static_cast<slong>(c.precision))) {
        retire();
        return CPS_STOP_NONE;
    }
    if (c.certificate && c.certificate->count == c.certificate->multiplicity && o.goal != CPS_CLUSTER_RADIUS &&
        c.precision < o.max_work_bits && c.precision < o.bits + 32) {
        s.raise(c, o);
        rotate();
        return CPS_STOP_NONE;
    }
    if (c.splitting && c.untested.empty() && c.split_index < c.cells.size()) {
        const auto &parent = c.cells[c.split_index];
        std::vector<Cell> children;
        for (int dx : {-1, 1})
            for (int dy : {-1, 1}) {
                Cell child;
                child.half = parent.half;
                arb_mul_2exp_si(child.half.value(), child.half.value(), -1);
                Float delta;
                arf_mul_si(delta.value(), arb_midref(child.half.value()), dx, ARF_PREC_EXACT, ARF_RND_NEAR);
                int inexact = arf_add(arb_midref(acb_realref(child.center.value())),
                                      arb_midref(acb_realref(parent.center.value())), delta.value(),
                                      static_cast<slong>(resources::precision_sum(c.precision, 16)), ARF_RND_NEAR);
                arf_mul_si(delta.value(), arb_midref(child.half.value()), dy, ARF_PREC_EXACT, ARF_RND_NEAR);
                inexact |= arf_add(arb_midref(acb_imagref(child.center.value())),
                                   arb_midref(acb_imagref(parent.center.value())), delta.value(),
                                   static_cast<slong>(resources::precision_sum(c.precision, 16)), ARF_RND_NEAR);
                if (inexact) {
                    if (!s.raise(c, o))
                        s.abandon_split(c);
                    rotate();
                    return CPS_STOP_NONE;
                }
                fmpz_mul_2exp(child.x.value(), parent.x.value(), 1);
                fmpz_add_si(child.x.value(), child.x.value(), dx);
                fmpz_mul_2exp(child.y.value(), parent.y.value(), 1);
                fmpz_add_si(child.y.value(), child.y.value(), dy);
                children.push_back(std::move(child));
            }
        for (auto &child : children)
            c.untested.push_back(std::move(child));
        ++c.split_index;
        ++s.statistics.boxes_split;
        rotate();
        return CPS_STOP_NONE;
    }
    if (c.splitting && c.untested.empty()) {
        auto components = connected(std::move(c.retained), c);
        s.active.pop_front();
        for (auto &child : components)
            s.active.push_back(std::move(child));
        return CPS_STOP_NONE;
    }
    ComplexPolynomial p;
    {
        StageTimer timer(enclosure_ns);
        enclose(p.value(), static_cast<slong>(c.working_precision()));
    }
    if (c.splitting) {
        auto &child = c.untested.front();
        RealBall radius(child.half);
        // 3/2 of the half-width strictly contains the square's corners.
        arb_mul_ui(radius.value(), radius.value(), 3, static_cast<slong>(c.precision));
        arb_mul_2exp_si(radius.value(), radius.value(), -1);
        // Any positive certified disk count ends the exclusion attempt. It
        // does not certify occupancy of the square: retain the child and let
        // component counting establish its root set. Only zero excludes it.
        const auto result = pellet_count(p.value(), child.center.value(), arb_midref(radius.value()),
                                         static_cast<slong>(c.precision), subdivision_graeffe_steps(acb_poly_degree(p.value())), -1, &s.counting);
        if (stopped())
            return CPS_STOP_NONE;
        if (result.count == 0) {
            c.untested.pop_front();
            ++s.statistics.boxes_excluded;
        } else if (result.count < 0 && (result.needs_precision || result.width_limited) && s.raise(c, o)) {
            rotate();
            return CPS_STOP_NONE;
        } else if (result.count < 0 && result.exclusion_limited) {
            // At the precision cap, the direct zero-exclusion comparison is
            // the stop heuristic that avoids multiplying width-obstructed
            // cells. A broader all-index precision request is not sufficient
            // to abandon the split: geometry may still separate the roots.
            // Likewise retain a cell whose direct test was decisive but
            // whose later Graeffe comparisons became inconclusive.
            s.abandon_split(c);
            rotate();
            return CPS_STOP_NONE;
        } else {
            c.retained.push_back(std::move(child));
            c.untested.pop_front();
        }
        rotate();
        return CPS_STOP_NONE;
    }
    if (c.certificate && meets_goal(*c.certificate, o, static_cast<slong>(c.precision))) {
        retire();
        return CPS_STOP_NONE;
    }
    if (!c.certificate && s.separated(c)) {
        const auto count = pellet_count(p.value(), c.center.value(), arb_midref(c.radius.value()),
                                        static_cast<slong>(c.precision), subdivision_graeffe_steps(acb_poly_degree(p.value())), -1, &s.counting);
        if (stopped())
            return CPS_STOP_NONE;
        if (count.count == 0) {
            s.active.pop_front();
            return CPS_STOP_NONE;
        }
        Root root;
        bool certificate_needs_precision = false;
        if (count.count > 0 &&
            certify_cluster_box(root, p.value(), c.center.value(), arb_midref(c.radius.value()), count.count,
                                s.multiplicity, static_cast<slong>(c.precision), subdivision_graeffe_steps(acb_poly_degree(p.value())),
                                &s.counting, &certificate_needs_precision)) {
            c.certificate = std::move(root);
            c.certificate_center = c.center;
            c.certificate_radius = c.radius;
            if (local_trial && o.goal == CPS_ROOT_ISOLATION &&
                meets_goal(*c.certificate, o, static_cast<slong>(c.precision)))
                retire();
            else
                rotate();
            return CPS_STOP_NONE;
        }
        if (stopped())
            return CPS_STOP_NONE;
        // The inner disk may already be decided while uncertainty still
        // obstructs the outer disk that certifies the exported square.
        if ((count.needs_precision || certificate_needs_precision) && s.raise(c, o)) {
            rotate();
            return CPS_STOP_NONE;
        }
    }
    if (c.certificate) {
        const bool temporary = local_trial && c.certificate->count / s.multiplicity > 1;
        if (temporary && !c.trial_precision)
            c.trial_precision = c.precision;
        bool needs_precision = false;
        if (s.newton(c, p.value(), o, needs_precision, accuracy, static_cast<slong>(c.working_precision()))) {
            if (c.trial_precision) {
                c.precision = *c.trial_precision;
                c.trial_precision.reset();
            }
            rotate();
            return CPS_STOP_NONE;
        }
        if (stopped())
            return CPS_STOP_NONE;
        if (o.cancelled())
            return CPS_STOP_CANCELLED;
        if (temporary) {
            if (needs_precision) {
                if (!s.raise_trial(c, o)) {
                    // A cap is not a refusal. Keep P_try, e, r and the old
                    // certificate for an honest stop and later resumption.
                    c.blocked = true;
                    c.state = CPS_COMPONENT_PRECISION_LIMIT;
                }
                rotate();
                return CPS_STOP_NONE;
            }
            c.trial_precision.reset(); // definitive failure: split at P_base
        } else {
            const uint64_t attempted = c.newton_bits;
            c.newton_bits = std::max<uint64_t>(1, c.newton_bits / 2);
            // Retry failed contractions with fresh coefficient enclosures.
            if (needs_precision && s.raise(c, o)) {
                rotate();
                return CPS_STOP_NONE;
            }
            // At the precision budget, a contraction that failed says nothing
            // about a gentler one. The factor 2^-newton_bits doubles its exponent
            // on every success, so the attempt that fails is the first to reach
            // below what the precision resolves, or below the diameter of the
            // component's own roots: retry at half the exponent before bisecting.
            // Below the budget a failed step is retried at a raised precision
            // (above) or bisected, and the children's tests raise the precision
            // when width obstructs them; at the budget those tests are abandoned
            // instead, and the component stops for precision short of a
            // resolution a gentler step certifies. The exponent halves on every
            // failure and is not retried below least_newton_bits (a gentler step
            // could undo a bisection), so the retries end, and then the bisection
            // decides. History: docs/code-history-20261021.md, subdivision.cpp 1.
            if (c.precision == o.max_work_bits && attempted > least_newton_bits) {
                rotate();
                return CPS_STOP_NONE;
            }
            // A failed contraction is not a failure of the search: bisection below
            // separates the component's roots whenever Newton cannot contract it
            // (roots spanning its disc, or a contraction the ceiling cannot
            // certify). The search stops only for its budgets: a child count or
            // centre that max_work_bits cannot decide or represent (the split is
            // abandoned, PRECISION_LIMIT), the cell allowance, or a requested
            // accuracy beyond the precision budget.
            if (o.goal != CPS_CLUSTER_RADIUS && o.bits > o.max_work_bits) {
                c.blocked = true;
                c.state = CPS_COMPONENT_PRECISION_LIMIT;
                rotate();
                return CPS_STOP_NONE;
            }
        }
    }
    c.splitting = true;
    c.before_split = std::move(c.certificate);
    c.certificate.reset();
    rotate();
    return CPS_STOP_NONE;
}
} // namespace certpolysolve::detail
