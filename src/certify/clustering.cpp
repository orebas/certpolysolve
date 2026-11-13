// Candidate-guided cluster recovery (recover_clusters, declared in counting.hpp, as this file has no header): from one
// candidate point per root of a factor, a disjoint cover of that factor by squares whose root counts are proved.
//
// The candidates only guide the search; every count is proved on the polynomial by certify_cluster_box (counting.cpp:
// Pellet's test, with up to four Graeffe steps, on a disk and on the disk of twice its radius). Each hint starts a
// group with a power-of-two radius (geometry), lowered to the distance to the nearest other hint. In each of at most
// 16 rounds, overlapping groups merge and add their hint counts (merge_overlaps), and each uncertified group of k hints
// is tested for exactly k roots; after a failure its radius grows to 4 max_{j<k} (|q_j|/|q_k|)^(1/(k-j)), q the Taylor
// coefficients at its centre (increase_radius_from_lower_terms), or else doubles, held to the larger of twice the
// distance to the nearest other centre and twice the old radius (limit_growth_to_neighbors). The cover is accepted
// when the squares are disjoint and their counts sum to the degree times factor_multiplicity. Below max_work_bits,
// recovery returns false while a footprint misses an approximation query's goal (approximation_met), and returns false
// with *needs_precision set when a test is undecided by arithmetic width, so that the caller raises precision first.
//
// State: none; `out` changes only on success. Work guards: the degree at most heuristics::cluster_recovery_degree, and
// at most min(4096, max(64, 8 x degree)) tests.
// Entry points: recover_clusters, called only by isolate_complex (classical_isolation.cpp) when its inclusion test
// leaves roots unisolated and its allow_clusters argument is set.
// Design notes: docs/mathematics.md ("Local counts and candidate-guided recovery"),
// docs/cluster-precision-findings-20260924.md (why an arithmetic-limited count returns to precision escalation).
#include "certify/counting.hpp"
#include "core/telemetry.hpp"
#include <algorithm>
#include <optional>

namespace certpolysolve::detail {
namespace {
struct Group {
    ComplexBall center, footprint;
    RealBall radius; // Always an exact positive power of two.
    slong hints = 1;
    std::optional<Root> certificate;
};

void footprint(Group &group) {
    acb_set(group.footprint.value(), group.center.value());
    acb_add_error_arf(group.footprint.value(), arb_midref(group.radius.value()));
}

bool geometry(Group &group, const acb_t box, slong precision) {
    acb_get_mid(group.center.value(), box);
    if (!acb_is_finite(group.center.value()))
        return false;
    Float radius, component, scale;
    if (acb_is_finite(box)) {
        arf_set_mag(radius.value(), arb_radref(acb_realref(box)));
        arf_set_mag(component.value(), arb_radref(acb_imagref(box)));
        if (arf_cmp(component.value(), radius.value()) > 0)
            arf_set(radius.value(), component.value());
    }
    // Give point/degenerate hints a scale-aware seed radius. This is just a
    // geometric search choice and is never used as a root error estimate.
    acb_get_abs_ubound_arf(scale.value(), group.center.value(), precision);
    if (arf_is_zero(scale.value()))
        arf_one(scale.value());
    arf_mul_2exp_si(scale.value(), scale.value(), -std::max<slong>(8, precision / 2));
    if (arf_cmp(scale.value(), radius.value()) > 0)
        arf_set(radius.value(), scale.value());
    Integer exponent;
    arf_abs_bound_le_2exp_fmpz(exponent.value(), radius.value());
    arb_one(group.radius.value());
    arb_mul_2exp_fmpz(group.radius.value(), group.radius.value(), exponent.value());
    footprint(group);
    group.certificate.reset();
    return true;
}

bool merge_overlaps(std::vector<Group> &groups, slong precision) {
    // Greedy merging is deliberately heuristic. A merged geometric envelope
    // must be recertified from the original polynomial; counts are not guessed.
    bool changed;
    do {
        changed = false;
        for (size_t i = 0; i < groups.size(); ++i) {
            for (size_t j = i + 1; j < groups.size();) {
                if (!acb_overlaps(groups[i].footprint.value(), groups[j].footprint.value())) {
                    ++j;
                    continue;
                }
                ComplexBall united;
                acb_union(united.value(), groups[i].footprint.value(), groups[j].footprint.value(), precision);
                groups[i].hints += groups[j].hints;
                if (!geometry(groups[i], united.value(), precision))
                    return false;
                groups.erase(groups.begin() + static_cast<std::ptrdiff_t>(j));
                changed = true;
            }
        }
    } while (changed);
    return true;
}

bool increase_radius_from_lower_terms(Group &group, const acb_poly_t polynomial, slong precision,
                                      CountContext *context) {
    // A Newton-polygon-style proposal: balance lower Taylor terms against the
    // requested kth term. This is a SEARCH HEURISTIC, not a count certificate;
    // higher Taylor terms and all uncertainty still face the full Pellet test.
    // It avoids dozens of blind doublings when fixed input uncertainty makes
    // an approximate multiple root's Newton inclusion radius meaningless.
    ComplexPolynomial shifted;
    {
        uint64_t unused = 0;
        StageTimer timer(context ? context->statistics.taylor_shift_ns : unused);
        if (context)
            ++context->statistics.taylor_shifts;
        acb_poly_taylor_shift(shifted.value(), polynomial, group.center.value(), precision);
    }
    if (group.hints > acb_poly_degree(shifted.value()))
        return false;
    Float denominator, numerator, upper, maximum;
    acb_get_abs_lbound_arf(denominator.value(), shifted.value()->coeffs + group.hints, precision);
    if (!arf_is_finite(denominator.value()) || arf_sgn(denominator.value()) <= 0)
        return false;
    RealBall quotient;
    for (slong j = 0; j < group.hints; ++j) {
        acb_get_abs_ubound_arf(numerator.value(), shifted.value()->coeffs + j, precision);
        if (!arf_is_finite(numerator.value()))
            return false;
        arf_div(upper.value(), numerator.value(), denominator.value(), precision, ARF_RND_CEIL);
        arb_set_arf(quotient.value(), upper.value());
        arb_root_ui(quotient.value(), quotient.value(), static_cast<ulong>(group.hints - j), precision);
        if (!arb_is_finite(quotient.value()))
            return false;
        arb_get_ubound_arf(upper.value(), quotient.value(), precision);
        if (arf_cmp(upper.value(), maximum.value()) > 0)
            arf_set(maximum.value(), upper.value());
    }
    arf_mul_2exp_si(maximum.value(), maximum.value(), 2);
    if (arf_cmp(maximum.value(), arb_midref(group.radius.value())) <= 0)
        return false;
    Integer exponent;
    arf_abs_bound_le_2exp_fmpz(exponent.value(), maximum.value());
    arb_one(group.radius.value());
    arb_mul_2exp_fmpz(group.radius.value(), group.radius.value(), exponent.value());
    footprint(group);
    return true;
}

void limit_growth_to_neighbors(Group &group, const std::vector<Group> &groups, const arf_t previous_radius,
                               slong precision) {
    Float nearest, distance, minimum_growth;
    arf_pos_inf(nearest.value());
    ComplexBall difference;
    for (const auto &other : groups) {
        if (&group == &other)
            continue;
        acb_sub(difference.value(), group.center.value(), other.center.value(), precision);
        acb_get_abs_ubound_arf(distance.value(), difference.value(), precision);
        if (arf_cmp(distance.value(), nearest.value()) < 0)
            arf_set(nearest.value(), distance.value());
    }
    if (!arf_is_finite(nearest.value()))
        return;
    arf_mul_2exp_si(nearest.value(), nearest.value(), 1);
    arf_mul_2exp_si(minimum_growth.value(), previous_radius, 1);
    if (arf_cmp(nearest.value(), minimum_growth.value()) < 0)
        arf_set(nearest.value(), minimum_growth.value());
    if (arf_cmp(nearest.value(), arb_midref(group.radius.value())) >= 0)
        return;
    Integer exponent;
    arf_abs_bound_le_2exp_fmpz(exponent.value(), nearest.value());
    arb_one(group.radius.value());
    arb_mul_2exp_fmpz(group.radius.value(), group.radius.value(), exponent.value());
    footprint(group);
}
} // namespace

bool recover_clusters(std::vector<Root> &out, const acb_poly_t polynomial, acb_srcptr hints, slong hint_count,
                      uint64_t factor_multiplicity, const cps_options &options, slong precision,
                      CountContext *context, bool *needs_precision) {
    if (needs_precision)
        *needs_precision = false;
    uint64_t unused = 0;
    ExclusiveTimer timer(context ? context->statistics.refinement_ns : unused,
                         context ? &context->statistics.counting_ns : nullptr,
                         context ? &context->statistics.taylor_shift_ns : nullptr);
    const slong degree = acb_poly_degree(polynomial);
    // Work guards for the first candidate-guided implementation. A skipped or
    // exhausted recovery leaves the caller's global coverage certificate intact.
    if (!hints || degree < 1 || hint_count != degree || degree > heuristics::cluster_recovery_degree)
        return false;
    const size_t max_tests = std::min<size_t>(4096, std::max<size_t>(64, 8 * static_cast<size_t>(degree)));
    size_t tests = 0;
    std::vector<Group> groups;
    groups.reserve(static_cast<size_t>(degree));
    for (slong j = 0; j < degree; ++j) {
        Group group;
        if (!geometry(group, hints + j, precision))
            return false;
        groups.push_back(std::move(group));
    }
    // Individual Newton inclusion radii may explode at an uncertain collision.
    // Cap initial search neighborhoods by nearby candidate spacing, while still
    // retaining a positive scale-aware seed for coincident points. This changes
    // only candidate geometry; no original input or certificate is narrowed.
    for (size_t i = 0; i < groups.size(); ++i) {
        Float nearest, distance;
        arf_pos_inf(nearest.value());
        ComplexBall difference;
        for (size_t j = 0; j < groups.size(); ++j) {
            if (i == j)
                continue;
            acb_sub(difference.value(), groups[i].center.value(), groups[j].center.value(), precision);
            acb_get_abs_ubound_arf(distance.value(), difference.value(), precision);
            if (arf_cmp(distance.value(), nearest.value()) < 0)
                arf_set(nearest.value(), distance.value());
        }
        if (arf_is_finite(nearest.value()) && arf_cmp(nearest.value(), arb_midref(groups[i].radius.value())) < 0) {
            ComplexBall bounded;
            acb_set(bounded.value(), groups[i].center.value());
            acb_add_error_arf(bounded.value(), nearest.value());
            if (!geometry(groups[i], bounded.value(), precision))
                return false;
        }
    }
    for (unsigned round = 0; round < 16; ++round) {
        if (options.cancelled() || !merge_overlaps(groups, precision))
            return false;
        // Early recovery is optional. Count only when the current footprints
        // are already on the requested approximation scale; otherwise retain
        // candidate iteration instead of spending on a coarse intermediate
        // cover. At the precision ceiling, still recover any useful coverage.
        if (options.goal == CPS_ROOT_APPROXIMATION &&
            precision < static_cast<slong>(options.max_work_bits))
            for (const auto &group : groups) {
                Root proposed;
                proposed.box = group.footprint;
                if (!approximation_met(proposed, options.bits, options.absolute_tolerance, precision))
                    return false;
            }
        // A footprint is only a proposal. Every accepted count below still
        // requires the original polynomial's Pellet certificate.
        bool complete = true;
        for (auto &group : groups) {
            if (group.certificate)
                continue;
            if (options.cancelled() || tests++ >= max_tests)
                return false;
            Root root;
            bool arithmetic_limited = false;
            if (certify_cluster_box(root, polynomial, group.center.value(), arb_midref(group.radius.value()),
                                    group.hints, factor_multiplicity, precision, default_graeffe_steps, context,
                                    &arithmetic_limited)) {
                acb_set(group.footprint.value(), root.box.value());
                group.certificate = std::move(root);
            } else {
                // A failed interval comparison is not evidence that the disk
                // needs to grow. Let the owner improve the source enclosure
                // before spending on more counts with the same uncertainty.
                // At the precision ceiling, retain geometric recovery.
                if (arithmetic_limited && precision < static_cast<slong>(options.max_work_bits)) {
                    if (needs_precision)
                        *needs_precision = true;
                    return false;
                }
                complete = false;
                Float previous_radius;
                arf_set(previous_radius.value(), arb_midref(group.radius.value()));
                if (!increase_radius_from_lower_terms(group, polynomial, precision, context)) {
                    arb_mul_2exp_si(group.radius.value(), group.radius.value(), 1);
                    footprint(group);
                }
                // Let nearby hints merge before treating an almost-zero first
                // derivative as evidence for a huge one-root neighborhood.
                limit_growth_to_neighbors(group, groups, previous_radius.value(), precision);
            }
        }
        if (complete) {
            std::vector<Root> certified;
            uint64_t count = 0;
            for (const auto &group : groups) {
                certified.push_back(*group.certificate);
                count += certified.back().count;
            }
            if (disjoint(certified) && count == static_cast<uint64_t>(degree) * factor_multiplicity) {
                for (auto &root : certified)
                    out.push_back(std::move(root));
                return true;
            }
            // Export may have enlarged a rectangle. Its actual footprint will
            // be merged at the next round; no pre-export disjointness is assumed.
        }
    }
    return false;
}
} // namespace certpolysolve::detail
