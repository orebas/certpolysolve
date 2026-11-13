// The certificate attempt the recipes share (machine_newton_cover): a binary64 Aberth-Ehrlich packet on the whole
// input (prepare_machine_engine) in quanta while it settles points, then a certificate at each point: Newton steps on
// rungs of doubling precision and Rouche's test of the linear Taylor model on a goal-sized rectangle about the point
// c or its Newton step (certify_at_step), the remainder bounded by the coefficient majorant of p'' (linear_remainder);
// where that is too coarse, the inclusion disk of radius n|p(c)/p'(c)|. Rounding refusals retry at doubled
// precision; repair rounds move refused points by multiprecision Aberth packets, certified ones held
// (candidate_packet). Success is n pairwise disjoint regions; a refusal can return its certificates
// (PartialNewtonCover). Its spending (newton_cover_costs) is capped by the caller's route's cost
// (fallback_operations). certify_newton_point certifies one point the same way on a monic source. Before the ball
// certificate, a point first tries the compensated tier (compensated_certificate.hpp): the same Rouché test from one
// evaluation in binary64 with error-free transformations, about a hundredth of the ball certificate's cost.
//
// State: none across calls; each call owns its packet, points and certificates.
// Entry points: machine_newton_cover (classical_isolation.cpp; secular_refinement.cpp, with `settled` and
// PartialNewtonCover; prepared_piecewise.cpp, before piecewise); certify_newton_point (tests and benchmark tools);
// newton_cover_costs (here).
// Not here: Rouche's test on a model (certify_root_model) and the rectangle's first half-width (goal_half_width),
// chart_certificate.cpp; the binary64 packet (machine_candidates.cpp); the multiprecision packet (candidates.cpp); the
// disk evaluators (disk_*_polynomial_evaluator.hpp).
// Design notes: docs/mathematics.md ("The certificate attempt"), docs/production-methods.md (the recipes).
// Moved verbatim from piecewise_roots.cpp in W8's S3 (docs/w8-plan-20261021.md).
#include "certify/attempts/newton_cover.hpp"
#include "core/engine_policy.hpp"
#include "core/engine_statistics.hpp"
#include "certify/attempts/newton_evaluation.hpp"
#include "certify/newton_point.hpp"
#include "numeric/candidates/candidates.hpp"
#include "numeric/arithmetic/arf_binary64.hpp"
#include "certify/attempts/compensated_certificate.hpp"
#include "numeric/evaluate/disk_sparse_polynomial_evaluator.hpp"
#include "core/dovetail.hpp"
#include "core/heuristics.hpp"
#include "numeric/arithmetic/machine_arithmetic.hpp"
#include "numeric/candidates/machine_candidates.hpp"
#include "certify/chart_certificate.hpp"
#include "numeric/arithmetic/modulus_bound.hpp"
#include "numeric/models/root_model.hpp"
#include "numeric/evaluate/structured_source.hpp"
#include "core/telemetry.hpp"
#include "core/wide_int.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <numbers>
#include <optional>
#include <span>
#include <vector>

namespace certpolysolve::detail {
namespace {
// What machine_newton_cover's work costs in the selector's operations at its first precision (heuristics.hpp):
// one certificate evaluation (value and slope) of a source of degree n, dense (sparse_support empty: 2(n + 1)
// block coefficients) or sparse (its support given, ascending: its terms' two dot products and its powers'
// products); and one sweep at one point of a repair round: n - 1 interactions and the evaluation of the
// coefficient enclosure (heuristics::newton_cover_interaction_weight, newton_cover_coefficient_weight).
struct NewtonCoverCosts {
    uint64_t evaluation = 0, round_sweep = 0; // counts of operations, exact
};
NewtonCoverCosts newton_cover_costs(slong n, std::span<const slong> sparse_support) {
    // Exact counts; for a degree of at least 1, as machine_newton_cover asks.
    const uint64_t coefficient = heuristics::certificate_coefficient_weight;
    const auto degree = resources::nonnegative(n);
    NewtonCoverCosts costs;
    costs.evaluation =
        sparse_support.empty()
            ? coefficient * 2 * (degree + 1)
            : coefficient * 2 * sparse_support.size() + heuristics::certificate_multiplication_weight *
                                                            DiskSparsePolynomialEvaluator::plan(sparse_support).multiplications;
    costs.round_sweep = heuristics::newton_cover_interaction_weight * (degree - 1) +
                        (sparse_support.empty() ? heuristics::newton_cover_coefficient_weight * 2 * (degree + 1)
                                                : costs.evaluation);
    return costs;
}

// log2 of a nonzero finite magnitude as an exponent and a fraction in [-1, 0), the fraction from the mantissa alone: an
// exact scaling of the input by 2^t then moves only the exponent. The solve's answers are equivariant under such
// scalings (test_equivalence_properties); mag_get_d_log2_approx rounds to the exponent beyond 2^20 and is not.
bool magnitude_log2(const mag_t m, slong &exponent, double &fraction) {
    if (mag_is_zero(m) || mag_is_inf(m) || COEFF_IS_MPZ(MAG_EXP(m)))
        return false;
    exponent = MAG_EXP(m);
    fraction = std::log2(std::ldexp(static_cast<double>(MAG_MAN(m)), -MAG_BITS));
    return true;
}
// Proposal geometry only: every pair of disks in a selected component meets.
// This does not prove that a component contains distinct roots or a root count.
// Keep admission, cached evaluation, union order and member order unchanged;
// the caller alone decides the work allowance, starts packets and commits proofs.
std::vector<std::vector<slong>> cluster_start_groups(ConstRootView points, std::span<const slong> pool,
                                                   const CompensatedCertificate &compensated,
                                                   std::span<const CompensatedEvaluation> seen,
                                                   const cps_options &options) {
    // Each refused point's disk radius n |p/p'|, an upper bound from the compensated tier's evaluation at the
    // binary64 point nearest it: the one the tier made when it refused the point there, or a new one where the
    // tier's admission (may_certify) would make it. A point the admission refuses takes no part: the
    // coefficients' own errors there exceed what the tier's test accepts (chebyshev-8192's degree-4,096
    // factor: all 4,096 points refused, and new evaluations at them, 0.45 s of 26, bounded no slope away
    // from 0); nor does one whose slope's bound holds 0. Ball Horner at the certificate's precision instead
    // cost about a tenth of orthogonal-laguerre-d032's solve on its 25 refused points.
    std::vector<slong> member;
    std::vector<double> reach, at_re, at_im;
    for (const slong i : pool) {
        const double re = arf_to_binary64(arb_midref(acb_realref(points[i])), ARF_RND_NEAR);
        const double im = arf_to_binary64(arb_midref(acb_imagref(points[i])), ARF_RND_NEAR);
        const CompensatedEvaluation &kept = seen[static_cast<size_t>(i)];
        const bool kept_here = kept.ok && kept.x_re == re && kept.x_im == im;
        if (!kept_here && !compensated.may_certify(re, im, options))
            continue;
        const CompensatedEvaluation e = kept_here ? kept : compensated.evaluate(re, im);
        const double slope = std::hypot(e.slope_re, e.slope_im) - e.slope_error;
        if (!e.ok || !(slope > 0))
            continue;
        const double radius = (std::hypot(e.value_re + e.correction_re, e.value_im + e.correction_im) +
                               e.value_error) / slope * static_cast<double>(points.size());
        if (!std::isfinite(radius))
            continue;
        member.push_back(i);
        reach.push_back(radius);
        at_re.push_back(re);
        at_im.push_back(im);
    }
    const size_t with_disk = member.size();
    if (with_disk < 2)
        return {};
    // Groups: single linkage on meeting disks.
    std::vector<size_t> group_of(with_disk);
    for (size_t a = 0; a < with_disk; ++a)
        group_of[a] = a;
    const auto find = [&](size_t a) {
        while (group_of[a] != a)
            a = group_of[a] = group_of[group_of[a]];
        return a;
    };
    // Two disks meet, at the binary64 points the radii were taken at (a grouping, nothing proved by it).
    const auto meet = [&](size_t a, size_t b) {
        const double dr = at_re[a] - at_re[b], di = at_im[a] - at_im[b], r = reach[a] + reach[b];
        return dr * dr + di * di <= r * r;
    };
    for (size_t a = 0; a < with_disk; ++a)
        for (size_t b = a + 1; b < with_disk; ++b)
            if (meet(a, b))
                group_of[find(a)] = find(b);
    std::vector<std::vector<size_t>> members(with_disk);
    for (size_t a = 0; a < with_disk; ++a)
        members[find(a)].push_back(a);
    // A group is a cluster when every two of its disks meet. Near a cluster of k roots about c, far from the
    // others, |p/p'| is about |z - c| / k, so a point's disk (radius n |p/p'|) reaches (n/k) |z - c| and two
    // points z, w are |z - w| <= |z - c| + |w - c| apart: their disks meet. A chain of separate roots whose
    // neighbours' disks meet is no cluster, and the circle about its centre no start (binary64's points on
    // orthogonal-laguerre-d032: 25 of 32 in one chain, 189 of its 300 pairs meeting).
    std::vector<std::vector<slong>> groups;
    for (const auto &group : members) {
        bool cluster = group.size() >= 2;
        for (size_t x = 0; cluster && x < group.size(); ++x)
            for (size_t y = x + 1; cluster && y < group.size(); ++y)
                cluster = meet(group[x], group[y]);
        if (!cluster)
            continue;
        groups.emplace_back();
        for (const size_t a : group)
            groups.back().push_back(member[a]);
    }
    return groups;
}

// Geometry for one proposed cluster start, not a root certificate. The caller
// supplies a group of at least two valid point indices admitted by the grouping
// policy. Source snapshots and point views are borrowed only during sizing.
struct ClusterStartCircle {
    ComplexBall centre;
    slong radius_exponent;
    double radius_fraction;
    double spread;
    slong packet_bits;
};

// Keep every weighted-work addition at its evaluation site. Returning a summed
// cost would regroup floating-point accounting. This cannot move proposals,
// publish certificates, or change the attempt's working precision.
std::optional<ClusterStartCircle> size_cluster_start_circle(
    PolynomialSource &source, ConstRootView points, std::span<const slong> group,
    const cps_options &options, slong precision, double &spent) {
    const slong k = static_cast<slong>(group.size());
    // The centre: Newton's iteration on p^(k-1) at the certificate's precision. Each step evaluates p^(k-1)
    // and p^(k), of degrees n - k + 1 and n - k, together no more than one of the attempt's evaluations (the
    // value and the slope of p), and is charged one; so is each pair p(c), p^(k)(c) below. The radius: p(c) is about
    // rho^k |q_k|, so until its enclosure has a few correct bits rho is the evaluation's noise; p(c) and
    // p^(k)(c) are evaluated at a precision doubling from the certificate's, the centre refined by one
    // Newton step at each (a triple 2^-60 apart needs about 3 x 60 bits beyond the leading ones).
    // At the certificate's precision as it now stands: a group certified before may have raised it.
    const auto enclosed_here = source.coefficient_enclosure(precision);
    ComplexBall centre;
    // rho = 2^(rho_exponent + rho_fraction), rho_fraction in [0, 1): a scaling of the input by 2^t moves
    // only the integer (magnitude_log2), so the seeds scale exactly.
    slong rho_exponent = 0;
    double rho_fraction = 0;
    bool sized = false;
    {
        const slong work = precision;
        ComplexPolynomial lower, upper_d;
        acb_poly_set(lower.value(), enclosed_here->polynomial.value());
        for (slong j = 0; j < k - 1; ++j)
            acb_poly_derivative(lower.value(), lower.value(), work);
        acb_poly_derivative(upper_d.value(), lower.value(), work);
        acb_zero(centre.value());
        for (const slong i : group)
            acb_add(centre.value(), centre.value(), points[i], work);
        acb_div_si(centre.value(), centre.value(), k, work);
        acb_get_mid(centre.value(), centre.value());
        // A group where a centred source evaluates the uncentred polynomial (PolynomialSource::origin_point)
        // stays with the routes that follow that rule: the attempt's certificates evaluate the centred
        // expansion, which there pays the translate's cancellation (test_coordinate_policy: 16 roots of
        // modulus 2^-40 beside 48 about 2, centred at about 3/2, 640 bits).
        ComplexBall original;
        if (source.origin_point(original.value(), centre.value()))
            return std::nullopt;
        double last = std::numeric_limits<double>::infinity();
        bool centred = true;
        for (unsigned step = 0; step < heuristics::newton_cover_steps(options.bits); ++step) {
            ComplexBall v, d, corr;
            acb_poly_evaluate(v.value(), lower.value(), centre.value(), work);
            acb_poly_evaluate(d.value(), upper_d.value(), centre.value(), work);
            spent += operation_weight(work);
            if (acb_contains_zero(d.value())) {
                centred = false;
                break;
            }
            acb_div(corr.value(), v.value(), d.value(), work);
            Magnitude m;
            acb_get_mag(m.value(), corr.value());
            const double size = mag_get_d(m.value());
            acb_sub(centre.value(), centre.value(), corr.value(), work);
            acb_get_mid(centre.value(), centre.value());
            if (!(size < last / 2) && step > 0)
                break;
            last = size;
        }
        for (slong at = work; centred && acb_is_finite(centre.value()) &&
                              at <= static_cast<slong>(options.max_work_bits);
             at = resources::precision_product(at, 2)) {
            const auto enclosed_at = at == work ? enclosed_here : source.coefficient_enclosure(at);
            ComplexPolynomial lower_at, upper_at;
            acb_poly_set(lower_at.value(), enclosed_at->polynomial.value());
            for (slong j = 0; j < k - 1; ++j)
                acb_poly_derivative(lower_at.value(), lower_at.value(), at);
            acb_poly_derivative(upper_at.value(), lower_at.value(), at);
            if (at > work) {
                ComplexBall v, d, corr;
                acb_poly_evaluate(v.value(), lower_at.value(), centre.value(), at);
                acb_poly_evaluate(d.value(), upper_at.value(), centre.value(), at);
                spent += operation_weight(at);
                if (acb_contains_zero(d.value()))
                    break;
                acb_div(corr.value(), v.value(), d.value(), at);
                acb_sub(centre.value(), centre.value(), corr.value(), at);
                acb_get_mid(centre.value(), centre.value());
            }
            ComplexBall q0, qk;
            acb_poly_evaluate(q0.value(), enclosed_at->polynomial.value(), centre.value(), at);
            acb_poly_evaluate(qk.value(), upper_at.value(), centre.value(), at);
            spent += operation_weight(at);
            if (acb_rel_accuracy_bits(q0.value()) < heuristics::newton_cover_correction_bits ||
                acb_contains_zero(qk.value()))
                continue;
            Magnitude m0, mk;
            acb_get_mag(m0.value(), q0.value());
            acb_get_mag_lower(mk.value(), qk.value());
            slong e0 = 0, ek = 0;
            double f0 = 0, fk = 0;
            if (!magnitude_log2(m0.value(), e0, f0) || !magnitude_log2(mk.value(), ek, fk))
                break;
            double log_factorial = 0;
            for (slong j = 2; j <= k; ++j)
                log_factorial += std::log2(static_cast<double>(j));
            // log2 rho = (e0 - ek + fraction) / k, the fraction free of any scaling; its whole part joins
            // the integers, which divide by k with a remainder.
            const double fraction = f0 - fk + log_factorial, whole = std::floor(fraction);
            const slong total = e0 - ek + static_cast<slong>(whole);
            rho_exponent = total >= 0 ? total / k : -((-total + k - 1) / k);
            rho_fraction = (static_cast<double>(total - k * rho_exponent) + (fraction - whole)) /
                           static_cast<double>(k);
            sized = true;
            break;
        }
    }
    if (!sized)
        return std::nullopt;
    Magnitude mc;
    acb_get_mag(mc.value(), centre.value());
    slong ec = 0;
    double fc = 0;
    if (!magnitude_log2(mc.value(), ec, fc))
        return std::nullopt;
    // The packet's precision: it resolves the cluster (2^(-P/k) below its relative size 2^-L: P > k L) and
    // brings each root to the goal (Aberth's error about 2^-P 2^((k-1) L), below 2^-(goal + margin)).
    // L = log2(|c| / rho), its integer part a difference of exponents that a scaling moves alike. A group
    // as wide as its centre's modulus has no local model.
    const double spread = static_cast<double>(ec - rho_exponent) + (fc - rho_fraction);
    if (!(spread > 0))
        return std::nullopt;
    const double wanted = std::max({static_cast<double>(precision),
                                    static_cast<double>(k) * spread,
                                    static_cast<double>(k - 1) * spread + static_cast<double>(options.bits) +
                                        heuristics::double_double_margin_bits});
    // A cluster that needs the caller's whole precision is left to the routes that work within it: at
    // that precision every refusal is the source's limit.
    const slong packet_bits = static_cast<slong>(std::ceil(wanted / FLINT_BITS)) * FLINT_BITS;
    if (static_cast<uint64_t>(packet_bits) >= options.max_work_bits)
        return std::nullopt;
    return ClusterStartCircle{std::move(centre), rho_exponent, rho_fraction, spread, packet_bits};
}

// Compute the next precision from this rounding refusal. This calculation has
// no retry, cancellation, cache or work-accounting effects; the caller decides
// whether the allowance permits the raise and commits it only at that point.
slong certificate_rounding_precision(slong bits, uint64_t goal_bits, uint64_t maximum_bits,
                                    const mag_t uncertainty, const acb_t derivative,
                                    const arb_t half_width) {
    // The raise: to one limb above the deficit the evaluation shows (the value's radius falls as 2^-bits and
    // the test needs it below the linear term's size, |p'(c)| times the half-width; the limb covers the
    // test's own factors), and at least heuristics::certificate_raise_ratio, at most a doubling, rounded up
    // to whole limbs. The doubling it replaces overshoots the need by up to 2 for every later point
    // (laguerre-4000: 10,880 bits where none needs more than 6,938). Raises of a limb at a time, where the
    // deficit is small, took a point needing about a doubling dozens of raises and re-evaluations
    // (unisolve_mand255 at 3,322 bits: 1.6 s to the time limit, 28 October 2026).
    slong raised = resources::precision_product(bits, 2);
    {
        Float lower, allowed;
        acb_get_abs_lbound_arf(lower.value(), derivative, bits);
        RealBall linear_term;
        arb_set_arf(linear_term.value(), lower.value());
        arb_mul(linear_term.value(), linear_term.value(), half_width, bits);
        arb_get_lbound_arf(allowed.value(), linear_term.value(), bits);
        Magnitude allowed_size;
        arf_get_mag_lower(allowed_size.value(), allowed.value());
        const bool measured = mag_is_finite(uncertainty) && !mag_is_zero(uncertainty) &&
                              arf_sgn(allowed.value()) > 0 && arf_is_finite(allowed.value()) &&
                              !mag_is_zero(allowed_size.value());
        if (measured) {
            const double deficit = std::max(
                0.0, std::ceil(mag_get_d_log2_approx(uncertainty) -
                               mag_get_d_log2_approx(allowed_size.value())));
            const double wanted =
                std::max(static_cast<double>(bits) * heuristics::certificate_raise_ratio,
                         static_cast<double>(bits) + deficit + FLINT_BITS);
            // Below two rungs of the Newton ladder (heuristics::newton_ladder_first_rung) a refusal without
            // a deficit keeps the doubling: a raise between the rungs leaves one below the certificate, and
            // the points take their Newton steps at the certificate's precision (Chebyshev and Legendre
            // polynomials of degree 320 raised from 170 to 256 bits: 13 to 23% slower, 28 October 2026).
            const bool ladder = bits >= 2 * heuristics::newton_ladder_first_rung || deficit > 0;
            // Where the goal is at least half the precision the doubling stays too: the points reach the
            // goal by Newton's steps at the certificate's precision, and the room above it (the precision
            // less the goal and the evaluation's loss) is what the doubling widens. Raised to the deficit
            // there, 3,322-bit jobs certified at 4,800 bits where they took 6,708, and three of five were
            // slower (bernoulli_256 1.35, chrmc_d171 1.18, chrmc343 1.09; 28 October 2026).
            const bool loss_bound = 2 * goal_bits <= static_cast<uint64_t>(bits);
            if (ladder && loss_bound && wanted < static_cast<double>(raised)) {
                const slong limbs = static_cast<slong>(std::ceil(wanted / FLINT_BITS));
                raised = std::min(raised, std::max<slong>(limbs * FLINT_BITS,
                                                          resources::precision_sum(bits, FLINT_BITS)));
            }
        }
    }
    raised = std::min<slong>(static_cast<slong>(maximum_bits), raised);
    return raised;
}

// The proof strength of a point result. RootInclusion proves presence only;
// exact singleton counts follow when a complete disjoint cover accounts for
// the degree. PartialNewtonCover retains only that presence guarantee.
enum class PointProof { Refused, CompensatedSingleton, Singleton, RootInclusion };

// Progress shared by the ordinary point attempts and their controller. Work
// is accumulated at each original evaluation (not regrouped on return). The
// source limit is a persistent latch; cluster starts save and restore it.
struct NewtonAttemptProgress {
    slong precision;
    double spent = 0;
    bool source_limited = false;
};

// Improve and prove one proposed point. No point array, retained cover, packet
// state or publication crosses this boundary. The caller owns rollback after
// a refusal; a failed preparation never permits evaluation at that precision.
PointProof attempt_ball_point(Root &root, const acb_t point, NewtonCoverEvaluation &evaluations,
                              const cps_options &options, slong degree, double raise_allowance,
                              NewtonAttemptProgress &progress, StudyStatistics *counts) {
    const double entry = progress.spent;
    Magnitude previous_uncertainty, previous_correction;
    mag_inf(previous_uncertainty.value());
    ComplexBall center;
    acb_get_mid(center.value(), point);
    unsigned steps = 0;
    // Newton at doubling precision. A step from a point with b correct bits
    // needs about 2b bits of evaluation to give 2b, so the steps below the
    // certificate precision evaluate on rungs doubling from
    // heuristics::newton_ladder_first_rung, and only the last doublings and
    // the certificate itself (below, unchanged) are paid at the certificate
    // precision: from a binary64 point to a 3,322-bit goal, five rung
    // evaluations cost about one at the goal, where six steps at the goal
    // would cost six. Nothing is certified on a rung; a correction a rung
    // leaves without enough bits moves to the next rung, and one that does
    // not contract fourfold ends the ladder (the steps below judge the
    // point as before).
    // History: docs/code-history-20261021.md, piecewise_roots.cpp 3.
    for (slong rung = heuristics::newton_ladder_first_rung;
         rung < progress.precision && steps < heuristics::newton_cover_steps(options.bits);
         rung = std::min<slong>(progress.precision, resources::precision_product(rung, 2))) {
        ComplexBall value, derivative, correction;
        if (options.cancelled() || !evaluations.rung_evaluate(rung, value.value(), derivative.value(), center.value()))
            break;
        progress.spent += operation_weight(rung);
        if (counts)
            ++counts->machine_cover_rung_evaluations;
        if (acb_contains_zero(derivative.value()))
            break;
        acb_div(correction.value(), value.value(), derivative.value(), rung);
        if (!acb_is_finite(correction.value()) ||
            acb_rel_accuracy_bits(correction.value()) < heuristics::newton_cover_correction_bits)
            continue;
        Magnitude size, bound;
        acb_get_mag(size.value(), correction.value());
        mag_mul_2exp_si(bound.value(), previous_correction.value(), heuristics::newton_cover_contraction_log2);
        if (steps && mag_cmp(size.value(), bound.value()) > 0)
            break;
        acb_sub(center.value(), center.value(), correction.value(), rung);
        acb_get_mid(center.value(), center.value());
        if (!acb_is_finite(center.value()) || acb_is_zero(center.value()))
            return PointProof::Refused;
        mag_set(previous_correction.value(), size.value());
        ++steps;
        if (counts)
            ++counts->machine_cover_newton_steps;
    }
    for (;;) {
        const slong bits = progress.precision;
        if (counts)
            counts->machine_cover_peak_bits =
                std::max(counts->machine_cover_peak_bits, static_cast<uint64_t>(bits));
        ComplexBall value, derivative, correction, linear, zero;
        if (options.cancelled() ||
            !evaluations.evaluate(value.value(), derivative.value(), center.value()))
            return PointProof::Refused;
        progress.spent += operation_weight(bits);
        if (counts && evaluations.has_disk_evaluator())
            ++counts->machine_cover_block_evaluations;
        RootModel model;
        model.center = center;
        RealBall half_width;
        goal_half_width(half_width.value(), model.center.value(), options, bits);
        // The raise (a rounding refusal): false where the point is refused, true where the certificate is tried
        // again at the raised precision.
        const auto raise = [&]() -> bool {
            // Retry the point at a higher precision: the refusal is rounding.
            // No root is accepted until a complete certificate passes.
            Magnitude uncertainty;
            mag_add(uncertainty.value(), arb_radref(acb_realref(value.value())),
                    arb_radref(acb_imagref(value.value())));
            if (static_cast<uint64_t>(bits) >= options.max_work_bits ||
                mag_cmp(uncertainty.value(), previous_uncertainty.value()) >= 0) {
                // Rounding that precision does not reduce (the source's own
                // uncertainty, or the caller's precision budget): no later
                // round of this attempt can certify the point either.
                progress.source_limited = true;
                return false;
            }
            const slong raised = certificate_rounding_precision(
                bits, options.bits, options.max_work_bits, uncertainty.value(), derivative.value(), half_width.value());
            if (progress.spent - entry + operation_weight(raised) > raise_allowance)
                return false;
            previous_uncertainty = uncertainty;
            progress.precision = raised;
            if (counts)
                ++counts->machine_cover_precision_retries;
            if (!evaluations.prepare(progress.precision))
                return false;
            return true;
        };
        // The Newton step this evaluation offers (heuristics.hpp has the
        // rules): the next round judges the stepped point by the same
        // certificate.
        bool steppable = !acb_contains_zero(derivative.value());
        if (steppable) {
            acb_div(correction.value(), value.value(), derivative.value(), bits);
            steppable = acb_is_finite(correction.value()) &&
                        acb_rel_accuracy_bits(correction.value()) >= heuristics::newton_cover_correction_bits;
        }
        const auto take_step = [&] {
            if (!steppable || steps == heuristics::newton_cover_steps(options.bits))
                return false;
            Magnitude size, bound;
            acb_get_mag(size.value(), correction.value());
            mag_mul_2exp_si(bound.value(), previous_correction.value(), heuristics::newton_cover_contraction_log2);
            if (steps && mag_cmp(size.value(), bound.value()) > 0)
                return false;
            acb_sub(center.value(), center.value(), correction.value(), bits);
            acb_get_mid(center.value(), center.value());
            if (!acb_is_finite(center.value()) || acb_is_zero(center.value()))
                return false;
            mag_set(previous_correction.value(), size.value());
            ++steps;
            if (counts)
                ++counts->machine_cover_newton_steps;
            return true;
        };
        // The certificate at the Newton point, from the evaluation in hand: the linear model at the centre,
        // with its remainder on a disk that holds the step and the goal's rectangle about the step's end
        // (certify_root_model takes the off-centre point as it takes piecewise's candidates). It costs the
        // remainder's majorant (the coefficients once, in magnitudes) where an evaluation at the stepped
        // point costs the blocks again. It is refused where that majorant on the step's disk exceeds what
        // the rectangle allows (a point more than one step from its rectangle, by the majorant's measure):
        // the point is then stepped and evaluated again. A point of a converged binary64 packet is one
        // step from its rectangle, so without this certificate every point would be evaluated twice.
        // History: docs/code-history-20261021.md, piecewise_roots.cpp 4.
        const auto certify_at_step = [&] {
            if (!steppable)
                return false;
            RootModel stepped;
            stepped.center = center;
            // The disk: the correction's size and two half-widths (the rectangle's corners lie sqrt(2) of
            // its half-width from the step's end, and the half-width there is this one to the goal's own
            // accuracy), as an exact radius. The correction's part is taken with a relative slack of
            // 2^-(modulus_bound_bits / 2): the certificate bounds the moduli of the step's end and of the
            // rectangle from above at modulus_bound_bits, and a disk that held the step exactly would be
            // refused for that rounding wherever the half-width is below 2^-modulus_bound_bits of the
            // correction (at a long goal). Half the bound's bits leave the rounding 2^-31 of the slack, and
            // the slack adds 2^-31 of the remainder.
            // History: docs/code-history-20261021.md, piecewise_roots.cpp 5.
            RealBall reach, slack;
            Float bound;
            upper(reach.value(), correction.value());
            arb_mul_2exp_si(slack.value(), reach.value(), -(modulus_bound_bits / 2));
            arb_add(reach.value(), reach.value(), slack.value(), bits);
            arb_addmul_si(reach.value(), half_width.value(), 2, bits);
            arb_get_ubound_arf(bound.value(), reach.value(), bits);
            arb_set_arf(stepped.radius.value(), bound.value());
            if (!arb_is_positive(stepped.radius.value()) || !arb_is_finite(stepped.radius.value()))
                return false;
            evaluations.linear_remainder(stepped.error.value(), stepped.center.value(), stepped.radius.value());
            ComplexBall scaled_slope, local;
            acb_poly_set_coeff_acb(stepped.polynomial.value(), 0, value.value());
            acb_mul_arb(scaled_slope.value(), derivative.value(), stepped.radius.value(), bits);
            acb_poly_set_coeff_acb(stepped.polynomial.value(), 1, scaled_slope.value());
            acb_div_arb(local.value(), correction.value(), stepped.radius.value(), bits);
            acb_neg(local.value(), local.value());
            Root proved;
            if (!certify_root_model(proved, stepped, local.value(), options, bits))
                return false;
            proved.certificate = "original-coefficient second-derivative Taylor remainder";
            root = std::move(proved);
            if (counts)
                ++counts->machine_cover_stepped_certificates;
            return true;
        };
        if (steppable) {
            // The rectangle's Rouche test needs |p(c)| below |p'(c)| times
            // the half-width, whatever the remainder is. A point that
            // misses that is certified at its Newton point, or stepped,
            // before the remainder about the point itself is computed.
            RealBall constant, slope;
            Float lower;
            upper(constant.value(), value.value());
            acb_get_abs_lbound_arf(lower.value(), derivative.value(), bits);
            arb_set_arf(slope.value(), lower.value());
            arb_mul(slope.value(), slope.value(), half_width.value(), bits);
            if (arb_ge(constant.value(), slope.value())) {
                if (certify_at_step())
                    return PointProof::Singleton;
                if (!take_step())
                    return PointProof::Refused;
                continue;
            }
        }
        // The model's disk is the smallest dyadic one around the rectangle
        // the certificate tries: its corners lie sqrt(2) half-widths from
        // the centre, so the radius is two half-widths. The remainder
        // majorant grows with the square of this radius.
        arb_mul_2exp_si(model.radius.value(), half_width.value(), 1);
        evaluations.linear_remainder(model.error.value(), model.center.value(), model.radius.value());
        acb_poly_set_coeff_acb(model.polynomial.value(), 0, value.value());
        acb_mul_arb(linear.value(), derivative.value(), model.radius.value(), bits);
        acb_poly_set_coeff_acb(model.polynomial.value(), 1, linear.value());
        ChartFailure failure = ChartFailure::Candidate;
        if (certify_root_model(root, model, zero.value(), options, bits, &failure)) {
            root.certificate = "original-coefficient second-derivative Taylor remainder";
            return PointProof::Singleton;
        }
        if (failure == ChartFailure::Candidate && steppable) {
            // Not accurate enough for the rectangle after all (the test's
            // margins): the Newton point, or a step, as above.
            if (certify_at_step())
                return PointProof::Singleton;
            if (!take_step())
                return PointProof::Refused;
            continue;
        }
        // Rounding: the certificate says so, or the value is too
        // inaccurate to step from.
        bool rounding = failure == ChartFailure::Rounding || failure == ChartFailure::Candidate;
        if (!rounding) {
            // The linear model cannot carry the rectangle (its remainder
            // majorant is a coefficient bound, far above |p''| where the
            // polynomial cancels). A disk needs no model: p'(c)/p(c) is the
            // sum of 1/(c - zeta) over the n roots, so some root lies
            // within n |p(c)/p'(c)| of c. The disk holds at least one
            // root; the attempt's closing test (n pairwise disjoint
            // regions of a polynomial of degree n) makes that exactly one.
            if (acb_contains_zero(derivative.value()))
                return PointProof::Refused;
            RealBall radius;
            Float lower;
            upper(radius.value(), value.value());
            acb_get_abs_lbound_arf(lower.value(), derivative.value(), bits);
            arb_div_arf(radius.value(), radius.value(), lower.value(), bits);
            arb_mul_si(radius.value(), radius.value(), degree, bits);
            Root disk;
            disk.box = center;
            arb_add_error(acb_realref(disk.box.value()), radius.value());
            arb_add_error(acb_imagref(disk.box.value()), radius.value());
            disk.count = disk.multiplicity = 1;
            if (acb_is_finite(disk.box.value())) {
                disk.publish(bits);
                if (meets_goal(disk, options, bits)) {
                    disk.certificate = "original-coefficient Newton inclusion disk in a complete disjoint cover";
                    root = std::move(disk);
                    if (counts)
                        ++counts->machine_cover_disks;
                    return PointProof::RootInclusion;
                }
            }
            // Too wide for the goal. A more accurate point narrows it;
            // when the value is already its own rounding, precision does.
            if (steppable) {
                if (!take_step())
                    return PointProof::Refused;
                continue;
            }
            rounding = true;
        }
        if (!raise())
            return PointProof::Refused;
    }
}
} // namespace

bool machine_newton_cover(std::vector<Root> &out, PolynomialSource &source, acb_srcptr initial,
                          const cps_options &options, slong precision, bool allow_underflow,
                          certpolysolve::detail::CandidateStatistics &statistics, bool coefficient_bound, const EnginePolicy *diagnostic,
                          std::unique_ptr<RootVector> *retained_proposals, bool *machine_converged,
                          slong *certificate_bits, bool scaled_packet, PartialNewtonCover *partial,
                          std::span<const uint8_t> settled, uint64_t fallback_operations) {
    if (machine_converged) *machine_converged = false;
    if (partial) *partial = {};
    const bool retry = !settled.empty();
    if (retry && (!partial || !coefficient_bound || !initial ||
                  settled.size() != static_cast<size_t>(source.degree())))
        return false;
    // A complete cover from the attempt is n pairwise disjoint regions of a
    // degree-n polynomial, each holding one root: every root isolated, as the
    // accuracy goal asks, and the callers give each region its factor's
    // multiplicity from the exact preparation. Either goal is served.
    if (options.real_only ||
        (engine_goal(options) != CPS_ROOT_APPROXIMATION && engine_goal(options) != CPS_ROOT_ACCURACY) || options.cancelled() ||
        precision < 2 || static_cast<uint64_t>(precision) > options.max_work_bits)
        return false;
    TraceScope shortcut(diagnostic,"machine_newton_cover",resources::nonnegative(precision),
                        resources::nonnegative(source.degree()),&statistics);
    auto *counts = diagnostic ? diagnostic->study_statistics.get() : nullptr;
    if (counts)
        counts->machine_cover_peak_bits = std::max(counts->machine_cover_peak_bits, static_cast<uint64_t>(precision));
    const auto coefficients = source.coefficient_enclosure(precision);
    const slong n = source.degree();
    if (n < 1 || acb_poly_degree(coefficients->polynomial.value()) != n)
        return false;
    RootVector points(n);
    if (initial)
        _acb_vec_set(points.data(), initial, n);
    else {
        // A sparse polynomial (at most a quarter of its coefficients nonzero, as candidate_jets.hpp decides for its
        // evaluation) starts from its Newton-polygon edges' binomial roots (initial_candidates).
        slong nonzero = 0;
        for (slong j = 0; j <= n; ++j)
            nonzero += !acb_is_zero(coefficients->polynomial.value()->coeffs + j);
        if (!initial_candidates(points.data(), coefficients->polynomial.value(), precision, 4 * nonzero <= n))
            return false;
    }
    // The attempt's prices, in the selector's complex coefficient operations (auto_selection.cpp): one
    // certificate evaluation (value and slope) and one sweep at one point of a repair round (costs), what the
    // attempt as a whole may spend (fallback: below, at the budget), and one quantum of its binary64 packet
    // at the packet's bound (heuristics::machine_cover_sweeps sweeps).
    std::vector<slong> support;
    if (source.sparse_available(precision))
        for (slong j = 0; j <= n; ++j)
            if (!acb_is_zero(coefficients->polynomial.value()->coeffs + j))
                support.push_back(j);
    const auto costs = newton_cover_costs(n, support);
    // The prices and the budget as exact counts, in 128 bits (at degree 2^20 the budget passes 2^48 and a quantum
    // 2^46; their sums stay far below 2^128): the packet's quantum at its bound is its operations at a point times the
    // degree (heuristics::machine_packet_point_operations, machine_cover_sweeps sweeps).
    const auto degree = static_cast<uint64_t>(n);
    const uint64_t terms = support.empty() ? degree + 1 : support.size();
    const uint128 fallback = fallback_operations > 0
        ? uint128{fallback_operations}
        : uint128{heuristics::initial_sweep_allowance} * degree * (degree + terms);
    const uint128 quantum_operations = uint128{heuristics::machine_packet_point_operations(degree, terms)} * degree;
    uint64_t quanta = 0;
    if (!retry) {
        // NumericalEngine::advance records candidate_ns itself. Binary64
        // first; a polynomial outside its range (no engine, or an arithmetic
        // failure in the packet) gets the same packet with a separate
        // exponent (ScaledBinary64), as the classical ladder's machine tier
        // does. Each packet is within the call's iteration budget; with none
        // left there is no attempt.
        bool advanced = false;
        for (const auto tier : {ArithmeticTier::Binary64, ArithmeticTier::ScaledBinary64}) {
            if (tier == ArithmeticTier::ScaledBinary64 && !scaled_packet)
                break;
            const unsigned sweeps = options.packet_sweeps(heuristics::machine_cover_sweeps, statistics);
            if (!sweeps)
                return false;
            auto engine = prepare_machine_engine(IterationMethod::Aberth, tier, coefficients->polynomial.value(),
                                                 points.view(), false, false, allow_underflow,
                                                 options.temporary_bytes, tier == ArithmeticTier::Binary64,
                                                 tier == ArithmeticTier::Binary64);
            if (!engine)
                continue;
            auto step = engine->advance(sweeps, options, statistics);
            if (step.candidates.cancelled)
                return false;
            if (step.arithmetic_failure)
                continue;
            ++quanta;
            // A quantum (heuristics::machine_cover_sweeps) schedules the packet; it does not end it. While a
            // quantum settles points, and the budget pays for another at its bound beside one certificate pass
            // (the certificate evaluations at every point), the packet goes on: a point that is still moving
            // costs a repair round below thousands of times more (n multiprecision interactions a sweep,
            // against the tree's few in binary64).
            // History: docs/code-history-20261021.md, piecewise_roots.cpp 2.
            const uint128 pass_operations = uint128{heuristics::certificate_evaluations} * degree * costs.evaluation;
            for (uint64_t settling = static_cast<uint64_t>(n);
                 !step.exhausted && !step.arithmetic_limited && !step.candidates.converged && step.active < settling &&
                 uint128{quanta + 1} * quantum_operations + pass_operations <= fallback;) {
                settling = step.active;
                const unsigned more = options.packet_sweeps(heuristics::machine_cover_sweeps, statistics);
                if (!more)
                    break;
                const auto next = engine->advance(more, options, statistics);
                ++quanta;
                if (next.candidates.cancelled)
                    return false;
                if (next.arithmetic_failure)
                    break; // the engine keeps the points this quantum began from
                step = next;
            }
            engine->export_candidates(points.data());
            // Converged: corrections reached the arithmetic's floor or became
            // small, rather than the sweep allowance running out.
            if (machine_converged)
                *machine_converged = !step.exhausted && (step.arithmetic_limited || step.candidates.converged);
            if (counts && tier == ArithmeticTier::ScaledBinary64)
                ++counts->scaled_packets;
            advanced = true;
            break;
        }
        if (!advanced)
            return false;
    }
    if (retained_proposals && !retry) {
        bool finite=true;
        for (slong i=0;i<n;++i) finite=finite && acb_is_finite(points[i]);
        if (finite) {
            auto proposed=std::make_unique<RootVector>(n);
            _acb_vec_set(proposed->data(),points.data(),n);
            *retained_proposals=std::move(proposed);
        }
    }
    std::optional<StageTimer> timer(std::in_place, statistics.validation_ns);
    NewtonAttemptProgress attempt_progress{precision};
    NewtonCoverEvaluation evaluations(source, n);
    if (coefficient_bound && !evaluations.prepare(attempt_progress.precision))
        return false;
    // The certificate's first tier (compensated_certificate.hpp), prepared once from the source's coefficients at the
    // two-word precision and a binary64 precision beyond it, so that the enclosure's own radius lies below the
    // words' residual. Where the attempt evaluates sparsely, a dense evaluation would cost the degree where the
    // sparse one costs the terms, and the ball certificate stays alone.
    std::optional<CompensatedCertificate> compensated;
    if (evaluations.available() && !evaluations.sparse() && (!diagnostic || diagnostic->compensated_certificates)) {
        const auto enclosed = source.coefficient_enclosure(machine::double_double_bits + machine::binary64_bits);
        compensated.emplace(enclosed->polynomial.value());
        if (!compensated->available())
            compensated.reset();
    }
    // What a point may still spend on precision raises once other points hold
    // certificates (set by the attempt below).
    double raise_allowance = std::numeric_limits<double>::infinity();
    // The precision at which the attempt's first pass has already asked the compensated tier about every point it
    // passes on, where they are (0 outside that pass): certify does not ask it again there (an evaluation of the
    // whole polynomial per point, 5% of a degree-3200 solve, 27 October). A raised precision asks again.
    slong compensated_refused_precision = 0;
    const auto certify = [&](Root &root, const acb_t point) -> PointProof {
        if (!evaluations.available())
            return certify_newton_point(root, source, point, options, precision)
                       ? PointProof::Singleton : PointProof::Refused;
        if (!acb_is_finite(point) || acb_is_zero(point) || options.cancelled())
            return PointProof::Refused;
        if (compensated && attempt_progress.precision != compensated_refused_precision &&
            compensated->certify(root, point, options, attempt_progress.precision)) {
            if (counts)
                ++counts->machine_cover_compensated_certificates;
            return PointProof::CompensatedSingleton;
        }
        return attempt_ball_point(root, point, evaluations, options, n, raise_allowance, attempt_progress, counts);
    };
    std::vector<Root> checked;
    for (slong i = 0; !evaluations.available() && i < n && !options.cancelled(); ++i) {
        ComplexBall point;
        acb_get_mid(point.value(), points[i]);
        // Without the coefficient bound the older certificate polishes first,
        // at the attempt's own precision, and the first refusal ends it.
        for (unsigned step = 0; step < 4 && !options.cancelled(); ++step) {
            ComplexBall value, derivative, correction;
            if (!source.evaluate(value.value(), derivative.value(), point.value(), precision) ||
                acb_contains_zero(derivative.value()))
                break;
            acb_div(correction.value(), value.value(), derivative.value(), precision);
            acb_sub(point.value(), point.value(), correction.value(), precision);
            if (!acb_is_finite(point.value()))
                break;
            acb_get_mid(point.value(), point.value());
            if (acb_is_zero(correction.value()) || acb_rel_accuracy_bits(correction.value()) < 4)
                break;
        }
        Root root;
        if (certify(root, point.value()) == PointProof::Refused)
            return false;
        checked.push_back(std::move(root));
    }
    // Whether the bounds pass below already found `checked` disjoint, so that
    // the final test need not compare its regions again.
    bool disjoint_known = false;
    // The highest precision a cluster start certified a group at (the certificate's precision returns to the attempt's
    // after it): the attempt reports the larger.
    slong start_bits = 0;
    if (evaluations.available()) {
        // Each accepted region contains at least one root independently of
        // the other points. Exact counts wait for the complete disjoint cover,
        // so a refusal need not end the
        // attempt: the refused points are iterated further, with the certified
        // ones held fixed, and tried again.
        checked.resize(static_cast<size_t>(n));
        std::vector<uint8_t> certified(static_cast<size_t>(n), 0);
        std::vector<slong> refused;
        // The certificate's precision when each refused point was refused in the pass (the re-ask below).
        std::vector<slong> refused_at(static_cast<size_t>(n), 0);
        // A refusal returns the certificates held (PartialNewtonCover): each
        // is a proof about the source whatever became of the other points.
        // Two certified regions that meet could hold the same root, and both
        // are given up; the rest are pairwise disjoint.
        bool pass_complete = false;
        const auto refuse = [&] {
            if (!partial || !std::count(certified.begin(), certified.end(), uint8_t(1)))
                return false;
            std::vector<size_t> index;
            std::vector<Root> held_regions;
            for (size_t i = 0; i < certified.size(); ++i)
                if (certified[i]) {
                    index.push_back(i);
                    held_regions.push_back(checked[i]);
                }
            const auto conflict = overlapping_roots(held_regions);
            for (size_t k = 0; k < index.size(); ++k)
                if (conflict[k])
                    certified[index[k]] = 0;
            partial->points = std::make_unique<RootVector>(n);
            _acb_vec_set(partial->points->data(), points.data(), n);
            partial->roots = std::move(checked);
            partial->certified = std::move(certified);
            partial->certificate_bits = std::max(attempt_progress.precision, start_bits);
            partial->pass_complete = pass_complete;
            return false;
        };
        // Refusals may cost one pass over the points at the attempt's first
        // precision and no more: that is what a single validation costs on
        // the path a failed attempt leaves the solve to, so an attempt that
        // cannot succeed adds at most one such pass.
        const double allowance = static_cast<double>(n) * operation_weight(precision);
        double refusals = 0;
        // The attempt as a whole spends at most what the caller's route costs
        // if it fails (fallback_operations, the selector's operations), so a
        // failed attempt at most doubles the solve. It is counted as the
        // certificates are, in evaluations weighted by operation_weight of
        // their precision; the selector's costs.evaluation operations are one
        // evaluation at the attempt's first precision (its weights were
        // measured there), which that count weighs operation_weight(precision).
        // The binary64 packet's quanta, each at its bound, come off the
        // budget first (none in a retry). What the rounds' packets below have
        // cost, in the same unit: packets.
        const uint128 packet_operations = uint128{quanta} * quantum_operations;
        // What a failed attempt leaves to its fallback, priced where that
        // fallback works. In a retry after the staged solve (settled) it is
        // the secular refinement of the refused points at the goal's
        // precision, a ladder that rebuilds the model of the whole
        // polynomial, at the precision the certificates run at: priced at
        // the attempt's first precision, half of it at a high goal, the
        // certificate pass alone would spend the budget and the round that
        // completes the cover would never run. A first attempt's fallback is
        // its caller's route from the attempt's first precision (the staged
        // solve, or piecewise's estimate): priced at the certificate
        // precision, a first attempt at a high goal would buy rounds at that
        // precision that the staged solve repairs at a far lower one.
        // History: docs/code-history-20261021.md, piecewise_roots.cpp 6.
        // A ratio, in binary64: the operations the packets left, exact (negative where they took more than the
        // budget), over an evaluation's.
        const double fallback_evaluations =
            static_cast<double>(static_cast<int128>(fallback) - static_cast<int128>(packet_operations)) /
            static_cast<double>(costs.evaluation);
        const auto budget = [&] {
            return fallback_evaluations * operation_weight(retry ? attempt_progress.precision : precision);
        };
        double packets = 0;
        bool exhausted = false;
        // The points that hold a ball certificate. A compensated certificate is proved in binary64 and says nothing
        // about the precision the ball certificate needs, so it is not counted: since the double-double repair
        // (proposal 1) proves every point the tier can first, counting those left the ball certificates of the
        // points after them no precision raise (rounded double roots at degree 64 stayed at 85 bits where they need
        // 170, 27 October 2026).
        size_t held = 0;
        // Whether the certificate is one of the cluster start's (below): paid from the budget, as the start's packet is,
        // it makes no raise, and a refusal spends none of the refusals' allowance and sets no source limit for the rest.
        bool starting = false;
        const auto attempt = [&](slong i) {
            const double before = attempt_progress.spent;
            const slong entry_precision = attempt_progress.precision;
            // While no point holds a ball certificate the precision is still
            // being found, as before; after that a point that asks for more
            // than sufficed for the others spends the refusals' allowance.
            raise_allowance = starting ? 0.0
                              : held   ? allowance - refusals
                                       : std::numeric_limits<double>::infinity();
            ComplexBall point;
            acb_get_mid(point.value(), points[i]);
            Root root;
            const PointProof proof = certify(root, point.value());
            if (proof == PointProof::Refused) {
                if (!starting)
                    refusals += attempt_progress.spent - before;
                // A precision raise is kept only when it produced a
                // certificate: a refused point asked for it in vain, and the
                // other points must not pay for it.
                if (attempt_progress.precision != entry_precision) {
                    attempt_progress.precision = entry_precision;
                    if (!evaluations.prepare(attempt_progress.precision))
                        exhausted = true;
                }
                exhausted = exhausted || (!starting && attempt_progress.source_limited) || refusals > allowance;
                return false;
            }
            // The certified centre: later packets hold it fixed.
            acb_get_mid(points[i], root.box.value());
            checked[static_cast<size_t>(i)] = std::move(root);
            certified[static_cast<size_t>(i)] = 1;
            if (proof != PointProof::CompensatedSingleton)
                ++held;
            return true;
        };
        // The compensated tier alone at point i (no ball certificate, nothing spent on a refusal): the double-double
        // repair's first and second tries (below). Its evaluation at a point it refuses is kept for the cluster start.
        std::vector<CompensatedEvaluation> compensated_seen(compensated ? static_cast<size_t>(n) : 0);
        const auto attempt_compensated = [&](slong i) {
            ComplexBall point;
            acb_get_mid(point.value(), points[i]);
            Root root;
            if (!acb_is_finite(point.value()) || acb_is_zero(point.value()) ||
                !compensated->certify(root, point.value(), options, attempt_progress.precision,
                                      &compensated_seen[static_cast<size_t>(i)]))
                return false;
            if (counts)
                ++counts->machine_cover_compensated_certificates;
            acb_get_mid(points[i], root.box.value());
            checked[static_cast<size_t>(i)] = std::move(root);
            certified[static_cast<size_t>(i)] = 1;
            return true;
        };
        // A round: a multiprecision Aberth packet on the refused points, at
        // the certificate's precision and within the call's iteration
        // budget, every other point held where it is; then the certificate
        // again. A round on two or more points must certify at least half
        // of them, or the attempt ends: the rounds on a set then cost at most
        // twice the first one, and all of them at most the attempt's budget.
        //
        // One point has no half. The packet is a quantum (heuristics.hpp): a
        // point may need two, so a round that leaves its one point refused
        // keeps it. What a point that stays can cost is bounded as every
        // repair is, by the balance below and by the refusals' allowance;
        // when it is the last, every other point is certified and held, the
        // polynomial deflated by them is linear, and a packet either brings
        // the point to its root or leaves a refusal that iterating cannot
        // change.
        //
        // A round takes the refused points, in order, that the attempt's
        // budget (above) still pays for; the others wait for the next round.
        // One sweep at one point costs costs.round_sweep operations: the
        // dense enclosure's value and slope and n interactions, which for a
        // sparse source is thousands of its certificate evaluations at high
        // degree.
        // History: docs/code-history-20261021.md, piecewise_roots.cpp 7.
        const double sweep_evaluations = static_cast<double>(costs.round_sweep) / static_cast<double>(costs.evaluation);
        CandidateWorkspace workspace(diagnostic);
        std::vector<uint8_t> fixed(static_cast<size_t>(n), 1);
        const auto round = [&] {
            const unsigned sweeps = options.packet_sweeps(heuristics::newton_cover_round_sweeps(options.bits), statistics);
            if (!sweeps)
                return false;
            // One refused point's share of the packet.
            const double share = static_cast<double>(sweeps) * sweep_evaluations * operation_weight(attempt_progress.precision);
            const size_t taken = std::min(
                refused.size(), static_cast<size_t>(std::max(0.0, std::floor((budget() - attempt_progress.spent - packets) / share))));
            if (!taken)
                return false;
            packets += static_cast<double>(taken) * share;
            const auto enclosed = source.coefficient_enclosure(attempt_progress.precision);
            if (acb_poly_degree(enclosed->polynomial.value()) != n)
                return false;
            const std::span<const slong> moving(refused.data(), taken);
            for (const slong i : moving)
                fixed[static_cast<size_t>(i)] = 0;
            RootVector moved(n);
            timer.reset(); // the packet keeps its own time (candidate_ns)
            const auto progress =
                candidate_packet(moved.data(), enclosed->polynomial.value(), points.data(), attempt_progress.precision,
                                 sweeps, CandidateEngine::Aberth, options, statistics, fixed, false, &workspace);
            timer.emplace(statistics.validation_ns);
            for (const slong i : moving)
                fixed[static_cast<size_t>(i)] = 1;
            if (progress.cancelled || !progress.finite)
                return false;
            if (counts)
                ++counts->machine_cover_rounds;
            std::vector<slong> left(refused.begin() + static_cast<std::ptrdiff_t>(taken), refused.end());
            size_t kept = 0;
            for (const slong i : moving) {
                if (options.cancelled())
                    return false;
                acb_get_mid(points[i], moved[i]);
                if (attempt(i))
                    continue;
                if (exhausted)
                    return false;
                left.push_back(i);
                ++kept;
            }
            if (taken > 1 && 2 * kept > taken)
                return false;
            refused = std::move(left);
            return true;
        };
        // The cluster start (research/d8-clusters): points the tiers refused whose Newton inclusion disks (n |p/p'| about each)
        // meet form a group of k; Newton's iteration on p^(k-1) from their centroid finds the cluster's centre (a
        // simple root of p^(k-1)), the ratio of the Taylor coefficients about it its radius rho = (|q_0| / |q_k|)^(1/k),
        // and an Aberth packet from k points on that circle, at a precision that resolves it (k log2(|c| / rho) above
        // the goal), brings them to the roots in 1 to 4 sweeps where Aberth from binary64's points needs 6 to 64 or
        // more (the probe's 27 degree-8 clusters, research/d8-clusters/README.md). Within the attempt's budget, as
        // the rounds; a group not certified goes back to where it was.
        // 0 when it ran (the points it certified leave the pool), 1 cancelled, 2 the refusals' allowance spent.
        // Points of groups the start resolved whose certified regions could not be disjoint at the goal (deferred): the
        // pass and the rounds leave them, and the attempt ends after the pass, returning the other points' certificates
        // where its caller keeps them (PartialNewtonCover); the caller's route finishes the deferred points.
        std::vector<uint8_t> deferred(static_cast<size_t>(n), 0);
        bool deferred_any = false;
        const auto cluster_start = [&](std::vector<slong> &pool) -> int {
            if (compensated && pool.size() >= 2 && pool.size() <= heuristics::cluster_start_points) {
                const auto groups = cluster_start_groups(points.view(), pool, *compensated,
                                                         compensated_seen, options);
                if (groups.empty())
                    return 0;
                std::vector<uint8_t> done(static_cast<size_t>(n), 0);
                // What the groups the start did not certify spent: kept apart, so that the start's failures together stay
                // within one budget while the rest of the attempt decides with what they gave back.
                double start_spent = 0;
                // The deferred groups and their points before the start, should their deferral be undone (below).
                std::vector<std::pair<std::vector<slong>, std::vector<ComplexBall>>> deferrals;
                // One group's start: Committed when every member is certified. Its packet's last correction, the
                // largest relative move of a point in the last sweep (log2), is kept for the deferral below.
                enum class Start { Declined, Committed, Deferred, Exhausted, Cancelled };
                double start_correction_log2 = 0;
                const auto start_group = [&](const std::vector<slong> &group, const slong k) -> Start {
                    const auto circle = size_cluster_start_circle(
                        source, points.view(), group, options, attempt_progress.precision, attempt_progress.spent);
                    if (!circle)
                        return Start::Declined;
                    const auto &centre = circle->centre;
                    const slong rho_exponent = circle->radius_exponent;
                    const double rho_fraction = circle->radius_fraction;
                    const double spread = circle->spread;
                    const slong packet_bits = circle->packet_bits;
                    // The packet's sweeps: a round's, as many as the budget still pays for at the packet's precision.
                    const double sweep_cost = static_cast<double>(k) * sweep_evaluations * operation_weight(packet_bits);
                    const double affordable = std::floor((budget() - attempt_progress.spent - packets - start_spent) / sweep_cost);
                    const unsigned sweeps = static_cast<unsigned>(std::clamp(
                        affordable, 0.0,
                        static_cast<double>(options.packet_sweeps(heuristics::newton_cover_round_sweeps(options.bits), statistics))));
                    if (!sweeps)
                        return Start::Declined;
                    packets += static_cast<double>(sweeps) * sweep_cost;
                    // The seeds: k points on the circle of radius rho about the centre (rho as a power of two and a factor in
                    // [1, 2): a cluster can be tighter than binary64's range), turned by heuristics::cluster_start_turn.
                    for (slong j = 0; j < k; ++j) {
                        const double angle = 2 * std::numbers::pi * static_cast<double>(j) / static_cast<double>(k) +
                                             heuristics::cluster_start_turn;
                        ComplexBall offset;
                        acb_set_d_d(offset.value(), std::cos(angle) * std::exp2(rho_fraction),
                                    std::sin(angle) * std::exp2(rho_fraction));
                        acb_mul_2exp_si(offset.value(), offset.value(), rho_exponent);
                        acb_add(points[group[static_cast<size_t>(j)]], centre.value(), offset.value(), packet_bits);
                        acb_get_mid(points[group[static_cast<size_t>(j)]], points[group[static_cast<size_t>(j)]]);
                    }
                    const auto at_packet = source.coefficient_enclosure(packet_bits);
                    for (const slong i : group)
                        fixed[static_cast<size_t>(i)] = 0;
                    RootVector moved(n);
                    timer.reset();
                    const auto progress = candidate_packet(moved.data(), at_packet->polynomial.value(), points.data(), packet_bits,
                                                           sweeps, CandidateEngine::Aberth, options, statistics, fixed, false,
                                                           &workspace);
                    timer.emplace(statistics.validation_ns);
                    for (const slong i : group)
                        fixed[static_cast<size_t>(i)] = 1;
                    if (progress.cancelled)
                        return Start::Cancelled;
                    // Resolved: the last sweep moved each point by less than the cluster's size. Otherwise the points have
                    // not separated, and their certificates would be spent for nothing.
                    start_correction_log2 = progress.correction_log2;
                    bool all = progress.finite && progress.correction_log2 < -spread;
                    // The group's certificates at the packet's precision at least: the cluster needed it, and below it
                    // the certificate's own Newton steps on cheaper rungs bring the resolved points back to their noise.
                    if (all && packet_bits > attempt_progress.precision) {
                        attempt_progress.precision = packet_bits;
                        if (!evaluations.prepare(attempt_progress.precision))
                            all = false;
                    }
                    // The group's certificates (starting): paid from the budget, a refusal spends none of the refusals'
                    // allowance, and no raise is made (the packet's precision is the cluster's need; a member that asks
                    // for more was not resolved by it).
                    starting = true;
                    if (all)
                        for (const slong i : group) {
                            acb_get_mid(points[i], moved[i]);
                            if (attempt(i))
                                done[static_cast<size_t>(i)] = 1;
                            else
                                all = false;
                            if (exhausted)
                                break;
                        }
                    starting = false;
                    if (exhausted)
                        return Start::Exhausted;
                    if (!all)
                        return Start::Declined;
                    // Committed only as part of a cover. Regions that meet cannot be one: a cluster tighter than the goal's
                    // regions (a pair 2^-60 apart at degree 32, regions of half-width about 2^-58), or two points at one
                    // root (the fuzz's axis case 13594533328001059812, a pair 2^-153 apart; its route asked again for
                    // every refused cover until its work ran out). The group is resolved but not certifiable here.
                    for (size_t m = 0; m < group.size(); ++m) {
                        const auto i = static_cast<size_t>(group[m]);
                        for (size_t j = 0; j < static_cast<size_t>(n); ++j)
                            if (j != i && certified[j] && acb_overlaps(checked[i].box.value(), checked[j].box.value()))
                                return Start::Deferred;
                    }
                    return Start::Committed;
                };
                for (const auto &group : groups) {
                    const slong k = static_cast<slong>(group.size());
                    if (k < 2 || k > heuristics::local_component_degree || options.cancelled())
                        continue;
                    // The packet costs at least a sweep at the certificate's precision: where the budget cannot pay that,
                    // the sizing would be spent for nothing.
                    if (budget() - attempt_progress.spent - packets - start_spent <
                        static_cast<double>(k) * sweep_evaluations * operation_weight(attempt_progress.precision))
                        continue;
                    // A group the start does not certify whole leaves the attempt as it found it: its points, the budget's
                    // accounts and the source limit (the refusals' allowance it never spends). After every group the
                    // certificate's precision is the attempt's again (a certified group's is reported as the attempt's
                    // highest) and the ball certificates' count (held) is the pass's: the group's certificates are at the
                    // cluster's precision, not one the pass has found.
                    const double spent_before = attempt_progress.spent, packets_before = packets;
                    const size_t held_before = held;
                    const bool limited_before = attempt_progress.source_limited;
                    const slong entry_bits = attempt_progress.precision;
                    std::vector<ComplexBall> group_before(group.size());
                    for (size_t m = 0; m < group.size(); ++m)
                        acb_set(group_before[m].value(), points[group[m]]);
                    const Start outcome = start_group(group, k);
                    if (outcome == Start::Cancelled)
                        return 1;
                    if (outcome == Start::Exhausted)
                        return 2;
                    held = held_before;
                    attempt_progress.source_limited = limited_before;
                    if (outcome == Start::Committed) {
                        start_bits = std::max(start_bits, attempt_progress.precision);
                        if (counts)
                            ++counts->machine_cover_cluster_starts;
                    } else {
                        // A deferred group leaves the pass and the rounds, its points back where they were: those would
                        // refuse them one by one and spend the budget on rounds that cannot certify them.
                        if (outcome == Start::Deferred) {
                            for (const slong i : group)
                                deferred[static_cast<size_t>(i)] = 1;
                            deferred_any = true;
                            deferrals.emplace_back(group, std::vector<ComplexBall>(group.size()));
                            for (size_t m = 0; m < group.size(); ++m)
                                acb_set(deferrals.back().second[m].value(), group_before[m].value());
                        }
                        start_spent += (attempt_progress.spent - spent_before) + (packets - packets_before);
                        attempt_progress.spent = spent_before;
                        packets = packets_before;
                        // A deferred group's resolved points stay where the caller keeps them as proposals (its route
                        // certifies them itself: perturbed-m2-k60-d032, the secular route then finishes the pair in 3 work
                        // units where it takes 8 from binary64's points), when every two are farther apart than their
                        // last moves together (twice the packet's last correction): closer, they are not known to sit at
                        // different roots, and two points at one root sent the route asking again until its work ran out
                        // (the fuzz's axis case 13594533328001059812). Its certificates go.
                        bool proposals_kept = outcome == Start::Deferred && partial;
                        for (size_t a = 0; proposals_kept && a < group.size(); ++a)
                            for (size_t b = a + 1; proposals_kept && b < group.size(); ++b) {
                                ComplexBall difference;
                                acb_sub(difference.value(), points[group[a]], points[group[b]], attempt_progress.precision);
                                Magnitude apart, size;
                                acb_get_mag_lower(apart.value(), difference.value());
                                acb_get_mag(size.value(), points[group[a]]);
                                slong ea = 0, es = 0;
                                double fa = 0, fs = 0;
                                proposals_kept = magnitude_log2(apart.value(), ea, fa) && magnitude_log2(size.value(), es, fs) &&
                                                 static_cast<double>(ea - es) + (fa - fs) > start_correction_log2 + 1;
                            }
                        for (size_t m = 0; m < group.size(); ++m) {
                            const auto i = static_cast<size_t>(group[m]);
                            if (done[i]) {
                                done[i] = 0;
                                certified[i] = 0;
                                checked[i] = Root{};
                            }
                            if (!proposals_kept)
                                acb_set(points[group[m]], group_before[m].value());
                        }
                    }
                    if (attempt_progress.precision != entry_bits) {
                        attempt_progress.precision = entry_bits;
                        if (!evaluations.prepare(attempt_progress.precision))
                            return 2;
                    }
                }
                // A deferral ends the attempt after its pass, which fits a pool of deferred clusters alone. With other
                // refused points left, the rounds take them together with the clusters' points (a round holds refused
                // neighbours together), and the deferred groups are declined: their points back where they were
                // (lsr1 at 768 bits of noise: a deferred triple beside two refused points took 33 work units where the
                // attempt's own path takes 14).
                if (deferred_any && std::any_of(pool.begin(), pool.end(), [&](slong i) {
                        return !done[static_cast<size_t>(i)] && !deferred[static_cast<size_t>(i)];
                    })) {
                    for (const auto &[group, saved] : deferrals)
                        for (size_t m = 0; m < group.size(); ++m) {
                            deferred[static_cast<size_t>(group[m])] = 0;
                            acb_set(points[group[m]], saved[m].value());
                        }
                    deferred_any = false;
                }
                std::vector<slong> still;
                for (const slong i : pool)
                    if (!done[static_cast<size_t>(i)] && !deferred[static_cast<size_t>(i)])
                        still.push_back(i);
                pool = std::move(still);
            }
            return 0;
        };
        // The rounds come after the pass over every point. A round holds
        // every point outside it where it is, and a point refused near an
        // ill-conditioned root has refused neighbours as far from their roots
        // as itself, so a round on some of them repairs against neighbours in
        // the wrong places; one round on all of them moves them together.
        // Rounds taken as the refusals come, in batches that double, would
        // end attempts that the rounds after the pass complete. An attempt
        // that cannot succeed costs its pass, within the refusals' allowance,
        // and its rounds, within the budget.
        // History: docs/code-history-20261021.md, piecewise_roots.cpp 8.
        std::vector<slong> order;
        for (slong i = 0; i < n; ++i)
            if (!(retry && settled[static_cast<size_t>(i)]))
                order.push_back(i);
        // Proposal 1 (dev/docs/cluster-double-double-design-20261027.md): the compensated tier alone on every point
        // first. A few points it refuses, the roots of small clusters that binary64 leaves at about u / s^(k-1), get
        // one repair in double-double (CompensatedCertificate::refine, every other point held where it is) and the
        // tier again, before any ball certificate; the points still refused go on below as before. The repair is a
        // proposal: every root it helps is proved by the compensated tier's own test.
        if (compensated && (!diagnostic || diagnostic->double_double_repair)) {
            std::vector<slong> waiting;
            for (const slong i : order) {
                if (options.cancelled())
                    return false;
                if (!attempt_compensated(i))
                    waiting.push_back(i);
            }
            // Paid from the attempt's budget, as a round is: a step at a point costs a round's sweep at a point at
            // double-double's weight, the steps are capped so that every point taking every one of them stays within
            // the repair's share of what is left, and the steps it takes are charged.
            const double point_step = sweep_evaluations * heuristics::double_double_operation_weight;
            const double affordable =
                waiting.empty() ? 0
                                : std::floor(heuristics::double_double_budget_share * (budget() - attempt_progress.spent - packets) /
                                             (static_cast<double>(waiting.size()) * point_step));
            const unsigned steps = static_cast<unsigned>(
                std::clamp(affordable, 0.0, static_cast<double>(heuristics::double_double_refine_steps)));
            if (steps && waiting.size() <= heuristics::double_double_refine_points) {
                // A point the repair does not let the tier prove goes on from where it was: the ball certificate and
                // its precision raises then see the points they see without the repair (a repaired point the tier
                // still refuses can need a higher precision that its new place does not reach: rounded double roots
                // of degree 64, 27 October).
                std::vector<ComplexBall> before(waiting.size());
                for (size_t m = 0; m < waiting.size(); ++m)
                    acb_set(before[m].value(), points[waiting[m]]);
                uint64_t evaluated = 0;
                const bool moved =
                    compensated->refine(points.data(), n, waiting, static_cast<int64_t>(options.bits), evaluated, steps);
                packets += static_cast<double>(evaluated) * point_step;
                if (moved) {
                    if (counts)
                        ++counts->machine_cover_double_double_repairs;
                    std::vector<slong> still;
                    for (size_t m = 0; m < waiting.size(); ++m) {
                        if (attempt_compensated(waiting[m]))
                            continue;
                        acb_set(points[waiting[m]], before[m].value());
                        still.push_back(waiting[m]);
                    }
                    waiting = std::move(still);
                }
            }
            order = std::move(waiting);
            // Before the ball certificate's pass: groups of the points the tiers refuse, clusters the double-double
            // repair could not reach, go to the cluster start, where the pass would refuse them one by one.
            // Not at a goal staged below the caller's (cps_options::staged_goal): a group certified there sent the caller's
            // goal down a path the stage's iteration does not (edge53 5276379880489863530 of the fuzz: the refinement then
            // stayed at 468 bits for its 100,000 units where main reaches 730 bits in 10).
            switch ((!diagnostic || diagnostic->cluster_start) && !options.staged_goal ? cluster_start(order) : 0) {
            case 1:
                return false;
            case 2:
                return refuse();
            default:
                break;
            }
            // Every point left was refused by the tier where it now is, at this precision.
            compensated_refused_precision = attempt_progress.precision;
        }
        for (const slong i : order) {
            if (options.cancelled())
                return false;
            if (attempt(i))
                continue;
            if (exhausted)
                return refuse();
            refused_at[static_cast<size_t>(i)] = attempt_progress.precision;
            refused.push_back(i);
        }
        // The rounds move the points: the tier is asked again.
        compensated_refused_precision = 0;
        pass_complete = true;
        // A point refused before a later point raised the attempt's precision (its value or slope without a correct
        // bit at the precision of its turn) is asked again at the raised one before any round: the points after it
        // were certified there directly, where a round spends a packet on each (Chebyshev 1,024: 238 points refused
        // at 85 bits, then rounds for 36 of them at 1,360 and none affordable after, the attempt failing; 28 October
        // 2026). Like the rounds, within the attempt's budget: an evaluation at the raised precision a point, while
        // what is left pays for one; the points left go to the rounds. Unbounded, it certified Chebyshev 8,192's 2,030
        // refused points at 6,656 bits one by one where the attempt had spent its budget 2.6 times over and the
        // secular route covers them for less (28 s to 51 s).
        {
            const double reask = operation_weight(attempt_progress.precision);
            std::vector<slong> still;
            for (const slong i : refused) {
                if (options.cancelled())
                    return false;
                if (refused_at[static_cast<size_t>(i)] < attempt_progress.precision &&
                    budget() - attempt_progress.spent - packets >= reask && attempt(i))
                    continue;
                if (exhausted)
                    return refuse();
                still.push_back(i);
            }
            refused = std::move(still);
        }
        while (!refused.empty()) {
            if (options.cancelled())
                return false;
            if (!round())
                return refuse();
        }
        // The deferred groups' points hold no certificate: the other points' are returned for the caller to keep.
        if (deferred_any)
            return refuse();
        if (retry || !disjoint(checked))
            return refuse();
        disjoint_known = true;
    }
    if (checked.size() != static_cast<size_t>(n) || (!disjoint_known && !disjoint(checked)))
        return false;
    out = std::move(checked);
    // The precision the certificates were computed at: the caller reports it.
    if (certificate_bits)
        *certificate_bits = std::max(attempt_progress.precision, start_bits);
    return true;
}
} // namespace certpolysolve::detail
