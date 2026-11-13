// The exact real engine (RealRefinement, real_refinement.hpp): it isolates the real roots of an exact square-free
// integer factor p by the Descartes method and refines each interval to the goal, one bounded action a call of step.
//
// Isolation (isolate_step) takes nodes from a stack (pending) that starts with (-B, B), B a power of two with every
// root strictly inside. A node holds an exact primitive integer polynomial q whose roots in (0, 1) are p's in its open
// interval (restrict_unit); variations bounds their number by Descartes' rule of signs on (1+s)^n q(1/(1+s)), a Taylor
// shift of q's reversal. No variation drops the node; one, with p nonzero at both ends (exact_sign), makes it a
// component. Otherwise, with k >= 2 variations at an even depth, contract tries the Newton step t - k q(t)/q'(t) from
// t = 1/2 and keeps the quarter of the interval about it only when Descartes' rule proves the discarded sides
// root-free and the new ends are not roots; failing that the node is bisected, the half with fewer variations first.
// Refinement (refine, on the unfinished components in turn) tries an interval Newton step with p' enclosed on the
// whole interval, then Abbott's quadratic interval refinement (refine_by_cells: the secant's zero placed among
// 2^grid_bits cells, the cell confirmed by signs), then bisection. Signs are certified (sign_at: balls at doubling
// precision, then exact_sign where none decides); an exact zero at a tested point becomes a rational component. A
// component at the precision cap waits (capped_at) while the others go on; if all meet the goal but overlap,
// extra_bits tightens the goal by 32 bits. cover returns the regions once isolation is done and they are disjoint.
//
// State: RealRefinement::Impl: the factor and its derivative, pending (IsolationNode), components (Component:
// interval, any exact rational root, precision, sign_precision, grid_bits, capped_at, and the published region and
// goal verdict, cached while they hold), the round-robin cursor, extra_bits and the cached goal (set_goal).
// Entry points: RealRefinement (step, cover, count, precision, cells), called from solver.cpp only (SolveState::Impl:
// factor_step's real-only branch keeps one per factor in FactorState::real_search; live_cells counts its cells) and
// from tests/test_refinement_state.cpp. A real-only Auto query on an exact real polynomial takes the complex route
// (real_only_route.hpp); this branch answers the real-only queries of native secular input and of D5 routing (tests).
// Not here: squarefree_factors (exact_preprocessing.cpp), the Sturm engine the portfolio policy can run instead
// (SturmRefinement, sturm.cpp), accuracy_met and disjoint (certification.cpp).
// Design notes: docs/mathematics.md, "Exact real isolation and refinement"; docs/real-roots-inventory-20261011.md.
#include "methods/real/real_refinement.hpp"
#include "core/telemetry.hpp"
#include <algorithm>
#include <bit>
#include <cmath>
#include <flint/arb_fmpz_poly.h>
#include <optional>
#include <stdexcept>

namespace certpolysolve::detail {
namespace {
struct Component {
    RealBall interval;
    std::optional<Rational> rational;
    slong precision;
    // The precision at which a ball last decided a sign near this root
    // (sign_at): where the next sign is tried first. It changes no answer.
    slong sign_precision;
    // The published region (Impl::published) and its goal verdict, kept
    // while the component stays as it was published: a refinement step
    // changes one component, and the goal test of every step reads all of
    // them. History: docs/code-history-20261021.md, real_refinement.cpp 1.
    mutable std::optional<Root> region;
    mutable RealBall region_interval;
    mutable slong region_precision = 0;
    mutable bool region_rational = false;
    mutable uint64_t verdict_generation = 0;
    mutable bool verdict = false;
    // log2 of the number of cells the next refinement step places the root
    // among (Impl::refine_by_cells); four to begin with, as in Abbott's QIR.
    ulong grid_bits = 2;
    // The cap at which this component's refinement last stopped for
    // precision (zero: never). It is passed over while the cap stays there,
    // so one root at the cap does not stop the refinement of the others.
    uint64_t capped_at = 0;
    explicit Component(slong p) : precision(p), sign_precision(p) {}
};
struct IsolationNode {
    RealBall interval;
    // A nonzero scalar multiple of the original p(a+(b-a)t).
    // Descartes counts below concern the OPEN unit interval; endpoint roots
    // are recorded separately and cannot be counted by two adjacent nodes.
    IntegerPolynomial polynomial;
    uint64_t depth = 0;
};
slong variations(const fmpz_poly_t polynomial) {
    IntegerPolynomial transformed;
    Integer one;
    fmpz_one(one.value());
    const slong length = fmpz_poly_length(polynomial);
    fmpz_poly_reverse(transformed.value(), polynomial, length);
    fmpz_poly_taylor_shift(transformed.value(), transformed.value(), one.value());
    // (t+1)^n p(1/(t+1)): positive t corresponds exactly to x in (0,1).
    slong changes = 0;
    int previous = 0;
    for (slong i = 0; i < fmpz_poly_length(transformed.value()); ++i) {
        const int sign = fmpz_sgn(transformed.value()->coeffs + i);
        if (sign) {
            changes += previous && sign != previous;
            previous = sign;
        }
    }
    return changes;
}
IntegerPolynomial restrict_unit(const fmpz_poly_t polynomial, ulong a, ulong b, ulong denominator) {
    // denominator^n p((a+(b-a)t)/denominator), with exact integer algebra.
    IntegerPolynomial scaled, affine, result;
    Integer power, coefficient;
    fmpz_one(power.value());
    for (slong j = fmpz_poly_degree(polynomial); j >= 0; --j) {
        fmpz_mul(coefficient.value(), polynomial->coeffs + j, power.value());
        fmpz_poly_set_coeff_fmpz(scaled.value(), j, coefficient.value());
        fmpz_mul_ui(power.value(), power.value(), denominator);
    }
    fmpz_poly_set_coeff_ui(affine.value(), 0, a);
    fmpz_poly_set_coeff_ui(affine.value(), 1, b - a);
    fmpz_poly_compose(result.value(), scaled.value(), affine.value());
    fmpz_poly_primitive_part(result.value(), result.value());
    return result;
}
// Isolation endpoints and bisections of an Arb ball have a radius that is
// exactly representable by mag. The generic interval constructor adds a guard
// ulp even in this case, which can move the endpoint outside its parent.
void exact_interval(arb_t out, const arf_t a, const arf_t b) {
    Float midpoint, radius;
    arf_add(midpoint.value(), a, b, ARF_PREC_EXACT, ARF_RND_NEAR);
    arf_mul_2exp_si(midpoint.value(), midpoint.value(), -1);
    arf_sub(radius.value(), b, a, ARF_PREC_EXACT, ARF_RND_NEAR);
    arf_mul_2exp_si(radius.value(), radius.value(), -1);
    arb_set_arf(out, midpoint.value());
    arf_get_mag_lower(arb_radref(out), radius.value());
    Float represented;
    arf_set_mag(represented.value(), arb_radref(out));
    if (!arf_equal(represented.value(), radius.value()))
        throw std::logic_error("bisection radius was not exactly representable");
}
int exact_sign(const fmpz_poly_t polynomial, const arf_t point) {
    Rational x, value;
    arf_get_fmpq(x.value(), point);
    fmpz_poly_evaluate_fmpq(value.value(), polynomial, x.value());
    return fmpq_sgn(value.value());
}
} // namespace
struct RealRefinement::Impl {
    IntegerPolynomial polynomial, derivative;
    std::vector<Component> components;
    std::vector<IsolationNode> pending;
    uint64_t multiplicity;
    certpolysolve::detail::SearchStatistics &statistics;
    certpolysolve::detail::CandidateStatistics &extended;
    size_t cursor = 0;
    uint64_t extra_bits = 0;
    slong initial_precision;
    Impl(const fmpz_poly_t source, uint64_t m, slong prec, certpolysolve::detail::SearchStatistics &s, certpolysolve::detail::CandidateStatistics &e)
        : multiplicity(m), statistics(s), extended(e), initial_precision(prec) {
        fmpz_poly_set(polynomial.value(), source);
        const slong degree = fmpz_poly_degree(source);
        if (degree > 0) {
            // At B=2^(h+1), every |a_j/a_n|/B^(n-j) < 2^(-(n-j));
            // their sum is <1. Thus all roots lie STRICTLY inside (-B,B).
            slong h = 0;
            const slong leading_bits = static_cast<slong>(fmpz_bits(source->coeffs + degree));
            for (slong j = 0; j < degree; ++j)
                if (!fmpz_is_zero(source->coeffs + j)) {
                    const slong bits = static_cast<slong>(fmpz_bits(source->coeffs + j)) - leading_bits + 1;
                    if (bits > 0)
                        h = std::max(h, (bits + degree - j - 1) / (degree - j));
                }
            Integer bound;
            fmpz_one(bound.value());
            fmpz_mul_2exp(bound.value(), bound.value(), static_cast<ulong>(h + 1));
            IntegerPolynomial affine;
            fmpz_poly_set_coeff_fmpz(affine.value(), 1, bound.value());
            fmpz_poly_scalar_mul_2exp(affine.value(), affine.value(), 1);
            fmpz_neg(bound.value(), bound.value());
            fmpz_poly_set_coeff_fmpz(affine.value(), 0, bound.value());
            IsolationNode first;
            fmpz_poly_compose(first.polynomial.value(), source, affine.value());
            fmpz_poly_primitive_part(first.polynomial.value(), first.polynomial.value());
            arb_zero(first.interval.value());
            mag_set_ui_2exp_si(arb_radref(first.interval.value()), 1, h + 1);
            pending.push_back(std::move(first));
        }
        fmpz_poly_derivative(derivative.value(), polynomial.value());
    }
    bool contract(IsolationNode &node, slong count, const cps_options &options) {
        if (count < 2 || node.depth % 2)
            return false;
        ++statistics.newton_attempts;
        // Modified Newton proposes a cluster center in normalized coordinates.
        // The proposal is accepted ONLY after exact Descartes exclusions prove
        // both discarded open sides empty and both new boundaries nonzero.
        IntegerPolynomial slope;
        fmpz_poly_derivative(slope.value(), node.polynomial.value());
        RealBall midpoint, value, derivative_value, next;
        arb_one(midpoint.value());
        arb_mul_2exp_si(midpoint.value(), midpoint.value(), -1);
        arb_fmpz_poly_evaluate_arb(value.value(), node.polynomial.value(), midpoint.value(), initial_precision);
        arb_fmpz_poly_evaluate_arb(derivative_value.value(), slope.value(), midpoint.value(), initial_precision);
        if (arb_contains_zero(derivative_value.value()))
            return false;
        arb_div(next.value(), value.value(), derivative_value.value(), initial_precision);
        arb_mul_si(next.value(), next.value(), count, initial_precision);
        arb_sub(next.value(), midpoint.value(), next.value(), initial_precision);
        const double proposal = arf_get_d(arb_midref(next.value()), ARF_RND_NEAR);
        if (!std::isfinite(proposal) || proposal < 0.125 || proposal > 0.875)
            return false;
        // Keep the retained width a power-of-two fraction of the parent.
        // Otherwise repeated contractions eventually exceed mag's exact
        // mantissa, and rounding a boundary could admit an excluded root.
        const ulong a = static_cast<ulong>(std::floor(16 * proposal)) - 2, b = a + 4;
        auto lower = restrict_unit(node.polynomial.value(), 0, a, 16);
        auto upper = restrict_unit(node.polynomial.value(), b, 16, 16);
        if (variations(lower.value()) || variations(upper.value()))
            return false;
        Rational boundary, at_boundary;
        for (ulong numerator : {a, b}) {
            fmpq_set_ui(boundary.value(), numerator, 16);
            fmpz_poly_evaluate_fmpq(at_boundary.value(), node.polynomial.value(), boundary.value());
            if (fmpq_is_zero(at_boundary.value()))
                return false;
        }
        if (options.cancelled())
            return false;
        auto inside = restrict_unit(node.polynomial.value(), a, b, 16);
        Float lo, hi, width, new_lo, new_hi;
        arb_get_interval_arf(lo.value(), hi.value(), node.interval.value(), ARF_PREC_EXACT);
        arf_sub(width.value(), hi.value(), lo.value(), ARF_PREC_EXACT, ARF_RND_NEAR);
        arf_mul_2exp_si(width.value(), width.value(), -4);
        arf_mul_ui(new_lo.value(), width.value(), a, ARF_PREC_EXACT, ARF_RND_NEAR);
        arf_mul_ui(new_hi.value(), width.value(), b, ARF_PREC_EXACT, ARF_RND_NEAR);
        arf_add(new_lo.value(), new_lo.value(), lo.value(), ARF_PREC_EXACT, ARF_RND_NEAR);
        arf_add(new_hi.value(), new_hi.value(), lo.value(), ARF_PREC_EXACT, ARF_RND_NEAR);
        exact_interval(node.interval.value(), new_lo.value(), new_hi.value());
        node.polynomial = std::move(inside);
        ++node.depth;
        ++statistics.newton_successes;
        return true;
    }
    cps_stop_reason isolate_step(const cps_options &options) {
        auto &node = pending.back();
        ++extended.real_refinement_steps;
        StageTimer timer(statistics.refinement_ns);
        const slong count = variations(node.polynomial.value());
        if (!count) {
            pending.pop_back();
            return CPS_STOP_NONE;
        }
        Float a, b;
        arb_get_interval_arf(a.value(), b.value(), node.interval.value(), ARF_PREC_EXACT);
        if (count == 1 && exact_sign(polynomial.value(), a.value()) && exact_sign(polynomial.value(), b.value())) {
            Component component(initial_precision);
            component.interval = node.interval;
            components.push_back(std::move(component));
            pending.pop_back();
            return CPS_STOP_NONE;
        }
        if (pending.size() + components.size() + 2 > options.live_cells)
            return CPS_STOP_CELL_BUDGET;
        const slong needed = resources::precision_sum(std::max(arf_bits(a.value()), arf_bits(b.value())), 5);
        if (needed > static_cast<slong>(options.max_work_bits))
            return CPS_STOP_PRECISION;
        if (contract(node, count, options))
            return CPS_STOP_NONE;
        if (options.cancelled())
            return CPS_STOP_CANCELLED;
        Float midpoint;
        arf_set(midpoint.value(), arb_midref(node.interval.value()));
        if (exact_sign(polynomial.value(), midpoint.value()) == 0) {
            Component exact(initial_precision);
            exact.rational.emplace();
            arf_get_fmpq(exact.rational->value(), midpoint.value());
            arb_set_arf(exact.interval.value(), midpoint.value());
            components.push_back(std::move(exact));
        }
        IsolationNode left, right;
        left.polynomial = restrict_unit(node.polynomial.value(), 0, 1, 2);
        right.polynomial = restrict_unit(node.polynomial.value(), 1, 2, 2);
        exact_interval(left.interval.value(), a.value(), midpoint.value());
        exact_interval(right.interval.value(), midpoint.value(), b.value());
        left.depth = right.depth = node.depth + 1;
        pending.pop_back();
        // The half with fewer sign variations is taken first, so a half that holds a cluster waits only for its
        // sibling's few steps. Taken left first, a pair that lies in the left halves would leave every right
        // half on the stack until the pair is split, each with its own exact polynomial whose coefficients grow
        // by the degree's bits a bisection. Descartes' count bounds what a half's isolation can cost, so the
        // stack holds one half a level of a cluster's descent. History:
        // docs/code-history-20261021.md, real_refinement.cpp 2.
        const bool left_first = variations(left.polynomial.value()) <= variations(right.polynomial.value());
        pending.push_back(std::move(left_first ? right : left));
        pending.push_back(std::move(left_first ? left : right));
        ++extended.real_bisections;
        return CPS_STOP_NONE;
    }
    const Root &published(const Component &component) const {
        if (component.region && component.region_precision == component.precision &&
            component.region_rational == component.rational.has_value() &&
            arb_equal(component.region_interval.value(), component.interval.value()))
            return *component.region;
        Root r;
        acb_set_arb(r.box.value(), component.interval.value());
        r.count = r.multiplicity = multiplicity;
        r.realness = CPS_REAL;
        r.certificate =
            component.rational ? "exact-real-rational-root" : "native-descartes-and-interval-refinement";
        r.publish(component.precision);
        component.region = std::move(r);
        arb_set(component.region_interval.value(), component.interval.value());
        component.region_precision = component.precision;
        component.region_rational = component.rational.has_value();
        component.verdict_generation = 0;
        return *component.region;
    }
    std::vector<Root> roots() const {
        std::vector<Root> result;
        for (const auto &c : components)
            result.push_back(published(c));
        return result;
    }
    // The goal every component is tested against: the requested bits plus
    // extra_bits and, without extra bits, the caller's absolute tolerance.
    // goal_generation changes whenever it does, which discards every kept
    // verdict.
    uint64_t goal_generation = 1, goal_bits = 0;
    Number goal_tolerance;
    void set_goal(const cps_options &options) {
        const uint64_t bits = std::min<uint64_t>(resources::precision_max, options.bits + extra_bits);
        const Number &tolerance = extra_bits ? Number() : options.absolute_tolerance;
        if (bits == goal_bits && goal_tolerance.radix == tolerance.radix &&
            fmpq_equal(goal_tolerance.significand.value(), tolerance.significand.value()) &&
            fmpz_equal(goal_tolerance.exponent.value(), tolerance.exponent.value()))
            return;
        goal_bits = bits;
        goal_tolerance = tolerance;
        ++goal_generation;
    }
    bool goal(const Component &c) const {
        // Once touching exports require separation, the user's loose absolute
        // tolerance must not make every component ready before it contracts.
        const Root &region = published(c);
        if (c.verdict_generation != goal_generation) {
            c.verdict = accuracy_met(region, goal_bits, goal_tolerance, c.precision);
            c.verdict_generation = goal_generation;
        }
        return c.verdict;
    }
    // Bits from the lowest to the highest bit the exact value of p can have at
    // a dyadic point x = M 2^s (M with b bits): every term c_j M^j 2^(js)
    // lies within bits(c) + n (b + |s|) bits of the lowest bit any term can
    // have, and n + 1 terms add bit_width(n + 1). A ball evaluation at this
    // precision rounds nothing, so it is the exact evaluation.
    slong exact_value_bits(const arf_t point) const {
        const slong bits = arf_bits(point);
        const slong scale = arf_abs_bound_lt_2exp_si(point) - bits;
        const auto degree = static_cast<uint64_t>(std::max<slong>(0, fmpz_poly_degree(polynomial.value())));
        const auto spread = static_cast<uint64_t>(bits) + static_cast<uint64_t>(scale < 0 ? -scale : scale);
        return resources::precision_sum(
            resources::precision_sum(std::abs(fmpz_poly_max_bits(polynomial.value())),
                                     static_cast<slong>(std::min<uint64_t>(resources::multiply(degree, spread),
                                                                           resources::precision_max))),
            std::bit_width(degree + 1));
    }
    // The sign of p at a dyadic point, as exact_sign gives it. A ball
    // evaluation decides it unless the ball contains zero; the precision then
    // doubles, up to the cap and up to the exact value's size (where the ball
    // is the exact evaluation), and the exact rational evaluation decides only
    // what no ball does: the value zero, or a sign the cap cannot resolve.
    // Polynomials with many real roots cancel about their degree in bits near
    // a root, so at a component's working precision a ball rarely decides, and
    // every refinement step would pay one or two exact evaluations. The balls
    // start at `from` or at `hint`, the precision of the last ball that
    // decided near this root, whichever is larger; `hint` is updated. History:
    // docs/code-history-20261021.md, real_refinement.cpp 3.
    int sign_at(const arf_t point, slong from, slong &hint, const cps_options &options) const {
        const slong cap = static_cast<slong>(options.max_work_bits), exact = exact_value_bits(point);
        RealBall at, value;
        arb_set_arf(at.value(), point);
        for (slong p = std::min(cap, std::max(from, hint)); p >= from;
             p = std::min(cap, resources::precision_sum(p, p))) {
            arb_fmpz_poly_evaluate_arb(value.value(), polynomial.value(), at.value(), p);
            if (arb_is_zero(value.value()))
                return 0;
            if (!arb_contains_zero(value.value())) {
                hint = p;
                return arb_is_positive(value.value()) ? 1 : -1;
            }
            if (p >= exact || p >= cap)
                break;
        }
        return exact_sign(polynomial.value(), point);
    }
    enum class Cells { None, Cell, Side, Root };
    // Quadratic interval refinement (J. Abbott, "Quadratic interval
    // refinement for real roots", 2006; as A. Kerber and M. Sagraloff use it
    // in "Efficient real root approximation", ISSAC 2011). [a, b] holds one
    // simple root, with p(a) of sign sign_a and p(b) of the other. The secant
    // through the ends places the root in one of N = 2^g equal cells, and the
    // signs at that cell's two ends confirm it: the interval shrinks N-fold.
    // Near a simple root the secant's error is quadratic in the width, so a
    // confirmed cell squares N for the next step and a refuted one halves g
    // (N = 2 is a bisection). The signs are exact (sign_at): a wrong guess
    // costs time, never an answer.
    //
    // Why: the interval Newton step encloses p' on the whole interval by
    // ball Horner, which overestimates by about the degree in bits where many
    // roots are real, so it certifies nothing until the interval is far below
    // the goal, and nearly every step would bisect. History:
    // docs/code-history-20261021.md, real_refinement.cpp 4.
    //
    // Cell: [a, b] is the confirmed cell; Side: the cell was refuted, and
    // [a, b] is the side of the tested point that holds the root; Root: the
    // component is an exact rational root; None: bisect instead (the new
    // ends would pass the cap, no ball short of the cap gave an end's value,
    // or a refuted cell left more than half the interval). sign_a is p's sign
    // at a when it was found (zero otherwise).
    Cells refine_by_cells(Component &c, arf_t a, arf_t b, int &sign_a, slong work, const cps_options &options) {
        const slong cap = static_cast<slong>(options.max_work_bits);
        const auto g = static_cast<slong>(c.grid_bits);
        // The new ends carry up to g more bits than the old (and two guard
        // bits, as the bisection's): beyond the cap, bisect instead.
        if (resources::precision_sum(std::max(arf_bits(a), arf_bits(b)), resources::precision_sum(g, 2)) > cap) {
            c.grid_bits = std::max<ulong>(1, c.grid_bits / 2);
            return Cells::None;
        }
        // lambda = p(a) / (p(a) - p(b)) puts the secant's zero at
        // a + lambda (b - a). Choosing among 2^g cells needs lambda to about
        // 2^-(g + 2), so values with that relative error: g + 2 bits beyond
        // the precision at which a ball last decided a sign here. p(a) and
        // p(b) have opposite signs, so their difference does not cancel.
        RealBall value_a, value_b, lambda;
        const auto exact_root = [&](const arf_t x) {
            c.rational.emplace();
            arf_get_fmpq(c.rational->value(), x);
            arb_set_arf(c.interval.value(), x);
            return Cells::Root;
        };
        // p at an end, to g + 2 bits beyond the precision at which a ball last
        // decided a sign here; where that ball still holds zero, the sign is
        // found first (sign_at raises that precision) and the value taken
        // again. Returns the sign (zero: an exact root); `known` is false when
        // no ball short of the cap decides it.
        bool known = true;
        const auto end_value = [&](const arf_t x, arb_t value) {
            RealBall at;
            arb_set_arf(at.value(), x);
            const auto evaluate = [&] {
                const slong precision = std::min(
                    cap, resources::precision_sum(std::max(work, c.sign_precision), resources::precision_sum(g, 2)));
                arb_fmpz_poly_evaluate_arb(value, polynomial.value(), at.value(), precision);
                return !arb_contains_zero(value);
            };
            if (evaluate())
                return arb_is_positive(value) ? 1 : -1;
            const int sign = sign_at(x, work, c.sign_precision, options);
            known = known && sign && evaluate();
            return sign;
        };
        sign_a = end_value(a, value_a.value());
        if (!sign_a)
            return exact_root(a);
        const int sign_b = end_value(b, value_b.value());
        if (!sign_b)
            return exact_root(b);
        if (sign_a == sign_b)
            throw std::logic_error("a refined real root's interval lost its sign change");
        if (!known)
            return Cells::None;
        // lambda to g + 2 bits: its two roundings stay below the values' own
        // errors at that many bits beyond the working precision.
        const slong precision = resources::precision_sum(work, resources::precision_sum(g, 2));
        arb_sub(lambda.value(), value_a.value(), value_b.value(), precision);
        arb_div(lambda.value(), value_a.value(), lambda.value(), precision);
        // The tested point a + j (b - a) 2^-g, j = round(lambda 2^g) kept in
        // [1, 2^g - 1] so that the cell on either side of it lies in [a, b].
        Integer j, last;
        Float cell, point, other;
        arf_mul_2exp_si(point.value(), arb_midref(lambda.value()), g);
        arf_get_fmpz(j.value(), point.value(), ARF_RND_NEAR);
        fmpz_one(last.value());
        fmpz_mul_2exp(last.value(), last.value(), static_cast<ulong>(g));
        fmpz_sub_ui(last.value(), last.value(), 1);
        if (fmpz_cmp_ui(j.value(), 1) < 0)
            fmpz_one(j.value());
        if (fmpz_cmp(j.value(), last.value()) > 0)
            fmpz_set(j.value(), last.value());
        arf_sub(cell.value(), b, a, ARF_PREC_EXACT, ARF_RND_DOWN);
        arf_mul_2exp_si(cell.value(), cell.value(), -g);
        arf_mul_fmpz(point.value(), cell.value(), j.value(), ARF_PREC_EXACT, ARF_RND_DOWN);
        arf_add(point.value(), point.value(), a, ARF_PREC_EXACT, ARF_RND_DOWN);
        const int at_point = sign_at(point.value(), work, c.sign_precision, options);
        if (at_point == 0)
            return exact_root(point.value());
        // The root is to the right of the point when p keeps a's sign there.
        const bool right = at_point == sign_a;
        if (right)
            arf_add(other.value(), point.value(), cell.value(), ARF_PREC_EXACT, ARF_RND_DOWN);
        else
            arf_sub(other.value(), point.value(), cell.value(), ARF_PREC_EXACT, ARF_RND_DOWN);
        const int at_other = sign_at(other.value(), work, c.sign_precision, options);
        if (at_other == 0)
            return exact_root(other.value());
        const bool confirmed = right ? at_other != sign_a : at_other == sign_a;
        if (confirmed) {
            // The root lies in the cell between the two points.
            arf_set(right ? a : b, point.value());
            arf_set(right ? b : a, other.value());
            c.grid_bits = std::min<ulong>(resources::precision_max, 2 * c.grid_bits);
            return Cells::Cell;
        }
        // Refuted: the root lies beyond the second point, on a side a whole
        // number of cells wide. The side is widened outward to the next power
        // of two of cells, which stays within [a, b]. Its ends are then on the cells' grid, as the ends of
        // a confirmed cell are, and its radius is an exact magnitude, as
        // exact_interval requires: an end of another width (any whole number
        // of cells) would carry its bits into every later point, and short
        // dyadic roots that the grid holds would never be met exactly.
        // History: docs/code-history-20261021.md, real_refinement.cpp 5.
        Integer count;
        if (right)
            fmpz_sub(count.value(), last.value(), j.value());
        else
            fmpz_sub_ui(count.value(), j.value(), 1);
        const auto widths = static_cast<slong>(fmpz_clog_ui(count.value(), 2));
        c.grid_bits = std::max<ulong>(1, c.grid_bits / 2);
        // A side of more than half the cells widens to the whole interval:
        // no progress. Bisect instead (one more sign, at the midpoint).
        if (widths >= g)
            return Cells::None;
        Float side;
        arf_mul_2exp_si(side.value(), cell.value(), widths);
        if (right)
            arf_sub(a, b, side.value(), ARF_PREC_EXACT, ARF_RND_DOWN);
        else
            arf_add(b, a, side.value(), ARF_PREC_EXACT, ARF_RND_DOWN);
        return Cells::Side;
    }
    cps_stop_reason step(const cps_options &options) {
        if (options.cancelled())
            return CPS_STOP_CANCELLED;
        if (!pending.empty())
            return isolate_step(options);
        if (components.empty())
            return CPS_STOP_NONE;
        set_goal(options);
        bool all_ready = std::all_of(components.begin(), components.end(),
                                     [&](const auto &c) { return goal(c); });
        if (all_ready) {
            if (disjoint(roots()))
                return CPS_STOP_NONE;
            if (extra_bits >= options.max_work_bits)
                return CPS_STOP_PRECISION;
            extra_bits = std::min<uint64_t>(options.max_work_bits, extra_bits + 32);
            // The tighter target may still be met by every individual box,
            // especially exact dyadic roots rounded to overlapping exports.
            // Continue tightening on a later controller step instead of
            // mistaking the absence of a selected component for a hard limit.
            return CPS_STOP_NONE;
        }
        Component *selected = nullptr;
        for (size_t i = 0; i < components.size(); ++i) {
            auto &candidate = components[cursor++ % components.size()];
            if (!goal(candidate) && candidate.capped_at != options.max_work_bits) {
                selected = &candidate;
                break;
            }
        }
        // Every unfinished component is at the cap.
        if (!selected)
            return CPS_STOP_PRECISION;
        // A component that cannot go on within the cap is set aside; the
        // others go on, and the stop comes when none can. Stopping at the
        // first one would leave its neighbours unrefined: two close roots
        // could still share an isolation endpoint, so the cover at the
        // stop would not be disjoint and would be lost. History:
        // docs/code-history-20261021.md, real_refinement.cpp 6.
        const auto reason = refine(*selected, options);
        if (reason != CPS_STOP_PRECISION)
            return reason;
        selected->capped_at = options.max_work_bits;
        return CPS_STOP_NONE;
    }
    // One refinement step of a component short of its goal.
    cps_stop_reason refine(Component &c, const cps_options &options) {
        ++extended.real_refinement_steps;
        StageTimer timer(statistics.refinement_ns);
        if (c.rational) {
            if (c.precision >= static_cast<slong>(options.max_work_bits))
                return CPS_STOP_PRECISION;
            c.precision = static_cast<slong>(std::min<uint64_t>(
                options.max_work_bits,
                std::max<uint64_t>(static_cast<uint64_t>(c.precision) * 2, static_cast<uint64_t>(resources::precision_sum(resources::precision_sum(options.bits, extra_bits), 32)))));
            arb_set_fmpq(c.interval.value(), c.rational->value(), c.precision);
            ++statistics.precision_increases;
            return CPS_STOP_NONE;
        }
        const slong work = c.precision;
        RealBall midpoint, value, slope, next, clipped;
        arb_get_mid_arb(midpoint.value(), c.interval.value());
        arb_fmpz_poly_evaluate_arb(value.value(), polynomial.value(), midpoint.value(), work);
        // The exact value zero is found exactly; never demand a positive
        // relative accuracy from it.
        const auto midpoint_root = [&] {
            c.rational.emplace();
            arf_get_fmpq(c.rational->value(), arb_midref(midpoint.value()));
            arb_set(c.interval.value(), midpoint.value());
            return CPS_STOP_NONE;
        };
        if (arb_is_zero(value.value()))
            return midpoint_root();
        arb_fmpz_poly_evaluate_arb(slope.value(), derivative.value(), c.interval.value(), work);
        if (!arb_contains_zero(slope.value())) {
            arb_div(next.value(), value.value(), slope.value(), work);
            arb_sub(next.value(), midpoint.value(), next.value(), work);
            if (arb_intersection(clipped.value(), next.value(), c.interval.value(), work) &&
                arb_contains(c.interval.value(), clipped.value()) &&
                mag_cmp(arb_radref(clipped.value()), arb_radref(c.interval.value())) < 0) {
                arb_set(c.interval.value(), clipped.value());
                return CPS_STOP_NONE;
            }
        }
        Float a, b;
        arb_get_interval_arf(a.value(), b.value(), c.interval.value(), ARF_PREC_EXACT);
        int sign_a = 0;
        const auto refined = refine_by_cells(c, a.value(), b.value(), sign_a, work, options);
        if (refined == Cells::Root)
            return CPS_STOP_NONE;
        if (refined == Cells::None) {
            // Bisection: the signs at the midpoint and at a.
            const int sign = arb_is_positive(value.value())   ? 1
                             : arb_is_negative(value.value()) ? -1
                                                              : sign_at(arb_midref(midpoint.value()),
                                                                        std::min(static_cast<slong>(options.max_work_bits),
                                                                                 resources::precision_sum(work, work)),
                                                                        c.sign_precision, options);
            if (sign == 0)
                return midpoint_root();
            if (!sign_a)
                sign_a = sign_at(a.value(), work, c.sign_precision, options);
            if (sign_a == 0) {
                c.rational.emplace();
                arf_get_fmpq(c.rational->value(), a.value());
                arb_set_arf(c.interval.value(), a.value());
                return CPS_STOP_NONE;
            }
            if (sign_a == sign)
                arf_set(a.value(), arb_midref(midpoint.value()));
            else
                arf_set(b.value(), arb_midref(midpoint.value()));
        }
        const slong endpoint_bits = resources::precision_sum(std::max(arf_bits(a.value()), arf_bits(b.value())), 2);
        if (endpoint_bits > work) {
            if (work == static_cast<slong>(options.max_work_bits))
                return CPS_STOP_PRECISION;
            c.precision = static_cast<slong>(std::min<uint64_t>(
                options.max_work_bits,
                std::max<uint64_t>(static_cast<uint64_t>(work) * 2, static_cast<uint64_t>(endpoint_bits))));
            ++statistics.precision_increases;
        }
        exact_interval(next.value(), a.value(), b.value());
        if (!arb_contains(c.interval.value(), next.value()) ||
            mag_cmp(arb_radref(next.value()), arb_radref(c.interval.value())) >= 0) {
            if (c.precision == static_cast<slong>(options.max_work_bits))
                return CPS_STOP_PRECISION;
            c.precision = static_cast<slong>(
                std::min<uint64_t>(options.max_work_bits, static_cast<uint64_t>(c.precision) * 2));
            ++statistics.precision_increases;
            return CPS_STOP_NONE;
        }
        arb_set(c.interval.value(), next.value());
        if (refined == Cells::None)
            ++extended.real_bisections;
        return CPS_STOP_NONE;
    }
};
RealRefinement::RealRefinement(const fmpz_poly_t p, uint64_t m, slong prec, certpolysolve::detail::SearchStatistics &s,
                               certpolysolve::detail::CandidateStatistics &e)
    : impl_(std::make_unique<Impl>(p, m, prec, s, e)) {}
RealRefinement::~RealRefinement() = default;
cps_stop_reason RealRefinement::step(const cps_options &o) { return impl_->step(o); }
bool RealRefinement::cover(std::vector<Root> &out) const {
    if (!impl_->pending.empty())
        return false;
    auto roots = impl_->roots();
    if (!disjoint(roots))
        return false;
    out = std::move(roots);
    return true;
}
uint64_t RealRefinement::count() const { return impl_->components.size() * impl_->multiplicity; }
uint64_t RealRefinement::precision() const {
    uint64_t p = static_cast<uint64_t>(impl_->initial_precision);
    for (const auto &c : impl_->components)
        p = std::max(p, static_cast<uint64_t>(c.precision));
    return p;
}
uint64_t RealRefinement::cells() const { return impl_->pending.size() + impl_->components.size(); }
} // namespace certpolysolve::detail
