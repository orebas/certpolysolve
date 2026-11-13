#pragma once
// Certified signs of a real polynomial at exact points, and the refinement of a bracket with proved opposite
// signs to a relative width, over an evaluator the caller supplies. These are the axis charts' routines (C1,
// piecewise_axis.cpp: source_refine over the dense input), which B6's sparse recursion
// (docs/b6-sparse-real-roots-note-20261015.md) runs over a term-by-term evaluator. History:
// docs/code-history-20261021.md, real_sign_refinement.hpp 1. The evaluator provides
//
//     void evaluate2(arb_t value, arb_t slope, const arb_t x, slong bits);
//
// enclosures of p(x) and p'(x) at `bits`, the exact coefficients enclosed at that precision. A sign is decided
// when the value's ball excludes zero.
#include "core/dovetail.hpp"
#include "core/query.hpp"
#include "core/raii.hpp"
#include "core/resource_limits.hpp"
#include <flint/arb.h>
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <limits>

namespace certpolysolve::detail {

// The precisions a sign evaluation may use, and the counter it charges.
struct RealSignBudget {
    slong start_bits; // the first evaluation's precision: the caller's working precision
    slong cap_bits;   // the precision cap, never exceeded
    slong guard_bits; // added to a goal's width for the bracket's own arithmetic (real_refine)
    uint64_t &evaluations;
    // The regime of B6's recursion, where a bracket may be wide and its ends may be balls: real_refine takes a
    // Newton step only where it shrinks the bracket by at least half (far from a root of a high-degree polynomial
    // Newton converges linearly, by a factor 1 - 1/n a step: x^4096 - 2 from [1, 2] needs 1,660 steps where the
    // pass bound allows 1,056, while bisection halves), and bisects at an exact point (a midpoint of ball ends
    // inherits their radius, and a ball holding the root has no sign at any precision). Off for the axis charts,
    // whose brackets start narrow with exact ends and whose arithmetic stays as promoted. History:
    // docs/code-history-20261021.md, real_sign_refinement.hpp 2. Also under this flag: a value that is exactly zero
    // at an evaluated point makes the point the root (Sagraloff 2014 returns the interval of width zero there): a
    // level's polynomial can have a dyadic root, which Newton lands on exactly, and an exact zero has no sign at
    // any precision. History: docs/code-history-20261021.md, real_sign_refinement.hpp 3. Also: a Newton step too
    // short to halve the bracket is the error estimate of a converged iterate, and the bracket is closed about that
    // iterate by probes at four times the step each side (the root lies within the step's order near a simple root,
    // so the probes bracket it and the bracket shrinks quadratically, where bisection would walk the far end down a
    // factor of two at a time), and a sign evaluation starts at half the precision the last one needed instead of
    // the budget's start. History: docs/code-history-20261021.md, real_sign_refinement.hpp 4.
    bool wide_brackets = false;
    mutable slong last_bits = 0; // (wide_brackets) the precision the last decided sign needed
    // (wide_brackets) The evaluations against a ceiling: past it every sign is undecided and over_limit is set.
    // Also, for evaluators that count multiplications, the work so far as counted multiplications weighted by
    // operation_weight(bits), a statistic.
    uint64_t evaluation_limit = std::numeric_limits<uint64_t>::max();
    mutable double weighted_work = 0;
    mutable bool over_limit = false;
};

// One evaluation charged to the budget (wide_brackets) where the evaluator counts its multiplications.
template <class Evaluator>
void charged_evaluate2(Evaluator &p, arb_t value, arb_t slope, const arb_t x, slong bits, const RealSignBudget &budget) {
    ++budget.evaluations;
    if (budget.wide_brackets && budget.evaluations > budget.evaluation_limit) {
        budget.over_limit = true;
        arb_zero(value);
        arb_zero(slope);
        return;
    }
    if constexpr (requires { p.multiplications(); }) {
        if (budget.wide_brackets) {
            const uint64_t before = p.multiplications();
            p.evaluate2(value, slope, x, bits);
            budget.weighted_work += static_cast<double>(p.multiplications() - before) * operation_weight(bits);
            return;
        }
    }
    p.evaluate2(value, slope, x, bits);
}

// The precision cap of a query as a precision: the cap, or half the largest slong where the cap exceeds it (the
// doubling below never overflows).
inline slong real_sign_cap(const cps_options &o) {
    return static_cast<slong>(std::min<uint64_t>(o.max_work_bits, static_cast<uint64_t>(std::numeric_limits<slong>::max() / 2)));
}

inline int real_sign_of(const arb_t x) { return arb_is_positive(x) ? 1 : arb_is_negative(x) ? -1 : 0; }

// Sign of p at a point x (exact), with p and p' there; the precision doubles from the budget's start while the
// sign is undecided, up to its cap. 0 when still undecided there, or cancelled; also 0 with `exact_zero` set when
// the value is exactly zero (wide_brackets): x is a root.
template <class Evaluator>
int real_sign(Evaluator &p, const arb_t x, arb_t value, arb_t slope, const RealSignBudget &budget, const cps_options &o,
              bool *exact_zero = nullptr) {
    if (exact_zero) *exact_zero = false;
    const slong first = budget.wide_brackets ? std::max(budget.start_bits, std::min(budget.last_bits / 2, budget.cap_bits))
                                             : budget.start_bits;
    for (slong bits = first;; bits = std::min(2 * bits, budget.cap_bits)) {
        charged_evaluate2(p, value, slope, x, bits, budget);
        if (budget.over_limit) return 0;
        if (const int sign = real_sign_of(value)) {
            if (budget.wide_brackets) budget.last_bits = bits;
            return sign;
        }
        if (budget.wide_brackets && arb_is_zero(value)) {
            if (exact_zero) *exact_zero = true;
            return 0;
        }
        if (bits >= budget.cap_bits || o.cancelled()) return 0;
    }
}

// Refines a one-root interval [lo, hi] (exact ends of proved opposite signs) to the goal's relative width: Newton
// on p from the middle, kept inside the bracket (a step leaving it is a bisection), the bracket shrinking to every
// evaluated point by its sign. Near the root a step shorter than a quarter of the goal's width is closed by two
// probes about the iterate. False where a sign stays undecided at the cap, or where the ends' signs are not
// opposite and decided. This is the axis charts' routine; under wide_brackets the call goes to real_refine_wide
// instead. History: docs/code-history-20261021.md, real_sign_refinement.hpp 5.
template <class Evaluator>
bool real_refine_narrow(Evaluator &p, RealBall &lo, RealBall &hi, slong goal_bits, const RealSignBudget &budget,
                        const cps_options &o) {
    // The bracket's arithmetic at the goal's width and a guard: a midpoint rounded to fewer bits than the
    // bracket's would not move it.
    const slong work = std::max(budget.start_bits, resources::precision_sum(goal_bits, budget.guard_bits));
    RealBall value, slope, x, width, limit, step, probe;
    const int s_lo = real_sign(p, lo.value(), value.value(), slope.value(), budget, o);
    if (!s_lo || real_sign(p, hi.value(), value.value(), slope.value(), budget, o) != -s_lo) return false;
    const auto reached = [&] {
        arb_sub(width.value(), hi.value(), lo.value(), work);
        RealBall a, b;
        arb_abs(a.value(), lo.value());
        arb_abs(b.value(), hi.value());
        arb_min(limit.value(), a.value(), b.value(), work);
        arb_mul_2exp_si(limit.value(), limit.value(), -goal_bits);
        return arb_le(width.value(), limit.value());
    };
    const auto midpoint = [&](RealBall &m) {
        arb_add(m.value(), lo.value(), hi.value(), work);
        arb_mul_2exp_si(m.value(), m.value(), -1);
    };
    const auto move = [&](const RealBall &point, int sign) {
        if (sign == s_lo) arb_set(lo.value(), point.value());
        else arb_set(hi.value(), point.value());
    };
    midpoint(x);
    // Each pass halves the bracket at least (a bisection) or converges quadratically; 4 bisections per bit of the
    // goal and its guard bound the passes generously.
    for (slong pass = 0; pass < 4 * (goal_bits + budget.guard_bits) && !o.cancelled(); ++pass) {
        if (reached()) return true;
        const int sign = real_sign(p, x.value(), value.value(), slope.value(), budget, o);
        if (!sign) return false;
        move(x, sign);
        if (reached()) return true;
        // Newton step x - p/p', exact midpoint; a bisection where it leaves (lo, hi)
        RealBall next;
        bool newton = arb_is_nonzero(slope.value());
        if (newton) {
            arb_div(step.value(), value.value(), slope.value(), work);
            arb_sub(next.value(), x.value(), step.value(), work);
            mag_zero(arb_radref(next.value()));
            newton = arb_gt(next.value(), lo.value()) && arb_lt(next.value(), hi.value());
        }
        if (!newton) {
            midpoint(x);
            continue;
        }
        // Close the bracket about a converged iterate: probes a quarter of the goal's width each side.
        RealBall quarter;
        arb_mul_2exp_si(quarter.value(), limit.value(), -2);
        arb_abs(step.value(), step.value());
        if (arb_lt(step.value(), quarter.value())) {
            for (const int side : {-1, 1}) {
                arb_set(probe.value(), next.value());
                if (side < 0) arb_sub(probe.value(), probe.value(), quarter.value(), work);
                else arb_add(probe.value(), probe.value(), quarter.value(), work);
                mag_zero(arb_radref(probe.value()));
                if (!arb_gt(probe.value(), lo.value()) || !arb_lt(probe.value(), hi.value())) continue;
                const int probe_sign = real_sign(p, probe.value(), value.value(), slope.value(), budget, o);
                if (!probe_sign) return false;
                move(probe, probe_sign);
            }
            if (reached()) return true;
        }
        arb_set(x.value(), next.value());
    }
    return false;
}

// The wide regime (B6: a bracket from the recursion, possibly 0 to Cauchy's bound, its ends points or tiny
// balls). Two phases. While the bracket is not local (its lower end zero, or its ends' exponents two or more
// apart), bisection in the exponent: at hi / 2^(2^j) while the lower end is zero, j counting those steps (their
// number grows as the logarithm of the binades crossed), else at the power of two halfway between the ends'
// exponents. Then, every pass: a bisection at an exact midpoint (the halving is guaranteed), then Newton
// from it with the bracket closed about the iterate by probes at four times the step each side (at least the
// goal's quarter width): near a simple root the probes bracket the root, so the bracket shrinks quadratically; a
// probe that does not bracket still moves its end. A value exactly zero at any evaluated point is the root (an
// interval of width zero). False where a sign is undecided at the cap, over the work ceiling, or at the pass
// bound; the ends' signs must be opposite and decided. History:
// docs/code-history-20261021.md, real_sign_refinement.hpp 6.
template <class Evaluator>
bool real_refine_wide(Evaluator &p, RealBall &lo, RealBall &hi, slong goal_bits, const RealSignBudget &budget,
                      const cps_options &o) {
    const slong work = std::max(budget.start_bits, resources::precision_sum(goal_bits, budget.guard_bits));
    RealBall value, slope, x, width, limit, step, probe, next, quarter, delta;
    bool exact = false;
    const auto exact_root = [&](const RealBall &point) {
        arb_set(lo.value(), point.value());
        arb_set(hi.value(), point.value());
        return true;
    };
    const int s_lo = real_sign(p, lo.value(), value.value(), slope.value(), budget, o, &exact);
    if (exact) return exact_root(lo);
    const int s_hi = real_sign(p, hi.value(), value.value(), slope.value(), budget, o, &exact);
    if (exact) return exact_root(hi);
    if (!s_lo || s_hi != -s_lo) return false;
    const auto reached = [&] {
        arb_sub(width.value(), hi.value(), lo.value(), work);
        RealBall a, b;
        arb_abs(a.value(), lo.value());
        arb_abs(b.value(), hi.value());
        arb_min(limit.value(), a.value(), b.value(), work);
        arb_mul_2exp_si(limit.value(), limit.value(), -goal_bits);
        return arb_le(width.value(), limit.value());
    };
    const auto move = [&](const RealBall &point, int sign) {
        if (sign == s_lo) arb_set(lo.value(), point.value());
        else arb_set(hi.value(), point.value());
    };
    // A sign at a point inside the bracket, moving its end: 1 decided, 0 undecided (or the ceiling), -1 exact root.
    const auto decide = [&](const RealBall &point) {
        const int sign = real_sign(p, point.value(), value.value(), slope.value(), budget, o, &exact);
        if (exact) return -1;
        if (!sign) return 0;
        move(point, sign);
        return 1;
    };
    // Exponents are compact integers too: neither the gallop nor the midpoint may clamp them to a machine word.
    Integer zero_step, e_lo, e_hi, gap, middle;
    fmpz_set_si(zero_step.value(), -1);
    // A research trace, compiled only with CERTPOLYSOLVE_DEBUG_PROBES: a release library reads no environment variable
    // and writes nothing to stderr (docs/thread-safety-audit-20261018.md, A3). History:
    // docs/code-history-20261021.md, real_sign_refinement.hpp 7.
#ifdef CERTPOLYSOLVE_DEBUG_PROBES
    static const bool trace = std::getenv("CERTPOLYSOLVE_TRACE_REFINE") != nullptr;
#else
    constexpr bool trace = false;
#endif
    for (slong pass = 0; pass < 4 * (goal_bits + budget.guard_bits) && !o.cancelled(); ++pass) {
        if (reached()) return true;
        if (trace && (pass < 40 || pass % 32 == 0))
            std::fprintf(stderr, "wide pass %ld: lo exp %ld hi exp %ld width exp %ld limit exp %ld\n", (long)pass,
                         (long)(arb_is_zero(lo.value()) ? -99999 : arf_abs_bound_lt_2exp_si(arb_midref(lo.value()))),
                         (long)arf_abs_bound_lt_2exp_si(arb_midref(hi.value())),
                         (long)arf_abs_bound_lt_2exp_si(arb_midref(width.value())),
                         (long)(arb_is_zero(limit.value()) ? -99999 : arf_abs_bound_lt_2exp_si(arb_midref(limit.value()))));
        // Phase one: not yet local.
        if (arb_is_zero(lo.value())) {
            arb_mul_2exp_fmpz(x.value(), hi.value(), zero_step.value());
            fmpz_mul_2exp(zero_step.value(), zero_step.value(), 1);
            mag_zero(arb_radref(x.value()));
            const int d = decide(x);
            if (d < 0) return exact_root(x);
            if (!d) return false;
            continue;
        }
        arf_abs_bound_lt_2exp_fmpz(e_lo.value(), arb_midref(lo.value()));
        arf_abs_bound_lt_2exp_fmpz(e_hi.value(), arb_midref(hi.value()));
        fmpz_sub(gap.value(), e_hi.value(), e_lo.value());
        if (fmpz_cmp_ui(gap.value(), 3) >= 0) {
            // lo < 2^e_lo <= 2^floor((e_lo + e_hi) / 2) <= 2^(e_hi - 2) < hi: the floor, not C++'s division toward
            // zero, which at negative exponents could land on an upper end that is a power of two (a stuck loop).
            fmpz_add(middle.value(), e_lo.value(), e_hi.value());
            fmpz_fdiv_q_2exp(middle.value(), middle.value(), 1);
            arb_one(x.value());
            arb_mul_2exp_fmpz(x.value(), x.value(), middle.value());
            if (arb_gt(x.value(), lo.value()) && arb_lt(x.value(), hi.value())) {
                const int d = decide(x);
                if (d < 0) return exact_root(x);
                if (!d) return false;
                continue;
            }
        }
        // Phase two: a bisection, then Newton closed by probes.
        arb_add(x.value(), lo.value(), hi.value(), work);
        arb_mul_2exp_si(x.value(), x.value(), -1);
        mag_zero(arb_radref(x.value()));
        const int d = decide(x);
        if (d < 0) return exact_root(x);
        if (!d) return false;
        if (reached()) return true;
        if (!arb_is_nonzero(slope.value())) continue;
        arb_div(step.value(), value.value(), slope.value(), work);
        arb_sub(next.value(), x.value(), step.value(), work);
        mag_zero(arb_radref(next.value()));
        if (!(arb_gt(next.value(), lo.value()) && arb_lt(next.value(), hi.value()))) continue;
        // The probes' distance: four times the step (the error's order near a simple root), at least the goal's
        // quarter width, at most a quarter of the bracket (so both probes lie inside it while the step is still
        // of the bracket's order).
        arb_abs(step.value(), step.value());
        arb_mul_2exp_si(quarter.value(), limit.value(), -2);
        arb_mul_2exp_si(delta.value(), step.value(), 2);
        RealBall bracket_quarter;
        arb_sub(bracket_quarter.value(), hi.value(), lo.value(), work);
        arb_mul_2exp_si(bracket_quarter.value(), bracket_quarter.value(), -2);
        if (arb_gt(delta.value(), bracket_quarter.value())) arb_set(delta.value(), bracket_quarter.value());
        if (arb_lt(delta.value(), quarter.value())) arb_set(delta.value(), quarter.value());
        for (const int side : {-1, 1}) {
            arb_set(probe.value(), next.value());
            if (side < 0) arb_sub(probe.value(), probe.value(), delta.value(), work);
            else arb_add(probe.value(), probe.value(), delta.value(), work);
            mag_zero(arb_radref(probe.value()));
            if (!arb_gt(probe.value(), lo.value()) || !arb_lt(probe.value(), hi.value())) continue;
            const int pd = decide(probe);
            if (pd < 0) return exact_root(probe);
            if (!pd) return false;
        }
    }
    return false;
}

template <class Evaluator>
bool real_refine(Evaluator &p, RealBall &lo, RealBall &hi, slong goal_bits, const RealSignBudget &budget,
                 const cps_options &o) {
    return budget.wide_brackets ? real_refine_wide(p, lo, hi, goal_bits, budget, o)
                                : real_refine_narrow(p, lo, hi, goal_bits, budget, o);
}

} // namespace certpolysolve::detail
