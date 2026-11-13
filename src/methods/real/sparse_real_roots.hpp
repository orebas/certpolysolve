#pragma once
// B6: the real roots of a real polynomial with few terms, isolated by the recursion through its fractional
// derivatives (Sagraloff 2014; docs/b6-sparse-real-roots-note-20261015.md), with certified sign evaluations only:
// every sign is ball arithmetic on exact coefficients of the exact derivative chain, evaluated term by term
// (sparse_real_evaluator.hpp), so the cost is counted in terms, never in the degree.
//
// On the positive axis, after stripping x^(i_0): p_0 = p and p_(j+1) = x^(-(e_1 - 1)) p_j', where e_1 is the
// lowest positive exponent of p_j, so p_(j+1) has one term fewer and a nonzero constant term; p_(k-1) is a
// constant. Bottom up, the positive roots of p_(j+1) (isolating intervals with exact ends and proved opposite
// signs), with 0 and a root bound R, are the breakpoints of level j: between two consecutive breakpoints p_j is
// strictly monotone, so it has exactly one root there if its certified signs at the two ends differ, none
// otherwise. The sign of p_j "at" a breakpoint is its sign over the breakpoint's whole interval; undecided, the
// interval is refined (real_refine on p_(j+1)) and the precision raised, in turn, up to the cap: a sign that
// stays undecided there is a possible common root of p_j and p_(j+1), and the isolation stops honestly.
//
// The negative roots are the positive roots of p(-x). Zero is a root of multiplicity i_0 and is reported apart.
#include "core/query.hpp"
#include "certify/real_sign_refinement.hpp"
#include "numeric/evaluate/sparse_real_evaluator.hpp"
#include <flint/arb.h>
#include <flint/fmpq.h>
#include <cstdint>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace certpolysolve::detail {

// A positive root of a level's polynomial: the ends of an interval holding exactly one root of it, with the
// polynomial's proved opposite signs at them (sign_low at `low`), or an exact dyadic root as an interval of width
// zero (low = high the root, sign_low 0). An end is a dyadic point or a tiny ball (a midpoint rounded at the
// working precision) with the sign proved over all of it, so no root lies in an end and the closed interval of
// the outer bounds holds the same single root.
struct SparseRootInterval {
    RealBall low, high;
    int sign_low = 0;
    // The level's monotone interval the root was found in (the consecutive breakpoints' facing ends): the proved
    // one-root region, kept as found while low and high are refined inside it.
    RealBall region_low, region_high;
};

struct SparseRealIsolation {
    bool complete = false;        // every level decided; the roots below are all of the real roots
    bool work_limited = false;    // the stop is the ceiling (Settings::evaluation_limit)
    std::string stop;             // why it is not complete
    uint64_t zero_multiplicity = 0;
    std::vector<SparseRootInterval> positive, negative; // ascending; the negative ones as roots of p(-x)
    uint64_t evaluations = 0, refinements = 0, levels = 0;
    uint64_t multiplications = 0; // the evaluators' counted real multiplications (the work meter's unit)
    double weighted_work = 0;     // the multiplications weighted by their precision (operation_weight)
    slong deepest_bits = 0;       // the highest precision a sign needed
};

class SparseRealIsolator {
  public:
    // The chain's precisions: the first evaluation's, the cap (never exceeded), and the guard added to a goal's
    // width in the bracket arithmetic (real_sign_refinement.hpp).
    struct Settings {
        slong start_bits = 64;
        slong cap_bits = 65536;
        slong guard_bits = 64;
        uint64_t evaluation_limit = std::numeric_limits<uint64_t>::max(); // the ceiling; see RealSignBudget
    };

    // `terms` ascending, distinct exponents, exact coefficients; the polynomial's nonzero part. Throws on an empty
    // or constant list.
    SparseRealIsolator(std::vector<SparseRealTerm> terms, Settings settings) : settings_(settings) {
        if (terms.size() < 2) throw std::invalid_argument("a sparse isolation needs at least two terms");
        zero_ = static_cast<uint64_t>(terms.front().exponent);
        for (auto &term : terms) term.exponent -= static_cast<slong>(zero_);
        terms_ = std::move(terms);
    }

    // The whole isolation, both half-axes.
    SparseRealIsolation isolate(const cps_options &o) {
        SparseRealIsolation out;
        out.zero_multiplicity = zero_;
        if (!half_axis(terms_, out.positive, out, o)) return out;
        std::vector<SparseRealTerm> reflected = terms_;
        for (auto &term : reflected)
            if (term.exponent & 1) fmpq_neg(term.coefficient.significand.value(), term.coefficient.significand.value());
        if (!half_axis(reflected, out.negative, out, o)) return out;
        out.complete = true;
        return out;
    }

  private:
    Settings settings_;
    uint64_t zero_ = 0;
    std::vector<SparseRealTerm> terms_;

    // p_(j+1) from p_j: the terms with positive exponents, each coefficient times its exponent, the exponents
    // lowered by the least of them (the derivative divided by x^(e_1 - 1)).
    static std::vector<SparseRealTerm> fractional_derivative(const std::vector<SparseRealTerm> &p) {
        std::vector<SparseRealTerm> next;
        slong lowest = 0;
        for (const auto &term : p) {
            if (term.exponent == 0) continue;
            if (next.empty()) lowest = term.exponent - 1;
            SparseRealTerm d;
            d.exponent = term.exponent - 1 - lowest;
            d.coefficient = term.coefficient;
            fmpq_mul_si(d.coefficient.significand.value(), d.coefficient.significand.value(), term.exponent);
            next.push_back(std::move(d));
        }
        return next;
    }

    // Cauchy's bound as a power of two: every root has |x| < 1 + max_(j<n) |c_j| / |c_n| <= R.
    static bool root_bound(RealBall &bound, const std::vector<SparseRealTerm> &p) {
        constexpr slong bound_bits = 64;
        RealBall lead, c, ratio, largest;
        p.back().coefficient.enclose(lead.value(), bound_bits);
        arb_abs(lead.value(), lead.value());
        for (size_t j = 0; j + 1 < p.size(); ++j) {
            p[j].coefficient.enclose(c.value(), bound_bits);
            arb_abs(c.value(), c.value());
            arb_div(ratio.value(), c.value(), lead.value(), bound_bits);
            arb_max(largest.value(), largest.value(), ratio.value(), bound_bits);
        }
        arb_add_ui(largest.value(), largest.value(), 1, bound_bits);
        Float upper;
        arb_get_ubound_arf(upper.value(), largest.value(), bound_bits);
        if (!arf_is_finite(upper.value()) || arf_sgn(upper.value()) <= 0) return false;
        // Keep the exponent exact: the machine-word variant saturates, which can put this bound BELOW a root
        // of compact input such as x^128 - 2^(128 k) with k > WORD_MAX.
        Integer exponent;
        arf_abs_bound_lt_2exp_fmpz(exponent.value(), upper.value());
        arb_one(bound.value());
        arb_mul_2exp_fmpz(bound.value(), bound.value(), exponent.value());
        return true;
    }

    // The sign of p over the ball `x` (a point, or a breakpoint's whole interval) at `bits`: one evaluation,
    // charged to the budget; 0 when undecided (or over the work ceiling).
    static int sign_over(const SparseRealEvaluator &p, const arb_t x, slong bits, const RealSignBudget &budget) {
        RealBall value, slope;
        charged_evaluate2(p, value.value(), slope.value(), x, bits, budget);
        return budget.over_limit ? 0 : real_sign_of(value.value());
    }

    // A breakpoint of level j: a point (0 or R, exact) or a root interval of p_(j+1).
    struct Breakpoint {
        RealBall low, high; // equal for a point
        bool is_root = false;
        size_t root_index = 0; // into the level above's roots, for refinement
    };

    // The sign of p_j over a breakpoint, refining the breakpoint's interval (on p_(j+1)) and raising the precision
    // in turn until it is decided; false at the cap.
    bool breakpoint_sign(const SparseRealEvaluator &p, const SparseRealEvaluator *above, Breakpoint &b,
                         std::vector<SparseRootInterval> &above_roots, int &sign, SparseRealIsolation &out,
                         const cps_options &o, RealSignBudget &budget) {
        slong bits = settings_.start_bits;
        slong goal = settings_.start_bits;
        for (;;) {
            RealBall x;
            if (b.is_root) arb_union(x.value(), b.low.value(), b.high.value(), bits);
            else arb_set(x.value(), b.low.value());
            sign = sign_over(p, x.value(), bits, budget);
            out.deepest_bits = std::max(out.deepest_bits, bits);
            if (sign) return true;
            if (budget.over_limit) {
                out.work_limited = true;
                out.stop = "the recursion exceeded its work ceiling";
                return false;
            }
            if (o.cancelled()) {
                out.stop = "cancelled";
                return false;
            }
            // Undecided: a root interval wider than the precision resolves is refined first, else the precision
            // doubles; both bounded by the cap.
            // The goal stays a guard below the cap: the bracket's arithmetic works at goal + guard bits, and the
            // signs that close it need that margin under the cap.
            const slong goal_cap = settings_.cap_bits - settings_.guard_bits;
            if (b.is_root && goal < goal_cap && above) {
                goal = std::min(2 * goal, goal_cap);
                ++out.refinements;
                auto &root = above_roots[b.root_index];
                if (!real_refine(*above, root.low, root.high, goal, budget, o)) {
                    out.work_limited = budget.over_limit;
                    out.stop = budget.over_limit ? "the recursion exceeded its work ceiling"
                               : o.cancelled()   ? "cancelled"
                                                 : "a breakpoint's root could not be refined (a possible common root "
                                                   "of consecutive derivatives, or the cap)";
                    return false;
                }
                arb_set(b.low.value(), root.low.value());
                arb_set(b.high.value(), root.high.value());
                if (bits < goal) bits = std::min(2 * bits, settings_.cap_bits);
                continue;
            }
            if (bits >= settings_.cap_bits) {
                out.stop = b.is_root ? "a sign at a breakpoint undecided at the cap (a possible common root of "
                                       "consecutive derivatives)"
                                     : "a sign at an end undecided at the cap";
                return false;
            }
            bits = std::min(2 * bits, settings_.cap_bits);
        }
    }

    // The positive roots of the chain's bottom level up to p_0 = `top`; false on a stop (out.stop says why).
    bool half_axis(const std::vector<SparseRealTerm> &top, std::vector<SparseRootInterval> &roots,
                   SparseRealIsolation &out, const cps_options &o) {
        // The chain p_0 .. p_(k-1), and an evaluator for each.
        std::vector<std::vector<SparseRealTerm>> chain{top};
        while (chain.back().size() > 1) chain.push_back(fractional_derivative(chain.back()));
        std::vector<std::unique_ptr<SparseRealEvaluator>> evaluators;
        for (const auto &level : chain) evaluators.push_back(std::make_unique<SparseRealEvaluator>(level));
        out.levels = std::max<uint64_t>(out.levels, chain.size());
        RealBall bound;
        if (!root_bound(bound, top)) {
            out.stop = "a finite Cauchy bound could not be established";
            return false;
        }
        RealSignBudget budget{settings_.start_bits, settings_.cap_bits, settings_.guard_bits, out.evaluations, true};
        budget.weighted_work = out.weighted_work; // the other half-axis's work is counted too
        budget.evaluation_limit = settings_.evaluation_limit; // one ceiling for both half-axes (out.evaluations is shared)
        struct CountOnExit {
            const std::vector<std::unique_ptr<SparseRealEvaluator>> &evaluators;
            const RealSignBudget &budget;
            SparseRealIsolation &out;
            ~CountOnExit() {
                for (const auto &evaluator : evaluators) out.multiplications += evaluator->multiplications();
                out.weighted_work = budget.weighted_work;
            }
        } counting{evaluators, budget, out};
        // Bottom up: the constant level has no root.
        std::vector<SparseRootInterval> above;
        for (size_t level = chain.size() - 1; level-- > 0;) {
            const auto &p = *evaluators[level];
            const SparseRealEvaluator *above_p = evaluators[level + 1].get();
            // The breakpoints: 0, the roots above, R.
            std::vector<Breakpoint> points;
            points.push_back({});
            arb_zero(points.back().low.value());
            arb_zero(points.back().high.value());
            for (size_t r = 0; r < above.size(); ++r) {
                Breakpoint b;
                arb_set(b.low.value(), above[r].low.value());
                arb_set(b.high.value(), above[r].high.value());
                b.is_root = true;
                b.root_index = r;
                points.push_back(std::move(b));
            }
            points.push_back({});
            arb_set(points.back().low.value(), bound.value());
            arb_set(points.back().high.value(), bound.value());
            std::vector<int> signs(points.size());
            for (size_t i = 0; i < points.size(); ++i)
                if (!breakpoint_sign(p, above_p, points[i], above, signs[i], out, o, budget)) return false;
            std::vector<SparseRootInterval> found;
            for (size_t i = 0; i + 1 < points.size(); ++i) {
                if (signs[i] == signs[i + 1]) continue;
                // p is monotone between the breakpoints and changes sign: exactly one root in (high_i, low_(i+1)).
                SparseRootInterval root;
                arb_set(root.low.value(), points[i].high.value());
                arb_set(root.high.value(), points[i + 1].low.value());
                arb_set(root.region_low.value(), root.low.value());
                arb_set(root.region_high.value(), root.high.value());
                root.sign_low = signs[i];
                found.push_back(std::move(root));
            }
            above = std::move(found);
        }
        roots = std::move(above);
        return true;
    }
};

} // namespace certpolysolve::detail
