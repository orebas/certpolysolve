#pragma once
// B6 in a session (docs/b6-sparse-real-roots-note-20261015.md): a real-only query's first attempt on an exact real
// input with few terms, before the axis charts and the complex route (real_only_route.cpp, which session.cpp
// chooses). The sparse recursion isolates the real roots (sparse_real_roots.hpp); each is refined by the input
// itself to the goal inside its proved one-root region and published (under the isolation goal, just far enough
// to lie inside it). A state answers one query: the route builds a fresh one for each query it answers. A stop of
// any kind (a possible common root of consecutive derivatives, the cap, the work ceiling, cancellation) leaves the
// query to the other routes, its work and statistics carried.
#include "core/heuristics.hpp"
#include "core/cover.hpp"
#include "input/coefficients/polynomial_input.hpp"
#include "core/query.hpp"
#include "methods/real/sparse_real_roots.hpp"
#include <limits>
#include <memory>
#include <string>
#include <vector>

namespace certpolysolve::detail {

class SparseRealAttempt {
  public:
    // `input` exact and real. The isolation's cap is the query's (real_sign_cap).
    explicit SparseRealAttempt(const DensePolynomial &input) {
        for (size_t j = 0; j < input.coefficients.size(); ++j)
            if (!input.coefficients[j].zero()) terms_.push_back({static_cast<slong>(j), input.coefficients[j].real});
        degree_ = static_cast<int64_t>(input.coefficients.size()) - 1;
        reflect();
    }
    // The same from the input object's terms (C40's second step): its terms are the nonzero coefficients in
    // increasing order, and its degree the last one's exponent (a session's coefficient list has no trailing zeros).
    explicit SparseRealAttempt(const cps_polynomial &input) {
        for (const auto &[exponent, coefficient] : input.terms)
            terms_.push_back({static_cast<slong>(exponent), coefficient.real});
        degree_ = input.terms.empty() ? -1 : static_cast<int64_t>(input.terms.rbegin()->first);
        reflect();
    }
    size_t support() const { return terms_.size(); }

  private:
    void reflect() {
        reflected_ = terms_;
        for (auto &term : reflected_)
            if (term.exponent & 1) fmpq_neg(term.coefficient.significand.value(), term.coefficient.significand.value());
    }

  public:
    cps_result advance(const cps_options &o) {
        cps_result result;
        result.degree = degree_;
        result.goal = engine_goal(o);
        const slong cap = real_sign_cap(o);
        SparseRealIsolator::Settings settings;
        settings.cap_bits = cap;
        evaluation_limit_ = static_cast<uint64_t>(heuristics::sparse_attempt_evaluation_ratio *
                                                  heuristics::sparse_recursion_evaluations(
                                                      static_cast<uint64_t>(degree_), terms_.size()));
        settings.evaluation_limit = evaluation_limit_;
        SparseRealIsolator isolator(terms_, settings);
        isolation_ = isolator.isolate(o);
        // A unit is a controller step: a level decided, a breakpoint refined, the preparation.
        result.progress.work_units_this_call = 1 + isolation_.levels + isolation_.refinements;
        total_evaluations_ += isolation_.evaluations;
        positive_ = std::make_unique<SparseRealEvaluator>(terms_);
        negative_ = std::make_unique<SparseRealEvaluator>(reflected_);
        if (!isolation_.complete) {
            result.termination = o.cancelled() ? CPS_CANCELLED : isolation_.work_limited ? CPS_WORK_LIMIT : CPS_PRECISION_LIMIT;
            result.progress.stop_reason = o.cancelled() ? CPS_STOP_CANCELLED
                                          : isolation_.work_limited ? CPS_STOP_WORK_BUDGET : CPS_STOP_PRECISION;
            result.diagnostics = "sparse recursion: " + isolation_.stop;
            fill_statistics(result);
            return result;
        }
        // The roots, each refined to the goal inside its region: the zero root, the negative ones ascending (as
        // roots of p(-x) in descending order), the positive ones.
        std::vector<Root> roots;
        slong work = std::max<slong>(64, isolation_.deepest_bits);
        if (isolation_.zero_multiplicity) {
            Root zero;
            acb_zero(zero.box.value());
            zero.count = zero.multiplicity = isolation_.zero_multiplicity;
            zero.realness = CPS_REAL;
            zero.certificate = "exact zero root: the lowest exponent of the input";
            zero.publish(64);
            roots.push_back(std::move(zero));
        }
        const char *failure = nullptr;
        for (size_t k = isolation_.negative.size(); k-- > 0 && !failure;)
            failure = finish_root(*negative_, isolation_.negative[k], true, o, roots, work);
        for (size_t k = 0; k < isolation_.positive.size() && !failure; ++k)
            failure = finish_root(*positive_, isolation_.positive[k], false, o, roots, work);
        fill_statistics(result);
        result.work_bits = static_cast<uint64_t>(work);
        result.progress.max_outer_bits = result.work_bits;
        if (failure) {
            result.termination = o.cancelled() ? CPS_CANCELLED : CPS_PRECISION_LIMIT;
            result.progress.stop_reason = o.cancelled() ? CPS_STOP_CANCELLED : CPS_STOP_PRECISION;
            result.diagnostics = std::string("sparse recursion: ") + failure;
            return result;
        }
        result.roots = std::move(roots);
        result.complete = result.goal_met = true;
        result.termination = CPS_SUCCESS;
        result.progress.stop_reason = CPS_STOP_NONE;
        result.diagnostics = "sparse recursion: " + std::to_string(terms_.size()) + " terms, " +
                             std::to_string(isolation_.levels) + " levels, " + std::to_string(total_evaluations_) +
                             " evaluations";
        return result;
    }

  private:
    std::vector<SparseRealTerm> terms_, reflected_;
    int64_t degree_ = -1;
    SparseRealIsolation isolation_;
    std::unique_ptr<SparseRealEvaluator> positive_, negative_;
    uint64_t total_evaluations_ = 0;
    // The ceiling on all of the attempt's evaluations, the isolation's and the publication's (a close pair's root
    // is published only inside its monotone interval, which needs the pair's separation and can cost more than
    // the complex route's whole solve). History: docs/code-history-20261021.md, sparse_real_attempt.hpp 1.
    uint64_t evaluation_limit_ = std::numeric_limits<uint64_t>::max();

    void fill_statistics(cps_result &result) const {
        result.statistics_v3.evaluations = total_evaluations_;
        result.statistics_v2.work_units = result.progress.work_units_this_call;
    }

    // Refines `root` to the goal (its bracket inside its region, by the input itself) and publishes it; the goal's
    // width first, then the goal doubled until the published box lies strictly inside the region (a ball's
    // representation widens a box by about 2^-30 of its width) and meets the goal by the library's own test.
    // History: docs/code-history-20261021.md, sparse_real_attempt.hpp 2. Doubling, not a fixed step: a root of a
    // close pair needs a box narrower than the pair's gap, which doubling reaches in a number of refinements
    // logarithmic in its bits, where a fixed step takes a number linear in them. History:
    // docs/code-history-20261021.md, sparse_real_attempt.hpp 3. An exact root (low = high) is published as its
    // point. Under the isolation goal no width is asked: the ladder starts at a 1-bit goal, a bracket within a
    // factor of two of the root, the least refinement that can move the bracket's ends off its region's ends
    // (where it starts). Returns the failure or null.
    const char *finish_root(const SparseRealEvaluator &p, SparseRootInterval &root, bool negative, const cps_options &o,
                            std::vector<Root> &roots, slong &work) {
        constexpr slong guard_bits = 64;
        constexpr slong isolation_first_goal = 1;
        const slong first_goal = o.goal == CPS_ROOT_ISOLATION ? isolation_first_goal : static_cast<slong>(o.bits) + 5;
        const slong cap = real_sign_cap(o);
        const bool exact = arb_equal(root.low.value(), root.high.value());
        for (slong goal = first_goal;; goal = resources::precision_product(goal, 2)) {
            if (goal > cap) return "a root's box does not fit its region at the precision cap";
            RealSignBudget budget{64, cap, guard_bits, total_evaluations_, true};
            budget.evaluation_limit = evaluation_limit_;
            if (!exact && !real_refine(p, root.low, root.high, goal, budget, o)) {
                if (budget.over_limit) return "the recursion exceeded its work ceiling while publishing a root";
                if (o.cancelled()) return "cancelled";
                return "a root not refined to the goal";
            }
            const slong precision = std::max(work, resources::precision_sum(goal, guard_bits));
            Root out;
            acb_zero(out.box.value());
            if (negative) {
                RealBall lo, hi;
                arb_neg(lo.value(), root.high.value());
                arb_neg(hi.value(), root.low.value());
                arb_union(acb_realref(out.box.value()), lo.value(), hi.value(), precision);
            } else {
                arb_union(acb_realref(out.box.value()), root.low.value(), root.high.value(), precision);
            }
            out.count = out.multiplicity = 1;
            out.realness = CPS_REAL;
            out.certificate = "sparse recursion: a sign change of the input inside an interval where it is monotone";
            out.publish(precision);
            bool inside = exact;
            if (!exact) {
                Float box_lo, box_hi, region_lo, region_hi;
                arb_get_lbound_arf(box_lo.value(), acb_realref(out.box.value()), precision);
                arb_get_ubound_arf(box_hi.value(), acb_realref(out.box.value()), precision);
                // The region's facing ends: their outer bounds, reflected for a negative root.
                if (negative) {
                    arb_get_ubound_arf(region_hi.value(), root.region_low.value(), precision);
                    arb_get_lbound_arf(region_lo.value(), root.region_high.value(), precision);
                    arf_neg(region_hi.value(), region_hi.value());
                    arf_neg(region_lo.value(), region_lo.value());
                } else {
                    arb_get_ubound_arf(region_lo.value(), root.region_low.value(), precision);
                    arb_get_lbound_arf(region_hi.value(), root.region_high.value(), precision);
                }
                inside = arf_cmp(box_lo.value(), region_lo.value()) > 0 && arf_cmp(box_hi.value(), region_hi.value()) < 0;
            }
            if (inside && meets_goal(out, o, precision)) {
                work = std::max(work, precision);
                roots.push_back(std::move(out));
                return nullptr;
            }
            if (exact) return "an exact root's point does not meet the goal"; // cannot happen: width zero
        }
    }
};

} // namespace certpolysolve::detail
