// Auto's direct route for exact polynomials of degree 1 to 4 (proposal C of dev/docs/easy-inputs-profile-20261025.md;
// degrees 3 and 4 since 26 October 2026, dev/docs/small-degree-closed-forms-20261026/): the closed forms in ball arithmetic,
// so that such an input skips the selector, the exact preparation, the candidate search and the per-root certificate,
// which cost it 80,000 to 400,000 instructions.
//
// direct_small_degree sets apart the exact zero roots (x^v, as preprocess does), then solves the nonzero part q:
// - degree 1: -q0 / q1;
// - degree 2: the discriminant D = q1^2 - 4 q2 q0 in exact Gaussian rationals. D = 0 is a double root -q1 / (2 q2),
//   one region of count and multiplicity 2. Otherwise the stable formula: the larger of -(q1 +- sqrt D) / 2 as Q, the
//   roots Q / q2 and q0 / Q (Vieta, no cancellation), as small_squarefree_factor in prepared_piecewise.cpp computes
//   them for piecewise's quadratic factors;
// - degree 3: q scaled to Gaussian integers by the least common multiple of its denominators, then ExactCubic
//   (small_degree_forms.hpp): its exact invariants decide square-freeness and give the multiple roots exactly, and its
//   forms give the boxes (Viete's form for three real roots, the real Cardano form for one, Cardano's form for complex
//   coefficients);
// - degree 4: the same scaling, then ExactQuartic: the exact square-free factors of a real quartic with a multiple root,
//   the biquadratic, or Descartes' two quadratic factors from a root of the resolvent (complex coefficients with a
//   multiple root are declined).
// Every value is a ball enclosing its exact expression, so each root's box holds it; the boxes are published and
// tested against the goal (meets_goal) and for disjointness, and the precision doubles from max(64, bits + 32) to the
// cap until both hold. Disjoint boxes for distinct roots hold one root each; a multiple root's box holds all its copies.
// At a spent work limit, a complete finite enumeration can instead retain disjoint unions of overlapping branches,
// after cancellation is checked. These are counted clusters, with no claim of an individual multiplicity.
//
// State: none. Entry point: direct_small_degree (SessionState::Impl::advance, session.cpp, for Auto's all-roots
// queries). Not here: piecewise's own small factors (prepared_piecewise.cpp), the general route (solver.cpp).
#include "methods/direct/small_degree.hpp"
#include "core/cover.hpp"
#include "certify/factor_cover.hpp"
#include "core/raii.hpp"
#include "core/resource_limits.hpp"
#include "methods/direct/small_degree_forms.hpp"
#include <algorithm>
#include <span>
#include <string>
#include <vector>

namespace certpolysolve::detail {
namespace {
struct Gaussian {
    Rational re, im;
};
// The same length-square, word-weighted admission as bulk exact preparation.
// This is a stage estimate, not a count of FLINT's internal operations. Charge
// the largest algebraic operand, including the invariants' growth below.
bool exact_effort_fits(size_t length, uint64_t operand_bits, const cps_options &options) {
    const auto words = resources::add(1, operand_bits / FLINT_BITS);
    return options.exact_effort_fits(resources::multiply(resources::multiply(length, length), words));
}
bool exact_operands_fit(size_t length, uint64_t operand_bits, const cps_options &options) {
    return resources::storage_fits(operand_bits, options.exact_storage_limit()) &&
           exact_effort_fits(length, operand_bits, options);
}
// The exact value of an exact coefficient within the exact storage limit, `expansion` saying whether the
// compact-exponent rule admits the list's expansion (as streamed_real_polynomial decides it); nullopt where refused.
std::optional<Gaussian> exact_value(const Coefficient &c, const cps_options &options, bool expansion) {
    Gaussian g;
    uint64_t budget = options.exact_storage_limit();
    if (!c.real.rational_with_budget(g.re.value(), budget, expansion) ||
        !c.imag.rational_with_budget(g.im.value(), budget, expansion))
        return std::nullopt;
    return g;
}
// The coefficients times the least common multiple of their denominators: Gaussian integers with the same roots.
std::optional<std::vector<GaussianInteger>> integer_coefficients(const std::vector<Gaussian> &q,
                                                                const cps_options &options) {
    Integer scale, factor;
    fmpz_one(scale.value());
    for (const auto &g : q)
        for (const auto *x : {&g.re, &g.im}) {
            const auto denominator = fmpq_denref(x->value());
            // lcm(s, d) divides s*d. Admit that uncancelled product before
            // calling FLINT; checking the final invariants would be too late.
            const auto product_bits = resources::add(fmpz_bits(scale.value()), fmpz_bits(denominator));
            if (!exact_operands_fit(q.size(), product_bits, options)) return std::nullopt;
            fmpz_lcm(scale.value(), scale.value(), denominator);
        }
    std::vector<GaussianInteger> out(q.size());
    for (size_t k = 0; k < q.size(); ++k) {
        for (const auto *x : {&q[k].re, &q[k].im}) {
            // The quotient s/d has at most bits(s) - bits(d) + 1 bits.
            const auto factor_bits = fmpz_bits(scale.value()) - fmpz_bits(fmpq_denref(x->value())) + 1;
            const auto product_bits = resources::add(fmpz_bits(fmpq_numref(x->value())), factor_bits);
            if (!exact_operands_fit(q.size(), product_bits, options)) return std::nullopt;
        }
        fmpz_divexact(factor.value(), scale.value(), fmpq_denref(q[k].re.value()));
        fmpz_mul(out[k].re.value(), fmpq_numref(q[k].re.value()), factor.value());
        fmpz_divexact(factor.value(), scale.value(), fmpq_denref(q[k].im.value()));
        fmpz_mul(out[k].im.value(), fmpq_numref(q[k].im.value()), factor.value());
    }
    return out;
}
// x y in Gaussian rationals.
Gaussian product(const Gaussian &x, const Gaussian &y) {
    Gaussian z;
    Rational t;
    fmpq_mul(z.re.value(), x.re.value(), y.re.value());
    fmpq_mul(t.value(), x.im.value(), y.im.value());
    fmpq_sub(z.re.value(), z.re.value(), t.value());
    fmpq_mul(z.im.value(), x.re.value(), y.im.value());
    fmpq_mul(t.value(), x.im.value(), y.re.value());
    fmpq_add(z.im.value(), z.im.value(), t.value());
    return z;
}
} // namespace

std::optional<cps_result> direct_small_degree(const DensePolynomial &input, const cps_options &options) {
    // The goals of one region a distinct root (as machine_newton_cover serves): the cluster-radius goal asks for the
    // coarsest regions the radius allows, which the general route finds.
    const auto goal = engine_goal(options);
    if (options.real_only || options.region.set() || (goal != CPS_ROOT_ACCURACY && goal != CPS_ROOT_APPROXIMATION))
        return std::nullopt;
    size_t length = input.coefficients.size();
    while (length && input.coefficients[length - 1].zero()) --length;
    constexpr size_t highest_degree = 4;
    if (length < 2 || length > highest_degree + 1) return std::nullopt;
    if (!std::all_of(input.coefficients.begin(), input.coefficients.begin() + static_cast<ptrdiff_t>(length),
                     [](const Coefficient &c) { return c.exact(); }))
        return std::nullopt;
    // The exact zero roots apart: q = p / x^v.
    size_t zeros = 0;
    while (zeros + 1 < length && input.coefficients[zeros].zero()) ++zeros;
    const auto *q = input.coefficients.data() + zeros;
    const size_t degree = length - 1 - zeros;
    if (degree == 0) return std::nullopt;
    // A linear formula uses ball division only. Every higher degree imports
    // exact coefficients and performs algebra before its numerical attempts.
    if (degree > 1 && !exact_effort_fits(degree + 1, 0, options)) return std::nullopt;
    bool double_root = false;
    if (degree == 2) {
        const bool expansion = options.allow_exponent_expansion ||
                               compact_expansion_fits(std::span<const Coefficient>(q, degree + 1),
                                                      options.exact_expansion_allowance());
        const auto a = exact_value(q[2], options, expansion), b = exact_value(q[1], options, expansion),
                   c = exact_value(q[0], options, expansion);
        if (!a || !b || !c) return std::nullopt;
        uint64_t coefficient_bits = 0;
        for (const auto *g : {&*a, &*b, &*c})
            for (const auto *x : {&g->re, &g->im})
                coefficient_bits = std::max<uint64_t>(coefficient_bits,
                    std::max(fmpz_bits(fmpq_numref(x->value())), fmpz_bits(fmpq_denref(x->value()))));
        // Each component of b² - 4ac sums four products of rationals.
        // Uncancelled denominators use at most 8B bits; scaling/sums add 4.
        const auto discriminant_bits = resources::add(resources::multiply(8, coefficient_bits), 4);
        if (!exact_operands_fit(degree + 1, discriminant_bits, options)) return std::nullopt;
        auto d = product(*b, *b);
        const auto ac = product(*a, *c);
        Rational four_ac;
        fmpq_mul_2exp(four_ac.value(), ac.re.value(), 2);
        fmpq_sub(d.re.value(), d.re.value(), four_ac.value());
        fmpq_mul_2exp(four_ac.value(), ac.im.value(), 2);
        fmpq_sub(d.im.value(), d.im.value(), four_ac.value());
        double_root = fmpq_is_zero(d.re.value()) && fmpq_is_zero(d.im.value());
    }
    std::optional<ExactCubic> cubic;
    std::optional<ExactQuartic> quartic;
    if (degree >= 3) {
        const bool expansion = options.allow_exponent_expansion ||
                               compact_expansion_fits(std::span<const Coefficient>(q, degree + 1),
                                                      options.exact_expansion_allowance());
        std::vector<Gaussian> exact;
        for (size_t k = 0; k <= degree; ++k) {
            auto value = exact_value(q[k], options, expansion);
            if (!value) return std::nullopt;
            exact.push_back(std::move(*value));
        }
        const auto integers = integer_coefficients(exact, options);
        if (!integers) return std::nullopt;
        // The invariants' size before they are formed: E is of degree 6 in the coefficients, its integer factors
        // (for the quartic 27^2, 4 27, 72^2 and their sums) add fewer than 32 bits; the exact storage limit bounds it
        // as it bounds the input.
        constexpr uint64_t invariant_degree = 6, constant_bits = 32;
        uint64_t coefficient_bits = 0;
        for (const auto &c : *integers)
            coefficient_bits = std::max<uint64_t>(coefficient_bits, std::max(fmpz_bits(c.re.value()), fmpz_bits(c.im.value())));
        const auto invariant_bits = resources::add(resources::multiply(invariant_degree, coefficient_bits), constant_bits);
        if (!exact_operands_fit(degree + 1, invariant_bits, options))
            return std::nullopt;
        if (degree == 3) {
            cubic.emplace(std::span<const GaussianInteger, 4>(integers->data(), 4));
        } else {
            quartic.emplace(std::span<const GaussianInteger, 5>(integers->data(), 5));
            if (quartic->declined()) return std::nullopt;
        }
    }
    const bool multiple = double_root || (cubic && !cubic->square_free()) || (quartic && !quartic->square_free());
    std::vector<Root> terminal_cover;
    cps_result result;
    result.degree = static_cast<int64_t>(length - 1);
    result.goal = engine_goal(options);
    const uint64_t requested = std::max<uint64_t>(64, options.bits + 32);
    uint64_t precision = std::min(requested, options.max_work_bits);
    while (precision && !options.cancelled() && result.progress.work_units_this_call < options.work_units) {
        ++result.progress.work_units_this_call;
        ++result.statistics_v2.work_units;
        result.work_bits = result.progress.max_outer_bits = precision;
        // Each precision is one attempt; every one after the first is a raise.
        if (result.statistics.precision_attempts++) ++result.statistics_v2.precision_increases;
        const auto bits = static_cast<slong>(precision);
        ComplexBall a, b, c, discriminant, radical, product_ac, big, small;
        q[0].enclose(c.value(), bits);
        q[1].enclose(b.value(), bits);
        std::vector<Root> roots;
        // A cubic whose forms leave a choice undecided at this precision has no boxes: the precision is raised.
        bool decided = true;
        if (cubic || quartic) {
            auto boxes = cubic ? cubic->roots(bits) : quartic->roots(bits);
            decided = boxes.has_value();
            if (decided) roots = std::move(*boxes);
        } else if (degree == 1) {
            roots.resize(1);
            acb_div(roots[0].box.value(), c.value(), b.value(), bits);
            acb_neg(roots[0].box.value(), roots[0].box.value());
            roots[0].count = roots[0].multiplicity = 1;
        } else if (double_root) {
            q[2].enclose(a.value(), bits);
            roots.resize(1);
            acb_div(roots[0].box.value(), b.value(), a.value(), bits);
            acb_mul_2exp_si(roots[0].box.value(), roots[0].box.value(), -1);
            acb_neg(roots[0].box.value(), roots[0].box.value());
            roots[0].count = roots[0].multiplicity = 2;
        } else {
            q[2].enclose(a.value(), bits);
            acb_mul(discriminant.value(), b.value(), b.value(), bits);
            acb_mul(product_ac.value(), a.value(), c.value(), bits);
            acb_mul_2exp_si(product_ac.value(), product_ac.value(), 2);
            acb_sub(discriminant.value(), discriminant.value(), product_ac.value(), bits);
            if (arb_is_negative(acb_realref(discriminant.value()))) {
                // Either branch of the square root gives the same pair; a left-half-plane discriminant is rotated away
                // from the cut, so that a tiny imaginary uncertainty does not enclose both signs.
                acb_neg(radical.value(), discriminant.value());
                acb_sqrt(radical.value(), radical.value(), bits);
                acb_mul_onei(radical.value(), radical.value());
            } else {
                acb_sqrt(radical.value(), discriminant.value(), bits);
            }
            // Q = -(b +- sqrt D) / 2, the sign giving the larger modulus; then Q / a and c / Q.
            acb_add(big.value(), b.value(), radical.value(), bits);
            acb_sub(small.value(), b.value(), radical.value(), bits);
            Float big_lower, small_lower;
            acb_get_abs_lbound_arf(big_lower.value(), big.value(), bits);
            acb_get_abs_lbound_arf(small_lower.value(), small.value(), bits);
            if (arf_cmp(small_lower.value(), big_lower.value()) > 0) acb_swap(big.value(), small.value());
            acb_neg(big.value(), big.value());
            acb_mul_2exp_si(big.value(), big.value(), -1);
            roots.resize(2);
            acb_div(roots[0].box.value(), big.value(), a.value(), bits);
            acb_div(roots[1].box.value(), c.value(), big.value(), bits);
            roots[0].count = roots[0].multiplicity = 1;
            roots[1].count = roots[1].multiplicity = 1;
        }
        bool valid = decided, goal_met = true;
        for (auto &root : roots) {
            if (!valid || !acb_is_finite(root.box.value())) {
                valid = false;
                break;
            }
            root.certificate = cubic     ? "closed form on the exact input (degree 3), the multiplicities decided exactly"
                               : quartic ? "closed form on the exact input (degree 4), the multiplicities decided exactly"
                                         : "closed form on the exact input (degree 1 or 2), the multiplicity decided exactly";
            if (arb_is_zero(acb_imagref(root.box.value())))
                root.realness = CPS_REAL;
            else if (!arb_contains_zero(acb_imagref(root.box.value())))
                root.realness = CPS_NONREAL;
            root.publish(bits);
            goal_met = goal_met && meets_goal(root, options, bits);
        }
        if (valid && zeros) roots.push_back(exact_zero_point(zeros));
        if (valid && disjoint(roots)) {
            sort_roots(roots);
            result.roots = std::move(roots);
            result.complete = true;
            result.goal_met = goal_met;
            result.squarefree = zeros > 1 || multiple ? CPS_SQUAREFREE_DISPROVED : CPS_SQUAREFREE_PROVED;
        }
        if (valid && !result.complete &&
            result.progress.work_units_this_call >= options.work_units) {
            uint64_t counted = 0;
            for (const auto &root : roots) counted += root.count;
            if (counted == length - 1 && !disjoint(roots))
                terminal_cover = std::move(roots);
        }
        if (result.goal_met) {
            result.termination = CPS_SUCCESS;
            result.diagnostics = std::string("closed form for degree ") + std::to_string(degree) +
                                 (zeros ? " after the exact zero roots" : "") +
                                 (double_root ? ", an exact double root" : (cubic || quartic) && multiple ? ", exact multiple roots" : "") +
                                 "; " + std::to_string(result.progress.work_units_this_call) + " work units; stop=0";
            return result;
        }
        if (precision >= options.max_work_bits) break;
        precision = precision > options.max_work_bits / 2 ? options.max_work_bits : 2 * precision;
    }
    if (!result.progress.work_units_this_call) return std::nullopt;
    const bool cancelled = options.cancelled();
    const bool spent = result.progress.work_units_this_call >= options.work_units;
    if (!cancelled && spent && !terminal_cover.empty()) {
        // Each finite branch encloses one member of the complete exact root
        // multiset, with its proved multiplicity. Disjoint unions therefore
        // contain exactly the sum of their assigned branch counts.
        merge_overlapping(terminal_cover, static_cast<slong>(result.work_bits));
        sort_roots(terminal_cover);
        result.goal_met = true;
        for (auto &root : terminal_cover) {
            if (root.certificate == "complete-factor-cover-union")
                root.certificate = "closed-form complete branch-cover union";
            result.goal_met = result.goal_met &&
                              meets_goal(root, options, static_cast<slong>(result.work_bits));
        }
        result.roots = std::move(terminal_cover);
        result.complete = true;
        result.squarefree = zeros > 1 || multiple ? CPS_SQUAREFREE_DISPROVED : CPS_SQUAREFREE_PROVED;
        if (result.goal_met) {
            result.termination = CPS_SUCCESS;
            result.diagnostics = "closed-form complete branch cover at the final work unit; " +
                                 std::to_string(result.progress.work_units_this_call) + " work units; stop=0";
            return result;
        }
    }
    result.termination = cancelled ? CPS_CANCELLED : spent ? CPS_WORK_LIMIT : CPS_PRECISION_LIMIT;
    result.progress.stop_reason = cancelled ? CPS_STOP_CANCELLED : spent ? CPS_STOP_WORK_BUDGET : CPS_STOP_PRECISION;
    result.diagnostics = "closed-form attempt stopped before the goal; " +
                         std::to_string(result.progress.work_units_this_call) + " work units; stop=" +
                         std::to_string(result.progress.stop_reason);
    return result;
}
} // namespace certpolysolve::detail
