// Piecewise's admission (admit_piecewise): whether the piecewise method takes a query, decided before any root work
// from the query, the input's structure and a work estimate; the verdict and its reason are a PiecewiseAdmission.
//
// The decisions, in their order, under a PiecewiseAdmissionPolicy. structural_decision refuses a real-only query or a
// goal other than accuracy (the default) and approximation, an unusable precision or a cancelled call, a constant or
// zero input, a degree below the policy's minimum_degree, inexact coefficients, a power composite (support_stride above
// one) and, under a nonzero minimum_degree, a support of at most a quarter of the degree; admission_by_structure makes
// only these and returns nothing when they all pass. admit_impl then encloses a DensePolynomial's coefficients at
// 96 bits, or admit_terms a cps_polynomial's terms (without the square-free proof), and admit_enclosed sums
// estimate_piecewise_models over the model ladder's rungs, 16, 32, ... up to bit_ceil(2 bits) model bits, refusing an
// envelope bend above the policy's cap, unavailable geometry, or an estimate that reaches the global-sweep cap while
// its rings are walked (a monotone test, so stopping early changes no verdict). The total must then pass the work-ratio
// cap on work / (degree x support) and the global-sweep cap on work / (degree x (degree + support)), both in exact
// integer quotients; under require_attempt_fits_call, twice the estimated models must be below the call's work units;
// and the modular square-free proof comes last (prove_squarefree: nonzero_part_squarefree_modular).
//
// State: none. Each call is a function of its arguments; its own time goes to PiecewiseAdmission::elapsed_ns.
// Entry points: admit_piecewise and admission_by_structure, from select_auto_route (auto_selection.cpp), and
// nonzero_part_squarefree_modular, from prepare_squarefree_input (exact_preprocessing.cpp).
// Not here: the estimate's ring walk (piecewise_models.cpp), the modular gcd (squarefree_modular in
// modular_squarefree.cpp), the input's structure (input_structure.hpp), the selector's policy (auto_selection.cpp).
// Design notes: docs/mathematics.md, "Routing: the selector" (the caps as Auto applies them);
// docs/piecewise-selection.md (the rule's structure and its modular proof; superseded otherwise, its status line);
// docs/selector-campaign-20260928.md (the sweep cap); docs/same-routes-every-platform-note-20261021.md (the caps in
// exact integers).
#include "methods/piecewise/piecewise_admission.hpp"
#include "numeric/exact/modular_squarefree.hpp"
#include "methods/piecewise/piecewise_estimate.hpp"
#include "core/telemetry.hpp"
#include <algorithm>
#include <bit>
#include <cmath>
#include <numeric>

namespace certpolysolve::detail {
// The admission's first decisions, which need only the query and the input's structure, in their order (the
// cancellation poll among them). One function, because the selector makes them from the input object's terms without
// its coefficient list (admission_by_structure, C40's second step). True when one of them decided, its reason and the
// fields it had set in `result`.
static bool structural_decision(PiecewiseAdmission &result, const cps_options &options,
                                const PiecewiseAdmissionPolicy &policy, const InputStructure &metadata) {
    if (options.real_only || (engine_goal(options) != CPS_ROOT_ACCURACY && engine_goal(options) != CPS_ROOT_APPROXIMATION))
        return true;
    if (options.bits < 2 || options.bits > options.max_work_bits ||
        options.bits > (static_cast<uint64_t>(std::numeric_limits<slong>::max()) - 64) / 8 ||
        options.cancelled()) {
        result.reason = "precision_or_cancellation";
        return true;
    }
    const auto start = metadata.start, length = metadata.length;
    result.zero_order = start;
    if (length <= start + 1) {
        result.reason = "constant_or_zero";
        return true;
    }
    result.degree = length - start - 1;
    result.support = metadata.support;
    result.support_stride = metadata.support_stride;
    if (result.degree < policy.minimum_degree) {
        result.reason = "small_degree";
        return true;
    }
    if (!metadata.exact) {
        result.reason = "uncertainty_or_cancellation";
        return true;
    }
    if (result.support_stride > 1) {
        result.reason = "power_composition";
        return true;
    }
    if (policy.minimum_degree && result.support <= result.degree / 4) {
        result.reason = "sparse_coefficients";
        return true;
    }
    return false;
}

// The admission after its coefficients are enclosed: the precision rule, the estimate over the model rungs
// (`estimator(model_bits, bend, beyond)`), the work rules and the square-free proof (`squarefree()`, asked only under
// policy.prove_squarefree). Shared by the coefficient list and the input object's terms (C40's second step).
template <class Estimator, class Squarefree>
static PiecewiseAdmission admit_enclosed(PiecewiseAdmission result, const cps_options &options,
                                         const PiecewiseAdmissionPolicy &policy, Estimator &&estimator,
                                         Squarefree &&squarefree) {
    if (options.bits > resources::precision_max / 4) {
        result.reason = "model_precision_not_representable";
        return result;
    }
    result.model_bits = std::max<uint64_t>(16, std::bit_ceil(2 * options.bits));
    uint64_t work = 0;
    // The caller's cap in global sweeps (see PiecewiseAdmissionPolicy), exact: total >= sweeps x degree x (degree +
    // support) is floor(floor(total / degree) / (degree + support)) >= sweeps, for these positive integers. The final
    // check uses it too, and it never turns false as the work grows, so stopping the estimate early cannot change the
    // answer.
    const uint64_t degree = result.degree, support = result.support;
    const auto above_global_sweeps = [&](uint64_t total) {
        return total / degree / (degree + support) >= *policy.maximum_global_sweeps;
    };
    const bool sweep_cap = policy.maximum_global_sweeps.has_value();
    const std::function<bool(uint64_t)> beyond =
        sweep_cap ? [&](uint64_t level) { return above_global_sweeps(resources::add(work, level)); }
                  : std::function<bool(uint64_t)>{};
    try {
        for (uint64_t model_bits = 16; model_bits <= result.model_bits; model_bits *= 2) {
            if (options.cancelled()) {
                result.reason = "cancelled";
                return result;
            }
            const auto estimate = estimator(static_cast<slong>(model_bits), policy.maximum_envelope_bend, beyond);
            result.envelope_bend = estimate.envelope_bend;
            if (estimate.bend_exceeded) {
                result.reason = "broad_radial_geometry";
                return result;
            }
            if (estimate.truncated) {
                result.reason = "estimate_above_global_sweeps";
                return result;
            }
            if (!estimate.available) {
                result.reason = "unavailable_geometry";
                return result;
            }
            if (estimate.models > (std::numeric_limits<uint64_t>::max() / 2) - result.estimated_models) {
                result.reason = "model_budget_overflow";
                return result;
            }
            result.estimated_models += estimate.models;
            work = resources::add(work, resources::add(estimate.construction_work, estimate.local_work));
        }
    } catch (const std::length_error &) {
        result.reason = "unavailable_geometry";
        return result;
    } catch (const std::invalid_argument &) {
        result.reason = "unavailable_geometry";
        return result;
    }
    result.work = work;
    result.work_ratio = static_cast<double>(work) / (static_cast<double>(degree) * static_cast<double>(support));
    // The cap exact: work / (degree x support) > ratio is ceil(ceil(work / degree) / support) > ratio.
    const auto ceil_quotient = [](uint64_t a, uint64_t b) { return a / b + (a % b != 0); };
    if (policy.maximum_work_ratio && ceil_quotient(ceil_quotient(work, degree), support) > *policy.maximum_work_ratio) {
        result.reason = "estimated_work";
        return result;
    }
    if (sweep_cap && above_global_sweeps(work)) {
        result.reason = "estimate_above_global_sweeps";
        return result;
    }
    result.estimated_attempt_work = 2 * result.estimated_models;
    // Leave at least one charged step for fallback. This is an admission
    // estimate, never a promise of completion or permission to exceed quota.
    if (policy.require_attempt_fits_call && result.estimated_attempt_work >= options.work_units) {
        result.reason = "remaining_call_work";
        return result;
    }
    result.squarefree_proved = policy.prove_squarefree && squarefree();
    if (policy.prove_squarefree && !result.squarefree_proved) {
        result.reason = "squarefree_not_proved";
        return result;
    }
    result.selected = true;
    result.reason = "bounded_local_work";
    return result;
}

static PiecewiseAdmission admit_impl(const DensePolynomial &input, const cps_options &options,
                                     PiecewiseAdmissionPolicy policy, const InputStructure &metadata) {
    PiecewiseAdmission result;
    if (structural_decision(result, options, policy, metadata))
        return result;
    const auto start = metadata.start, length = metadata.length;
    ComplexPolynomial polynomial;
    ComplexBall coefficient;
    for (size_t j = start; j < length; ++j) {
        const auto &value = input.coefficients[j];
        if (!value.exact() || options.cancelled()) {
            result.reason = "uncertainty_or_cancellation";
            return result;
        }
        value.enclose(coefficient.value(), 96);
        acb_poly_set_coeff_acb(polynomial.value(), static_cast<slong>(j - start), coefficient.value());
    }
    return admit_enclosed(
        std::move(result), options, policy,
        [&](slong bits, std::optional<uint64_t> bend, const std::function<bool(uint64_t)> &beyond) {
            return estimate_piecewise_models(polynomial.value(), bits, bend, beyond);
        },
        [&] { return nonzero_part_squarefree_modular(input, options); });
}

// The same on the input object's terms, without the square-free proof (policy.prove_squarefree false): the real-only
// route's selection for its complex query. The terms are enclosed as the coefficient loop encloses them, polling the
// cancellation once a term where the loop polls once a coefficient.
static PiecewiseAdmission admit_terms(const cps_polynomial &input, const cps_options &options,
                                      PiecewiseAdmissionPolicy policy, const InputStructure &metadata) {
    if (policy.prove_squarefree)
        throw std::logic_error("the terms admission makes no square-free proof");
    PiecewiseAdmission result;
    if (structural_decision(result, options, policy, metadata))
        return result;
    const auto start = static_cast<uint64_t>(metadata.start);
    std::vector<std::pair<slong, ComplexBall>> terms;
    terms.reserve(input.terms.size());
    for (const auto &[exponent, value] : input.terms) {
        if (!value.exact() || options.cancelled()) {
            result.reason = "uncertainty_or_cancellation";
            return result;
        }
        terms.emplace_back(static_cast<slong>(exponent - start), ComplexBall{});
        value.enclose(terms.back().second.value(), 96);
    }
    return admit_enclosed(
        std::move(result), options, policy,
        [&](slong bits, std::optional<uint64_t> bend, const std::function<bool(uint64_t)> &beyond) {
            return estimate_piecewise_models(std::span<const std::pair<slong, ComplexBall>>(terms), bits, bend, beyond);
        },
        [] { return false; });
}

PiecewiseAdmission admit_piecewise(const cps_polynomial &input, const cps_options &options,
                                   PiecewiseAdmissionPolicy policy) {
    PiecewiseAdmission result;
    uint64_t elapsed = 0;
    {
        StageTimer timer(elapsed);
        result = admit_terms(input, options, policy, input_structure_of_terms(input.terms));
    }
    result.elapsed_ns = elapsed;
    return result;
}

std::optional<PiecewiseAdmission> admission_by_structure(const InputStructure &metadata, const cps_options &options,
                                                         PiecewiseAdmissionPolicy policy) {
    PiecewiseAdmission result;
    uint64_t elapsed = 0;
    bool decided = false;
    {
        StageTimer timer(elapsed);
        decided = structural_decision(result, options, policy, metadata);
    }
    if (!decided)
        return std::nullopt;
    result.elapsed_ns = elapsed;
    return result;
}

PiecewiseAdmission admit_piecewise(const DensePolynomial &input, const cps_options &options,
                                   PiecewiseAdmissionPolicy policy) {
    PiecewiseAdmission result;
    uint64_t elapsed = 0;
    {
        StageTimer timer(elapsed);
        result = admit_impl(input, options, policy, input_structure(input.coefficients));
    }
    result.elapsed_ns = elapsed;
    return result;
}
} // namespace certpolysolve::detail
