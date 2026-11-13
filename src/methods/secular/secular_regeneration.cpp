// The secular regeneration: the residues of the secular model at its nodes, each to the accuracy its node's iteration
// needs (regenerate_residues): exact nodes and their products of differences (product_derivative), the disk values
// of a polynomial source (DiskValues), a residue's adequacy (adequate, adequate_accuracy_bits), the budget-limited
// nodes a step leaves unsteered (unsteered_nodes). Part of the secular refinement (secular_refinement.cpp), split
// from it verbatim in W8's S5 (docs/w8-plan-20261021.md); state: none, every function works on its arguments.
// Entry points: regenerate_residues and unsteered_nodes (SecularRefinement::step).
#include "methods/secular/secular_refinement.hpp"
#include "core/query.hpp"
#include "methods/secular/secular_internal.hpp"
#include "numeric/evaluate/secular_disk_values.hpp"
#include "methods/secular/secular_geometry.hpp"
#include "core/dovetail.hpp"
#include "numeric/evaluate/node_products.hpp"
#include "core/telemetry.hpp"
#include <algorithm>
#include <cmath>
#include <complex>
#include <functional>
#include <limits>
#include <numeric>
#include <optional>

namespace certpolysolve::detail {
namespace {
bool product_derivative(acb_t product, mag_t separation, std::span<const ComplexBall> nodes, size_t i,
                        slong precision, const cps_options &options) {
    std::vector<ComplexBall> terms;
    terms.reserve(nodes.size());
    mag_inf(separation);
    Magnitude distance;
    for (size_t j = 0; j < nodes.size(); ++j)
        if (i != j) {
            if ((j & 63u) == 0 && options.cancelled())
                return false;
            terms.emplace_back();
            acb_sub(terms.back().value(), nodes[i].value(), nodes[j].value(), precision);
            if (acb_contains_zero(terms.back().value()))
                return false;
            acb_get_mag_lower(distance.value(), terms.back().value());
            mag_min(separation, separation, distance.value());
        }
    if (terms.empty()) {
        acb_one(product);
        return true;
    }
    size_t length = terms.size();
    while (length > 1) {
        size_t next = 0;
        for (size_t j = 0; j < length; j += 2) {
            if (j + 1 < length)
                acb_mul(terms[next].value(), terms[j].value(), terms[j + 1].value(), precision);
            else
                acb_set(terms[next].value(), terms[j].value());
            ++next;
        }
        length = next;
    }
    acb_set(product, terms[0].value());
    return !acb_contains_zero(product);
}
} // namespace
// The relative accuracy that makes a residue adequate whatever its size (adequate()'s first clause): eight bits
// beyond the node's iteration bits, and for an iteration model no more than sixteen beyond the goal's.
slong adequate_accuracy_bits(slong bits, uint64_t goal, bool iteration_model) {
    const slong own = resources::precision_sum(bits, 8);
    return iteration_model ? std::min<slong>(own, resources::precision_sum(goal, 16)) : own;
}
namespace {
// A residue may be a rigorous enclosure yet have an uninformative midpoint.
// The regenerated candidate model needs reliable coefficients; its requested
// accuracy is independent of the precision needed to resolve cancellation in
// the original polynomial. The resource-ceiling exit below still retains wide
// valid enclosures for certification without promising candidate convergence.
bool adequate(const acb_t residue, const acb_t node, const mag_t separation, size_t n, slong bits,
              uint64_t goal, bool iteration_model = false, const acb_struct *origin_offset = nullptr) {
    if (iteration_model) {
        if (acb_is_zero(residue) || acb_rel_accuracy_bits(residue) >= adequate_accuracy_bits(bits, goal, true))
            return true;
        // A residue whose whole enclosure, times n (the inclusion radius), lies below the node's resolution
        // and below its separation's moves nothing: the step it gives its node is below what the node's bits
        // represent (for a centred source, at the node's distance from the original origin when that is the
        // nearer), and its term in the model is below the model's resolution at every other node (they are
        // at least the separation away). A node exactly on a root has residue exactly zero, whose enclosure
        // gains no relative accuracy short of exact arithmetic: with coefficients that are not dyadic, the
        // enclosure is rounding noise at every precision, and without this clause the precision would double
        // to the cap. History: docs/code-history-20261021.md, secular_refinement.cpp 3.
        Magnitude bound, resolution;
        acb_get_mag(bound.value(), residue);
        mag_mul_ui(bound.value(), bound.value(), static_cast<ulong>(n));
        acb_get_mag(resolution.value(), node);
        mag_min(resolution.value(), resolution.value(), separation);
        if (origin_offset) {
            ComplexBall original;
            Magnitude from_origin;
            acb_add(original.value(), node, origin_offset, ARF_PREC_EXACT);
            acb_get_mag(from_origin.value(), original.value());
            mag_min(resolution.value(), resolution.value(), from_origin.value());
        }
        mag_mul_2exp_si(resolution.value(), resolution.value(), -bits);
        return !mag_is_zero(resolution.value()) && mag_cmp(bound.value(), resolution.value()) <= 0;
    }
    if (acb_rel_accuracy_bits(residue) >= adequate_accuracy_bits(bits, goal, false))
        return true;
    Magnitude error, target, local;
    mag_add(error.value(), arb_radref(acb_realref(residue)), arb_radref(acb_imagref(residue)));
    mag_mul_ui(error.value(), error.value(), static_cast<ulong>(n));
    acb_get_mag(target.value(), node);
    mag_mul_2exp_si(target.value(), target.value(), -static_cast<slong>(goal + 16));
    mag_mul_2exp_si(local.value(), separation, -5);
    mag_min(target.value(), target.value(), local.value());
    return mag_cmp(error.value(), target.value()) <= 0;
}
} // namespace
// The share of the nodes that may have moved for an unmoved node's residue to be updated in ball arithmetic
// (its factors one by one at the evaluation's precision) and not rebuilt: a quarter.
inline constexpr size_t update_quarter = 4;
bool regenerate_residues(std::vector<ComplexBall> &out, MonicSource &source,
                         std::span<const ComplexBall> nodes, std::span<const ComplexBall> old_nodes,
                         std::span<const ComplexBall> old_residues, std::span<const slong> iteration_bits,
                         const cps_options &options, certpolysolve::detail::CandidateStatistics &statistics, uint64_t &maximum_precision,
                         const EnginePolicy *policy, bool partial, std::span<slong> source_precision_hints,
                         RegenerationWork *work, size_t *budget_limited, std::vector<uint8_t> *input_limited,
                         const acb_struct *origin_offset, std::span<const uint8_t> estimated) {
    StageTimer timer(statistics.regeneration_ns);
    size_t limited = 0;
    const bool inexact_source = source.metadata().input == InputSemantics::CoefficientEnclosures;
    if (input_limited) input_limited->assign(nodes.size(), 0);
    TraceScope trace(policy, "regeneration", maximum_precision, nodes.size(), &statistics);
    const size_t n = nodes.size();
    const bool iteration_model = conditioned_secular(source, policy) ||
        (policy && policy->model_stopping == SecularModelPolicy::Continuation);
    const RegenerationPolicy experiment = policy ? policy->regeneration : RegenerationPolicy{};
    const bool remember = experiment.remember_precision && source_precision_hints.size() == n;
    std::vector<slong> proposed_hints;
    if (remember) proposed_hints.assign(source_precision_hints.begin(), source_precision_hints.end());
    auto *counts = policy ? policy->study_statistics.get() : nullptr;
    if (n != static_cast<size_t>(source.degree()) || iteration_bits.size() != n || !exact_nodes(nodes))
        return false;
    std::vector<size_t> changed;
    if (old_nodes.size() == n && old_residues.size() == n) {
        for (size_t i = 0; i < n; ++i)
            if (!acb_equal(nodes[i].value(), old_nodes[i].value()))
                changed.push_back(i);
    } else
        partial = false;
    // The products of nodes of at most two limbs a part (node_products.hpp), where there is more than one node
    // and every one is; and the nodes before, where an unmoved node's update takes its factors there too.
    std::optional<ShortNodeProducts> node_products, old_node_products;
    if (experiment.integer_node_products && experiment.reuse_products) {
        node_products.emplace(nodes);
        if (!node_products->usable())
            node_products.reset();
        else if (partial) {
            old_node_products.emplace(old_nodes);
            if (!old_node_products->usable()) old_node_products.reset();
        }
    }
    // An unmoved node's residue is updated by the moved nodes' factors while that costs less than rebuilding
    // it. With the factors in ball arithmetic: where at most a quarter of the nodes moved. In integer mantissas
    // an update is two factors a moved node, against a factor a node and an evaluation of the source for the
    // rebuild: less whatever moved, so every unmoved node is updated.
    const size_t update_limit = old_node_products ? n - 1 : n / update_quarter;
    partial = partial && !changed.empty() && changed.size() <= update_limit;
    std::vector<ComplexBall> proposed(n);
    ComplexBall product, value, derivative, numerator, denominator;
#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
    ComplexBall alternative; // the exact residue at the budget's (research build)
#endif
    Magnitude separation;
    // The meter's operand sizes (dovetail.hpp, horner_weight): an evaluation
    // multiplies by its own node; a node product by differences of nodes,
    // taken as the longest node.
    const auto node_bits = [&](size_t i) {
        return static_cast<int64_t>(std::max(arf_bits(arb_midref(acb_realref(nodes[i].value()))),
                                             arf_bits(arb_midref(acb_imagref(nodes[i].value())))));
    };
    int64_t longest_node = 0;
    if (work)
        for (size_t i = 0; i < n; ++i)
            longest_node = std::max(longest_node, node_bits(i));
#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
    auto *polynomial = dynamic_cast<PolynomialSource *>(&source);
    const bool exact_at_budget = polynomial && polynomial->exact_available();
#endif
    DiskValues disk(source, experiment.limb_horner);
    for (size_t i = 0; i < n; ++i) {
        if (options.cancelled())
            return false;
        slong prec = std::min<slong>(static_cast<slong>(options.max_work_bits),
                                     resources::precision_sum(iteration_bits[i], regeneration_guard_bits));
        const slong ceiling = static_cast<slong>(options.max_work_bits);
        const slong product_floor = prec;
        if (remember && source_precision_hints[i] > prec) {
            prec = std::min(ceiling, source_precision_hints[i]);
            if (counts) ++counts->source_hint_starts;
        }
        const slong residue_bits = iteration_model
            ? std::min<slong>(resources::precision_sum(iteration_bits[i], 8), resources::precision_sum(options.bits, 16))
            : resources::precision_sum(iteration_bits[i], 8);
        // The product depends only on this node geometry. Keep its enclosure
        // during source cancellation retries; no cache survives a model change.
        slong product_precision = product_floor;
        bool cached_product = false;
        bool updating = partial && acb_equal(nodes[i].value(), old_nodes[i].value());
        Magnitude previous_error; // of the last full evaluation on an inexact source
        mag_inf(previous_error.value());
        slong previous_precision = 0; // the precision of that evaluation
        // An estimated start that falls short takes its measured shortfall once (below).
        bool shortfall_step = estimated.size() == n && estimated[i] && remember && source_precision_hints[i] == prec;
        for (;;) {
            maximum_precision = std::max(maximum_precision, static_cast<uint64_t>(prec));
            bool valid = true, product_ready = false;
            if (updating) {
                acb_set(proposed[i].value(), old_residues[i].value());
                mag_inf(separation.value());
                if (work)
                    work->operations += static_cast<double>(changed.size()) * horner_weight(prec, longest_node);
                // The residue changes by prod (b_i - old b_j) / (b_i - b_j) over the moved nodes j. Between
                // short nodes the two products come from integer mantissas (node_products.hpp), where their
                // bound serves the residue as it must for a whole product; the factors one by one in ball
                // arithmetic otherwise, at the evaluation's precision. History:
                // docs/code-history-20261021.md, secular_refinement.cpp 11.
                if (old_node_products && ShortNodeProducts::accuracy_bits(changed.size() + 1) >=
                                               resources::precision_sum(residue_bits, 8)) {
                    if (work) work->short_product_factors += 2 * changed.size();
                    valid = node_products->product_over(numerator.value(), i, changed, *old_node_products) &&
                            node_products->product_over(denominator.value(), i, changed, *node_products);
                    if (valid) {
                        if (work) work->service_operations += operation_weight(product_floor) + operation_weight(prec);
                        acb_div(numerator.value(), numerator.value(), denominator.value(), product_floor);
                        acb_mul(proposed[i].value(), proposed[i].value(), numerator.value(), prec);
                    }
                    if (counts) ++counts->integer_node_updates;
                } else if (changed.size() > n / update_quarter) {
                    // More moved than a ball update is worth: rebuild.
                    updating = false;
                    continue;
                } else {
                    if (work) work->service_operations += static_cast<double>(changed.size()) * horner_weight(prec, longest_node);
                    for (auto j : changed) {
                        acb_sub(numerator.value(), nodes[i].value(), old_nodes[j].value(), prec);
                        acb_sub(denominator.value(), nodes[i].value(), nodes[j].value(), prec);
                        if (acb_contains_zero(denominator.value())) {
                            valid = false;
                            break;
                        }
                        acb_div(numerator.value(), numerator.value(), denominator.value(), prec);
                        acb_mul(proposed[i].value(), proposed[i].value(), numerator.value(), prec);
                    }
                }
            } else {
                if (!experiment.reuse_products || !cached_product) {
                    if (!experiment.reuse_products)
                        product_precision = prec;
                    // The integer product serves when its accuracy is what the retry below asks of a product
                    // (the residue's bits and eight more); a retry at a higher precision is a ball product.
                    const bool integers = node_products && product_precision == product_floor &&
                        ShortNodeProducts::accuracy_bits(n) >= resources::precision_sum(residue_bits, 8);
                    product_ready = integers
                        ? node_products->product(product.value(), separation.value(), i)
                        : product_derivative(product.value(), separation.value(), nodes, i, product_precision,
                                             options);
                    if (integers && counts) ++counts->integer_node_products;
                    if (work)
                        work->operations += static_cast<double>(n) * horner_weight(product_precision, longest_node);
                    if (work) {
                        if (integers) work->short_product_factors += n - 1;
                        else work->service_operations += static_cast<double>(n) * horner_weight(product_precision, longest_node);
                    }
                    cached_product = product_ready;
                    maximum_precision = std::max(maximum_precision, static_cast<uint64_t>(product_precision));
                    if (counts) {
                        ++counts->product_builds;
                        counts->product_peak_bits = std::max(counts->product_peak_bits,
                                                             static_cast<uint64_t>(product_precision));
                    }
                } else {
                    product_ready = true;
                    if (counts) ++counts->product_reuses;
                }
                valid = product_ready;
                if (valid) {
                    StageTimer evaluation_timer(statistics.evaluation_ns);
                    ++statistics.evaluations;
                    statistics.evaluation_terms += source.terms();
                    if (work) {
                        work->operations += static_cast<double>(source.terms()) * horner_weight(prec, node_bits(i));
                        work->source_terms += source.terms();
                    }
                    if (counts) {
                        counts->source_peak_bits = std::max(counts->source_peak_bits, static_cast<uint64_t>(prec));
                        if (experiment.value_only) ++counts->source_values;
                        else ++counts->source_joints;
                    }
                    if (auto prepared = disk.prepare(prec, nodes[i], experiment.value_only)) {
                        if (work) {
                            if (prepared->uses_limb_horner()) work->limb_horner_work += prepared->research_limb_work(source.terms());
                            else work->service_operations += static_cast<double>(source.terms()) * horner_weight(prec, node_bits(i));
                        }
                        if (counts && prepared->uses_limb_horner()) ++counts->limb_horner_evaluations;
                        valid = prepared->evaluate(value.value(), derivative.value());
                    } else if (experiment.value_only) {
                        if (work) work->service_operations += static_cast<double>(source.terms()) * horner_weight(prec, node_bits(i));
                        auto evaluated = source.value(nodes[i], prec);
                        valid = evaluated.status == EvaluationStatus::Enclosed;
                        value = std::move(evaluated.value);
                    } else {
                        if (work) work->service_operations += static_cast<double>(source.terms()) * horner_weight(prec, node_bits(i));
                        valid = source.evaluate(value.value(), derivative.value(), nodes[i].value(), prec);
                    }
                }
                if (valid) {
                    acb_div(proposed[i].value(), value.value(), product.value(), prec);
                    acb_neg(proposed[i].value(), proposed[i].value());
                }
            }
            valid = valid && acb_is_finite(proposed[i].value());
            if (valid && adequate(proposed[i].value(), nodes[i].value(), separation.value(), n, iteration_bits[i],
                                  options.bits, iteration_model, origin_offset)) {
                if (remember && !updating) proposed_hints[i] = prec;
                break;
            }
            if (updating) {
                updating = false;
                continue;
            }
            if (valid && inexact_source) {
                // The input's part of the error: a rung that failed to halve
                // the bound (see the declaration). Accepted as the best the
                // input allows, the node marked for the growth, which raises
                // no precision for it.
                Magnitude error;
                mag_add(error.value(), arb_radref(acb_realref(proposed[i].value())),
                        arb_radref(acb_imagref(proposed[i].value())));
                if (!mag_is_inf(previous_error.value())) {
                    Magnitude half;
                    mag_mul_2exp_si(half.value(), previous_error.value(), -1);
                    if (mag_cmp(error.value(), half.value()) > 0) {
                        if (input_limited) (*input_limited)[i] = 1;
                        // The input's part already dominated at the lower
                        // rung: that is the precision to start from next
                        // time. (Remembering this rung instead would raise
                        // the start by one rung per regeneration, the same
                        // climb to the budget one step at a time.)
                        if (remember) proposed_hints[i] = previous_precision;
                        break;
                    }
                }
                mag_set(previous_error.value(), error.value());
                previous_precision = prec;
            }
            if (options.cancelled())
                return false;
            if (experiment.reuse_products && product_precision < ceiling &&
                (!product_ready || acb_rel_accuracy_bits(product.value()) < resources::precision_sum(residue_bits, 8))) {
                product_precision = std::min(ceiling, resources::precision_product(product_precision, 2));
                cached_product = false;
                if (counts) ++counts->product_precision_retries;
                continue;
            }
            if (prec == static_cast<slong>(options.max_work_bits)) {
                // The research build may first take the exact residue at the budget (exact_residues_at_budget).
#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
                // The ladder has reached the precision budget with a residue
                // the plain evaluation leaves inadequate. The budget bounds the
                // working precision, not the exact source: its value at the
                // exact node, computed in Gaussian-integer arithmetic and
                // rounded once to the budget's precision, has no cancellation
                // beyond that rounding, at about degree x (coefficient bits +
                // degree x node bits) bit operations, charged to the exact
                // budgets (exact_fits). Without it a residue whose evaluation
                // cancels beyond the budget has no correct bit, where the exact
                // value has most of the budget's bits (the old default route's
                // exact local model, an exact Taylor shift per coarse centre,
                // resolves such cancellation at a cost quadratic in the
                // degree). Tried where it can serve: when the demanded accuracy
                // is within the precision (no evaluation at prec bits meets a
                // demand above them), or when the plain residue has no correct
                // bit (any correct bit then steers the node). History:
                // docs/code-history-20261021.md, secular_refinement.cpp 12.
                double exact_operations = 0;
                if (exact_at_budget && product_ready &&
                    (!valid || residue_bits <= prec || residue_unresolved(proposed[i].value())) &&
                    polynomial->exact_fits(nodes[i].value(), options, exact_operations)) {
                    if (options.cancelled()) return false;
                    if (counts) ++counts->exact_retries;
                    TraceScope retry(policy, "exact_regeneration", static_cast<uint64_t>(prec), 1, &statistics);
                    bool enclosed;
                    {
                        StageTimer evaluation_timer(statistics.evaluation_ns);
                        ++statistics.evaluations;
                        statistics.evaluation_terms += source.terms();
                        if (work) {
                            work->operations += exact_operations;
                            work->service_operations += exact_operations;
                            work->source_terms += source.terms();
                        }
                        enclosed = polynomial->evaluate_exact(value.value(), derivative.value(), nodes[i].value(),
                                                              prec, !experiment.value_only);
                    }
                    if (enclosed) {
                        acb_div(alternative.value(), value.value(), product.value(), prec);
                        acb_neg(alternative.value(), alternative.value());
                        if (acb_is_finite(alternative.value())) {
                            // Both are enclosures of the same residue; the
                            // tighter one is kept whether or not it is
                            // adequate (it steers the node where a residue
                            // with no correct bit would settle it).
                            Magnitude plain_error, exact_error;
                            mag_add(plain_error.value(), arb_radref(acb_realref(proposed[i].value())),
                                    arb_radref(acb_imagref(proposed[i].value())));
                            mag_add(exact_error.value(), arb_radref(acb_realref(alternative.value())),
                                    arb_radref(acb_imagref(alternative.value())));
                            if (!valid || mag_cmp(exact_error.value(), plain_error.value()) < 0) {
                                acb_swap(proposed[i].value(), alternative.value());
                                valid = true;
                                if (counts) ++counts->exact_accepted;
                            }
                            if (adequate(proposed[i].value(), nodes[i].value(), separation.value(), n, iteration_bits[i],
                                         options.bits, iteration_model)) {
                                if (remember) proposed_hints[i] = prec;
                                break;
                            }
                        }
                    }
                }
#endif
                if (!valid)
                    return false;
                ++limited;
                break; // A wide valid enclosure still proves a conservative cover.
            }
            if (counts) ++counts->source_precision_retries;
            if (shortfall_step && valid) {
                // The start was an estimate (a restart's: the precision its component's model showed, and the
                // accuracy wanted). The enclosure's relative accuracy says by how many bits it fell short: its
                // radius is the evaluation's bound, which every bit of precision halves. Those bits and the
                // guard, in whole limbs, where that is less than the ladder's doubling. History:
                // docs/code-history-20261021.md, secular_refinement.cpp 13.
                shortfall_step = false;
                const slong accuracy = acb_rel_accuracy_bits(proposed[i].value());
                const slong wanted = adequate_accuracy_bits(iteration_bits[i], options.bits, iteration_model);
                if (accuracy < wanted && accuracy > -prec) {
                    const slong more = resources::precision_sum(wanted - accuracy, regeneration_guard_bits);
                    const slong stepped = resources::whole_limbs(resources::precision_sum(prec, more));
                    if (stepped < resources::precision_product(prec, 2)) {
                        prec = std::min(ceiling, stepped);
                        continue;
                    }
                }
            }
            prec = std::min(ceiling, resources::precision_product(prec, 2));
        }
    }
    out = std::move(proposed);
    if (budget_limited) *budget_limited = limited;
    if (remember) std::copy(proposed_hints.begin(), proposed_hints.end(), source_precision_hints.begin());
    ++statistics.secular_regenerations;
    trace.progress(maximum_precision, n);
    return true;
}
std::vector<uint8_t> unsteered_nodes(std::span<const ComplexBall> nodes, std::span<const ComplexBall> residues,
                                     std::span<const uint8_t> inactive) {
    // One word of working precision: the question is whether the step keeps a
    // single correct bit, and the radii of the residue enclosures decide it.
    // The nodes are exact, so their differences lose nothing to cancellation.
    const slong bits = FLINT_BITS;
    const size_t n = nodes.size();
    std::vector<uint8_t> unsteered(n, 0);
    ComplexBall difference, term, denominator, step;
    for (size_t k = 0; k < n; ++k) {
        if (inactive[k])
            continue;
        acb_zero(denominator.value());
        for (size_t j = 0; j < n; ++j)
            if (j != k) {
                acb_sub(difference.value(), nodes[k].value(), nodes[j].value(), bits);
                acb_add(term.value(), residues[j].value(), residues[k].value(), bits);
                acb_div(term.value(), term.value(), difference.value(), bits);
                acb_add(denominator.value(), denominator.value(), term.value(), bits);
            }
        acb_sub_ui(denominator.value(), denominator.value(), 1, bits);
        acb_div(step.value(), residues[k].value(), denominator.value(), bits);
        unsteered[k] = !acb_is_zero(step.value()) &&
                       (!acb_is_finite(step.value()) || acb_rel_accuracy_bits(step.value()) <= 0);
    }
    return unsteered;
}
} // namespace certpolysolve::detail
