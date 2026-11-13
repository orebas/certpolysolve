// Exact preparation of coefficient input: exact real coefficients as one integer polynomial
// (streamed_real_polynomial), its square-free factors over Z (squarefree_factors, costed by squarefree_effort), and
// the piecewise route's whole preparation (prepare_squarefree_input).
//
// The conversion clears the denominators' least common multiple, charging each expansion of a compressed coefficient to
// the exact storage limit before it is made: one coefficient at a time, in two passes, after the effort,
// temporary-storage and exponent-expansion limits (the compact-exponent rule: compact_expansion_fits with one
// precision cap a number, cps_options::exact_expansion_allowance), and it can report why it refused
// (PreparationDiagnostic).
// squarefree_factors is FLINT's fmpz_poly_factor_squarefree, which cannot be cancelled, so callers first fit
// squarefree_effort to the effort allowance. prepare_squarefree_input sets aside x^zero_order; with power_composition,
// writes the rest as q(x^k), q(0) != 0 (prepare_power_composition), and if q's preparation is available multiplies its
// factors' powers by k (x -> x^k is unramified over q's roots); else proves the rest square-free by degree <= 1 or by
// nonzero_part_squarefree_modular (modular_squarefree.cpp), returning no factors; else factors it over Z when it is real
// and within budget, or over Q(i) (exact_complex_factors, exact_complex.cpp); and with power_composition it then
// reduces each factor that is itself a power composition.
//
// State: none.
// Entry points: streamed_real_polynomial (SolveState::Impl::preprocess in solver.cpp, and every_root_seeded in
// real_line_routing.cpp); squarefree_effort and squarefree_factors (preprocess); prepare_squarefree_input
// (PreparedPiecewiseState::Impl::prepare in prepared_piecewise.cpp).
// Design notes: docs/mathematics.md ("Individual root certificates"), docs/resource-policy.md (the allowances).
#include "numeric/exact/exact_preprocessing.hpp"
#include "numeric/exact/modular_squarefree.hpp"
#include "numeric/exact/power_composition.hpp"
#include <algorithm>

namespace certpolysolve::detail {
bool streamed_real_polynomial(IntegerPolynomial &out, std::span<const Coefficient> coefficients,
                             const cps_options &options, PreparationDiagnostic *diagnostic) {
    if (diagnostic) *diagnostic = {};
    const auto refusal = [&](const char *reason) {
        if (diagnostic) diagnostic->reject(reason, options.exact_storage_limit());
        return false;
    };
    if (!options.exact_storage_limit() ||
        !options.exact_effort_fits(resources::multiply(coefficients.size(), coefficients.size())))
        return refusal("exact-preparation-effort-or-storage");
    const bool expansion =
        options.allow_exponent_expansion || compact_expansion_fits(coefficients, options.exact_expansion_allowance());
    if (!expansion) return refusal("compact-exponent-expansion-needs-explicit-allowance");
    if (!resources::storage_fits(resources::multiply(coefficients.size(), sizeof(fmpz)), options.temporary_bytes))
        return refusal("exact-preparation-temporary-storage");
    const auto charge = [](uint64_t &available, uint64_t bits) {
        if (bits > available) return false;
        available -= bits;
        return true;
    };
    Integer denominator;
    fmpz_one(denominator.value());
    uint64_t expanded_budget = options.exact_storage_limit();
    for (const auto &coefficient : coefficients) {
        Rational value;
        if (options.cancelled() || !coefficient.exact() || !coefficient.is_real()) {
            if (diagnostic) diagnostic->reason = options.cancelled() ? "cancelled" : "not-exact-real";
            return false;
        }
        if (!coefficient.real.rational_with_budget(value.value(), expanded_budget, expansion))
            return refusal("coefficient-expansion-budget");
        if (diagnostic) diagnostic->input_bits = options.exact_storage_limit() - expanded_budget;
        // The common denominator is the only object retained between values.
        // Bound both live operands and the uncancelled LCM product before GMP
        // materializes it. These are conservative bit guards, not RAM quotas.
        uint64_t live = options.exact_storage_limit();
        const auto previous = static_cast<uint64_t>(fmpz_bits(denominator.value()));
        const auto divisor = static_cast<uint64_t>(fmpz_bits(fmpq_denref(value.value())));
        if (!charge(live, previous) || !charge(live, previous) ||
            !charge(live, divisor) || !charge(live, divisor) ||
            !charge(live, static_cast<uint64_t>(fmpz_bits(fmpq_numref(value.value())))))
            return refusal("denominator-lcm-budget");
        fmpz_lcm(denominator.value(), denominator.value(), fmpq_denref(value.value()));
    }
    const auto denominator_bits = static_cast<uint64_t>(fmpz_bits(denominator.value()));
    uint64_t retained = denominator_bits;
    IntegerPolynomial proposed;
    for (size_t i = 0; i < coefficients.size(); ++i) {
        if (options.cancelled()) {
            if (diagnostic) diagnostic->reason = "cancelled";
            return false;
        }
        if (retained > options.exact_storage_limit()) return refusal("retained-coefficients-budget");
        uint64_t available = options.exact_storage_limit() - retained;
        Rational value;
        if (!coefficients[i].real.rational_with_budget(value.value(), available))
            return refusal("denominator-clearing-input-budget");
        if (fmpq_is_zero(value.value())) continue;
        const auto numerator_bits = static_cast<uint64_t>(fmpz_bits(fmpq_numref(value.value())));
        const auto divisor_bits = static_cast<uint64_t>(fmpz_bits(fmpq_denref(value.value())));
        if (divisor_bits > denominator_bits) {
            if (diagnostic) diagnostic->reason = "invalid-common-denominator";
            return false;
        }
        const auto multiplier_bits = denominator_bits - divisor_bits + 1;
        // Charge the multiplier and two product copies separately, avoiding
        // overflowing an aggregate estimate. The second copy is the stored
        // output coefficient while its local product still exists.
        if (!charge(available, multiplier_bits) ||
            !charge(available, numerator_bits) || !charge(available, multiplier_bits) ||
            !charge(available, numerator_bits) || !charge(available, multiplier_bits))
            return refusal("denominator-clearing-scratch-budget");
        Integer multiplier, integer;
        fmpz_divexact(multiplier.value(), denominator.value(), fmpq_denref(value.value()));
        fmpz_mul(integer.value(), multiplier.value(), fmpq_numref(value.value()));
        fmpz_poly_set_coeff_fmpz(proposed.value(), static_cast<slong>(i), integer.value());
        retained += static_cast<uint64_t>(fmpz_bits(integer.value()));
    }
    out = std::move(proposed);
    if (diagnostic) {
        diagnostic->reason = "accepted";
        diagnostic->accepted = true;
        diagnostic->output_bits = retained;
        diagnostic->required_bits_lower_bound = std::max(diagnostic->input_bits, retained);
    }
    return true;
}

uint64_t squarefree_effort(const IntegerPolynomial &polynomial) {
    const auto length = static_cast<uint64_t>(fmpz_poly_length(polynomial.value()));
    const auto words = 1 + static_cast<uint64_t>(FLINT_ABS(fmpz_poly_max_bits(polynomial.value()))) / FLINT_BITS;
    return resources::multiply(resources::multiply(length, length), words);
}

std::vector<IntegerFactor> squarefree_factors(const IntegerPolynomial &polynomial) {
    Factorization decomposition;
    fmpz_poly_factor_squarefree(decomposition.value(), polynomial.value());
    std::vector<IntegerFactor> result;
    for (slong i = 0; i < decomposition.value()->num; ++i) {
        IntegerFactor factor;
        fmpz_poly_set(factor.polynomial.value(), decomposition.value()->p + i);
        factor.multiplicity = static_cast<uint64_t>(decomposition.value()->exp[i]);
        result.push_back(std::move(factor));
    }
    return result;
}

SquarefreePreparation prepare_squarefree_input(const DensePolynomial &input, const cps_options &options,
                                               bool power_composition) {
    SquarefreePreparation result;
    size_t length = input.coefficients.size();
    while (length && input.coefficients[length - 1].zero()) --length;
    while (result.zero_order + 1 < length && input.coefficients[result.zero_order].zero())
        ++result.zero_order;
    if (options.cancelled()) {
        result.cancelled = true;
        result.status = "square-free preparation cancelled";
        return result;
    }
    if (!length || std::any_of(input.coefficients.begin(), input.coefficients.begin() +
            static_cast<ptrdiff_t>(length), [](const Coefficient &c) { return !c.exact(); })) {
        result.status = "square-free preparation requires a nonzero exact polynomial";
        return result;
    }
    const auto coefficients = std::span(input.coefficients).subspan(
        static_cast<size_t>(result.zero_order), length - static_cast<size_t>(result.zero_order));
    // Already elementary inputs have outward closed forms. Preserve that
    // cheaper path (including exact realness) without an unnecessary root map.
    if (power_composition && coefficients.size() > 3) {
        if (auto composition = prepare_power_composition(coefficients)) {
            auto reduced = prepare_squarefree_input(composition->polynomial, options, true);
            if (reduced.cancelled) return reduced;
            if (reduced.available()) {
                // q(0) != 0, so x -> x^k is unramified at every root. Square-
                // freeness, factor multiplicities and coprimality are preserved.
                if (reduced.nonzero_squarefree) {
                    reduced.factors.push_back({std::move(composition->polynomial), 1, 1});
                    reduced.nonzero_squarefree = false;
                }
                for (auto &factor : reduced.factors) factor.power *= composition->power;
                reduced.zero_order = result.zero_order;
                reduced.status += "; exact coefficient-index power reduction";
                return reduced;
            }
            // A power identity alone is not a square-freeness proof. Preserve
            // the established unavailable/raw path when preparation is refused.
        }
    }
    if (length <= result.zero_order + 2 || nonzero_part_squarefree_modular(input, options)) {
        result.nonzero_squarefree = true;
        result.status = "nonzero part square-free by degree or modular proof";
        return result;
    }
    IntegerPolynomial real;
    if (streamed_real_polynomial(real, coefficients, options) &&
        options.exact_effort_fits(squarefree_effort(real))) {
        for (auto &factor : squarefree_factors(real)) {
            ExactInputFactor owned;
            owned.multiplicity = factor.multiplicity;
            const auto size = fmpz_poly_length(factor.polynomial.value());
            owned.polynomial.coefficients.resize(static_cast<size_t>(size));
            for (slong j = 0; j < size; ++j)
                fmpq_set_fmpz(owned.polynomial.coefficients[static_cast<size_t>(j)].real.significand.value(),
                              factor.polynomial.value()->coeffs + j);
            result.factors.push_back(std::move(owned));
        }
        result.status = "budgeted exact real square-free preparation";
    } else {
        auto decomposition = exact_complex_factors(coefficients, options);
        result.cancelled = decomposition.cancelled;
        result.status = decomposition.status;
        for (auto &factor : decomposition.factors) {
            ExactInputFactor owned;
            owned.multiplicity = factor.multiplicity;
            owned.polynomial.coefficients.resize(factor.polynomial.size());
            for (size_t j = 0; j < factor.polynomial.size(); ++j) {
                owned.polynomial.coefficients[j].real.significand = std::move(factor.polynomial[j].real);
                owned.polynomial.coefficients[j].imag.significand = std::move(factor.polynomial[j].imag);
            }
            result.factors.push_back(std::move(owned));
        }
        if (result.factors.empty() && !result.cancelled)
            result.status = "square-free preparation unavailable within exact algebra budget";
    }
    if (options.cancelled() || result.cancelled) {
        result.factors.clear();
        result.cancelled = true;
        result.status = "square-free preparation cancelled";
    }
    if (power_composition && !result.cancelled) {
        for (auto &factor : result.factors) {
            if (factor.polynomial.coefficients.size() <= 3) continue;
            if (auto composition = prepare_power_composition(factor.polynomial.coefficients)) {
                factor.polynomial = std::move(composition->polynomial);
                factor.power *= composition->power;
            }
        }
    }
    return result;
}
} // namespace certpolysolve::detail
