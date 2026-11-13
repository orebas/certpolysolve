// The sources that iterations and certificates evaluate: MonicSource (a monic polynomial of known degree, evaluated as
// enclosures of p and p' over a query ball), PolynomialSource on coefficients, SecularSource on native secular input.
//
// PolynomialSource evaluates its monic polynomial from its encloser's coefficient balls at the precision asked
// (evaluate_impl). A centred source evaluates the uncentred polynomial at z + c where that expansion's largest term is
// more than a bit below the centred one's (origin_point, from TermHull's upper hulls of (k, log2 |a_k|)), as the
// accuracy goals are relative to the original origin. With exact translations (exact_taylor), below the precision
// budget (exact_precision) the coefficient value is kept where it excludes zero; otherwise the exact translation to a
// coarse dyadic centre (cached) is evaluated at z - centre (evaluate_translated). The secular restarts take local
// models about an exact centre, by the same origin rule: taylor, the whole expansion; taylor_pair, coefficients k - 1
// and k from one joint evaluation of p^(k-1) / (k-1)!; taylor_head, the first coefficients by Horner's rule on exact
// accumulators with one disk radius of error each, or term by term where cheaper (cost: taylor_head_passes);
// taylor_tail_bound, a bound on what the head drops, from magnitudes alone. Other evaluators take its coefficients
// (coefficient_enclosure, origin_enclosure), dense or sparse (rectangular_available, sparse_available).
//
// SecularSource evaluates Q (1 - sum a_i / (x - b_i)), Q = prod (x - b_i), and its derivative by streaming Q and the
// numerator, never dividing by x - b_i, on poles and residues enclosed at the precision (enclose_terms). The seeds of
// secular_block_seeds are eigenvalue guesses for diag(b) + a 1^T, whose characteristic polynomial this is: b_i + a_i,
// or, for greedily matched pairs closer than a few times sqrt|a_i a_j| (a power-of-two score test), their 2-by-2
// principal block's eigenvalues.
//
// State: PolynomialSource owns its encloser and options and three caches of at most eight entries, filled on demand:
// prepared_ and origin_prepared_ by precision (lowest evicted first), translated_ by centre and precision (oldest
// first), whose Coefficients handles outlive eviction. SecularSource: poles_ and residues_ at precision_.
// Entry points: MonicSource, which secular_refinement.cpp, candidates.cpp, piecewise_roots.cpp, program_seeds.cpp and
// program_certificate.cpp take. PolynomialSource, built in solver.cpp and, for machine_newton_cover, in
// classical_isolation.cpp and prepared_piecewise.cpp; secular_refinement.cpp and candidates.cpp downcast to it.
// SecularSource and secular_block_seeds: solver.cpp. The research build's exact residues (exact_fits,
// evaluate_exact): research/solvers/exact_residues.cpp.
// Not here: the evaluators (CoefficientEvaluator: evaluator.cpp), the exact translation (ExactTaylorCache,
// centering.cpp), each factor's enclosers and centre (factor_state.cpp), ProgramSource (program_source.cpp).
// Design notes: docs/mathematics.md (Secular inputs and evaluation groundwork; C2 source and representation contracts).
#include "numeric/evaluate/structured_source.hpp"
#include "core/service_meter.hpp"
#include "core/heuristics.hpp"
#include "core/query.hpp"
#include <climits>
#include <limits>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>

namespace certpolysolve::detail {
ValueEvaluation MonicSource::value(const ComplexBall &z, slong precision) {
    if (precision < 2)
        throw std::invalid_argument("evaluation precision must be at least two bits");
    ValueEvaluation result;
    result.precision = precision;
    if (evaluate_value(result.value.value(), z.value(), precision))
        result.status = EvaluationStatus::Enclosed;
    return result;
}
void PolynomialSource::TermHull::build(const std::vector<double> &log2_coefficients) {
    vertices.clear();
    available = false;
    if (log2_coefficients.empty() || !std::isfinite(log2_coefficients.back()))
        return;
    for (size_t k = 0; k < log2_coefficients.size(); ++k) {
        if (!std::isfinite(log2_coefficients[k]))
            continue;
        const std::pair<double, double> q{static_cast<double>(k), log2_coefficients[k]};
        while (vertices.size() >= 2) {
            const auto &a = vertices[vertices.size() - 2], &b = vertices.back();
            if ((b.second - a.second) * (q.first - b.first) > (q.second - b.second) * (b.first - a.first))
                break;
            vertices.pop_back();
        }
        vertices.push_back(q);
    }
    leading = log2_coefficients.back();
    available = true;
}
double PolynomialSource::TermHull::largest_log2(double t) const {
    if (t == -HUGE_VAL)
        return vertices.front().first == 0 ? vertices.front().second - leading : -HUGE_VAL;
    // The hull's slopes fall along it: the largest term at 2^t is at the first vertex whose next edge
    // rises by less than -t a step.
    size_t low = 0, high = vertices.size() - 1;
    while (low < high) {
        const size_t mid = (low + high) / 2;
        const auto &a = vertices[mid], &b = vertices[mid + 1];
        if ((b.second - a.second) / (b.first - a.first) + t > 0)
            low = mid + 1;
        else
            high = mid;
    }
    return vertices[low].second + vertices[low].first * t - leading;
}
bool PolynomialSource::origin_point(acb_t w, const acb_t z) {
    if (!options_.origin || !centred_terms_.available || !origin_terms_.available || !acb_is_finite(z))
        return false;
    Magnitude size;
    acb_get_mag(size.value(), z);
    const double centred_t = mag_is_zero(size.value()) ? -HUGE_VAL : mag_get_d_log2_approx(size.value());
    // z + c exactly: the point's distance from p's origin is what the comparison and the evaluation need,
    // and a rounded sum would lose exactly that. Points far from both origins compare by magnitude alone.
    acb_add(w, z, options_.origin_offset.value(), ARF_PREC_EXACT);
    acb_get_mag(size.value(), w);
    const double origin_t = mag_is_zero(size.value()) ? -HUGE_VAL : mag_get_d_log2_approx(size.value());
    // A bit of hysteresis: where the two expansions are within a bit, the centred one stays.
    return origin_terms_.largest_log2(origin_t) + 1 < centred_terms_.largest_log2(centred_t);
}
PolynomialSource::Coefficients PolynomialSource::origin_enclosure(slong prec) {
    if (!origin_prepared_.contains(prec)) {
        if (origin_prepared_.size() >= 8)
            origin_prepared_.erase(origin_prepared_.begin());
        ComplexPolynomial polynomial;
        options_.origin(polynomial.value(), prec);
        auto entry = std::make_shared<Prepared>(std::move(polynomial), prec);
        entry->evaluator = std::make_unique<CoefficientEvaluator>(entry->polynomial);
        origin_prepared_[prec] = std::move(entry);
    }
    return origin_prepared_.at(prec);
}
PolynomialSource::Coefficients PolynomialSource::coefficient_enclosure(slong prec) {
    if (!prepared_.contains(prec)) {
        if (prepared_.size() >= 8)
            prepared_.erase(prepared_.begin());
        ComplexPolynomial polynomial;
        enclose_(polynomial.value(), prec);
        auto entry = std::make_shared<Prepared>(std::move(polynomial), prec);
#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
        if (local())
            entry->evaluator = std::make_unique<LocalPolynomialEvaluator>(entry->polynomial);
        else
#endif
            entry->evaluator = std::make_unique<CoefficientEvaluator>(entry->polynomial);
        prepared_[prec] = std::move(entry);
    }
    return prepared_.at(prec);
}
bool PolynomialSource::evaluate(acb_t value, acb_t derivative, const acb_t z, slong prec) {
    return evaluate_impl(value, derivative, z, prec, true);
}
bool PolynomialSource::evaluate_value(acb_t value, const acb_t z, slong prec) {
    ComplexBall derivative;
    return evaluate_impl(value, derivative.value(), z, prec, false);
}
PolynomialSource::Translated PolynomialSource::evaluate_translated(const ExactTaylorEncloser &encloser, acb_t value,
                                                                    acb_t derivative, const acb_t z, slong prec,
                                                                    bool joint) {
    if (!encloser || !acb_is_finite(z))
        return Translated::Unavailable;
    // Coarse dyadic centres coalesce nearby queries into a reusable exact
    // translation. All coefficients are retained, so the model is valid
    // on every query ball; this is a cache choice, not a proof threshold.
    constexpr slong exact_center_bits = heuristics::local_center_bits;
    constexpr size_t translated_cache_entries = 8;
    ComplexBall center, delta;
    const slong scale = std::max(arf_abs_bound_lt_2exp_si(arb_midref(acb_realref(z))),
                                 arf_abs_bound_lt_2exp_si(arb_midref(acb_imagref(z))));
    Float scaled;
    Integer grid;
    for (bool imaginary : {false, true})
        if (scale < WORD_MAX && scale > WORD_MIN + exact_center_bits) {
            const auto coordinate = imaginary ? acb_imagref(z) : acb_realref(z);
            const auto target = imaginary ? acb_imagref(center.value()) : acb_realref(center.value());
            arf_mul_2exp_si(scaled.value(), arb_midref(coordinate), exact_center_bits - scale);
            arf_get_fmpz(grid.value(), scaled.value(), ARF_RND_NEAR);
            arb_set_fmpz(target, grid.value());
            arb_mul_2exp_si(target, target, scale - exact_center_bits);
        }
    Coefficients translated;
    for (const auto &[cached_center, entry] : translated_)
        if (entry->precision == prec && acb_equal(cached_center.value(), center.value())) {
            translated = entry;
            break;
        }
    if (!translated && !acb_is_zero(center.value())) {
        ComplexPolynomial polynomial;
        if (encloser(polynomial.value(), center.value(), prec) &&
            acb_poly_degree(polynomial.value()) == degree_ &&
            !acb_contains_zero(polynomial.value()->coeffs + degree_)) {
            auto entry = std::make_shared<Prepared>(std::move(polynomial), prec);
            entry->evaluator = std::make_unique<CoefficientEvaluator>(entry->polynomial);
            translated = entry;
            if (translated_.size() >= translated_cache_entries)
                translated_.erase(translated_.begin());
            translated_.emplace_back(center, std::move(entry));
        }
    }
    if (!translated)
        return Translated::Unavailable;
    acb_sub(delta.value(), z, center.value(), prec);
    const bool evaluated = joint ? translated->evaluator->evaluate(value, derivative, delta.value(), prec)
                                 : translated->evaluator->evaluate_value(value, delta.value(), prec);
    if (!evaluated)
        return Translated::Unavailable;
    const acb_struct *leading = translated->polynomial.value()->coeffs + degree_;
    acb_div(value, value, leading, prec);
    if (joint)
        acb_div(derivative, derivative, leading, prec);
    return acb_is_finite(value) && (!joint || acb_is_finite(derivative)) ? Translated::Evaluated : Translated::Failed;
}
bool PolynomialSource::evaluate_impl(acb_t value, acb_t derivative, const acb_t z, slong prec, bool joint) {
    const auto evaluate = [&](const Evaluator &evaluator, const acb_t point) {
        return joint ? evaluator.evaluate(value, derivative, point, prec)
                     : evaluator.evaluate_value(value, point, prec);
    };
    const auto normalize = [&](const acb_t leading) {
        acb_div(value, value, leading, prec);
        if (joint)
            acb_div(derivative, derivative, leading, prec);
        return acb_is_finite(value) && (!joint || acb_is_finite(derivative));
    };
    // The value from the coefficients, as every source gives it.
    const auto plain = [&] {
        ComplexBall shifted;
        if (origin_point(shifted.value(), z)) {
            // The monic centred polynomial at z is the monic uncentred one at z + c, and so is its derivative.
            const auto origin = origin_enclosure(prec);
            const auto *polynomial = origin->polynomial.value();
            if (acb_poly_degree(polynomial) != degree_ || acb_contains_zero(polynomial->coeffs + degree_) ||
                !evaluate(*origin->evaluator, shifted.value()))
                return false;
            return normalize(polynomial->coeffs + degree_);
        }
        const auto handle = coefficient_enclosure(prec);
        const auto &entry = *handle;
        if (acb_poly_degree(entry.polynomial.value()) != degree_ ||
            acb_contains_zero(entry.polynomial.value()->coeffs + degree_) ||
            !evaluate(*entry.evaluator, z))
            return false;
        return normalize(entry.polynomial.value()->coeffs + degree_);
    };
    if (!options_.exact_taylor)
        return plain();
    // A source with exact translations (small exact factors). An exact translation gives the value rounded
    // once, whatever cancels, and costs an exact Taylor shift, n^2 exact operations, for every point not about
    // a centre in the cache; Horner's rule costs n operations at the precision asked, and a caller that needs
    // more of the value's bits asks again with more precision. So below the precision budget
    // (exact_precision) the value from the coefficients is taken wherever it has a correct bit, and the
    // translation where it has none (a point on a root is an exact zero only there) and at the budget, where
    // no more precision can be asked for. Taken for every point, the translations would be most of such
    // inputs' time. History: docs/code-history-20261021.md, structured_source.cpp 1.
    ComplexBall plain_value, plain_derivative;
    bool known = false;
    if (prec < options_.exact_precision) {
        known = plain();
        if (known && !acb_contains_zero(value))
            return true;
        if (known) {
            acb_set(plain_value.value(), value);
            if (joint)
                acb_set(plain_derivative.value(), derivative);
        }
    }
    switch (evaluate_translated(options_.exact_taylor, value, derivative, z, prec, joint)) {
    case Translated::Evaluated: return true;
    case Translated::Failed: return false;
    case Translated::Unavailable: break;
    }
    if (prec >= options_.exact_precision)
        return plain();
    if (known) {
        acb_set(value, plain_value.value());
        if (joint)
            acb_set(derivative, plain_derivative.value());
    }
    return known;
}
// Whether the source is evaluated from its coefficients (not an exact
// translation or the research build's local model), and then whether
// its evaluator treats them as sparse.
std::optional<bool> PolynomialSource::coefficient_sparsity(slong prec) {
    if (options_.exact_taylor || local())
        return std::nullopt;
    const auto handle = coefficient_enclosure(prec);
    const auto *evaluator = dynamic_cast<const CoefficientEvaluator *>(handle->evaluator.get());
    if (!evaluator)
        return std::nullopt;
    return evaluator->sparse();
}
bool PolynomialSource::rectangular_available(slong prec) {
    // Keep exact translations, sparse evaluation and experimental
    // local models intact. This probe changes only a dense coefficient schedule.
    return coefficient_sparsity(prec) == false;
}
bool PolynomialSource::sparse_available(slong prec) { return coefficient_sparsity(prec) == true; }
bool PolynomialSource::taylor_pair(acb_t lower, acb_t upper, const acb_t center, slong prec, slong k) {
    if (!enclose_ || !acb_is_exact(center) || !acb_is_finite(center) || k < 1 || k > degree_)
        return false;
    ComplexBall shifted;
    ComplexPolynomial original;
    const bool uncentred = origin_point(shifted.value(), center);
    if (uncentred)
        options_.origin(original.value(), prec);
    else if (options_.exact_taylor) {
        // An exact local model (small inputs): its own expansion, whole.
        ComplexPolynomial whole;
        if (!taylor(whole.value(), center, prec) || acb_poly_length(whole.value()) <= k)
            return false;
        acb_set(lower, whole.value()->coeffs + k - 1);
        acb_set(upper, whole.value()->coeffs + k);
        return true;
    } else
        enclose_(original.value(), prec);
    if (acb_poly_degree(original.value()) != degree_ || acb_contains_zero(original.value()->coeffs + degree_))
        return false;
    // q = p^(k-1) / (k-1)! = sum over i >= k - 1 of C(i, k - 1) a_i x^(i - k + 1): coefficient k - 1 of the
    // expansion is q(center), coefficient k is q'(center) / k.
    ComplexPolynomial derivative;
    ComplexBall term;
    Integer binomial;
    fmpz_one(binomial.value());
    for (slong i = k - 1; i <= degree_; ++i) {
        acb_mul_fmpz(term.value(), original.value()->coeffs + i, binomial.value(), prec);
        acb_poly_set_coeff_acb(derivative.value(), i - k + 1, term.value());
        // C(i + 1, k - 1) = C(i, k - 1) (i + 1) / (i - k + 2), exactly.
        fmpz_mul_ui(binomial.value(), binomial.value(), static_cast<ulong>(i + 1));
        fmpz_divexact_ui(binomial.value(), binomial.value(), static_cast<ulong>(i - k + 2));
    }
    acb_poly_evaluate2(lower, upper, derivative.value(), uncentred ? shifted.value() : center, prec);
    acb_div_ui(upper, upper, static_cast<ulong>(k), prec);
    acb_div(lower, lower, original.value()->coeffs + degree_, prec);
    acb_div(upper, upper, original.value()->coeffs + degree_, prec);
    ++taylor_pairs_;
    return acb_is_finite(lower) && acb_is_finite(upper);
}

namespace {
// The multiplications of an expansion's head summed term by term: for each of the `support` terms a power of
// the centre by repeated squaring, at most two multiplications a bit of the exponent, and one product for each
// of the head's coefficients.
double termwise_head_multiplications(slong support, slong degree, slong length) {
    return static_cast<double>(support) *
           static_cast<double>(2 * static_cast<slong>(FLINT_BIT_COUNT(static_cast<ulong>(degree))) + length);
}
slong support_of(const acb_poly_struct *polynomial) {
    slong support = 0;
    for (slong i = 0; i < acb_poly_length(polynomial); ++i)
        support += !acb_is_zero(polynomial->coeffs + i);
    return support;
}
} // namespace
std::optional<double> PolynomialSource::taylor_head_passes(slong prec, slong length) {
    if (!coefficient_sparsity(prec).has_value() || !enclose_ || length < 1 || length > degree_ + 1)
        return std::nullopt;
    const auto handle = coefficient_enclosure(prec);
    const double pass = static_cast<double>(degree_ + 1);
    return std::min(static_cast<double>(length),
                    termwise_head_multiplications(support_of(handle->polynomial.value()), degree_, length) / pass);
}

bool PolynomialSource::taylor_tail_bound(mag_t out, const acb_t center, const acb_t offset, slong prec,
                                         slong length) {
    if (options_.exact_taylor || !enclose_ || !acb_is_exact(center) || !acb_is_finite(center) ||
        !acb_is_finite(offset) || length < 1 || length > degree_ + 1)
        return false;
    ComplexBall shifted;
    const bool uncentred = origin_point(shifted.value(), center);
    const auto handle = uncentred ? origin_enclosure(prec) : coefficient_enclosure(prec);
    const auto original = handle->polynomial.value();
    const acb_struct *point = uncentred ? shifted.value() : center;
    if (acb_poly_degree(original) != degree_ || acb_contains_zero(original->coeffs + degree_))
        return false;
    // A term a_i x^i contributes sum over j from L = length to i of C(i, j) a_i p^(i - j) t^j at the point p,
    // with t the offset. Consecutive terms of that sum differ by the factor (i - j) / (j + 1) |t| / |p|, at
    // most r = (i - L) / (L + 1) |t| / |p| from j = L on; where r <= 1/2 the sum is at most twice its first
    // term (a geometric series), and otherwise at most |a_i| (|p| + |t|)^i, the whole binomial sum.
    Magnitude above, below, step, coefficient, first, part, ratio, whole;
    acb_get_mag(above.value(), point);
    acb_get_mag_lower(below.value(), point);
    acb_get_mag(step.value(), offset);
    mag_zero(out);
    for (slong i = length; i <= degree_; ++i) {
        if (acb_is_zero(original->coeffs + i))
            continue;
        acb_get_mag(coefficient.value(), original->coeffs + i);
        const auto lowest = static_cast<ulong>(length);
        const auto top = static_cast<ulong>(i);
        bool geometric = !mag_is_zero(below.value());
        if (geometric) {
            mag_mul_ui(ratio.value(), step.value(), top - lowest);
            mag_div(ratio.value(), ratio.value(), below.value());
            mag_div_ui(ratio.value(), ratio.value(), lowest + 1);
            geometric = mag_cmp_2exp_si(ratio.value(), -1) <= 0;
        }
        if (geometric) {
            // C(i, L) as the product of (i - m) / (L - m) over m < L: FLINT's mag_bin_uiui bounds it through
            // Stirling's formula, (e i / L)^L, about 4.5 times C(257, 3).
            mag_one(first.value());
            for (ulong m = 0; m < lowest; ++m) {
                mag_mul_ui(first.value(), first.value(), top - m);
                mag_div_ui(first.value(), first.value(), lowest - m);
            }
            mag_pow_ui(part.value(), above.value(), top - lowest);
            mag_mul(first.value(), first.value(), part.value());
            mag_pow_ui(part.value(), step.value(), lowest);
            mag_mul(first.value(), first.value(), part.value());
            mag_mul_2exp_si(part.value(), first.value(), 1);
        } else {
            mag_add(whole.value(), above.value(), step.value());
            mag_pow_ui(part.value(), whole.value(), top);
        }
        mag_mul(part.value(), part.value(), coefficient.value());
        mag_add(out, out, part.value());
    }
    // taylor()'s normalisation: the monic polynomial.
    acb_get_mag_lower(part.value(), original->coeffs + degree_);
    mag_div(out, out, part.value());
    return mag_is_finite(out);
}

bool PolynomialSource::taylor_head(acb_poly_t out, const acb_t center, slong prec, slong length) {
    if (options_.exact_taylor || !enclose_ || !acb_is_exact(center) || !acb_is_finite(center) ||
        length < 1 || length > degree_ + 1)
        return false;
    ComplexBall shifted;
    const bool uncentred = origin_point(shifted.value(), center);
    const auto handle = uncentred ? origin_enclosure(prec) : coefficient_enclosure(prec);
    const auto original = handle->polynomial.value();
    const acb_struct *point = uncentred ? shifted.value() : center;
    if (acb_poly_degree(original) != degree_ || acb_contains_zero(original->coeffs + degree_) || !acb_is_exact(point))
        return false;
    if (termwise_head_multiplications(support_of(original), degree_, length) <
        static_cast<double>(degree_ + 1) * static_cast<double>(length)) {
        // Term by term, in ball arithmetic. A power's chain of at most 2 bits(degree) multiplications and the
        // products down the head widen a rectangle by at most half a bit each at a complex point; that many
        // bits beyond the precision keep the sum's terms at the precision asked for.
        const slong chain = 2 * static_cast<slong>(FLINT_BIT_COUNT(static_cast<ulong>(degree_))) + length;
        const slong inner = resources::precision_sum(prec, chain);
        std::vector<ComplexBall> sums(static_cast<size_t>(length));
        ComplexBall power, term;
        Integer binomial;
        for (slong i = 0; i <= degree_; ++i) {
            if (acb_is_zero(original->coeffs + i))
                continue;
            const slong top = std::min(length - 1, i);
            acb_pow_ui(power.value(), point, static_cast<ulong>(i - top), inner);
            for (slong j = top; j >= 0; --j) {
                fmpz_bin_uiui(binomial.value(), static_cast<ulong>(i), static_cast<ulong>(j));
                acb_mul(term.value(), original->coeffs + i, power.value(), inner);
                acb_mul_fmpz(term.value(), term.value(), binomial.value(), inner);
                acb_add(sums[static_cast<size_t>(j)].value(), sums[static_cast<size_t>(j)].value(), term.value(),
                        inner);
                if (j)
                    acb_mul(power.value(), power.value(), point, inner);
            }
        }
        ++taylor_heads_;
        acb_poly_zero(out);
        for (slong j = 0; j < length; ++j) {
            const auto sum = sums[static_cast<size_t>(j)].value();
            acb_div(sum, sum, original->coeffs + degree_, prec);
            if (!acb_is_finite(sum))
                return false;
            acb_poly_set_coeff_acb(out, j, sum);
        }
        return true;
    }
    // errors[j] bounds the modulus of accumulator j's error: it is multiplied by |x| with the accumulator, takes
    // the error of the accumulator below (or the coefficient's radii), and the rounding of the step's two
    // operations, whose radii the accumulator then gives up.
    std::vector<ComplexBall> rows(static_cast<size_t>(length));
    std::vector<Magnitude> errors(static_cast<size_t>(length));
    Magnitude norm, rounding;
    acb_get_mag(norm.value(), point);
    for (slong i = degree_; i >= 0; --i)
        // The accumulators above the terms taken so far are still zero.
        for (slong j = std::min(length - 1, degree_ - i); j >= 0; --j) {
            const auto row = rows[static_cast<size_t>(j)].value();
            const auto error = errors[static_cast<size_t>(j)].value();
            acb_mul(row, row, point, prec);
            acb_add(row, row, j ? rows[static_cast<size_t>(j - 1)].value() : original->coeffs + i, prec);
            mag_mul(error, error, norm.value());
            if (j)
                mag_add(error, error, errors[static_cast<size_t>(j - 1)].value());
            mag_add(rounding.value(), arb_radref(acb_realref(row)), arb_radref(acb_imagref(row)));
            mag_add(error, error, rounding.value());
            acb_get_mid(row, row);
        }
    ++taylor_heads_;
    acb_poly_zero(out);
    for (slong j = 0; j < length; ++j) {
        const auto row = rows[static_cast<size_t>(j)].value();
        acb_add_error_mag(row, errors[static_cast<size_t>(j)].value());
        acb_div(row, row, original->coeffs + degree_, prec);
        if (!acb_is_finite(row))
            return false;
        acb_poly_set_coeff_acb(out, j, row);
    }
    return true;
}

bool PolynomialSource::taylor(acb_poly_t out, const acb_t center, slong prec) {
    last_taylor_quote_ = {};
    if (!enclose_ || !acb_is_exact(center) || !acb_is_finite(center))
        return false;
    ComplexBall shifted;
    if (origin_point(shifted.value(), center)) {
        // The local model about a centre nearer p's origin: p shifted by center + c (origin_point).
        ComplexPolynomial uncentred;
        options_.origin(uncentred.value(), prec);
        if (acb_poly_degree(uncentred.value()) != degree_ || acb_contains_zero(uncentred.value()->coeffs + degree_))
            return false;
        if (service_meter_) last_taylor_quote_ = service_taylor_quote(uncentred.value(), shifted.value(), prec);
        acb_poly_taylor_shift(out, uncentred.value(), shifted.value(), prec);
        ++taylor_shifts_;
        for (slong j = 0; j < acb_poly_length(out); ++j)
            acb_div(out->coeffs + j, out->coeffs + j, uncentred.value()->coeffs + degree_, prec);
        for (slong j = 0; j < acb_poly_length(out); ++j)
            if (!acb_is_finite(out->coeffs + j))
                return false;
        return true;
    }
    ComplexPolynomial original;
    if (options_.exact_taylor && options_.exact_taylor(original.value(), center, prec)) {
        // The owned cache pays its exact construction and this enclosure.
        // It is not an Arb shift of the already translated coefficients.
        if (service_meter_ && !options_.exact_taylor_work)
            last_taylor_quote_ = {ServiceOperation::Taylor, service_taylor_base(original.value(), center, prec)};
        ++taylor_shifts_;
        if (acb_poly_degree(original.value()) != degree_ || acb_contains_zero(original.value()->coeffs + degree_))
            return false;
        ComplexBall leading;
        acb_set(leading.value(), original.value()->coeffs + degree_);
        for (slong j = 0; j <= degree_; ++j)
            acb_div(original.value()->coeffs + j, original.value()->coeffs + j, leading.value(), prec);
        acb_poly_set(out, original.value());
        return true;
    }
    enclose_(original.value(), prec);
    if (acb_poly_degree(original.value()) != degree_ || acb_contains_zero(original.value()->coeffs + degree_))
        return false;
    if (service_meter_) last_taylor_quote_ = service_taylor_quote(original.value(), center, prec);
    acb_poly_taylor_shift(out, original.value(), center, prec);
    ++taylor_shifts_;
    for (slong j = 0; j < acb_poly_length(out); ++j)
        acb_div(out->coeffs + j, out->coeffs + j, original.value()->coeffs + degree_, prec);
    for (slong j = 0; j < acb_poly_length(out); ++j)
        if (!acb_is_finite(out->coeffs + j))
            return false;
    return true;
}

SecularSource::SecularSource(const cps_secular &input) : input_(input) {
    for (const auto &term : input_.terms)
        real_ = real_ && fmpq_is_zero(term.residue.imag.value()) && fmpq_is_zero(term.pole.imag.value());
}
void SecularSource::enclose_terms(slong prec) {
    if (precision_ == prec)
        return;
    poles_.resize(input_.terms.size());
    residues_.resize(input_.terms.size());
    for (size_t i = 0; i < input_.terms.size(); ++i) {
        const auto &term = input_.terms[i];
        arb_set_fmpq(acb_realref(poles_[i].value()), term.pole.real.value(), prec);
        arb_set_fmpq(acb_imagref(poles_[i].value()), term.pole.imag.value(), prec);
        arb_set_fmpq(acb_realref(residues_[i].value()), term.residue.real.value(), prec);
        arb_set_fmpq(acb_imagref(residues_[i].value()), term.residue.imag.value(), prec);
    }
    precision_ = prec;
}
bool SecularSource::evaluate(acb_t value, acb_t derivative, const acb_t z, slong prec) {
    enclose_terms(prec);
    // Streaming Q and the numerator sum avoids every division by x-b_i,
    // including on a disk meeting several poles. Its arithmetic error is
    // enclosed uniformly; callers may increase precision when cancellation wins.
    ComplexBall product, dproduct, numerator, dnumerator, difference, next, dnext, term;
    acb_one(product.value());
    for (size_t i = 0; i < poles_.size(); ++i) {
        acb_sub(difference.value(), z, poles_[i].value(), prec);
        acb_mul(next.value(), numerator.value(), difference.value(), prec);
        acb_mul(term.value(), residues_[i].value(), product.value(), prec);
        acb_add(next.value(), next.value(), term.value(), prec);
        acb_mul(dnext.value(), dnumerator.value(), difference.value(), prec);
        acb_add(dnext.value(), dnext.value(), numerator.value(), prec);
        acb_mul(term.value(), residues_[i].value(), dproduct.value(), prec);
        acb_add(dnext.value(), dnext.value(), term.value(), prec);
        acb_set(numerator.value(), next.value());
        acb_set(dnumerator.value(), dnext.value());
        acb_mul(dproduct.value(), dproduct.value(), difference.value(), prec);
        acb_add(dproduct.value(), dproduct.value(), product.value(), prec);
        acb_mul(product.value(), product.value(), difference.value(), prec);
    }
    acb_sub(value, product.value(), numerator.value(), prec);
    acb_sub(derivative, dproduct.value(), dnumerator.value(), prec);
    return acb_is_finite(value) && acb_is_finite(derivative);
}
bool secular_block_seeds(acb_ptr out, SecularSource &source, slong prec, const cps_options &options) {
    const size_t n = static_cast<size_t>(source.degree());
    if (n < 16 || n > heuristics::secular_block_seed_degree)
        return false;
    source.enclose_terms(prec);
    const auto &a = source.residues();
    const auto &b = source.poles();
    std::vector<ComplexBall> diagonal(n);
    std::vector<double> heights(n);
    Float bound;
    for (size_t i = 0; i < n; ++i) {
        acb_add(diagonal[i].value(), b[i].value(), a[i].value(), prec);
        acb_get_abs_ubound_arf(bound.value(), a[i].value(), 64);
        if (!arf_is_finite(bound.value()) || arf_is_zero(bound.value()))
            return false;
        heights[i] = static_cast<double>(arf_abs_bound_lt_2exp_si(bound.value()));
        acb_get_mid(out + static_cast<slong>(i), diagonal[i].value());
    }
    struct Edge {
        double score;
        size_t i, j;
    };
    std::vector<Edge> edges;
    ComplexBall difference;
    for (size_t i = 0; i < n; ++i) {
        if (options.cancelled())
            return false;
        for (size_t j = i + 1; j < n; ++j) {
            acb_sub(difference.value(), diagonal[i].value(), diagonal[j].value(), prec);
            acb_get_abs_ubound_arf(bound.value(), difference.value(), 64);
            if (!arf_is_finite(bound.value()))
                continue;
            // Coincident entries have no distance; -1e12 makes their score
            // smaller than any real pair's, so they are always linked.
            const double distance = arf_is_zero(bound.value())
                                        ? -1.e12
                                        : static_cast<double>(arf_abs_bound_lt_2exp_si(bound.value()));
            const double score = distance - (heights[i] + heights[j]) / 2;
            if (score <= 2)
                edges.push_back({score, i, j});
        }
    }
    std::sort(edges.begin(), edges.end(), [](const Edge &x, const Edge &y) {
        if (x.score != y.score)
            return x.score < y.score;
        if (x.i != y.i)
            return x.i < y.i;
        return x.j < y.j;
    });
    std::vector<uint8_t> used(n, 0);
    ComplexBall center, discriminant, product, left, right;
    bool paired = false;
    for (const auto &edge : edges)
        if (!used[edge.i] && !used[edge.j]) {
            if (options.cancelled())
                return false;
            acb_add(center.value(), diagonal[edge.i].value(), diagonal[edge.j].value(), prec);
            acb_mul_2exp_si(center.value(), center.value(), -1);
            acb_sub(discriminant.value(), diagonal[edge.i].value(), diagonal[edge.j].value(), prec);
            acb_mul_2exp_si(discriminant.value(), discriminant.value(), -1);
            acb_sqr(discriminant.value(), discriminant.value(), prec);
            acb_mul(product.value(), a[edge.i].value(), a[edge.j].value(), prec);
            acb_add(discriminant.value(), discriminant.value(), product.value(), prec);
            acb_sqrt(discriminant.value(), discriminant.value(), prec);
            acb_add(left.value(), center.value(), discriminant.value(), prec);
            acb_sub(right.value(), center.value(), discriminant.value(), prec);
            if (!acb_is_finite(left.value()) || !acb_is_finite(right.value()))
                continue;
            acb_get_mid(out + static_cast<slong>(edge.i), left.value());
            acb_get_mid(out + static_cast<slong>(edge.j), right.value());
            used[edge.i] = used[edge.j] = 1;
            paired = true;
        }
    return paired;
}
} // namespace certpolysolve::detail
