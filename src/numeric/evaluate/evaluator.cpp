// An implementation of the Evaluator contract (evaluator.hpp: p(z) and p'(z) enclosed over the ball z and the
// coefficient balls; false is unresolved, never a statement about roots): CoefficientEvaluator, from the coefficients.
//
// CoefficientEvaluator applies Horner's rule for p and p' across the gaps between the exponents it lists (those of the
// nonzero coefficients, and 0) when they are fewer than a quarter of degree + 1, and acb_poly_evaluate2 otherwise.
//
// State: none at file scope; an evaluator keeps its coefficient snapshot and the exponents it lists.
// Entry points: CoefficientEvaluator in structured_source.cpp (PolynomialSource), candidates.cpp (regenerate_secular).
// The local Taylor cache (LocalPolynomialEvaluator, the evaluator of CandidateEngine::Series, an engine nothing in src/
// selects) is the research build's since W8's S2 (research/solvers/local_polynomial_evaluator.cpp); the secular
// function's evaluator, which no library code called, is test support (tests/secular_evaluator.hpp).
// Not here: the native secular input's evaluation (SecularSource, structured_source.cpp).
// Design notes: docs/mathematics.md (Secular inputs and evaluation groundwork; Separate C1 experiments).
#include "numeric/evaluate/evaluator.hpp"
#include <stdexcept>

namespace certpolysolve::detail {
CoefficientEvaluator::CoefficientEvaluator(PolynomialSnapshot polynomial, bool sparse)
    : polynomial_(std::move(polynomial)), sparse_(sparse) {
    const auto p = polynomial_.value();
    for (slong k = acb_poly_degree(p); k >= 0; --k)
        if (!acb_is_zero(p->coeffs + k) || k == 0)
            indices_.push_back(k);
    sparse_ = sparse_ && indices_.size() * 4 < static_cast<size_t>(acb_poly_length(p));
}
bool CoefficientEvaluator::evaluate(acb_t value, acb_t derivative, const acb_t z, slong prec) const {
    if (!sparse_) {
        acb_poly_evaluate2(value, derivative, polynomial_.value(), z, prec);
    } else {
        ComplexBall power, dpower, temporary;
        acb_set(value, polynomial_.value()->coeffs + indices_[0]);
        acb_zero(derivative);
        for (size_t j = 1; j < indices_.size(); ++j) {
            const auto gap = static_cast<ulong>(indices_[j - 1] - indices_[j]);
            acb_pow_ui(dpower.value(), z, gap - 1, prec);
            acb_mul(power.value(), dpower.value(), z, prec);
            acb_mul_ui(dpower.value(), dpower.value(), gap, prec);
            acb_mul(derivative, derivative, power.value(), prec);
            acb_mul(temporary.value(), value, dpower.value(), prec);
            acb_add(derivative, derivative, temporary.value(), prec);
            acb_mul(value, value, power.value(), prec);
            acb_add(value, value, polynomial_.value()->coeffs + indices_[j], prec);
        }
    }
    return acb_is_finite(value) && acb_is_finite(derivative);
}

bool CoefficientEvaluator::evaluate_value(acb_t value, const acb_t z, slong prec) const {
    if (!sparse_) {
        acb_poly_evaluate(value, polynomial_.value(), z, prec);
    } else {
        ComplexBall power;
        acb_set(value, polynomial_.value()->coeffs + indices_[0]);
        for (size_t j = 1; j < indices_.size(); ++j) {
            const auto gap = static_cast<ulong>(indices_[j - 1] - indices_[j]);
            acb_pow_ui(power.value(), z, gap, prec);
            acb_mul(value, value, power.value(), prec);
            acb_add(value, value, polynomial_.value()->coeffs + indices_[j], prec);
        }
    }
    return acb_is_finite(value);
}
} // namespace certpolysolve::detail
