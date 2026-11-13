#pragma once
// The secular function S(z) = -1 + sum_j a_j/(z - b_j) and its derivative, enclosed: the oracle test_representations'
// secular_identity checks regenerate_secular's residues against. Test support: until W8's S2 (23 October 2026) it was
// src/evaluator.cpp's SecularEvaluator, which no library code called; moved here unchanged but for `inline`.
#include "numeric/evaluate/evaluator.hpp"
#include <span>
#include <stdexcept>
#include <vector>

namespace certpolysolve::detail {
class SecularEvaluator final : public Evaluator {
    std::vector<ComplexBall> residues_, poles_;

  public:
    SecularEvaluator(std::span<const ComplexBall> residues, std::span<const ComplexBall> poles);
    slong degree() const override { return static_cast<slong>(poles_.size()); }
    bool evaluate(acb_t value, acb_t derivative, const acb_t z, slong precision) const override;
};
inline SecularEvaluator::SecularEvaluator(std::span<const ComplexBall> a, std::span<const ComplexBall> b)
    : residues_(a.begin(), a.end()), poles_(b.begin(), b.end()) {
    if (a.size() != b.size())
        throw std::invalid_argument("secular vector size mismatch");
}
inline bool SecularEvaluator::evaluate(acb_t value, acb_t derivative, const acb_t z, slong prec) const {
    acb_set_si(value, -1);
    acb_zero(derivative);
    ComplexBall difference, term;
    for (size_t j = 0; j < poles_.size(); ++j) {
        acb_sub(difference.value(), z, poles_[j].value(), prec);
        if (acb_contains_zero(difference.value())) {
            acb_indeterminate(value);
            acb_indeterminate(derivative);
            return false;
        }
        acb_div(term.value(), residues_[j].value(), difference.value(), prec);
        acb_add(value, value, term.value(), prec);
        acb_div(term.value(), term.value(), difference.value(), prec);
        acb_sub(derivative, derivative, term.value(), prec);
    }
    return acb_is_finite(value) && acb_is_finite(derivative);
}
} // namespace certpolysolve::detail
