#pragma once
#include "input/coefficients/polynomial_input.hpp"
#include "numeric/exact/exact_complex.hpp"
#include "numeric/exact/preparation_diagnostics.hpp"

namespace certpolysolve::detail {
struct IntegerFactor {
    IntegerPolynomial polynomial;
    uint64_t multiplicity = 1;
};
// Shared with the regular driver; the conversion preserves exact input and
// charges denominator expansion against the caller's existing algebra budget.
bool streamed_real_polynomial(IntegerPolynomial &, std::span<const Coefficient>, const cps_options &,
                              PreparationDiagnostic * = nullptr);
std::vector<IntegerFactor> squarefree_factors(const IntegerPolynomial &);
// Word operations charged to exact_operations for squarefree_factors: a
// gcd of p and p' costs about (n + 1)^2 products of coefficient-sized words.
// FLINT's factorisation cannot be cancelled, so this charge is what bounds it
// (principled audit, fix 14).
uint64_t squarefree_effort(const IntegerPolynomial &);

struct ExactInputFactor {
    DensePolynomial polynomial;
    uint64_t multiplicity = 1;
    // The original nonzero factor is polynomial(x^power)^multiplicity.
    // A nontrivial power always has a proved nonzero constant coefficient.
    ulong power = 1;
};
struct SquarefreePreparation {
    std::vector<ExactInputFactor> factors;
    uint64_t zero_order = 0;
    bool nonzero_squarefree = false;
    bool cancelled = false;
    std::string status;
    bool available() const { return nonzero_squarefree || !factors.empty(); }
};
// A successful modular test returns the original nonzero polynomial implicitly,
// without expanding compressed coefficients or allocating an owned copy.
// Otherwise return exact square-free factors of the nonzero part, or an honest
// unavailable result. A failed modular test does not establish repeated roots.
SquarefreePreparation prepare_squarefree_input(const DensePolynomial &, const cps_options &,
                                               bool power_composition = false);
} // namespace certpolysolve::detail
