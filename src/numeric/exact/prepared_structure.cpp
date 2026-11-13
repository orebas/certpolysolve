// Two overloads of input_structure (the third, on coefficients, is inline in input_structure.hpp): the InputStructure
// of an exact factor that preparation retained, an IntegerPolynomial (every coefficient exact and real) or Gaussian
// rationals (real where the imaginary part is zero), observed in increasing index order (InputStructure::observe).
// State: none. Caller of both: SolveState::prepared_structure (solver.cpp), itself called by tests and benchmarks.
#include "numeric/exact/prepared_structure.hpp"

namespace certpolysolve::detail {
InputStructure input_structure(const IntegerPolynomial &polynomial) {
    InputStructure result;
    for (slong j = 0; j < fmpz_poly_length(polynomial.value()); ++j)
        result.observe(static_cast<size_t>(j), true, true, fmpz_is_zero(polynomial.value()->coeffs + j));
    return result;
}
InputStructure input_structure(std::span<const GaussianRational> polynomial) {
    InputStructure result;
    for (size_t j = 0; j < polynomial.size(); ++j)
        result.observe(j, true, fmpq_is_zero(polynomial[j].imag.value()), polynomial[j].zero());
    return result;
}
} // namespace certpolysolve::detail
