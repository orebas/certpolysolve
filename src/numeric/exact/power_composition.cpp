// Exact power compositions p(x) = q(x^k), k > 1, q(0) != 0, found from which coefficients vanish (coefficient_power,
// prepare_power_composition, reduce_power_composition).
//
// k is the gcd of the exponents of the nonzero terms, sought only when the constant term is nonzero, so that zero,
// where x -> x^k ramifies, is no root. coefficient_power and prepare_power_composition read the stored, exact
// coefficients without expanding them; reduce_power_composition rewrites an integer or Gaussian-rational factor.
// Lifting a proved cover of q's roots to p's belongs to certification (certify/power_cover.cpp).
//
// State: none.
// Entry points: prepare_power_composition (SolveState::Impl::preprocess in solver.cpp; prepare_squarefree_input in
// exact_preprocessing.cpp); reduce_power_composition (preprocess); coefficient_power (only from here).
// Design notes: docs/mathematics.md ("C2 source and representation contracts": the p(z) = q(z^k) reduction).
#include "numeric/exact/power_composition.hpp"
#include <algorithm>
#include <numeric>

namespace certpolysolve::detail {
ulong coefficient_power(std::span<const Coefficient> coefficients) {
    if (coefficients.size() < 3 || coefficients.front().zero() ||
        std::any_of(coefficients.begin(), coefficients.end(), [](const auto &c) { return !c.exact(); }))
        return 1;
    ulong power = static_cast<ulong>(coefficients.size() - 1);
    for (size_t i = 1; i + 1 < coefficients.size() && power > 1; ++i)
        if (!coefficients[i].zero())
            power = std::gcd(power, static_cast<ulong>(i));
    return power;
}

std::optional<ExactPowerComposition> prepare_power_composition(std::span<const Coefficient> coefficients) {
    const auto power = coefficient_power(coefficients);
    if (power == 1) return std::nullopt;
    ExactPowerComposition result;
    result.power = power;
    result.polynomial.coefficients.reserve((coefficients.size() - 1) / power + 1);
    for (size_t i = 0; i < coefficients.size(); i += power)
        result.polynomial.coefficients.push_back(coefficients[i]);
    return result;
}

ulong reduce_power_composition(IntegerPolynomial &polynomial) {
    const slong degree = fmpz_poly_degree(polynomial.value());
    if (degree < 2 || fmpz_is_zero(polynomial.value()->coeffs))
        return 1;
    ulong power = static_cast<ulong>(degree);
    for (slong i = 1; i < degree && power > 1; ++i)
        if (!fmpz_is_zero(polynomial.value()->coeffs + i))
            power = std::gcd(power, static_cast<ulong>(i));
    if (power == 1)
        return 1;
    IntegerPolynomial reduced;
    for (slong j = 0; j <= degree / static_cast<slong>(power); ++j)
        fmpz_poly_set_coeff_fmpz(reduced.value(), j, polynomial.value()->coeffs + j * static_cast<slong>(power));
    polynomial = std::move(reduced);
    return power;
}

ulong reduce_power_composition(GaussianPolynomial &polynomial) {
    if (polynomial.size() < 3 ||
        (fmpq_is_zero(polynomial[0].real.value()) && fmpq_is_zero(polynomial[0].imag.value())))
        return 1;
    const auto degree = polynomial.size() - 1;
    size_t power = degree;
    for (size_t i = 1; i < degree && power > 1; ++i)
        if (!fmpq_is_zero(polynomial[i].real.value()) || !fmpq_is_zero(polynomial[i].imag.value()))
            power = std::gcd(power, i);
    if (power == 1)
        return 1;
    GaussianPolynomial reduced(degree / power + 1);
    for (size_t j = 0; j < reduced.size(); ++j)
        reduced[j] = polynomial[j * power];
    polynomial = std::move(reduced);
    return static_cast<ulong>(power);
}

} // namespace certpolysolve::detail
