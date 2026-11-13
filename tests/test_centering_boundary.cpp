// Exact centering must preserve the original polynomial on either side of
// its degree crossover. Rebuild the translated polynomial from rational
// planted roots with the independent test arithmetic; do not use a second
// Taylor shift as the oracle. This exercises preparation, not a root solve.
#include "solve/factor_state.hpp"
#include "solve/policies.hpp"
#include "numeric/exact/centering.hpp"
#include "winding_math.hpp"
#include <certpolysolve/certpolysolve.hpp>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace d = certpolysolve::detail;
using reference::Complex;
using reference::ComplexPolynomial;
using reference::Rational;
using reference::require;

namespace {
Rational rational(const d::Rational &value) {
    return Rational(d::integer_string(fmpq_numref(value.value()))) /
           Rational(d::integer_string(fmpq_denref(value.value())));
}

Complex value(const d::GaussianRational &coefficient) {
    return {rational(coefficient.real), rational(coefficient.imag)};
}

ComplexPolynomial coefficients(const d::FactorState &factor, bool translated) {
    ComplexPolynomial result;
    if (factor.kind == d::FactorState::Real) {
        const auto &polynomial = translated ? factor.centered : factor.integer;
        d::Integer coefficient;
        for (slong j = 0; j <= fmpz_poly_degree(polynomial.value()); ++j) {
            fmpz_poly_get_coeff_fmpz(coefficient.value(), polynomial.value(), j);
            result.emplace_back(Rational(d::integer_string(coefficient.value())));
        }
    } else {
        const auto &polynomial = translated ? factor.centered_gaussian : factor.gaussian;
        for (const auto &coefficient : polynomial) result.push_back(value(coefficient));
    }
    return result;
}

void import_source(d::FactorState &factor, const ComplexPolynomial &polynomial, bool gaussian) {
    factor.precision = 128;
    if (gaussian) {
        factor.kind = d::FactorState::Gaussian;
        factor.gaussian.resize(polynomial.size());
        for (size_t j = 0; j < polynomial.size(); ++j) {
            require(!fmpq_set_str(factor.gaussian[j].real.value(), polynomial[j].real.str().c_str(), 10),
                    "real coefficient import failed");
            require(!fmpq_set_str(factor.gaussian[j].imag.value(), polynomial[j].imag.str().c_str(), 10),
                    "imaginary coefficient import failed");
        }
    } else {
        factor.kind = d::FactorState::Real;
        reference::Polynomial primitive;
        for (const auto &coefficient : polynomial) {
            require(coefficient.imag == 0, "real fixture has an imaginary coefficient");
            primitive.push_back(coefficient.real);
        }
        reference::make_primitive(primitive);
        d::Integer coefficient;
        for (size_t j = 0; j < primitive.size(); ++j) {
            require(boost::multiprecision::denominator(primitive[j]) == 1,
                    "reference primitive coefficient is not an integer");
            require(!fmpz_set_str(coefficient.value(), primitive[j].str().c_str(), 10),
                    "integer coefficient import failed");
            fmpz_poly_set_coeff_fmpz(factor.integer.value(), static_cast<slong>(j), coefficient.value());
        }
    }
}

void check(slong degree, bool gaussian) {
    std::vector<Complex> roots;
    // A separated rational cloud away from the original origin. Non-dyadic
    // translations exercise the rounded proposal followed by an exact shift.
    const Complex center(Rational(5) / 3, gaussian ? Rational(-7) / 5 : Rational(0));
    for (slong j = 0; j < degree; ++j) {
        const Rational real_offset = Rational(2 * j + 1 - degree) / 4096;
        const Rational imag_offset = gaussian ? Rational(j % 3 - 1) / 8192 : Rational(0);
        roots.push_back(center + Complex(real_offset, imag_offset));
    }
    const auto original = reference::from_roots(roots);
    d::FactorState factor;
    import_source(factor, original, gaussian);
    const auto retained_source = coefficients(factor, false);
    d::ComplexPolynomial enclosure;
    if (gaussian) d::enclose_gaussian(enclosure.value(), factor.gaussian, factor.precision);
    else acb_poly_set_fmpz_poly(enclosure.value(), factor.integer.value(), factor.precision);

    cps_options options;
    options.bits = 53;
    options.max_work_bits = 512;
    auto policy = d::maintained_method_policy(CPS_METHOD_SECULAR);
    d::SearchStatistics statistics{};
    d::prepare_secular_center(factor, enclosure.value(), options, statistics, policy);

    const auto after_source = coefficients(factor, false);
    require(after_source.size() == retained_source.size(), "centering changed the source degree");
    for (size_t j = 0; j < retained_source.size(); ++j)
        require(reference::zero(after_source[j] - retained_source[j]),
                "centering changed an original source coefficient");

    if (degree == 63) {
        require(!factor.have_center && !factor.center_rejected && statistics.taylor_shifts == 0,
                "below-boundary fixture no longer exercises skipped centering");
        return;
    }
    require(factor.have_center && statistics.taylor_shifts > 0,
            "at/above-boundary fixture did not exercise accepted centering");
    const auto actual_center = value(factor.center);
    require(!reference::zero(actual_center), "accepted centering has a zero translation");
    for (auto &root : roots) root = root - actual_center;
    const auto expected = reference::from_roots(roots);
    const auto actual = coefficients(factor, true);
    require(actual.size() == expected.size() && !reference::zero(actual.back()),
            "centering changed the translated polynomial degree");
    // Exact preparation may remove a common scalar. Comparing every
    // coefficient relative to the leading one preserves that freedom.
    for (size_t j = 0; j < expected.size(); ++j)
        require(reference::zero(actual[j] * expected.back() - expected[j] * actual.back()),
                "centering disagrees with the independent translated-root expansion");
}
} // namespace

int main() {
    certpolysolve::ThreadScope lifetime;
    try {
        for (const slong degree : {63, 64, 65})
            for (const bool gaussian : {false, true}) check(degree, gaussian);
        std::cout << "centering degree-boundary exact coefficient checks passed\n";
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
