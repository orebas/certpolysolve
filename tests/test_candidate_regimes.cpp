// The sparse candidate evaluator must compute the same polynomial and its
// derivatives on both sides of its support crossover. At 0 and the fourth
// roots of unity, small Gaussian-integer coefficients give exact integer
// answers: no floating-point error tolerance or solver trajectory is needed.
#include "numeric/candidates/candidate_jets.hpp"
#include <array>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace m = certpolysolve::detail::machine;
using Candidate = m::Complex<double>;

namespace {
struct Integer {
    int64_t real = 0, imag = 0;
    Integer &operator+=(Integer other) { real += other.real; imag += other.imag; return *this; }
    friend Integer operator*(Integer a, Integer b) {
        return {a.real * b.real - a.imag * b.imag, a.real * b.imag + a.imag * b.real};
    }
    friend Integer operator*(Integer a, int64_t scale) { return {a.real * scale, a.imag * scale}; }
};
void require(bool value, const char *message) {
    if (!value) throw std::runtime_error(message);
}

Integer power(Integer z, size_t exponent) {
    Integer result{1, 0};
    for (size_t i = 0; i < exponent; ++i) result = result * z;
    return result;
}

// Differentiate monomials, then sum their values. This reference uses neither
// a Horner recurrence nor the evaluator's sparse gap/product-rule recurrence.
Integer exact_derivative(const std::vector<Integer> &coefficients, Integer z, unsigned order) {
    Integer result{};
    for (size_t exponent = order; exponent < coefficients.size(); ++exponent) {
        int64_t factor = 1;
        for (unsigned k = 0; k < order; ++k) factor *= static_cast<int64_t>(exponent - k);
        result += coefficients[exponent] * factor * power(z, exponent - order);
    }
    return result;
}

void check(const std::vector<Integer> &coefficients) {
    std::vector<Candidate> converted;
    for (const auto &coefficient : coefficients)
        converted.emplace_back(static_cast<double>(coefficient.real), static_cast<double>(coefficient.imag));
    const m::PreparedJets<double, 2> full(converted);
    const m::PreparedJets<double, 1> first(converted);
    for (const Integer z : {Integer{0, 0}, Integer{1, 0}, Integer{-1, 0}, Integer{0, 1}, Integer{0, -1}}) {
        const Candidate point{static_cast<double>(z.real), static_cast<double>(z.imag)};
        const auto value = full.evaluate(point);
        const auto reduced = first.evaluate(point);
        const std::array<Candidate, 3> actual{value.value, value.derivative, value.second_derivative};
        for (unsigned order = 0; order < actual.size(); ++order) {
            const auto expected = exact_derivative(coefficients, z, order);
            require(actual[order].re == static_cast<double>(expected.real) &&
                        actual[order].im == static_cast<double>(expected.imag),
                    "candidate jet disagrees with the exact monomial derivative sum");
        }
        require(reduced.value == value.value && reduced.derivative == value.derivative,
                "first-order jet changed the polynomial or its derivative");
    }
}
} // namespace

int main() {
    try {
        check({});
        check({Integer{3, -2}});
        // Adjacent powers-of-two degrees also cross the long-gap power-table
        // sizes and the integer rounding in the quarter-support decision.
        for (const size_t degree : {31, 32, 33, 63, 64, 65}) {
            for (const int offset : {-1, 0, 1}) {
                const size_t support = static_cast<size_t>(static_cast<int>(degree / 4) + offset);
                for (const bool constant : {false, true}) {
                    for (const bool long_gap : {false, true}) {
                        std::vector<Integer> coefficients(degree + 1);
                        for (size_t j = 0; j < support; ++j) {
                            const size_t first_exponent = constant ? 0 : 1;
                            const size_t exponent = long_gap
                                ? (j == 0 ? first_exponent : degree - (support - 1 - j))
                                : first_exponent + j * (degree - first_exponent) / (support - 1);
                            const int64_t sign = j % 2 ? -1 : 1;
                            coefficients[exponent] = {sign * (1 + static_cast<int64_t>(j % 3)),
                                                      static_cast<int64_t>(j % 5) - 2};
                        }
                        // All intermediate integer components are far below
                        // 2^53 at these degrees, including the second derivative.
                        check(coefficients);
                    }
                }
            }
        }
        std::cout << "candidate support-regime exact derivative checks passed\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
