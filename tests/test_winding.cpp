#include "exact_winding_api.hpp"
#include <iostream>
#include <random>

using namespace reference;
namespace {
void checker_examples() {
    require(cauchy_index({1}, {-1, 2}) == 1, "simple positive Cauchy index failed");
    require(cauchy_index({-1}, {-1, 2}) == -1, "negative Cauchy index failed");
    require(cauchy_index({1}, {1, -4, 4}) == 0, "even pole has a nonzero Cauchy index");
    require(cauchy_index({1}, {-1, 6, -12, 8}) == 1, "odd high-order pole lost its index");
    const auto polynomial = from_roots({{-1}, {1}, {0, 1}, {0, -1}, {2}, {2}});
    require(rectangle_count(polynomial, -3, 3, -3, 3) == 6, "global winding count failed");
    require(rectangle_count(polynomial, Rational(-3) / 2, Rational(3) / 2, Rational(-1) / 2, Rational(1) / 2) == 2,
            "local winding count failed");
    require(rectangle_count(polynomial, Rational(-1) / 2, Rational(1) / 2, Rational(1) / 2, Rational(3) / 2) == 1,
            "complex root winding count failed");
    require(rectangle_count(polynomial, 4, 5, 4, 5) == 0, "empty rectangle count failed");
    require(!rectangle_count(polynomial, -2, 2, -2, 2), "boundary root accepted");
    require(!rectangle_count(polynomial, 2, 3, 0, 1), "vertex root accepted");
    require(rectangle_count(polynomial, -3, 3, 0, 0) == 4, "closed real segment count failed");
    require(rectangle_count(polynomial, 0, 0, -2, 2) == 2, "closed imaginary segment count failed");
    require(rectangle_count(polynomial, 2, 2, 0, 0) == 2, "exact point multiplicity failed");
    require(rectangle_count(polynomial, 3, 3, 0, 0) == 0, "empty point count failed");
    require(rectangle_count(ComplexPolynomial{{3, 7}}, -1, 1, -1, 1) == 0, "constant polynomial count failed");
}
void known_roots() {
    std::mt19937 generator(20260906);
    for (int sample = 0; sample < 80; ++sample) {
        std::vector<Complex> roots;
        for (int j = 0; j < 1 + sample % 8; ++j)
            roots.emplace_back(Rational(static_cast<int>(generator() % 33) - 16) / 4,
                               Rational(static_cast<int>(generator() % 33) - 16) / 4);
        auto p = from_roots(roots);
        if (sample % 2) for (auto& coefficient : p) coefficient = coefficient * Complex(-3, 7);
        const Rational left = Rational(static_cast<int>(generator() % 17) - 16) / 4;
        const Rational bottom = Rational(static_cast<int>(generator() % 17) - 16) / 4;
        const Rational right = left + Rational(1 + generator() % 24) / 4;
        const Rational top = bottom + Rational(1 + generator() % 24) / 4;
        int expected = 0;
        bool boundary = false;
        for (const auto& z : roots) {
            const bool inside = z.real >= left && z.real <= right && z.imag >= bottom && z.imag <= top;
            if (inside) {
                ++expected;
                boundary = boundary || z.real == left || z.real == right || z.imag == bottom || z.imag == top;
            }
        }
        const auto count = rectangle_count(p, left, right, bottom, top);
        require(boundary ? !count.has_value() : count == expected, "constructed-root winding oracle mismatch");
    }
}
void validate_solver_rectangles() {
    for (int degree : {3, 4, 5, 6}) {
        ComplexPolynomial p(static_cast<size_t>(degree) + 1);
        p[0] = {-2, -1}; p.back() = {1};
        certpolysolve::Options options;
        options.bits = 40;
        options.max_bits = 256;
        auto result = certpolysolve::solve(input(p), options);
        require(result.goal_met(), "complex winding regression did not solve");
        const auto& roots = result.roots();
        check_disjoint(roots);
        for (const auto& root : roots) {
            const auto count = rectangle_count(p, root);
            require(count && *count == static_cast<int>(root.count), "independent exact winding rejected exported root count");
        }
    }
    // Complex coefficients exclude the current exact-real preprocessing path.
    const auto repeated = from_roots({{1, 1}, {1, 1}, {-2, 1}, {-2, 1}});
    certpolysolve::Options options;
    options.goal = certpolysolve::Goal::clusters(certpolysolve::Number("1e-6"));
    options.bits = 40;
    options.max_bits = 256;
    auto result = certpolysolve::solve(input(repeated), options);
    require(result.goal_met() && result.roots().size() == 2, "complex repeated clusters failed");
    for (const auto& root : result.roots())
        require(rectangle_count(repeated, root) == static_cast<int>(root.count), "independent winding rejected repeated cluster");
}
}
int main() {
    certpolysolve::ThreadScope thread_lifetime;
    try {
        checker_examples();
        known_roots();
        validate_solver_rectangles();
        std::cout << "Independent exact Sturm--Cauchy rectangle counts passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
