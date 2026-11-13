// Programs: give the polynomial by the steps that compute it, not by its coefficients.
#include <certpolysolve/certpolysolve.hpp>
#include <algorithm>
#include <cmath>
#include <numbers>
#include <print>
#include <string_view>
#include <vector>

void summary(std::string_view title, const certpolysolve::Result &result) {
    std::println("{}: degree {}, {} boxes, {}", title, result.degree(), result.roots().size(),
                 result.ok() ? "goal met" : "stopped");
}

int main() {
    // The Mandelbrot polynomial p_8(c): p_1 = c, p_(k+1) = p_k^2 + c. Degree 128.
    auto mandelbrot = certpolysolve::Program::from([](auto c) {
        auto p = c;
        for (int k = 1; k < 8; ++k) p = p * p + c;
        return p;
    });
    summary("Mandelbrot p_8", certpolysolve::solve(mandelbrot));

    // Chebyshev's T_50 by its recurrence T_(k+1) = 2x T_k - T_(k-1).
    auto chebyshev = certpolysolve::Program::from([](auto x) {
        auto previous = x, current = 2 * x * x - 1; // T_1 and T_2
        for (int k = 2; k < 50; ++k) {
            auto next = 2 * x * current - previous;
            previous = current;
            current = next;
        }
        return current;
    });
    summary("Chebyshev T_50", certpolysolve::solve(chebyshev));

    // A tridiagonal matrix's characteristic polynomial det(T - x I), by the three-term recurrence. Here the 1D
    // Laplacian of size 40 (2 on the diagonal, -1 beside it), whose eigenvalues are 2 - 2 cos(k pi / 41).
    const int n = 40;
    const std::vector<double> diagonal(n, 2.0), beside(n - 1, -1.0);
    auto laplacian = certpolysolve::Program::from([&](auto x) {
        auto previous = diagonal[0] - x;                                     // the leading 1 x 1 minor
        auto current = (diagonal[1] - x) * previous - beside[0] * beside[0]; // the leading 2 x 2 minor
        for (int k = 2; k < n; ++k) {
            auto next = (diagonal[k] - x) * current - beside[k - 1] * beside[k - 1] * previous;
            previous = current;
            current = next;
        }
        return current;
    });
    const certpolysolve::Result eigenvalues = certpolysolve::solve(laplacian);
    summary("Laplacian of size 40", eigenvalues);
    std::vector<double> computed;
    for (const certpolysolve::Root &root : eigenvalues.roots())
        if (root.value) computed.push_back(root.value->real());
    std::ranges::sort(computed);
    double error = 0;
    for (int k = 1; k <= n; ++k)
        error = std::max(error, std::abs(computed[k - 1] - (2 - 2 * std::cos(k * std::numbers::pi / (n + 1)))));
    std::println("  largest difference from 2 - 2 cos(k pi / 41): {:.3g}", error);

    // Wilkinson's polynomial (x - 1)(x - 2)...(x - 30) as a product: its coefficients, up to 30!, never appear.
    auto wilkinson = certpolysolve::Program::from([](auto x) {
        auto p = x - 1;
        for (int k = 2; k <= 30; ++k) p *= x - k;
        return p;
    });
    // every root in its own box (a program's default allows close roots to share one)
    summary("Wilkinson W_30", certpolysolve::solve(wilkinson, {.goal = certpolysolve::Goal::accuracy()}));

    // Powers and exact fractions: (x^2 + x + 1)^5 - 1/3.
    auto power = certpolysolve::Program::from([](auto x) { return pow(x * x + x + 1, 5) - "1/3"; });
    summary("(x^2 + x + 1)^5 - 1/3", certpolysolve::solve(power));
}
