// Coefficients known only within bounds: about x^2 - 3x + 2, each coefficient within 1e-9.
#include <certpolysolve/certpolysolve.hpp>
#include <print>

int main() {
    using certpolysolve::ball;
    certpolysolve::Polynomial p{ball(2, "1e-9"), ball(-3, "1e-9"), ball(1, "1e-9")}; // value, then radius
    std::println("{}", p);

    // The boxes hold the roots of EVERY polynomial within those bounds, so ask for boxes of a size, not for digits.
    certpolysolve::Result result = certpolysolve::solve(p, {.goal = certpolysolve::Goal::clusters("1e-6")});
    std::println("{}", result.ok() ? "goal met" : "stopped");
    for (const certpolysolve::Root &root : result.roots())
        std::println("  real part in {:.12}, count {}", root.box.real, root.count);
}
