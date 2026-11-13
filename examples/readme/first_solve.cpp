// The README's first solve: every root of x^2 - 2, printed as numbers.
#include <certpolysolve/certpolysolve.hpp>
#include <print>

int main() {
    const certpolysolve::Polynomial p{-2, 0, 1}; // constant term first: x^2 - 2
    std::println("p(x) = {}", p);
    const certpolysolve::Result result = certpolysolve::solve(p);
    if (!result.ok()) {
        std::println(stderr, "{}: {}", result.stop(), result.message());
        return 1;
    }
    for (const certpolysolve::Root &root : result.roots())
        std::println("{}", root); // the centre of a box proved to hold one root
}
