// The first solve: every root of x^3 - 2x - 5.
#include <certpolysolve/certpolysolve.hpp>
#include <print>

int main() {
    certpolysolve::Polynomial p{-5, -2, 0, 1}; // constant term first: -5 - 2x + 0x^2 + x^3
    certpolysolve::Result result = certpolysolve::solve(p);

    for (const certpolysolve::Root &root : result.roots()) std::println("{}", root);
}
