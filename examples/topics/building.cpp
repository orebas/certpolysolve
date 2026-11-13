// Building a polynomial: with exact arithmetic, from its roots, by a recurrence, term by term; and reading it back.
#include <certpolysolve/certpolysolve.hpp>
#include <print>
#include <vector>

int main() {
    using certpolysolve::Polynomial;
    const Polynomial x = Polynomial::x();

    // Arithmetic is exact: +, -, * and pow, with numbers mixed in. A fraction is text: 1/3 in C++ is 0.
    const Polynomial p = pow(x - 1, 3) * (x + 2) - "1/3";
    std::println("p = {}", p);
    std::println("p' = {}", p.derivative());

    // Wilkinson's polynomial (x - 1)(x - 2)...(x - 20), as a product and from its roots.
    Polynomial w = 1;
    std::vector<int> roots;
    for (int k = 1; k <= 20; ++k) {
        w *= x - k;
        roots.push_back(k);
    }
    std::println("the same both ways: {}", w == Polynomial::from_roots(roots));
    std::println("degree {}; the coefficient of x^19 is {}, the constant term {}", w.degree(), w[19], w[0]);

    // Chebyshev's T_6 by its recurrence T_(k+1) = 2x T_k - T_(k-1).
    Polynomial previous = 1, current = x;
    for (int k = 1; k < 6; ++k) {
        Polynomial next = 2 * x * current - previous;
        previous = current;
        current = next;
    }
    std::println("T_6 = {}", current);
    std::print("its coefficients, constant term first:");
    for (const certpolysolve::Coefficient &c : current.coefficients()) std::print(" {}", c);
    std::println("");

    // Term by term: coefficients you do not set are zero, and only the nonzero terms are stored, at any degree.
    Polynomial sparse;
    sparse.set(1000, 1).set(7, -3).set(0, "1/2");
    std::print("{}: degree {}, terms", sparse, sparse.degree());
    for (const Polynomial::Term &term : sparse.terms()) std::print(" ({}, {})", term.exponent, term.coefficient);
    std::println("");
}
