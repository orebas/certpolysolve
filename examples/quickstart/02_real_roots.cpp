// Real roots: all of them, those in an interval, the positive ones, and only how many.
#include <certpolysolve/certpolysolve.hpp>
#include <print>
#include <string_view>

using certpolysolve::Rounding;

void print(std::string_view title, const certpolysolve::Result &result) {
    std::println("{}: {}", title, result.roots().size());
    for (const certpolysolve::Root &root : result.roots()) // the interval's ends as doubles, rounded outward
        std::println("  [{}, {}]", root.box.real.lo.to_double(Rounding::down), root.box.real.hi.to_double(Rounding::up));
}

int main() {
    certpolysolve::Polynomial quintic{-1, -1, 0, 0, 0, 1}; // x^5 - x - 1: one real root of five
    print("real roots of x^5 - x - 1", certpolysolve::real_roots(quintic));

    certpolysolve::Polynomial t5{-0.5, 5, 0, -20, 0, 16}; // 16x^5 - 20x^3 + 5x - 1/2 (the double -0.5 is exactly -1/2)
    print("positive roots of T5(x) - 1/2", certpolysolve::real_roots(t5, certpolysolve::at_least(0)));
    print("roots in [-1/2, 1/2]", certpolysolve::real_roots(t5, certpolysolve::between("-1/2", "1/2"))); // closed: ends count

    // each root alone in an interval, not refined: the cheapest query
    const certpolysolve::Result count = certpolysolve::real_roots(t5, {.goal = certpolysolve::Goal::isolation()});
    std::println("T5(x) - 1/2 has {} distinct real roots", count.roots().size());
}
