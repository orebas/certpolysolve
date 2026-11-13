// The first solve of docs/api.md, section 1, in C++. The suite builds and runs it (guide_first_solve_cpp).
#include <certpolysolve/certpolysolve.hpp>
#include <chrono>
#include <print>

using namespace std::chrono_literals;

int main() {
    const certpolysolve::Polynomial p{-2, 0, 1}; // x^2 - 2, ascending
    const certpolysolve::Result result = certpolysolve::solve(p, {.bits = 128, .time_limit = 10s});
    if (!result.ok()) { // an error, or a stop with proved regions short of the goal
        std::println(stderr, "{}: {}", result.stop(), result.message());
        return 1;
    }
    for (const certpolysolve::Root &root : result.roots())
        std::println("{} holds {} root(s), centre {:.40}", root.box.real, root.count, root);
}
