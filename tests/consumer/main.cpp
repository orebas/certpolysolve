// The installed package from C++ (the 0.5 interface): a double root, a complex double root, a refined solve, a
// secular equation every way, and each method forced.
#include <certpolysolve/certpolysolve.hpp>

int main() {
    certpolysolve::ThreadScope thread_lifetime;
    const auto result = certpolysolve::solve(certpolysolve::Polynomial{1, -2, 1});
    if (!result.ok() || result.roots()[0].multiplicity != 2) return 1;
    certpolysolve::Polynomial complex;
    complex.set(0, -1).set(1, certpolysolve::Number::complex(0, -2)).set(2, 1); // (x - i)^2
    const auto repeated = certpolysolve::solve(complex);
    if (!repeated.ok() || repeated.roots()[0].multiplicity != 2 || repeated.engine_counters().precision_attempts < 1 ||
        repeated.statistics().regions != 1 || repeated.statistics().roots != 2) return 1;
    const auto initial = certpolysolve::solve(certpolysolve::Polynomial{-2, 0, 1});
    const auto refined = certpolysolve::solve(certpolysolve::Polynomial{-2, 0, 1}, {.bits = 256});
    const certpolysolve::SecularEquation secular{{.residue = "2/3", .pole = 0}, {.residue = "-2/3", .pole = 3}};
    const auto secular_result = certpolysolve::solve(secular);
    if (!secular_result.ok() || secular_result.roots().size() != 2) return 1;
    if (secular_result.statistics().total_ns == 0) return 1;
    (void)secular_result.components();
    if (!certpolysolve::real_roots(secular).ok() ||
        !certpolysolve::solve(secular, {.goal = certpolysolve::Goal::clusters("1e-8")}).ok())
        return 1;
    for (const auto method : {certpolysolve::Method::classical, certpolysolve::Method::secular, certpolysolve::Method::piecewise}) {
        certpolysolve::Options selected;
        selected.method = method;
        selected.goal = certpolysolve::Goal::approximation();
        selected.bits = 80;
        const auto solved = certpolysolve::solve(certpolysolve::Polynomial{-1, 0, 0, 1}, selected);
        if (!solved.ok() || solved.roots().size() != 3) return 1;
    }
    return initial.ok() && refined.ok() && refined.complete() ? 0 : 1;
}
