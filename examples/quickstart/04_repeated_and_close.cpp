// Repeated roots and roots very close together.
#include <certpolysolve/certpolysolve.hpp>
#include <print>
#include <string_view>

void print(std::string_view title, const certpolysolve::Result &result) {
    std::println("{} ({})", title, result.ok() ? "goal met" : "stopped");
    for (const certpolysolve::Root &root : result.roots()) {
        std::print("  real {}, imaginary {}: count {}, ", root.box.real.mid().format(12), root.box.imag.mid().format(3),
                   root.count);
        if (root.multiplicity) std::println("multiplicity {}", *root.multiplicity);
        else std::println("multiplicity not known");
    }
}

int main() {
    // (x - 1)^3 (x + 2) with exact coefficients: the library proves that 1 is a triple root.
    certpolysolve::Polynomial exact{-2, 5, -3, -1, 1};
    print("(x - 1)^3 (x + 2)", certpolysolve::solve(exact));

    // (x - 0.1)^3 with its coefficients rounded to doubles: no longer a cube, but three distinct roots close to 0.1.
    const double a = 0.1;
    certpolysolve::Polynomial rounded{-a * a * a, 3 * a * a, -3 * a, 1};
    print("(x - 0.1)^3 in doubles, every root alone", certpolysolve::solve(rounded));

    // Often one box and its count is the useful answer: clusters of radius at most 1e-4.
    print("the same, as clusters of radius 1e-4", certpolysolve::solve(rounded, {.goal = certpolysolve::Goal::clusters("1e-4")}));
}
