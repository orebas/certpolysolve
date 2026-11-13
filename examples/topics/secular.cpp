// Eigenvalues of D + rho z z^T, D diagonal: the roots of the secular equation sum_j a_j / (x - b_j) - 1 = 0,
// with a_j = rho z_j^2 and b_j = d_j. Here D = diag(1, 2, 3, 4), z = (1, 1/2, 1/3, 1/4), rho = 1.
#include <certpolysolve/certpolysolve.hpp>
#include <print>
#include <string>
#include <vector>

int main() {
    const std::vector<std::string> d = {"1", "2", "3", "4"};
    const std::vector<std::string> a = {"1", "1/4", "1/9", "1/16"}; // rho z_j^2

    std::vector<certpolysolve::SecularEquation::Term> terms;
    for (size_t j = 0; j < d.size(); ++j) terms.push_back({.residue = a[j], .pole = d[j]});
    const certpolysolve::SecularEquation s(terms);

    const certpolysolve::Result result = certpolysolve::real_roots(s);
    std::println("{}", result.ok() ? "goal met" : "stopped");
    for (const certpolysolve::Root &root : result.roots()) std::println("  {:.17}", root);
}
