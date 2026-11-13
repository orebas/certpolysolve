#include "numeric/arithmetic/polynomial_arithmetic.hpp"
#include <certpolysolve/certpolysolve.hpp>
#include "numeric/exact/input_structure.hpp"
#include "internal.hpp"
#include <iostream>
#include <stdexcept>
#include <type_traits>

namespace d = certpolysolve::detail;
static void require(bool value, const char *message) {
    if (!value) throw std::runtime_error(message);
}
static_assert(std::is_trivially_copyable_v<d::InputStructure>);

int main() {
    certpolysolve::ThreadScope lifetime;
    try {
        certpolysolve::Polynomial zero;
        zero.set(80, certpolysolve::Number("0"));
        const auto z = d::input_structure(certpolysolve::detail::dense(*zero.get()).coefficients);
        require(z.start == 0 && z.length == 0 && z.support == 0 && z.reduced_degree() == 0,
                "zero polynomial has fictitious support");

        // x^5 (2 + 3 i x^12 - x^36): zero lifting and power lifting are distinct.
        certpolysolve::Polynomial p;
        p.set(5, certpolysolve::Number("2"))
            .set(17, certpolysolve::Number::complex(certpolysolve::Number("0"), certpolysolve::Number("3")))
            .set(41, certpolysolve::Number("-1"))
            .set(90, certpolysolve::Number("0"));
        const auto original = d::input_structure(certpolysolve::detail::dense(*p.get()).coefficients);
        require(original.start == 5 && original.length == 42 && original.support == 3 &&
                    original.exact && !original.real && original.support_stride == 12 &&
                    original.degree() == 36 && original.reduced_degree() == 3 &&
                    original.unshifted_power() == 1,
                "zero and power factors were conflated");
        // Metadata owns no references to the caller's mutable input.
        p.set(6, certpolysolve::Number("1"));
        require(original.support_stride == 12 && d::input_structure(certpolysolve::detail::dense(*p.get()).coefficients).support_stride == 1,
                "input mutation changed an existing value snapshot");

        certpolysolve::Polynomial giant_power;
        giant_power.set(0, certpolysolve::Number("-1")).set(65537, certpolysolve::Number("1"));
        const auto power = d::input_structure(certpolysolve::detail::dense(*giant_power.get()).coefficients);
        require(power.unshifted_power() == 65537 && power.reduced_degree() == 1 && power.real,
                "large binomial was treated as a large irreducible numerical problem");

        certpolysolve::Polynomial compressed;
        compressed.set(0, certpolysolve::Number("1@-1000000000000000000000000000000000000000"))
                  .set(8, certpolysolve::Number("1e1000000000000000000000000000000000000000"));
        const auto huge = d::input_structure(certpolysolve::detail::dense(*compressed.get()).coefficients);
        require(huge.exact && huge.real && huge.unshifted_power() == 8,
                "compressed exponents changed algebraic support");
        compressed.set(12, certpolysolve::ball(certpolysolve::Number("0"), certpolysolve::Number("0"), certpolysolve::Number("1")));
        const auto uncertain = d::input_structure(certpolysolve::detail::dense(*compressed.get()).coefficients);
        require(!uncertain.exact && !uncertain.real && uncertain.length == 13 &&
                    uncertain.support == 3 && uncertain.unshifted_power() == 1,
                "uncertain zero midpoint was trimmed or treated as exact");

        certpolysolve::Polynomial constant;
        constant.set(0, certpolysolve::Number("7")).set(40, certpolysolve::Number("0"));
        const auto c = d::input_structure(certpolysolve::detail::dense(*constant.get()).coefficients);
        require(c.length == 1 && c.degree() == 0 && c.unshifted_power() == 1,
                "constant was treated as a power reduction");

        // Exercise the cached identity through both maintained global methods.
        // The original polynomial has twelve unit roots, not just one root of q.
        certpolysolve::Polynomial roots;
        roots.set(0, certpolysolve::Number("-1")).set(12, certpolysolve::Number("1")).set(20, certpolysolve::Number("0"));
        for (auto method : {certpolysolve::Method::classical, certpolysolve::Method::secular}) {
            certpolysolve::Options options;
            options.method = method;
            options.goal = certpolysolve::Goal::approximation();
            options.bits = 53;
            auto answer = certpolysolve::solve(roots, options);
            require(answer.goal_met() && answer.complete() && answer.roots().size() == 12,
                    "cached power metadata lost original roots");
        }
        std::cout << "input structure checks passed\n";
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
