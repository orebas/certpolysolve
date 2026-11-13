// A coefficient ball denotes every polynomial in that box, not just its
// midpoint. Check exact corner and interior witnesses with independent rational
// argument counts. Finite sampling strengthens tests; it is not the interval
// certificate's proof for the whole family.
#include "exact_winding_api.hpp"
#include <iostream>

using namespace reference;

namespace {
struct CoefficientBox {
    Complex midpoint;
    Rational real_radius, imaginary_radius;
};
struct Family {
    std::vector<CoefficientBox> coefficients;
    certpolysolve::Polynomial input() const {
        certpolysolve::Polynomial result;
        for (size_t j = 0; j < coefficients.size(); ++j) {
            const auto &c = coefficients[j];
            result.set(j, certpolysolve::ball(
                certpolysolve::Number::complex(certpolysolve::Number(c.midpoint.real.str()),
                                          certpolysolve::Number(c.midpoint.imag.str())),
                certpolysolve::Number(c.real_radius.str()), certpolysolve::Number(c.imaginary_radius.str())));
        }
        require(result.ok(), "exact coefficient-family input was refused");
        return result;
    }
    std::vector<ComplexPolynomial> witnesses() const {
        ComplexPolynomial midpoint;
        struct Axis { size_t coefficient; bool imaginary; Rational radius; };
        std::vector<Axis> axes;
        for (size_t j = 0; j < coefficients.size(); ++j) {
            const auto &c = coefficients[j];
            midpoint.push_back(c.midpoint);
            if (c.real_radius != 0) axes.push_back({j, false, c.real_radius});
            if (c.imaginary_radius != 0) axes.push_back({j, true, c.imaginary_radius});
        }
        require(axes.size() <= 6, "test family has too many corners for its bounded oracle");
        const auto move = [](ComplexPolynomial &p, const Axis &axis, const Rational &amount) {
            auto &c = p[axis.coefficient];
            (axis.imaginary ? c.imag : c.real) += amount;
        };
        std::vector<ComplexPolynomial> result{midpoint};
        for (size_t corner = 0; corner < (size_t{1} << axes.size()); ++corner) {
            auto witness = midpoint;
            for (size_t j = 0; j < axes.size(); ++j)
                move(witness, axes[j], (corner & (size_t{1} << j)) ? axes[j].radius : -axes[j].radius);
            result.push_back(std::move(witness));
        }
        for (const auto &axis : axes)
            for (const int sign : {-1, 1}) {
                auto witness = midpoint;
                move(witness, axis, sign * axis.radius / 2);
                result.push_back(std::move(witness));
            }
        return result;
    }
};

static void check_family(const Family &family, size_t &checked) {
    const auto input = family.input();
    const auto witnesses = family.witnesses();
    const uint64_t degree = family.coefficients.size() - 1;
    for (const uint64_t bits : {12, 80}) {
        for (const uint64_t work : {1, 128, 2048}) {
            certpolysolve::Options options;
            options.goal = certpolysolve::Goal::approximation();
            options.bits = bits;
            options.max_bits = 256;
            options.work_budget = work;
            const auto result = certpolysolve::solve(input, options);
            require(result.get() && result.stop() != certpolysolve::Stop::error, "family solve returned an API error");
            require(result.degree() == static_cast<int64_t>(degree), "nonzero leading family lost its degree");
            check_disjoint(result.roots());
            uint64_t total = 0;
            for (const auto &root : result.roots()) {
                for (const auto &witness : witnesses) {
                    const auto count = rectangle_count(witness, root);
                    require(count.has_value(), "independent witness contour count was inconclusive");
                    require(*count == static_cast<int>(root.count),
                            "a published family region miscounted an exact coefficient witness");
                    ++checked;
                }
                if (result.goal_met()) check_accuracy(root, static_cast<int>(bits));
                total += root.count;
            }
            if (result.complete()) require(total == degree, "complete coefficient-family cover lost roots");
            if (bits == 12 && work == 2048)
                require(result.complete() && result.goal_met(), "loose family goal did not finish");
        }
    }
}
} // namespace

int main() {
    certpolysolve::ThreadScope lifetime;
    try {
        const Rational error = power2(-40);
        size_t checked = 0;
        // A complex linear family: both the numerator and leading coefficient vary.
        check_family({{{Complex(-2, 1), error, error}, {Complex(1), error, error}}}, checked);
        // The exact midpoint has roots -1+i and 2-i. Every coefficient varies independently.
        check_family({{{Complex(-1, 3), error, error}, {Complex(-1), error, error},
                       {Complex(1), error, error}}}, checked);
        // The family includes a double root and nearby split roots: a counted cluster is required.
        check_family({{{Complex(1), error, error}, {Complex(-2), 0, 0}, {Complex(1), 0, 0}}}, checked);
        certpolysolve::Polynomial changing_degree;
        changing_degree.set(0, 1).set(1, certpolysolve::ball(0, 1));
        require(changing_degree.ok(), "degree-changing family fixture was refused");
        certpolysolve::Options options;
        options.max_bits = 64;
        const auto refused = certpolysolve::solve(changing_degree, options);
        require(refused.get() && cps_result_termination(refused.get()) == CPS_INDETERMINATE_DEGREE &&
                !refused.complete() && !refused.goal_met() && refused.roots().empty(),
                "a family with varying degree was given a finite complete root cover");
        require(checked > 0, "the coefficient-family witness oracle never ran");
        std::cout << checked << " exact coefficient-witness region counts passed\n";
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
