#include "solve/prepared_piecewise.hpp"
#include "solve/policies.hpp"
#include "numeric/arithmetic/polynomial_arithmetic.hpp"
#include "internal.hpp"
#include "exact_winding_api.hpp"
#include "methods/piecewise/piecewise_roots.hpp"
#include <iostream>

namespace d = certpolysolve::detail;
using namespace reference;

static void check_snapshot(const d::PreparedStructure &s, uint64_t degree,
                           uint64_t zeros, uint64_t numerical_degree) {
    require(s.input_degree == static_cast<int64_t>(degree) && s.accounted_degree() == degree,
            "prepared factors do not account for original degree");
    if (s.zero_order != zeros || s.numerical_degree() != numerical_degree)
        throw std::runtime_error("prepared factor mismatch: zero order " + std::to_string(s.zero_order) +
            " (expected " + std::to_string(zeros) + "), numerical degree " +
            std::to_string(s.numerical_degree()) + " (expected " + std::to_string(numerical_degree) + ")");
}

int main() {
    certpolysolve::ThreadScope lifetime;
    try {
        // Independently construct x^5 (x^4-(1+2i))^3 (x^4-(-2+i))^2.
        // Its two nonzero square-free factors are coprime degree-four factors;
        // all maintained methods additionally reduce each to degree one.
        const auto q = from_roots({Complex(1,2), Complex(1,2), Complex(1,2),
                                   Complex(-2,1), Complex(-2,1)});
        ComplexPolynomial p(26);
        for (size_t j = 0; j < q.size(); ++j) p[5 + 4*j] = q[j];
        auto input = reference::input(p);
        for (auto method : {CPS_METHOD_CLASSICAL, CPS_METHOD_SECULAR}) {
            d::SolveState state(certpolysolve::detail::dense(*input.get()), d::maintained_method_policy(method));
            require(!state.prepared_structure(), "unprepared state exposed factors");
            cps_options options;
            options.goal = CPS_ROOT_APPROXIMATION;
            options.work_units = 2; // Degree proof, then exact preparation.
            state.advance(&options);
            const auto saved = state.prepared_structure();
            require(saved.has_value(), "preparation did not commit");
            check_snapshot(*saved, 25, 5, 2);
            require(saved->factors.size() == 2, "square-free factor count lost");
            uint64_t multiplicities = 0;
            for (const auto &factor : saved->factors) {
                require(factor.power == 4 && factor.squarefree_proved &&
                        factor.representation == d::PreparedRepresentation::Gaussian,
                        "actual factor representation or power was misreported");
                multiplicities += factor.multiplicity;
            }
            require(multiplicities == 5, "exact multiplicities lost");
            // The whole solve is a fresh state (one advance a solve state since 18 October 2026).
            options.work_units = 10000;
            d::SolveState whole(certpolysolve::detail::dense(*input.get()), d::maintained_method_policy(method));
            auto answer = whole.advance(&options);
            require(answer.goal_met, "the solve after an inspected preparation failed");
            check_snapshot(*saved, 25, 5, 2); // Owned values survive later work.
        }
        d::PreparedPiecewiseState piecewise(certpolysolve::detail::dense(*input.get()));
        require(!piecewise.prepared_structure(), "unprepared piecewise exposed factors");
        cps_options options;
        options.goal = CPS_ROOT_APPROXIMATION;
        options.work_units = 1;
        piecewise.advance(options);
        const auto local = piecewise.prepared_structure();
        require(local.has_value(), "piecewise preparation did not commit");
        check_snapshot(*local, 25, 5, 2);
        for (const auto &factor : local->factors)
            require(factor.power == 4 && factor.shape.nonzero_part_power() == 1 && factor.squarefree_proved,
                    "applied piecewise power substitution was lost");

        // Large apparent degree can remain reduced even when exact algebra is
        // refused. The snapshot must not invent a square-freeness proof.
        certpolysolve::Polynomial binomial;
        binomial.set(0, certpolysolve::Number("-2")).set(120, certpolysolve::Number("1"));
        d::SolveState limited(certpolysolve::detail::dense(*binomial.get()), d::maintained_method_policy(CPS_METHOD_CLASSICAL));
        options.exact_budget = 0;
        options.work_units = 2;
        limited.advance(&options);
        const auto reduced = limited.prepared_structure();
        require(reduced.has_value(), "budgeted raw factor unavailable");
        check_snapshot(*reduced, 120, 0, 1);
        require(reduced->factors.size() == 1 && reduced->factors[0].power == 120 &&
                !reduced->factors[0].squarefree_proved, "budget refusal fabricated proof");
        auto raw_config = d::maintained_piecewise_config();
        raw_config.power_composition = false;
        d::PreparedPiecewiseState raw(certpolysolve::detail::dense(*binomial.get()), raw_config);
        options.work_units = 1;
        raw.advance(options);
        const auto unreduced = raw.prepared_structure();
        require(unreduced.has_value(), "raw piecewise preparation absent");
        check_snapshot(*unreduced, 120, 0, 120);
        require(unreduced->factors[0].squarefree_proved && unreduced->factors[0].power == 1,
                "the exponent-free modular proof was incorrectly tied to expansion storage");
        d::PreparedPiecewiseState refused(certpolysolve::detail::dense(*binomial.get()), raw_config);
        options.exact_operations = 0;
        refused.advance(options);
        const auto unproved = refused.prepared_structure();
        require(unproved.has_value() && !unproved->factors[0].squarefree_proved,
                "refused modular work fabricated a square-freeness proof");
        d::PreparedPiecewiseState elementary(certpolysolve::detail::dense(*binomial.get()));
        elementary.advance(options);
        const auto power_reduced = elementary.prepared_structure();
        require(power_reduced.has_value(), "power preparation requires unnecessary exact expansion");
        check_snapshot(*power_reduced, 120, 0, 1);
        require(power_reduced->factors[0].squarefree_proved && power_reduced->factors[0].power == 120,
                "a nonzero linear factor did not retain its degree-based square-free proof");

        // Cancellation never exposes a half-built preparation. Mutating the
        // caller's input after construction cannot change the owned snapshot.
        d::SolveState cancelled(certpolysolve::detail::dense(*input.get()), d::maintained_method_policy(CPS_METHOD_CLASSICAL));
        d::PreparedPiecewiseState cancelled_local(certpolysolve::detail::dense(*input.get()));
        // The preparation after the cancellation is a fresh state, made before the caller's input changes (one
        // advance a solve state since 18 October 2026).
        d::SolveState prepared(certpolysolve::detail::dense(*input.get()), d::maintained_method_policy(CPS_METHOD_CLASSICAL));
        input.set(25, certpolysolve::Number("7"));
        options = cps_options{};
        options.cancel = [](void *) { return 1; };
        cancelled.advance(&options);
        cancelled_local.advance(options);
        require(!cancelled.prepared_structure() && !cancelled_local.prepared_structure(),
                "cancellation published uncommitted factors");
        options.cancel = nullptr;
        options.work_units = 2;
        prepared.advance(&options);
        check_snapshot(*prepared.prepared_structure(), 25, 5, 2);
        check_snapshot(*local, 25, 5, 2);
        std::cout << "prepared snapshots preserve ownership, degree, multiplicity, power and refusal semantics\n";
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
