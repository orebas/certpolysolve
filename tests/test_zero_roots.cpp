#include "solve/policies.hpp"
// S8's B1 and B2 in the global driver (SolveState::Impl::preprocess, solver.cpp; the S8 map's rules 4 and 5).
// B1 sets the exact zero roots apart: the point 0, counted v times, joins every
// assembly beside the factors of the nonzero part q = p / x^v. Checked under forced classical, forced secular and
// Auto at the accuracy and approximation goals, on x (x^2 - 2), which FLINT kept whole as x^3 - 2x; the power
// composite x^2 (x^4 - 3); x^3 (x - 5); 7 x^5, with zero roots only; and x^2 (x^2 + 1)^2 (x - 3). The goal is met,
// the regions are disjoint, each region's count is the input's own count in it (rectangle_count, exact) and the
// counts sum to the degree; exactly one region is the point 0, with count and multiplicity v, real, certified as an
// exact zero monomial factor; and the square-free verdict is disproved where v is 2 or more. At work budgets of 2, 3
// and 4 units (C3's stop covers) under both forced methods the answer is complete and exact in its counts: a stop
// cover that meets the point 0 merges with it into one counted region, as overlapping factor covers merge (rule 1). B2 proves q square-free by the modular proof before FLINT: on a square-free input the preparation's trace
// reports the modular proof, on an input with a repeated factor FLINT's factorisation; and the factor the proof stands
// for is FLINT's own (fmpz_poly_factor_squarefree on 400 random square-free integer polynomials of degree 2 to 40:
// one factor, the primitive part with a positive leading coefficient, multiplicity 1).
#include "numeric/arithmetic/polynomial_arithmetic.hpp"
#include "exact_winding_api.hpp"
#include "api_test_support.hpp"
#include "numeric/exact/exact_preprocessing.hpp"
#include "internal.hpp"
#include <flint/fmpz_poly.h>
#include <iostream>
#include <random>
#include <string>

using namespace reference;
namespace detail = certpolysolve::detail;

namespace {
const std::string zero_certificate = "exact zero monomial factor before square-free preparation";

bool zero_point(const certpolysolve::Root &r) {
    return endpoint(r.box.real.lo) == Rational(0) && endpoint(r.box.real.hi) == Rational(0) &&
           endpoint(r.box.imag.lo) == Rational(0) && endpoint(r.box.imag.hi) == Rational(0);
}

ComplexPolynomial integers(std::initializer_list<int> coefficients) {
    ComplexPolynomial p;
    for (const int c : coefficients) p.push_back(Complex(Rational(c), Rational(0)));
    return p;
}

struct Case {
    std::string name;
    ComplexPolynomial polynomial;
    uint64_t degree, zeros;
};

// goal: the solve must meet its goal, and then the point 0 is a region of its own.
void check(const Case &c, const certpolysolve::Result &result, const std::string &what, bool goal) {
    require(result.complete(), (what + ": not complete: " + result.message()).c_str());
    if (goal) require(result.goal_met(), (what + ": the goal is not met: " + result.message()).c_str());
    const auto regions = result.roots();
    check_disjoint(regions);
    uint64_t total = 0, points = 0;
    for (const auto &region : regions) {
        const auto count = rectangle_count(c.polynomial, region);
        require(count && static_cast<uint64_t>(*count) == region.count,
                (what + ": a region's count is not the input's").c_str());
        total += region.count;
        if (zero_point(region)) {
            ++points;
            require(region.count == c.zeros && region.multiplicity == c.zeros &&
                        region.realness == certpolysolve::Realness::real && region.certificate == zero_certificate,
                    (what + ": the point 0 is not the zero roots set apart: certificate " + region.certificate).c_str());
        }
    }
    require(total == c.degree, (what + ": the counts do not sum to the degree").c_str());
    if (goal) require(points == 1, (what + ": no region is the point 0").c_str());
    if (c.zeros >= 2)
        require(result.squarefree() == certpolysolve::Squarefree::disproved, (what + ": a double zero root not reported").c_str());
}

void zero_roots_apart() {
    const std::vector<Case> cases{
        {"x (x^2 - 2)", integers({0, -2, 0, 1}), 3, 1},
        {"x^2 (x^4 - 3)", integers({0, 0, -3, 0, 0, 0, 1}), 6, 2},
        {"x^3 (x - 5)", integers({0, 0, 0, -5, 1}), 4, 3},
        {"7 x^5", integers({0, 0, 0, 0, 0, 7}), 5, 5},
        // x^2 (x^2 + 1)^2 (x - 3) = x^7 - 3 x^6 + 2 x^5 - 6 x^4 + x^3 - 3 x^2
        {"x^2 (x^2 + 1)^2 (x - 3)", integers({0, 0, -3, 1, -6, 2, -3, 1}), 7, 2},
    };
    for (const auto &c : cases)
        for (const certpolysolve::Method method : {certpolysolve::Method::classical, certpolysolve::Method::secular, certpolysolve::Method::automatic})
            for (const auto &goal : {certpolysolve::Goal::accuracy(), certpolysolve::Goal::approximation()}) {
                certpolysolve::Options options;
                options.method = method;
                options.goal = goal;
                const std::string what = c.name + ", method " + std::to_string(static_cast<int>(method)) + ", goal " +
                                         std::to_string(static_cast<int>(goal.value()));
                check(c, certpolysolve::solve(input(c.polynomial), options), what, true);
            }
    // At a stop short of the goal (C3): every factor without a step completes with its stop cover, and the point 0
    // is published beside the factors.
    for (const auto &c : cases)
        for (const certpolysolve::Method method : {certpolysolve::Method::classical, certpolysolve::Method::secular})
            for (const uint64_t budget : {uint64_t{2}, uint64_t{3}, uint64_t{4}}) {
                certpolysolve::Options options;
                options.method = method;
                options.work_budget = budget;
                const std::string what = c.name + ", method " + std::to_string(static_cast<int>(method)) + ", budget " +
                                         std::to_string(budget);
                check(c, certpolysolve::solve(input(c.polynomial), options), what, false);
            }
}

// The preparation stage the trace reports for the global driver's factorisation ("exact-squarefree"): its reason.
std::vector<std::string> preparation_reasons;
std::vector<std::string> traced_preparation(const ComplexPolynomial &p) {
    auto policy = detail::production_policy();
    policy.trace = [](void *, const detail::TraceEvent &event) {
        if (std::string(event.stage) == "exact-squarefree") preparation_reasons.emplace_back(event.reason);
    };
    preparation_reasons.clear();
    auto polynomial = input(p);
    detail::SolveState state(detail::dense(*polynomial.get()), policy);
    cps_options options;
    options.method = CPS_METHOD_CLASSICAL;
    state.advance(&options);
    return preparation_reasons;
}

void modular_proof_first() {
    const auto square_free = traced_preparation(integers({-6, 11, -6, 1}));          // (x - 1)(x - 2)(x - 3)
    const auto repeated = traced_preparation(integers({-12, 16, -7, 1}));            // (x - 2)^2 (x - 3)
    const auto zeros_apart = traced_preparation(integers({0, 0, 0, -6, 11, -6, 1})); // x^3 (x - 1)(x - 2)(x - 3)
    require(square_free == std::vector<std::string>{"modular-squarefree-proof"},
            "a square-free input was factored by FLINT");
    require(repeated == std::vector<std::string>{"exact-squarefree-factors"},
            "an input with a repeated factor was not factored by FLINT");
    require(zeros_apart == std::vector<std::string>{"modular-squarefree-proof"},
            "the nonzero part of x^3 q, q square-free, was factored by FLINT");
    // The one factor the proof stands for is FLINT's own on square-free input.
    std::mt19937_64 random(20261023);
    int checked = 0;
    while (checked < 400) {
        detail::IntegerPolynomial p;
        const slong degree = 2 + static_cast<slong>(random() % 39);
        for (slong i = 0; i <= degree; ++i)
            fmpz_poly_set_coeff_si(p.value(), i, static_cast<slong>(random() % 2001) - 1000);
        if (fmpz_poly_degree(p.value()) < 2 || !fmpz_poly_is_squarefree(p.value())) continue;
        const auto factors = detail::squarefree_factors(p);
        detail::IntegerPolynomial primitive;
        fmpz_poly_primitive_part(primitive.value(), p.value());
        require(factors.size() == 1 && factors.front().multiplicity == 1 &&
                    fmpz_poly_equal(factors.front().polynomial.value(), primitive.value()),
                "FLINT's square-free factor of a square-free input is not its primitive part");
        ++checked;
    }
}
} // namespace

int main() {
    certpolysolve::ThreadScope lifetime;
    try {
        zero_roots_apart();
        modular_proof_first();
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
    std::cout << "zero roots tests passed\n";
    return 0;
}
