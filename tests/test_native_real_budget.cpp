// Native secular real queries share one live-cell allowance across every exact
// factor. Check the public route and its partial answers with independent Sturm
// arithmetic, including multiplicities and a zero root of a separate factor.
#include "exact_reference_api.hpp"
#include <iostream>

using namespace reference;

static certpolysolve::SecularEquation secular(const Polynomial &polynomial) {
    const size_t degree = polynomial.size() - 1;
    std::vector<certpolysolve::SecularEquation::Term> terms;
    Polynomial product{1};
    // These poles lie outside every fixture's roots, so every residue is nonzero.
    constexpr int first_pole = 10;
    for (size_t j = 0; j < degree; ++j)
        product = multiply(product, {-Rational(first_pole + j), 1});
    Polynomial reconstructed = product;
    for (auto &coefficient : reconstructed) coefficient = -coefficient;
    for (size_t j = 0; j < degree; ++j) {
        const Rational pole = first_pole + j;
        const auto [quotient, remainder] = divide(product, {-pole, 1});
        require(remainder.empty(), "reference pole is not a factor of its product");
        const Rational residue = -evaluate(polynomial, pole) / evaluate(quotient, pole);
        require(residue != 0, "test secular residue is zero");
        terms.push_back({.residue = certpolysolve::Number(residue.str()), .pole = certpolysolve::Number(pole.str())});
        for (size_t k = 0; k < quotient.size(); ++k) reconstructed[k] += residue * quotient[k];
    }
    // sum r_j/(x-b_j) - 1 = -p(x)/product_j(x-b_j), exactly.
    for (size_t j = 0; j < polynomial.size(); ++j)
        require(reconstructed[j] == -polynomial[j], "reference secular numerator changed the polynomial");
    certpolysolve::SecularEquation source(terms);
    require(source.ok(), "exact secular test input was refused");
    return source;
}

static void verify(const Polynomial &polynomial, const certpolysolve::Result &result) {
    require(result.get() && result.stop() != certpolysolve::Stop::error,
            "native secular real query returned an API error");
    check_disjoint(result.roots());
    const auto factors = squarefree(polynomial);
    uint64_t total = 0;
    for (const auto &root : result.roots()) {
        const auto lower = endpoint(root.box.real.lo), upper = endpoint(root.box.real.hi);
        require(root.realness == certpolysolve::Realness::real && endpoint(root.box.imag.lo) == 0 &&
                    endpoint(root.box.imag.hi) == 0, "real query did not publish a proved real interval");
        require(roots_closed(polynomial, lower, upper) == static_cast<int>(root.count),
                "independent exact Sturm count disagrees with a returned interval");
        int distinct = 0;
        uint64_t multiplicity = 0;
        for (const auto &[factor, exponent] : factors) {
            const int count = distinct_roots_closed(factor, lower, upper);
            distinct += count;
            if (count) multiplicity = static_cast<uint64_t>(exponent);
        }
        require(distinct == 1 && root.multiplicity == multiplicity && root.count == multiplicity,
                "a real interval lost its independently established multiplicity");
        if (result.goal_met()) check_accuracy(root, 53);
        total += root.count;
    }
    if (result.complete())
        require(total == polynomial.size() - 1, "complete real query lost roots of an all-real polynomial");
}

static void budgets(const Polynomial &polynomial) {
    const auto source = secular(polynomial);
    for (const bool accuracy : {false, true}) {
        for (const uint64_t cells : {1, 2, 3, 4, 5, 6, 8, 16, 65536}) {
            for (const uint64_t work : {1, 3, 8, 4096}) {
                certpolysolve::Options options;
                options.bits = 53;
                options.max_bits = 256;
                options.goal = accuracy ? certpolysolve::Goal::accuracy() : certpolysolve::Goal::approximation();
                options.cell_budget = cells;
                options.work_budget = work;
                const auto result = certpolysolve::real_roots(source, options);
                verify(polynomial, result);
                const auto counters = result.engine_counters();
                if (counters.live_cells > cells)
                    throw std::runtime_error("native real query exceeded aggregate cell allowance: " +
                        std::to_string(counters.live_cells) + " live, budget " + std::to_string(cells) +
                        ", work " + std::to_string(work));
                if (work == 4096) {
                    require(counters.real_refinement_steps > 0, "fixture missed the real refinement route");
                    if (cells >= 16)
                        require(result.complete() && result.goal_met(), "ample cell allowance did not finish");
                    if (!result.goal_met())
                        require(cps_result_stop_reason(result.get()) == CPS_STOP_CELL_BUDGET,
                                "cell-limited real query reported an unrelated stop");
                }
            }
        }
    }
}

int main() {
    try {
        const Polynomial two_factors = multiply({-2, 0, 1}, multiply({-3, 0, 1}, {-3, 0, 1}));
        budgets(two_factors);
        budgets(multiply(two_factors, {0, 0, 0, 1}));
        budgets(multiply(two_factors, multiply({-5, 0, 1}, multiply({-5, 0, 1}, {-5, 0, 1}))));
        std::cout << "Native real cell budgets and independent interval counts passed\n";
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
