// Check an actual Aberth update against an enclosed direct formula on both
// sides of the interaction-tree crossovers. The reference differentiates
// p(z) = z^n - 1 directly and sums reciprocal distances at higher precision;
// it uses neither the candidate evaluator nor an interaction tree.
#include "internal.hpp"
#include "core/heuristics.hpp"
#include "numeric/candidates/candidates.hpp"
#include "numeric/candidates/machine_candidates.hpp"
#include <certpolysolve/certpolysolve.hpp>
#include <cmath>
#include <complex>
#include <iostream>
#include <numbers>
#include <stdexcept>
#include <vector>

namespace d = certpolysolve::detail;
namespace {
constexpr slong reference_bits = 256;
// The existing machine_tree regression asks for 40 agreeing bits between
// candidate trajectories. Use the same accuracy here against a separate
// enclosed formula, not as a claimed forward-error bound for arbitrary data.
constexpr slong checked_bits = 40;

void require(bool condition, const char *message) {
    if (!condition) throw std::runtime_error(message);
}

void initialize(d::RootVector &points, slong degree, slong rotation) {
    // A visible displacement from the unit roots keeps every point moving
    // during this single sweep, including at the machine-tree crossover.
    constexpr double radius = 1 + 0x1p-10;
    for (slong j = 0; j < degree; ++j) {
        const double angle = 2 * std::numbers::pi * static_cast<double>((j + rotation) % degree) /
                             static_cast<double>(degree);
        const auto point = std::polar(radius, angle);
        acb_set_d_d(points[j], point.real(), point.imag());
    }
}

void check_update(acb_srcptr initial, acb_srcptr updated, slong degree, slong target, bool gauss_seidel) {
    d::ComplexBall value, derivative, sum, difference, term, denominator, expected;
    const auto z = initial + target;
    acb_pow_ui(value.value(), z, static_cast<ulong>(degree), reference_bits);
    acb_sub_ui(value.value(), value.value(), 1, reference_bits);
    acb_pow_ui(derivative.value(), z, static_cast<ulong>(degree - 1), reference_bits);
    acb_mul_ui(derivative.value(), derivative.value(), static_cast<ulong>(degree), reference_bits);
    for (slong j = 0; j < degree; ++j) {
        if (j == target) continue;
        // Binary64 follows earlier updates within the sweep. Using the
        // exported exact binary64 neighbours isolates the update at target
        // from errors already made at earlier targets.
        const auto neighbour = gauss_seidel && j < target ? updated + j : initial + j;
        acb_sub(difference.value(), z, neighbour, reference_bits);
        acb_inv(term.value(), difference.value(), reference_bits);
        acb_add(sum.value(), sum.value(), term.value(), reference_bits);
    }
    // z - p(z) / (p'(z) - p(z) sum_{j != target} 1/(z-z_j)).
    acb_mul(denominator.value(), value.value(), sum.value(), reference_bits);
    acb_sub(denominator.value(), derivative.value(), denominator.value(), reference_bits);
    require(!acb_contains_zero(denominator.value()), "reference Aberth denominator is unresolved");
    acb_div(term.value(), value.value(), denominator.value(), reference_bits);
    acb_sub(expected.value(), z, term.value(), reference_bits);
    require(acb_is_finite(expected.value()) && acb_rel_accuracy_bits(expected.value()) > 3 * checked_bits,
            "reference enclosure is too wide for the comparison");
    require(!acb_overlaps(expected.value(), z), "fixture did not propose a genuine update");
    arb_add_error_2exp_si(acb_realref(expected.value()), -checked_bits);
    arb_add_error_2exp_si(acb_imagref(expected.value()), -checked_bits);
    require(acb_contains(expected.value(), updated + target),
            "candidate update disagrees with the enclosed direct Aberth formula");
}

void multiprecision_boundary(slong degree) {
    d::ComplexPolynomial polynomial;
    acb_poly_set_coeff_si(polynomial.value(), 0, -1);
    acb_poly_set_coeff_si(polynomial.value(), degree, 1);
    for (const slong rotation : {slong{0}, degree / 3}) {
        d::RootVector initial(degree), output(degree);
        initialize(initial, degree, rotation);
        // Only the first point moves: a classical packet's far fields and
        // near fields then refer to the same positions, despite Gauss-Seidel.
        std::vector<uint8_t> inactive(static_cast<size_t>(degree), 1);
        inactive[0] = 0;
        cps_options options;
        d::CandidateStatistics statistics{};
        const auto progress = d::candidate_packet(output.data(), polynomial.value(), initial.data(), 128, 1,
                                                  d::CandidateEngine::Aberth, options, statistics, inactive, true);
        require(progress.finite && !progress.cancelled && statistics.candidate_sweeps == 1,
                "multiprecision boundary packet was not accepted");
        check_update(initial.data(), output.data(), degree, 0, false);
        for (slong j = 1; j < degree; ++j)
            require(acb_equal(output[j], initial[j]), "packet changed an inactive proposal");
        // Branch witnesses supplement the mathematical check; they keep a
        // silent loss of tree coverage from turning the regression vacuous.
        const auto direct_terms = static_cast<uint64_t>(degree - 1);
        require(degree < d::heuristics::interaction_tree_degree
                    ? statistics.root_interactions == direct_terms
                    : statistics.root_interactions < direct_terms,
                "fixture did not exercise the intended direct/tree interaction path");
        std::cout << "multiprecision degree " << degree << " rotation " << rotation
                  << ": " << statistics.root_interactions << " direct interaction terms\n";
    }
}

void machine_boundary(slong degree) {
    d::ComplexPolynomial polynomial;
    acb_poly_set_coeff_si(polynomial.value(), 0, -1);
    acb_poly_set_coeff_si(polynomial.value(), degree, 1);
    d::RootVector initial(degree), output(degree);
    initialize(initial, degree, 0);
    auto engine = d::prepare_machine_engine(d::IterationMethod::Aberth, d::ArithmeticTier::Binary64,
                                            polynomial.value(), initial.view(), false, false, false,
                                            UINT64_MAX, false, true);
    require(engine != nullptr, "binary64 boundary engine was refused");
    cps_options options;
    d::CandidateStatistics statistics{};
    const auto progress = engine->advance(1, options, statistics);
    require(!progress.arithmetic_failure && !progress.candidates.cancelled && statistics.candidate_sweeps == 1,
            "binary64 boundary packet was not accepted");
    engine->export_candidates(output.data());
    for (const slong target : {slong{0}, degree / 4, degree / 2, degree - 1})
        check_update(initial.data(), output.data(), degree, target, true);
    const auto direct_terms = static_cast<uint64_t>(degree) * static_cast<uint64_t>(degree - 1);
    require(static_cast<size_t>(degree) < d::heuristics::machine_interaction_tree_degree
                ? statistics.root_interactions == direct_terms
                : statistics.root_interactions < direct_terms,
            "fixture did not exercise the intended machine direct/tree path");
    std::cout << "binary64 degree " << degree << ": " << statistics.root_interactions
              << " charged interaction terms\n";
}
} // namespace

int main() {
    certpolysolve::ThreadScope lifetime;
    try {
        for (const slong delta : {-1, 0, 1}) {
            multiprecision_boundary(d::heuristics::interaction_tree_degree + delta);
            machine_boundary(static_cast<slong>(d::heuristics::machine_interaction_tree_degree) + delta);
        }
        std::cout << "interaction crossover candidate updates match the direct mathematical reference\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
