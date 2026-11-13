#include "solve/policies.hpp"
#include "numeric/arithmetic/polynomial_arithmetic.hpp"
#include "internal.hpp"
#include "exact_winding_api.hpp"
#include "api_test_support.hpp"
#include "numeric/candidates/candidate_jets.hpp"
#include "numeric/candidates/machine_adapter.hpp"
#include "numeric/candidates/study_candidates.hpp"
#include "numeric/evaluate/structured_source.hpp"
#include "methods/secular/secular_refinement.hpp"
#include <iostream>
#include <certpolysolve/certpolysolve.hpp>

using namespace certpolysolve::detail;
namespace {
void require(bool value, const char *message) {
    if (!value) throw std::runtime_error(message);
}
// The D5 policy with its private counters attached (observation only).
EnginePolicy observed() {
    auto p = production_policy();
    p.study_statistics = std::make_shared<StudyStatistics>();
    return p;
}
void check(const cps_result &raw, const std::vector<reference::Complex> &known, uint64_t bits) {
    const auto result = certpolysolve_test::published(raw);
    require(result.goal_met() && result.complete(), "D3 did not complete the known-root input");
    const auto &roots = result.roots();
    reference::check_disjoint(roots);
    uint64_t total = 0;
    for (const auto &r : roots) {
        uint64_t count = 0;
        for (const auto &z : known)
            count += reference::endpoint(r.box.real.lo) <= z.real && z.real <= reference::endpoint(r.box.real.hi) &&
                     reference::endpoint(r.box.imag.lo) <= z.imag && z.imag <= reference::endpoint(r.box.imag.hi);
        require(count == r.count, "D3 enclosure disagrees with exact roots");
        total += count;
        reference::check_accuracy(r, bits);
    }
    require(total == known.size(), "D3 lost roots or multiplicities");
}
void uncertainty() {
    certpolysolve::Polynomial input;
    input.set(0, certpolysolve::ball(certpolysolve::Number("-1"), certpolysolve::Number("1/16")));
    input.set(2, certpolysolve::Number("1"));
    SolveState state(certpolysolve::detail::dense(*input.get()), production_policy());
    cps_options o;
    o.bits = 53;
    o.max_work_bits = 128;
    o.work_units = 64;
    require(!state.advance(&o).goal_met, "arithmetic continuation erased coefficient uncertainty");
}
void representations_and_tighter_precision() {
    const std::vector<reference::Complex> known{{1}, {2}, {4}};
    const auto expanded = reference::input(reference::from_roots(known));
    const certpolysolve::SecularEquation native{{.residue = certpolysolve::Number("36/5"), .pole = certpolysolve::Number("-2")},
                                           {.residue = certpolysolve::Number("-4/3"), .pole = certpolysolve::Number("0")},
                                           {.residue = certpolysolve::Number("2/15"), .pole = certpolysolve::Number("3")}};
    const auto policy = production_policy();
    // Each representation at each goal: a cancelled solve says so, a one-unit solve stops honestly, and a solve
    // with the 128 units the one-unit calls of a session used to share meets the goal (one advance a solve state
    // since 18 October 2026; until then one state was cancelled, then resumed a unit at a time).
    for (const bool native_input : {false, true}) {
        const auto fresh = [&] {
            // Two returns rather than ?: because MSVC copies the operands of ?:, and SolveState can be neither
            // copied nor moved (C2280); each return of a prvalue is elided on every compiler.
            if (native_input) return SolveState(*native.get(), policy);
            return SolveState(certpolysolve::detail::dense(*expanded.get()), policy);
        };
        for (uint64_t bits : {53, 128, 256}) {
            cps_options options;
            options.max_work_bits = 2048;
            options.bits = bits;
            options.work_units = 1;
            options.cancel = [](void *) { return 1; };
            auto cancelled = fresh();
            require(cancelled.advance(&options).termination == CPS_CANCELLED, "a representation ignored cancellation");
            options.cancel = nullptr;
            auto unit = fresh();
            const auto stopped = unit.advance(&options);
            require(stopped.goal_met || stopped.termination == CPS_WORK_LIMIT,
                    "a one-unit solve of a representation stopped for another reason");
            options.work_units = 128;
            auto whole = fresh();
            check(whole.advance(&options), known, bits);
        }
    }
}
// The double-double machine arithmetic (machine_arithmetic.hpp), which the
// research arithmetic tiers still use; the D3 bridge that drove it is gone.
void double_double_ratio() {
    // An independent ball calculation checks the ~106-bit reciprocal ratio,
    // including coefficients with exact zero factors and a very large forward
    // polynomial value. A binary64-rounded result fails this tolerance.
    using R = machine::DoubleDouble;
    using C = machine::Complex<R>;
    std::vector<C> coefficients(2001, C(R(1)));
    coefficients[0] = {};
    machine::PreparedNewton<R> evaluator(coefficients);
    const C point{R(1.25), R(.25)};
    const auto ratio = evaluator.evaluate(point).ratio();
    require(ratio.has_value(), "double-double reciprocal evaluation failed");
    ComplexPolynomial polynomial;
    for (slong j = 1; j <= 2000; ++j) acb_poly_set_coeff_si(polynomial.value(), j, 1);
    ComplexBall z, value, derivative, expected, actual, error;
    to_ball(z.value(), point);
    acb_poly_evaluate2_rectangular(value.value(), derivative.value(), polynomial.value(), z.value(), 2048);
    acb_div(expected.value(), value.value(), derivative.value(), 2048);
    to_ball(actual.value(), *ratio);
    acb_sub(error.value(), actual.value(), expected.value(), 2048);
    Magnitude bound;
    acb_get_mag(bound.value(), error.value());
    require(mag_cmp_2exp_si(bound.value(), -90) < 0, "double-double Newton ratio lost significand accuracy");
}
// A residue below 2^-2b of its node is a root of the binary64 model, not a
// range failure: a certified node 2^40 away with a residue of 2^-1000 made
// every term it contributes subnormal, the packet raised underflow and was
// refused for the multiprecision engine.
void secular_residue_below_resolution() {
    std::vector<ComplexBall> nodes(4), residues(4);
    acb_set_si(nodes[0].value(), -1);
    acb_set_si(nodes[1].value(), 1);
    acb_onei(nodes[2].value());
    acb_one(nodes[3].value());
    acb_mul_2exp_si(nodes[3].value(), nodes[3].value(), 40);
    for (size_t i = 0; i < 3; ++i) {
        acb_one(residues[i].value());
        acb_mul_2exp_si(residues[i].value(), residues[i].value(), -4);
    }
    acb_one(residues[3].value());
    acb_mul_2exp_si(residues[3].value(), residues[3].value(), -1000);
    RootVector initial(4), out(4);
    for (slong i = 0; i < 4; ++i) acb_set(initial.data() + i, nodes[static_cast<size_t>(i)].value());
    for (slong i = 0; i < 3; ++i) {
        ComplexBall step;
        acb_set_d_d(step.value(), 0.125, 0.0625);
        acb_add(initial.data() + i, initial.data() + i, step.value(), 64);
    }
    const std::vector<uint8_t> inactive{0, 0, 0, 1};
    const std::vector<slong> bits(4, 85);
    cps_options o;
    o.bits = 53;
    o.max_work_bits = 256;
    certpolysolve::detail::CandidateStatistics statistics{};
    auto p = observed();
    SecularMachineWorkspace workspace;
    CandidateProgress progress;
    require(workspace.advance(progress, out.data(), initial.data(), nodes, residues, inactive, bits, 0, 1, o,
                              statistics, p) && progress.finite,
            "a residue below its node's resolution refused the binary64 packet");
    require(p.study_statistics->continuations == 0, "the binary64 packet handed off a resolved model");
}
void contains_exact(const acb_t ball, const reference::Complex &value) {
    Rational real, imag;
    require(fmpq_set_str(real.value(), value.real.str().c_str(), 10) == 0 &&
            fmpq_set_str(imag.value(), value.imag.str().c_str(), 10) == 0, "exact reference conversion failed");
    require(arb_contains_fmpq(acb_realref(ball), real.value()) &&
            arb_contains_fmpq(acb_imagref(ball), imag.value()), "ball evaluation omitted an exact reference value");
}
void rectangular_enclosures() {
    // The rectangular probe (rectangular_available) leaves a sparse source to its sparse evaluation. Until W8's S2
    // this also checked PolynomialSource::evaluate_rectangular, a leftover of the D3 retry with no caller.
    PolynomialSource sparse([](acb_poly_t out, slong) {
        acb_poly_set_coeff_si(out, 0, 1);
        acb_poly_set_coeff_si(out, 31, 1);
    }, 31, true);
    require(!sparse.rectangular_available(128), "rectangular probe replaced sparse evaluation");
}
// Residue regeneration (the default policy) against exact residues, and a
// cancelled regeneration keeping the published ones. Until 7 October 2026
// this also screened the D3 rectangular retry, since removed.
void residue_regeneration() {
    constexpr size_t degree = 128;
    PolynomialSource source([](acb_poly_t out, slong) {
        for (slong j = 0; j <= static_cast<slong>(degree); ++j) acb_poly_set_coeff_si(out, j, 1);
    }, degree, true);
    std::vector<ComplexBall> nodes(degree), residues;
    std::vector<slong> bits(degree, 32);
    for (size_t j = 0; j < degree; ++j) {
        arb_set_ui(acb_realref(nodes[j].value()), 160 + static_cast<ulong>(j));
        arb_mul_2exp_si(acb_realref(nodes[j].value()), acb_realref(nodes[j].value()), -7);
        arb_set_d(acb_imagref(nodes[j].value()), .25);
    }
    cps_options o;
    o.bits = 24;
    o.max_work_bits = 256;
    certpolysolve::detail::CandidateStatistics statistics{};
    uint64_t maximum = 0;
    require(regenerate_residues(residues, source, nodes, {}, {}, bits, o, statistics, maximum),
            "reference residue regeneration failed");
    const reference::ComplexPolynomial polynomial(degree + 1, reference::Complex(1));
    for (const size_t i : {size_t(0), degree/2, degree-1}) {
        const reference::Complex node{reference::Rational(160 + i, 128), reference::Rational(1, 4)};
        reference::Rational product(1);
        for (size_t j = 0; j < degree; ++j)
            if (j != i) product *= reference::Rational(static_cast<long>(i) - static_cast<long>(j), 128);
        const auto value = reference::evaluate(polynomial, node);
        contains_exact(residues[i].value(), {-value.real/product, -value.imag/product});
    }
    // A cancelled new evaluation must retain the previously published residues.
    const auto saved = residues;
    o.cancel = [](void *) { return 1; };
    require(!regenerate_residues(residues, source, nodes, {}, {}, bits, o, statistics, maximum),
            "regeneration ignored cancellation");
    for (size_t i = 0; i < degree; ++i)
        require(acb_equal(residues[i].value(), saved[i].value()), "cancelled regeneration changed retained enclosures");
}
}
int main() {
    certpolysolve::ThreadScope scope;
    try {
        uncertainty();
        representations_and_tighter_precision();
        double_double_ratio();
        secular_residue_below_resolution();
        rectangular_enclosures();
        residue_regeneration();
        std::cout << "D3 independent checks passed\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
