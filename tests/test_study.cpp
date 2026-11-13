#include "solve/policies.hpp"
#include "numeric/arithmetic/polynomial_arithmetic.hpp"
#include "internal.hpp"
#include "api_test_support.hpp"
#include "numeric/candidates/candidate_jets.hpp"
#include "numeric/candidates/candidates.hpp"
#include "core/heuristics.hpp"
#include "exact_winding_api.hpp"
#include "numeric/candidates/study_candidates.hpp"
#include <cfenv>
#include <cmath>
#include <iostream>

using namespace certpolysolve::detail;
void
require(bool value, const char *message) {
    if (!value) throw std::runtime_error(message);
}
EnginePolicy
policy(StudyMethod method, bool machine) {
    EnginePolicy p;
    p.study = method;
    p.study_machine = machine;
    p.refinement = method == StudyMethod::None ? RefinementMethod::Automatic :
        method == StudyMethod::Aberth ? RefinementMethod::Aberth : RefinementMethod::Secular;
    p.study_statistics = std::make_shared<StudyStatistics>();
    return p;
}
void
check(const cps_result &raw, const std::vector<reference::Complex> &known, uint64_t bits, bool require_goal = true) {
    const auto result = certpolysolve_test::published(raw);
    require(result.complete(), "study did not retain a complete root cover");
    const auto &roots = result.roots();
    reference::check_disjoint(roots);
    uint64_t total = 0;
    for (const auto &r : roots) {
        uint64_t count = 0;
        for (const auto &z : known)
            count += reference::endpoint(r.box.real.lo) <= z.real && z.real <= reference::endpoint(r.box.real.hi) &&
                     reference::endpoint(r.box.imag.lo) <= z.imag && z.imag <= reference::endpoint(r.box.imag.hi);
        require(count == r.count, "study certificate disagrees with independent exact roots");
        total += count;
        if (result.goal_met()) reference::check_accuracy(r, bits);
    }
    require(total == known.size(), "study lost root coverage");
    if (require_goal) require(result.goal_met(), "study did not complete known-root input");
}
void
known_roots() {
    std::vector<std::vector<reference::Complex>> examples{
        { { -3 }, { -1 }, { 2 }, { 4 }, { 7 }, { 11 }, { 13 }, { 17 } },
        { { -3 }, { 0 }, { 0 }, { 0 }, { 1 }, { 1 }, { 2 } },
        { { -5 },
          { -2 },
          { 1 },
          { reference::Rational(1) + reference::Rational(1) / (reference::Integer(1) << 40) },
          { 7 } },
        { { 1, 2 }, { -2, 1 }, { 3, -1 }, { 5, 4 }, { -4, -3 } }
    };
    std::vector<reference::Complex> translated;
    for (int j = 0; j < 64; ++j) translated.push_back({ reference::Rational(1048576) + reference::Rational(j) / 16 });
    examples.push_back(translated);
    for (auto method : { StudyMethod::None, StudyMethod::Aberth, StudyMethod::Secular })
        for (bool machine : { false, true }) {
            for (const auto &known : examples) {
                auto input = reference::input(reference::from_roots(known));
                auto p = policy(method, machine);
                cps_options o;
                o.bits = 32;
                o.max_work_bits = 4096;
                o.work_units = 2048;
                o.exact_budget = 8000000;
                // This clustered real-root construction checks translated covers
                // even when a fixed method hits its limit. Completion is measured
                // in the campaign, not asserted for every method in a unit test.
                if (known.size() == 64) o.max_work_bits = 256;
                SolveState state(certpolysolve::detail::dense(*input.get()), p);
                auto result = state.advance(&o);
                check(result, known, o.bits, known.size() != 64);
                if (known.size() == 8) {
                    if (method == StudyMethod::None && machine)
                        require(p.study_statistics->machine_packets > 0,
                                "C2 integration did not enter ordinary machine arithmetic");
                    // A tighter goal is a fresh solve (one advance a solve state since 18 October 2026).
                    o.bits = 96;
                    SolveState tighter(certpolysolve::detail::dense(*input.get()), p);
                    result = tighter.advance(&o);
                    check(result, known, o.bits);
                }
            }
        }
}
void
secular_packet() {
    ComplexPolynomial p;
    acb_poly_set_coeff_si(p.value(), 0, -8);
    acb_poly_set_coeff_si(p.value(), 1, 14);
    acb_poly_set_coeff_si(p.value(), 2, -7);
    acb_poly_set_coeff_si(p.value(), 3, 1);
    std::vector<ComplexBall> nodes(3), residues;
    acb_set_si(nodes[0].value(), -2);
    acb_zero(nodes[1].value());
    acb_set_si(nodes[2].value(), 3);
    RootVector initial(3), ordinary(3), adaptive(3);
    for (slong j = 0; j < 3; ++j) acb_set(initial.data() + j, nodes[static_cast<size_t>(j)].value());
    cps_options o;
    o.bits = 53;
    o.max_work_bits = 512;
    certpolysolve::detail::CandidateStatistics s{};
    require(regenerate_secular(residues, p.value(), nodes, 128, o, s), "test residue regeneration");
    auto control = policy(StudyMethod::Secular, true);
    SecularMachineWorkspace cache;
    const std::vector<uint8_t> inactive{ 0, 1, 0 };
    const std::vector<slong> bits(3, 85);
    for (uint64_t revision : { 1, 1, 2 }) {
        CandidateProgress progress;
        require(cache.advance(
                  progress, ordinary.data(), initial.data(), nodes, residues, inactive, bits, revision, 1, o, s, control),
                "machine secular packet rejected");
        const auto exact =
          secular_candidate_packet(adaptive.data(), initial.data(), nodes, residues, inactive, bits, 1, o, s);
        require(exact.finite && progress.finite, "nonfinite secular proposal");
        for (slong j = 0; j < 3; ++j) {
            ComplexBall difference;
            acb_sub(difference.value(), ordinary.data() + j, adaptive.data() + j, 128);
            Magnitude error;
            acb_get_mag(error.value(), difference.value());
            require(mag_get_d(error.value()) < 1e-11, "machine/adaptive secular correction disagreement");
        }
        require(acb_equal(ordinary.data() + 1, initial.data() + 1), "machine moved certified inactive root");
    }
    require(control.study_statistics->preparations == 2 && control.study_statistics->reuses == 1,
            "secular revision cache contract");
    CandidateProgress progress;
    o.max_work_bits = 32;
    require(!cache.advance(progress, ordinary.data(), initial.data(), nodes, residues, inactive, bits, 3, 1, o, s, control),
            "machine exceeded precision ceiling");
    o.max_work_bits = 512;
    std::fesetround(FE_DOWNWARD);
    require(!cache.advance(progress, ordinary.data(), initial.data(), nodes, residues, inactive, bits, 3, 1, o, s, control),
            "machine admitted unsupported rounding mode");
    require(std::fegetround() == FE_DOWNWARD, "machine changed caller rounding mode");
    std::fesetround(FE_TONEAREST);
}
// The binary64 packet's sums in lanes (EnginePolicy::packet_lane_sums; PacketSums in study_candidates.cpp) against
// the one-accumulator loop: the same model, the same points, a few sweeps each way. A packet is rejected on any
// floating-point exception, so "accepted" is the test that the lanes raise none: points on their own nodes (where
// a packet starts, the model's poles), a point on another node, points beside the nodes. The proposals agree
// to rounding. Parts outside the lanes' range leave the packet to the general loop, bit for bit; two equal
// points are refused either way.
void
lane_sums() {
    constexpr slong n = 23;
    auto unit = [state = uint64_t{88172645463325252ULL}]() mutable {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        return static_cast<double>(state >> 11) * 0x1p-53;
    };
    RootVector roots(n);
    std::vector<ComplexBall> nodes(n), residues;
    for (slong j = 0; j < n; ++j) {
        const double re = 2 * unit() - 1, im = 2 * unit() - 1;
        acb_set_d_d(roots[j], re, im);
        // Nodes a little off the roots, so that the residues are small and the packet has steps to take.
        acb_set_d_d(nodes[static_cast<size_t>(j)].value(), re + 0x1p-12 * (2 * unit() - 1), im + 0x1p-12 * (2 * unit() - 1));
    }
    ComplexPolynomial p;
    acb_poly_product_roots(p.value(), roots.data(), n, 2048);
    cps_options o;
    o.bits = 53;
    o.max_work_bits = 512;
    certpolysolve::detail::CandidateStatistics s{};
    require(regenerate_secular(residues, p.value(), nodes, 256, o, s), "lane test residue regeneration");
    const std::vector<uint8_t> inactive(static_cast<size_t>(n), 0);
    const std::vector<slong> bits(static_cast<size_t>(n), 85);
    const auto run = [&](RootVector &out, const RootVector &initial, const std::vector<ComplexBall> &at,
                         const std::vector<ComplexBall> &weights, bool lanes) {
        auto control = policy(StudyMethod::Secular, true);
        control.packet_lane_sums = lanes;
        SecularMachineWorkspace cache;
        CandidateProgress progress;
        certpolysolve::detail::CandidateStatistics statistics{};
        return cache.advance(progress, out.data(), initial.data(), at, weights, inactive, bits, 1, 4, o, statistics,
                             control) && progress.finite;
    };
    const auto agree = [&](const RootVector &a, const RootVector &b, double tolerance) {
        for (slong j = 0; j < n; ++j) {
            ComplexBall difference;
            Magnitude error, size;
            acb_sub(difference.value(), a[j], b[j], 128);
            acb_get_mag(error.value(), difference.value());
            acb_get_mag(size.value(), a[j]);
            if (mag_get_d(error.value()) > tolerance * (1 + mag_get_d(size.value()))) return false;
        }
        return true;
    };
    RootVector initial(n), plain(n), lanes(n);
    // On their own nodes: every point of the first sweep is at a pole of the model.
    for (slong j = 0; j < n; ++j) acb_set(initial[j], nodes[static_cast<size_t>(j)].value());
    require(run(plain, initial, nodes, residues, false), "the one-accumulator packet was rejected on its nodes");
    require(run(lanes, initial, nodes, residues, true), "the lane packet was rejected on its nodes");
    require(agree(plain, lanes, 0x1p-36), "the lane packet disagrees on its nodes");
    // One point on another point's node (the two then coincide after no step: taken beside it instead), and the
    // others beside their nodes.
    for (slong j = 0; j < n; ++j) {
        acb_set(initial[j], nodes[static_cast<size_t>(j)].value());
        ComplexBall offset;
        acb_set_d_d(offset.value(), 0x1p-20 * (2 * unit() - 1), 0x1p-20 * (2 * unit() - 1));
        acb_add(initial[j], initial[j], offset.value(), 128);
        acb_get_mid(initial[j], initial[j]);
    }
    acb_set(initial[3], nodes[7].value());
    require(run(plain, initial, nodes, residues, false), "the one-accumulator packet was rejected beside its nodes");
    require(run(lanes, initial, nodes, residues, true), "the lane packet was rejected beside its nodes");
    require(agree(plain, lanes, 0x1p-36), "the lane packet disagrees beside its nodes");
    // Two equal points: refused either way.
    acb_set(initial[5], initial[11]);
    require(!run(plain, initial, nodes, residues, false) && !run(lanes, initial, nodes, residues, true),
            "a packet with two equal points was taken");
    // Parts beyond the lanes' range (the model scaled by 2^450): the general loop, bit for bit.
    std::vector<ComplexBall> far_nodes(n), far_residues(n);
    for (slong j = 0; j < n; ++j) {
        acb_mul_2exp_si(far_nodes[static_cast<size_t>(j)].value(), nodes[static_cast<size_t>(j)].value(), 450);
        acb_mul_2exp_si(far_residues[static_cast<size_t>(j)].value(), residues[static_cast<size_t>(j)].value(), 450);
        acb_set(initial[j], far_nodes[static_cast<size_t>(j)].value());
    }
    const bool plain_far = run(plain, initial, far_nodes, far_residues, false);
    const bool lanes_far = run(lanes, initial, far_nodes, far_residues, true);
    require(plain_far == lanes_far, "the lanes changed whether a packet far from one is taken");
    if (plain_far)
        for (slong j = 0; j < n; ++j) require(acb_equal(plain[j], lanes[j]), "the lanes ran outside their range");
}
void
uncertainty_and_cancel() {
    certpolysolve::Polynomial input;
    input.set(0, certpolysolve::ball(certpolysolve::Number("-1"), certpolysolve::Number("1/16")));
    input.set(2, certpolysolve::Number("1"));
    require(input.ok(), "uncertain test input was refused");
    for (auto method : { StudyMethod::None, StudyMethod::Aberth, StudyMethod::Secular }) {
        auto p = policy(method, true);
        SolveState state(certpolysolve::detail::dense(*input.get()), p);
        cps_options o;
        o.bits = 32;
        o.max_work_bits = 128;
        o.work_units = 64;
        const auto uncertain = state.advance(&o);
        require(!uncertain.goal_met, "arithmetic erased input uncertainty");
    }
    const std::vector<reference::Complex> known{ { -3 }, { -1 }, { 2 }, { 4 }, { 7 } };
    auto polynomial = reference::input(reference::from_roots(known));
    auto p = policy(StudyMethod::Secular, true);
    SolveState state(certpolysolve::detail::dense(*polynomial.get()), p);
    cps_options o;
    o.bits = 24;
    o.max_work_bits = 1024;
    o.cancel = [](void *) { return 1; };
    require(state.advance(&o).termination == CPS_CANCELLED, "study ignored cancellation");
    o.cancel = nullptr;
    SolveState uncancelled(certpolysolve::detail::dense(*polynomial.get()), p);
    check(uncancelled.advance(&o), known, o.bits);
}
void
machine_precision_handoff(bool translate) {
    ComplexPolynomial polynomial;
    acb_poly_set_coeff_si(polynomial.value(), 0, 1);
    acb_poly_set_coeff_si(polynomial.value(), 1, 1);
    acb_poly_set_coeff_si(polynomial.value(), 5, 1);
    RootVector initial(5), first(5), second(5);
    require(initial_candidates(initial.data(), polynomial.value(), 85), "handoff initialization failed");
    auto p = policy(StudyMethod::Aberth, true);
    StudyCandidates workspace(p);
    cps_options options;
    options.max_work_bits = 512;
    certpolysolve::detail::CandidateStatistics statistics{};
    const auto machine = workspace.advance(first.data(), polynomial.value(), initial.data(), 85, 1,
                                             options, statistics, {}, false);
    require(machine.finite && p.study_statistics->machine_packets > 0 &&
                p.study_statistics->adaptive_packets == 0,
            "handoff test did not exercise a live machine stage");
    const auto machine_packets = p.study_statistics->machine_packets;
    if (translate) {
        ComplexBall shift;
        acb_one(shift.value());
        acb_poly_taylor_shift(polynomial.value(), polynomial.value(), shift.value(), 170);
        for (slong j = 0; j < first.size(); ++j)
            acb_sub_ui(first.data() + j, first.data() + j, 1, 170);
        workspace.coordinates_changed();
    }
    const auto adaptive = workspace.advance(second.data(), polynomial.value(), first.data(), 170, 1,
                                              options, statistics, {}, false);
    require(adaptive.finite && p.study_statistics->machine_packets == machine_packets &&
                p.study_statistics->adaptive_packets == 1,
            "higher precision restarted machine arithmetic instead of continuing adaptively");
    // Coordinate invalidation also preserves an already completed handoff at
    // the same precision; only the numerical preparation is discarded.
    workspace.coordinates_changed();
    workspace.advance(first.data(), polynomial.value(), second.data(), 170, 1,
                      options, statistics, {}, false);
    require(p.study_statistics->machine_packets == machine_packets &&
                p.study_statistics->adaptive_packets == 2,
            "coordinate invalidation restarted a completed machine stage");
}
void
reciprocal_newton_evaluation() {
    using C = machine::Complex<double>;
    constexpr slong degree = 2000;
    // Rectangular complex-ball Horner can accumulate wrapping even at exact
    // input points. Give this independent reference enough precision to remain
    // finite; a wide reference interval cannot check the candidate ratio.
    constexpr slong reference_bits = degree + 256;
    std::vector<C> coefficients(static_cast<size_t>(degree + 1), C(1));
    machine::PreparedNewton<double> evaluator(coefficients);
    ComplexPolynomial polynomial;
    for (slong i = 0; i <= degree; ++i)
        acb_poly_set_coeff_si(polynomial.value(), i, 1);
    for (C point : {C(0), C(.5), C(2), C(1.25, .25)}) {
        const auto sample = evaluator.evaluate(point);
        const auto ratio = sample.ratio();
        require(bool(ratio), "forward/reciprocal Newton ratio failed");
        ComplexBall z, value, derivative, expected, actual, difference;
        acb_set_d_d(z.value(), point.re, point.im);
        acb_poly_evaluate2(value.value(), derivative.value(), polynomial.value(), z.value(), reference_bits);
        acb_div(expected.value(), value.value(), derivative.value(), reference_bits);
        require(acb_is_finite(expected.value()), "independent Newton reference needs more precision");
        acb_set_d_d(actual.value(), ratio->re, ratio->im);
        acb_sub(difference.value(), actual.value(), expected.value(), reference_bits);
        Magnitude error, scale;
        acb_get_mag(error.value(), difference.value());
        acb_get_mag(scale.value(), expected.value());
        if (!(mag_get_d(error.value()) < 1e-12 * mag_get_d(scale.value())))
            std::cerr << "Newton point " << point.re << ',' << point.im << " ratio "
                      << ratio->re << ',' << ratio->im << " error " << mag_get_d(error.value())
                      << " scale " << mag_get_d(scale.value()) << '\n';
        require(mag_get_d(error.value()) < 1e-12 * mag_get_d(scale.value()),
                "Newton ratio disagrees with independent adaptive evaluation");
    }
    const auto forward = machine::jet<double, 1>(coefficients, C(2));
    require(!machine::finite(forward.value), "range example did not exceed forward binary64 evaluation");
    std::feclearexcept(FE_ALL_EXCEPT);
    std::vector<C> monomial(9);
    monomial.back() = C(1);
    machine::PreparedNewton<double> zero_factor(monomial);
    const auto ratio = zero_factor.evaluate(C(2)).ratio();
    require(ratio && machine::norm(*ratio - C(.25)) < 1e-14,
            "reversal lost the original degree of an exact zero factor");
    require(!zero_factor.evaluate(C(0)).ratio(), "undefined Newton quotient was accepted");
}
void
per_root_machine_stop() {
    ComplexPolynomial polynomial;
    acb_poly_set_coeff_si(polynomial.value(), 0, -8);
    acb_poly_set_coeff_si(polynomial.value(), 1, 14);
    acb_poly_set_coeff_si(polynomial.value(), 2, -7);
    acb_poly_set_coeff_si(polynomial.value(), 3, 1);
    RootVector seeds(3), output(3);
    acb_one(seeds.data());
    acb_set_d(seeds.data() + 1, 2.125);
    acb_set_d(seeds.data() + 2, 4.125);
    auto engine = prepare_machine_engine(IterationMethod::Aberth, ArithmeticTier::Binary64,
                                         polynomial.value(), seeds.view());
    require(bool(engine), "mixed-progress machine preparation failed");
    cps_options options;
    certpolysolve::detail::CandidateStatistics statistics{};
    const auto progress = engine->advance(8, options, statistics);
    require(progress.candidates.finite, "mixed-progress packet lost finite proposals");
    require(statistics.candidate_sweeps > 1, "one exact seed stopped the entire machine packet");
    engine->export_candidates(output.data());
    require(acb_equal(output.data(), seeds.data()), "exact seed moved after its arithmetic stop");
    for (slong i = 1; i < 3; ++i) {
        ComplexBall difference;
        acb_sub_ui(difference.value(), output.data() + i, i == 1 ? 2 : 4, 128);
        Magnitude error;
        acb_get_mag(error.value(), difference.value());
        require(mag_get_d(error.value()) < 1e-11, "remaining machine roots did not improve independently");
    }
}
// Aberth's correction where binary64 cannot hold the Newton ratio (machine_candidates.cpp). At
// p = x^2200 - 3x^1100 + 1 the value is about 1 inside the unit circle, and p'(z) = 2200 z^2199 - 3300 z^1099
// underflows: to exactly zero at z = 1/2 (2^-1099 is below binary64's least subnormal), and to a subnormal at
// z = 33/64, where p/p' (about 2^1038) overflows. Aberth's correction 1 / (p'/p - S) is defined at both. One
// starting point at either place, the others where the attempt starts them, must not end the packet for
// every point (the zero) or fail it (the overflow): the packet converges.
void
unrepresentable_newton_ratio() {
    constexpr slong n = 2200, k = 1100;
    ComplexPolynomial polynomial;
    acb_poly_set_coeff_si(polynomial.value(), 0, 1);
    acb_poly_set_coeff_si(polynomial.value(), k, -3);
    acb_poly_set_coeff_si(polynomial.value(), n, 1);
    for (const double start : {0.5, 33.0 / 64.0}) {
        RootVector seeds(n);
        require(initial_candidates(seeds.data(), polynomial.value(), 64), "no starting points");
        acb_set_d(seeds.data(), start);
        // Underflow allowed, as in the certificate attempt's packet; otherwise the zero is an arithmetic failure.
        auto engine = prepare_machine_engine(IterationMethod::Aberth, ArithmeticTier::Binary64,
                                             polynomial.value(), seeds.view(), false, false, true);
        require(bool(engine), "binary64 preparation failed");
        cps_options options;
        certpolysolve::detail::CandidateStatistics statistics{};
        const auto step = engine->advance(heuristics::machine_cover_sweeps, options, statistics);
        require(!step.exhausted, "an unrepresentable Newton ratio ended the packet");
        require(!step.arithmetic_failure, "an unrepresentable Newton ratio failed the packet");
        require(step.candidates.converged, "the packet did not converge");
        RootVector output(n);
        engine->export_candidates(output.data());
        RealBall modulus;
        acb_abs(modulus.value(), output.data(), 64);
        require(std::abs(arf_get_d(arb_midref(modulus.value()), ARF_RND_NEAR) - 1) < 0x1p-8,
                "the point did not reach the roots' circles");
    }
}
// Class E of the equivalence screen: p(3x) with p = x^3200 + x^1600 + 1.
// Rounding log2 3 to a power-of-two variable scale leaves 0.75^3200, below
// binary64, so the dyadic binary64 engine is refused. The Newton cover's
// opt-in non-dyadic scale prepares it; where the dyadic scale fits (p(2x))
// the opt-in changes nothing.
void
nondyadic_machine_scaling() {
    const slong n = 3200;
    for (const int factor : {3, 2}) {
        ComplexPolynomial polynomial;
        ComplexBall coefficient;
        for (const slong k : {slong(0), n / 2, n}) {
            acb_set_ui(coefficient.value(), static_cast<ulong>(factor));
            acb_pow_ui(coefficient.value(), coefficient.value(), static_cast<ulong>(k), 20000);
            acb_poly_set_coeff_acb(polynomial.value(), k, coefficient.value());
        }
        // Seeds on a circle a little inside the roots' modulus 1/factor.
        RootVector seeds(n);
        for (slong k = 0; k < n; ++k) {
            const double angle = 6.283185307179586 * (static_cast<double>(k) + 0.25) / static_cast<double>(n);
            acb_set_d_d(seeds.data() + k, 0.9 * std::cos(angle) / factor, 0.9 * std::sin(angle) / factor);
        }
        auto dyadic = prepare_machine_engine(IterationMethod::Aberth, ArithmeticTier::Binary64,
                                             polynomial.value(), seeds.view());
        auto scaled = prepare_machine_engine(IterationMethod::Aberth, ArithmeticTier::Binary64,
                                             polynomial.value(), seeds.view(), false, false, false, UINT64_MAX, true);
        require(bool(scaled), "non-dyadic machine scaling refused p(kx)");
        require(bool(dyadic) == (factor == 2), "dyadic machine scaling range classification changed");
        cps_options options;
        certpolysolve::detail::CandidateStatistics statistics{};
        const auto progress = scaled->advance(4, options, statistics);
        require(progress.candidates.finite, "non-dyadic machine packet lost finite proposals");
        RootVector output(n);
        scaled->export_candidates(output.data());
        for (slong k = 0; k < n; ++k) {
            RealBall modulus;
            acb_abs(modulus.value(), output.data() + k, 64);
            const double r = arf_get_d(arb_midref(modulus.value()), ARF_RND_NEAR) * factor;
            require(r > 0.5 && r < 2, "non-dyadic machine candidates left the root annulus");
        }
        if (dyadic) {
            certpolysolve::detail::CandidateStatistics other{};
            dyadic->advance(4, options, other);
            RootVector reference(n);
            dyadic->export_candidates(reference.data());
            for (slong k = 0; k < n; ++k)
                require(acb_equal(reference.data() + k, output.data() + k),
                        "opting in changed a dyadic machine preparation");
        }
    }
}
void
secular_input_equivalence() {
    // These residues give the numerator (z-1)(z-2)(z-4).
    const std::vector<certpolysolve::SecularEquation::Term> terms{
        { .residue = certpolysolve::Number("36/5"), .pole = certpolysolve::Number("-2") },
        { .residue = certpolysolve::Number("-4/3"), .pole = certpolysolve::Number("0") },
        { .residue = certpolysolve::Number("2/15"), .pole = certpolysolve::Number("3") } };
    const std::vector<reference::Complex> known{ { 1 }, { 2 }, { 4 } };
    certpolysolve::SecularEquation secular(terms);
    require(secular.ok(), "secular test input was refused");
    const auto coefficients = reference::input(reference::from_roots(known));
    for (auto method : { StudyMethod::None, StudyMethod::Aberth, StudyMethod::Secular })
        for (bool machine : { false, true }) {
            cps_options options;
            options.bits = 53;
            options.max_work_bits = 1024;
            SolveState native(*secular.get(), policy(method, machine));
            SolveState expanded(certpolysolve::detail::dense(*coefficients.get()), policy(method, machine));
            check(native.advance(&options), known, options.bits);
            check(expanded.advance(&options), known, options.bits);
        }
}
// The public entry's route for each input: native secular input keeps the
// old default driver; coefficient input goes to the selector, which sends this
// small polynomial to the maintained secular recipe (whose certificate attempt
// now serves the default accuracy goal too).
// What a solver state is made from: a polynomial's coefficient list (the input object keeps terms, C40), or a
// secular equation as it is.
certpolysolve::detail::DensePolynomial solver_input(const cps_polynomial &p) { return certpolysolve::detail::dense(p); }
const cps_secular &solver_input(const cps_secular &s) { return s; }

template<class Input>
void
public_arithmetic_route(const Input &input, const std::vector<reference::Complex> &known, EnginePolicy selected) {
    auto options = certpolysolve_test::raw_options();
    require(cps_options_set_precision(options.get(), 53) == CPS_OK && cps_options_set_max_work_bits(options.get(), 1024) == CPS_OK,
            "test options refused the precision");
    require(certpolysolve_test::set_exact_allowance(options.get(), 8000000) == CPS_OK, "test options refused the exact allowance");
    SolveState expected(solver_input(*input.get()), selected);
    const auto reference_result = expected.advance(options.get());
    require(selected.study_statistics->machine_packets > 0 || reference_result.statistics_v3.candidate_sweeps > 0,
            "public routing regression did not exercise machine arithmetic");
    check(reference_result, known, 53);
    const auto verify = [&](const certpolysolve::Result &result) {
        check(*result.get(), known, 53);
        const auto actual = result.engine_counters();
        const auto sweeps = result.statistics().sweeps;
        const auto &reference = reference_result.statistics_v3;
        // Work counters are deterministic; elapsed-time counters are excluded.
        // Check the public entry takes the qualified private path, in addition
        // to independently checking its exact roots and certificates.
        require(sweeps == reference.candidate_sweeps &&
                    actual.evaluations == reference.evaluations &&
                    actual.root_interactions == reference.root_interactions &&
                    actual.secular_regenerations == reference.secular_regenerations,
                "public entry differs from qualified candidate arithmetic");
    };
    verify(certpolysolve_test::solve(input, options.get()));
}
void
public_arithmetic_routes() {
    const std::vector<reference::Complex> known{ { -3 }, { -1 }, { 2 }, { 4 },
                                               { 7 }, { 11 }, { 13 }, { 17 } };
    auto maintained = maintained_method_policy(CPS_METHOD_SECULAR);
    maintained.study_statistics = std::make_shared<StudyStatistics>();
    public_arithmetic_route(reference::input(reference::from_roots(known)), known, maintained);
    const certpolysolve::SecularEquation secular{
        { .residue = certpolysolve::Number("36/5"), .pole = certpolysolve::Number("-2") },
        { .residue = certpolysolve::Number("-4/3"), .pole = certpolysolve::Number("0") },
        { .residue = certpolysolve::Number("2/15"), .pole = certpolysolve::Number("3") } };
    require(secular.ok(), "secular test input was refused");
    public_arithmetic_route(secular, { { 1 }, { 2 }, { 4 } }, policy(StudyMethod::None, true));
}
int
main() {
    certpolysolve::ThreadScope scope;
    try {
        secular_packet();
        lane_sums();
        known_roots();
        uncertainty_and_cancel();
        machine_precision_handoff(false);
        machine_precision_handoff(true);
        reciprocal_newton_evaluation();
        per_root_machine_stop();
        unrepresentable_newton_ratio();
        nondyadic_machine_scaling();
        secular_input_equivalence();
        public_arithmetic_routes();
        std::cout << "D1 independent mathematical checks passed\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
