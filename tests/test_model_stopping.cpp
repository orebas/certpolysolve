#include "solve/policies.hpp"
#include "internal.hpp"
#include <cmath>
#include <flint/fmpq_poly.h>
#include "certpolysolve/certpolysolve.hpp"
#include "numeric/candidates/study_candidates.hpp"
#include <iostream>
#include <stdexcept>
#include <string_view>

using namespace certpolysolve::detail;
static void require(bool ok, const char *message) {
    if (!ok)
        throw std::runtime_error(message);
}
static int cancelled(void *) { return 1; }

static void check_classical_floor() {
    // At z=1+2^-8, p=(z-1)^8-2^-80 is about 2^-64. Its expanded
    // coefficients suffer severe cancellation in 53-bit Horner arithmetic.
    // A paused proposal must resume once the caller provides enough precision.
    ComplexPolynomial linear, polynomial;
    acb_poly_set_coeff_si(linear.value(), 0, -1);
    acb_poly_set_coeff_si(linear.value(), 1, 1);
    acb_poly_pow_ui(polynomial.value(), linear.value(), 8, 192);
    ComplexBall epsilon;
    acb_one(epsilon.value());
    acb_mul_2exp_si(epsilon.value(), epsilon.value(), -80);
    acb_sub(polynomial.value()->coeffs, polynomial.value()->coeffs, epsilon.value(), 192);
    RootVector initial(8), low(8), high(8);
    acb_one(initial[0]);
    acb_one(epsilon.value());
    acb_mul_2exp_si(epsilon.value(), epsilon.value(), -8);
    acb_add(initial[0], initial[0], epsilon.value(), 192);
    for (slong j = 1; j < 8; ++j) {
        acb_set_si(initial[j], j + 3);
        arb_set_si(acb_imagref(initial[j]), j);
    }
    std::vector<uint8_t> inactive(8, 1);
    inactive[0] = 0;
    EnginePolicy policy;
    policy.classical_error_stopping = true;
    CandidateWorkspace workspace(&policy);
    cps_options options;
    certpolysolve::detail::CandidateStatistics statistics{};
    const auto floor = candidate_packet(low.data(), polynomial.value(), initial.data(), 53, 8,
        CandidateEngine::Aberth, options, statistics, inactive, false, &workspace);
    require(floor.finite && !floor.converged &&
                floor.continuation == CandidateContinuation::PrecisionDecision,
            "a cancellation floor must request precision, not claim convergence");
    require(acb_equal(low[0], initial[0]), "unreliable Newton evaluation moved a proposal");
    const auto evaluations = statistics.evaluations;
    const auto retained = candidate_packet(low.data(), polynomial.value(), initial.data(), 53, 8,
        CandidateEngine::Aberth, options, statistics, inactive, false, &workspace);
    require(retained.continuation == CandidateContinuation::PrecisionDecision &&
                statistics.evaluations == evaluations,
            "unchanged paused classical point was evaluated again");
    // Replacing a coordinate at the same precision invalidates the temporary
    // mask. In particular, root reordering by inclusion cannot reuse it blindly.
    acb_set_si(initial[0], -2);
    const auto replaced = candidate_packet(high.data(), polynomial.value(), initial.data(), 53, 1,
        CandidateEngine::Aberth, options, statistics, inactive, false, &workspace);
    require(replaced.finite && !acb_equal(high[0], initial[0]), "replacement reused a stale arithmetic mask");
    acb_set(initial[0], low[0]);
    const auto resumed = candidate_packet(high.data(), polynomial.value(), low.data(), 192, 1,
        CandidateEngine::Aberth, options, statistics, inactive, false, &workspace);
    require(resumed.finite && !acb_equal(high[0], low[0]), "higher precision did not resume a paused proposal");
    for (slong j = 1; j < 8; ++j)
        require(acb_equal(high[j], initial[j]), "classical arithmetic stop moved an inactive point");
    options.cancel = cancelled;
    const auto stopped = candidate_packet(high.data(), polynomial.value(), low.data(), 192, 1,
        CandidateEngine::Aberth, options, statistics, inactive, false, &workspace);
    require(stopped.cancelled, "classical arithmetic stop hid cancellation");
}

// Asked to (a mixed round's binary64 sweeps, EnginePolicy::machine_keeps_unmoved), a node the binary64 packet
// does not move keeps its entry point. The model's root 1 + 2^-80 is 1 in binary64: the node settles at its first
// evaluation, and its 81 bits must survive the packet; an ordinary packet writes back its rounding.
static void check_unmoved_entry() {
    const slong precision = 256;
    ComplexBall root, offset;
    acb_one(offset.value());
    acb_mul_2exp_si(offset.value(), offset.value(), -80);
    acb_one(root.value());
    acb_add(root.value(), root.value(), offset.value(), precision);
    // (x-r)(x-2) at the nodes 0 and 3: a_j = -p(b_j) / prod_{k != j} (b_j - b_k).
    std::vector<ComplexBall> nodes(2), residues(2);
    acb_zero(nodes[0].value());
    acb_set_si(nodes[1].value(), 3);
    acb_mul_2exp_si(residues[0].value(), root.value(), 1);
    acb_div_ui(residues[0].value(), residues[0].value(), 3, precision);
    acb_sub_ui(residues[1].value(), root.value(), 3, precision);
    acb_div_ui(residues[1].value(), residues[1].value(), 3, precision);
    RootVector initial(2), out(2);
    acb_set(initial[0], root.value());
    acb_set_si(initial[1], 2);
    std::vector<uint8_t> certified(2, 0);
    std::vector<slong> bits(2, 85);
    cps_options options;
    options.bits = 53;
    options.max_work_bits = 512;
    EnginePolicy policy = production_policy();
    policy.model_stopping = SecularModelPolicy::Continuation;
    for (bool keep : {true, false}) {
        certpolysolve::detail::CandidateStatistics statistics{};
        ModelStoppingState model;
        SecularMachineWorkspace machine;
        CandidateProgress progress;
        require(machine.advance(progress, out.data(), initial.data(), nodes, residues, certified, bits, 1, 8,
                                options, statistics, policy, &model, keep) &&
                    model.settled(0),
                "the model's root did not settle in binary64");
        if (keep)
            require(acb_equal(out[0], initial[0]), "a node the binary64 packet did not move lost its entry bits");
        else
            require(!acb_equal(out[0], initial[0]), "the entry fits binary64: the check above proves nothing");
    }
}

// A converged cubic step needs no confirming sweep (candidates.cpp). On the
// model of (x-1)(x-2) at nodes 0 and 3, from points 2^-20 off the roots at 160
// bits, a point's first Aberth step leaves it about 2^-60 off and its second
// about 2^-180, below the arithmetic's resolution (2^-158). The second step's
// order (60 / 20 = 3) predicts that, so each point settles after it: two
// evaluations a point, where a confirming third sweep made three.
static void cubic_steps_settle() {
    std::vector<ComplexBall> nodes(2), residues(2);
    acb_zero(nodes[0].value());
    acb_set_si(nodes[1].value(), 3);
    acb_set_si(residues[0].value(), 2);
    acb_div_ui(residues[0].value(), residues[0].value(), 3, 192);
    acb_neg(residues[1].value(), residues[0].value());
    RootVector initial(2), out(2);
    acb_set_d(initial[0], 1.0 + std::ldexp(1.0, -20)); // exact in binary64
    acb_set_d(initial[1], 2.0 - std::ldexp(1.0, -20));
    std::vector<uint8_t> certified(2, 0);
    std::vector<slong> bits(2, 160);
    cps_options options;
    options.bits = 53;
    options.max_work_bits = 512;
    certpolysolve::detail::CandidateStatistics statistics{};
    EnginePolicy policy = production_policy();
    policy.model_stopping = SecularModelPolicy::Continuation;
    ModelStoppingState model;
    const auto progress = secular_candidate_packet(out.data(), initial.data(), nodes, residues, certified, bits, 8,
                                                   options, statistics, nullptr, false, true, nullptr, 1, &policy,
                                                   &model);
    require(progress.finite && progress.model_exhausted, "the cubic model roots did not settle");
    require(statistics.evaluations == 4, "a converged cubic step was confirmed by another sweep");
    for (slong k = 0; k < 2; ++k) {
        ComplexBall error;
        acb_sub_si(error.value(), out[k], k + 1, 192);
        mag_t bound;
        mag_init(bound);
        acb_get_mag(bound, error.value());
        require(mag_cmp_2exp_si(bound, -150) < 0, "a settled point is not at its root");
        mag_clear(bound);
    }
}
// The steps a binary64 secular packet cannot take, taken at a fresh model's floor (EnginePolicy::machine_node_offsets,
// NodeOffsets), on the model of (x - 1/3)(x - 2/3), whose roots binary64 cannot hold (no packet lands on them by
// rounding). Nodes 2^-30 off the roots: the packet moves them, there is no floor, and its output is the ordinary
// packet's, bit for bit. Nodes at binary64's floor (1/3 and 2/3 rounded): the ordinary packet settles each at its
// first evaluation without a step (the floor that sends the solve up the precision ladder); with the policy the
// sub-floor steps follow, and the points come out about 2^-107 off with the floor still signalled.
static void node_offsets_refine() {
    const slong precision = 256;
    std::vector<ComplexBall> roots(2), nodes(2), residues(2);
    for (int j = 0; j < 2; ++j) {
        acb_set_si(roots[j].value(), j + 1);
        acb_div_ui(roots[j].value(), roots[j].value(), 3, precision);
    }
    RootVector initial(2), out(2), ordinary(2);
    const auto model_of = [&](double offset) {
        for (int j = 0; j < 2; ++j) {
            acb_set_d(nodes[j].value(), offset);
            if (offset == 0) acb_set_d(nodes[j].value(), (j + 1) / 3.0); // binary64's floor
            else acb_add(nodes[j].value(), roots[j].value(), nodes[j].value(), precision);
        }
        for (int j = 0; j < 2; ++j) {
            // w_j = -p(b_j) / (b_j - b_k) for p(x) = (x - r_0)(x - r_1).
            ComplexBall value, factor, between;
            acb_sub(value.value(), nodes[j].value(), roots[0].value(), precision);
            acb_sub(factor.value(), nodes[j].value(), roots[1].value(), precision);
            acb_mul(value.value(), value.value(), factor.value(), precision);
            acb_sub(between.value(), nodes[j].value(), nodes[1 - j].value(), precision);
            acb_div(residues[j].value(), value.value(), between.value(), precision);
            acb_neg(residues[j].value(), residues[j].value());
            acb_set(initial[j], nodes[j].value());
        }
    };
    std::vector<uint8_t> certified(2, 0);
    std::vector<slong> bits(2, SecularMachineWorkspace::stage_bits(53));
    cps_options options;
    options.bits = 53;
    options.max_work_bits = 512;
    const auto error_log2 = [&](int j) {
        ComplexBall error;
        acb_sub(error.value(), out[j], roots[j].value(), precision);
        mag_t bound;
        mag_init(bound);
        acb_get_mag(bound, error.value());
        const double result = mag_is_zero(bound) ? -1e9 : mag_get_d_log2_approx(bound);
        mag_clear(bound);
        return result;
    };
    const auto packet = [&](bool offsets, CandidateProgress &progress) {
        EnginePolicy policy = production_policy();
        policy.model_stopping = SecularModelPolicy::Continuation;
        policy.machine_node_offsets = offsets;
        certpolysolve::detail::CandidateStatistics statistics{};
        ModelStoppingState model;
        SecularMachineWorkspace machine;
        require(machine.advance(progress, out.data(), initial.data(), nodes, residues, certified, bits, 1, 8,
                                options, statistics, policy, &model, true, policy.machine_node_offsets) &&
                    progress.finite,
                "the binary64 secular packet refused the model");
    };
    model_of(std::ldexp(1.0, -30)); // the second node is placed at r_1 + 2^-30
    CandidateProgress progress;
    packet(false, progress);
    require(!progress.model_immediate && std::max(error_log2(0), error_log2(1)) > -60,
            "nodes 2^-30 off: the ordinary packet should move them and stop at binary64's floor");
    for (int j = 0; j < 2; ++j) acb_set(ordinary[j], out[j]);
    packet(true, progress);
    for (int j = 0; j < 2; ++j)
        require(acb_equal(out[j], ordinary[j]), "a packet that does not end at the floor changed with the policy");
    model_of(0);
    for (bool offsets : {false, true}) {
        packet(offsets, progress);
        require(progress.model_exhausted && progress.model_immediate, "nodes at binary64's floor did not settle at once");
        const double worst = std::max(error_log2(0), error_log2(1));
        if (offsets)
            require(worst < -100, "the sub-floor steps at the floor were not taken");
        else
            require(worst > -60, "the ordinary packet moved past binary64's floor: the check above proves nothing");
    }
}
// The floor steps (EnginePolicy::machine_node_offsets) are a default-goal measure: a higher goal's staged default
// solve keeps its multiprecision packet, whose points start that goal's certificate attempt (floor-stepped points made
// nektarios at 3,322 bits 11.8 -> 15.2 s). On Wilkinson's 20 plus one, the default goal takes floor steps; at 3,322
// bits the staged solve reaches the same floor, takes none, and returns exactly what it returns without the policy.
static void floor_steps_default_goal_only() {
    // prod_{k=1}^{20} (x - k) + 1, coefficients ascending.
    IntegerPolynomial wilkinson, factor;
    fmpz_poly_one(wilkinson.value());
    for (slong k = 1; k <= 20; ++k) {
        fmpz_poly_zero(factor.value());
        fmpz_poly_set_coeff_si(factor.value(), 1, 1);
        fmpz_poly_set_coeff_si(factor.value(), 0, -k);
        fmpz_poly_mul(wilkinson.value(), wilkinson.value(), factor.value());
    }
    Integer constant;
    fmpz_poly_get_coeff_fmpz(constant.value(), wilkinson.value(), 0);
    fmpz_add_ui(constant.value(), constant.value(), 1);
    fmpz_poly_set_coeff_fmpz(wilkinson.value(), 0, constant.value());
    struct Run {
        cps_result result;
        unsigned floor_steps = 0;
    };
    const auto solve = [&](bool offsets, uint64_t bits) {
        certpolysolve::detail::DensePolynomial input;
        input.coefficients.resize(static_cast<size_t>(fmpz_poly_length(wilkinson.value())));
        for (slong k = 0; k < fmpz_poly_length(wilkinson.value()); ++k)
            fmpq_set_fmpz(input.coefficients[static_cast<size_t>(k)].real.significand.value(),
                          wilkinson.value()->coeffs + k);
        auto policy = maintained_method_policy(CPS_METHOD_SECULAR);
        require(policy.machine_node_offsets, "the secular recipe no longer takes floor steps");
        policy.machine_node_offsets = offsets;
        Run run;
        policy.trace_data = &run.floor_steps;
        policy.trace = [](void *data, const TraceEvent &event) {
            if (std::string_view(event.stage) == "study_secular_candidates" &&
                std::string_view(event.reason) == "floor-steps")
                ++*static_cast<unsigned *>(data);
        };
        SolveState state(std::move(input), policy);
        cps_options options;
        options.method = CPS_METHOD_SECULAR;
        options.goal = CPS_ROOT_APPROXIMATION;
        options.bits = bits;
        options.max_work_bits = 65536;
        for (unsigned calls = 0; calls < 200 && !run.result.goal_met; ++calls) run.result = state.advance(&options);
        require(run.result.goal_met, "Wilkinson's 20 plus one did not finish");
        return run;
    };
    require(solve(true, cps_options{}.bits).floor_steps > 0, "the default goal took no floor steps: the fixture proves nothing");
    const auto with = solve(true, 3322), without = solve(false, 3322);
    require(with.floor_steps == 0, "a higher goal's staged default solve took floor steps");
    require(with.result.roots.size() == without.result.roots.size(), "the policy changed a higher goal's cover");
    for (size_t i = 0; i < with.result.roots.size(); ++i)
        require(acb_equal(with.result.roots[i].box.value(), without.result.roots[i].box.value()) &&
                    with.result.roots[i].count == without.result.roots[i].count,
                "the policy changed a higher goal's result");
}
// A node a validation returns to the active set takes the precision growths it missed while certified
// (SecularRefinement::missed_growths_). The fixture is the matrix input jacobi-characteristic-3-affine-cluster: the
// characteristic polynomial of a 256 x 256 symmetric tridiagonal matrix (diagonal ((11 i mod 19) - 9)/8, off-diagonal
// products (((i mod 3) + 1)/5)^2; benchmarks/consolidation/engineering.py), at 256 x - 256. At the accuracy goal, the
// floor steps certify nodes before the first growth and a later validation returns one; left at binary64's stage
// among nodes grown past it, it kept every later round a mixed one (25 binary64 packets after the growth; 1.12x).
struct RationalPolynomial { // a local owner for the fixture's exact recurrence
    fmpq_poly_t poly;
    RationalPolynomial() { fmpq_poly_init(poly); }
    ~RationalPolynomial() { fmpq_poly_clear(poly); }
    RationalPolynomial(const RationalPolynomial &) = delete;
    RationalPolynomial &operator=(const RationalPolynomial &) = delete;
    fmpq_poly_struct *value() { return poly; }
};
static void returned_nodes_take_missed_growths() {
    const slong n = 256;
    RationalPolynomial previous, current, next, linear;
    fmpq_poly_one(previous.value());
    Rational value;
    const auto diagonal = [&](slong i) {
        fmpq_set_si(value.value(), (i * 11) % 19 - 9, 8);
    };
    diagonal(0);
    fmpq_poly_set_coeff_si(current.value(), 1, 1);
    fmpq_neg(value.value(), value.value());
    fmpq_poly_set_coeff_fmpq(current.value(), 0, value.value());
    for (slong k = 1; k < n; ++k) {
        // next = (x - d_k) current - b_{k-1} previous
        fmpq_poly_zero(linear.value());
        fmpq_poly_set_coeff_si(linear.value(), 1, 1);
        diagonal(k);
        fmpq_neg(value.value(), value.value());
        fmpq_poly_set_coeff_fmpq(linear.value(), 0, value.value());
        fmpq_poly_mul(next.value(), linear.value(), current.value());
        fmpq_set_si(value.value(), ((k - 1) % 3) + 1, 5);
        fmpq_mul(value.value(), value.value(), value.value());
        fmpq_poly_scalar_mul_fmpq(previous.value(), previous.value(), value.value());
        fmpq_poly_sub(next.value(), next.value(), previous.value());
        fmpq_poly_swap(previous.value(), current.value());
        fmpq_poly_swap(current.value(), next.value());
    }
    fmpq_poly_zero(linear.value());
    fmpq_poly_set_coeff_si(linear.value(), 1, 256);
    fmpq_poly_set_coeff_si(linear.value(), 0, -256);
    fmpq_poly_compose(next.value(), current.value(), linear.value());
    certpolysolve::detail::DensePolynomial input;
    input.coefficients.resize(static_cast<size_t>(fmpq_poly_length(next.value())));
    for (slong k = 0; k < fmpq_poly_length(next.value()); ++k)
        fmpq_poly_get_coeff_fmpq(input.coefficients[static_cast<size_t>(k)].real.significand.value(), next.value(), k);
    auto policy = maintained_method_policy(CPS_METHOD_SECULAR);
    struct Seen {
        bool grown = false;
        unsigned binary64_after_growth = 0;
    } seen;
    policy.trace_data = &seen;
    policy.trace = [](void *data, const TraceEvent &event) {
        auto &s = *static_cast<Seen *>(data);
        if (std::string_view(event.stage) == "secular_precision_growth") s.grown = true;
        if (s.grown && std::string_view(event.stage) == "study_secular_candidates" &&
            std::string_view(event.tier) == "binary64")
            ++s.binary64_after_growth;
    };
    SolveState state(std::move(input), policy);
    cps_options options;
    options.method = CPS_METHOD_SECULAR;
    options.goal = CPS_ROOT_ACCURACY;
    options.bits = cps_options{}.bits;
    options.max_work_bits = 65536;
    cps_result result;
    for (unsigned calls = 0; calls < 400 && !result.goal_met; ++calls) result = state.advance(&options);
    require(result.goal_met, "the clustered Jacobi characteristic polynomial did not finish");
    require(seen.grown, "the fixture no longer grows precision: the check below proves nothing");
    require(seen.binary64_after_growth == 0, "a node returned to the active set below its peers' precision");
}
static void check() {
    check_classical_floor();
    check_unmoved_entry();
    cubic_steps_settle();
    node_offsets_refine();
    floor_steps_default_goal_only();
    returned_nodes_take_missed_growths();
    const double infinity = std::numeric_limits<double>::infinity();
    require(model_stop(-60, 0, -10, 0, 2, 53) == ModelStop::EvaluationFloor,
            "cancellation must stop model work without claiming a small Newton correction");
    require(model_stop(NAN, NAN, -60, 0, 2, 53) == ModelStop::SmallCorrection,
            "pole limit must admit correction stopping");
    require(model_stop(NAN, NAN, -60, -infinity, 2, 53) == ModelStop::Active,
            "a small absolute correction at zero is not a relative floor");

    std::vector<ComplexBall> nodes(2), residues(2);
    acb_zero(nodes[0].value());
    acb_set_si(nodes[1].value(), 3);
    acb_set_si(residues[0].value(), 2);
    acb_div_ui(residues[0].value(), residues[0].value(), 3, 192);
    acb_neg(residues[1].value(), residues[0].value());
    // These residues represent (x-1)(x-2) at nodes 0 and 3. Its known roots
    // test stopping on the rounded model independently of any solver output.
    RootVector initial(2), out(2);
    acb_one(initial[0]);
    acb_set_si(initial[1], 2);
    std::vector<uint8_t> certified(2, 0);
    std::vector<slong> bits(2, 85);
    cps_options options;
    options.bits = 53;
    options.max_work_bits = 512;
    certpolysolve::detail::CandidateStatistics statistics{};
    EnginePolicy policy = production_policy();
    policy.model_stopping = SecularModelPolicy::Continuation;
    ModelStoppingState model;
    SecularMachineWorkspace machine;
    CandidateProgress progress;
    require(machine.advance(progress, out.data(), initial.data(), nodes, residues, certified, bits, 1, 8, options,
                            statistics, policy, &model),
            "machine model rejected");
    require(progress.model_exhausted && model.exhausted(certified), "known model roots did not settle");
    require(certified == std::vector<uint8_t>(2, 0), "temporary stopping changed certificate flags");
    const auto evaluations = statistics.evaluations;
    require(machine.advance(progress, out.data(), out.data(), nodes, residues, certified, bits, 1, 8, options,
                            statistics, policy, &model),
            "retained model rejected");
    require(statistics.evaluations == evaluations, "settled model points were evaluated again");
    require(machine.advance(progress, out.data(), initial.data(), nodes, residues, certified, bits, 2, 8, options,
                            statistics, policy, &model),
            "new model rejected");
    require(statistics.evaluations > evaluations, "model revision did not reset temporary stopping");

    // Changing arithmetic must discard the machine mask. The adaptive kernel
    // starts at interpolation poles, exercising the removable Newton limits.
    acb_set(initial[0], nodes[0].value());
    acb_set(initial[1], nodes[1].value());
    bits.assign(2, 160);
    const auto before = statistics.evaluations;
    progress = secular_candidate_packet(out.data(), initial.data(), nodes, residues, certified, bits, 8, options,
                                        statistics, nullptr, false, true, nullptr,2, &policy, &model);
    require(progress.finite && statistics.evaluations > before, "arithmetic transition reused stale mask");

    // A certified point is immutable even if it is not the numerical model's
    // root. The kernel cannot change a proof owner's restart coordinates.
    certified[0] = 1;
    acb_set_si(initial[0], 7);
    progress = secular_candidate_packet(out.data(), initial.data(), nodes, residues, certified, bits, 3, options,
                                        statistics, nullptr, false, true, nullptr,3, &policy, &model);
    require(acb_equal(out[0], initial[0]), "certified point moved");
    options.cancel = cancelled;
    progress = secular_candidate_packet(out.data(), initial.data(), nodes, residues, certified, bits, 3, options,
                                        statistics, nullptr, false, true, nullptr,4, &policy, &model);
    require(progress.cancelled, "model stopping hid cancellation");
}
int main() {
    certpolysolve::ThreadScope scope;
    try {
        check();
        std::cout << "temporary model stopping and ownership passed\n";
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
