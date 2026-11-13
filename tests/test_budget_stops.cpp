#include "solve/policies.hpp"
// Honest stops at the precision budget: a secular model that cannot steer its
// nodes ends the solve for precision instead of moving them by noise.
#include "exact_winding_api.hpp"
#include "certify/counting.hpp"
#include "internal.hpp"
#include "methods/secular/secular_refinement.hpp"
#include "numeric/candidates/study_candidates.hpp"
#include <iostream>
#include <string>

namespace d = certpolysolve::detail;
using namespace reference;

namespace {
uint64_t total(const certpolysolve::Result &result) {
    uint64_t sum = 0;
    for (const auto &region : result.roots()) sum += region.count;
    return sum;
}
// A ball mid * 2^mid_exponent with radius 2^radius_exponent on the real axis.
void real_ball(acb_t out, slong mid, slong mid_exponent, slong radius_exponent) {
    acb_set_si(out, mid);
    acb_mul_2exp_si(out, out, mid_exponent);
    mag_set_ui_2exp_si(arb_radref(acb_realref(out)), 1, radius_exponent);
}
// The steering test itself. Nodes 0 and 2^-30 share a neighbourhood, node 1
// stands apart. The residue at 2^-30 has no correct bit, so the step at node
// 0, a_0 / (sum_j (a_j + a_0)/(b_0 - b_j) - 1), has an unbounded enclosure
// although a_0 is known to 30 bits. Node 1 is steered: its neighbours'
// uncertainty is 2^-28 of its denominator.
void steering() {
    std::vector<d::ComplexBall> nodes(3), residues(3);
    acb_zero(nodes[0].value());
    acb_one(nodes[1].value());
    acb_one(nodes[2].value());
    acb_mul_2exp_si(nodes[2].value(), nodes[2].value(), -30);
    real_ball(residues[0].value(), 1, -30, -60);
    real_ball(residues[1].value(), 1, -2, -60);
    real_ball(residues[2].value(), 1, -30, -28);
    const std::vector<uint8_t> active(3, 0);
    require(d::unsteered_nodes(nodes, residues, active) == std::vector<uint8_t>{1, 0, 1},
            "an unresolved neighbour did not leave its cluster unsteered, or a separate node lost its step");
    require(d::unsteered_nodes(nodes, residues, std::vector<uint8_t>{1, 0, 0}) == std::vector<uint8_t>{0, 0, 1},
            "an inactive node was given a steering verdict");
    // The same geometry with every residue resolved steers every node.
    real_ball(residues[2].value(), 1, -30, -60);
    require(d::unsteered_nodes(nodes, residues, active) == std::vector<uint8_t>(3, 0),
            "a model of accurate residues was called unsteered");
    // An exact zero residue is an exact root, not an unsteered node.
    acb_zero(residues[0].value());
    require(d::unsteered_nodes(nodes, residues, active)[0] == 0, "an exact root was called unsteered");

    // Both packet engines settle an unsteered node at its first evaluation and
    // leave it where it is; without the verdict they step it.
    real_ball(residues[0].value(), 1, -30, -60);
    cps_options options{};
    options.bits = 53;
    options.max_work_bits = 4096;
    auto policy = d::maintained_method_policy(CPS_METHOD_SECULAR);
    const std::vector<slong> bits(3, 64);
    d::RootVector initial(3), out(3);
    for (slong i = 0; i < 3; ++i) acb_set(initial.data() + i, nodes[static_cast<size_t>(i)].value());
    for (const bool machine : {false, true})
        for (const bool verdict : {false, true}) {
            d::ModelStoppingState model;
            const uint64_t revision = 7;
            if (verdict) model.steering(revision, {1, 0, 0});
            certpolysolve::detail::CandidateStatistics statistics{};
            d::CandidateProgress progress;
            if (machine) {
                d::SecularMachineWorkspace workspace;
                require(workspace.advance(progress, out.data(), initial.data(), nodes, residues, active, bits, revision,
                                          1, options, statistics, policy, &model),
                        "the binary64 packet rejected the model");
            } else
                progress = d::secular_candidate_packet(out.data(), initial.data(), nodes, residues, active, bits, 1,
                                                       options, statistics, nullptr, false, true, nullptr,
                                                       revision, &policy, &model);
            require(progress.finite, "the packet was not finite");
            require(model.settled(0) == verdict && acb_equal(out.data(), initial.data()) == verdict,
                    "a packet engine ignored the steering verdict, or stopped a steered node");
        }
}
// The step at a node whose Newton denominator D_k rounds to exactly zero
// (NewtonPole in model_stopping.hpp). Nodes 2, 3, 4; at node 0, D_0 =
// -(a_1 + a_0) - (a_2 + a_0)/2 - 1 and A_0 = -3/2, so the Aberth limit -1/A
// moves it from 2 to 4/3. Three residue sets, all dyadic:
//  (a) a = (1/4, -3/4, -5/4): D_0 = 0 exactly and a_0 A_0 stands far above
//      D_0's rounding error. The node must move by -1/A_0; the packet used
//      to be rejected as a whole (binary64) or the node left unmoved
//      (multiprecision).
//  (b) a = (2^-80, -3/4, -1/2 - 2^-69): D_0 = 2^-70 - 3 2^-81 rounds to zero
//      at 53 and 64 bits and a_0 is below that rounding error. No step is
//      determined: the node settles at the evaluation floor. At 192 bits D_0
//      is resolved and the node steps (by about 2^-10).
//  (c) a = (1/4, -3/4, -5/4 + 2^-48): D_0 = -2^-49, tiny but not zero, the
//      shape of a node inside a curve of roots where p is nearly constant
//      (sendra320). Newton's correction is huge and the Aberth step is -1/A_0
//      to 2^-47; the node must take it. (An error floor on D_0 settles it
//      instead: that floor was tried and withdrawn, 30 September.)
// Node 1 has an ordinary step in every case and must take it.
void limit_step_pole() {
    cps_options options{};
    options.bits = 53;
    options.max_work_bits = 4096;
    auto policy = d::maintained_method_policy(CPS_METHOD_SECULAR);
    const std::vector<uint8_t> active(3, 0);
    struct Case { const char *name; slong a0_mid, a0_exp, a2_mid, a2_exp, a2_offset_mid, a2_offset_exp; };
    const Case cases[] = {{"tiny nonzero D", 1, -2, -5, -2, 1, -48}, {"resolved pole", 1, -2, -5, -2, 0, 0},
                          {"unresolved pole", 1, -80, -1, -1, -1, -69}};
    d::ComplexBall target; // 4/3
    acb_set_ui(target.value(), 4);
    acb_div_ui(target.value(), target.value(), 3, 256);
    d::RealBall distance;
    for (const auto &c : cases) {
        std::vector<d::ComplexBall> nodes(3), residues(3);
        acb_set_si(nodes[0].value(), 2);
        acb_set_si(nodes[1].value(), 3);
        acb_set_si(nodes[2].value(), 4);
        real_ball(residues[0].value(), c.a0_mid, c.a0_exp, -300);
        real_ball(residues[1].value(), -3, -2, -300);
        real_ball(residues[2].value(), c.a2_mid, c.a2_exp, -300);
        if (c.a2_offset_mid) {
            d::ComplexBall offset;
            acb_set_si(offset.value(), c.a2_offset_mid);
            acb_mul_2exp_si(offset.value(), offset.value(), c.a2_offset_exp);
            acb_add(residues[2].value(), residues[2].value(), offset.value(), 256);
        }
        d::RootVector initial(3), out(3);
        for (slong i = 0; i < 3; ++i) acb_set(initial.data() + i, nodes[static_cast<size_t>(i)].value());
        // (engine, working bits, does node 0 step)
        struct Arithmetic { bool machine; slong bits; bool steps; };
        const bool unresolved = std::string(c.name) == "unresolved pole";
        for (const auto arithmetic : {Arithmetic{true, 64, !unresolved}, Arithmetic{false, 64, !unresolved},
                                      Arithmetic{false, 192, true}}) {
            d::ModelStoppingState model;
            const std::vector<slong> bits(3, arithmetic.bits);
            certpolysolve::detail::CandidateStatistics statistics{};
            d::CandidateProgress progress;
            if (arithmetic.machine) {
                d::SecularMachineWorkspace workspace;
                require(workspace.advance(progress, out.data(), initial.data(), nodes, residues, active, bits, 3, 1,
                                          options, statistics, policy, &model),
                        (std::string(c.name) + ": the binary64 packet was rejected as a whole").c_str());
            } else
                progress = d::secular_candidate_packet(out.data(), initial.data(), nodes, residues, active, bits, 1,
                                                       options, statistics, nullptr, false, true, nullptr, 3,
                                                       &policy, &model);
            require(progress.finite, (std::string(c.name) + ": the packet was not finite").c_str());
            const bool stayed = acb_equal(out.data(), initial.data()) != 0;
            require(model.settled(0) == !arithmetic.steps && stayed == !arithmetic.steps,
                    (std::string(c.name) + (arithmetic.steps ? ": node 0 did not take its step"
                                                             : ": node 0 stepped where no step is determined")).c_str());
            if (arithmetic.steps && !unresolved) {
                acb_sub(out.data(), out.data(), target.value(), 256);
                acb_abs(distance.value(), out.data(), 256);
                d::RealBall tolerance;
                arb_one(tolerance.value());
                arb_mul_2exp_si(tolerance.value(), tolerance.value(), -40);
                require(arb_lt(distance.value(), tolerance.value()),
                        (std::string(c.name) + ": node 0 did not take the Aberth step -1/A").c_str());
            }
            require(!model.settled(1) && !acb_equal(out.data() + 1, initial.data() + 1),
                    (std::string(c.name) + ": node 1's ordinary step was not taken").c_str());
        }
    }
}
// Wilkinson-41 scaled by 2^-6 under a 32-bit budget (fuzz seed
// 10794763200463370038): 32 bits resolve no residue, 18 of the 41 nodes moved
// by noise every round, and the solve ran to the fuzzer's 30 s deadline
// (20,000 units). It must stop for precision, with its cover, long before the
// work budget.
void wilkinson_below_its_precision() {
    std::vector<Complex> known;
    for (int k = 1; k <= 41; ++k) known.emplace_back(Rational(k) / 64);
    const auto p = input(from_roots(known));
    certpolysolve::Options o;
    const uint64_t budget = 2000;
    o.method = certpolysolve::Method::secular;
    o.goal = certpolysolve::Goal::accuracy();
    o.bits = 53;
    o.max_bits = 32;
    o.work_budget = budget;
    o.cell_budget = CPS_DEFAULT_LIVE_CELLS;
    const auto result = certpolysolve::solve(p, o);
    require(cps_result_termination(result.get()) == CPS_PRECISION_LIMIT &&
                cps_result_stop_reason(result.get()) == CPS_STOP_PRECISION,
            "secular below its precision did not stop for precision");
    require(result.complete() && total(result) == 41, "the precision stop lost its cover");
    require(result.statistics().work_units < budget / 10, "the precision stop came late");
}
// One root and five more 2^-100 apart (the shape of fuzz seed
// 12828861162913203822), 100-bit accuracy. A 96-bit budget cannot separate
// the five (nor can a 192-bit one, unless the private control
// exact_residues_at_budget is on; at 192 bits the solve once ran out its work
// budget, 100,000 units in 7.8 s). It must stop for precision at 96 bits and
// succeed at 384.
void cluster_below_its_precision() {
    std::vector<Complex> known{Complex(Rational(7) / 8, Rational(-1) / 11)};
    for (int j = 0; j < 5; ++j)
        known.emplace_back(Rational(-5) / 8 + j * power2(-100), Rational(3) / 5 + (j % 2) * power2(-101));
    const auto p = input(from_roots(known));
    const uint64_t budget = 2000;
    for (const auto method : {certpolysolve::Method::secular, certpolysolve::Method::automatic}) {
        certpolysolve::Options o;
        o.method = method;
        o.goal = certpolysolve::Goal::accuracy();
        o.bits = 100;
        o.max_bits = 96;
        o.work_budget = budget;
        o.cell_budget = CPS_DEFAULT_LIVE_CELLS;
        const auto limited = certpolysolve::solve(p, o);
        require(cps_result_termination(limited.get()) == CPS_PRECISION_LIMIT &&
                    cps_result_stop_reason(limited.get()) == CPS_STOP_PRECISION &&
                    limited.complete() && total(limited) == 6,
                "a cluster below the precision budget did not stop for precision with its cover");
        require(limited.statistics().work_units < budget / 10, "the cluster's precision stop came late");
        o.max_bits = 384;
        const auto enough = certpolysolve::solve(p, o);
        require(enough.goal_met() && enough.roots().size() == 6, "the cluster failed with enough precision");
    }
}
// Subdivision decides a width-limited count on the rotated real shift before
// it asks for precision (audit fix 13). Fuzz case 11793351552132886930: 38
// roots in general position, 100-bit accuracy under a 128-bit budget. The
// direct Taylor shift about an off-axis centre re-boxes the error at every
// step and two components could not be counted: the solve stopped for
// precision after 588 work units. With the retry it meets its goal (543).
void rotated_retry_at_the_budget() {
    const std::vector<std::pair<const char *, const char *>> planted{
            {"1/4", "1159/2048"}, {"1/4", "-1/2"}, {"853401/1048576", "-4987/8192"}, {"62/35", "-5445/16384"},
            {"61/512", "1571/4096"}, {"-239/897", "869/4096"}, {"-3/2048", "-1807/2048"},
            {"-93/128", "-34/39"}, {"7/128", "-51187/131072"}, {"-619/596", "-1"}, {"-1/2", "-300/233"},
            {"-3/8", "-4783/16384"}, {"-1/2", "992/531"}, {"13/32", "-5029/65536"},
            {"74095/131072", "-2269/16384"}, {"3/4", "4119/32768"}, {"-79/363", "229/1024"},
            {"1867/4096", "3/4"}, {"237/782", "8043/8192"}, {"321/512", "-315/512"},
            {"1410/983", "377913/524288"}, {"25279/262144", "13/32"}, {"-977/950", "-797/2048"},
            {"-75867/131072", "-231207/262144"}, {"915/668", "45061/131072"}, {"453/16384", "-1519/2048"},
            {"-1/2", "438679/1048576"}, {"-5/16", "-9/5"}, {"1/4", "127/128"}, {"344/495", "12/299"},
            {"149/79", "283/512"}, {"1459/2048", "1225/723"}, {"5/16", "175/97"}, {"-31/32", "-425/1024"},
            {"420409/4194304", "7879/32768"}, {"820/577", "2421609/4194304"}, {"-1/2", "-43809/262144"},
            {"-17/64", "-68/175"}};
    std::vector<Complex> known;
    // The input is the product of (d x - d z) over the roots, as the fuzzer
    // builds it: Gaussian-integer coefficients.
    Rational denominators(1);
    for (const auto &[real, imag] : planted) {
        known.emplace_back(Rational(real), Rational(imag));
        denominators *= Rational(boost::multiprecision::lcm(denominator(known.back().real),
                                                          denominator(known.back().imag)));
    }
    auto polynomial = from_roots(known);
    for (auto &coefficient : polynomial) coefficient = coefficient * Complex(denominators);
    const auto p = input(polynomial);
    certpolysolve::Options o;
    o.method = certpolysolve::Method::classical;
    o.goal = certpolysolve::Goal::accuracy();
    o.bits = 100;
    o.max_bits = 128;
    o.work_budget = 1000;
    o.cell_budget = CPS_DEFAULT_LIVE_CELLS;
    const auto result = certpolysolve::solve(p, o);
    require(result.goal_met() && result.roots().size() == known.size() && total(result) == known.size(),
            "a width-limited count at the precision budget was not decided on the rotated shift");
}
// The retry is for a comparison that arithmetic width obstructed, and costs a
// second Taylor shift. A test left undecided by the roots' positions is not
// retried: 32 roots j/8 + i (j^2 mod 7)/8 and a disk whose circle passes
// within 2^-10 of one of them. At 64 bits the direct comparison is clean
// (not width-limited) and all nine Graeffe steps stay undecided; the last
// step's own flag is set only because width obstructs every comparison at
// that depth. One shift, not two.
void retry_only_when_width_obstructs() {
    d::ComplexPolynomial p;
    acb_poly_one(p.value());
    for (int j = 1; j <= 32; ++j) {
        d::ComplexPolynomial factor;
        d::ComplexBall root;
        acb_set_d_d(root.value(), -j / 8.0, -((j * j) % 7) / 8.0);
        acb_poly_set_coeff_acb(factor.value(), 0, root.value());
        acb_poly_set_coeff_si(factor.value(), 1, 1);
        acb_poly_mul(p.value(), p.value(), factor.value(), 4096); // exact: dyadic roots
    }
    d::ComplexBall centre;
    acb_set_d_d(centre.value(), 0.4375, 0.1875);
    d::Float radius; // the distance to the root 5/8 + 4i/8, less one part in 2^10
    arf_set_d(radius.value(), std::hypot(5 / 8.0 - 0.4375, 4 / 8.0 - 0.1875) * (1 - std::ldexp(1.0, -10)));
    certpolysolve::detail::SearchStatistics statistics{};
    d::CountContext context(statistics);
    context.real_shift = true;
    const auto count = d::pellet_count(p.value(), centre.value(), radius.value(), 64, 9, -1, &context);
    require(count.count < 0 && count.needs_precision && !count.width_limited,
            "the disk test is no longer the geometrically undecided case this check needs");
    require(statistics.taylor_shifts == 1, "a geometrically undecided test was retried on the rotated shift");
}
} // namespace

// With an argument, run only the named check (for mutation checks).
int main(int argc, char **argv) {
    certpolysolve::ThreadScope lifetime;
    try {
        const std::string only = argc > 1 ? argv[1] : "";
        if (only.empty() || only == "steering") steering();
        if (only.empty() || only == "pole") limit_step_pole();
        if (only.empty() || only == "wilkinson") wilkinson_below_its_precision();
        if (only.empty() || only == "cluster") cluster_below_its_precision();
        if (only.empty() || only == "retry") rotated_retry_at_the_budget();
        if (only.empty() || only == "retry-cost") retry_only_when_width_obstructs();
        std::cout << "unsteered models stop for precision; steered ones keep their steps; a width-limited count is "
                     "retried on the rotated shift\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
