#include "support/cost-progress/prepared_isolation.hpp"
#include "exact_winding_api.hpp"
#include "numeric/arithmetic/polynomial_arithmetic.hpp"
#include "platform/deadline.hpp"
#include "solve/policies.hpp"
#include "solve/solve_state.hpp"
#include <filesystem>
#include <iostream>

namespace d = certpolysolve::detail;
namespace ir = isolation_research;
using reference::require;

namespace {
constexpr uint64_t turn_limit = 30000; // test deadline, not an algorithmic threshold

reference::Rational endpoint(const d::Endpoint &point) {
    const auto text = point.text();
    return reference::dyadic(text.mantissa, text.exponent);
}
void check(const std::vector<d::Root> &cover, const reference::ComplexPolynomial &input, int bits) {
    std::vector<reference::Rectangle> rectangles;
    uint64_t total = 0;
    for (const auto &root : cover) {
        reference::Rectangle box{endpoint(root.endpoints[0]), endpoint(root.endpoints[1]),
                                 endpoint(root.endpoints[2]), endpoint(root.endpoints[3]), root.count};
        const auto count = reference::rectangle_count(input, box.real_lower, box.real_upper,
                                                       box.imag_lower, box.imag_upper);
        require(count && *count > 0 && static_cast<uint64_t>(*count) == root.count,
                "independent exact rectangle count differs");
        require(root.count == root.multiplicity, "refinement must retain individual-root multiplicities");
        if (bits) reference::check_accuracy(box, bits);
        total += root.count;
        rectangles.push_back(std::move(box));
    }
    require(total == input.size() - 1, "refinement lost root coverage");
    reference::check_disjoint(rectangles);
}
reference::ComplexPolynomial exact_input(const std::filesystem::path &path) {
    reference::ComplexPolynomial polynomial;
    for (const auto &row : scaling_io::read(path.string()))
        polynomial.emplace_back(reference::Rational(reference::Integer(row[0])),
                                reference::Rational(reference::Integer(row[1])));
    return polynomial;
}
void finish(ir::Run &run) {
    while (!run.goal_met && run.turns < turn_limit)
        require(run.tick() == CPS_STOP_NONE, "prepared continuation unexpectedly stopped");
    require(run.goal_met, "prepared continuation exhausted its test allowance");
}
void continuation(const std::filesystem::path &path) {
    const auto input = exact_input(path);
    auto options = ir::isolation_options(4096);
    options.deadline = d::deadline_after(30, std::chrono::steady_clock::now());
    ir::Run run(ir::IntegerFactor(path.string()), 64, options, d::SubdivisionPrecision::TrialLocal);
    finish(run);
    check(run.cover, input, 0);
    auto *const engine = &run.engine();
    const auto isolation_turns = run.turns;
    for (const unsigned bits : {53u, 160u, 512u}) {
        const auto deliveries = run.enclosure_calls, turns = run.turns;
        run.refine(bits);
        require(&run.engine() == engine && run.enclosure_calls == deliveries && run.turns == turns,
                "refinement restarted the prepared solve");
        finish(run);
        check(run.cover, input, static_cast<int>(bits));
        std::cout << path.stem().string() << " bits=" << bits << " isolation_turns=" << isolation_turns
                  << " additional_turns=" << run.turns - turns << " peak=" << run.peak_enclosure_bits << '\n';
    }
}
void cap_then_resume(const std::filesystem::path &path) {
    const auto input = exact_input(path);
    auto options = ir::isolation_options(64);
    options.deadline = d::deadline_after(30, std::chrono::steady_clock::now());
    ir::Run run(ir::IntegerFactor(path.string()), 64, options, d::SubdivisionPrecision::TrialLocal);
    finish(run);
    run.refine(160);
    cps_stop_reason reason = CPS_STOP_NONE;
    while (!run.goal_met && reason == CPS_STOP_NONE && run.turns < turn_limit) reason = run.tick();
    require(reason == CPS_STOP_PRECISION && !run.goal_met, "refinement must expose an inadequate cap");
    check(run.cover, input, 0);
    run.options.max_work_bits = 1024;
    run.engine().reopen(run.options);
    finish(run);
    check(run.cover, input, 160);
}
void whole_input() {
    // x (x²-2) (x-1)²: zero extraction, a repeated factor, and irrational
    // roots. The oracle counts against this original polynomial, not factors.
    const reference::ComplexPolynomial polynomial{{0}, {-2}, {4}, {-1}, {-2}, {1}};
    auto input = reference::input(polynomial);
    cps_options options;
    options.goal = CPS_ROOT_ACCURACY;
    options.bits = 160;
    options.max_work_bits = 4096;
    options.work_units = turn_limit;
    options.deadline = d::deadline_after(30, std::chrono::steady_clock::now());
    d::SolveState state(d::dense(*input.get()), d::standalone_subdivision_policy());
    const auto result = state.advance(&options);
    require(result.termination == CPS_SUCCESS && result.goal_met && result.complete,
            "whole-input standalone refinement did not complete");
    require(!result.statistics_v3.candidate_sweeps && !result.statistics_v3.secular_regenerations,
            "standalone refinement ran a numerical candidate engine");
    check(result.roots, polynomial, 160);
}
} // namespace

int main(int argc, char **argv) {
    certpolysolve::ThreadScope lifetime;
    try {
        require(argc == 2, "expected refinement-fixture directory");
        const std::filesystem::path directory(argv[1]);
        for (const auto *name : {"sqrt2", "sqrt_i", "translated_sqrt2", "near_zero"})
            continuation(directory / (std::string(name) + ".ps"));
        cap_then_resume(directory / "sqrt2.ps");
        whole_input();
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
