#include "exact_winding_api.hpp"
#include "solve/complex_query.hpp"
#include "solve/query_execution.hpp"
#include "solve/local_model_solver.hpp"
#include "methods/secular/secular_refinement.hpp"
#include <array>
#include <iostream>
#include <type_traits>
#include <utility>

namespace d = certpolysolve::detail;
using namespace reference;
static_assert(!std::is_copy_constructible_v<d::QueryExecution>);
static_assert(!std::is_move_constructible_v<d::QueryExecution>);
static_assert(!std::is_move_constructible_v<d::QueryInput>);

static void known_cover(const cps_result &result) {
    require(result.complete && result.goal_met && result.degree == 5 && result.roots.size() == 5,
            "fresh complex query lost a known simple root");
    require(d::disjoint(result.roots), "fresh complex query has overlapping regions");
    constexpr std::array<std::pair<int, int>, 5> known{{{0, 0}, {-2, 0}, {1, 0}, {3, 1}, {3, -1}}};
    for (const auto &[re, im] : known) {
        d::ComplexBall exact;
        acb_set_d_d(exact.value(), re, im);
        size_t containing = 0;
        for (const auto &root : result.roots) {
            require(root.count == 1, "fresh complex query changed a known multiplicity");
            containing += acb_contains(root.box.value(), exact.value()) != 0;
        }
        require(containing == 1, "a planted root is not in exactly one published rectangle");
    }
}

static void queries() {
    const auto object = input(from_roots({Complex(0), Complex(-2), Complex(1), Complex(3, 1), Complex(3, -1)}));
    const auto polynomial = d::dense(*object.get());
    cps_options options;
    options.max_work_bits = 4096;
    const auto automatic = d::solve_complex_query(polynomial, options);
    known_cover(automatic);
    require(automatic.diagnostics.find("heuristic-packets=") != std::string::npos,
            "direct exact coefficient Auto did not install its shared complete allowance");
    for (const auto method : {CPS_METHOD_CLASSICAL, CPS_METHOD_SECULAR}) {
        auto forced = options;
        forced.method = method;
        const auto result = d::solve_complex_query(polynomial, forced);
        known_cover(result);
        require(result.diagnostics.find("heuristic-packets=") == std::string::npos,
                "Auto's complete policy leaked into an explicitly forced method");
    }
    auto staged = options;
    staged.staged_goal = true;
    const auto stage = d::solve_complex_query(polynomial, staged);
    known_cover(stage);
    require(stage.diagnostics.find("heuristic-packets=") == std::string::npos,
            "a wrapper's staged query silently acquired the direct Auto account");
    auto isolation = options;
    isolation.goal = CPS_ROOT_ISOLATION;
    known_cover(d::solve_complex_query(polynomial, isolation));
    for (int invalid = 0; invalid < 2; ++invalid) {
        auto query = options;
        if (invalid == 0) query.real_only = true;
        if (invalid == 1) query.region.real_lower = d::Number("0");
        bool refused = false;
        try { (void)d::solve_complex_query(polynomial, query); }
        catch (const std::invalid_argument &) { refused = true; }
        require(refused, "complex entry admitted a query requiring real routing");
    }

    // A nonmonotonic caller callback must stay cancelled through a fresh child
    // query, and a later independent query must not retain that borrowed state.
    int polls = 0;
    auto cancelled = options;
    cancelled.cancel_data = &polls;
    cancelled.cancel = [](void *p) -> int { return ++*static_cast<int *>(p) == 1; };
    {
        d::QueryInput source(polynomial);
        d::QueryExecution outer(source, &cancelled);
        require(!outer.refusal && outer.next.cancelled(), "outer query did not latch cancellation");
        const auto result = d::solve_complex_query(polynomial, outer.next);
        require(result.termination == CPS_CANCELLED && !result.goal_met, "fresh child resumed after cancellation");
        require(outer.next.cancelled() && polls == 1, "nested query repolled an already latched callback");
        auto local_options = d::local_model_options(outer.next, 80, {}, {});
        require(local_options.has_value(), "the cancellation fixture spent its work or sweep allowance");
        const auto local = d::solve_local_model(d::DensePolynomial(polynomial), *local_options);
        require(local.termination == CPS_CANCELLED && !local.goal_met && polls == 1,
                "local-model adapter lost the caller's latched cancellation");
    }
    known_cover(d::solve_complex_query(polynomial, options));

    d::QueryInput source(polynomial);
    d::ComplexQueryState state(source);
    d::QueryExecution query(source, &options);
    known_cover(query.run([&] { return state.advance(query.next, query.degree); }));
    bool refused = false;
    try { (void)state.advance(query.next, query.degree); }
    catch (const std::logic_error &) { refused = true; }
    require(refused, "ordinary engine state allowed a second advance on consumed input");
}

int main() {
    certpolysolve::ThreadScope lifetime;
    try { queries(); }
    catch (const std::exception &error) { std::cerr << error.what() << '\n'; return 1; }
}
