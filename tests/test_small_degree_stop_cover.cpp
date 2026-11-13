// A complete algebraic enumeration can retain a counted cluster when its
// final work unit cannot separate individual branch boxes. Exact reference
// counts check the entire published rectangles, including zeros/multiplicity.
#include "exact_winding_api.hpp"
#include <certpolysolve/certpolysolve_testing.h>
#include <iostream>
#include <limits>
#include <memory>

using namespace reference;

struct Fixture {
    const char *name;
    std::vector<Complex> points;
    bool approximation_met;
};

static void retained_cover(const Fixture &fixture, bool accuracy) {
    auto polynomial = from_roots(fixture.points);
    if (std::string(fixture.name) == "quartic-complex")
        for (auto &coefficient : polynomial) coefficient = coefficient * Complex(Rational(5) / 7, Rational(8) / 11);
    certpolysolve::Options options;
    options.bits = 53;
    options.max_bits = 1024;
    options.work_budget = 1;
    options.goal = accuracy ? certpolysolve::Goal::accuracy() : certpolysolve::Goal::approximation();
    const auto result = certpolysolve::solve(input(polynomial), options);
    require(result.get() && result.complete(), "final algebraic attempt discarded a complete counted cover");
    const bool expected_goal = !accuracy && fixture.approximation_met;
    require(result.goal_met() == expected_goal, "counted cluster attained the wrong goal");
    require(result.stop() == (expected_goal ? certpolysolve::Stop::none : certpolysolve::Stop::work_budget),
            "terminal counted cover has the wrong stop reason");
    require(result.statistics().work_units == 1 && result.engine_counters().precision_attempts == 1 &&
                result.engine_counters().isolation_passes == 0 && result.statistics().exact_factors == 0,
            "cover recovery spent another attempt or launched fallback");
    const auto roots = result.roots();
    check_disjoint(roots);
    uint64_t total = 0;
    bool merged = false;
    for (const auto &root : roots) {
        const auto exact = rectangle_count(polynomial, root);
        require(exact && *exact == static_cast<int>(root.count), "independent whole-box count disagrees");
        total += root.count;
        if (root.certificate.find("branch-cover union") != std::string::npos) {
            merged = true;
            require(root.count > 1 && !root.multiplicity, "merged branches retained an individual multiplicity");
            require(!accuracy || !root.value, "unresolved accuracy cluster exposes an approximate point");
        }
        // Every root is planted, so realness labels can be checked without
        // trusting the production formulas or their reflection argument.
        const auto box = rectangle(root);
        for (const auto &point : fixture.points) {
            if (point.real < box.real_lower || point.real > box.real_upper || point.imag < box.imag_lower || point.imag > box.imag_upper) continue;
            if (root.realness == certpolysolve::Realness::real) require(point.imag == 0, "nonreal root labelled real");
            if (root.realness == certpolysolve::Realness::not_real) require(point.imag != 0, "real root labelled nonreal");
        }
        if (expected_goal) check_accuracy(root, 53);
    }
    require(merged && total == polynomial.size() - 1, "terminal union lost its complete degree ledger");
}

static void cancellation_precedes_recovery() {
    const auto source = input(from_roots({{1}, {1 + power2(-200)}}));
    cps_options *created = nullptr;
    require(cps_options_create(&created) == CPS_OK, "create cancellation options");
    std::unique_ptr<cps_options, decltype(&cps_options_destroy)> options(created, cps_options_destroy);
    require(cps_options_set_work_budget(options.get(), 1) == CPS_OK &&
                cps_options_set_max_work_bits(options.get(), 1024) == CPS_OK &&
                cps_options_set_goal(options.get(), CPS_ROOT_ACCURACY) == CPS_OK,
            "configure cancellation options");
    const auto solve = [&](uint64_t checkpoint) {
        require(cps_testing_options_stop_at_poll(options.get(), checkpoint) == CPS_OK, "set cancellation checkpoint");
        cps_result *raw = nullptr;
        const auto code = cps_solve(source.get(), options.get(), &raw);
        std::unique_ptr<cps_result, decltype(&cps_result_destroy)> result(raw, cps_result_destroy);
        require(code == CPS_OK && result, "cancellation solve failed");
        return result;
    };
    const auto counted = solve(std::numeric_limits<uint64_t>::max());
    const auto polls = cps_testing_result_polls(counted.get());
    require(polls > 0 && cps_result_complete(counted.get()), "uncancelled control missed its counted cover");
    for (uint64_t checkpoint = 1; checkpoint <= polls; ++checkpoint) {
        const auto cancelled = solve(checkpoint);
        require(cps_result_termination(cancelled.get()) == CPS_CANCELLED &&
                    cps_result_stop_reason(cancelled.get()) == CPS_STOP_CANCELLED &&
                    !cps_result_goal_met(cancelled.get()) && !cps_result_complete(cancelled.get()) &&
                    cps_result_size(cancelled.get()) == 0,
                "cover recovery ran after cancellation");
        require(cps_testing_result_polls(cancelled.get()) == checkpoint, "cancellation was polled again");
    }
}

// The real route consumes the algebraic cover through its nested complex
// query. At a work stop it may keep a counted interval crossing the requested
// region's end; the count is for the whole interval, and the answer stays
// incomplete. This is the documented partial-result contract.
static void real_partial_regions() {
    const Rational a = power2(15), scale = power2(-45);
    const ComplexPolynomial polynomial{{-2 * scale}, {4 * a * scale}, {-2 * a * a * scale}, {scale}};
    const auto p = input(polynomial);
    certpolysolve::Options options;
    options.bits = 100;
    options.max_bits = 64;
    options.work_budget = 1;
    options.goal = certpolysolve::Goal::isolation();
    const auto check = [&](const certpolysolve::Result &result, uint64_t expected) {
        require(result.get() && !result.complete() && !result.goal_met() &&
                    result.stop() == certpolysolve::Stop::work_budget,
                "partial real cover has the wrong termination");
        require(result.statistics().work_units == 1 && result.engine_counters().isolation_passes == 0,
                "partial real cover spent fallback work");
        check_disjoint(result.roots());
        uint64_t count = 0;
        for (const auto &root : result.roots()) {
            require(root.realness == certpolysolve::Realness::real && endpoint(root.box.imag.lo) == 0 &&
                        endpoint(root.box.imag.hi) == 0,
                    "partial real cover contains a nonreal interval");
            const auto exact = rectangle_count(polynomial, root);
            require(exact && *exact == static_cast<int>(root.count), "partial real interval count is wrong");
            count += root.count;
        }
        require(count == expected, "partial real count ledger is wrong");
    };
    check(certpolysolve::real_roots(p, options), 3);
    const Rational lower = Rational(1) / 32768, upper = Rational(1) / 16384;
    const auto requested_count = rectangle_count(polynomial, lower, upper, Rational(0), Rational(0));
    require(requested_count && *requested_count == 1, "independent requested-interval count");
    const auto crossing = certpolysolve::real_roots(
        p, certpolysolve::between(certpolysolve::Number(lower.str()), certpolysolve::Number(upper.str())), options);
    check(crossing, 2);
    const auto roots = crossing.roots();
    require(roots.size() == 1 && !roots[0].multiplicity && endpoint(roots[0].box.real.lo) < lower &&
                endpoint(roots[0].box.real.hi) > lower,
            "stopped cluster was collapsed at the query endpoint");
    // The distant individual root is retained when the close pair is outside.
    options.goal = certpolysolve::Goal::accuracy();
    check(certpolysolve::real_roots(p, certpolysolve::at_least(certpolysolve::Number(9)), options), 1);
}

int main() {
    try {
        real_partial_regions();
        const auto d = power2(-200);
        const std::vector<Fixture> fixtures{
            {"quadratic-real", {{1}, {1 + d}}, false},
            {"quadratic-complex", {{1, 1}, {1 + d, 1}}, false},
            {"cubic-complex", {{-power2(-139), Rational(1) / 2}, {power2(-139), Rational(1) / 2},
                               {0, Rational(107) / 1024}}, true},
            {"cubic-zero", {{0}, {1}, {1 + d}}, false},
            {"cubic-repeated", {{1}, {1}, {1 + d}}, true},
            {"quartic-real", {{-2}, {3}, {1}, {1 + d}}, true},
            {"quartic-zero", {{0}, {-2}, {1}, {1 + d}}, true},
            {"quartic-repeated", {{1}, {1}, {1 + d}, {1 + d}}, true},
            {"quartic-complex", {{1, power2(-103)}, {1, -power2(-103)},
                                 {Rational(39) / 64, power2(-99)}, {Rational(39) / 64, -power2(-99)}}, false}};
        for (const auto &fixture : fixtures)
            for (const bool accuracy : {true, false}) retained_cover(fixture, accuracy);
        cancellation_precedes_recovery();
        std::cout << "Terminal algebraic covers: twenty-one exact-count queries and every cancellation checkpoint passed\n";
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
