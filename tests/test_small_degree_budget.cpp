// Failed closed-form attempts are work too. Exercise the public Auto entry and
// check its retained/handed-off covers with independent exact Sturm arithmetic.
#include "exact_winding_api.hpp"
#include <iostream>

using namespace reference;

static void verify(const Polynomial &polynomial, const certpolysolve::Result &result, uint64_t budget) {
    require(result.get() && result.stop() != certpolysolve::Stop::error, "small-degree solve returned an API error");
    require(result.statistics().work_units <= budget, "closed form and fallback exceeded their shared allowance");
    const auto roots = result.roots();
    check_disjoint(roots);
    uint64_t total = 0;
    // Every fixture is all-real. Thus exact Sturm counts on each rectangle's
    // real projection, with zero in its imaginary projection, count its roots.
    for (const auto &root : roots) {
        require(endpoint(root.box.imag.lo) <= 0 && endpoint(root.box.imag.hi) >= 0,
                "returned rectangle misses the real axis");
        require(roots_closed(polynomial, endpoint(root.box.real.lo), endpoint(root.box.real.hi)) ==
                    static_cast<int>(root.count), "independent count disagrees with a returned rectangle");
        if (root.multiplicity)
            require(root.multiplicity == 1 && root.count == 1, "square-free fixture gained a multiple root");
        total += root.count;
    }
    if (result.complete()) require(total == polynomial.size() - 1, "complete cover lost roots");
    if (result.goal_met()) require(result.complete(), "successful goal lacks a complete cover");
}

static void capped_attempt(bool accuracy) {
    const Polynomial polynomial{-2, 0, 1};
    const auto source = input(polynomial);
    certpolysolve::Options options;
    options.bits = 200;
    options.goal = accuracy ? certpolysolve::Goal::accuracy() : certpolysolve::Goal::approximation();
    for (const uint64_t cap : {64, 128}) {
        options.max_bits = cap;
        options.work_budget = 1;
        const auto stopped = certpolysolve::solve(source, options);
        verify(polynomial, stopped, 1);
        require(!stopped.goal_met() && stopped.stop() == certpolysolve::Stop::work_budget,
                "spent closed-form allowance did not stop as work-limited");
        require(stopped.statistics().work_units == 1 && stopped.engine_counters().precision_attempts == 1 &&
                    stopped.engine_counters().isolation_passes == 0 && stopped.statistics().exact_factors == 0,
                "spent closed-form allowance launched a general solve");
        require(stopped.complete() && stopped.roots().size() == 2,
                "work-limited closed form discarded its proved disjoint cover");

        // At this cap the direct attempt spends one unit. The general route
        // must see exactly the seven remaining units, and both engines count.
        options.work_budget = 8;
        const auto automatic = certpolysolve::solve(source, options);
        auto forced_options = options;
        forced_options.method = certpolysolve::Method::secular;
        forced_options.work_budget = 7;
        const auto forced = certpolysolve::solve(source, forced_options);
        verify(polynomial, automatic, 8);
        verify(polynomial, forced, 7);
        require(automatic.statistics().work_units == forced.statistics().work_units + 1 &&
                    automatic.engine_counters().precision_attempts == forced.engine_counters().precision_attempts + 1,
                "fallback omitted the closed-form attempt's work or precision attempt");
        require(automatic.statistics().highest_work_bits >= cap, "fallback lost the earlier working precision");
        require(!automatic.goal_met() && automatic.stop() == certpolysolve::Stop::precision_cap,
                "insufficient precision was reported as success or an unrelated stop");

        options.work_budget = 2;
        const auto short_fallback = certpolysolve::solve(source, options);
        verify(polynomial, short_fallback, 2);
        require(short_fallback.statistics().work_units == 2 &&
                    short_fallback.engine_counters().precision_attempts == 2 &&
                    short_fallback.stop() == certpolysolve::Stop::work_budget,
                "fallback did not honor the one remaining work unit");
    }
    // Enough precision: the original successful direct path still costs one.
    options.max_bits = 256;
    options.work_budget = 1;
    const auto success = certpolysolve::solve(source, options);
    verify(polynomial, success, 1);
    require(success.goal_met() && success.statistics().work_units == 1 &&
                success.engine_counters().precision_attempts == 1, "successful closed form changed its work");
    if (accuracy) for (const auto &root : success.roots()) check_accuracy(root, 200);
}

static void repeated_attempts() {
    // Distinct roots 1 and 1 + 2^-200. The discriminant needs more precision
    // than the first three attempts, so each exhausted allowance is observable.
    const Polynomial polynomial = multiply(Polynomial{-1, 1}, Polynomial{-1 - power2(-200), 1});
    const auto source = input(polynomial);
    certpolysolve::Options options;
    options.bits = 53;
    options.max_bits = 1024;
    options.goal = certpolysolve::Goal::accuracy();
    for (const uint64_t budget : {1, 2, 3, 8}) {
        options.work_budget = budget;
        const auto result = certpolysolve::solve(source, options);
        verify(polynomial, result, budget);
        const auto statistics = result.statistics();
        const auto counters = result.engine_counters();
        require(counters.isolation_passes == 0 && statistics.exact_factors == 0,
                "closed-form work exhaustion launched another engine");
        require(counters.precision_attempts == statistics.work_units &&
                    statistics.precision_increases + 1 == statistics.work_units,
                "closed-form precision attempts were not all charged");
        if (budget <= 3) {
            require(!result.goal_met() && result.stop() == certpolysolve::Stop::work_budget &&
                        statistics.work_units == budget, "multiple direct attempts exceeded their allowance");
        } else {
            require(result.goal_met(), "ample-budget close roots no longer finish directly");
            for (const auto &root : result.roots()) check_accuracy(root, 53);
        }
    }
    // Two direct attempts reach the cap; one remaining general-route unit
    // starts at a lower precision. The high-water mark must include both.
    options.max_bits = 128;
    options.work_budget = 3;
    const auto handed_off = certpolysolve::solve(source, options);
    verify(polynomial, handed_off, 3);
    require(handed_off.stop() == certpolysolve::Stop::work_budget && handed_off.statistics().work_units == 3 &&
                handed_off.engine_counters().precision_attempts == 3,
            "multi-attempt fallback did not retain all spent work");
    require(handed_off.statistics().highest_work_bits == 128 && handed_off.statistics().final_work_bits < 128,
            "fallback lost the earlier, higher working precision");
}

static void stopped_cubic_quartic_covers() {
    // Both real and complex formulas can prove a full cover before they reach
    // 200-bit accuracy. Audit the published rectangles, including the nonreal
    // roots, with exact Sturm--Cauchy counts independent of the formulas.
    const std::vector<ComplexPolynomial> polynomials{
        {{1}, {-3}, {0}, {1}}, {{6}, {0}, {-5}, {0}, {1}},
        {{-2}, {0}, {0}, {1}}, {{1}, {1}, {0}, {0}, {1}},
        from_roots({{Rational(1) / 3, Rational(1) / 7}, {Rational(2) / 5, -Rational(1) / 11},
                    {-Rational(7) / 3, Rational(1) / 13}}),
        from_roots({{Rational(1) / 3, Rational(1) / 7}, {Rational(2) / 5, -Rational(1) / 11},
                    {-Rational(7) / 3, Rational(1) / 13}, {Rational(1) / 17, -Rational(2) / 19}})};
    for (const auto &polynomial : polynomials) {
        certpolysolve::Options options;
        options.bits = 200;
        options.max_bits = 64;
        options.work_budget = 1;
        options.goal = certpolysolve::Goal::accuracy();
        const auto result = certpolysolve::solve(input(polynomial), options);
        require(result.get() && result.complete() && !result.goal_met() &&
                    result.stop() == certpolysolve::Stop::work_budget && result.statistics().work_units == 1,
                "cubic/quartic stop lost a cover or exceeded its allowance");
        require(result.engine_counters().isolation_passes == 0 && result.statistics().exact_factors == 0,
                "cubic/quartic budget exhaustion launched the general route");
        const auto roots = result.roots();
        check_disjoint(roots);
        uint64_t total = 0;
        for (const auto &root : roots) {
            const auto count = rectangle_count(polynomial, root);
            require(count && *count == static_cast<int>(root.count),
                    "exact argument count disagrees with a stopped formula's cover");
            require(root.count == 1 && root.multiplicity == 1, "square-free formula lost individual roots");
            total += root.count;
        }
        require(total == polynomial.size() - 1, "stopped formula cover lost roots");
    }
}

int main() {
    try {
        capped_attempt(false);
        capped_attempt(true);
        repeated_attempts();
        stopped_cubic_quartic_covers();
        std::cout << "Closed-form work budgets and independent cover counts passed\n";
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
