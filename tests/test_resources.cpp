#include "numeric/arithmetic/polynomial_arithmetic.hpp"
#include <certpolysolve/certpolysolve.hpp>
#include "api_test_support.hpp"
#include "numeric/exact/exact_preprocessing.hpp"
#include "numeric/candidates/machine_candidates.hpp"
#include "methods/piecewise/piecewise_roots.hpp"
#include "numeric/exact/centering.hpp"
#include "methods/direct/small_degree.hpp"
#include "exact_winding_api.hpp"
#include "platform/deadline.hpp"
#include <chrono>
#include <cmath>
#include <limits>
#include <iostream>
#include <stdexcept>

namespace d = certpolysolve::detail;
static void require(bool ok, const char *why) { if (!ok) throw std::runtime_error(why); }

static void deadline_boundaries() {
    using Clock = std::chrono::steady_clock;
    using Duration = Clock::duration;
    const auto second = std::chrono::duration_cast<Duration>(std::chrono::seconds(1));
    const Clock::time_point now{10 * second};
    const auto last = Clock::time_point::max();
    require(d::deadline_after(0, now) == now, "zero delay changed the clock");
    require(d::deadline_after(60, now) == now + 60 * second, "ordinary deadline changed");
    require(d::deadline_after(2, Clock::time_point{-second}) == Clock::time_point{second},
            "deadline assumed a nonnegative clock epoch");
    const double tick = std::chrono::duration<double>(Duration{1}).count();
    require(d::deadline_after(tick / 2, now) == now, "sub-tick limit did not truncate");
    require(d::deadline_after(5, last - 5 * second) == last, "exact final deadline was lost");
    require(d::deadline_after(10, last - 5 * second) == last, "deadline addition overflowed");
    require(d::deadline_after(0.5, last - second) == last - second / 2,
            "representable deadline near the end was clamped");
    const double maximum_seconds = std::chrono::duration<double>(Duration::max()).count();
    for (double seconds : {maximum_seconds, std::nextafter(maximum_seconds, INFINITY),
                           1e100, std::numeric_limits<double>::max()})
        require(d::deadline_after(seconds, now) == last, "large finite deadline wrapped or expired");
}

static void closed_form_effort() {
    // (3x - 1)^n has one known real root of multiplicity n. Exact algebra
    // must honor the same effort policy through Auto's closed-form entry.
    reference::Polynomial polynomial{1};
    for (uint64_t degree = 1; degree <= 4; ++degree) {
        polynomial = reference::multiply(polynomial, reference::Polynomial{-1, 3});
        const auto source = reference::input(polynomial);
        cps_options direct_options;
        direct_options.goal = CPS_ROOT_APPROXIMATION;
        direct_options.max_work_bits = 256;
        direct_options.work_units = 128;
        for (const uint64_t allowance : {UINT64_C(0), UINT64_C(1), UINT64_C(100000)}) {
            direct_options.exact_operations = allowance;
            const auto direct = d::direct_small_degree(d::dense(*source.get()), direct_options);
            if (degree > 1 && allowance <= 1)
                require(!direct, "closed-form exact preparation ignored its effort allowance");
            else
                require(direct && direct->complete && direct->goal_met,
                        "linear or amply admitted closed form was refused");

            certpolysolve::Options options;
            options.goal = certpolysolve::Goal::approximation();
            options.max_bits = 256;
            options.work_budget = 128;
            certpolysolve::ResourceLimits limits;
            limits.exact_operations = allowance;
            options.resources = limits;
            const auto result = certpolysolve::solve(source, options);
            require(result.get() && result.stop() != certpolysolve::Stop::error, "effort-limited Auto query failed");
            if (degree > 1 && allowance <= 1)
                require(result.preparation_limited(), "Auto hid its exact-preparation refusal");
            const auto roots = result.roots();
            reference::check_disjoint(roots);
            uint64_t total = 0;
            for (const auto &root : roots) {
                require(reference::endpoint(root.box.imag.lo) <= 0 && reference::endpoint(root.box.imag.hi) >= 0 &&
                            reference::roots_closed(polynomial, reference::endpoint(root.box.real.lo),
                                                    reference::endpoint(root.box.real.hi)) == int(root.count),
                        "effort-limited result disagrees with independent exact root count");
                total += root.count;
            }
            if (result.complete()) require(total == degree, "effort-limited complete cover lost roots");
        }
        if (degree == 1) continue;
        // The operation estimate must grow with exact operand size, even
        // though scaling by 2^4096 changes none of these roots.
        auto scaled = polynomial;
        for (auto &coefficient : scaled) coefficient *= reference::power2(4096);
        const auto scaled_source = reference::input(scaled);
        direct_options.exact_operations = 64;
        require(!d::direct_small_degree(d::dense(*scaled_source.get()), direct_options),
                "closed-form exact effort ignored multiword operands");
        direct_options.exact_operations = UINT64_MAX;
        const auto unlimited = d::direct_small_degree(d::dense(*scaled_source.get()), direct_options);
        require(unlimited && unlimited->goal_met,
                "unlimited exact effort refused a scaled small-degree input");
    }
}

namespace {
// Observe GMP without changing its allocator or refusing any request. The
// large-denominator regression checks that refusal precedes growing the LCM.
size_t largest_request = 0;
void *(*original_allocate)(size_t);
void *(*original_reallocate)(void *, size_t, size_t);
void (*original_release)(void *, size_t);
void *record_allocate(size_t n) {
    largest_request = std::max(largest_request, n);
    return original_allocate(n);
}
void *record_reallocate(void *p, size_t old_n, size_t n) {
    largest_request = std::max(largest_request, n);
    return original_reallocate(p, old_n, n);
}
struct RecordAllocations {
    RecordAllocations() {
        largest_request = 0;
        mp_get_memory_functions(&original_allocate, &original_reallocate, &original_release);
        mp_set_memory_functions(record_allocate, record_reallocate, original_release);
    }
    ~RecordAllocations() { mp_set_memory_functions(original_allocate, original_reallocate, original_release); }
};
}

static void closed_form_storage() {
    const char *integers[] = {"9223372036854775783", "9223372036854775805", "9223372036854775643",
                              "9223372036854775549", "9223372036854775507"};
    for (size_t degree = 2; degree <= 4; ++degree) {
        reference::ComplexPolynomial polynomial;
        for (size_t i = 0; i <= degree; ++i) {
            const reference::Rational value(integers[i]);
            polynomial.emplace_back(degree == 2 ? value : 1 / value);
        }
        const auto source = reference::input(polynomial);
        cps_options direct_options;
        direct_options.goal = CPS_ROOT_APPROXIMATION;
        direct_options.max_work_bits = 256;
        direct_options.work_units = 128;
        direct_options.exact_budget = 70;
        require(!d::direct_small_degree(d::dense(*source.get()), direct_options),
                "closed form ignored exact operand storage");
        direct_options.exact_budget = UINT64_MAX;
        direct_options.temporary_bytes = 9; // 72 exact bits, still below the products/LCM.
        require(!d::direct_small_degree(d::dense(*source.get()), direct_options),
                "closed form ignored temporary storage for exact operands");
        direct_options.temporary_bytes = UINT64_MAX;
        direct_options.exact_budget = 4096;
        const auto direct = d::direct_small_degree(d::dense(*source.get()), direct_options);
        require(direct && direct->goal_met, "ample exact storage refused a small closed form");

        for (const uint64_t allowance : {UINT64_C(70), UINT64_C(4096)}) {
            certpolysolve::Options options;
            options.max_bits = 256;
            options.work_budget = 128;
            options.goal = certpolysolve::Goal::approximation();
            options.resources.emplace();
            options.resources->exact_materialization_bits = allowance;
            const auto result = certpolysolve::solve(source, options);
            require(result.goal_met() && result.complete(), "storage refusal lost the numerical fallback");
            const auto roots = result.roots();
            reference::check_disjoint(roots);
            uint64_t count = 0;
            for (const auto &root : roots) {
                const auto independent = reference::rectangle_count(polynomial, root);
                require(independent && *independent == int(root.count),
                        "storage-limited small-degree answer failed independent exact count");
                reference::check_accuracy(root, 53);
                count += root.count;
            }
            require(count == degree, "storage-limited complete cover lost roots");
        }
    }
}

static void closed_form_denominator_storage() {
    constexpr size_t coefficient_bits = 8192;
    reference::Polynomial polynomial;
    for (const int offset : {17, 19, 23, 29, 31})
        polynomial.push_back(1 / (reference::power2(coefficient_bits) - offset));
    for (const size_t degree : {3, 4}) {
        const auto source = reference::input(reference::Polynomial(polynomial.begin(), polynomial.begin() + degree + 1));
        cps_options options;
        options.goal = CPS_ROOT_APPROXIMATION;
        options.exact_budget = coefficient_bits + 64;
        options.max_work_bits = 256;
        options.work_units = 128;
        RecordAllocations recording;
        require(!d::direct_small_degree(d::dense(*source.get()), options),
                "oversized exact denominator was admitted");
        // Enough room for two input-sized operands plus four limb-rounding
        // words, but not the three/four-input LCM the old route constructed.
        require(largest_request <= 2 * coefficient_bits / 8 + 4 * sizeof(mp_limb_t),
                "closed form constructed its oversized denominator before refusal");
    }
}

int main(int argc, char **argv) {
    try {
        deadline_boundaries();
        if (argc == 2) {
            if (std::string(argv[1]) == "storage") closed_form_storage();
            else if (std::string(argv[1]) == "denominator-storage") closed_form_denominator_storage();
            else return 2;
            return 0;
        }
        closed_form_storage();
        closed_form_denominator_storage();
        closed_form_effort();
        {
            auto options = certpolysolve_test::raw_options();
            cps_options *raw = options.get();
            certpolysolve::ResourceLimits limits;
            require(cps_options_get_resource_limits(raw, &limits) == CPS_OK, "resource limits unreadable");
            require(limits.exact_materialization_bits == UINT64_MAX && limits.exact_operations == UINT64_MAX &&
                    limits.temporary_bytes == UINT64_MAX && limits.flags == 0, "new resource defaults");
            require(cps_options_set_precision(raw, 1000001) == CPS_OK && cps_options_set_max_work_bits(raw, 2000002) == CPS_OK,
                    "large precision refused");
            require(cps_options_set_max_iterations(raw, 2000002) == CPS_OK, "large iteration budget refused");
            require(cps_options_set_work_budget(raw, 0) == CPS_OK && cps_options_set_live_cells(raw, 10000001) == CPS_OK,
                    "large search budget refused");
            require(certpolysolve_test::set_exact_allowance(raw, UINT64_C(2000000000)) == CPS_OK, "exact allowance refused");
            require(cps_options_get_resource_limits(raw, &limits) == CPS_OK, "resource limits unreadable");
            require(limits.exact_materialization_bits == UINT64_C(2000000000), "a large exact allowance was capped");
            require(certpolysolve_test::set_exact_allowance(raw, 0) == CPS_OK, "zero exact allowance refused");
            require(cps_options_get_resource_limits(raw, &limits) == CPS_OK, "resource limits unreadable");
            require(limits.exact_materialization_bits == 0, "zero must disable exact preparation");
            require(cps_options_set_precision(raw, UINT64_MAX) == CPS_OUT_OF_RANGE &&
                        cps_options_set_max_work_bits(raw, UINT64_MAX) == CPS_OUT_OF_RANGE,
                    "precision overflow accepted");
            limits.flags = UINT64_C(2);
            require(cps_options_set_resource_limits(raw, &limits) == CPS_INVALID_ARGUMENT, "unknown flags accepted");
        }
        {
            certpolysolve::Polynomial polynomial;
            polynomial.set(1000001, 1.0);
            require(polynomial.ok(), "a large exponent refused");
            require(cps_polynomial_set_coefficient_d(const_cast<cps_polynomial *>(polynomial.get()),
                    UINT64_MAX, 1, 0) == CPS_OUT_OF_RANGE, "index overflow accepted");
            d::Number long_literal(std::string(1000001, '0'));
            require(long_literal.zero(), "old literal-length cap remains");
        }
        {
            certpolysolve::detail::DensePolynomial compact;
            compact.coefficients.resize(2);
            compact.coefficients[0].real = d::Number("1@1000000000000");
            compact.coefficients[1].real = d::Number("1");
            d::IntegerPolynomial exact;
            d::PreparationDiagnostic diagnostic;
            cps_options options;
            require(!d::streamed_real_polynomial(exact, compact.coefficients, options, &diagnostic), "huge exponent expanded implicitly");
            require(std::string(diagnostic.reason) == "compact-exponent-expansion-needs-explicit-allowance", "compact refusal unexplained");
            compact.coefficients[0].real = d::Number("1@4096");
            options.allow_exponent_expansion = true;
            options.exact_budget = 20000;
            require(d::streamed_real_polynomial(exact, compact.coefficients, options), "explicit exponent allowance ignored");
            options.temporary_bytes = 1;
            require(!d::streamed_real_polynomial(exact, compact.coefficients, options), "exact temporary allowance ignored");
            options.temporary_bytes = UINT64_MAX;
            options.exact_operations = 0;
            require(!d::streamed_real_polynomial(exact, compact.coefficients, options), "exact effort zero ignored");
            // A linear nonzero part is square-free by degree, which needs no
            // exact effort. A quadratic needs the modular proof or exact
            // algebra, and both must respect the zero allowance.
            require(d::prepare_squarefree_input(compact, options).available(), "square-free by degree refused");
            certpolysolve::detail::DensePolynomial quadratic;
            quadratic.coefficients.resize(3);
            quadratic.coefficients[0].real = d::Number("1@4096");
            quadratic.coefficients[2].real = d::Number("1");
            require(!d::prepare_squarefree_input(quadratic, options).available(), "modular preparation ignored zero effort");
        }
        {
            cps_secular_term_view term{"1@8000001", "0", "0", "0"};
            cps_secular *source = nullptr;
            require(cps_secular_create(&term, 1, nullptr, &source) == CPS_OUT_OF_RANGE && !source,
                    "compact secular input expanded without permission");
            certpolysolve::ResourceLimits limits;
            limits.flags = CPS_RESOURCE_ALLOW_EXPONENT_EXPANSION;
            limits.exact_materialization_bits = 9000000;
            require(cps_secular_create(&term, 1, &limits, &source) == CPS_OK, "secular exponent allowance ignored");
            cps_secular_destroy(source);
            limits.temporary_bytes = 1;
            require(cps_secular_create(&term, 1, &limits, &source) == CPS_OUT_OF_RANGE && !source,
                    "secular temporary allowance ignored");
        }
        {
            d::ComplexPolynomial p;
            acb_poly_set_coeff_si(p.value(), 0, -1);
            acb_poly_set_coeff_si(p.value(), 2, 1);
            d::RootVector seeds(2);
            acb_one(seeds[0]); acb_set_si(seeds[1], -1);
            require(!d::prepare_machine_engine(d::IterationMethod::Aberth, d::ArithmeticTier::Binary64,
                p.value(), seeds.view(), false, false, false, 1), "machine preallocation allowance ignored");
            d::PiecewiseRootStatistics stats;
            d::PiecewiseModelCursor cursor(p.value(), 16, 128, stats);
            cps_options options; options.temporary_bytes = 1;
            auto step = cursor.advance(options);
            if (step == d::PiecewiseModelStep::Work) step = cursor.advance(options);
            require(step == d::PiecewiseModelStep::Refused, "piecewise preallocation allowance ignored");
        }
        require(d::resources::precision_product(d::resources::precision_max, 2) > 0, "precision doubled through signed zero");
        require(!d::resources::storage_fits(UINT64_MAX, UINT64_MAX), "saturated estimate admitted");
        require(!cps_options{}.exact_effort_fits(UINT64_MAX), "saturated exact effort estimate admitted");
        {
            d::ComplexBall center;
            acb_one(center.value());
            acb_mul_2exp_si(center.value(), center.value(), 1000000000000L);
            require(!d::center_materialization_fits(center.value(), 1024), "center expanded before size admission");
            acb_set_si(center.value(), 3);
            require(d::center_materialization_fits(center.value(), 1024), "ordinary center refused");
        }
        cps_cleanup_thread();
        std::cout << "resource policy checks passed\n";
    } catch (const std::exception &e) { std::cerr << e.what() << '\n'; return 1; }
}
