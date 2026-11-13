// D1 (docs/real-interval-contract-note-20261015.md): the isolation goal of a real-only query (Part 1) and the real
// roots in a closed region (Part 2). Each answer is checked independently of the solver: an interval's real roots are counted by Sturm's theorem on the exact input
// (exact_reference_api.hpp) where the degree allows, and otherwise by the input's certified signs at the interval's
// dyadic ends against a count known by construction.
#include "api_test_support.hpp"
#include "exact_reference_api.hpp"
#include "internal.hpp"
#include "certify/real_region.hpp"
#include "numeric/evaluate/sparse_real_evaluator.hpp"
#include <flint/fmpq.h>
#include <chrono>
#include <iostream>
#include <tuple>

using namespace reference;
static void require(bool condition, const std::string &message) { reference::require(condition, message.c_str()); }

namespace {
// Count requests only; no allocation is denied and the original GMP allocator always performs it.
// This executable is single-threaded, and the hooks are restored before any later test group.
size_t largest_gmp_request = 0;
void *(*gmp_allocate)(size_t);
void *(*gmp_reallocate)(void *, size_t, size_t);
void (*gmp_release)(void *, size_t);
void *record_allocate(size_t n) {
    largest_gmp_request = std::max(largest_gmp_request, n);
    return gmp_allocate(n);
}
void *record_reallocate(void *p, size_t old_n, size_t n) {
    largest_gmp_request = std::max(largest_gmp_request, n);
    return gmp_reallocate(p, old_n, n);
}
struct RecordAllocations {
    RecordAllocations() {
        largest_gmp_request = 0;
        mp_get_memory_functions(&gmp_allocate, &gmp_reallocate, &gmp_release);
        mp_set_memory_functions(record_allocate, record_reallocate, gmp_release);
    }
    ~RecordAllocations() { mp_set_memory_functions(gmp_allocate, gmp_reallocate, gmp_release); }
};
}

static void region_compact_exponents() {
    RecordAllocations recording;
    for (const long exponent : {8000000, -8000000}) {
        for (const int sign : {-1, 1}) {
            certpolysolve::Polynomial p;
            p.set(0, (sign > 0 ? "-1@" : "1@") + std::to_string(exponent)).set(1, 1);
            const certpolysolve::Dyadic known(std::to_string(sign), std::to_string(exponent));
            certpolysolve::Options options;
            options.max_bits = 256;
            options.resources.emplace();
            options.resources->exact_materialization_bits = 4096;
            for (int region = 0; region < 3; ++region) {
                const auto answer = certpolysolve::real_roots(p, region == 0 ? certpolysolve::at_least(0)
                                                          : region == 1 ? certpolysolve::at_most(0)
                                                                        : certpolysolve::between(0, 1), options);
                const bool inside = region == 0 ? sign > 0 : region == 1 ? sign < 0 : sign > 0 && exponent < 0;
                require(answer.complete() && answer.goal_met() && answer.roots().size() == size_t(inside),
                        "compact interval query: " + answer.message());
                for (const auto &root : answer.roots()) {
                    require(root.count == 1 && root.multiplicity.value_or(0) == 1 &&
                                root.realness == certpolysolve::Realness::real && root.box.real.contains(known) &&
                                root.box.imag.is_point() && root.box.imag.contains(certpolysolve::Dyadic(0)),
                            "compact interval query lost its known exact root");
                }
            }
        }
    }
    // A valid coarse enclosure of the one root 1/3 of 3x - 1. Crossing the region end requires placement;
    // a wholly contained interval needs no such test, even when its other endpoint is enormous.
    namespace d = certpolysolve::detail;
    for (const bool crossing : {false, true}) {
        d::Float lo, hi;
        arf_one(hi.value());
        arf_mul_2exp_si(hi.value(), hi.value(), 8000000);
        if (crossing) arf_neg(lo.value(), hi.value());
        d::Root root;
        root.endpoints = {d::Endpoint(lo.value()), d::Endpoint(hi.value()), d::Endpoint{}, d::Endpoint{}};
        arb_set_interval_arf(acb_realref(root.box.value()), lo.value(), hi.value(), 128);
        root.count = root.multiplicity = 1;
        root.realness = CPS_REAL;
        std::vector<d::Root> roots;
        roots.push_back(std::move(root));
        cps_options options;
        options.exact_budget = 4096;
        options.max_work_bits = 256;
        options.region.real_lower = d::Number("0");
        d::RealRegionFilter filter(options, {{0, d::Number("-1")}, {1, d::Number("3")}});
        const auto outcome = filter.keep(roots, options, 128);
        require(roots.size() == 1 && outcome.undecided == size_t(crossing) && bool(outcome.refusal) == crossing,
                "compact boundary placement misreported its exact-storage refusal");
        require(arf_equal(roots[0].endpoints[0].value(), lo.value()) &&
                    arf_equal(roots[0].endpoints[1].value(), hi.value()),
                "compact boundary placement changed its retained enclosure");
    }
    // Generous headroom for backend bookkeeping, far below the million-byte request made by expanding 2^8000000.
    require(largest_gmp_request < 64 * 1024, "region filtering expanded a compact endpoint into a large rational");
}

static bool dyadic_end(const Rational &x) {
    const Integer d = boost::multiprecision::denominator(x);
    return (d & (d - 1)) == 0;
}
static bool answered_by(const certpolysolve::Result &r, const char *route) {
    return r.message().find(std::string("real roots by the ") + route) != std::string::npos;
}
// Every real root alone in its interval: each enclosure is real with its multiplicity established, its interval
// holds exactly its count of real roots by Sturm (so no other real root), the intervals are disjoint, and the
// counts make the input's real roots. `bits`, if positive, is an accuracy the enclosures must also meet.
static void verify_isolated(const Polynomial &p, const certpolysolve::Result &r, int count, int bits,
                            const std::string &label) {
    require(r.complete() && r.goal_met(), label + ": " + r.message());
    check_disjoint(r.roots());
    int total = 0;
    for (const auto &root : r.roots()) {
        require(root.realness == certpolysolve::Realness::real && root.multiplicity.value_or(0) == root.count, label + ": realness or multiplicity");
        require(roots_closed(p, endpoint(root.box.real.lo), endpoint(root.box.real.hi)) == static_cast<int>(root.count),
                label + ": an interval does not hold exactly its root");
        if (bits > 0) check_accuracy(root, bits);
        total += static_cast<int>(root.count);
    }
    require(total == count, label + ": " + std::to_string(total) + " real roots");
}

// (3x - 1)(5x - 7)(7x + 17)(x^20 + 2): three real roots that are not dyadic, and twenty roots far off the axis,
// where the selector chooses the global recipes (no attempt before the complex route).
static Polynomial dense_input() {
    Polynomial far(21, 0);
    far[0] = 2;
    far[20] = 1;
    return multiply(multiply(multiply({-1, 3}, {-7, 5}), {17, 7}), far);
}

// The complex route: the isolation goal, and the accuracy and approximation goals on the same roots. (Until 18
// October 2026 one session tightened isolation to accuracy; the 0.5 API has one-shot solves only.)
static void isolation_dense() {
    const auto p = dense_input();
    certpolysolve::Options options;
    options.goal = certpolysolve::Goal::isolation();
    options.bits = 53;
    options.max_bits = 4096;
    auto result = certpolysolve::real_roots(input(p), options);
    require(cps_result_goal(result.get()) == CPS_ROOT_ISOLATION && answered_by(result, "complex route"),
            "isolation on a dense input: " + result.message());
    verify_isolated(p, result, 3, 0, "dense isolation");
    options.goal = certpolysolve::Goal::accuracy();
    options.bits = 53;
    options.max_bits = 4096;
    result = certpolysolve::real_roots(input(p), options);
    require(cps_result_goal(result.get()) == CPS_ROOT_ACCURACY, "the accuracy solve's goal");
    verify_isolated(p, result, 3, 53, "dense accuracy");
    options.goal = certpolysolve::Goal::approximation();
    verify_isolated(p, certpolysolve::real_roots(input(p), options), 3, 53, "dense approximation");
}

// The sparse recursion on x^4096 - 2 (roots +-2^(1/4096)): its isolating intervals are published as they are,
// each a certified sign change of the input (Sturm is out of reach at this degree), and not refined to the
// default 53 bits; under the accuracy goal they are refined.
static void isolation_sparse() {
    Polynomial p(4097, 0);
    p[0] = -2;
    p[4096] = 1;
    const certpolysolve::detail::SparseRealEvaluator q({{0, certpolysolve::detail::Number("-2")}, {4096, certpolysolve::detail::Number("1")}});
    const auto sign_at = [&](const Rational &x) {
        certpolysolve::detail::RealBall point, value, slope;
        fmpq_t exact;
        fmpq_init(exact);
        fmpq_set_str(exact, x.str().c_str(), 10);
        arb_set_fmpq(point.value(), exact, 4096);
        fmpq_clear(exact);
        require(arb_is_exact(point.value()), "an endpoint is not a short dyadic");
        q.evaluate2(value.value(), slope.value(), point.value(), 512);
        return arb_is_positive(value.value()) ? 1 : arb_is_negative(value.value()) ? -1 : 0;
    };
    certpolysolve::Options options;
    options.goal = certpolysolve::Goal::isolation();
    options.bits = 53;
    options.max_bits = 4096;
    auto result = certpolysolve::real_roots(input(p), options);
    require(result.complete() && result.goal_met() && answered_by(result, "sparse recursion"),
            "isolation by the sparse recursion: " + result.message());
    require(result.roots().size() == 2, "x^4096 - 2 has two real roots");
    check_disjoint(result.roots());
    bool wide = false;
    for (const auto &root : result.roots()) {
        const Rational a = endpoint(root.box.real.lo), b = endpoint(root.box.real.hi);
        require(root.realness == certpolysolve::Realness::real && root.count == 1 && root.multiplicity.value_or(0) == 1, "sparse realness or count");
        require(sign_at(a) * sign_at(b) < 0, "an isolating interval without a sign change of the input");
        const Rational magnitude = a > 0 ? a : Rational(-b);
        wide = wide || (b - a) > power2(-53) * magnitude;
    }
    require(wide, "the isolation goal refined every interval to the default accuracy");
    options.goal = certpolysolve::Goal::accuracy();
    options.bits = 53;
    options.max_bits = 4096;
    result = certpolysolve::real_roots(input(p), options);
    require(result.complete() && result.goal_met() && answered_by(result, "sparse recursion"),
            "the accuracy goal by the sparse recursion: " + result.message());
    for (const auto &root : result.roots()) {
        check_accuracy(root, 53);
        require(sign_at(endpoint(root.box.real.lo)) * sign_at(endpoint(root.box.real.hi)) < 0,
                "a refined interval without a sign change of the input");
    }
}

// The axis charts (the selector's piecewise route): (3x - 1)(5x - 7)(7x + 17) prod_{i=1..126} (x^2 + 9^i), as in
// test_real_refinement; its real roots are known by construction (Sturm would take too long at degree 255).
static void isolation_axis() {
    auto p = multiply(multiply({-1, 3}, {-7, 5}), {17, 7});
    Rational scale = 1;
    for (int i = 1; i <= 126; ++i) {
        scale *= 9;
        p = multiply(p, {scale, 0, 1});
    }
    certpolysolve::Options options;
    options.goal = certpolysolve::Goal::isolation();
    const auto result = certpolysolve::real_roots(input(p), options);
    require(result.complete() && result.goal_met() && answered_by(result, "axis charts"),
            "isolation by the axis charts: " + result.message());
    const std::vector<Rational> known{Rational(1) / 3, Rational(7) / 5, Rational(-17) / 7};
    require(result.roots().size() == known.size(), "the axis charts' isolation lost or added a root");
    check_disjoint(result.roots());
    for (const auto &root : result.roots()) {
        int inside = 0;
        for (const auto &x : known) inside += endpoint(root.box.real.lo) <= x && x <= endpoint(root.box.real.hi);
        require(root.realness == certpolysolve::Realness::real && root.count == 1 && inside == 1, "an axis interval without its root");
    }
}

// Both branches of the axis mode's convexity test (axis_convex, piecewise_axis.cpp), which no other test and no fuzz
// case reaches: isolation_axis's input times a pair 3/5 +- i 2^-12
// beside the axis, whose part the test proves root-free in one step, and times two real roots 3/5 and 3/5 + 2^-8,
// which it counts by the part's ends. Under the debug probes (CERTPOLYSOLVE_PROBE_AXIS) the published rounds report
// "free 1" and "counted 2" (22 October 2026); the release library reports no counters, so this checks the answers
// exactly (each planted real root alone in its interval, nothing else) and that the axis charts gave them.
static void axis_convex_branches() {
    auto p = multiply(multiply({-1, 3}, {-7, 5}), {17, 7});
    Rational scale = 1;
    for (int i = 1; i <= 126; ++i) {
        scale *= 9;
        p = multiply(p, {scale, 0, 1});
    }
    const std::vector<Rational> base{Rational(1) / 3, Rational(7) / 5, Rational(-17) / 7};
    const Rational shift = Rational(1) / (1 << 8), pair = Rational(1 << 12);
    // 2^24 (5x - 3)^2 + 25: the pair 3/5 +- i 2^-12; (5x - 3)(2^8 (5x - 3) - 5): the roots 3/5 and 3/5 + 2^-8.
    const auto pair_input = multiply(p, {9 * pair * pair + 25, -30 * pair * pair, 25 * pair * pair});
    const auto close_input = multiply(p, multiply({-3, 5}, {-3 * Rational(1 << 8) - 5, 5 * Rational(1 << 8)}));
    auto close_roots = base;
    close_roots.push_back(Rational(3) / 5);
    close_roots.push_back(Rational(3) / 5 + shift);
    for (const auto &[label, q, known] : {std::tuple{"a pair beside the axis (root-free)", pair_input, base},
                                          std::tuple{"a close real pair (counted)", close_input, close_roots}}) {
        certpolysolve::Options options;
        options.goal = certpolysolve::Goal::approximation();
        const auto result = certpolysolve::real_roots(input(q), options);
        require(result.complete() && result.goal_met() && answered_by(result, "axis charts"),
                std::string(label) + ": " + result.message());
        require(result.roots().size() == known.size(), std::string(label) + ": a real root lost or added");
        check_disjoint(result.roots());
        for (const auto &root : result.roots()) {
            int inside = 0;
            for (const auto &x : known) inside += endpoint(root.box.real.lo) <= x && x <= endpoint(root.box.real.hi);
            require(root.realness == certpolysolve::Realness::real && root.count == 1 && inside == 1,
                    std::string(label) + ": an interval without exactly its root");
        }
    }
}

// Isolation with an explicit method or native secular input is a usage error
// (C7: the solve's error, no result); on inexact input it is an input-dependent refusal.
static void isolation_unsupported() {
    const auto p = input(dense_input());
    const auto refused = [](const certpolysolve::Result &r, const std::string &label) {
        require(r.stop() == certpolysolve::Stop::unsupported && !r.complete(), label + ": " + r.message());
    };
    const auto usage = [](const certpolysolve::Result &r, const std::string &label) {
        require(r.stop() == certpolysolve::Stop::error && !r.get() && !r.message().empty(), label + ": " + r.message());
    };
    certpolysolve::Options classical;
    classical.goal = certpolysolve::Goal::isolation();
    classical.method = certpolysolve::Method::classical;
    usage(certpolysolve::real_roots(p, classical), "isolation with an explicit method");
    certpolysolve::Polynomial ball;
    ball.set(0, certpolysolve::Number("-2"));
    ball.set(2, certpolysolve::ball(certpolysolve::Number("1"), certpolysolve::Number("1/1000000")));
    require(ball.ok(), "the inexact input was refused: " + ball.status());
    certpolysolve::Options real;
    real.goal = certpolysolve::Goal::isolation();
    refused(certpolysolve::real_roots(ball, real), "isolation on inexact input");
    const certpolysolve::SecularEquation secular{{.residue = 1, .pole = 1}, {.residue = 1, .pole = 2}};
    require(secular.ok(), "the secular input was refused: " + secular.status());
    usage(certpolysolve::real_roots(secular, real), "isolation on native secular input");
}

// ---- Part 2: the real roots in a closed region ----

// A bound on the moduli of p's roots (Cauchy's): the reference count of an unbounded side.
static Rational cauchy(const Polynomial &p) {
    Rational bound = 1;
    for (size_t k = 0; k + 1 < p.size(); ++k) bound = std::max(bound, Rational(1 + absolute(p[k] / p.back())));
    return bound;
}
static Rational parse(const std::string &text) {
    return text.find('/') == std::string::npos ? Rational(Integer(text)) : Rational(text);
}
// The answer holds exactly the real roots in [lower, upper] (an empty string: unbounded), by Sturm: each enclosure is
// real and alone in its interval (count = multiplicity, Sturm's count of the interval), the part of its interval in
// the region holds its root, and the counts make the region's count. A root at a dyadic end is that exact point.
static void verify_region(const Polynomial &p, const certpolysolve::Result &r, const std::string &lower,
                          const std::string &upper, int bits, const std::string &label) {
    require(r.complete() && r.goal_met(), label + ": " + r.message());
    check_disjoint(r.roots());
    const Rational a = lower.empty() ? Rational(-cauchy(p)) : parse(lower);
    const Rational b = upper.empty() ? cauchy(p) : parse(upper);
    int total = 0;
    for (const auto &root : r.roots()) {
        const Rational lo = endpoint(root.box.real.lo), hi = endpoint(root.box.real.hi);
        require(root.realness == certpolysolve::Realness::real && root.multiplicity.value_or(0) == root.count, label + ": realness or multiplicity");
        require(roots_closed(p, lo, hi) == static_cast<int>(root.count), label + ": an interval without its one root");
        require(lo <= b && hi >= a && roots_closed(p, std::max(lo, a), std::min(hi, b)) == static_cast<int>(root.count),
                label + ": a root outside the region");
        bool at_end = false;
        for (const auto &end : {a, b})
            if (lo <= end && end <= hi && evaluate(p, end) == 0) {
                at_end = true;
                if (dyadic_end(end)) require(lo == end && hi == end, label + ": a root at a dyadic end is not its point");
            }
        // A root that is not an end was refined until its interval lies in the region.
        require(at_end || (lo >= a && hi <= b), label + ": an interval still holds an end that is not its root");
        if (bits > 0 && lo != hi) check_accuracy(root, bits);
        total += static_cast<int>(root.count);
    }
    require(total == roots_closed(p, a, b), label + ": " + std::to_string(total) + " roots in the region");
}
// A region from its ends in text; an empty string leaves a side unbounded.
static certpolysolve::Region end_region(const std::string &lower, const std::string &upper) {
    certpolysolve::Region region;
    if (!lower.empty()) region.lower = certpolysolve::Number(lower);
    if (!upper.empty()) region.upper = certpolysolve::Number(upper);
    return region;
}
static certpolysolve::Result in_region(const Polynomial &p, const std::string &lower, const std::string &upper,
                                       certpolysolve::Goal goal = certpolysolve::Goal::accuracy()) {
    certpolysolve::Options options;
    options.goal = goal;
    options.bits = 53;
    options.max_bits = 4096;
    return certpolysolve::real_roots(input(p), end_region(lower, upper), options);
}

// Roots at the ends, dyadic (1/2, 2) and not (1/3), with closed ends; regions just missing roots; a point region;
// unbounded sides.
static void region_ends() {
    const auto p = multiply(multiply(multiply({-1, 2}, {-1, 3}), {-2, 1}), {1, 0, 1}); // 1/2, 1/3, 2; +-i
    const std::vector<std::tuple<std::string, std::string, size_t>> cases{
        {"1/3", "1/2", 2}, {"0", "1/3", 1}, {"", "0", 0}, {"1/2", "", 2}, {"1001/3000", "499/1000", 0},
        {"2", "2", 1}, {"-1", "10", 3}, {"", "", 3}, {"1/3", "1/3", 1}};
    for (const auto &[lower, upper, count] : cases) {
        const std::string label = "[" + lower + ", " + upper + "]";
        const auto result = in_region(p, lower, upper);
        verify_region(p, result, lower, upper, 53, label);
        require(result.roots().size() == count, label + ": " + std::to_string(result.roots().size()) + " roots");
    }
}
// A root 2^-1000 from an end: inside [1, 2] and outside [1 + 2^-999, 2]. The input's signs at the end and at the
// root's interval decide it, and the interval is cut to the region's side, without a refinement (18 October 2026;
// until then the route was solved again at doubled accuracy until about 1,000 bits).
static void region_near_end() {
    const Rational close = 1 + power2(-1000);
    const auto p = multiply(multiply({-close, 1}, {-2, 3}), {1, 0, 1});
    const auto inside = in_region(p, "1", "2");
    verify_region(p, inside, "1", "2", 53, "[1, 2] with a root 2^-1000 above 1");
    require(inside.work_bits() < 1000, "deciding a root 2^-1000 from an end refined the route to " +
                                           std::to_string(inside.work_bits()) + " bits");
    require(inside.message().find("fresh solves") == std::string::npos,
            "deciding a root by the input's signs solved the route again: " + inside.message());
    const auto beyond = Rational(1 + power2(-999)).str();
    const auto result = in_region(p, beyond, "2");
    verify_region(p, result, beyond, "2", 53, "[1 + 2^-999, 2]");
    require(result.roots().empty(), "a root 2^-1001 below the lower end was kept");
}
// The solve's sweep budget covers the refinements a region needs: under max_iterations 2 a solve of [1, 2], whose
// double root lies 2^-1000 above 1, makes at most 2 sweeps, and one that stops for the budget has spent it; without
// the budget the answer is right. (The D1 fuzz family's first campaign: each refinement had the whole budget, one
// call made 99 sweeps under 64; seed 2473157081156417207 and seven more.) The root is double because the input does
// not change sign at a root of even multiplicity: its side of the end is not read from signs, and the route is
// solved again at doubled accuracy.
static void region_iteration_budget() {
    const Rational close = 1 + power2(-1000);
    const auto p = multiply(multiply(multiply({-close, 1}, {-close, 1}), {-2, 3}), {1, 0, 1});
    certpolysolve::Options options;
    constexpr uint64_t budget = 2;
    options.bits = 53;
    options.max_bits = 4096;
    options.iteration_budget = budget;
    const auto budgeted = certpolysolve::real_roots(input(p), certpolysolve::between("1", "2"), options);
    const uint64_t sweeps = budgeted.statistics().sweeps;
    require(sweeps <= budget, "a solve made " + std::to_string(sweeps) + " sweeps under a budget of 2");
    require(cps_result_stop_reason(budgeted.get()) != CPS_STOP_ITERATION_BUDGET || sweeps == budget,
            "a solve stopped for the sweep budget without spending it");
    if (budgeted.goal_met()) verify_region(p, budgeted, "1", "2", 53, "[1, 2] under a sweep budget");
    options.iteration_budget = 0;
    const auto refined = certpolysolve::real_roots(input(p), certpolysolve::between("1", "2"), options);
    verify_region(p, refined, "1", "2", 53, "[1, 2]");
    // Deciding the double root 2^-1000 above the end takes the route solved again at doubled accuracy, each a fresh
    // solve (18 October 2026): the answer says so, and its counters cover every one of them (a fresh solve counts from
    // zero, so without the superseded solves' statistics the work counter would fall below the call's work).
    require(refined.message().find("fresh solves at doubled accuracy") != std::string::npos,
            "the refinement of [1, 2] did not take fresh solves: " + refined.message());
    // The result keeps two work counters (internal.hpp; until stage 3b-3 they were the pre-0.5 blocks' fields); the
    // 0.5 statistics publish one.
    const auto &counters = refined.get()->statistics_v2;
    const auto &progress = refined.get()->progress;
    require(counters.work_units == progress.work_units_this_call &&
                refined.statistics().work_units == progress.work_units_this_call,
            "the counters do not cover every fresh solve: work units " + std::to_string(counters.work_units) +
                " against the call's " + std::to_string(progress.work_units_this_call) + ", published " +
                std::to_string(refined.statistics().work_units));
}
// Roots placed against an end by the input's signs (RealRegionFilter::place, 18 October 2026), where refinement within
// the cap could not have decided them:
// - the D1 fuzz family's seed 94282258552801974: (2/7)(x^2 + 6x + 13)(x^2 + 1)(x + 2)(x^2 + 4) in [-2 - 2^-85, -2]
//   at a 16-bit goal under a 96-bit cap and 1,000 work units; the lower end lies 2^-85 below the root -2, the upper
//   end is the root. The fresh solves at doubled accuracy reached the cap (96 bits of accuracy under a 96-bit cap,
//   which no solve certifies) and stopped; the signs decide it in the first solve;
// - a root 2^-200 above the non-dyadic end 1/3, under a 96-bit cap: inside [1/3, 1], outside [1/3 + 2^-199, 1];
// - both ends in the root's interval, 2^-300 each side of the root: the interval is cut twice;
// - a triple root (odd multiplicity: the input changes sign there) 2^-1000 above the end 1;
// - the non-dyadic end 1/3 under a cap of 2^40 bits: the end is held at each evaluation's precision, not the cap's
//   (an end set at the cap's precision would take 128 GB);
// - a cut point that is the root: 5/16 is the first multiple of 2^-k above the end 5/16 - 2^-200, and the input
//   vanishes there, so the root is published as that point.
static void region_signs_place() {
    {
        const Polynomial p = multiply(multiply(multiply(multiply({13, 6, 1}, {1, 0, 1}), {2, 1}), {4, 0, 1}),
                                      {Rational(2) / 7});
        const auto lower = Rational(-2 - power2(-85)).str();
        certpolysolve::Options options;
        options.bits = 16;
        options.max_bits = 96;
        options.work_budget = 1000;
        const auto result = certpolysolve::real_roots(input(p), certpolysolve::between(lower, "-2"), options);
        verify_region(p, result, lower, "-2", 16, "seed 94282258552801974's region");
        require(result.roots().size() == 1, "seed 94282258552801974: the root -2 was not answered");
    }
    const Rational root = Rational(1) / 3 + power2(-200);
    const auto p = multiply(multiply({-root, 1}, {-2, 1}), {1, 0, 1});
    const auto in = [&](const std::string &lower, const std::string &upper) {
        certpolysolve::Options options;
        options.bits = 53;
        options.max_bits = 96;
        return certpolysolve::real_roots(input(p), end_region(lower, upper), options);
    };
    verify_region(p, in("1/3", "1"), "1/3", "1", 53, "[1/3, 1], a root 2^-200 above 1/3, a 96-bit cap");
    const auto beyond = Rational(Rational(1) / 3 + power2(-199)).str();
    const auto outside = in(beyond, "1");
    verify_region(p, outside, beyond, "1", 53, "[1/3 + 2^-199, 1], a 96-bit cap");
    require(outside.roots().empty(), "a root 2^-200 below the lower end was kept");
    const auto a = Rational(root - power2(-300)).str(), b = Rational(root + power2(-300)).str();
    certpolysolve::Options wide;
    wide.bits = 53;
    wide.max_bits = 4096;
    const auto both = certpolysolve::real_roots(input(p), certpolysolve::between(a, b), wide);
    verify_region(p, both, a, b, 53, "both ends 2^-300 from the root");
    require(both.roots().size() == 1, "both ends 2^-300 from the root: the root was not answered");
    const Rational close = 1 + power2(-1000);
    const auto triple = multiply(multiply(multiply({-close, 1}, {-close, 1}), {-close, 1}), {1, 0, 1});
    const auto result = in_region(triple, "1", "2");
    verify_region(triple, result, "1", "2", 53, "[1, 2] with a triple root 2^-1000 above 1");
    require(result.message().find("fresh solves") == std::string::npos,
            "a triple root was refined instead of placed by signs: " + result.message());
    {
        certpolysolve::Options huge_cap;
        huge_cap.bits = 53;
        huge_cap.max_bits = uint64_t{1} << 40;
        verify_region(p, certpolysolve::real_roots(input(p), certpolysolve::between("1/3", "1"), huge_cap), "1/3", "1", 53,
                      "[1/3, 1] under a cap of 2^40 bits");
    }
    const Rational dyadic_root = Rational(5) / 16;
    const auto q = multiply(multiply({-dyadic_root, 1}, {-2, 1}), {1, 0, 1});
    const auto below_root = Rational(dyadic_root - power2(-200)).str();
    const auto at_cut = in_region(q, below_root, "1");
    verify_region(q, at_cut, below_root, "1", 53, "[5/16 - 2^-200, 1] with the root 5/16");
    require(at_cut.roots().size() == 1 && endpoint(at_cut.roots()[0].box.real.lo) == dyadic_root &&
                endpoint(at_cut.roots()[0].box.real.hi) == dyadic_root,
            "the root 5/16 at a cut point is not published as its point");
}
// A counted cluster at a stop keeps its interval: (x + 1)^2 (x + 1 - 2^-159) under a 72-bit cap, the double root
// -1 exactly at the lower end of [-1, 1], the third root 2^-159 above it. The answer stops at the cap with one
// region of count 3 whose interval must still hold all three roots; it is not the end's point with that count
// (the D1 fuzz family, 16 October, seed 5753768794565409950: a false count on a published region).
static void region_cluster_at_end() {
    const auto p = multiply(multiply({1, 1}, {1, 1}), {1 - power2(-159), 1});
    for (const certpolysolve::Goal goal : {certpolysolve::Goal::isolation(), certpolysolve::Goal::accuracy()}) {
        certpolysolve::Options options;
        options.goal = goal;
        options.bits = 16;
        options.max_bits = 72;
        const auto result = certpolysolve::real_roots(input(p), certpolysolve::between("-1", "1"), options);
        require(!result.goal_met() && cps_result_termination(result.get()) == CPS_PRECISION_LIMIT,
                "the cluster at the end was resolved under a 72-bit cap: " + result.message());
        require(!result.roots().empty(), "the answer at the stop dropped the cluster at the end");
        for (const auto &root : result.roots()) {
            const Rational lo = endpoint(root.box.real.lo), hi = endpoint(root.box.real.hi);
            require(roots_closed(p, lo, hi) == static_cast<int>(root.count),
                    "a region at the stop does not hold exactly its count of roots");
            require(root.count == 1 || lo < hi, "a cluster was published as a point");
        }
    }
}
// The zero root is in the answer exactly when lower <= 0 <= upper: [0, +inf) for "the positive roots", which the
// caller filters.
static void region_zero_root() {
    const auto p = multiply(multiply(multiply(multiply({0, 1}, {-1, 1}), {1, 1}), {-3, 1}), {2, 0, 1});
    for (const auto &[lower, upper, count] : std::vector<std::tuple<std::string, std::string, size_t>>{
             {"0", "", 3}, {"", "0", 2}, {"0", "0", 1}, {"1/2", "", 2}}) {
        const auto result = in_region(p, lower, upper);
        verify_region(p, result, lower, upper, 53, "zero root, [" + lower + ", " + upper + "]");
        require(result.roots().size() == count, "zero root: the region's count");
    }
}
// The sparse recursion's answer filtered: x^4096 - 2, whose roots are +-2^(1/4096) = +-1.000169..., checked by the
// input's signs at the interval's ends (Sturm is out of reach).
static void region_sparse() {
    Polynomial p(4097, 0);
    p[0] = -2;
    p[4096] = 1;
    const certpolysolve::detail::SparseRealEvaluator q({{0, certpolysolve::detail::Number("-2")}, {4096, certpolysolve::detail::Number("1")}});
    const auto sign_at = [&](const Rational &x) {
        certpolysolve::detail::RealBall point, value, slope;
        fmpq_t exact;
        fmpq_init(exact);
        fmpq_set_str(exact, x.str().c_str(), 10);
        arb_set_fmpq(point.value(), exact, 4096);
        fmpq_clear(exact);
        q.evaluate2(value.value(), slope.value(), point.value(), 512);
        return arb_is_positive(value.value()) ? 1 : arb_is_negative(value.value()) ? -1 : 0;
    };
    for (const auto &[lower, upper, count] : std::vector<std::tuple<std::string, std::string, size_t>>{
             {"10001/10000", "10002/10000", 1}, {"10002/10000", "", 0}, {"", "-10002/10000", 0},
             {"-10002/10000", "-10001/10000", 1}, {"-2", "2", 2}}) {
        const std::string label = "x^4096 - 2 in [" + lower + ", " + upper + "]";
        const auto result = in_region(p, lower, upper);
        require(result.complete() && result.goal_met() && answered_by(result, "sparse recursion"),
                label + ": " + result.message());
        require(result.roots().size() == count, label + ": the count");
        for (const auto &root : result.roots()) {
            const Rational lo = endpoint(root.box.real.lo), hi = endpoint(root.box.real.hi);
            require(sign_at(lo) * sign_at(hi) < 0, label + ": an interval without a sign change");
            require((lower.empty() || lo >= parse(lower)) && (upper.empty() || hi <= parse(upper)),
                    label + ": an interval reaching out of the region");
            check_accuracy(root, 53);
        }
    }
}
// The isolation goal in a region, and the accuracy goal there: the same roots.
static void region_isolation() {
    const auto p = dense_input(); // 1/3, 7/5, -17/7
    certpolysolve::Options options;
    options.goal = certpolysolve::Goal::isolation();
    options.bits = 53;
    options.max_bits = 4096;
    auto result = certpolysolve::real_roots(input(p), certpolysolve::between("0", "2"), options);
    verify_region(p, result, "0", "2", 0, "isolation in [0, 2]");
    require(result.roots().size() == 2 && cps_result_goal(result.get()) == CPS_ROOT_ISOLATION, "isolation in [0, 2]: count");
    options.goal = certpolysolve::Goal::accuracy();
    options.bits = 53;
    options.max_bits = 4096;
    result = certpolysolve::real_roots(input(p), certpolysolve::between("0", "2"), options);
    verify_region(p, result, "0", "2", 53, "accuracy in [0, 2]");
}
// A rectangle that misses the real axis: an empty, complete answer, without work.
static void region_off_axis() {
    // Imaginary bounds: C options (the 0.5 Region is real).
    const auto options = certpolysolve_test::raw_options();
    require(cps_options_set_real_only(options.get(), 1) == CPS_OK &&
                cps_options_set_region(options.get(), "", "", "1", "2") == CPS_OK,
            "the off-axis options were refused");
    const auto result = certpolysolve_test::solve(input(dense_input()), options.get());
    require(result.complete() && result.goal_met() && result.roots().empty(), "off the axis: " + result.message());
    require(result.get()->progress.work_units_this_call == 0, "an off-axis region did work");
}
// No benchmark luck: an anchor's images under maps of the line, x -> 3x and x -> x + 1/3, with the mapped roots as
// the ends (21/5 and 26/15 are not dyadic).
static void region_maps() {
    const auto anchor = dense_input();
    const auto compose = [&](const Rational &scale, const Rational &shift) { // anchor((x - shift) / scale)
        Polynomial result{anchor.back()};
        for (size_t k = anchor.size() - 1; k-- > 0;) {
            result = multiply(result, {-shift / scale, 1 / scale});
            result[0] += anchor[k];
        }
        return result;
    };
    const auto scaled = compose(3, 0);                   // roots 1, 21/5, -51/7
    const auto shifted = compose(1, Rational(1) / 3); // roots 2/3, 26/15, -44/21
    verify_region(scaled, in_region(scaled, "1", "21/5"), "1", "21/5", 53, "scaled anchor in [1, 21/5]");
    verify_region(shifted, in_region(shifted, "2/3", "26/15"), "2/3", "26/15", 53, "shifted anchor in [2/3, 26/15]");
    verify_region(shifted, in_region(shifted, "-44/21", "2/3"), "-44/21", "2/3", 53, "shifted anchor in [-44/21, 2/3]");
}
// Usage errors: a malformed end at the setter; a lower end above the upper one, a region without real_only or with an
// explicit method, and an end the default exact storage policy will not expand, at the solve (C7, C8).
static void region_refused() {
    const auto p = input(dense_input());
    const auto reversed = certpolysolve::real_roots(p, certpolysolve::between("2", "1"));
    require(reversed.stop() == certpolysolve::Stop::error && !reversed.message().empty(),
            "a region with its ends reversed was accepted");
    const auto malformed = certpolysolve::real_roots(p, certpolysolve::between("one", "2"));
    require(malformed.stop() == certpolysolve::Stop::error && !malformed.message().empty(),
            "a malformed end was accepted");
    // A region without real_only: C options (the 0.5 header takes a region only through real_roots).
    const auto complex = certpolysolve_test::raw_options();
    require(cps_options_set_region(complex.get(), "0", "1", "", "") == CPS_OK, "the region options were refused");
    require(certpolysolve_test::solve(p, complex.get()).stop() == certpolysolve::Stop::error, "a region without real_only");
    certpolysolve::Options classical;
    classical.method = certpolysolve::Method::classical;
    require(certpolysolve::real_roots(p, certpolysolve::between("0", "1"), classical).stop() == certpolysolve::Stop::error,
            "a region with an explicit method");
    const auto beyond = certpolysolve::real_roots(p, certpolysolve::at_least("1@100000"));
    require(beyond.stop() == certpolysolve::Stop::error && beyond.message().find("exact storage policy") != std::string::npos,
            "an end beyond the exact storage policy: " + beyond.message());
}

int main() {
    certpolysolve::ThreadScope lifetime;
    // Each group's time, so a run shows that every route was exercised.
    const auto timed = [](const char *name, void (*group)()) {
        const auto begin = std::chrono::steady_clock::now();
        group();
        std::cout << name << " " << std::chrono::duration<double>(std::chrono::steady_clock::now() - begin).count()
                  << " s\n";
    };
    try {
        timed("region_compact_exponents", region_compact_exponents);
        timed("isolation_unsupported", isolation_unsupported);
        timed("isolation_dense", isolation_dense);
        timed("isolation_sparse", isolation_sparse);
        timed("isolation_axis", isolation_axis);
        timed("axis_convex_branches", axis_convex_branches);
        timed("region_refused", region_refused);
        timed("region_off_axis", region_off_axis);
        timed("region_ends", region_ends);
        timed("region_zero_root", region_zero_root);
        timed("region_near_end", region_near_end);
        timed("region_iteration_budget", region_iteration_budget);
        timed("region_signs_place", region_signs_place);
        timed("region_cluster_at_end", region_cluster_at_end);
        timed("region_maps", region_maps);
        timed("region_isolation", region_isolation);
        timed("region_sparse", region_sparse);
    } catch (const std::exception &error) {
        std::cerr << error.what() << "\n";
        return 1;
    }
    std::cout << "real interval tests passed\n";
    return 0;
}
