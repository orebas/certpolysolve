// The invariants that tie a result's termination, stop reason and goal together (docs/api-spec-0.5.md section 6.2:
// I1 goal_met implies complete; I2 CPS_SUCCESS if and only if goal_met; I3 CPS_STOP_NONE if and only if CPS_SUCCESS;
// I4's input verdicts; change C25), on solves cancelled at chosen checkpoints: a cancellation callback that answers
// yes from its N-th poll on. Until 19 October 2026 each route polled once more after its goal was met, and a
// cancellation seen there reported CPS_CANCELLED with goal_met set (the program fuzz family's first trial found it):
// on that code the sweep's last checkpoint breaks I2. Each case sweeps the last checkpoints and a spread of earlier
// ones: every method on a coefficient polynomial, a real-only query, a program and a native secular input. Then the
// zero polynomial, with and without a requested token, keeps its own verdict (CPS_ZERO_POLYNOMIAL, CPS_STOP_INPUT).
// Every result's regions also come in increasing order of their exact midpoints, real part then imaginary part
// (C26), the cancelled solves' partial answers included. The other invariants of section 6.2 (C27) are checked on
// every result too: I4's table of termination, coverage, stop reason and regions; I5 a proved count in every region;
// I6 pairwise disjoint regions; I7 a complete answer's counts summing to the roots asked for; I8 a multiplicity 0 or
// the count; I9 a real region on the axis. Last, section 5.2's tables: every combination of input kind, method,
// real-only, goal and region on a few inputs is answered, refused or a usage error exactly as the tables say.
#include "api_test_support.hpp"
#include "exact_reference_api.hpp"
#include <algorithm>
#include <functional>
#include <iostream>
#include <string>

using namespace reference;

namespace {
// C options with the precision and method of the case; the sweep adds the checkpoint at which the solve is cancelled.
using RawOptions = certpolysolve_test::RawOptions;

void ok(cps_error code, const char *what) { require(code == CPS_OK, what); }

// `expected_roots`: what a complete answer's counts sum to (the degree for all complex roots, the real roots for a
// real-only query), or -1 not to check I7.
void invariants(const certpolysolve::Result &r, const std::string &what, int64_t expected_roots = -1) {
    require(r.get() != nullptr, (what + ": the solve gave an error, not a result: " + r.message()).c_str());
    const auto stop = cps_result_stop_reason(r.get());
    const auto term = cps_result_termination(r.get());
    const std::string at = what + " (termination " + std::to_string(term) + ", stop " + std::to_string(stop) +
                           ", goal_met " + std::to_string(r.goal_met()) + "): ";
    require(!r.goal_met() || r.complete(), (at + "goal_met without complete (I1)").c_str());
    require((term == CPS_SUCCESS) == r.goal_met(), (at + "CPS_SUCCESS and goal_met disagree (I2)").c_str());
    require((stop == CPS_STOP_NONE) == (term == CPS_SUCCESS), (at + "CPS_STOP_NONE and CPS_SUCCESS disagree (I3)").c_str());
    if (term == CPS_CANCELLED)
        require(stop == CPS_STOP_CANCELLED || stop == CPS_STOP_TIME_LIMIT, (at + "a cancellation's stop reason (I4)").c_str());
    // I4's other rows.
    if (term == CPS_PRECISION_LIMIT)
        require(stop == CPS_STOP_PRECISION || stop == CPS_STOP_INPUT, (at + "a precision stop's reason (I4)").c_str());
    if (term == CPS_WORK_LIMIT)
        require(stop == CPS_STOP_WORK_BUDGET || stop == CPS_STOP_CELL_BUDGET || stop == CPS_STOP_ITERATION_BUDGET,
                (at + "a work limit's reason (I4)").c_str());
    if (term == CPS_UNSUPPORTED_QUERY || term == CPS_ZERO_POLYNOMIAL || term == CPS_INDETERMINATE_DEGREE)
        require(!r.complete() && stop == CPS_STOP_INPUT, (at + "an input verdict's coverage or reason (I4)").c_str());
    if (term == CPS_ZERO_POLYNOMIAL || term == CPS_INDETERMINATE_DEGREE)
        require(r.roots().empty() && (term != CPS_ZERO_POLYNOMIAL || r.degree() == -1), (at + "an input verdict's regions (I4)").c_str());
    int64_t counted = 0;
    for (const auto &root : r.roots()) {
        require(root.count >= 1, (at + "a region without a proved count (I5)").c_str());
        require(!root.multiplicity || *root.multiplicity == root.count, (at + "a multiplicity that is not the count (I8)").c_str());
        if (root.realness == certpolysolve::Realness::real)
            require(endpoint(root.box.imag.lo) == 0 && endpoint(root.box.imag.hi) == 0, (at + "a real region off the axis (I9)").c_str());
        counted += static_cast<int64_t>(root.count);
    }
    if (expected_roots >= 0 && r.complete())
        require(counted == expected_roots, (at + "a complete answer counts " + std::to_string(counted) + " roots of " +
                                            std::to_string(expected_roots) + " (I7)").c_str());
    // I6: pairwise disjoint as closed rectangles.
    const auto &all = r.roots();
    for (size_t i = 0; i < all.size(); ++i)
        for (size_t j = i + 1; j < all.size(); ++j) {
            const auto &a = all[i].box, &b = all[j].box;
            const bool apart = endpoint(a.real.hi) < endpoint(b.real.lo) || endpoint(b.real.hi) < endpoint(a.real.lo) ||
                               endpoint(a.imag.hi) < endpoint(b.imag.lo) || endpoint(b.imag.hi) < endpoint(a.imag.lo);
            require(apart, (at + "regions " + std::to_string(i) + " and " + std::to_string(j) + " meet (I6)").c_str());
        }
    // C26: the sum of a region's ends on an axis orders as its midpoint does; exact.
    const auto &roots = r.roots();
    for (size_t i = 1; i < roots.size(); ++i) {
        const auto &a = roots[i - 1].box, &b = roots[i].box;
        const Rational ar = endpoint(a.real.lo) + endpoint(a.real.hi), br = endpoint(b.real.lo) + endpoint(b.real.hi);
        const Rational ai = endpoint(a.imag.lo) + endpoint(a.imag.hi), bi = endpoint(b.imag.lo) + endpoint(b.imag.hi);
        require(ar < br || (ar == br && ai < bi),
                (at + "regions " + std::to_string(i - 1) + " and " + std::to_string(i) +
                 " out of exact midpoint order (C26)").c_str());
    }
}

// Solves with the solve cancelled at every one of the last `tail` checkpoints and at `spread` earlier ones (the
// testing counter of certpolysolve_testing.h answers "cancelled" from the n-th poll on, as a callback would), each checked
// against the invariants; the solve with no cancellation (counting only) must meet its goal.
using Solve = std::function<certpolysolve::Result(const cps_options *)>;
void sweep(const std::string &what, const std::function<RawOptions()> &options_for, const Solve &solve, int64_t roots,
           uint64_t tail = 48, uint64_t spread = 24) {
    auto options = options_for();
    ok(cps_testing_options_stop_at_poll(options.get(), UINT64_MAX), "counting the polls was refused");
    const auto full = solve(options.get());
    invariants(full, what + ", not cancelled", roots);
    require(full.goal_met(), (what + ": the solve without cancellation did not meet its goal").c_str());
    const uint64_t polls = cps_testing_result_polls(full.get());
    require(polls > 0, (what + ": the solve polled no checkpoint").c_str());
    std::vector<uint64_t> points;
    for (uint64_t k = 0; k < tail && k < polls; ++k) points.push_back(polls - k);
    for (uint64_t k = 1; k <= spread; ++k) points.push_back(std::max<uint64_t>(1, polls * k / (spread + 1)));
    for (const uint64_t n : points) {
        auto cancelled = options_for();
        ok(cps_testing_options_stop_at_poll(cancelled.get(), n), "stopping at a poll was refused");
        const auto result = solve(cancelled.get());
        invariants(result, what + ", cancelled at checkpoint " + std::to_string(n) + " of " + std::to_string(polls), roots);
    }
}

// Integer coefficients in a fixed pattern in [-500, 499] (test_threads' input).
certpolysolve::Polynomial patterned(uint64_t degree) {
    certpolysolve::Polynomial p;
    for (uint64_t k = 0; k <= degree; ++k)
        p.set(k, static_cast<double>(static_cast<int>((k * 7919 + 13) % 1000) - 500));
    return p;
}

void zero_polynomial() {
    const certpolysolve::Polynomial zero; // every coefficient 0
    for (const bool requested : {false, true}) {
        std::stop_source token;
        if (requested) token.request_stop();
        certpolysolve::Options options;
        options.stop = token.get_token();
        const auto result = certpolysolve::solve(zero, options);
        const std::string what = std::string("the zero polynomial") + (requested ? " with a requested token" : "");
        require(cps_result_termination(result.get()) == CPS_ZERO_POLYNOMIAL, (what + ": not CPS_ZERO_POLYNOMIAL").c_str());
        require(cps_result_stop_reason(result.get()) == CPS_STOP_INPUT, (what + ": its stop reason is not CPS_STOP_INPUT").c_str());
        invariants(result, what);
    }
}
// Section 5.2's tables (C7, C27), encoded from the spec, not from the library's own check: for every input kind, method,
// real-only, goal and region, the solve is a usage error (CPS_INVALID_ARGUMENT, no result), a refusal (CPS_OK,
// CPS_UNSUPPORTED_QUERY, CPS_STOP_INPUT) or an answer, as the tables say.
enum class Outcome { answered, refusal, usage };
enum class Kind { exact_real, exact_nonreal, ball, secular, program };
enum class Goal { by_default, accuracy, approximation, cluster, isolation };
Outcome expected(Kind kind, cps_method method, bool real_only, Goal goal, bool region) {
    const bool forced = method != CPS_METHOD_AUTO, cluster = goal == Goal::cluster, isolation = goal == Goal::isolation;
    switch (kind) {
    case Kind::exact_real:
    case Kind::exact_nonreal:
    case Kind::ball:
        if (region && !real_only) return Outcome::usage;
        if (forced && (real_only || cluster || isolation)) return Outcome::usage;
        if (real_only && cluster) return Outcome::usage;
        if ((forced || isolation) && kind == Kind::ball) return Outcome::refusal;
        if (real_only && kind != Kind::exact_real) return Outcome::refusal;
        return Outcome::answered;
    case Kind::secular:
        if (method == CPS_METHOD_CLASSICAL || method == CPS_METHOD_PIECEWISE) return Outcome::usage;
        if (isolation || region) return Outcome::usage;
        if (forced && (real_only || cluster)) return Outcome::usage;
        if (real_only && cluster) return Outcome::usage;
        return Outcome::answered; // real terms: the real-only query's numerator is real and small
    case Kind::program:
        if (method == CPS_METHOD_PIECEWISE || real_only || isolation || region) return Outcome::usage;
        return Outcome::answered;
    }
    return Outcome::usage;
}
void query_table() {
    // (x - 1)(x - 2)(x^2 + 1), x - i, and the first with a ball coefficient; a real secular equation (roots 1, 2, 4);
    // a program, (x - 1)(x - 2)(x^2 + 1) as a product.
    const auto exact = input(multiply(multiply({-1, 1}, {-2, 1}), {1, 0, 1}));
    certpolysolve::Polynomial nonreal;
    nonreal.set(0, certpolysolve::Number::complex(certpolysolve::Number("0"), certpolysolve::Number("-1"))).set(1, certpolysolve::Number("1"));
    auto ball = exact;
    ball.set(0, certpolysolve::ball(certpolysolve::Number("2"), certpolysolve::Number("1@-60")));
    const certpolysolve::SecularEquation secular{{.residue = "36/5", .pole = "-2"}, {.residue = "-4/3", .pole = "0"},
                                             {.residue = "2/15", .pole = "3"}};
    const auto program = certpolysolve::Program::from([](auto x) { return (x - 1) * (x - 2) * (x * x + 1); });
    require(ball.ok() && nonreal.ok() && secular.ok() && program.ok(), "the table's inputs were refused");
    uint64_t solves = 0, usage = 0, refusals = 0;
    for (const Kind kind : {Kind::exact_real, Kind::exact_nonreal, Kind::ball, Kind::secular, Kind::program})
        for (const cps_method method : {CPS_METHOD_AUTO, CPS_METHOD_CLASSICAL, CPS_METHOD_SECULAR, CPS_METHOD_PIECEWISE})
            for (const bool real_only : {false, true})
                for (const Goal goal : {Goal::by_default, Goal::accuracy, Goal::approximation, Goal::cluster, Goal::isolation})
                    for (const bool region : {false, true}) {
                        auto o = certpolysolve_test::raw_options();
                        ok(cps_options_set_method(o.get(), method), "the method was refused");
                        ok(cps_options_set_real_only(o.get(), real_only), "real-only was refused");
                        if (goal == Goal::accuracy) ok(cps_options_set_goal(o.get(), CPS_ROOT_ACCURACY), "the goal was refused");
                        if (goal == Goal::approximation) ok(cps_options_set_goal(o.get(), CPS_ROOT_APPROXIMATION), "the goal was refused");
                        if (goal == Goal::isolation) ok(cps_options_set_goal(o.get(), CPS_ROOT_ISOLATION), "the goal was refused");
                        if (goal == Goal::cluster) ok(cps_options_set_cluster_radius(o.get(), "1e-3"), "the radius was refused");
                        if (region) ok(cps_options_set_region(o.get(), "0", "3", nullptr, nullptr), "the region was refused");
                        const certpolysolve::Result r = kind == Kind::exact_real    ? certpolysolve_test::solve(exact, o.get())
                                                    : kind == Kind::exact_nonreal ? certpolysolve_test::solve(nonreal, o.get())
                                                    : kind == Kind::ball         ? certpolysolve_test::solve(ball, o.get())
                                                    : kind == Kind::secular      ? certpolysolve_test::solve(secular, o.get())
                                                                                 : certpolysolve_test::solve(program, o.get());
                        const Outcome want = expected(kind, method, real_only, goal, region);
                        const std::string what = "kind " + std::to_string(static_cast<int>(kind)) + ", method " +
                                                 std::to_string(method) + ", real-only " + std::to_string(real_only) +
                                                 ", goal " + std::to_string(static_cast<int>(goal)) + ", region " +
                                                 std::to_string(region);
                        ++solves;
                        if (want == Outcome::usage) {
                            ++usage;
                            require(r.stop() == certpolysolve::Stop::error && !r.get() && !r.message().empty(),
                                    (what + ": not a usage error: " + r.message()).c_str());
                            continue;
                        }
                        require(r.get() != nullptr, (what + ": an error, not an answer: " + r.message()).c_str());
                        invariants(r, what);
                        const bool refused = cps_result_termination(r.get()) == CPS_UNSUPPORTED_QUERY;
                        if (want == Outcome::refusal) ++refusals;
                        require(refused == (want == Outcome::refusal),
                                (what + (refused ? ": refused, but the table answers it: " : ": answered, but the table refuses it: ") +
                                 r.message()).c_str());
                    }
    std::cout << "section 5.2: " << solves << " queries, " << usage << " usage errors, " << refusals << " refusals\n";
}
} // namespace

int main() {
    certpolysolve::ThreadScope lifetime;
    try {
        const auto p = patterned(24);
        for (const cps_method method : {CPS_METHOD_AUTO, CPS_METHOD_CLASSICAL, CPS_METHOD_SECULAR, CPS_METHOD_PIECEWISE})
            sweep("method " + std::to_string(method), [method] {
                auto o = certpolysolve_test::raw_options();
                ok(cps_options_set_precision(o.get(), 64), "precision was refused");
                ok(cps_options_set_max_work_bits(o.get(), 4096), "the cap was refused");
                ok(cps_options_set_method(o.get(), method), "the method was refused");
                return o;
            }, [&](const cps_options *o) { return certpolysolve_test::solve(p, o); }, 24);
        // (x^2 - 2)(x^2 + 1)(3x - 1): two real roots beside a complex pair, real-only
        const auto real_input = input(multiply(multiply({-2, 0, 1}, {1, 0, 1}), {-1, 3}));
        sweep("a real-only query", [] {
            auto o = certpolysolve_test::raw_options();
            ok(cps_options_set_precision(o.get(), 64), "precision was refused");
            ok(cps_options_set_max_work_bits(o.get(), 4096), "the cap was refused");
            ok(cps_options_set_real_only(o.get(), 1), "real-only was refused");
            return o;
        }, [&](const cps_options *o) { return certpolysolve_test::solve(real_input, o); }, 3);
        // Mandelbrot's p_3 times the Chebyshev polynomial T_4: a product, so the product lift runs
        const auto program = certpolysolve::Program::from([](auto x) {
            auto mandelbrot = x * 0 + 1;
            for (int j = 0; j < 3; ++j) mandelbrot = x * mandelbrot * mandelbrot + 1;
            auto previous = x * 0 + 1, chebyshev = x;
            for (int k = 1; k < 4; ++k) {
                auto next = 2 * x * chebyshev - previous;
                previous = chebyshev;
                chebyshev = next;
            }
            return mandelbrot * chebyshev;
        });
        sweep("a program", [] {
            auto o = certpolysolve_test::raw_options();
            ok(cps_options_set_precision(o.get(), 64), "precision was refused");
            ok(cps_options_set_max_work_bits(o.get(), 4096), "the cap was refused");
            return o;
        }, [&](const cps_options *o) { return certpolysolve_test::solve(program, o); }, 11);
        // test_d3's native secular equation, whose roots are 1, 2 and 4
        const certpolysolve::SecularEquation secular{{.residue = "36/5", .pole = "-2"},
                                                 {.residue = "-4/3", .pole = "0"},
                                                 {.residue = "2/15", .pole = "3"}};
        require(secular.ok(), "the secular equation was refused");
        sweep("a native secular equation", [] {
            auto o = certpolysolve_test::raw_options();
            ok(cps_options_set_precision(o.get(), 64), "precision was refused");
            ok(cps_options_set_max_work_bits(o.get(), 4096), "the cap was refused");
            return o;
        }, [&](const cps_options *o) { return certpolysolve_test::solve(secular, o); }, 3);
        zero_polynomial();
        query_table();
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
    std::cout << "result invariant tests passed\n";
    return 0;
}
