#include "solve/prepared_piecewise.hpp"
#include "numeric/arithmetic/polynomial_arithmetic.hpp"
#include "exact_winding_api.hpp"
#include "api_test_support.hpp"
#include "numeric/exact/exact_preprocessing.hpp"
#include "core/heuristics.hpp"
#include "methods/piecewise/piecewise_roots.hpp"
#include <chrono>
#include <vector>
#include <numbers>
#include <cmath>
#include <iostream>
#include <thread>

namespace d = certpolysolve::detail;
using namespace reference;
static Rational endpoint(const d::Endpoint &value) {
    return reference::endpoint(certpolysolve::Dyadic{value.mantissa_text(), value.exponent_text()});
}
static void known_roots(const std::vector<Complex> &known, const char *absolute_tolerance = "0",
                        cps_goal goal = CPS_ROOT_APPROXIMATION) {
    const auto polynomial = from_roots(known);
    auto input = reference::input(polynomial);
    cps_options options;
    options.goal = goal;
    options.max_work_bits = 2048;
    options.absolute_tolerance = d::Number(absolute_tolerance);
    d::PiecewiseRootStatistics statistics;
    const auto result = d::solve_prepared_piecewise_roots(certpolysolve::detail::dense(*input.get()), options, statistics);
    require(result.complete && result.goal_met && d::disjoint(result.roots), "prepared cover incomplete/overlapping");
    uint64_t total = 0;
    for (const auto &root : result.roots) {
        uint64_t count = 0;
        for (const auto &z : known)
            if (endpoint(root.endpoints[0]) <= z.real && z.real <= endpoint(root.endpoints[1]) &&
                endpoint(root.endpoints[2]) <= z.imag && z.imag <= endpoint(root.endpoints[3])) ++count;
        require(count == root.count, "known exact roots disagree with the restored count");
        // A region of one factor's root carries its exact multiplicity. Overlapping factor covers are merged at once
        // (W8's S8, rule 1), and a merged region is a counted cluster with no multiplicity (0), which only the
        // approximation goal accepts; under the accuracy goal every region is one root with its multiplicity.
        require(root.multiplicity == 0 ? goal == CPS_ROOT_APPROXIMATION : count == root.multiplicity,
                "known exact roots disagree with the restored multiplicity");
        require(d::meets_goal(root, options, static_cast<slong>(result.work_bits)), "published accuracy lost");
        total += count;
    }
    require(total == known.size(), "prepared cover lost roots");
}
static void quadratic_cancellation() {
    for (int exponent : {3000, 40000}) {
        const Rational a = power2(exponent), t = power2(-exponent);
        const Polynomial factor{1, -a, 1};
        const auto input = reference::input(multiply(factor, factor));
        // Independent dyadic brackets: f(t)>0, f(t+2t³)<0, and f is
        // decreasing there. The other root follows from their exact sum a.
        const Rational lower = t, upper = t + 2 * t * t * t;
        require(evaluate(factor, lower) > 0 && evaluate(factor, upper) < 0 && upper < a / 2,
                "invalid independent quadratic brackets");
        for (const auto goal : {CPS_ROOT_APPROXIMATION, CPS_ROOT_ACCURACY}) {
            cps_options options;
            options.goal = goal;
            options.max_work_bits = 256;
            d::PiecewiseRootStatistics statistics;
            const auto result = d::solve_prepared_piecewise_roots(certpolysolve::detail::dense(*input.get()), options, statistics);
            require(result.complete && result.goal_met && result.roots.size() == 2 && d::disjoint(result.roots),
                    "quadratic subtraction lost the small root within its precision allowance");
            unsigned small = 0, large = 0;
            for (const auto &root : result.roots) {
                const Rational lo = endpoint(root.endpoints[0]), hi = endpoint(root.endpoints[1]);
                require(root.count == 2 && root.multiplicity == 2 && root.realness == CPS_REAL &&
                            endpoint(root.endpoints[2]) == 0 && endpoint(root.endpoints[3]) == 0,
                        "quadratic formula lost realness or multiplicity");
                if (hi < a / 2) {
                    require(lo <= lower && upper <= hi && evaluate(factor, lo) > 0 && evaluate(factor, hi) < 0,
                            "small root failed exact containment/sign checks");
                    ++small;
                } else {
                    require(lo > a / 2 && lo <= a - upper && a - lower <= hi &&
                                evaluate(factor, lo) < 0 && evaluate(factor, hi) > 0,
                            "large root failed exact containment/sign checks");
                    ++large;
                }
            }
            require(small == 1 && large == 1, "quadratic cover duplicated a root");
            const auto wrapped = certpolysolve_test::published(result);
            for (const auto &root : wrapped.roots())
                check_accuracy(root, 53);
        }
    }
}
// Piecewise decides by counted work, never by the clock (fuzz finding F6, seed
// 4088092304110591797: 22 real roots and a four-root cluster 2^-38 wide, 100
// bits under a 256-bit cap). The regime exit and the source finisher's
// budgets compared elapsed time, so the same solve ended differently from
// run to run (complete through the exit in 6 of 36 idle runs and 15 of 36
// loaded ones). A slow observer stands in for a loaded machine here: every
// decision, and with it every counter of work, must be the same with it.
static void decisions_do_not_read_the_clock() {
    std::vector<Complex> known;
    const auto at = [&](long numerator, long denominator) { known.emplace_back(Rational(numerator) / denominator); };
    at(199, 256); at(-37, 772); at(-105, 233); at(751, 1024); at(247, 256); at(-20, 17); at(214, 109);
    at(-53875, 65536); at(25, 32); at(1189, 2048); at(27489, 131072); at(27839, 32768); at(-373, 512);
    at(51, 64); at(-1, 1); at(-42467, 65536); at(28, 45); at(-42175, 65536); at(0, 1); at(-27, 128);
    at(83, 512); at(-1415, 2048); at(3, 4);
    for (const int offset : {3, 1, -1}) known.emplace_back(Rational(3) / 4 + offset * power2(-38));
    const auto input = reference::input(from_roots(known));
    struct Outcome {
        cps_termination termination;
        size_t regions;
        uint64_t regime_exits, rounds, certificates, finisher_attempts, finisher_certificates, finisher_refusals;
        double work;
        bool operator==(const Outcome &) const = default;
    };
    const auto solve = [&](std::chrono::milliseconds pause) {
        auto config = d::maintained_piecewise_config();
        config.observe = [pause](const d::PiecewiseRootStatistics &, const char *) {
            std::this_thread::sleep_for(pause);
        };
        cps_options options;
        options.goal = CPS_ROOT_APPROXIMATION;
        options.bits = 100;
        options.max_work_bits = 256;
        options.work_units = 100000;
        d::PiecewiseRootStatistics statistics;
        const auto result = d::solve_prepared_piecewise_roots(certpolysolve::detail::dense(*input.get()), options, statistics, config);
        return Outcome{result.termination, result.roots.size(), statistics.regime_exits, statistics.rounds,
                       statistics.certificates, statistics.finisher_attempts, statistics.finisher_certificates,
                       statistics.finisher_refusals, d::counted_work(statistics.work)};
    };
    const auto idle = solve(std::chrono::milliseconds(0));
    require(idle.work > 0, "explicit piecewise counted no work");
    require(solve(std::chrono::milliseconds(100)) == idle, "a slow observer changed piecewise's decisions");
    require(solve(std::chrono::milliseconds(0)) == idle, "two identical piecewise solves differ");
}
// The precision ladder's last rung is the cap (queue item 3.2). Twenty integer
// roots under a 160-bit cap: the ladder's base rung for 53 bits is 149, the
// next rung (181) would pass the cap, so the solve used to stop for precision
// at 149 bits; the rung is now cut to the cap, where the solve finishes. The
// regime exit is off so piecewise finishes by itself.
static void ladder_ends_at_the_cap() {
    std::vector<Complex> known;
    for (int k = 1; k <= 20; ++k) known.emplace_back(k);
    const auto input = reference::input(from_roots(known));
    auto config = d::maintained_piecewise_config();
    config.regime_exit = false;
    cps_options options;
    options.goal = CPS_ROOT_APPROXIMATION;
    options.bits = 53;
    options.max_work_bits = 160;
    options.work_units = 100000;
    require(d::piecewise_start_precision(options.bits) == 149, "the ladder's base rung for 53 bits is not 149");
    d::PiecewiseRootStatistics statistics;
    const auto result = d::solve_prepared_piecewise_roots(certpolysolve::detail::dense(*input.get()), options, statistics, config);
    require(result.goal_met && result.work_bits == 160,
            "the ladder did not reach the 160-bit cap, or did not finish there");
}
// Auto's exit (PiecewiseRootConfig::regime_exit_estimate): where the selector
// chose piecewise by its estimate, the state leaves once its counted work
// passes the estimate at its overhead plus the global budget,
// initial_sweep_allowance x degree x (degree + support) operations, all in
// binary64 operations. Twenty integer roots are far outside piecewise's
// regime; given an estimate of a quarter of the budget the state leaves
// unfinished, just past that bound. Asked for by the caller (no estimate),
// the same solve finishes inside piecewise, its counted work many times the
// bound.
static void selected_route_leaves_at_the_selectors_bound() {
    std::vector<Complex> known;
    for (int k = 1; k <= 20; ++k) known.emplace_back(k);
    const auto input = reference::input(from_roots(known));
    const uint64_t budget = d::heuristics::initial_sweep_allowance * 20 * (20 + 21);
    const uint64_t estimate = budget / 4;
    const double bound = d::binary64_operation_work *
                         (d::estimate_operation_overhead * static_cast<double>(estimate) + static_cast<double>(budget));
    cps_options options;
    options.goal = CPS_ROOT_APPROXIMATION;
    options.bits = 53;
    options.work_units = 100000;
    const auto solve = [&](uint64_t selected, d::PiecewiseRootStatistics &statistics) {
        auto config = d::maintained_piecewise_config();
        config.regime_exit_estimate = selected;
        return d::solve_prepared_piecewise_roots(certpolysolve::detail::dense(*input.get()), options, statistics, config);
    };
    d::PiecewiseRootStatistics asked, chosen;
    const auto finished = solve(0, asked);
    require(finished.goal_met && !asked.regime_exits, "the caller's piecewise solve of twenty integer roots left");
    require(d::counted_work(asked.work) > 4 * bound, "twenty integer roots are not far outside the selector's bound");
    const auto left = solve(estimate, chosen);
    require(!left.goal_met && chosen.regime_exits == 1, "the selector's route did not leave at its bound");
    require(d::counted_work(chosen.work) > bound && d::counted_work(chosen.work) < 2 * bound,
            "the selector's route left away from its bound");
}
// The round for the last roots (PiecewiseRootConfig::seeds_before_recovery).
// On the Swinnerton-Dyer polynomial of six primes (degree 64: the roots are
// the sums of +-sqrt(p) over p = 2, 3, 5, 7, 11, 13) a round leaves roots
// undiscovered, the deflated discovery seeds them, and in that same round a
// discovery recovery falls due. Tried first, the seeds finish the solve with
// no recovery round; with the switch off the recovery round runs first, as
// before. (On resultant_30, degree 900, that round was 23 of 26 s for one
// missing root.)
static void seeds_are_tried_before_a_recovery() {
    // P(x) -> P(x - sqrt(p)) P(x + sqrt(p)) = A^2 - p B^2, with P(x + y) = sum_m y^m P^(m)(x) / m!,
    // A the even m (y^m = p^(m/2)) and B the odd ones.
    d::IntegerPolynomial polynomial;
    fmpz_poly_set_coeff_si(polynomial.value(), 1, 1);
    for (const ulong prime : {2u, 3u, 5u, 7u, 11u, 13u}) {
        d::IntegerPolynomial even, odd, term, square;
        d::Integer power;
        fmpz_one(power.value());
        fmpz_poly_set(term.value(), polynomial.value());
        for (slong m = 0; !fmpz_poly_is_zero(term.value()); ++m) {
            d::IntegerPolynomial scaled;
            fmpz_poly_scalar_mul_fmpz(scaled.value(), term.value(), power.value());
            fmpz_poly_add(m % 2 ? odd.value() : even.value(), m % 2 ? odd.value() : even.value(), scaled.value());
            if (m % 2) fmpz_mul_ui(power.value(), power.value(), prime);
            fmpz_poly_derivative(term.value(), term.value());
            fmpz_poly_scalar_divexact_si(term.value(), term.value(), m + 1);
        }
        fmpz_poly_mul(square.value(), odd.value(), odd.value());
        fmpz_poly_scalar_mul_ui(square.value(), square.value(), prime);
        fmpz_poly_mul(polynomial.value(), even.value(), even.value());
        fmpz_poly_sub(polynomial.value(), polynomial.value(), square.value());
    }
    require(fmpz_poly_degree(polynomial.value()) == 64, "the Swinnerton-Dyer fixture is not of degree 64");
    certpolysolve::detail::DensePolynomial input;
    input.coefficients.resize(65);
    for (slong k = 0; k <= 64; ++k)
        fmpq_set_fmpz(input.coefficients[static_cast<size_t>(k)].real.significand.value(), polynomial.value()->coeffs + k);
    cps_options options;
    options.goal = CPS_ROOT_APPROXIMATION;
    options.bits = 53;
    options.max_work_bits = 65536;
    options.work_units = 100000;
    const auto solve = [&](bool seeds_first, d::PiecewiseRootStatistics &statistics) {
        auto config = d::maintained_piecewise_config();
        config.seeds_before_recovery = seeds_first;
        return d::solve_prepared_piecewise_roots(input, options, statistics, config);
    };
    d::PiecewiseRootStatistics before, now;
    const auto recovered = solve(false, before);
    require(recovered.goal_met && before.deflation_seeds && before.recovery_rounds == 1,
            "the fixture no longer seeds missing roots in a round where a recovery falls due");
    const auto seeded = solve(true, now);
    require(seeded.goal_met && seeded.roots.size() == 64, "the seeds did not finish the solve");
    require(now.deflation_seeds && now.recovery_rounds == 0, "a recovery round ran before the seeds were tried");
}
// The certificate attempt before piecewise (PiecewiseRootConfig::
// certificate_attempt): the global recipes' first step, made once by the
// prepared state on the squarefree nonzero part.
static std::vector<Complex> two_circles() {
    // 100 rational points on each of the circles of radius 1 and 5/4 (the
    // fixture of test_certificate_attempt: simple roots whose certificates
    // need more than the first precision), two real roots and a triple zero.
    std::vector<Complex> planted;
    for (const auto &radius : {Rational(1), Rational(5) / 4})
        for (int k = 1; k <= 25; ++k) {
            const Rational t = Rational(k) / 26, scale = radius / (1 + t * t);
            const Rational x = (1 - t * t) * scale, y = 2 * t * scale;
            for (const auto &z : {Complex(x, y), Complex(-x, y), Complex(x, -y), Complex(-x, -y)})
                planted.push_back(z);
        }
    planted.push_back(Complex(Rational(3) / 2));
    planted.push_back(Complex(Rational(-7) / 4));
    return planted;
}
// Every region holds exactly the planted roots its count says; returns the regions that are real.
static uint64_t check_planted(const cps_result &result, const std::vector<Complex> &planted, const char *what) {
    require(result.complete && d::disjoint(result.roots), what);
    uint64_t total = 0, real = 0;
    for (const auto &root : result.roots) {
        uint64_t count = 0;
        for (const auto &z : planted)
            if (endpoint(root.endpoints[0]) <= z.real && z.real <= endpoint(root.endpoints[1]) &&
                endpoint(root.endpoints[2]) <= z.imag && z.imag <= endpoint(root.endpoints[3])) ++count;
        require(count == root.count, what);
        total += count;
        real += root.realness == CPS_REAL;
    }
    require(total == planted.size(), what);
    return real;
}
static void attempt_before_piecewise() {
    cps_options options;
    options.goal = CPS_ROOT_APPROXIMATION;
    options.bits = 53;
    options.max_work_bits = 4096;
    auto with_attempt = d::maintained_piecewise_config();
    with_attempt.certificate_attempt = true;
    {
        // The attempt finishes the solve: no model is built, the zero root
        // and the realness of the two real roots are reported, and the
        // precision is the certificates'.
        auto planted = two_circles();
        planted.insert(planted.end(), 3, Complex(0));
        const auto input = reference::input(from_roots(planted));
        d::PreparedPiecewiseState state(certpolysolve::detail::dense(*input.get()), with_attempt);
        const auto result = state.advance(options);
        require(result.goal_met && result.termination == CPS_SUCCESS, "the attempt before piecewise did not finish an easy input");
        require(check_planted(result, planted, "the attempt's cover does not hold the planted roots") == 3,
                "the attempt's cover does not mark exactly the two real roots and the zero root real");
        const auto &statistics = state.statistics();
        require(statistics.certificate_attempts == 1 && statistics.certificate_attempt_finishes == 1 &&
                    statistics.models == 0 && statistics.rounds == 0,
                "the solve was not finished by the attempt alone");
        require(result.work_bits > 85 && result.statistics_v3.validation_calls == 0 &&
                    result.statistics_v3.candidate_sweeps != 0,
                "the attempt's precision or observations are missing from the result");
        for (const auto &root : result.roots)
            require(d::meets_goal(root, options, static_cast<slong>(result.work_bits)), "a region of the attempt misses the goal");
        // The same query again is answered from the cover, at no work.
        const auto again = state.advance(options);
        require(again.goal_met && again.progress.work_units_this_call == 0 && again.roots.size() == result.roots.size() &&
                    state.statistics().certificate_attempts == 1,
                "a repeated query did not reuse the attempt's cover");
        // One work unit a call: the preparation, then the attempt.
        d::PreparedPiecewiseState stepped(certpolysolve::detail::dense(*input.get()), with_attempt);
        auto unit = options;
        unit.work_units = 1;
        const auto first = stepped.advance(unit);
        require(!first.goal_met && stepped.statistics().certificate_attempts == 0, "the preparation's unit also made the attempt");
        const auto second = stepped.advance(unit);
        require(second.goal_met && second.progress.work_units_this_call == 1 &&
                    stepped.statistics().certificate_attempt_finishes == 1,
                "the second unit did not make the attempt");
        // Without the flag the same input is piecewise's, as before.
        d::PreparedPiecewiseState plain(certpolysolve::detail::dense(*input.get()), d::maintained_piecewise_config());
        const auto piecewise = plain.advance(options);
        require(piecewise.goal_met && plain.statistics().certificate_attempts == 0 && plain.statistics().models != 0,
                "a state without the flag made the attempt");
        check_planted(piecewise, planted, "piecewise's cover does not hold the planted roots");
    }
    {
        // A refusal leaves the solve to piecewise. Two sweeps of the call's
        // iteration budget are not enough for the attempt's packet (its
        // rounds then get none), and piecewise's local refinement does not
        // count against that budget. The attempt's observations stay in the
        // result, and the attempt is not made again.
        const auto planted = two_circles();
        const auto input = reference::input(from_roots(planted));
        d::PreparedPiecewiseState state(certpolysolve::detail::dense(*input.get()), with_attempt);
        auto budgeted = options;
        budgeted.max_iterations = 2;
        const auto result = state.advance(budgeted);
        require(result.goal_met, "piecewise did not finish after a refused attempt");
        check_planted(result, planted, "the cover after a refused attempt does not hold the planted roots");
        require(state.statistics().certificate_attempts == 1 && state.statistics().certificate_attempt_finishes == 0 &&
                    state.statistics().models != 0 && result.statistics_v3.candidate_sweeps == 2,
                "the refused attempt was not followed by piecewise, or its observations were dropped");
        const auto again = state.advance(options);
        require(again.goal_met && state.statistics().certificate_attempts == 1, "a refused attempt was made again");
    }
    {
        // Outside binary64's range (roots of modulus 2^-300 and 2^300) the
        // attempt is refused before any packet: the scaled packet is not
        // tried here.
        std::vector<Complex> planted;
        constexpr long long grid = 1ll << 20;
        constexpr int twelfth = 12;
        for (const auto &radius : {power2(-300), power2(300)})
            for (int k = 0; k < twelfth; ++k) {
                const double angle = std::numbers::pi / 2 * (k + 0.5) / twelfth;
                const Rational t = Rational(std::llround(std::tan(angle / 2) * grid)) / grid;
                const Rational x = radius * (1 - t * t) / (1 + t * t), y = radius * 2 * t / (1 + t * t);
                // Conjugate pairs in the right half plane only: with their
                // mirrors in the imaginary axis the polynomial is one in x^2,
                // which the preparation reduces before any attempt.
                for (const auto &z : {Complex(x, y), Complex(x, -y)})
                    planted.push_back(z);
            }
        const auto input = reference::input(from_roots(planted));
        d::PreparedPiecewiseState state(certpolysolve::detail::dense(*input.get()), with_attempt);
        const auto result = state.advance(options);
        require(state.statistics().certificate_attempts == 1 && state.statistics().certificate_attempt_finishes == 0 &&
                    result.statistics_v3.candidate_sweeps == 0,
                "an input outside binary64's range was given a packet before piecewise");
        require(result.goal_met, "piecewise did not finish an input outside binary64's range");
        check_planted(result, planted, "the out-of-range cover does not hold the planted roots");
    }
}

// W8's S8, rule 1: overlapping covers of two factors are merged at once, as the global driver's aggregate merges.
// The double root 1/3 and the cubic's root 2^-200 from it meet a loose tolerance as one counted cluster (count 3, no
// multiplicity, certificate complete-factor-cover-union) without the separation rounds the old rule spent.
static void overlapping_factors_merge_at_once() {
    const std::vector<Complex> known{Complex(Rational(1) / 3), Complex(Rational(1) / 3),
                                     Complex(Rational(1) / 3 + power2(-200)), Complex(2), Complex(3)};
    auto input = reference::input(from_roots(known));
    cps_options options;
    options.goal = CPS_ROOT_APPROXIMATION;
    options.max_work_bits = 2048;
    options.absolute_tolerance = d::Number("1@-100");
    d::PiecewiseRootStatistics statistics;
    const auto result = d::solve_prepared_piecewise_roots(certpolysolve::detail::dense(*input.get()), options, statistics);
    require(result.complete && result.goal_met, "the merged cover does not meet the loose goal");
    const auto merged = std::count_if(result.roots.begin(), result.roots.end(), [](const d::Root &root) {
        return root.count == 3 && root.multiplicity == 0 && root.certificate == "complete-factor-cover-union";
    });
    require(merged == 1 && result.roots.size() == 3, "the overlapping factors were not merged into one cluster");
}

// W8's S8, rule 2: when two factors' covers meet, only the factor whose region is the wider gets a tighter goal, as in
// the global driver; a factor that meets no other keeps the caller's goal. The simple roots 1/3 + 2^-80 and 2 and the
// double root 1/3 are refined until their covers separate under the accuracy goal; the triple root -22/3, its own
// factor 3x + 22, is solved once at the caller's goal (the closed form at 53 + 32 bits: -22/3 is no dyadic, so its
// ball keeps a radius), so its region stays far wider than the separated ones: under the old rule every factor's goal
// doubled with each round.
static void only_the_meeting_factor_is_tightened() {
    std::vector<Complex> known{Complex(Rational(1) / 3), Complex(Rational(1) / 3),
                               Complex(Rational(1) / 3 + power2(-80)), Complex(2)};
    const Rational far = Rational(-22) / 3;
    known.insert(known.end(), 3, Complex(far));
    auto input = reference::input(from_roots(known));
    cps_options options;
    options.goal = CPS_ROOT_ACCURACY;
    options.max_work_bits = 2048;
    d::PiecewiseRootStatistics statistics;
    const auto result = d::solve_prepared_piecewise_roots(certpolysolve::detail::dense(*input.get()), options, statistics);
    require(result.complete && result.goal_met && result.roots.size() == 4, "the separated cover is not complete");
    bool seen = false;
    for (const auto &root : result.roots) {
        const auto width = endpoint(root.endpoints[1]) - endpoint(root.endpoints[0]);
        if (endpoint(root.endpoints[0]) <= far && far <= endpoint(root.endpoints[1])) {
            seen = root.count == 3 && root.multiplicity == 3;
            require(width > power2(-110), "a factor that meets no other was tightened");
        } else if (root.count == 1 || root.multiplicity == 2) {
            require(width < power2(-80), "the meeting factors were not separated");
        }
    }
    require(seen, "the triple root -22/3 is not one region of multiplicity 3");
}

int main() {
    certpolysolve::ThreadScope lifetime;
    try {
        attempt_before_piecewise();
        overlapping_factors_merge_at_once();
        only_the_meeting_factor_is_tightened();
        quadratic_cancellation();
        ladder_ends_at_the_cap();
        seeds_are_tried_before_a_recovery();
        selected_route_leaves_at_the_selectors_bound();
        decisions_do_not_read_the_clock();
        // Public engineering fast path: an already square-free elementary
        // nonzero part also uses the proved formulas, including stripped zeros.
        for (auto known : {std::vector<Complex>{{1}, {2}}, std::vector<Complex>{{0}, {0}, {1}, {2}}}) {
            auto elementary = reference::input(from_roots(known));
            cps_options o;
            d::PiecewiseRootStatistics stats;
            auto r = d::solve_prepared_piecewise_roots(certpolysolve::detail::dense(*elementary.get()), o, stats);
            require(r.goal_met && stats.models == 0 && stats.rounds == 0,
                    "elementary source unnecessarily constructed local models");
            known_roots(known);
        }
        // The cubic's completed box initially overlaps the separate repeated root. Under the loose tolerance the
        // merged region (the double root and the cubic's root 2^-200 from it, a counted cluster) meets the goal;
        // under the accuracy goal the factors are refined until their covers separate.
        for (const auto goal : {CPS_ROOT_APPROXIMATION, CPS_ROOT_ACCURACY})
            known_roots({Complex(Rational(1) / 3), Complex(Rational(1) / 3),
                         Complex(Rational(1) / 3 + power2(-200)), Complex(2), Complex(3)}, "1@-100", goal);
        // (36x² + 12(1+i)x + 9+2i)² has these exact roots. Its
        // discriminant is negative real, but rational coefficient balls can
        // straddle the principal square-root branch cut at every precision.
        const Complex upper(Rational(-1) / 6, Rational(1) / 3);
        const Complex lower(Rational(-1) / 6, Rational(-2) / 3);
        for (const auto goal : {CPS_ROOT_APPROXIMATION, CPS_ROOT_ACCURACY})
            for (const auto &rotation : {Complex(1), Complex(-1), Complex(0, 1), Complex(3, 4)})
                known_roots({upper * rotation, upper * rotation, lower * rotation, lower * rotation}, "0", goal);
        known_roots(std::vector<Complex>(128, Complex(1)));
        std::vector<Complex> mixed(17, Complex(0));
        mixed.insert(mixed.end(), 64, Complex(1));
        mixed.insert(mixed.end(), 31, Complex(-2));
        mixed.insert(mixed.end(), 13, Complex(0, 1));
        mixed.insert(mixed.end(), 13, Complex(0, -1));
        known_roots(mixed);
        mixed.assign(64, Complex(1, 2));
        mixed.insert(mixed.end(), 17, Complex(-2, 1));
        mixed.insert(mixed.end(), 3, Complex(Rational(1)/3, Rational(-1)/5));
        known_roots(mixed);
        // Different factor exponents and roots closer than the initial output accuracy: a counted cluster under the
        // approximation goal, the combined-cover refinement under the accuracy goal.
        mixed.assign(4, Complex(Rational(1)/3));
        mixed.insert(mixed.end(), 3, Complex(Rational(1)/3 + power2(-160)));
        known_roots(mixed);
        known_roots(mixed, "0", CPS_ROOT_ACCURACY);
        // A bad reduction prime is inconclusive, not proof of repeated roots.
        known_roots({Complex(1), Complex(998244354)});
        // Irreducible quadratic repeated factor, checked against an independent
        // exact winding calculation on the original expanded polynomial.
        const ComplexPolynomial quadratic{{1}, {-1}, {1}};
        const auto expanded = multiply(quadratic, quadratic);
        auto input = reference::input(expanded);
        cps_options options;
        options.goal = CPS_ROOT_APPROXIMATION;
        d::PiecewiseRootStatistics statistics;
        auto result = d::solve_prepared_piecewise_roots(certpolysolve::detail::dense(*input.get()), options, statistics);
        require(result.complete && result.roots.size() == 2, "quadratic factor completion failed");
        for (const auto &root : result.roots)
            require(root.multiplicity == 2 && rectangle_count(expanded,
                endpoint(root.endpoints[0]), endpoint(root.endpoints[1]),
                endpoint(root.endpoints[2]), endpoint(root.endpoints[3])) == root.count,
                "quadratic formula count failed independent winding");
        auto low_precision = options;
        low_precision.max_work_bits = 32;
        const auto limited = d::solve_prepared_piecewise_roots(certpolysolve::detail::dense(*input.get()), low_precision, statistics);
        require(!limited.goal_met && limited.termination == CPS_PRECISION_LIMIT,
                "quadratic formula ignored insufficient precision");
        options.cancel = [](void *) { return 1; };
        result = d::solve_prepared_piecewise_roots(certpolysolve::detail::dense(*input.get()), options, statistics);
        require(!result.complete && result.termination == CPS_CANCELLED, "cancelled preparation fabricated completion");
        options.cancel = nullptr;
        options.work_units = 0;
        result = d::solve_prepared_piecewise_roots(certpolysolve::detail::dense(*input.get()), options, statistics);
        require(!result.goal_met && result.termination == CPS_WORK_LIMIT && result.progress.work_units_this_call == 0,
                "prepared factors ignored work allowance");
        options.exact_budget = 0;
        const auto unavailable = d::prepare_squarefree_input(certpolysolve::detail::dense(*input.get()), options);
        require(!unavailable.available(), "exhausted algebra allowance invented a decomposition");
        // A completed factor's box can contain an unresolved root of another
        // factor. Partial factor covers must not undercount the original input.
        std::vector<Complex> near{Complex(Rational(1)/3)};
        for (const auto &z : {Complex(Rational(1)/3 + power2(-160)), Complex(2), Complex(3)}) {
            near.push_back(z);
            near.push_back(z);
        }
        auto partial_input = reference::input(from_roots(near));
        options = cps_options{};
        options.goal = CPS_ROOT_APPROXIMATION;
        options.work_units = 1;
        result = d::solve_prepared_piecewise_roots(certpolysolve::detail::dense(*partial_input.get()), options, statistics);
        // One work unit certifies nothing: the stop returns the counted region
        // that holds every root (queue item 3.2), never an undercount.
        uint64_t stop_count = 0;
        for (const auto &root : result.roots) stop_count += root.count;
        require(result.complete && !result.goal_met && stop_count == near.size() && d::disjoint(result.roots),
                "partial factor cover understated original-input root count");
        // The saved tridiagonal degree-40 counterexample used to revisit the
        // same candidates at 181 bits until its external deadline. Reconstruct
        // it algebraically so this regression needs no benchmark artifacts.
        ComplexPolynomial stalled{{2}, {-5}, {5}, {-11}, {2}, {13}, {10}, {9}, {-30},
                                  {15}, {-9}, {21}, {-15}, {10}, {-6}, {3}, {-3}, {1}};
        stalled = multiply(multiply(stalled, quadratic), quadratic);
        std::vector<Complex> repeated(6, Complex(0));
        repeated.insert(repeated.end(), 6, Complex(1));
        repeated.insert(repeated.end(), 7, Complex(-1));
        stalled = multiply(stalled, from_roots(repeated));
        auto stalled_input = reference::input(stalled);
        options = cps_options{};
        options.goal = CPS_ROOT_APPROXIMATION;
        // The raw worker cannot certify the repeated roots at any precision;
        // the cap is where it must stop. 400 bits: the ladder's rungs for 53
        // bits are 149, 181, 245, 373, and the last is cut to the cap (queue
        // item 3.2); a 1024-bit cap would spend a minute on a 1024-bit rung.
        options.max_work_bits = 400;
        options.work_units = 10000;
        d::PiecewiseRootConfig config;
        config.reuse_guards = config.reuse_polish = true;
        config.conjugate_symmetry = config.recover_discovery = config.blocked_polish = true;
        result = d::solve_piecewise_roots(certpolysolve::detail::dense(*stalled_input.get()), options, statistics, config);
        if (result.goal_met || result.termination != CPS_PRECISION_LIMIT)
            std::cerr << "raw progress: termination=" << result.termination << " complete=" << result.complete
                      << " rounds=" << statistics.rounds << " model_bits=" << statistics.model_bits
                      << " bits=" << statistics.arithmetic_bits << " stationary=" << statistics.stationary_rounds
                      << " recovery=" << statistics.recovery_rounds << " tentative=" << statistics.unproved_seeds_retained
                      << " improvements=" << statistics.candidate_improvements << " rejected=" << statistics.seed_rejected
                      << " candidate_limited=" << statistics.candidate_limited << " model_limited=" << statistics.model_limited
                      << " rounding_limited=" << statistics.rounding_limited << '\n';
        require(!result.goal_met && result.termination == CPS_PRECISION_LIMIT,
                "raw repeated-root worker spent its work allowance cycling at unchanged precision");
        result = d::solve_prepared_piecewise_roots(certpolysolve::detail::dense(*stalled_input.get()), options, statistics, config);
        require(result.complete && result.goal_met && d::disjoint(result.roots),
                "prepared tridiagonal control failed to complete");
        std::cout << "prepared piecewise exact multiplicity, overlap, quadratic and budget checks passed\n";
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
