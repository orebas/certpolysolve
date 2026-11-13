#include "solve/policies.hpp"
#include "numeric/arithmetic/polynomial_arithmetic.hpp"
#include <certpolysolve/certpolysolve.hpp>
#include "address_space_cap.hpp"
#include "solve/auto_selection.hpp"
#include "numeric/evaluate/disk_sparse_polynomial_evaluator.hpp"
#include "exact_winding_api.hpp"
#include "core/heuristics.hpp"
#include "internal.hpp"
#include "methods/piecewise/piecewise_admission.hpp"
#include "methods/piecewise/piecewise_roots.hpp"
#include "core/wide_int.hpp"
#include <cmath>
#include <limits>
#include <numbers>
#include <iostream>
#include <string>

namespace d = certpolysolve::detail;
using namespace reference;

namespace {
// A flat envelope: |c_j| ~ 1/sqrt(j!), with deterministic signs. Its roots
// spread over a disk that grows like sqrt(degree), the regime of the
// flat/elliptic/hyperbolic families where piecewise is the fast method.
Polynomial flat(int degree) {
    Polynomial p(static_cast<size_t>(degree) + 1);
    for (int j = 0; j <= degree; ++j) {
        const auto exponent = static_cast<int>(std::llround(0.5 * std::lgamma(j + 1.0) / std::log(2.0)));
        const bool negative = ((j * 2654435761u) >> 7) & 1u;
        p[static_cast<size_t>(j)] = negative ? Rational(-power2(-exponent)) : power2(-exponent);
    }
    return p;
}

cps_options approximation() {
    cps_options options;
    options.goal = CPS_ROOT_APPROXIMATION;
    options.bits = 53;
    return options;
}
} // namespace

int main() {
    certpolysolve::ThreadScope lifetime;
    try {
        // Queries outside the explicit methods keep the D5 route.
        {
            const auto input = reference::input(flat(24));
            auto options = approximation();
            options.real_only = true;
            require(d::select_auto_route(certpolysolve::detail::dense(*input.get()), options).route == d::AutoRoute::D5,
                    "real-only query left the D5 route");
            // The cluster-radius goal is routed like an accuracy goal since
            // 1 October 2026; a real-only cluster query keeps D5.
            options = approximation();
            options.goal = CPS_CLUSTER_RADIUS;
            require(d::select_auto_route(certpolysolve::detail::dense(*input.get()), options).route == d::AutoRoute::Secular,
                    "cluster-radius query did not use the maintained recipes");
            options.real_only = true;
            require(d::select_auto_route(certpolysolve::detail::dense(*input.get()), options).route == d::AutoRoute::D5,
                    "real-only cluster-radius query left the D5 route");
            certpolysolve::Polynomial uncertain;
            uncertain.set(0, certpolysolve::ball(certpolysolve::Number("1"), certpolysolve::Number("1/1024")))
                .set(1, certpolysolve::ball(certpolysolve::Number("0"), certpolysolve::Number("0")))
                .set(2, certpolysolve::ball(certpolysolve::Number("1"), certpolysolve::Number("0")));
            require(uncertain.ok(), "the uncertain input recorded a bad value");
            require(d::select_auto_route(certpolysolve::detail::dense(*uncertain.get()), approximation()).route == d::AutoRoute::D5,
                    "uncertain coefficients left the D5 route");
        }
        // Moderate degree: the global methods' budget is cheaper than
        // piecewise's estimated model and local work.
        require(d::select_auto_route(certpolysolve::detail::dense(*reference::input(flat(256)).get()), approximation()).route ==
                    d::AutoRoute::Secular,
                "degree-256 flat input did not use secular");
        // A power composite q(x^2): piecewise refuses it, secular reduces it.
        {
            const auto base = flat(64);
            Polynomial composite(2 * base.size() - 1);
            for (size_t j = 0; j < base.size(); ++j) composite[2 * j] = base[j];
            require(d::select_auto_route(certpolysolve::detail::dense(*reference::input(composite).get()), approximation()).route ==
                        d::AutoRoute::Secular,
                    "power composite was routed to piecewise");
        }
        // Large flat degree: piecewise's estimate falls below the global
        // stage's minimum sweep budget.
        const auto large = d::select_auto_route(certpolysolve::detail::dense(*reference::input(flat(8192)).get()), approximation());
        require(large.route == d::AutoRoute::Piecewise, "degree-8192 flat input did not use piecewise");
        // A precision cap below piecewise's first working precision (its
        // ladder's base rung for the goal) is not routed to piecewise, which
        // would refuse it before any work (queue item 3.2).
        {
            auto capped = approximation();
            capped.max_work_bits = static_cast<uint64_t>(d::piecewise_start_precision(capped.bits)) - 1;
            const auto below = d::select_auto_route(certpolysolve::detail::dense(*reference::input(flat(8192)).get()), capped);
            require(below.route == d::AutoRoute::Secular &&
                        std::string(below.reason) == "precision-cap-below-piecewise-start",
                    "a precision cap below piecewise's start was routed to piecewise");
            capped.max_work_bits += 1;
            require(d::select_auto_route(certpolysolve::detail::dense(*reference::input(flat(8192)).get()), capped).route == d::AutoRoute::Piecewise,
                    "a precision cap at piecewise's start was refused");
        }
        // The piecewise estimate walks one ring per octave of the root window,
        // whatever the degree: about 1.6e9 rings for x^3 + 2^(2^30 + 6) x - 1.
        // The selector stops it once it passes the global sweep budget. (It
        // ran for minutes and then failed to allocate; Auto returned an error
        // on an input secular solves in three work units.)
        uint64_t wide_ns = 0;
        {
            const AddressSpaceCap cap(4ull << 30);
            certpolysolve::Polynomial wide;
            wide.set(0, certpolysolve::Number("-1")).set(1, certpolysolve::Number("1@1073741830")).set(3, certpolysolve::Number("1"));
            require(wide.ok(), "the wide input recorded a bad value");
            const auto options = approximation();
            const auto selection = d::select_auto_route(certpolysolve::detail::dense(*wide.get()), options);
            wide_ns = selection.elapsed_ns;
            require(selection.route == d::AutoRoute::Secular, "a wide root window was not routed to secular");
            d::PiecewiseAdmissionPolicy capped;
            capped.maximum_envelope_bend = capped.maximum_work_ratio = std::nullopt;
            capped.maximum_global_sweeps = d::heuristics::initial_sweep_allowance;
            const auto admission = d::admit_piecewise(certpolysolve::detail::dense(*wide.get()), options, capped);
            require(!admission.selected && std::string(admission.reason) == "estimate_above_global_sweeps",
                    "the estimate of a wide root window was not stopped at the sweep budget");
            certpolysolve::Options query;
            query.goal = certpolysolve::Goal::approximation();
            const auto solved = certpolysolve::solve(wide, query);
            require(solved.goal_met() && solved.roots().size() == 3, "Auto did not solve the wide input");
        }
        // The cap changes no decision: the capped admission selects exactly
        // when the full estimate satisfies the selector's inequality, on
        // inputs on both sides of it.
        for (const int degree : {16, 256, 1024, 8192}) {
            const auto input = reference::input(flat(degree));
            const auto options = approximation();
            d::PiecewiseAdmissionPolicy lifted;
            lifted.maximum_envelope_bend = lifted.maximum_work_ratio = std::nullopt;
            const auto full = d::admit_piecewise(certpolysolve::detail::dense(*input.get()), options, lifted);
            // Below the budget: the work under budget sweeps of degree x (degree + support) operations each.
            const auto budget = d::heuristics::initial_sweep_allowance;
            const bool below = full.selected && full.work < budget * full.degree * (full.degree + full.support);
            lifted.maximum_global_sweeps = budget;
            require(d::admit_piecewise(certpolysolve::detail::dense(*input.get()), options, lifted).selected == below,
                    "the sweep cap changed a piecewise selection");
            require((d::select_auto_route(certpolysolve::detail::dense(*input.get()), options).route == d::AutoRoute::Piecewise) == below,
                    "the selector no longer follows the sweep budget");
        }
        // The certificate attempt before piecewise. Its cost in global
        // sweeps: the binary64 packet's 32, and two evaluations of value and
        // slope at 17 operations a coefficient (heuristics.hpp): 66 for a
        // dense polynomial, 32 for a sparse one. Counted exactly at one point
        // (certificate_attempt_point_operations): those sweeps times degree + support.
        {
            // Degree 3,200 is past the interaction tree's crossover: the
            // packet's sweep is an evaluation of 3,201 terms and twelve levels
            // of the tree at each point.
            const uint64_t dense = d::heuristics::certificate_attempt_point_operations(3200, 3201);
            const uint64_t level = d::heuristics::machine_tree_level_weight;
            const uint64_t dense_expected = 32 * (3201 + level * 12) + 2 * 2 * 17 * 3201;
            // A sparse source's passes also build its terms' powers: nine
            // terms of degree 6,400 at exponents 0, 800, ..., 6,400 (gaps of
            // 800 = 0b1100100000, so the plan's products depend only on it).
            std::vector<slong> exponents;
            for (slong e = 0; e <= 6400; e += 800) exponents.push_back(e);
            const auto products = d::DiskSparsePolynomialEvaluator::plan(exponents).multiplications;
            const uint64_t sparse = d::heuristics::certificate_attempt_point_operations(6400, 9, products);
            // Its packet evaluates the nine terms at ten operations each.
            const uint64_t packet = 32 * (10 * 9 + level * 13);
            const uint64_t passes = 2 * (2 * 17 * 9 + 70 * products);
            require(dense == dense_expected && products > 0 && sparse == packet + passes,
                    "the attempt's cost is not the packet plus two evaluation passes");
            // Below the tree's crossover a dense packet sweep is one global sweep: an evaluation of 1,001 terms and
            // 999 interactions at a point.
            require(d::heuristics::machine_packet_point_operations(1000, 1001) == 32 * 2000,
                    "a direct packet sweep is not priced as a global sweep");
            // The order rule: the attempt first exactly when piecewise's
            // estimate is at least the golden ratio times that cost, on
            // inputs on both sides of it.
            bool first = false, later = false;
            for (const int degree : {2048, 4096, 8192, 16384}) {
                const auto input = reference::input(flat(degree));
                const auto options = approximation();
                d::PiecewiseAdmissionPolicy lifted;
                lifted.maximum_envelope_bend = lifted.maximum_work_ratio = std::nullopt;
                lifted.maximum_global_sweeps = d::heuristics::initial_sweep_allowance;
                const auto admission = d::admit_piecewise(certpolysolve::detail::dense(*input.get()), options, lifted);
                const auto selection = d::select_auto_route(certpolysolve::detail::dense(*input.get()), options);
                require(admission.selected == (selection.route == d::AutoRoute::Piecewise),
                        "the selector does not follow piecewise's admission");
                if (!admission.selected) {
                    require(!selection.attempt_first, "an attempt before piecewise without piecewise");
                    continue;
                }
                // The rule restated exactly: in global sweeps the estimate W is work / (n (n + s)) and the attempt's cost
                // A is its operations at a point over n + s, so W >= phi A is x >= phi for x = work / (A' n), A' those
                // operations; and for x > 0 that is x^2 >= x + 1 (phi is the positive root of x^2 - x - 1):
                // work^2 >= work A' n + (A' n)^2, in 128 bits.
                const d::uint128 work = admission.work;
                const d::uint128 attempt =
                    d::uint128{d::heuristics::certificate_attempt_point_operations(admission.degree, admission.support)} *
                    admission.degree;
                const bool expected = work * work >= work * attempt + attempt * attempt;
                require(selection.attempt_first == expected, "the selector does not follow the order rule");
                // The attempt before piecewise may spend piecewise's estimate (machine_newton_cover's budget).
                require(selection.piecewise_operations == admission.work,
                        "the selection does not carry piecewise's estimate in operations");
                require(std::string(selection.reason) == (expected ? "certificate-attempt-before-piecewise"
                                                                   : "piecewise-estimate-below-global-sweep-budget"),
                        "the selection's reason does not name the order");
                (expected ? first : later) = true;
            }
            require(first && later, "the order rule's fixtures are all on one side of it");
        }
        // A goal below the default is routed as the default goal is: the
        // selector's estimate drops piecewise's upper model rungs as the goal
        // falls, while neither route's cost falls with it (auto_selection.cpp).
        // At least one fixture is one whose estimate at the lower goal would
        // have passed under the global sweep budget.
        {
            bool would_flip = false;
            for (const int degree : {256, 1024, 2048, 4096}) {
                const auto input = reference::input(flat(degree));
                const auto at_default = d::select_auto_route(certpolysolve::detail::dense(*input.get()), approximation());
                for (const uint64_t bits : {10, 24}) {
                    auto options = approximation();
                    options.bits = bits;
                    const auto lower = d::select_auto_route(certpolysolve::detail::dense(*input.get()), options);
                    require(lower.route == at_default.route && lower.attempt_first == at_default.attempt_first,
                            "a goal below the default was routed differently from the default goal");
                    d::PiecewiseAdmissionPolicy lifted;
                    lifted.maximum_envelope_bend = lifted.maximum_work_ratio = std::nullopt;
                    lifted.maximum_global_sweeps = d::heuristics::initial_sweep_allowance;
                    would_flip = would_flip || (d::admit_piecewise(certpolysolve::detail::dense(*input.get()), options, lifted).selected &&
                                                at_default.route != d::AutoRoute::Piecewise);
                }
            }
            require(would_flip, "no fixture's lower-goal estimate passes the budget where the default's does not");
        }
        // End to end through the session: Auto finishes a dense polynomial of
        // degree 2,400 with coefficients 1, 2, ..., 2,401 by the attempt (it
        // would be piecewise's), and says so.
        {
            Polynomial easy(2401);
            for (size_t j = 0; j < easy.size(); ++j) easy[j] = Rational(static_cast<long long>(j) + 1);
            const auto input = reference::input(easy);
            auto options = approximation();
            const auto selection = d::select_auto_route(certpolysolve::detail::dense(*input.get()), options);
            require(selection.route == d::AutoRoute::Piecewise && selection.attempt_first,
                    "the easy fixture is not routed to the attempt before piecewise: it tests nothing");
            d::SessionState session(certpolysolve::detail::dense(*input.get()));
            const auto result = session.advance(&options);
            require(result.complete && result.goal_met && d::disjoint(result.roots),
                    "Auto did not finish the easy fixture");
            require(result.diagnostics.rfind("auto-route=piecewise (certificate-attempt-before-piecewise); "
                                             "certificate attempt before piecewise: complete cover", 0) == 0,
                    "the attempt before piecewise is missing from the diagnostics");
            uint64_t total = 0;
            for (const auto &root : result.roots) total += root.count;
            require(total == 2400 && result.roots.size() == 2400 && result.statistics_v3.validation_calls == 0,
                    "the attempt before piecewise miscounted the roots");
        }
        // End to end through the session: a complex-rotated input that D5
        // handles poorly completes on the secular route, with the route in
        // the diagnostics and every planted root accounted for.
        {
            std::vector<Complex> roots;
            const Complex unit(Rational(3, 5), Rational(4, 5));
            for (int j = 1; j <= 12; ++j)
                roots.push_back(unit * Complex(Rational(j, 4), Rational(j * j % 7 - 3, 8)));
            const auto input = reference::input(from_roots(roots));
            d::SessionState session(certpolysolve::detail::dense(*input.get()));
            auto options = approximation();
            const auto result = session.advance(&options);
            require(result.complete && result.goal_met && d::disjoint(result.roots), "selector route did not complete");
            require(result.diagnostics.rfind("auto-route=secular", 0) == 0, "route missing from diagnostics");
            uint64_t total = 0;
            for (const auto &root : result.roots) total += root.count;
            require(total == roots.size(), "selector route miscounted the roots");
            // The public entry routes through the selector too.
            d::SessionState shipped(certpolysolve::detail::dense(*input.get()));
            const auto published = shipped.advance(&options);
            require(published.diagnostics.find("auto-route=") != std::string::npos,
                    "public Auto did not take the selector's route");
        }
        // The D5 recipe's local cluster certificates (production_policy(), the recipe Auto's selector sends inexact
        // coefficients and unsupported queries to: session.cpp): repeated roots with exact preparation disabled end in
        // three Pellet-counted clusters. Driven directly, as the session drives it (until W8's S2 through a session
        // that skipped the selector).
        {
            const auto input = reference::input(from_roots({Complex(1), Complex(1), Complex(-2), Complex(-2), Complex(4)}));
            d::SolveState d5(certpolysolve::detail::dense(*input.get()), d::maintained_method_policy(CPS_METHOD_AUTO));
            cps_options options;
            options.bits = 128;
            options.max_work_bits = 128;
            options.exact_budget = 0;
            const auto result = d5.advance(&options);
            require(result.complete && !result.goal_met && result.roots.size() == 3,
                    "D5 did not recover three local clusters");
            for (const auto &root : result.roots)
                require(root.certificate.rfind("pellet", 0) == 0, "D5 local cluster lacks a Pellet count certificate");
        }
        std::cout << "auto selection: D5 for unsupported queries, secular below the piecewise estimate, piecewise at "
                     "degree 8192 ("
                  << large.elapsed_ns / 1000 << " us to decide); a root window of 1.6e9 octaves decided in "
                  << wide_ns / 1000 << " us\n";
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
