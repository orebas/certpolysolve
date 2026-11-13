#include "numeric/arithmetic/polynomial_arithmetic.hpp"
#include "exact_winding_api.hpp"
#include "methods/piecewise/piecewise_roots.hpp"
#include <iostream>
namespace d = certpolysolve::detail;
using reference::require;
int
main() {
    certpolysolve::ThreadScope lifetime;
    try {
        auto input =
          reference::input(reference::from_roots({ { 1, 2 }, { -3, 1 }, { 2, -1 }, { 5 }, { 0, 3 }, { -2, -3 } }));
        cps_options options;
        options.goal = CPS_ROOT_APPROXIMATION;
        options.max_work_bits = 2048;
        options.work_units = UINT64_MAX;
        auto config = d::maintained_piecewise_config();
        // Diagnostics must also work with source finishing disabled. Compare
        // the same explicit recipe in the control, observed and lazy runs.
        config.source_finisher = false;
        config.reuse_polish_bound = false;
        d::PiecewiseRootStatistics baseline_stats;
        const auto baseline = d::solve_piecewise_roots(certpolysolve::detail::dense(*input.get()), options, baseline_stats, config);
        require(baseline.goal_met, "diagnostic control did not complete");
        std::vector<d::PiecewiseRootSnapshot> snapshots;
        std::vector<d::PiecewiseTransformCapture> captures;
        size_t polls = 0;
        uint64_t observed_evaluations = 0, observed_guards = 0;
        std::vector<d::PiecewisePopulationSnapshot> populations;
        size_t selection_count = 0, attempt_count = 0;
        config.observe_selection = [&](const auto &e) {
            require(e.candidate != UINT64_MAX, "selection lost its candidate identity");
            require(std::isfinite(e.radius_log2), "selection lost chart scale");
            ++selection_count;
        };
        config.observe_candidate = [&](const auto &e) {
            require(e.raw_failure == e.failure, "unrequested failure reclassification");
            ++attempt_count;
        };
        config.observe_population = [&](auto s) { populations.push_back(std::move(s)); };
        config.observe_polish = [&](const auto &e) { observed_evaluations += e.evaluations; };
        config.observe_guard = [&](const auto &, const auto &, const auto &, slong, bool) { ++observed_guards; };
        config.snapshot_discovery = config.snapshot_rounds = true;
        config.snapshot_due = [&] { return ++polls % 32 == 0; };
        config.observe_snapshot = [&](auto s) { snapshots.push_back(std::move(s)); };
        config.want_transform = [&](const auto &) { return captures.empty(); };
        config.capture_transform = [&](auto c) { captures.push_back(std::move(c)); };
        cps_result result;
        {
            d::PiecewiseRootState state(certpolysolve::detail::dense(*input.get()), config);
            auto partial = options;
            partial.work_units = 100;
            state.advance(partial);
            result = state.advance(options);
            const auto &stats = state.statistics();
            require(stats.proposals == baseline_stats.proposals && stats.models == baseline_stats.models &&
                      stats.polish_evaluations == baseline_stats.polish_evaluations,
                    "diagnostics changed numerical work");
            require(stats.finisher_attempts == 0, "diagnostics accidentally enabled source finishing");
        }
        require(result.goal_met && result.roots.size() == baseline.roots.size(), "diagnostic cover changed");
        require(observed_evaluations == baseline_stats.polish_evaluations && observed_guards > 0,
                "per-candidate observation missed numerical work");
        require(attempt_count == baseline_stats.proposals,
                "population observation missed numerical work");
        require(!populations.empty() && populations.back().terminal && populations.back().accurate == 6,
                "population observation missed completion");
        for (const auto &p : populations.back().points)
            require(acb_is_finite(p.point.value()), "population point lost its owned lifetime");
        std::vector<d::ComplexBall> supplied_points;
        for (const auto &r : baseline.roots) {
            d::ComplexBall point;
            acb_get_mid(point.value(), r.box.value());
            supplied_points.push_back(std::move(point));
        }
        d::PiecewiseRootStatistics supplied_stats;
        const auto supplied_result = d::certify_piecewise_points(certpolysolve::detail::dense(*input.get()), options,
            supplied_stats, supplied_points, config);
        require(supplied_result.goal_met && selection_count > 0,
                "population observer missed carried-point selection");
        for (size_t i = 0; i < result.roots.size(); ++i)
            require(acb_equal(result.roots[i].box.value(), baseline.roots[i].box.value()),
                    "diagnostics changed an output rectangle");
        require(!snapshots.empty() && snapshots.back().terminal && snapshots.back().unfinished.empty(),
                "completion snapshot retains finished owners");
        require(snapshots.back().statistics.final_accurate_count == 6, "incorrect snapshot progress");
        require(!captures.empty() && captures[0].input.size() == captures[0].event.length,
                "transform capture is absent or malformed");
        // Both the source and captured balls remain usable after the state and
        // all of its coroutine/model scratch have been destroyed.
        require(snapshots.front().source->coefficients.size() == certpolysolve::detail::dense(*input.get()).coefficients.size(),
                "snapshot source lifetime lost");
        for (const auto &z : captures[0].output) require(acb_is_finite(z.value()), "dangling FFT output");
        auto lazy = d::maintained_piecewise_config();
        lazy.source_finisher = false;
        lazy.reuse_polish_bound = true;
        d::PiecewiseRootStatistics lazy_stats;
        const auto reused = d::solve_piecewise_roots(certpolysolve::detail::dense(*input.get()), options, lazy_stats, lazy);
        require(reused.goal_met && reused.roots.size() == baseline.roots.size(), "lazy bound lost completion");
        require(lazy_stats.polish_bounds_skipped > 0 &&
                  lazy_stats.polish_evaluations == baseline_stats.polish_evaluations &&
                  lazy_stats.polish_model_stops == baseline_stats.polish_model_stops &&
                  lazy_stats.polish_rounding_stops == baseline_stats.polish_rounding_stops &&
                  lazy_stats.polish_step_stops == baseline_stats.polish_step_stops,
                "lazy bound failed to preserve proposal work");
        for (size_t i = 0; i < reused.roots.size(); ++i)
            require(acb_equal(reused.roots[i].box.value(), baseline.roots[i].box.value()),
                    "lazy bound changed a rectangle");
        std::cout << "owned diagnostic equivalence and continuation passed\n";
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
