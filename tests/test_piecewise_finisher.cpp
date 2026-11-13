#include "numeric/arithmetic/polynomial_arithmetic.hpp"
#include "exact_winding_api.hpp"
#include "methods/piecewise/piecewise_roots.hpp"
#include <iostream>

namespace d = certpolysolve::detail;
using namespace reference;
static void verify(const cps_result &result, const std::vector<Complex> &known, unsigned bits = 0) {
    require(d::disjoint(result.roots), "finisher published overlapping balls");
    std::vector<reference::Rectangle> exported;
    uint64_t total = 0;
    for (const auto &root : result.roots) {
        const auto endpoint = [](const d::Endpoint &e) {
            return reference::endpoint(certpolysolve::Dyadic{e.mantissa_text(), e.exponent_text()});
        };
        const reference::Rectangle rectangle{endpoint(root.endpoints[0]), endpoint(root.endpoints[1]),
                                             endpoint(root.endpoints[2]), endpoint(root.endpoints[3]), root.count};
        exported.push_back(rectangle);
        uint64_t count = 0;
        const Complex *first = nullptr;
        bool one_location = true;
        for (const auto &z : known)
            if (rectangle.real_lower <= z.real && z.real <= rectangle.real_upper &&
                rectangle.imag_lower <= z.imag && z.imag <= rectangle.imag_upper) {
                ++count;
                if (!first) first = &z;
                else one_location = one_location && first->real == z.real && first->imag == z.imag;
                if (root.realness == CPS_REAL) require(z.imag == 0, "finisher labelled a nonreal root real");
                if (root.realness == CPS_NONREAL) require(z.imag != 0, "finisher labelled a real root nonreal");
                if (bits) {
                    const Rational re = (rectangle.real_lower + rectangle.real_upper) / 2 - z.real;
                    const Rational im = (rectangle.imag_lower + rectangle.imag_upper) / 2 - z.imag;
                    const Rational factor = Rational(Integer(1) << (2 * bits));
                    require((re * re + im * im) * factor <= z.real * z.real + z.imag * z.imag,
                            "exported midpoint missed exact relative accuracy");
                }
            }
        require(count == root.count, "original-source finisher lost root identity");
        if (root.multiplicity)
            require(root.multiplicity == count && one_location, "incorrect isolated-root multiplicity");
        total += count;
    }
    // This check uses Boost rational endpoint comparisons, independently of
    // the production ball/index predicates and the finisher's admission helper.
    reference::check_disjoint(exported);
    if (result.complete) require(total == known.size(), "finisher declared incomplete coverage complete");
}
int main(int argc,char **argv) {
    certpolysolve::ThreadScope lifetime;
    try {
        std::vector<Complex> known;
        for (int i = 1; i <= 64; ++i) known.emplace_back(i);
        reference::Polynomial coefficients{1};
        for (const auto &z : known) coefficients = reference::multiply(coefficients, {-z.real, 1});
        const auto input = reference::input(coefficients);
        d::PiecewiseRootConfig config;
        config.reuse_guards = config.reuse_polish = config.blocked_polish = true;
        config.conjugate_symmetry = config.recover_discovery = config.durable_ownership = true;
        config.local_band_error = config.aligned_polish = config.source_finisher = true;
        for (int i = 1; i < argc; ++i) {
            const std::string argument = argv[i];
            if (argument == "owned-only") config.owned_only_after_discovery = true;
            else throw std::invalid_argument("unknown finisher test variant");
        }
        d::PiecewiseRootState state(certpolysolve::detail::dense(*input.get()), config);
        cps_options options;
        options.goal = CPS_ROOT_APPROXIMATION;
        options.max_work_bits = 2048;
        options.work_units = 1;
        cps_result partial;
        // Stop immediately after admission of a resumable finisher, before
        // allowing its Newton/certificate step to consume the next work unit.
        do {
            partial = state.advance(options);
        } while (!partial.goal_met && !state.statistics().finisher_attempts &&
                 partial.progress.work_units_this_call);
        verify(partial, known);
        const auto before = state.certified_regions();
        options.cancel = [](void *) { return 1; };
        const auto cancelled = state.advance(options);
        require(cancelled.roots.size() == before, "cancellation discarded finisher proofs");
        options.cancel = nullptr;
        // Resume only through first source work, then cancel with prepared
        // caches retained. Cancelled calls must neither rebuild nor evaluate.
        do {
            partial = state.advance(options);
        } while (!partial.goal_met && !state.statistics().finisher_evaluations &&
                 partial.progress.work_units_this_call);
        const auto prepared = state.statistics().finisher_preparation_ns;
        const auto evaluations = state.statistics().finisher_evaluations;
        require(prepared && evaluations, "control missed prepared source state");
        const auto certified = state.certified_regions();
        options.cancel = [](void *) { return 1; };
        const auto prepared_cancelled = state.advance(options);
        verify(prepared_cancelled, known);
        require(prepared_cancelled.roots.size() == certified &&
                state.statistics().finisher_preparation_ns == prepared &&
                state.statistics().finisher_evaluations == evaluations,
                "cancelled prepared finisher changed proofs or performed work");
        options.cancel = nullptr;
        options.work_units = 100000;
        const auto result = state.advance(options);
        verify(result, known, 53);
        require(result.goal_met, "original-source finisher control did not finish");
        require(state.statistics().finisher_attempts > 0, "test missed the finisher admission path");
        require(state.statistics().finisher_certificates > 0, "test missed original-source publication");
        require(state.statistics().finisher_block_evaluations == state.statistics().finisher_evaluations,
                "test bypassed block/disk source evaluation");
        // A tighter query must not inherit success or a stale cost calibration
        // from the previous goal. Verify the result against exact known roots.
        auto tighter = options;
        tighter.bits = 69;
        const auto refined = state.advance(tighter);
        verify(refined, known, 69);
        require(refined.goal_met, "tighter finisher query did not complete");
        {
            // The outer state owns an exact double zero while its finisher
            // certifies the other roots against the full original polynomial.
            auto with_zero = known;
            with_zero.emplace_back(0);
            with_zero.emplace_back(0);
            const auto zero_input = reference::input(from_roots(with_zero));
            d::PiecewiseRootStatistics zero_statistics;
            const auto zero_result = d::solve_piecewise_roots(d::dense(*zero_input.get()), options,
                                                             zero_statistics, config);
            verify(zero_result, with_zero, 53);
            require(zero_result.goal_met && zero_statistics.finisher_certificates > 0,
                    "repeated-zero control missed a completed source-finisher solve");
            bool exact_zero = false;
            for (const auto &root : zero_result.roots)
                if (root.count == 2 && root.multiplicity == 2 && root.realness == CPS_REAL) {
                    exact_zero = true;
                    for (const auto &end : root.endpoints)
                        require(reference::endpoint(certpolysolve::Dyadic{end.mantissa_text(), end.exponent_text()}) == 0,
                                "the retained exact zero acquired a nonzero endpoint");
                }
            require(exact_zero, "finishing lost the exact zero's count, multiplicity or realness");
        }
        {
            // Independent exact complex roots exercise the original-source
            // certificate after rotation and scaling, not just real intervals.
            std::vector<Complex> rotated;
            const Complex scale(Rational(3) / 4, Rational(1) / 2);
            for (const auto &z : known) rotated.push_back(z * scale + Complex(1, -2));
            auto complex_coefficients = from_roots(rotated);
            // The roots stay exact, but a nonreal nonunit leading coefficient
            // exercises original versus monic source preparation as well.
            const Complex coefficient_scale(Rational(3)/8,Rational(5)/16);
            for (auto &coefficient : complex_coefficients) coefficient = coefficient * coefficient_scale;
            const auto complex_input = reference::input(complex_coefficients);
            d::PiecewiseRootStatistics complex_statistics;
            const auto complex_result = d::solve_piecewise_roots(certpolysolve::detail::dense(*complex_input.get()), options,
                                                                 complex_statistics, config);
            verify(complex_result, rotated, 53);
            require(complex_result.goal_met, "rotated complex finisher did not complete");
            require(complex_statistics.finisher_certificates > 0,
                    "complex root test missed source-finisher publication");
            require(complex_statistics.finisher_block_evaluations == complex_statistics.finisher_evaluations,
                    "complex source finisher bypassed block/disk evaluation");
            for (const auto &root : complex_result.roots)
                require(d::meets_goal(root, options, static_cast<slong>(complex_result.work_bits)),
                        "rotated complex finisher lost requested accuracy");
        }
        {
            // Real coefficients with conjugate pairs: under conjugate symmetry
            // the lower half-plane's models are skipped, so every publication
            // path, the original-source finisher included, must publish the
            // mirror in the same step. Check closure after every work unit.
            std::vector<Complex> pairs;
            for (int k = 1; k <= 32; ++k) {
                pairs.emplace_back(Rational(k), Rational(k, 2));
                pairs.emplace_back(Rational(k), Rational(-k, 2));
            }
            const auto paired = reference::input(from_roots(pairs));
            d::PiecewiseRootState mirrored(certpolysolve::detail::dense(*paired.get()), config);
            auto step = options;
            step.work_units = 1;
            const auto closed = [](const cps_result &cover) {
                auto endpoint = [](const d::Endpoint &e) {
                    return reference::endpoint(certpolysolve::Dyadic{e.mantissa_text(), e.exponent_text()});
                };
                for (const auto &root : cover.roots) {
                    const auto lo = endpoint(root.endpoints[2]), hi = endpoint(root.endpoints[3]);
                    if (lo <= Rational(0) && Rational(0) <= hi) continue;
                    bool found = false;
                    for (const auto &other : cover.roots)
                        found = found || (endpoint(other.endpoints[0]) <= endpoint(root.endpoints[1]) &&
                                          endpoint(root.endpoints[0]) <= endpoint(other.endpoints[1]) &&
                                          endpoint(other.endpoints[2]) <= -lo && -hi <= endpoint(other.endpoints[3]));
                    if (!found) return false;
                }
                return true;
            };
            cps_result current;
            for (int calls = 0; calls < 100000 && !current.goal_met; ++calls) {
                current = mirrored.advance(step);
                require(closed(current), "a published nonreal root lacks its conjugate");
                if (!current.progress.work_units_this_call) break;
            }
            verify(current, pairs, 53);
            require(current.goal_met, "conjugate-pair finisher control did not finish");
        }
        std::cout << "finisher attempts=" << state.statistics().finisher_attempts
                  << " certificates=" << state.statistics().finisher_certificates
                  << " evaluations=" << state.statistics().finisher_evaluations << '\n';
    } catch (const std::exception &e) { std::cerr << e.what() << '\n'; return 1; }
}
