#include "numeric/evaluate/structured_source.hpp"
#include "solve/policies.hpp"
#include "solve/local_model_solver.hpp"
#include "numeric/arithmetic/polynomial_arithmetic.hpp"
#include "internal.hpp"
#include <certpolysolve/certpolysolve.hpp>
#include "exact_winding_api.hpp"
#include "api_test_support.hpp"
#include "methods/secular/secular_refinement.hpp"
#include "numeric/exact/centering.hpp"
#include "numeric/exact/power_composition.hpp"
#include <iostream>

using namespace reference;
namespace d = certpolysolve::detail;

static void enclose(acb_t out, const Complex &value, slong bits) {
    d::Number(value.real.str()).enclose(acb_realref(out), bits);
    d::Number(value.imag.str()).enclose(acb_imagref(out), bits);
}

static bool contains(const acb_t ball, const Complex &value) {
    d::Rational real, imag;
    fmpq_set_str(real.value(), value.real.str().c_str(), 10);
    fmpq_set_str(imag.value(), value.imag.str().c_str(), 10);
    return arb_contains_fmpq(acb_realref(ball), real.value()) &&
           arb_contains_fmpq(acb_imagref(ball), imag.value());
}

static Complex quotient(const Complex &a, const Complex &b) {
    const auto norm = b.real * b.real + b.imag * b.imag;
    return {(a.real * b.real + a.imag * b.imag) / norm,
            (a.imag * b.real - a.real * b.imag) / norm};
}

static std::unique_ptr<d::PolynomialSource> exact_source(const ComplexPolynomial &p) {
    d::GaussianPolynomial exact(p.size());
    for (size_t j = 0; j < p.size(); ++j) {
        fmpq_set_str(exact[j].real.value(), p[j].real.str().c_str(), 10);
        fmpq_set_str(exact[j].imag.value(), p[j].imag.str().c_str(), 10);
    }
    d::PolynomialSourceOptions options;
    options.input = d::InputSemantics::Exact;
    options.candidate_direct = true;
    options.exact_taylor = [exact](acb_poly_t out, const acb_t center, slong bits) {
        d::GaussianRational shift;
        arf_get_fmpq(shift.real.value(), arb_midref(acb_realref(center)));
        arf_get_fmpq(shift.imag.value(), arb_midref(acb_imagref(center)));
        d::GaussianPolynomial translated;
        if (!d::center_gaussian(translated, exact, shift, 8000000)) return false;
        d::enclose_gaussian(out, translated, bits);
        return true;
    };
    return std::make_unique<d::PolynomialSource>(
        [exact](acb_poly_t out, slong bits) { d::enclose_gaussian(out, exact, bits); },
        static_cast<slong>(p.size() - 1), false, options);
}

static std::vector<Complex> tight_complex_roots() {
    std::vector<Complex> roots;
    for (int j = 0; j < 4; ++j) {
        // These centers are deliberately not on the exact-model cache grid.
        // Resolving their much smaller internal spacing requires cancellation
        // beyond the coordinate precision used for interpolation updates.
        roots.emplace_back(Rational(2) + power2(-73) + j * power2(-180),
                           Rational(1) + power2(-79) + j * power2(-179));
        roots.emplace_back(Rational(-3) + power2(-83) + j * power2(-170),
                           Rational(2) + power2(-81) - j * power2(-172));
        roots.emplace_back(j + 5, j + 3);
    }
    return roots;
}

static void rigorous_regeneration() {
    const auto roots = tight_complex_roots();
    const auto polynomial = from_roots(roots);
    auto source = exact_source(polynomial);
    std::vector<Complex> exact_nodes;
    std::vector<d::ComplexBall> nodes(roots.size()), residues;
    for (size_t j = 0; j < roots.size(); ++j) {
        exact_nodes.push_back(roots[j] + Complex(power2(-190), power2(-192)));
        enclose(nodes[j].value(), exact_nodes.back(), 256);
        require(acb_is_exact(nodes[j].value()), "test interpolation node was rounded");
    }
    std::vector<slong> iteration_bits(roots.size(), 80);
    cps_options options;
    options.bits = 53;
    options.max_work_bits = 4096;
    for (bool prepared : {false, true}) {
        d::EnginePolicy policy;
        policy.prepared_source_refinement = prepared;
        certpolysolve::detail::CandidateStatistics statistics{};
        uint64_t maximum = 0;
        require(d::regenerate_residues(residues, *source, nodes, {}, {}, iteration_bits,
                    options, statistics, maximum, &policy), "exact-model regeneration failed");
        for (size_t j = 0; j < roots.size(); ++j) {
            Complex product = polynomial.back();
            for (size_t k = 0; k < roots.size(); ++k)
                if (j != k) product = product * (exact_nodes[j] - exact_nodes[k]);
            const auto expected = quotient(Complex(0) - evaluate(polynomial, exact_nodes[j]), product);
            require(contains(residues[j].value(), expected), "regenerated residue lost its exact identity");
            if (prepared)
                require(acb_rel_accuracy_bits(residues[j].value()) >= 69,
                        "interpolation model accepted an unresolved nonzero residue");
        }
        if (prepared)
            require(maximum > 112, "conditioning test did not require separate coefficient precision");
    }
}

// The answer of a solve state (the internal result), as the 0.5 header's Result.
static certpolysolve::Result wrap(cps_result &&result) {
    return certpolysolve_test::published(std::move(result));
}

static void known_cover(const certpolysolve::Result &result, const std::vector<Complex> &roots,
                        int bits, cps_goal goal) {
    require(result.goal_met() && result.complete() && cps_result_goal(result.get()) == goal,
            "prepared-source solve did not meet its contract");
    const auto regions = result.roots();
    check_disjoint(regions);
    size_t total = 0;
    for (const auto &region : regions) {
        const auto lo = endpoint(region.box.real.lo), hi = endpoint(region.box.real.hi);
        const auto ilo = endpoint(region.box.imag.lo), ihi = endpoint(region.box.imag.hi);
        const auto count = std::count_if(roots.begin(), roots.end(), [&](const auto &root) {
            return lo <= root.real && root.real <= hi && ilo <= root.imag && root.imag <= ihi;
        });
        require(static_cast<uint64_t>(count) == region.count, "known-root region count disagrees");
        total += static_cast<size_t>(count);
        if (goal == CPS_ROOT_ACCURACY) {
            require(region.count == 1 && region.multiplicity.value_or(0) == 1, "individual solve retained a cluster");
            check_accuracy(region, bits);
        } else {
            const Rational re = lo > 0 ? lo : hi < 0 ? -hi : Rational(0);
            const Rational im = ilo > 0 ? ilo : ihi < 0 ? -ihi : Rational(0);
            require(((hi - lo) * (hi - lo) + (ihi - ilo) * (ihi - ilo)) * power2(2 * bits) <=
                        4 * (re * re + im * im), "counted region failed exact relative accuracy");
            require(region.multiplicity.value_or(0) == 0 || (region.multiplicity.value_or(0) == 1 && count == 1),
                    "cluster count became an invented multiplicity");
        }
    }
    require(total == roots.size(), "prepared-source cover omitted a root");
}

static void known_complex_solve() {
    const auto roots = tight_complex_roots();
    const auto polynomial = input(from_roots(roots));
    for (auto goal : {CPS_ROOT_ACCURACY, CPS_ROOT_APPROXIMATION}) {
        auto policy = d::production_policy();
        policy.prepared_source_refinement = true;
        unsigned regenerations = 0;
        policy.trace_data = &regenerations;
        policy.trace = [](void *data, const d::TraceEvent &event) {
            if (std::string_view(event.stage) == "regeneration" && std::string_view(event.phase) == "begin")
                ++*static_cast<unsigned *>(data);
        };
        d::SolveState state(certpolysolve::detail::dense(*polynomial.get()), policy);
        cps_options options;
        options.bits = 53;
        options.max_work_bits = 8192;
        options.goal = goal;
        auto result = state.advance(&options);
        if (!result.goal_met) std::cerr << result.diagnostics << '\n';
        known_cover(wrap(std::move(result)), roots, 53, goal);
        if (goal == CPS_ROOT_ACCURACY)
            require(regenerations > 0, "known complex test did not enter prepared refinement");
        // 160 bits: a fresh solve (one advance a solve state since 18 October 2026).
        options.bits = 160;
        d::SolveState tighter(certpolysolve::detail::dense(*polynomial.get()), policy);
        result = tighter.advance(&options);
        known_cover(wrap(std::move(result)), roots, 160, goal);
    }
}

static void sparse_source_counts() {
    // Exercise the sparse-update source directly, so exact power reduction in
    // the outer solver cannot bypass it. The roots of (x^4-1)(x^4-16) are exact
    // Gaussian integers, giving independent counts at negligible test cost.
    const std::vector<Complex> roots{{1}, {-1}, {0, 1}, {0, -1},
                                     {2}, {-2}, {0, 2}, {0, -2}};
    const auto exact = from_roots(roots);
    auto encloser = [exact](acb_poly_t out, slong bits) {
        d::ComplexBall value;
        for (size_t j = 0; j < exact.size(); ++j) {
            enclose(value.value(), exact[j], bits);
            acb_poly_set_coeff_acb(out, static_cast<slong>(j), value.value());
        }
    };
    d::ComplexPolynomial coefficients;
    encloser(coefficients.value(), 128);
    d::RootVector seeds(8);
    require(d::initial_candidates(seeds.data(), coefficients.value(), 128), "sparse seed construction failed");
    d::PolynomialSourceOptions source_options;
    source_options.input = d::InputSemantics::Exact;
    source_options.sparse_coefficients = true;
    for (bool prepared : {false, true}) {
        auto policy = d::production_policy();
        policy.prepared_source_refinement = prepared;
        auto source = std::make_unique<d::PolynomialSource>(encloser, 8, true, source_options);
        d::SecularRefinement refinement(std::move(source), seeds.data(), 1, 128, policy, d::solve_local_model);
        cps_options options;
        options.bits = 53;
        options.max_work_bits = 2048;
        cps_result result;
        for (unsigned steps = 0; steps < 512 && refinement.active(); ++steps)
            require(refinement.step(options, result.statistics_v2, result.statistics_v3) != CPS_STOP_PRECISION,
                    "sparse refinement exhausted precision");
        require(!refinement.active() && refinement.cover(result.roots), "sparse refinement did not complete");
        result.complete = true;
        result.goal_met = std::all_of(result.roots.begin(), result.roots.end(), [&](const auto &root) {
            return d::meets_goal(root, options, static_cast<slong>(refinement.precision()));
        });
        known_cover(wrap(std::move(result)), roots, 53, CPS_ROOT_ACCURACY);
    }
}

static void composition_without_denominator_expansion() {
    const std::vector<Complex> roots{{1}, {-1}, {0, 1}, {0, -1},
                                     {2}, {-2}, {0, 2}, {0, -2}};
    certpolysolve::Polynomial polynomial;
    polynomial.set(0, certpolysolve::Number("16@-1000000")).set(4, certpolysolve::Number("-17@-1000000")).set(8, certpolysolve::Number("1@-1000000"));
    require(polynomial.ok(), polynomial.status().c_str());
    require(d::coefficient_power(certpolysolve::detail::dense(*polynomial.get()).coefficients) == 4,
            "stored exact coefficients lost power identity");
    for (auto goal : {CPS_ROOT_ACCURACY, CPS_ROOT_APPROXIMATION}) {
        auto policy = d::production_policy();
        policy.early_power_composition = true;
        // One solve with the 512 units the one-unit calls of a session used to share (one advance a solve state
        // since 18 October 2026).
        d::SolveState state(certpolysolve::detail::dense(*polynomial.get()), policy);
        cps_options options;
        options.exact_budget = 64;
        options.bits = 53;
        options.max_work_bits = 2048;
        options.goal = goal;
        options.work_units = 512;
        cps_result result = state.advance(&options);
        // The write-out is refused (no exact factor), but the modular test proves the composed part square-free on
        // its coefficients as written (the compact-exponent rule, docs/resource-policy.md), so the preparation is
        // not reported as limited.
        require(result.statistics.exact_factors == 0, "test did not exceed exact expansion budget");
        require(result.squarefree == CPS_SQUAREFREE_PROVED && !result.progress.preprocessing_limited,
                "the square-free proof on the coefficients as written is lost");
        require(std::all_of(result.roots.begin(), result.roots.end(), [](const auto &r) {
                    return r.certificate.find("exact-power-composition") != std::string::npos;
                }), "coefficient power identity was not used");
        known_cover(wrap(std::move(result)), roots, 53, goal);
        options.bits = 160;
        options.work_units = 1024;
        d::SolveState tighter(certpolysolve::detail::dense(*polynomial.get()), policy);
        result = tighter.advance(&options);
        known_cover(wrap(std::move(result)), roots, 160, goal);
    }
    polynomial.set(4, certpolysolve::ball(certpolysolve::Number("-17@-1000000"), certpolysolve::Number("1@-1000020")));
    require(d::coefficient_power(certpolysolve::detail::dense(*polynomial.get()).coefficients) == 1,
            "early exact reduction admitted coefficient uncertainty");
    polynomial.set(4, certpolysolve::Number("-17@-1000000")).set(0, certpolysolve::Number("0"));
    require(d::coefficient_power(certpolysolve::detail::dense(*polynomial.get()).coefficients) == 1,
            "early reduction admitted the ramified zero branch");
}

static void bounded_rational_temporaries() {
    const std::vector<Complex> roots{{1}, {-1}, {0, 1}, {0, -1},
                                     {2}, {-2}, {0, 2}, {0, -2}};
    certpolysolve::Polynomial polynomial;
    polynomial.set(0, certpolysolve::Number("16@-2048")).set(4, certpolysolve::Number("-17@-2048")).set(8, certpolysolve::Number("1@-2048"));
    require(polynomial.ok(), polynomial.status().c_str());
    auto policy = d::production_policy();
    policy.early_power_composition = true;
    for (uint64_t budget : {uint64_t(10000), uint64_t(64)}) {
        cps_options options;
        options.bits = 80;
        options.max_work_bits = 2048;
        options.exact_budget = budget; // the exact allowance, with expansion allowed
        options.allow_exponent_expansion = true;
        d::SolveState state(certpolysolve::detail::dense(*polynomial.get()), policy);
        auto result = state.advance(&options);
        // At 64 bits the streamed conversion refuses (no exact factor); the modular test still proves the composed
        // part square-free on its coefficients as written, so the preparation is not reported as limited.
        require((result.statistics.exact_factors == 0) == (budget == 64),
                "streamed conversion did not respect its bit budget");
        require(result.squarefree == CPS_SQUAREFREE_PROVED && !result.progress.preprocessing_limited,
                "the square-free proof is lost");
        known_cover(wrap(std::move(result)), roots, 80, CPS_ROOT_ACCURACY);
    }
    // Unrelated rational denominators exercise the LCM, rather than only a
    // shared binary scale. Independent exact roots establish the conversion.
    const std::vector<Complex> rational_roots{{Rational(1) / 3}, {Rational(2) / 5},
        {Rational(4) / 7}, {Rational(8) / 11}, {Rational(16) / 13}};
    auto rational = input(from_roots(rational_roots));
    d::SolveState state(certpolysolve::detail::dense(*rational.get()), policy);
    cps_options options;
    options.bits = 80;
    auto result = state.advance(&options);
    known_cover(wrap(std::move(result)), rational_roots, 80, CPS_ROOT_ACCURACY);
}

static void composition_cancellation_resume() {
    const std::vector<Complex> roots{{1}, {-1}, {0, 1}, {0, -1},
                                     {2}, {-2}, {0, 2}, {0, -2}};
    certpolysolve::Polynomial polynomial;
    polynomial.set(0, certpolysolve::Number("16@-2048")).set(4, certpolysolve::Number("-17@-2048")).set(8, certpolysolve::Number("1@-2048"));
    require(polynomial.ok(), polynomial.status().c_str());
    for (const auto goal : {CPS_ROOT_ACCURACY, CPS_ROOT_APPROXIMATION})
        for (const unsigned checkpoint : {2u, 5u, 7u}) {
            struct Cancellation {
                bool armed = false;
                unsigned calls = 0;
                unsigned checkpoint;
            } cancellation{false, 0, checkpoint};
            auto policy = d::production_policy();
            policy.early_power_composition = true;
            policy.trace_data = &cancellation;
            policy.trace = [](void *data, const d::TraceEvent &event) {
                if (std::string_view(event.stage) == "preprocessing" &&
                    std::string_view(event.phase) == "begin")
                    static_cast<Cancellation *>(data)->armed = true;
            };
            cps_options options;
            options.bits = 80;
            options.max_work_bits = 2048;
            options.exact_budget = 10000; // the exact allowance, with expansion allowed
            options.allow_exponent_expansion = true;
            options.goal = goal;
            options.cancel_data = &cancellation;
            options.cancel = [](void *data) -> int {
                auto &request = *static_cast<Cancellation *>(data);
                return request.armed && ++request.calls >= request.checkpoint;
            };
            d::SolveState state(certpolysolve::detail::dense(*polynomial.get()), policy);
            auto interrupted = state.advance(&options);
            require(interrupted.termination == CPS_CANCELLED && interrupted.degree == 8,
                    "preprocessing cancellation lost the original input degree");
            require(interrupted.statistics.exact_factors == 0,
                    "cancelled preprocessing committed a partial factorization");
            // Without the cancellation, fresh solves at 80 and 160 bits hold the original eight roots (until 18
            // October 2026 the interrupted state was resumed and tightened, which had to apply the power identity
            // exactly once across calls).
            options.cancel = nullptr;
            d::SolveState whole(certpolysolve::detail::dense(*polynomial.get()), policy);
            auto solved = whole.advance(&options);
            known_cover(wrap(std::move(solved)), roots, 80, goal);
            options.bits = 160;
            d::SolveState tighter(certpolysolve::detail::dense(*polynomial.get()), policy);
            auto tightened = tighter.advance(&options);
            known_cover(wrap(std::move(tightened)), roots, 160, goal);
        }
}

int main() {
    certpolysolve::ThreadScope lifetime;
    try {
        require(!d::EnginePolicy{}.prepared_source_refinement, "explicit C2 control changed");
        rigorous_regeneration();
        known_complex_solve();
        sparse_source_counts();
        composition_without_denominator_expansion();
        bounded_rational_temporaries();
        composition_cancellation_resume();
        std::cout << "prepared evaluation sources preserve exact identities and certified roots\n";
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
