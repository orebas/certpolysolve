#include "solve/local_model_solver.hpp"
// Factor preparation, numerical steps and cover assembly for the global driver.
// Query-level control is in solver_advance.cpp; native and program routes have
// their own controllers. compare_numbers and solve are here too.
//
// On coefficient input, advance first proves the degree, enclosing the input at doubling precision until the leading
// coefficient excludes zero, and publishes Cauchy's square (cauchy_cluster). preprocess sets the exact zero roots
// apart (zero_order: the point 0, which aggregate adds to every cover) and prepares the nonzero part: square-free
// by the modular proof (squarefree_modular), else factored exactly over Z (squarefree_factors) or else over Q(i)
// (exact_complex_factors), or kept as one factor of its coefficient balls; a factor g(x^k) with g(0) nonzero is
// solved as g with power k (reduce_power_composition). advance
// steps the factors round-robin, and factor_step dispatches, in order: a power composite whose lifted cover misses the
// goal its reduced cover meets (more separation_bits), a blocked factor, a real-only query (real_only_step:
// RealRefinement), the research build's Portfolio engine (portfolio_subdivision_step), a pending local refinement, a
// factor Subdivision owns (subdivision_step), a structured factor (dovetail_turn: a Subdivision turn when its metered
// share is due, DovetailSchedule; else structured_step: SecularRefinement::step), and else candidate_step: a factor
// the policy sends to that iteration (built on a PolynomialSource), or the classical pass: exact centring
// (center_integer, center_gaussian), isolate_complex, the precision doubled. At the precision cap a complex factor
// passes to Subdivision, the complete method (candidate_handoff, structured_step).
//
// aggregate assembles the cover once every factor is ready: it lifts power covers (lift_power_cover), checks the
// counts against the degree, merges overlapping regions, proves regions real by reflection when the coefficients are
// real (prove_real_by_reflection) and publishes the cover if it meets the goal or prefer_cover ranks it ahead. Where
// finished factors' covers overlap, advance raises the accuracy goal of the wider region's factor: its separation_bits
// become the larger of 16 bits more and the goal's bits plus twice their old value. A full pass of finished or blocked
// factors ends the call with stall_reason (iteration budget, then cells, precision, input); a precision stop on integer
// factors reports the cap Mahler's bound shows sufficient (separation_bound). At a stop short of the goal that is not a
// cancellation, a factor that has had no step yet contributes its stop cover (stop_cover: complete_cover_at_stop on
// its own polynomial in x) and aggregate publishes the factors' covers, as piecewise completes an unfinished factor.
//
// native_advance (solver_native.cpp) runs SecularRefinement on a SecularSource seeded at the poles and falls back to
// the coefficients (secular_polynomial) at its precision stop. ProgramRouteState (solver_program.cpp) proves the
// degree, leading coefficient and a root bound (prepare_program), solves a product of distinct factors by nested
// solves and a power c q^k through q, places seeds by recurrence_seeds or program_seeds, runs the binary64 stage
// and ball packets of Aberth sweeps until the points converge or stall, then SeparatedRootsPass, and
// SecularRefinement on a ProgramSource where that pass fails.
//
// State: SolveState::Impl (the private solve_state_impl.hpp) for one advance (configure throws on a second): input,
// policy, settings, result; factors (FactorState, factor_state.hpp; dovetail_search holds the dovetailed
// Subdivision); degree_proved, preprocessed, terminal, initial_precision, cursor; native_secular;
// program_route owns the separate program controller's state.
// Entry points: SolveState (advance, hint_real_line_seeds), built by session.cpp (the global and program routes,
// piecewise's hand-off to secular), SecularRefinement::group_model (secular_refinement.cpp) and the program route;
// prepared_structure, with no caller in src/; compare_numbers, by validate_query (c_api.cpp) for a region's ends;
// solve, which advances a SessionState once, by cps_solve (c_api.cpp), in its cps_polynomial form only.
// Not here: routing and recipes (session.cpp, auto_selection.cpp), the methods (classical_isolation.cpp,
// secular_refinement.cpp, subdivision.cpp, dovetail.hpp, prepared_piecewise.cpp), centring (factor_state.cpp), covers
// and goals (certification.cpp), exact preparation (exact_preprocessing.cpp, exact_complex.cpp, power_composition.cpp).
// Design notes: docs/mathematics.md (subdivision, Cauchy's bound, export), docs/principled-audit-20260929.md
// (dovetailing), docs/session-state-design-note-20261017.md (one advance), docs/slp-input-spec-20261015.md,
// docs/slp-second-version-note-20261017.md and docs/slp-structural-seeding-note-20261018.md (program input).
#include "certify/power_cover.hpp"
#include "numeric/candidates/candidates.hpp"
#include "methods/classical/classical_isolation.hpp"
#include "numeric/exact/centering.hpp"
#include "certify/counting.hpp"
#include "core/dovetail.hpp"
#include "numeric/exact/exact_complex.hpp"
#include "numeric/exact/exact_preprocessing.hpp"
#include "certify/factor_cover.hpp"
#include "solve/factor_state.hpp"
#include "numeric/exact/modular_squarefree.hpp"
#include "core/cover.hpp"
#include "input/coefficients/polynomial_input.hpp"
#include "core/query.hpp"
#include "core/engine_policy.hpp"
#include "solve/solve_state.hpp"
#include "solve/solve_state_impl.hpp"
#include "core/service_meter.hpp"
#include "methods/real/real_refinement.hpp"
#include "input/secular/secular.hpp"
#include "methods/secular/secular_refinement.hpp"
#include "certify/separation.hpp"
#include "methods/subdivision/subdivision.hpp"
#include "core/telemetry.hpp"
#include "numeric/exact/power_composition.hpp"
#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
#include "portfolio.hpp"
#endif
#include "solve/program_route_state.hpp"
#include "solve/real_line_routing.hpp"
#include "numeric/candidates/real_line_seeds.hpp"
#include "numeric/candidates/study_candidates.hpp"
#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
#include "sturm.hpp"
#endif
#include <algorithm>
#include <chrono>
#include <cmath>
#include <flint/arb_fmpz_poly.h>
#include <limits>
#include <numeric>
#include <optional>
#include <span>
#include <locale>
#include <sstream>

namespace certpolysolve::detail {
namespace {
// The exact preparation's verdict on the input (C14, cps_result_squarefree): the factors are the square-free
// decomposition of the input's nonzero part, each a polynomial in x^power with an exact multiplicity, beside the
// zero roots set apart, counted zero_order times. A multiplicity above one or a zero root of order two or more
// proves a repeated root, as does a factor in x^power with power above one that vanishes at 0 (a real-only query
// keeps x among its factors: its root 0 then has order power); otherwise every root is simple. A ball placeholder
// (no exact preparation) proves nothing beyond the zero roots.
cps_squarefree squarefree_verdict(const std::vector<std::unique_ptr<FactorState>> &factors, uint64_t zero_order) {
    // Exact zero roots set apart: a zero of order 2 or more is a repeated root, whatever the other factors are, and
    // c x is square-free.
    if (zero_order > 1) return CPS_SQUAREFREE_DISPROVED;
    if (factors.empty()) return zero_order ? CPS_SQUAREFREE_PROVED : CPS_SQUAREFREE_UNKNOWN;
    for (const auto &f : factors)
        if (f->kind == FactorState::Ball && !f->squarefree_proved) return CPS_SQUAREFREE_UNKNOWN;
    for (const auto &f : factors) {
        const bool vanishes_at_zero = f->kind == FactorState::Real
                                          ? fmpz_poly_length(f->integer.value()) > 0 && fmpz_is_zero(f->integer.value()->coeffs)
                                          : !f->gaussian.empty() && f->gaussian.front().zero();
        if (f->multiplicity > 1 || (f->power > 1 && vanishes_at_zero)) return CPS_SQUAREFREE_DISPROVED;
    }
    return CPS_SQUAREFREE_PROVED;
}

bool same_number(const Number &a, const Number &b) {
    return (a.zero() && b.zero()) || (a.radix == b.radix && fmpz_equal(a.exponent.value(), b.exponent.value()) &&
                                      fmpq_equal(a.significand.value(), b.significand.value()));
}

// The all-real seed pass at the start of a factor's iteration (cps_options::real_line_seeds): computed once a factor,
// on the polynomial the start seeds from (its centred form when the study centres), its seeds replace the first points
// of the circles' start and the circles' other points stay (design note, "Fewer seeds than roots"). A pass that stops
// gives no seeds and leaves the start as it was.
void seed_from_real_line(RootVector &start, FactorState &f, bool centred, const cps_options &o,
                         const EnginePolicy &policy, CandidateStatistics &statistics) {
    auto &state = f.real_line_seeds;
    const fmpz_poly_struct *polynomial = centred ? f.centered.value() : f.integer.value();
    if (fmpz_poly_degree(polynomial) != start.size()) return;
    if (!state.computed || state.centred != centred) {
        state = {};
        state.computed = true;
        state.centred = centred;
        TraceScope trace(&policy, "real_line_seeds", f.precision, static_cast<uint64_t>(start.size()), &statistics);
        auto pass = real_line_seeds(polynomial, real_line_seed_limits(o));
        state.seeds = std::move(pass.seeds);
        state.statistics = pass.statistics;
    }
    const slong k = std::min<slong>(static_cast<slong>(state.seeds.size()), start.size());
    for (slong i = 0; i < k; ++i) acb_set_d(start[i], state.seeds[static_cast<size_t>(i)]);
}

// The test of a cover that factor_step's steps share (attained): every root meets o's goal at precision.
auto goal_check(const cps_options &o) {
    return [&o](const std::vector<Root> &roots, uint64_t precision) {
        return std::all_of(roots.begin(), roots.end(), [&](const Root &root) {
            return meets_goal(root, o, static_cast<slong>(precision));
        });
    };
}
uint64_t prepared_height(const FactorState &f) {
    if (f.kind == FactorState::Real) {
        const slong bits = fmpz_poly_max_bits(f.integer.value());
        return static_cast<uint64_t>(bits < 0 ? -bits : bits);
    }
    // Clear one common denominator, then remove integer content. This is
    // exact preparation, done once; no floating-point height estimate enters Q.
    Integer denominator;
    fmpz_one(denominator.value());
    for (const auto &c : f.gaussian)
        for (const auto *part : {&c.real, &c.imag})
            fmpz_lcm(denominator.value(), denominator.value(), fmpq_denref(part->value()));
    std::vector<Integer> coefficients;
    Integer content;
    for (const auto &c : f.gaussian)
        for (const auto *part : {&c.real, &c.imag}) {
            auto &a = coefficients.emplace_back();
            fmpz_divexact(a.value(), denominator.value(), fmpq_denref(part->value()));
            fmpz_mul(a.value(), a.value(), fmpq_numref(part->value()));
            fmpz_gcd(content.value(), content.value(), a.value());
        }
    uint64_t height = 0;
    for (auto &a : coefficients) {
        if (!fmpz_is_zero(content.value())) fmpz_divexact(a.value(), a.value(), content.value());
        height = std::max<uint64_t>(height, fmpz_bits(a.value()));
    }
    return height;
}
} // namespace
std::optional<int> compare_numbers(const Number &a, const Number &b, const cps_options &options) {
    if (same_number(a, b)) return 0;
    Rational x, y;
    uint64_t budget = options.exact_storage_limit();
    if (a.rational_with_budget(x.value(), budget, options.allow_exponent_expansion) &&
        b.rational_with_budget(y.value(), budget, options.allow_exponent_expansion)) {
        const int order = fmpq_cmp(x.value(), y.value());
        return (order > 0) - (order < 0);
    }
    const auto ceiling = static_cast<slong>(options.max_work_bits);
    for (slong precision = std::min<slong>(128, ceiling); !options.cancelled();) {
        RealBall xball, yball;
        a.enclose(xball.value(), precision);
        b.enclose(yball.value(), precision);
        if (arb_lt(xball.value(), yball.value())) return -1;
        if (arb_gt(xball.value(), yball.value())) return 1;
        if (arb_is_exact(xball.value()) && arb_is_exact(yball.value()) && arb_equal(xball.value(), yball.value()))
            return 0;
        if (precision == ceiling) break;
        precision = std::min(ceiling, resources::precision_product(precision, 2));
    }
    return std::nullopt;
}
SolveState::Impl::Impl(DensePolynomial p, EnginePolicy e) : input(std::move(p)), policy(e) {
    source_structure = input_structure(input.coefficients);
    input.coefficients.resize(source_structure->length);
    result.degree = static_cast<int64_t>(input.coefficients.size()) - 1;
    result.progress.input_uncertainty = !source_structure->exact;
}
SolveState::Impl::Impl(const cps_secular &p, EnginePolicy e) : secular_input(std::make_unique<cps_secular>(p)), policy(e) {
    result.degree = static_cast<int64_t>(p.terms.size());
}
SolveState::Impl::Impl(const cps_program &p, EnginePolicy e)
    : program_route(std::make_unique<ProgramRouteState>(p)), policy(e) {
    result.degree = program_route->degree_bound();
}
SolveState::Impl::~Impl() = default;

void SolveState::Impl::configure(const cps_options *proposed) {
    if (configured)
        throw std::logic_error("a solve state is advanced once");
    settings = proposed ? *proposed : cps_options{};
    settings = with_resolved_goal(settings, program_route != nullptr);
    settings.cancel = nullptr;
    settings.cancel_data = nullptr;
    initial_precision = std::min(settings.max_work_bits, std::max<uint64_t>(64, settings.bits + 32));
    if (policy.prepared_subdivision)
        initial_precision = std::min<uint64_t>(settings.max_work_bits, 64);
    configured = true;
    result.goal = settings.goal;
}
cps_options SolveState::Impl::factor_options(const FactorState &f, const cps_options &base) const {
    auto options = separated_options(base, f.separation_bits);
    // Isolation has no width goal. Once two factors' regions meet, adding
    // bits alone would therefore do no work: refine the affected factor to
    // an internal accuracy goal while keeping the exported isolation goal.
    if ((policy.prepared_subdivision || f.complete_fallback) && f.separation_bits && options.goal == CPS_ROOT_ISOLATION)
        options.goal = CPS_ROOT_ACCURACY;
    if ((policy.prepared_subdivision || f.complete_fallback) && !f.subdivision_isolated && !f.separation_bits) {
        options.goal = CPS_ROOT_ISOLATION;
        options.bits = 0;
        options.absolute_tolerance = Number();
    }
    // Every recomputation, including power-lift/cross-factor tightening,
    // retains the same shared allowance. Complete service restores caller caps.
    options.heuristic_budget = f.complete_fallback ? nullptr : settings.heuristic_budget;
    if (options.heuristic_budget) {
        options.max_work_bits = std::min(options.max_work_bits, options.heuristic_budget->ceiling());
        options.exact_budget = std::min(options.exact_budget, options.heuristic_budget->exact_limit());
        options.exact_operations = std::min(options.exact_operations, options.heuristic_budget->exact_limit());
    }
    return options;
}
bool SolveState::Impl::factor_goal(const FactorState &f) const {
    if (f.power > 1) {
        if (f.lifted_cover.empty() || !disjoint(f.lifted_cover) ||
            !goal(f.lifted_cover, std::max(f.precision, result.progress.max_outer_bits)))
            return false;
        // The lifted cover answers the query. After the overlap branch has raised the factor's separation bits,
        // its reduced cover must also meet that tighter goal: otherwise the factor would count as finished at
        // once and never be stepped, and the branch would raise the goal again at no work up to its
        // representable limit, which would then stop the solve "at the precision cap". The regression test:
        // (z^2 - 2)(z - c)^2 with c within 2^-130 of sqrt 2 (tests/test_power_overlap.cpp). History:
        // docs/code-history-20261021.md, solver.cpp 3.
        // factor_step's own tightening of a power composite (the inverse map's magnification) keeps the query's
        // goal as the test: it raises the bits only once the reduced cover meets the factor's goal already.
        if (!f.overlap_tightened)
            return true;
        const auto o = factor_options(f, settings);
        return std::all_of(f.cover.begin(), f.cover.end(), [&](const Root &r) {
            return meets_goal(r, o, static_cast<slong>(std::max(f.precision, result.progress.max_outer_bits)));
        });
    }
    if (!f.separation_bits)
        return goal(f.cover, std::max(f.precision, result.progress.max_outer_bits));
    const auto o = factor_options(f, settings);
    return std::all_of(f.cover.begin(), f.cover.end(), [&](const Root &r) {
        return meets_goal(r, o,
                          static_cast<slong>(std::max(f.precision, result.progress.max_outer_bits)));
    });
}
void SolveState::Impl::enclose_input(acb_poly_t p, slong precision, bool observe) {
    StageTimer timer(result.statistics.input_enclosure_ns);
    ComplexBall coefficient;
    for (size_t i = 0; i < input.coefficients.size(); ++i) {
        input.coefficients[i].enclose(coefficient.value(), precision);
        acb_poly_set_coeff_acb(p, static_cast<slong>(i), coefficient.value());
    }
    if (observe)
        observe_precision(precision);
}
void SolveState::Impl::enclose_factor(FactorState &f, acb_poly_t p, slong precision, bool observe) {
    if (f.kind == FactorState::Ball) {
        enclose_input(p, precision, observe);
        return;
    }
    StageTimer timer(result.statistics.input_enclosure_ns);
    if (f.kind == FactorState::Real)
        acb_poly_set_fmpz_poly(p, f.integer.value(), precision);
    else
        enclose_gaussian(p, f.gaussian, precision);
    if (observe)
        observe_precision(precision);
}
void SolveState::Impl::preprocess(const cps_options &o) {
    StageTimer timer(result.statistics.preprocessing_ns);
    // The exact zero roots apart (S8's B1), as prepare_squarefree_input sets them apart for piecewise: the v lowest
    // coefficients that are exact zeros (InputStructure::start; an enclosure that merely contains 0 is not one) leave
    // the point 0, counted v times, beside the factors in every assembly (aggregate), and all that follows prepares
    // the nonzero part q = p / x^v. Not for a real-only query, whose engine counts x among its factors.
    const size_t zeros = source_structure && !o.real_only ? source_structure->start : 0;
    const std::span<const Coefficient> nonzero = std::span<const Coefficient>(input.coefficients).subspan(zeros);
    // A coefficient-index identity needs no denominator expansion. Preserve
    // it even when exact squarefree preprocessing exceeds its own budget.
    // Commit the smaller owned input only after preprocessing succeeds, so
    // a cancelled or refused preprocessing leaves the input as it was.
    std::optional<ExactPowerComposition> composed;
    if (policy.early_power_composition && policy.backend && !o.real_only &&
        (!source_structure || source_structure->nonzero_part_power() > 1))
        composed = prepare_power_composition(nonzero);
    const ulong source_power = composed ? composed->power : 1;
    const std::span<const Coefficient> preparation =
        composed ? std::span<const Coefficient>(composed->polynomial.coefficients) : nonzero;
    IntegerPolynomial exact;
    PreparationDiagnostic diagnostic;
    const auto import_started = std::chrono::steady_clock::now();
    const bool imported = streamed_real_polynomial(exact, preparation, o, &diagnostic);
    // The modular square-free proof first (S8's B2), the proof piecewise and the selector's admission run: when it
    // proves q square-free, q's one factor is the import's primitive part with a positive leading coefficient, the
    // factor FLINT's square-free factorisation returns for a square-free input, and FLINT is not called. Any other
    // outcome says nothing about q, and FLINT decides as before.
    bool proved_squarefree = false;
    if (imported && preparation.size() >= 2) {
        const auto proof = squarefree_modular(preparation, o);
        if (proof == ModularSquarefreeResult::Cancelled) {
            cancel_observed = true;
            return;
        }
        proved_squarefree = proof == ModularSquarefreeResult::Proved;
    }
    // A refused import leaves the coefficients as written. The modular test reduces them directly (each m@e number
    // modulo the prime, its exponent modulo p - 1), so on real input it still runs, and a pass is kept on the ball
    // factor below. Complex input takes the same test inside exact_complex_factors. A zero exact allowance disables
    // exact preparation, this test with it.
    bool proved_as_written = false;
    if (!imported && preparation.size() >= 2 && o.exact_storage_limit() &&
        std::all_of(preparation.begin(), preparation.end(),
                    [](const Coefficient &c) { return c.exact() && c.is_real(); })) {
        const auto proof = squarefree_modular(preparation, o);
        if (proof == ModularSquarefreeResult::Cancelled) {
            cancel_observed = true;
            return;
        }
        proved_as_written = proof == ModularSquarefreeResult::Proved;
    }
    // The square-free factorisation is charged too: it cannot be
    // cancelled once started. A refusal takes the refused-import route.
    const bool have_real = imported && (proved_squarefree || o.exact_effort_fits(squarefree_effort(exact)));
    if (imported && !have_real) {
        diagnostic.reason = "squarefree-effort";
        diagnostic.accepted = false;
    }
    observe_preparation(policy, "exact-real-import", diagnostic, o.exact_budget, import_started, preparation.size());
    std::vector<std::unique_ptr<FactorState>> prepared;
    if (have_real) {
        const auto factor_started = std::chrono::steady_clock::now();
        std::vector<IntegerFactor> found;
        if (proved_squarefree) {
            found.emplace_back();
            fmpz_poly_primitive_part(found.back().polynomial.value(), exact.value());
        } else {
            found = squarefree_factors(exact);
        }
        for (auto &factor : found) {
            auto f = std::make_unique<FactorState>();
            f->kind = FactorState::Real;
            f->integer = std::move(factor.polynomial);
            f->multiplicity = factor.multiplicity;
            // The pass already ran on this polynomial before routing: its outcome is the factor's, uncentred
            // (seed_from_real_line runs it again only on a centred start).
            if (real_line_hint && fmpz_poly_equal(f->integer.value(), real_line_hint->polynomial.value())) {
                f->real_line_seeds.computed = true;
                f->real_line_seeds.centred = false;
                f->real_line_seeds.seeds = real_line_hint->pass.seeds;
                f->real_line_seeds.statistics = real_line_hint->pass.statistics;
            }
            prepared.push_back(std::move(f));
        }
        diagnostic.reason = proved_squarefree ? "modular-squarefree-proof" : "exact-squarefree-factors";
        diagnostic.accepted = true;
        observe_preparation(policy, "exact-squarefree", diagnostic, o.exact_budget,
                            factor_started, prepared.size());
        preprocessing_status = "budgeted exact real squarefree preprocessing";
    } else if (!o.real_only) {
        auto decomposition = exact_complex_factors(preparation, o, &policy);
        if (decomposition.cancelled) {
            cancel_observed = true;
            return;
        }
        preprocessing_status = decomposition.status == "Q(i) preprocessing not applicable"
            ? diagnostic.reason : decomposition.status;
        if (proved_as_written)
            preprocessing_status += "; squarefreeness proved modulo 998244353 on the coefficients as written";
        proved_as_written = proved_as_written || decomposition.squarefree_proved;
        for (auto &factor : decomposition.factors) {
            auto f = std::make_unique<FactorState>();
            f->kind = FactorState::Gaussian;
            f->gaussian = std::move(factor.polynomial);
            f->multiplicity = factor.multiplicity;
            prepared.push_back(std::move(f));
        }
    }
    if (o.cancelled())
        return;
    result.statistics.exact_factors = prepared.size();
    // Limited: a refused exact path left the multiplicities unknown. A square-free proof on the coefficients as
    // written makes them known (every root simple), as piecewise's preparation counts a proof without factors
    // (SquarefreePreparation::available).
    result.progress.preprocessing_limited =
        !have_real && prepared.empty() && !result.progress.input_uncertainty && !proved_as_written;
    if (o.real_only && !have_real) {
        result.termination = CPS_UNSUPPORTED_QUERY;
        result.progress.stop_reason = CPS_STOP_INPUT;
        terminal = true;
        preprocessing_status =
            "real_only requires exact real input within the exact preprocessing budget";
        return;
    }
    if (policy.prepared_subdivision && prepared.empty() && preparation.size() > 1) {
        // A modular proof on compressed coefficients alone is not the
        // prepared-factor service. Keep the initial counted cover, but do
        // not silently substitute the uncertain/legacy factor route.
        result.termination = CPS_UNSUPPORTED_QUERY;
        result.progress.stop_reason = CPS_STOP_INPUT;
        result.progress.preprocessing_limited = true;
        // Public input refusals are not complete answers (I4). Previously
        // certified regions may remain visible, including the initial bound.
        result.complete = result.goal_met = false;
        terminal = true;
        preprocessing_status += "; standalone subdivision requires exact squarefree factors within the preprocessing budget";
        return;
    }
    // A nonzero part of degree 0 (p = c x^v) leaves no factor: the point 0 is the whole answer (aggregate, below).
    // Otherwise a preparation without factors keeps one factor of the coefficient balls.
    if (prepared.empty() && preparation.size() > 1) {
        prepared.push_back(std::make_unique<FactorState>());
        prepared.back()->squarefree_proved = proved_as_written;
    }
    uint64_t factor_id = 0;
    for (auto &f : prepared) {
        f->id = factor_id++;
        f->precision = initial_precision;
        if (policy.backend && !o.real_only) {
            if (f->kind == FactorState::Real)
                f->power = reduce_power_composition(f->integer);
            else if (f->kind == FactorState::Gaussian)
                f->power = reduce_power_composition(f->gaussian);
            f->power *= source_power;
        }
    }
    // The smaller owned input, which a factor of the coefficient balls encloses: the composed nonzero part, or the
    // nonzero part.
    if (composed)
        input = std::move(composed->polynomial);
    else if (zeros)
        input.coefficients.erase(input.coefficients.begin(), input.coefficients.begin() + static_cast<ptrdiff_t>(zeros));
    zero_order = zeros;
    factors = std::move(prepared);
    if (policy.complete_fallback &&
        std::any_of(factors.begin(), factors.end(), [](const auto &f) { return f->kind == FactorState::Ball; })) {
        // Exact input can still exceed the caller's preparation allowance.
        // The delayed warm schedule requires exact prepared factors; applying
        // it to the ball placeholder would suppress the maintained early
        // interleaver without installing complete service. Decide eligibility
        // here, before any factor work, and keep that existing counted path.
        policy.complete_fallback = false;
        policy.research_integrated_dovetail = false;
        policy.research_service_meter = false;
        policy.research_identity_shift_meter = false;
        policy.research_local_packet_meter = false;
    }
    if (policy.complete_fallback && !settings.real_only && source_structure && source_structure->exact &&
        !settings.heuristic_budget && !factors.empty() &&
        std::all_of(factors.begin(), factors.end(), [](const auto &f) { return f->kind != FactorState::Ball; })) {
        uint64_t height = 0;
        for (auto &f : factors) {
            f->exact_height = prepared_height(*f);
            height = std::max(height, f->exact_height);
        }
        heuristic_budget.emplace(static_cast<uint64_t>(result.degree), height,
            settings.goal == CPS_ROOT_ISOLATION ? 0 : settings.bits,
            policy.heuristic_packets, policy.heuristic_precision_multiplier);
        settings.heuristic_budget = &*heuristic_budget;
    }
    result.squarefree = squarefree_verdict(factors, zero_order);
    preprocessed = true;
    if (factors.empty())
        aggregate();
}
uint64_t SolveState::Impl::subdivision_start(const cps_options &o) const {
    return std::min<uint64_t>(
        o.max_work_bits,
        std::max<uint64_t>({initial_precision, 64, static_cast<uint64_t>(resources::precision_sum(o.bits, 32))}));
}
void SolveState::Impl::candidate_handoff(FactorState &f, const cps_options &o) {
    f.search = std::make_unique<Subdivision>(f.cover, f.multiplicity, subdivision_start(o), result.statistics_v2);
    if (policy.trace) {
        TraceEvent event{"candidate_subdivision_handoff", "transition", f.precision,
                         static_cast<uint64_t>(f.cover.size()), &result.statistics_v3, &result.statistics_v2};
        event.engine = "subdivision";
        event.reason = "precision-budget";
        policy.trace(policy.trace_data, event);
    }
}
uint64_t SolveState::Impl::live_cells() const {
    uint64_t n = 0;
    for (const auto &f : factors) {
        if (f->search)
            n += f->search->cells();
        else if (f->real_search)
            n += f->real_search->cells();
        if (f->dovetail_search) n+=f->dovetail_search->cells();
#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
        if (f->sturm) n+=f->sturm->cells();
#endif
    }
    return n;
}
bool SolveState::Impl::stop_cover(FactorState &f, const cps_options &o, std::vector<Root> &cover) {
    const auto precision = static_cast<slong>(std::max<uint64_t>(64, f.precision));
    ComplexPolynomial stored, own;
    enclose_factor(f, stored.value(), precision, false);
    if (f.power > 1) {
        const auto power = static_cast<slong>(f.power);
        for (slong i = 0; i <= acb_poly_degree(stored.value()); ++i)
            acb_poly_set_coeff_acb(own.value(), i * power, stored.value()->coeffs + i);
    } else {
        acb_poly_swap(own.value(), stored.value());
    }
    cover.clear();
    if (!complete_cover_at_stop(cover, own.value(), precision, &o) || o.cancelled())
        return false;
    for (auto &root : cover) {
        root.count *= f.multiplicity;
        root.multiplicity *= f.multiplicity;
    }
    return true;
}
void SolveState::Impl::aggregate(const cps_options *stop) {
    if (!preprocessed || (factors.empty() && !zero_order))
        return;
    std::vector<Root> cover;
    uint64_t expected = 0;
    if (zero_order)
        cover.push_back(exact_zero_point(zero_order)); // the exact zero roots set apart (preprocess)
    for (const auto &f : factors) {
        if (!f->ready) {
            // A factor that has had no step yet. At a stop short of the goal it contributes its stop cover, so the
            // answer holds the other factors' regions and a counted square for its roots, as piecewise's factor loop
            // completes an unfinished factor; otherwise nothing is published until every factor has a cover (the
            // degree proof's Cauchy square stays). Real-only queries publish real roots only and take no stop cover.
            std::vector<Root> own;
            if (!stop || settings.real_only || !stop_cover(*f, *stop, own))
                return;
            cover.insert(cover.end(), own.begin(), own.end());
            continue;
        }
        if (f->power > 1) {
            f->precision =
                std::min(settings.max_work_bits, std::max({f->precision, result.progress.max_outer_bits,
                                                           static_cast<uint64_t>(resources::precision_sum(resources::precision_sum(settings.bits, f->separation_bits), 32))}));
            result.progress.max_outer_bits = std::max(result.progress.max_outer_bits, f->precision);
            result.work_bits = f->precision;
            TraceScope trace(&policy, "power_lifting", f->precision, f->cover.size(),
                             &result.statistics_v3);
            if (!lift_power_cover(f->lifted_cover, f->cover, f->power,
                                  static_cast<slong>(f->precision))) {
                f->lifted_cover.clear();
                return;
            }
        }
        const auto &output = f->output_cover();
        cover.insert(cover.end(), output.begin(), output.end());
        if (settings.real_only)
            expected += f->real_count;
    }
    if (!settings.real_only)
        expected = static_cast<uint64_t>(result.degree);
    assemble_factor_covers(cover, expected, static_cast<slong>(result.progress.max_outer_bits));
    if (std::all_of(input.coefficients.begin(), input.coefficients.end(),
                    [](const Coefficient &c) { return c.is_real(); })) {
        prove_real_by_reflection(cover, static_cast<slong>(result.progress.max_outer_bits));
    }
    const bool attained = goal(cover, result.progress.max_outer_bits);
    const bool published =
        publish_assembly(result, cover, attained, static_cast<slong>(result.progress.max_outer_bits));
    // A published assembly is result.roots now, the same roots at the same precision: its test is `attained`.
    result.goal_met = result.complete && (published ? attained : goal(result.roots, result.progress.max_outer_bits));
}
cps_stop_reason SolveState::Impl::factor_step(FactorState &f, const cps_options &base) {
    TraceScope trace(&policy, f.search ? "subdivision" : "factor", f.precision,
                     static_cast<uint64_t>(std::max<slong>(0, f.candidates.unresolved_count)),
                     &result.statistics_v3, &result.statistics_v2);
    auto o = factor_options(f, base);
    if (heuristic_budget && heuristic_budget->exhausted() && !f.complete_fallback)
        return start_complete_fallback(f, base);
    if (f.power > 1 && f.ready && !factor_goal(f) &&
        std::all_of(f.cover.begin(), f.cover.end(),
                    [&](const Root &r) { return meets_goal(r, o, static_cast<slong>(f.precision)); })) {
        // The inverse map may magnify an absolute error near zero, or two
        // lifted rectangles may overlap. Tighten only this reduced factor;
        // do not confuse its attained goal with the original query's goal.
        if (f.precision >= base.max_work_bits ||
            f.separation_bits >= resources::precision_max - base.bits)
            return CPS_STOP_PRECISION;
        f.separation_bits = std::min<uint64_t>(resources::precision_max - base.bits, f.separation_bits + 32);
        o = factor_options(f, base);
        if (f.dovetail_search) f.dovetail_search->reopen(o);
        if (f.search)
            f.search->reopen(o);
        else if (!f.structured)
            f.local_refinement = true;
    }
    if (f.blocked)
        return f.input_limited ? CPS_STOP_INPUT : CPS_STOP_PRECISION;
    if (o.real_only) return real_only_step(f, o);
#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
    if (policy.candidates==CandidateEngine::Portfolio && f.ready && !f.search &&
        ++f.portfolio_steps%controller_config::subdivision_period==0) return portfolio_subdivision_step(f, o);
#endif
    if (f.local_refinement && f.dovetail_search) {
        f.search=std::move(f.dovetail_search);f.search->reopen(o);f.local_refinement=false;
        return CPS_STOP_NONE;
    }
    if (f.local_refinement) {
#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
        if (policy.research_no_subdivision)
            return research_handoff(CPS_STOP_PRECISION);
#endif
        if (live_cells() + f.cover.size() > o.live_cells)
            return CPS_STOP_CELL_BUDGET;
        f.search = std::make_unique<Subdivision>(f.cover, f.multiplicity, subdivision_start(o),
                                                 result.statistics_v2,
                                                 policy.prepared_subdivision ? SubdivisionPrecision::TrialLocal
                                                                             : SubdivisionPrecision::Retained);
        f.local_refinement = false;
        return CPS_STOP_NONE;
    }
    if (f.search) {
        const auto reason = subdivision_step(f, o);
        if (reason == CPS_STOP_PRECISION && policy.research_integrated_dovetail &&
            !policy.prepared_subdivision && !f.complete_fallback &&
            f.search->precision_policy() == SubdivisionPrecision::TrialLocal &&
            !factor_goal(f)) {
            // The warm worker can hit its trial guard after EA has handed
            // over. Try the maintained finite-cap path once, from the best
            // certificate. Its Retained policy prevents another such restart.
            // subdivision_step has already settled the old worker's meters.
            if (o.heuristic_budget) {
                // Record the retiring worker's peak before replacement can
                // hide it from advance_factor's ordinary precision check.
                o.heuristic_budget->observe_precision(f.search->max_bits());
                if (o.heuristic_budget->exhausted())
                    return heuristic_budget ? CPS_STOP_NONE : CPS_STOP_WORK_BUDGET;
            }
            return start_retained_subdivision(f, o);
        }
        return reason;
    }
    if (policy.prepared_subdivision) return start_prepared_subdivision(f, o);
    if (const auto stop = dovetail_turn(f, o)) return *stop;
    if (f.structured)
        return structured_step(f, o);
    return candidate_step(f, o);
}

cps_stop_reason SolveState::Impl::real_only_step(FactorState &f, const cps_options &o) {
#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
    const auto attained = goal_check(o);
#endif
    StageTimer timer(result.statistics.isolation_ns);
#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
    // The research build's Portfolio engine: its Sturm step and its share of the cells.
    const bool portfolio=policy.candidates==CandidateEngine::Portfolio;
    if (portfolio && policy.real_method!=RealMethod::Descartes && !f.sturm_unavailable &&
        (policy.real_method==RealMethod::Sturm ||
         f.real_stagnant>=controller_config::stagnant_packets*controller_config::subdivision_period) &&
        (policy.real_method==RealMethod::Sturm ||
         ++f.portfolio_steps%controller_config::subdivision_period==0)) {
        TraceScope trace(&policy,"sturm",f.precision,0,&result.statistics_v3,&result.statistics_v2);
        if (!f.sturm) f.sturm=std::make_unique<SturmRefinement>(f.integer.value(),f.multiplicity,
            static_cast<slong>(f.precision),result.statistics_v2);
        auto sturm_options=o;
        const auto other_cells=live_cells()-f.sturm->cells();
        sturm_options.live_cells=o.live_cells>other_cells?o.live_cells-other_cells:0;
        const auto reason=f.sturm->step(sturm_options);
        if (reason==CPS_STOP_INPUT) {
            // An optional exact sequence exceeding its growth budget
            // is an engine rejection, not an unsupported polynomial.
            f.sturm_unavailable=true;f.sturm.reset();return CPS_STOP_NONE;
        }
        f.precision=std::max(f.precision,f.sturm->precision());
        result.work_bits=f.precision;
        result.progress.max_outer_bits=std::max(result.progress.max_outer_bits,f.precision);
        std::vector<Root> cover;
        if (f.sturm->cover(cover) && (!f.ready || attained(cover,f.precision) ||
            prefer_cover(cover,f.cover,static_cast<slong>(f.precision)))) {
            f.cover=std::move(cover);f.real_count=f.sturm->count();f.ready=true;
        }
        return reason;
    }
#endif
    // Each real factor starts with one pending interval. The allowance belongs
    // to the whole solve, including every other factor and research engine.
    if (!f.real_search && live_cells() >= o.live_cells)
        return CPS_STOP_CELL_BUDGET;
    if (!f.real_search)
        f.real_search = std::make_unique<RealRefinement>(f.integer.value(), f.multiplicity,
                                                         static_cast<slong>(f.precision),
                                                         result.statistics_v2, result.statistics_v3);
#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
    const auto previous_count=f.real_search->count();
#endif
    auto real_options=o;
    const auto other_cells = live_cells() - f.real_search->cells();
    real_options.live_cells = o.live_cells > other_cells ? o.live_cells - other_cells : 0;
    const auto reason = f.real_search->step(real_options);
#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
    f.real_stagnant=f.real_search->count()>previous_count?0:f.real_stagnant+1;
#endif
    f.precision = std::max(f.precision, f.real_search->precision());
    result.work_bits = f.precision;
    result.progress.max_outer_bits = std::max(result.progress.max_outer_bits, f.precision);
    std::vector<Root> cover;
    if (f.real_search->cover(cover)) {
        f.cover = std::move(cover);
        f.real_count = f.real_search->count();
        f.ready = true;
    }
    return reason;
}

#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
cps_stop_reason SolveState::Impl::portfolio_subdivision_step(FactorState &f, const cps_options &o) {
    const auto attained = goal_check(o);
    if (!f.dovetail_search) {
        if (live_cells()+f.cover.size()>o.live_cells) return CPS_STOP_CELL_BUDGET;
        f.dovetail_search=std::make_unique<Subdivision>(f.cover,f.multiplicity,f.precision,
            result.statistics_v2);
    }
    TraceScope trace(&policy,"portfolio_subdivision",f.precision,f.dovetail_search->pending(),
                      &result.statistics_v3,&result.statistics_v2);
    const auto used=live_cells();
    const auto reason=f.dovetail_search->step(o,[&](acb_poly_t p,slong bits){enclose_factor(f,p,bits);},
                                               o.live_cells>used?o.live_cells-used:0);
    std::vector<Root> cover;
    if (f.dovetail_search->complete_cover(cover) &&
        (attained(cover,f.precision) || prefer_cover(cover,f.cover,static_cast<slong>(f.precision))))
        f.cover=std::move(cover);
    // A blocked subdivision packet does not block eligible candidates.
    return reason==CPS_STOP_CANCELLED?reason:CPS_STOP_NONE;
}
#endif

namespace {
template<class Function> struct ResearchOnReturn {
    Function function;
    ~ResearchOnReturn() { function(); }
};
uint64_t research_elapsed(std::chrono::steady_clock::time_point started) {
    return static_cast<uint64_t>(std::chrono::duration_cast<std::chrono::nanoseconds>(
        std::chrono::steady_clock::now() - started).count());
}
}
std::unique_ptr<Subdivision> SolveState::Impl::make_research_warm_search(const FactorState &factor,
                                                                      const cps_options &options) {
    if (!factor.ready || factor.kind == FactorState::Ball || factor.power != 1)
        throw std::logic_error("warm subdivision requires an exact unmapped certified factor");
    return std::make_unique<Subdivision>(factor.cover, factor.multiplicity, subdivision_start(options),
                                         result.statistics_v2, SubdivisionPrecision::TrialLocal);
}

cps_stop_reason SolveState::Impl::subdivision_step(FactorState &f, const cps_options &o) {
    if (o.heuristic_budget) o.heuristic_budget->linear();
    auto *metrics = policy.research_scheduling_metrics;
    const double work_before = metrics ? f.search->work() : 0;
    const auto began = metrics ? std::chrono::steady_clock::now() : std::chrono::steady_clock::time_point{};
    ResearchOnReturn account{[&] {
        if (!metrics) return;
        metrics->subdivision_work += std::max(0., f.search->work() - work_before);
        metrics->subdivision_ns += research_elapsed(began);
        ++metrics->subdivision_turns;
    }};
    const auto attained = goal_check(o);
    ExclusiveTimer timer(result.statistics.isolation_ns, &result.statistics.input_enclosure_ns);
    const auto used = live_cells();
    double allowance = std::numeric_limits<double>::infinity();
#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
    if (research) allowance = research->subdivision_allowance;
#endif
    f.search->research_identity_shift_meter(policy.research_identity_shift_meter);
    f.search->research_service_meter(&policy);
    const auto reason = f.search->step(
        o, [&](acb_poly_t p, slong precision) { enclose_factor(f, p, precision); },
        o.live_cells > used ? o.live_cells - used : 0, allowance,
        // TrialLocal workers select PreparedExact internally. A retained
        // terminal handoff keeps its established finite-cap Newton behavior;
        // the integrated scheduler alone does not change that contract.
        (policy.prepared_subdivision || f.complete_fallback || policy.research_warm_subdivision)
            ? NewtonAccuracy::PreparedExact : NewtonAccuracy::Legacy);
#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
    if (research) research->subdivision_turn = {
        f.search->suspended(), f.search->idle(), f.search->pending_work()};
#endif
    std::vector<Root> cover;
    if (f.search->complete_cover(cover) &&
        (attained(cover, result.progress.max_outer_bits) ||
         prefer_cover(cover, f.cover, static_cast<slong>(result.progress.max_outer_bits))))
        f.cover = std::move(cover);
    if ((policy.prepared_subdivision || f.complete_fallback) && !f.subdivision_isolated && o.goal == CPS_ROOT_ISOLATION &&
        attained(f.cover, result.progress.max_outer_bits)) {
        f.subdivision_isolated = true;
        // Reopen the owned isolating geometry, keeping certificates and
        // component precision. No new Cauchy start or exact preparation.
        f.search->reopen(factor_options(f, settings));
    }
    return reason;
}

cps_stop_reason SolveState::Impl::start_retained_subdivision(FactorState &f, const cps_options &o) {
    const auto other_cells = live_cells() - (f.search ? f.search->cells() : 0);
    if (other_cells + f.cover.size() > o.live_cells) return CPS_STOP_CELL_BUDGET;
    if (o.heuristic_budget) o.heuristic_budget->quadratic();
    f.search = std::make_unique<Subdivision>(f.cover, f.multiplicity, subdivision_start(o), result.statistics_v2);
    return CPS_STOP_NONE;
}

cps_stop_reason SolveState::Impl::start_complete_fallback(FactorState &f, const cps_options &o) {
    if (f.kind == FactorState::Ball)
        throw std::logic_error("complete service requires an exact prepared factor");
    if (o.max_work_bits < 64) return CPS_STOP_PRECISION;
    // Keep an already admissible owned worker. Earlier heuristic warm forests
    // have no complete-entry provenance: coalesce their localization once,
    // retaining the fine certified cover for publication throughout the move.
    if (f.dovetail_search && f.dovetail_search->complete_entry())
        f.search = std::move(f.dovetail_search);
    if (!f.search || !f.search->complete_entry()) {
        const auto old_cells = (f.search ? f.search->cells() : 0) +
                               (f.dovetail_search ? f.dovetail_search->cells() : 0);
        if (live_cells() - old_cells >= o.live_cells) return CPS_STOP_CELL_BUDGET;
        ComplexPolynomial polynomial;
        enclose_factor(f, polynomial.value(), 64);
        Root global;
        if (!cauchy_cluster(global, polynomial.value(), 64)) return CPS_STOP_PRECISION;
        global.count *= f.multiplicity;
        global.multiplicity *= f.multiplicity;
        // Cauchy's radius is at most 2^(height+1). Degree covers factor
        // height growth; four additional bits cover rounding and padding.
        const uint64_t scale = resources::add(resources::add(std::max<uint64_t>(1, f.exact_height),
            static_cast<uint64_t>(acb_poly_degree(polynomial.value()))), 4);
        auto next = Subdivision::complete_start(global, f.cover, f.multiplicity, scale,
                                                result.statistics_v2);
        if (!next) return CPS_STOP_PRECISION;
        f.search = std::move(next);
        f.dovetail_search.reset();
        if (!f.ready) f.cover = {std::move(global)};
    }
    f.ready = true;
    f.complete_fallback = true;
    f.blocked = f.local_refinement = false;
    f.subdivision_isolated = false;
    f.search->reopen(factor_options(f, o));
    if (policy.trace) {
        TraceEvent event{"complete_fallback", "transition", 64,
                         static_cast<uint64_t>(f.cover.size()), &result.statistics_v3,
                         &result.statistics_v2};
        event.reason = heuristic_budget->stop() == HeuristicBudget::Stop::Precision ? "heuristic-precision" : "heuristic-packets";
        event.algorithm_cutoff = heuristic_budget->ticks();
        event.exact_budget = heuristic_budget->ceiling();
        policy.trace(policy.trace_data, event);
    }
    return CPS_STOP_NONE;
}

cps_stop_reason SolveState::Impl::start_prepared_subdivision(FactorState &f, const cps_options &o) {
    if (f.kind == FactorState::Ball || f.power != 1)
        throw std::logic_error("standalone subdivision requires an exact prepared factor without a power map");
    if (!f.ready) {
        ComplexPolynomial polynomial;
        enclose_factor(f, polynomial.value(), static_cast<slong>(f.precision));
        Root global;
        if (!cauchy_cluster(global, polynomial.value(), static_cast<slong>(f.precision))) {
            if (f.precision == o.max_work_bits) {
                f.blocked = true;
                return CPS_STOP_PRECISION;
            }
            f.precision = std::min(o.max_work_bits, f.precision * 2);
            ++result.statistics_v2.precision_increases;
            return CPS_STOP_NONE;
        }
        global.count *= f.multiplicity;
        global.multiplicity *= f.multiplicity;
        f.cover = {std::move(global)};
        f.ready = true;
    }
    if (live_cells() + f.cover.size() > o.live_cells)
        return CPS_STOP_CELL_BUDGET;
    f.search = std::make_unique<Subdivision>(f.cover, f.multiplicity, f.precision,
                                             result.statistics_v2, SubdivisionPrecision::TrialLocal);
    return CPS_STOP_NONE;
}

cps_stop_reason SolveState::Impl::structured_step(FactorState &f, const cps_options &o) {
    if (o.heuristic_budget) o.heuristic_budget->quadratic();
    const auto attained = goal_check(o);
    const double iteration_before = f.structured->work();
    auto *metrics = policy.research_scheduling_metrics;
    const auto began = metrics ? std::chrono::steady_clock::now() : std::chrono::steady_clock::time_point{};
    ResearchOnReturn account{[&] {
        if (!metrics) return;
        metrics->numerical_work += std::max(0., f.structured->work() - iteration_before);
        metrics->numerical_ns += research_elapsed(began);
    }};
    const auto precision_increases_before = result.statistics_v2.precision_increases;
    const uint64_t units_before = result.statistics_v2.work_units;
    const auto reason = f.structured->step(o, result.statistics_v2, result.statistics_v3);
    f.research_numerical_precision_increased |= result.statistics_v2.precision_increases > precision_increases_before;
    // A geometric restart may run a nested controller. Its charged units also
    // consume this call's allowance, beyond the outer atomic step already paid.
    result.progress.work_units_this_call = resources::add(
        result.progress.work_units_this_call, result.statistics_v2.work_units - units_before);
    f.precision = std::max(f.precision, f.structured->precision());
    result.work_bits = f.precision;
    result.progress.max_outer_bits = std::max(result.progress.max_outer_bits, f.precision);
    std::vector<Root> cover;
    const bool covered = f.structured->cover(cover);
    if (f.dovetail_engaged) {
        // Records are secular's own: Subdivision's progress must not
        // reset the share that lets it make that progress.
        f.schedule.charge_iteration(f.structured->work() - iteration_before);
        f.schedule.observe(covered ? cover_structure(cover) : DovetailSchedule::Structure{~uint64_t{0}, ~uint64_t{0}, 0},
                           covered ? unresolved_radius_log2(cover) : std::numeric_limits<double>::infinity(),
                           f.structured->last_correction_log2());
    }
    if (covered &&
        (attained(cover, f.precision) || !f.ready ||
         prefer_cover(cover, f.cover, static_cast<slong>(f.precision)))) {
        f.cover = std::move(cover);
        f.ready = true;
    }
    if ((reason != CPS_STOP_PRECISION && reason != CPS_STOP_INPUT) || o.real_only)
        return reason;
#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
    if (policy.research_no_subdivision)
        return research_handoff(reason);
#endif
    // The secular iteration has no global convergence proof, and it
    // has spent the precision budget, or reached the floor an inexact
    // input sets for its residues (CPS_STOP_INPUT: its disks are then
    // within a factor of the degree of what the input allows, and
    // the complete method's certificates can be tighter). Continue
    // with the complete method at the same budget, from the certified
    // cover (a Cauchy disc before the first one), as the classical
    // route does.
    if (!f.ready) {
        ComplexPolynomial p;
        enclose_factor(f, p.value(), static_cast<slong>(f.precision));
        Root disc;
        if (!cauchy_cluster(disc, p.value(), static_cast<slong>(f.precision))) {
            f.blocked = true;
            return reason;
        }
        disc.count *= f.multiplicity;
        if (disc.multiplicity)
            disc.multiplicity *= f.multiplicity;
        f.cover = {std::move(disc)};
        f.ready = true;
    }
    if (f.dovetail_search) {
        // Secular has stopped, so choosing the start once cannot starve
        // Subdivision: keep the dovetailed search unless it has no
        // complete cover yet (mid-split, which on a small share is
        // usual) or secular's final cover is ahead of it, by structure
        // or by the largest region still unresolved. Judged by
        // structure alone, a search whose unresolved cluster is still
        // wider than secular's would be kept, and its precision stop
        // would cost more than a start from secular's cover. History:
        // docs/code-history-20261021.md, solver.cpp 4.
        std::vector<Root> own;
        // This is the final numerical handoff, not a recurring race reset.
        // A delayed search may still hold its first, now coarser cover.
        // Use the maintained one-time choice even for the integrated race;
        // the shared heuristic allowance is neither replaced nor renewed.
        // The delayed event may only have frozen a cover, with no metered
        // arithmetic yet (unmetered bookkeeping may exist). A TrialLocal worker
        // stopped at its precision guard can leave room for ordinary retained
        // refinement under the same caller cap. Both cases use the existing
        // terminal constructor once, from the best certified cover.
        const bool continuable = !policy.research_integrated_dovetail ||
                                 (f.dovetail_search->work() > 0 && f.dovetail_stop != CPS_STOP_PRECISION);
        const bool keep = continuable && (!f.dovetail_engaged ||
                          (f.dovetail_search->complete_cover(own) &&
                           !cover_ahead(f.cover, own)));
        if (keep) {
            f.search = std::move(f.dovetail_search);
            f.search->reopen(o);
            return CPS_STOP_NONE;
        }
        f.dovetail_search.reset();
    }
    const auto start = start_retained_subdivision(f, o);
    if (start != CPS_STOP_NONE) return start;
    if (policy.trace) {
        TraceEvent event{"secular_subdivision_handoff", "transition", f.precision,
                         static_cast<uint64_t>(f.cover.size()), &result.statistics_v3, &result.statistics_v2};
        event.engine = "subdivision";
        event.reason = reason == CPS_STOP_INPUT ? "input-floor" : "precision-budget";
        policy.trace(policy.trace_data, event);
    }
    return CPS_STOP_NONE;
}

std::optional<cps_stop_reason> SolveState::Impl::dovetail_turn(FactorState &f, const cps_options &o) {
    const auto attained = goal_check(o);
    if (f.structured && policy.dovetail_subdivision && !portfolio_engine(policy.candidates) &&
        f.ready && !o.real_only) {
        // Dovetailing (dovetail.hpp): Subdivision runs from the certified
        // cover on its metered share, with persistent state that secular's
        // improvements never reset. At secular's precision stop the same
        // search continues alone (the handoff below moves it).
        if (!f.dovetail_engaged) {
            if (policy.research_integrated_dovetail) {
                using Start = EnginePolicy::DovetailStart;
                if (policy.research_dovetail_start == Start::Never) return std::nullopt;
                // The delayed policy only uses the same exact, unmapped
                // factors as the qualified warm fork. Other factors retain
                // the maintained sequential fallback.
                if (f.kind == FactorState::Ball || f.power != 1) return std::nullopt;
                if (policy.research_dovetail_start == Start::PrecisionCluster &&
                    (!f.research_numerical_precision_increased || !result.statistics_v3.candidate_sweeps ||
                     !std::get<0>(cover_structure(f.cover)))) return std::nullopt;
            }
            if (policy.research_scheduling_metrics) ++policy.research_scheduling_metrics->engagements;
            f.dovetail_engaged = true;
            // Subdivision's first step is a disk test of the factor at its
            // start precision: a Taylor shift and its Graeffe iterations.
            // It is predicted at degree^2 operations each, FLINT's worst
            // case, not by the factor's own block structure: the test
            // shifts to a cluster's centre and squares the result, whose
            // profile can be far worse than the input's. Later steps are
            // predicted from the measured charge of the last. History:
            // docs/code-history-20261021.md, solver.cpp 5.
            uint64_t roots = 0;
            for (const auto &root : f.cover)
                roots += root.count;
            const double degree = static_cast<double>(f.multiplicity ? roots / f.multiplicity : roots);
            double first = (1 + default_graeffe_steps) * degree * degree *
                           operation_weight(static_cast<int64_t>(subdivision_start(o)));
            unsigned share = policy.dovetail_share_log2;
            if (policy.research_integrated_dovetail) share = policy.research_dovetail_share_log2;
            if (policy.research_service_meter)
                first = service_graeffe_price({first / operation_weight(static_cast<int64_t>(subdivision_start(o))), 0},
                                              static_cast<slong>(subdivision_start(o)));
            f.schedule = DovetailSchedule(share, f.structured->work(), first);
            if (policy.research_integrated_dovetail) {
                // Freeze the event's cover immediately, as research_warm_fork
                // did. Admission delays work, not selection of the start cover.
                if (live_cells() + f.cover.size() > o.live_cells)
                    f.dovetail_stop = CPS_STOP_CELL_BUDGET;
                else {
                    const auto began = std::chrono::steady_clock::now();
                    f.dovetail_search = make_research_warm_search(f, o);
                    if (policy.research_scheduling_metrics)
                        policy.research_scheduling_metrics->subdivision_ns += research_elapsed(began);
                }
            }
        }
        // The factor's precision is the budget: the iteration's model is
        // final, and only Subdivision can decide whether the budget
        // suffices, so it runs as the iteration's equal partner there
        // (dovetail.hpp, at_budget). History:
        // docs/code-history-20261021.md, solver.cpp 6.
        if (!policy.research_integrated_dovetail)
        f.schedule.at_budget(f.precision >= o.max_work_bits);
        // Subdivision has spent the precision budget on every pending
        // component, so at this budget its share guarantees nothing more.
        // Secular keeps the factor while its records keep coming. Once it
        // has gone longer without one than all of its work before it (the
        // event that would raise Subdivision's share), the precision
        // budget is what holds the factor. Both states are kept for a
        // raised budget.
        // A capped warm worker does not establish that numerical iteration
        // has exhausted its own precision. Let it reach its terminal handoff
        // (or the shared K/Q boundary), where the better cover is chosen once.
        if (!policy.research_integrated_dovetail &&
            f.dovetail_stop == CPS_STOP_PRECISION && f.schedule.overdue()) {
            f.blocked = true;
            return CPS_STOP_PRECISION;
        }
        if (f.dovetail_stop == CPS_STOP_NONE && f.schedule.subdivision_due() &&
            (!f.dovetail_search || f.dovetail_search->pending())) {
            if (!f.dovetail_search) {
                if (live_cells() + f.cover.size() > o.live_cells)
                    f.dovetail_stop = CPS_STOP_CELL_BUDGET;
                else {
                    if (policy.research_integrated_dovetail)
                        f.dovetail_search = make_research_warm_search(f, o);
                    else
                    f.dovetail_search = std::make_unique<Subdivision>(
                        f.cover, f.multiplicity, subdivision_start(o), result.statistics_v2);
                    if (policy.trace) {
                        TraceEvent event{"dovetail_engaged", "transition", subdivision_start(o),
                                         static_cast<uint64_t>(f.cover.size()), &result.statistics_v3,
                                         &result.statistics_v2};
                        event.engine = "subdivision";
                        policy.trace(policy.trace_data, event);
                    }
                }
            }
            if (f.dovetail_search) {
                f.dovetail_search->research_identity_shift_meter(policy.research_identity_shift_meter);
                f.dovetail_search->research_service_meter(&policy);
                auto *metrics = policy.research_scheduling_metrics;
                const double work_before = metrics ? f.dovetail_search->work() : 0;
                const auto began = metrics ? std::chrono::steady_clock::now() : std::chrono::steady_clock::time_point{};
                bool settled = false;
                const auto settle = [&] {
                    if (settled) return;
                    settled = true;
                    if (!metrics) return;
                    metrics->subdivision_work += std::max(0., f.dovetail_search->work() - work_before);
                    metrics->subdivision_ns += research_elapsed(began);
                    ++metrics->subdivision_turns;
                };
                ResearchOnReturn account{settle};
                TraceScope trace(&policy, "dovetail_subdivision", f.dovetail_search->max_bits(),
                                 f.dovetail_search->pending(), &result.statistics_v3, &result.statistics_v2);
                // A turn costs the caller one work unit, like an iteration
                // step, but Subdivision's own steps are often far smaller.
                // A turn therefore runs as many of them as one average
                // iteration step is worth, while the share still covers
                // them: work units follow the share as metered work does,
                // and the driver's bookkeeping is paid once per turn.
                const double turn = f.schedule.iteration_step();
                double spent = 0;
                bool idle = false;
                cps_stop_reason reason = CPS_STOP_NONE;
                // Each step runs under the share's allowance: a disk test
                // that needs a costlier Taylor shift or Graeffe step than
                // the share covers stops before it, and the turn ends
                // until the share covers that operation. A turn whose
                // first step stopped before computing anything leaves its
                // work unit to the iteration (below): learning that the
                // share does not cover the next operation is no step of
                // the search, and charging it would spend work units on
                // turns that search nothing. History:
                // docs/code-history-20261021.md, solver.cpp 7.
                do {
                    if (o.heuristic_budget && o.heuristic_budget->exhausted()) break;
                    const auto used = live_cells();
                    const double before = f.dovetail_search->work();
                    if (o.heuristic_budget) o.heuristic_budget->linear();
                    reason = f.dovetail_search->step(
                        o, [&](acb_poly_t p, slong bits) { enclose_factor(f, p, bits); },
                        o.live_cells > used ? o.live_cells - used : 0, f.schedule.available()
                        , policy.research_integrated_dovetail ? NewtonAccuracy::PreparedExact : NewtonAccuracy::Legacy
                        );
                    if (o.heuristic_budget)
                        o.heuristic_budget->observe_precision(f.dovetail_search->max_bits());
                    if (f.dovetail_search->suspended()) {
                        f.schedule.wait_for(f.dovetail_search->pending_work());
                        idle = spent == 0 && f.dovetail_search->idle();
                    }
                    if (idle)
                        break;
                    const double work = std::max(1.0, f.dovetail_search->work() - before);
                    if (policy.research_integrated_dovetail)
                        f.schedule.charge_subdivision(work, policy.research_dovetail_growth);
                    else
                    f.schedule.charge_subdivision(work);
                    spent += work;
                    if (policy.research_integrated_dovetail && f.dovetail_search->suspended())
                        f.schedule.wait_for(f.dovetail_search->pending_work());
                    if (f.dovetail_search->suspended())
                        break;
                } while (reason == CPS_STOP_NONE && spent < turn && f.dovetail_search->pending() &&
                         f.schedule.subdivision_due());
                // Diagnostic: the schedule's metered totals (iteration,
                // Subdivision) in the turn's end event, in two spare fields.
                trace.event.preparation_bytes = static_cast<uint64_t>(f.schedule.iteration());
                trace.event.arithmetic_bits = static_cast<uint64_t>(f.schedule.subdivision());
                if (idle && reason == CPS_STOP_NONE) {
                    // Numerical continuation can transfer or discard this
                    // search. Settle its exclusive cost while it is owned here.
                    settle();
                    return structured_step(f, o);
                }
                std::vector<Root> cover;
                if (f.dovetail_search->complete_cover(cover) &&
                    (attained(cover, f.precision) || prefer_cover(cover, f.cover, static_cast<slong>(f.precision))))
                {
                    if (policy.research_scheduling_metrics && attained(cover, f.precision))
                        ++policy.research_scheduling_metrics->fallback_goals;
                    f.cover = std::move(cover);
                }
                else
                    step_kept_covers = true;
                if (reason == CPS_STOP_CANCELLED)
                    return reason;
                // Every component awaits a budget: secular keeps the factor
                // until its own stop hands the search the whole budget, or
                // until the rule above stops both.
                if (reason == CPS_STOP_PRECISION || reason == CPS_STOP_CELL_BUDGET)
                    f.dovetail_stop = reason;
                return CPS_STOP_NONE;
            }
        }
    }
    return std::nullopt;
}

cps_stop_reason SolveState::Impl::candidate_step(FactorState &f, const cps_options &o) {
    if (o.heuristic_budget) o.heuristic_budget->quadratic();
    auto *metrics = policy.research_scheduling_metrics;
    bool candidate_charged = false;
    const uint64_t ns_before = metrics ? metrics->numerical_ns : 0;
    const auto stats_before = result.statistics_v3;
    const auto precision_increases_before = result.statistics_v2.precision_increases;
    const auto began = metrics ? std::chrono::steady_clock::now() : std::chrono::steady_clock::time_point{};
    ResearchOnReturn account{[&] {
        f.research_numerical_precision_increased |= result.statistics_v2.precision_increases > precision_increases_before;
        if (!metrics) return;
        if (!candidate_charged) {
            const auto &after = result.statistics_v3;
            metrics->numerical_work += static_cast<double>(after.evaluation_terms - stats_before.evaluation_terms) +
                                      static_cast<double>(after.root_interactions - stats_before.root_interactions);
        }
        const uint64_t elapsed = research_elapsed(began);
        metrics->numerical_ns += elapsed - std::min(elapsed, metrics->numerical_ns - ns_before);
    }};
    const auto attained = goal_check(o);
    ComplexPolynomial p;
    enclose_factor(f, p.value(), static_cast<slong>(f.precision));
    const auto degree = acb_poly_degree(p.value());
    if (!f.degree_policy_reported) {
        trace_degree_policy(policy, degree);
        f.degree_policy_reported = true;
    }
    const bool adaptive = policy.candidates == CandidateEngine::Adaptive;
    slong support = 0;
    if (adaptive && degree > 3)
        for (slong j = 0; j <= degree; ++j)
            support += !acb_is_zero(p.value()->coeffs + j);
    if (policy.study == StudyMethod::Secular && support > degree / 4)
        prepare_secular_center(f, p.value(), o, result.statistics_v2, policy);
    const bool adaptive_direct = adaptive && degree > heuristics::secular_reconstruction_degree &&
                                 support > 2 && support * 4 <= degree;
    // Routing preferences below (degree above 3, two candidate attempts,
    // the heuristics:: degrees): which factors try the secular path.
    // They decide speed only. A factor sent there keeps the complete
    // method: Subdivision runs beside secular on its share and takes the
    // factor when secular stops (audit fix 4).
    const bool stalled_small_exact = adaptive && f.kind == FactorState::Real && degree > 3 &&
                                     degree <= heuristics::small_complex_degree &&
                                     f.candidates.candidate_attempts >= 2 &&
                                     f.precision > initial_precision;
    const bool adaptive_structured =
        policy.study==StudyMethod::Secular ||
        (adaptive && (adaptive_direct || stalled_small_exact ||
                     (degree > 3 && degree <= heuristics::secular_reconstruction_degree &&
                      f.candidates.candidate_attempts >= 2 &&
                      (support * 4 <= degree ||
                       (degree <= heuristics::small_complex_degree && f.kind == FactorState::Gaussian)))));
    if (policy.backend &&
        (test_structured_engine(policy.candidates) || adaptive_structured ||
         (policy.explicit_source_policy() && policy.study==StudyMethod::None)) &&
        degree > 1 && !acb_contains_zero(p.value()->coeffs)) {
        RootVector seeds(degree);
        ComplexPolynomial shifted_source;
        const bool study_center=policy.study!=StudyMethod::None && f.have_center &&
                                (f.kind == FactorState::Real || f.kind == FactorState::Gaussian);
        if (study_center) {
            if (f.kind == FactorState::Real) acb_poly_set_fmpz_poly(shifted_source.value(),f.centered.value(),static_cast<slong>(f.precision));
            else enclose_gaussian(shifted_source.value(),f.centered_gaussian,static_cast<slong>(f.precision));
        }
        if (initial_candidates(seeds.data(), study_center?shifted_source.value():p.value(), static_cast<slong>(f.precision))) {
            if (o.real_line_seeds && f.kind == FactorState::Real && degree >= heuristics::real_line_seed_min_degree)
                seed_from_real_line(seeds, f, study_center, o, policy, result.statistics_v3);
            if (adaptive && !adaptive_direct && degree >= 128 && f.candidates.centers &&
                f.candidates.centers->size() == degree) {
                bool reusable = true;
                for (slong i = 0; i < degree; ++i) {
                    reusable = reusable && acb_is_finite(f.candidates.centers->data() + i);
                    for (slong j = 0; j < i; ++j)
                        reusable = reusable && !acb_equal(f.candidates.centers->data() + i,
                                                          f.candidates.centers->data() + j);
                }
                if (reusable)
                    for (slong i = 0; i < degree; ++i)
                        acb_get_mid(seeds.data() + i, f.candidates.centers->data() + i);
            }
            PolynomialSourceOptions source_options;
            source_options.sparse_coefficients = support > 2 && support * 4 <= degree;
#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
            source_options.evaluation = policy.candidates == CandidateEngine::Series
                                            ? EvaluationMethod::TaylorCache
                                            : EvaluationMethod::Coefficients;
#endif
            source_options.candidate_direct = direct_engine(policy.candidates) ||
                                              adaptive_direct ||
                                              (adaptive && degree <= heuristics::small_complex_degree &&
                                               f.kind == FactorState::Gaussian);
            source_options.input = result.progress.input_uncertainty
                                       ? InputSemantics::CoefficientEnclosures
                                       : InputSemantics::Exact;
            if (stalled_small_exact ||
                (f.kind == FactorState::Gaussian && degree <= heuristics::small_complex_degree))
                source_options.exact_taylor = owned_exact_taylor(f, o.exact_storage_limit(), o.allow_exponent_expansion, o.exact_operations
                    , policy.research_service_meter ? &source_options.exact_taylor_work : nullptr
                    );
            source_options.exact_precision = static_cast<slong>(std::min<uint64_t>(o.max_work_bits, static_cast<uint64_t>(WORD_MAX)));
            // The exact polynomial this source evaluates (the centred
            // factor under a centred study), one shared copy for its
            // encloser and, in the research build, for the exact value at a
            // node that residue regeneration takes at the precision budget
            // (regenerate_residues, exact_residues_at_budget).
            const auto shared = shared_factor_polynomial(f, study_center);
            auto encloser = shared ? shared->encloser() : owned_encloser(f, input);
            if (study_center)
                source_options.exact_taylor = {};
            if (study_center)
                source_options.exact_taylor_work = {};
#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
            if (policy.exact_residues_at_budget && shared) {
                source_options.exact_integer = shared->integer;
                source_options.exact_gaussian = shared->gaussian;
            }
#endif
            if (policy.centred_origin_basis && study_center)
                attach_origin(source_options, f);
            auto source = std::make_unique<PolynomialSource>(std::move(encloser), degree,
                                                             acb_poly_is_real(study_center ? shifted_source.value() : p.value()), source_options);
            auto continuation = policy;
            continuation.snapshot_factor = f.id;
            continuation.snapshot_power = f.power;
            if (study_center) continuation.snapshot_offset = f.center;
            f.structured = std::make_unique<SecularRefinement>(
                std::move(source), seeds.data(), f.multiplicity, static_cast<slong>(f.precision), continuation,
                policy.research_integrated_dovetail ? solve_local_model_priced :
#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
                    policy.research_no_subdivision ? solve_local_model_numerical :
#endif
                solve_local_model,
                study_center ? RootCoordinates(f.center) : RootCoordinates{});
            if (metrics) {
                const auto &after = result.statistics_v3;
                metrics->numerical_work += static_cast<double>(after.evaluation_terms - stats_before.evaluation_terms) +
                                          static_cast<double>(after.root_interactions - stats_before.root_interactions);
                candidate_charged = true;
            }
            return structured_step(f, o);
        }
    }
    if (!f.ready && !o.real_only) {
        Root fallback;
        if (!cauchy_cluster(fallback, p.value(), static_cast<slong>(f.precision))) {
            if (f.precision == o.max_work_bits) {
                f.blocked = true;
                return CPS_STOP_PRECISION;
            }
            f.precision = std::min(o.max_work_bits, f.precision * 2);
            return CPS_STOP_NONE;
        }
        fallback.count *= f.multiplicity;
        if (fallback.multiplicity)
            fallback.multiplicity *= f.multiplicity;
        f.cover = {std::move(fallback)};
        f.ready = true;
    }
    std::vector<Root> candidates;
    bool isolated = false;
    ComplexPolynomial candidate_polynomial;
    ComplexBall translation;
    bool translated = false;
    auto engine = policy.candidates;
    if (engine == CandidateEngine::Adaptive) {
        bool dense_translation = (f.kind == FactorState::Real || (policy.complex_centering && f.kind == FactorState::Gaussian)) && !f.center_rejected &&
                                       degree >= heuristics::centering_min_degree &&
                                       support * 4 > degree && !acb_is_zero(p.value()->coeffs + degree - 1);
        if (dense_translation && !f.have_center &&
            !(f.kind == FactorState::Real ? center_height_allows_gain(f.integer.value()) :
              center_height_allows_gain(f.gaussian))) {
            f.center_rejected=true;
            dense_translation=false;
            PreparationDiagnostic diagnostic;
            diagnostic.reason="centering-height-gain-impossible";
            observe_preparation(policy,"exact-centering",diagnostic,o.exact_budget,
                                std::chrono::steady_clock::now(),resources::nonnegative(degree));
        }
        engine = dense_translation ? CandidateEngine::Centered : CandidateEngine::Aberth;
    }
#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
    if (engine == CandidateEngine::Hybrid)
        engine = acb_poly_degree(p.value()) >= 128 ? CandidateEngine::Aberth
                 : acb_poly_degree(p.value()) >= 8 && f.candidates.candidate_attempts >= 2
                     ? CandidateEngine::Centered
                     : CandidateEngine::Legacy;
#endif
    if (engine == CandidateEngine::Centered && f.kind != FactorState::Ball &&
        acb_poly_degree(p.value()) > 1) {
        if (f.center_precision != f.precision) {
            f.center_precision = f.precision;
            GaussianRational proposed;
            ComplexBall target, old_center;
            arb_set_fmpq(acb_realref(old_center.value()), f.center.real.value(),
                         static_cast<slong>(f.precision));
            arb_set_fmpq(acb_imagref(old_center.value()), f.center.imag.value(),
                         static_cast<slong>(f.precision));
            if (!f.have_center) {
                const auto n = acb_poly_degree(p.value());
                acb_div(target.value(), p.value()->coeffs + n - 1, p.value()->coeffs + n, 64);
                acb_div_ui(target.value(), target.value(), static_cast<ulong>(n), 64);
                acb_neg(target.value(), target.value());
                // A low-denominator dyadic centre can remove most global
                // cancellation while keeping the exact translation cheap.
                // Its accuracy is heuristic; the shifted coefficients are
                // exact and the final enclosures are translated outward.
                round_center(target.value(), heuristics::initial_center_bits);
            } else if (f.candidates.unresolved_count >= 2 &&
                       f.candidates.unresolved_count < acb_poly_degree(p.value()) &&
                       acb_is_finite(f.candidates.unresolved_center.value())) {
                acb_add(target.value(), old_center.value(), f.candidates.unresolved_center.value(),
                        static_cast<slong>(
                            std::min<uint64_t>(heuristics::global_centering_bits, f.precision / 2)));
            } else
                acb_set(target.value(), old_center.value());
            if (f.kind == FactorState::Real)
                arb_zero(acb_imagref(target.value()));
            const auto owned = f.kind == FactorState::Real ? center_owned_bits(f.integer.value()) : center_owned_bits(f.gaussian);
            const auto center_limit = o.allow_exponent_expansion ? o.exact_storage_limit() : std::min(o.exact_storage_limit(), owned);
            if (center_materialization_fits(target.value(), center_limit) &&
                o.exact_effort_fits(resources::multiply(resources::nonnegative(degree + 1), resources::nonnegative(degree + 1)))) {
                arf_get_fmpq(proposed.real.value(), arb_midref(acb_realref(target.value())));
                arf_get_fmpq(proposed.imag.value(), arb_midref(acb_imagref(target.value())));
                if (!f.have_center || !fmpq_equal(proposed.real.value(), f.center.real.value()) ||
                    !fmpq_equal(proposed.imag.value(), f.center.imag.value())) {
                    IntegerPolynomial shifted;
                    GaussianPolynomial shifted_gaussian;
                    StageTimer timer(result.statistics_v2.taylor_shift_ns);
                    PreparationDiagnostic diagnostic;
                    const auto center_started = std::chrono::steady_clock::now();
                    TraceScope construction(&policy,"exact-centering-construction",f.precision,resources::nonnegative(degree));
                    // The first centring of a real factor is kept only if the height rule below
                    // accepts it. When the translate's constant term alone shows that it cannot
                    // (translated_height_lower_bound; the same prediction as in prepare_secular_center's
                    // first centring), the translate is not built: the outcome is the refusal the
                    // rule would give after construction.
                    bool predicted_refusal = false;
                    if (adaptive && !f.have_center && f.kind == FactorState::Real)
                        if (const auto floor_height = translated_height_lower_bound(f.integer.value(), proposed.real.value())) {
                            const slong before = FLINT_ABS(fmpz_poly_max_bits(f.integer.value())) -
                                                 static_cast<slong>(fmpz_bits(fmpz_poly_lead(f.integer.value())));
                            predicted_refusal = before - *floor_height <
                                                static_cast<slong>(FLINT_BIT_COUNT(static_cast<ulong>(degree)));
                        }
                    bool success = !predicted_refusal &&
                        (f.kind == FactorState::Real
                            ? center_integer(shifted.value(), f.integer.value(), proposed.real.value(),
                                             o.exact_storage_limit(), policy.trace ? &diagnostic : nullptr)
                            : center_gaussian(shifted_gaussian, f.gaussian, proposed, o.exact_storage_limit()));
                    const bool converted = success;
                    if (adaptive && success && !f.have_center && f.kind == FactorState::Real) {
                        // A centroid shift can destroy useful radial scaling
                        // (for example, geometrically spaced roots). Accept
                        // this heuristic only if the exact monic coefficient
                        // height falls by about log2(degree) bits. This is a
                        // candidate-method choice, not a root certificate.
                        const auto height = [](const fmpz_poly_t polynomial) {
                            const slong maximum = FLINT_ABS(fmpz_poly_max_bits(polynomial));
                            return maximum - static_cast<slong>(fmpz_bits(fmpz_poly_lead(polynomial)));
                        };
                        const auto before = height(f.integer.value()), after = height(shifted.value());
                        success = before - after >= static_cast<slong>(FLINT_BIT_COUNT(static_cast<ulong>(degree)));
                    }
                    if (adaptive && success && !f.have_center && f.kind == FactorState::Gaussian)
                        success = gaussian_center_improves_height(f.gaussian, shifted_gaussian);
                    if (converted)
                        diagnostic.reason = success ? "centering-selected" : "centering-no-height-improvement";
                    else if (predicted_refusal)
                        diagnostic.reason = "centering-constant-term-bound";
                    else if (f.kind != FactorState::Real)
                        diagnostic.reason = "gaussian-centering-unavailable";
                    diagnostic.accepted = success;
                    observe_preparation(policy, "exact-centering", diagnostic, o.exact_budget,
                                        center_started, static_cast<uint64_t>(degree));
                    if (adaptive && !success && !f.have_center)
                        f.center_rejected = true;
                    if (success) {
                        ComplexBall delta;
                        Rational difference;
                        fmpq_sub(difference.value(), f.center.real.value(), proposed.real.value());
                        arb_set_fmpq(acb_realref(delta.value()), difference.value(),
                                     static_cast<slong>(f.precision));
                        fmpq_sub(difference.value(), f.center.imag.value(), proposed.imag.value());
                        arb_set_fmpq(acb_imagref(delta.value()), difference.value(),
                                     static_cast<slong>(f.precision));
                        if (f.candidates.centers)
                            for (slong i = 0; i < f.candidates.centers->size(); ++i) {
                                acb_add(f.candidates.centers->data() + i, f.candidates.centers->data() + i,
                                        delta.value(), static_cast<slong>(f.precision));
                                acb_get_mid(f.candidates.centers->data() + i,
                                            f.candidates.centers->data() + i);
                            }
                        const bool first = !f.have_center;
                        f.center = std::move(proposed);
                        f.centered = std::move(shifted);
                        f.centered_gaussian = std::move(shifted_gaussian);
                        f.have_center = true;
                        if (f.candidates.study) f.candidates.study->coordinates_changed();
                        ++result.statistics_v2.taylor_shifts;
                        if (first && policy.centering_refinement)
                            refine_center(f, o, result.statistics_v2, policy);
                    }
                }
            }
        }
        if (f.have_center) {
            if (f.kind == FactorState::Real)
                acb_poly_set_fmpz_poly(candidate_polynomial.value(), f.centered.value(),
                                       static_cast<slong>(f.precision));
            else
                enclose_gaussian(candidate_polynomial.value(), f.centered_gaussian,
                                 static_cast<slong>(f.precision));
            arb_set_fmpq(acb_realref(translation.value()), f.center.real.value(),
                         static_cast<slong>(f.precision));
            arb_set_fmpq(acb_imagref(translation.value()), f.center.imag.value(),
                         static_cast<slong>(f.precision));
            translated = true;
        }
    }
    {
        StageTimer timer(result.statistics.isolation_ns);
        ++result.statistics.isolation_passes;
        if (policy.backend || acb_poly_degree(p.value()) == 1) {
            std::optional<EnginePolicy> observation;
            if (policy.snapshot || policy.trace) {
                observation = policy;
                observation->snapshot_factor = f.id;
                observation->snapshot_power = f.power;
                if (translated) observation->snapshot_offset = f.center;
            }
            // The certificate attempt (isolate_complex) reads the factor
            // at any precision it needs, in this pass's coordinates.
            PolynomialEncloser factor_source;
            if (policy.newton_point_refinement && !f.candidates.newton_attempted) {
                // The same polynomial as this pass's (enclose_factor, or
                // the centred factor), without the pass's bookkeeping.
                factor_source = [&](acb_poly_t out, slong bits) {
                    if (translated && f.kind == FactorState::Real)
                        acb_poly_set_fmpz_poly(out, f.centered.value(), bits);
                    else if (translated)
                        enclose_gaussian(out, f.centered_gaussian, bits);
                    else if (f.kind == FactorState::Ball)
                        enclose_input(out, bits);
                    else if (f.kind == FactorState::Real)
                        acb_poly_set_fmpz_poly(out, f.integer.value(), bits);
                    else
                        enclose_gaussian(out, f.gaussian, bits);
                };
            }
            isolated =
                isolate_complex(candidates, translated ? candidate_polynomial.value() : p.value(),
                                f.multiplicity, f.candidates, o, static_cast<slong>(f.precision),
                                o.goal != CPS_ROOT_ACCURACY || f.precision == o.max_work_bits,
                                result.statistics_v2, result.statistics_v3, engine,
                                hybrid_engine(policy.candidates) ||
                                    (adaptive && degree >= 128 && support == 2),
                                observation ? &*observation : &policy,
                                factor_source ? &factor_source : nullptr);
            if (isolated && f.candidates.cover_bits > static_cast<slong>(result.work_bits)) {
                // The attempt certified at a precision above this pass's:
                // the precision reported is the precision used.
                result.work_bits = static_cast<uint64_t>(f.candidates.cover_bits);
                result.progress.max_outer_bits = std::max(result.progress.max_outer_bits, result.work_bits);
            }
        }
    }
    if (f.candidates.budget_stop) {
        // No whole packet ran (the solve's iteration budget was spent
        // before or during it): nothing was validated, and the factor
        // waits for the budget (iteration_wait) until the solve stops.
        f.candidates.budget_stop = false;
        return CPS_STOP_ITERATION_BUDGET;
    }
    if (translated)
        for (auto &root : candidates) {
            acb_add(root.box.value(), root.box.value(), translation.value(), static_cast<slong>(f.precision));
            root.realness = arb_is_zero(acb_imagref(root.box.value()))          ? CPS_REAL
                            : !arb_contains_zero(acb_imagref(root.box.value())) ? CPS_NONREAL
                                                                              : CPS_REAL_UNKNOWN;
            root.certificate += ";exact-rational-variable-shift";
            root.publish(static_cast<slong>(f.precision));
        }
    uint64_t count = 0;
    for (const auto &root : candidates)
        count += root.count;
    const uint64_t expected =
        o.real_only ? f.real_count : static_cast<uint64_t>(acb_poly_degree(p.value())) * f.multiplicity;
    if (isolated && count == expected && disjoint(candidates)) {
        if (!f.ready || attained(candidates, f.precision) ||
            prefer_cover(candidates, f.cover, static_cast<slong>(f.precision)))
            f.cover = std::move(candidates);
        f.ready = true;
        if (attained(f.cover, f.precision))
            return CPS_STOP_NONE;
        // Complex factors continue with independently scheduled certified components.
        if (!o.real_only && acb_poly_degree(p.value()) > 1 && f.precision == o.max_work_bits) {
#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
            if (policy.research_no_subdivision)
                return research_handoff(CPS_STOP_PRECISION);
#endif
            if (f.dovetail_search) {f.search=std::move(f.dovetail_search);f.search->reopen(o);return CPS_STOP_NONE;}
            if (live_cells() + f.cover.size() > o.live_cells)
                return CPS_STOP_CELL_BUDGET;
            candidate_handoff(f, o);
            return CPS_STOP_NONE;
        }
    }
    if (!policy.backend && !o.real_only && f.ready) {
#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
        if (policy.research_no_subdivision)
            return research_handoff(CPS_STOP_PRECISION);
#endif
        if (live_cells() + f.cover.size() > o.live_cells)
            return CPS_STOP_CELL_BUDGET;
        f.search = std::make_unique<Subdivision>(f.cover, f.multiplicity, f.precision, result.statistics_v2);
    } else if (!o.real_only && f.candidates.repeat) {
        return CPS_STOP_NONE;
    } else if (f.precision < o.max_work_bits) {
        f.precision = std::min(o.max_work_bits, f.precision * 2);
        ++result.statistics_v2.precision_increases;
    } else if (!o.real_only && acb_poly_degree(p.value()) > 1) {
#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
        if (policy.research_no_subdivision)
            return research_handoff(CPS_STOP_PRECISION);
#endif
        if (f.dovetail_search) {f.search=std::move(f.dovetail_search);f.search->reopen(o);return CPS_STOP_NONE;}
        if (live_cells() + f.cover.size() > o.live_cells)
            return CPS_STOP_CELL_BUDGET;
        candidate_handoff(f, o);
    } else
        f.blocked = true;
    return CPS_STOP_NONE;
}
SolveState::SolveState(const DensePolynomial &p, EnginePolicy policy)
    : impl_(std::make_unique<Impl>(p, policy)) {}
SolveState::SolveState(DensePolynomial &&p, EnginePolicy policy)
    : impl_(std::make_unique<Impl>(std::move(p), policy)) {}
SolveState::SolveState(const cps_secular &p, EnginePolicy policy) : impl_(std::make_unique<Impl>(p, policy)) {}
SolveState::SolveState(const cps_program &p, EnginePolicy policy) : impl_(std::make_unique<Impl>(p, policy)) {}
SolveState::~SolveState() = default;
cps_result SolveState::advance(const cps_options *options) { return impl_->advance(options); }
void SolveState::hint_real_line_seeds(RealLineSeedHint &&hint) { impl_->real_line_hint = std::move(hint); }
std::optional<PreparedStructure> SolveState::prepared_structure() const {
    if (!impl_->preprocessed) return std::nullopt;
    PreparedStructure out;
    out.input_degree = impl_->result.degree;
    // The exact zero roots set apart before the factors (preprocess); append adds any a factor holds.
    out.zero_order = impl_->zero_order;
    for (const auto &factor : impl_->factors) {
        PreparedFactorStructure shape;
        shape.multiplicity = factor->multiplicity;
        shape.power = factor->power;
        shape.squarefree_proved = factor->kind != FactorState::Ball || factor->squarefree_proved;
        if (factor->kind == FactorState::Real) {
            shape.representation = PreparedRepresentation::Integer;
            shape.shape = input_structure(factor->integer);
        } else if (factor->kind == FactorState::Gaussian) {
            shape.representation = PreparedRepresentation::Gaussian;
            shape.shape = input_structure(std::span(factor->gaussian));
        } else {
            // This owned input may already have undergone power reduction;
            // the original source metadata is not the right snapshot here.
            shape.shape = input_structure(impl_->input.coefficients);
        }
        out.append(std::move(shape));
    }
    return out;
}
cps_result solve(const DensePolynomial &polynomial, const cps_options &options) {
    SessionState state(polynomial);
    return state.advance(&options);
}
cps_result solve(DensePolynomial &&polynomial, const cps_options &options) {
    SessionState state(std::move(polynomial));
    return state.advance(&options);
}
cps_result solve(const cps_polynomial &borrowed, const cps_options &options) {
    SessionState state(borrowed);
    return state.advance(&options);
}
} // namespace certpolysolve::detail
