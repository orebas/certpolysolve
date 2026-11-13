// The classical route's isolation of one factor's roots at one precision: isolate_complex turns one pass's candidates
// into certified regions, or reports that this precision does not settle them; sort_roots orders a result's regions.
//
// In order: an exact factor x^k is split off as one region at 0 (no numerical deflation); a linear factor is solved by
// ball division, a binomial a + b x^n (with analytic_binomial) by an nth root of -a/b times the nth roots of unity
// (binomial_enclosures). Otherwise the certificate attempt (machine_newton_cover, once, under newton_point_refinement)
// may certify every root; if not, a packet of at most heuristics::classical_packet_sweeps sweeps runs from the
// retained centres, if any (StudyCandidates::advance in the maintained recipes: the Aberth iteration, binary64 first,
// unless the study is StudyMethod::Secular; else candidate_packet). FLINT's _acb_poly_validate_roots gives each
// candidate z a disk of radius n |p(z)/p'(z)|, which holds a root; all n pairwise disjoint certify the vector, one
// root each. Short of that, recover_clusters (with allow_clusters) tries a cover by disk counts, and
// CandidateState::repeat asks for another packet at this precision before a higher one, as the packets' progress
// within max(heuristics::initial_sweep_allowance, 2n) sweeps and the probe candidate_arithmetic_adequate decide.
// Regions are classified (classify: real where acb_poly_validate_real_roots proves it) and published (Root::publish).
//
// State: none of its own. A factor's CandidateState (classical_isolation.hpp), held in FactorState::candidates
// (factor_state.hpp), carries its centres, engines and packet accounting (packets, charged, budget_stop) across passes;
// factor_state.cpp and solver.cpp move the centres when they recentre the factor.
// Entry points: isolate_complex, called by solver.cpp's SolveState::Impl::factor_step (which translates centred regions
// back, then accepts the cover, repeats or doubles the precision, or hands the factor to Subdivision); sort_roots (by
// the midpoints' real, then imaginary parts), by solver.cpp's advance, native_advance (solver_native.cpp) and
// program_advance (solver_program.cpp).
// Not here: the attempt's body (piecewise_roots.cpp), packet iterations (study_candidates.cpp, candidates.cpp),
// cluster counts (clustering.cpp, counting.cpp), centring (centering.cpp, factor_state.cpp).
// Design notes: docs/mathematics.md ("Individual root certificates", "Local counts and candidate-guided recovery"),
// docs/max-iterations-note-20261001.md (the iteration budget and the packet it cuts short).
#include "methods/classical/classical_isolation.hpp"
#include "certify/attempts/newton_cover.hpp"
#include "core/engine_policy.hpp"
#include "core/engine_statistics.hpp"
#include "certify/counting.hpp"
#include "core/telemetry.hpp"

namespace certpolysolve::detail {
namespace {
bool binomial_enclosures(acb_ptr roots, const acb_poly_t p, slong precision, const cps_options &options) {
    const auto n = acb_poly_degree(p);
    if (n < 2)
        return false;
    for (slong i = 1; i < n; ++i)
        if (!acb_is_zero(p->coeffs + i))
            return false;
    ComplexBall ratio, principal, rotation;
    acb_div(ratio.value(), p->coeffs, p->coeffs + n, precision);
    acb_neg(ratio.value(), ratio.value());
    if (!acb_is_finite(ratio.value()) || acb_contains_zero(ratio.value()))
        return false;
    const bool rotate = arb_is_negative(acb_realref(ratio.value()));
    acb_set(principal.value(), ratio.value());
    if (rotate)
        acb_neg(principal.value(), principal.value());
    acb_root_ui(principal.value(), principal.value(), static_cast<ulong>(n), precision);
    if (rotate) {
        acb_unit_root(rotation.value(), static_cast<ulong>(2 * n), precision);
        acb_mul(principal.value(), principal.value(), rotation.value(), precision);
    }
    _acb_vec_unit_roots(roots, n, n, precision);
    for (slong i = 0; i < n; ++i) {
        if (options.cancelled())
            return false;
        acb_mul(roots + i, roots + i, principal.value(), precision);
        if (!acb_is_finite(roots + i))
            return false;
        if (acb_is_real(ratio.value()) && ((!rotate && (i == 0 || (n % 2 == 0 && i == n / 2))) ||
                                         (rotate && n % 2 == 1 && i == (n - 1) / 2)))
            arb_zero(acb_imagref(roots + i));
    }
    // For each nonzero ratio, one nth root times all nth roots of unity
    // exhausts exactly n simple roots. Final outward disjointness is checked
    // by the ordinary cover ledger; no numerical agreement establishes counts.
    return true;
}
bool candidate_arithmetic_adequate(const acb_poly_t polynomial, ConstRootView points, slong first_unresolved,
                                   slong precision) {
    // Stagnation is not itself evidence that more bits will help. Sample
    // unresolved Newton evaluations with outward arithmetic. Well-resolved
    // values/derivatives justify more bounded iteration at this precision.
    const slong remaining = points.size() - first_unresolved;
    if (remaining <= 0)
        return false;
    constexpr slong samples = 8, reliable_bits = 16;
    ComplexBall center, value, derivative;
    for (slong j = 0; j < std::min(samples, remaining); ++j) {
        const slong i = first_unresolved + j * remaining / std::min(samples, remaining);
        acb_get_mid(center.value(), points[i]);
        acb_poly_evaluate2(value.value(), derivative.value(), polynomial, center.value(), precision);
        if (!acb_is_finite(value.value()) || !acb_is_finite(derivative.value()) ||
            acb_contains_zero(derivative.value()) ||
            (!acb_is_zero(value.value()) && acb_rel_accuracy_bits(value.value()) < reliable_bits) ||
            acb_rel_accuracy_bits(derivative.value()) < reliable_bits)
            return false;
    }
    return true;
}


void classify(Root &root, bool real_polynomial, bool real_roots_validated) {
    auto imag = acb_imagref(root.box.value());
    if (real_polynomial && real_roots_validated && arb_contains_zero(imag))
        arb_zero(imag);
    if (arb_is_zero(imag))
        root.realness = CPS_REAL;
    else if (!arb_contains_zero(imag))
        root.realness = CPS_NONREAL;
}
} // namespace

bool isolate_complex(std::vector<Root> &result, const acb_poly_t polynomial, uint64_t multiplicity,
                     CandidateState &state, const cps_options &options, slong precision, bool allow_clusters,
                     certpolysolve::detail::SearchStatistics &statistics, certpolysolve::detail::CandidateStatistics &candidate_statistics,
                     CandidateEngine engine, bool analytic_binomial, const EnginePolicy *policy,
                     const PolynomialEncloser *exact) {
    const slong degree = acb_poly_degree(polynomial);
    if (degree < 1)
        return true;
    // Exact zero deflation is safe; numerical deflation is not used.
    slong zeros = 0;
    while (zeros < degree && acb_is_zero(polynomial->coeffs + zeros))
        ++zeros;
    if (zeros) {
        Root root;
        root.count = root.multiplicity = static_cast<uint64_t>(zeros) * multiplicity;
        root.realness = CPS_REAL;
        root.certificate = "exact-zero-factor";
        root.publish(precision);
        result.push_back(std::move(root));
    }
    if (zeros == degree)
        return true;
    ComplexPolynomial reduced;
    acb_poly_shift_right(reduced.value(), polynomial, zeros);
    const slong reduced_degree = degree - zeros;
    RootVector roots(reduced_degree);
    bool binomial = false;
    if (analytic_binomial && reduced_degree > 1) {
        StageTimer timer(candidate_statistics.validation_ns);
        binomial = binomial_enclosures(roots.data(), reduced.value(), precision, options);
        if (binomial)
            ++candidate_statistics.validation_calls;
    }
    if (reduced_degree == 1) {
        acb_div(roots.data(), reduced.value()->coeffs, reduced.value()->coeffs + 1, precision);
        acb_neg(roots.data(), roots.data());
        if (!acb_is_finite(roots.data()))
            return false;
    } else if (!binomial) {
        if (policy && policy->newton_point_refinement && !state.newton_attempted) {
            state.newton_attempted = true;
            // The attempt may raise its own precision: it reads the factor
            // itself when the caller supplies its encloser (`polynomial` is
            // that encloser's value at `precision`), and this call's
            // enclosure otherwise, whose rounding no precision reduces.
            PolynomialSource source(exact ? PolynomialEncloser([exact, zeros](acb_poly_t p, slong bits) {
                                                ComplexPolynomial whole;
                                                (*exact)(whole.value(), bits);
                                                acb_poly_shift_right(p, whole.value(), zeros);
                                            })
                                          : PolynomialEncloser([&reduced](acb_poly_t p, slong) {
                                                acb_poly_set(p, reduced.value());
                                            }),
                                    reduced_degree, acb_poly_is_real(reduced.value()));
            std::vector<Root> checked;
            // Keep finite machine proposals when the shortcut's proof fails.
            // Ordinary refinement still owes the full original-input proof.
            if (machine_newton_cover(checked, source, nullptr,
                                      options, precision,
                                      policy->candidate_allow_underflow, candidate_statistics,policy->coefficient_newton_bound,policy,
                                      &state.centers, nullptr, &state.cover_bits) &&
                (!zeros || std::none_of(checked.begin(), checked.end(), [](const auto &r) {
                    return acb_contains_zero(r.box.value());
                }))) {
                auto retained = std::make_unique<RootVector>(reduced_degree);
                for (slong i = 0; i < reduced_degree; ++i)
                    acb_get_mid((*retained)[i], checked[static_cast<size_t>(i)].box.value());
                state.centers = std::move(retained);
                for (auto &root : checked) root.count = root.multiplicity = multiplicity;
                result.insert(result.end(), std::make_move_iterator(checked.begin()),
                              std::make_move_iterator(checked.end()));
                return true;
            }
        }
        // Restarting a short iteration from scratch at every precision can
        // NEVER resolve a tight cluster, even after raising precision thousands
        // of bits. Preserve finite, distinct candidates and grow the default
        // iteration allowance with precision. These are heuristic decisions;
        // only the independent inclusion/coverage predicates certify results.
        acb_srcptr initial =
            state.centers && state.centers->size() == reduced_degree ? state.centers->data() : nullptr;
        slong isolated;
        state.repeat = false;
#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
        if (engine == CandidateEngine::DoubleSeed && !state.seeded && !initial) {
            StageTimer timer(candidate_statistics.candidate_ns);
            ++candidate_statistics.double_seed_calls;
            auto seeds = std::make_unique<RootVector>(reduced_degree);
            if (binary64_seeds(seeds->data(), reduced.value(), options, candidate_statistics)) {
                state.centers = std::move(seeds);
                initial = state.centers->data();
            }
            state.seeded = true;
        }
#endif
        // A packet's length: the engine's own; then what the solve's
        // iteration budget leaves of it.
        const unsigned packet_length = heuristics::classical_packet_sweeps;
        const unsigned packet_size = options.packet_sweeps(packet_length, candidate_statistics);
        if (!packet_size) {
            state.budget_stop = true;
            return false;
        }
        {
            TraceScope trace(policy, "native_candidates", static_cast<uint64_t>(precision),
                             static_cast<uint64_t>(reduced_degree), &candidate_statistics, &statistics);
            ++state.candidate_attempts;
            if (state.last_goal_bits != options.bits) {
                state.inactive.clear();
                state.last_goal_bits = options.bits;
            }
            if (state.precision != precision) {
                state.precision = precision;
                state.packets = 0;
                state.charged = 0;
                state.correction_log2 = std::numeric_limits<double>::infinity();
            }
            CandidateProgress progress;
            if (policy && (policy->study != StudyMethod::None ||
                           (policy->study_machine && policy->candidates == CandidateEngine::Adaptive))) {
                if (!state.study) {
                    auto arithmetic_policy = *policy;
                    // The surrounding C2 policy retains source selection and
                    // recovery. This owner supplies arithmetic only for the
                    // ordinary Aberth packet already selected at this point.
                    if (arithmetic_policy.study == StudyMethod::None)
                        arithmetic_policy.study = StudyMethod::Aberth;
                    state.study=std::make_unique<StudyCandidates>(arithmetic_policy);
                }
                progress=state.study->advance(roots.data(),reduced.value(),initial,precision,packet_size,
                    options,candidate_statistics,state.inactive,
                    reduced_degree>=heuristics::interaction_tree_degree);
#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
            } else if (engine == CandidateEngine::Portfolio && policy) {
                if (!state.portfolio)
                    state.portfolio=std::make_unique<CandidatePortfolio>(*policy);
                const auto packet=state.portfolio->advance(roots.data(),reduced.value(),initial,precision,options,candidate_statistics);
                progress=packet.candidates; state.repeat=packet.repeat;
#endif
            } else progress = candidate_packet(roots.data(), reduced.value(), initial, precision, packet_size, engine, options,
                                 candidate_statistics, state.inactive,
                                 policy && reduced_degree >= heuristics::interaction_tree_degree);
            if (progress.cancelled)
                return false;
            // The solve's iteration budget cut the packet short (the budget is
            // spent), and the packet did not end by itself: no verdict on this
            // precision, and the caller stops the solve for the budget.
            if (packet_size < packet_length && !options.sweeps_left(candidate_statistics) && progress.finite &&
                !progress.converged &&
                progress.continuation == CandidateContinuation::Ordinary && !portfolio_engine(engine)) {
                state.budget_stop = true;
                return false;
            }
            ++state.packets;
            state.charged += packet_length;
            const uint64_t allowance = std::max<uint64_t>(heuristics::initial_sweep_allowance,
                                                          2 * static_cast<uint64_t>(reduced_degree));
            const bool continue_precision =
                progress.continuation == CandidateContinuation::NextPacketAtRequestedPrecision;
            const bool stage_stopped = progress.continuation != CandidateContinuation::Ordinary;
            if (!portfolio_engine(engine)) {
                const uint64_t charged = state.charged;
                state.repeat = progress.finite &&
                    (continue_precision ? charged <= allowance && allowance - charged >= packet_length
                             : !stage_stopped && !progress.converged && charged < allowance &&
                               (state.packets < 4 || progress.correction_log2 < state.correction_log2 - 0.5));
            }
            if (!state.repeat && !progress.converged && !stage_stopped)
                ++candidate_statistics.candidate_stagnations;
            state.correction_log2 = progress.correction_log2;
            isolated = 0;
            if (progress.finite) {
                if (policy && policy->snapshot) {
                    NumericalSnapshot event;
                    event.stage = "ordinary_before_validation";
                    event.factor = policy->snapshot_factor;
                    event.power = policy->snapshot_power;
                    event.offset = policy->snapshot_offset;
                    event.multiplicity = multiplicity;
                    event.zero_count = static_cast<uint64_t>(zeros) * multiplicity;
                    event.source_precision = static_cast<uint64_t>(precision);
                    event.coordinate_precision.assign(static_cast<size_t>(reduced_degree), precision);
                    event.inactive = state.inactive;
                    for (slong i = 0; i < reduced_degree; ++i) {
                        event.points.emplace_back();
                        acb_get_mid(event.points.back().value(), roots[i]);
                    }
                    policy->snapshot(std::move(event));
                }
                // Save proposals BEFORE validation changes radii and ordering.
                state.centers = std::make_unique<RootVector>(reduced_degree);
                _acb_vec_set(state.centers->data(), roots.data(), reduced_degree);
                if (policy && policy->defer_classical_validation && state.repeat &&
                    progress.continuation == CandidateContinuation::NextPacketAtRequestedPrecision) {
                    // Active classical candidates explicitly request another
                    // packet at this precision. Preserve their exact ordering
                    // and certified inactivity; validate when that work stops,
                    // before any precision decision or successful return.
                    if (policy->study_statistics) ++policy->study_statistics->deferred_validations;
                    return false;
                }
                StageTimer timer(candidate_statistics.validation_ns);
                ++candidate_statistics.validation_calls;
                {
                    TraceScope validation_trace(policy,"original_input_validation",static_cast<uint64_t>(precision),
                        static_cast<uint64_t>(reduced_degree),&candidate_statistics,&statistics);
                    isolated = _acb_poly_validate_roots(roots.data(), reduced.value()->coeffs, reduced_degree + 1, precision);
                }
#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
                if (state.portfolio) state.portfolio->certified(static_cast<uint64_t>(isolated));
#endif
                if (!portfolio_engine(engine) && !state.repeat && !progress.converged && !stage_stopped &&
                    reduced_degree >= heuristics::interaction_tree_degree &&
                    state.charged < allowance) {
                    TraceScope probe(policy, "precision_probe", static_cast<uint64_t>(precision),
                                     static_cast<uint64_t>(reduced_degree - isolated), &candidate_statistics);
                    state.repeat =
                        candidate_arithmetic_adequate(reduced.value(), roots.view(), isolated, precision);
                }
                state.unresolved_count = reduced_degree - isolated;
                acb_zero(state.unresolved_center.value());
                for (slong i = isolated; i < reduced_degree; ++i)
                    acb_add(state.unresolved_center.value(), state.unresolved_center.value(), roots.data() + i,
                            precision);
                if (state.unresolved_count)
                    acb_div_ui(state.unresolved_center.value(), state.unresolved_center.value(),
                               static_cast<ulong>(state.unresolved_count), precision);
                acb_get_mid(state.unresolved_center.value(), state.unresolved_center.value());
            }
        }
        bool reusable = true;
        auto next = std::make_unique<RootVector>(reduced_degree);
        for (slong i = 0; i < reduced_degree; ++i) {
            acb_get_mid(next->data() + i, roots.data() + i);
            if (!acb_is_finite(next->data() + i))
                reusable = false;
            for (slong j = 0; j < i; ++j)
                if (acb_equal(next->data() + i, next->data() + j))
                    reusable = false;
        }
        if (legacy_centers_engine(engine))
            state.centers = reusable ? std::move(next) : nullptr;
        else if (reusable) {
            state.centers = std::move(next);
            state.inactive.assign(static_cast<size_t>(reduced_degree), 0);
            for (slong i = 0; i < isolated; ++i)
                state.inactive[static_cast<size_t>(i)] =
                    acb_rel_accuracy_bits(roots.data() + i) >= resources::precision_sum(options.bits, 16);
        }
        // The entire vector is accepted only after inclusion checks certify ALL
        // simple roots. Unvalidated centers are never promoted to certificates.
        if (isolated != reduced_degree) {
            CountContext context(statistics, precision < static_cast<slong>(options.max_work_bits));
            context.real_shift = policy && policy->pellet_real_shift;
            bool needs_precision = false;
            const bool recovered = allow_clusters &&
                recover_clusters(result, reduced.value(), roots.data(), reduced_degree, multiplicity, options,
                                 precision, &context, &needs_precision);
            if (needs_precision) {
                state.repeat = false;
                if (policy && policy->trace) {
                    TraceEvent event{"cluster_recovery", "precision", static_cast<uint64_t>(precision),
                                     static_cast<uint64_t>(reduced_degree - isolated), &candidate_statistics, &statistics};
                    event.reason = "counting-arithmetic-limited";
                    policy->trace(policy->trace_data, event);
                }
            }
            return recovered;
        }
    }
    const bool is_real = acb_poly_is_real(reduced.value());
    const bool real_validated =
        !binomial && is_real && acb_poly_validate_real_roots(roots.data(), reduced.value(), precision);
    for (slong i = 0; i < reduced_degree; ++i) {
        Root root;
        acb_set(root.box.value(), roots.data() + i);
        root.count = root.multiplicity = multiplicity;
        root.certificate = reduced_degree == 1 ? "linear-ball-division"
                           : binomial          ? "binomial-nth-root-and-unity-cover"
                                               : "native-candidates-with-acb-inclusion";
        classify(root, is_real, real_validated);
        root.publish(precision);
        result.push_back(std::move(root));
    }
    return true;
}

} // namespace certpolysolve::detail
