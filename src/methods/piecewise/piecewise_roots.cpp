// The piecewise method as a resumable state (PiecewiseRootState). Its certificates are Rouche's test on a chart's
// model (certify_root_model) with the chart's root count (count_root_model), in chart_certificate.cpp; the
// certificate attempt the recipes share (machine_newton_cover) is in certificate_attempt.cpp.
//
// PiecewiseRootState::Impl::advance spends one work unit a step. A round (start_round) streams the nonzero part's
// models at model_bits (PiecewiseModelCursor); each model (prepare_model) takes its count, the carried points in it
// and, where discovery is owed, fresh proposals, which are polished and certified (process_pending, process_candidate:
// a coarse guard, kept in PiecewiseOwnership, then the final certificate; for real input the conjugate too). At a batch
// boundary the source finisher may refine owned roots by Newton steps on the input itself (begin_finishing,
// advance_finisher). At a round's end deflated discovery may seed missing roots (seed_missing_roots: Aberth's
// correction against owned roots, from Fujiwara's bound); then a discovery recovery, another round, or the next rung:
// doubled precision where rounding refusals outnumber model refusals, else doubled model_bits. regime_exit_due leaves
// when counted work passes a global solve's predicted cost; a stop short of the goal returns a complete cover if it can
// (cover_at_stop). Axis mode (config.axis_only) builds only the charts centred on the real axis, whose real models
// isolate the real roots by Bernstein bisection (axis_classify, axis_settle, axis_solve).
//
// State: PiecewiseRootState::Impl, across advance calls: input, config, stats, result (the certified regions, with
// index and guards), model_bits and precision, carried (the candidate pool) and ownership, models and pending,
// finishing and the axis round. Stable proof identities live in piecewise_ownership.cpp. The free functions keep none.
// The class is in piecewise_state.hpp (W8's S4). Axis mode has its own PiecewiseAxisPass; the source finisher
// is in piecewise_finisher.cpp, the rounds in piecewise_rounds.cpp, the rest (the constructor, advance,
// configure, snapshot, the guards and the observers) here.
// Entry points: PiecewiseRootState (prepared_piecewise.cpp: one per square-free factor, or one for the whole input).
// solve_piecewise_roots, certify_piecewise_points,
// polish_piecewise_candidates and propose_piecewise_candidates only by tests and benchmark tools.
// Not here: the chart certificate (chart_certificate.cpp), the certificate attempt (certificate_attempt.cpp), the
// models and root_model_error (piecewise_models.cpp), proposals and polish (piecewise_proposals.cpp), square-free
// factors (prepared_piecewise.cpp, exact_preprocessing.cpp), the handoff to secular (session.cpp).
// Design notes: docs/piecewise-root-certification.md (models, the rectangle certificate, ownership, the finisher),
// docs/piecewise-ownership-invariants-20260920.md, docs/c1-axis-charts-note-20261015.md (axis mode),
// docs/production-methods.md (the recipes), docs/same-routes-every-platform-note-20261021.md (the budgets' arithmetic).
#include "methods/piecewise/piecewise_roots.hpp"
#include "methods/piecewise/piecewise_state.hpp"
#include "methods/piecewise/finishing/piecewise_finisher_source.hpp"
#include "methods/piecewise/piecewise_proposals.hpp"
#include <cstdio>
#include <cstdlib>
#include <chrono>
#include "core/dovetail.hpp"
#include "methods/piecewise/piecewise_ownership.hpp"
#include "methods/piecewise/root_index.hpp"
#include "methods/piecewise/piecewise_points.hpp"
#include "core/telemetry.hpp"
#include <algorithm>
#include <cmath>
#include <optional>
#include <limits>
#include <complex>
#include <numbers>
#include <tuple>
#include <set>
#include <span>

namespace certpolysolve::detail {

void polish_piecewise_candidates(RootVector &points, const RootModel &model, const cps_options &options,
                                 slong precision, PiecewiseRootStatistics &stats, PiecewiseRootConfig config) {
    std::vector<bool> stable(static_cast<size_t>(points.size()), false);
    polish(points, points.size(), model, options, precision, stats, &stable, config);
}

slong propose_piecewise_candidates(RootVector &points, const RootModel &model, const cps_options &options,
                                  slong precision, PiecewiseRootStatistics &stats, PiecewiseRootConfig config) {
    const slong degree = acb_poly_degree(model.polynomial.value());
    if (degree < 1 || points.size() != degree)
        throw std::invalid_argument("proposal replay requires one slot per local polynomial degree");
    std::vector<bool> stable(static_cast<size_t>(degree), false);
    certpolysolve::detail::CandidateStatistics machine_statistics{};
    const slong count = proposals(points, model, options, machine_statistics, precision, stats, config);
    polish(points, count, model, options, precision, stats, &stable, config);
    return count;
}

PiecewiseRootState::Impl::~Impl() = default;

PiecewiseRootState::Impl::Impl(const DensePolynomial &p, PiecewiseRootConfig c, const std::vector<ComplexBall> *points)
    : input(p), config(c), supplied(points != nullptr) {
    length = input.coefficients.size();
    while (length && input.coefficients[length - 1].zero()) --length;
    result.degree = static_cast<int64_t>(length) - 1;
    while (zeros + 1 < length && input.coefficients[zeros].zero()) ++zeros;
    degree = length ? static_cast<slong>(length - zeros - 1) : 0;
    if (config.durable_ownership && !supplied)
        ownership = std::make_unique<PiecewiseOwnership>(static_cast<uint64_t>(degree));
    reflect = config.conjugate_symmetry && std::all_of(input.coefficients.begin(),
        input.coefficients.end(), [](const Coefficient &c) { return c.is_real(); });
    if (points) carried = *points;
    if (zeros) {
        Root root;
        root.count = root.multiplicity = zeros;
        root.realness = CPS_REAL;
        root.certificate = "exact zero monomial factor";
        root.publish(64);
        index.insert(root, 0);
        result.roots.push_back(std::move(root));
        covered = zeros;
    }
}
bool PiecewiseRootState::Impl::valid(const cps_options &o) const {
    return !o.real_only && (engine_goal(o) == CPS_ROOT_ACCURACY || engine_goal(o) == CPS_ROOT_APPROXIMATION) &&
        std::all_of(input.coefficients.begin(), input.coefficients.begin() + static_cast<ptrdiff_t>(length),
                    [](const Coefficient &c) { return c.exact(); });
}
slong PiecewiseRootState::Impl::clamped_rung(slong next, const cps_options &o) const {
    const auto cap = static_cast<slong>(std::min<uint64_t>(o.max_work_bits, static_cast<uint64_t>(std::numeric_limits<slong>::max())));
    return next > cap && precision < cap ? cap : next;
}
std::vector<Root> PiecewiseRootState::Impl::cover_at_stop(const cps_options &o) const {
    std::vector<Root> cover = result.roots;
    if (!length || finished()) return cover;
    ComplexPolynomial polynomial;
    ComplexBall coefficient;
    for (size_t j = 0; j < length; ++j) {
        input.coefficients[j].enclose(coefficient.value(), precision);
        acb_poly_set_coeff_acb(polynomial.value(), static_cast<slong>(j), coefficient.value());
    }
    complete_cover_at_stop(cover, polynomial.value(), precision, &o);
    return cover;
}
void PiecewiseRootState::Impl::remember_ownership(const Root &proof, size_t candidate, bool mirror) {
    if (!ownership) return;
    OptionalStageTimer timer(stats.ownership_ns, config.detailed_accounting && mirror);
    const auto admitted = ownership->admit(proof);
    if (!admitted) { ++stats.ownership_conflicts; return; }
    const size_t id = *admitted;
    if (config.source_finisher) {
        accurate_owner.resize(std::max(accurate_owner.size(), id + 1), false);
        if (ownership->proof(id).count == 1 && !accurate_owner[id]) unfinished_owners.insert(id);
    }
    owner_candidate.resize(std::max(owner_candidate.size(), id + 1), new_owner);
    candidate_owner.resize(carried.size(), PiecewiseOwnership::none);
    if (owner_candidate[id] == new_owner) {
        if (candidate == new_owner || candidate >= carried.size() ||
            !ownership->contains(id, carried[candidate].value()) || owner_of(candidate) != PiecewiseOwnership::none) {
            candidate = carried.size();
            ComplexBall seed;
            acb_get_mid(seed.value(), ownership->proof(id).box.value());
            carried.push_back(std::move(seed));
            candidate_owner.push_back(PiecewiseOwnership::none);
        }
        candidate_owner[candidate] = id;
        owner_candidate[id] = candidate;
    }
    if (mirror && reflect && strictly_nonreal(proof))
        remember_ownership(conjugate_proof(proof), new_owner, false);
}
void PiecewiseRootState::Impl::snapshot(bool terminal, const cps_options &o, bool force) {
    if (!config.observe_snapshot) return;
    const bool first_discovery = config.snapshot_discovery && !discovery_snapshot_emitted && discovered();
    if (!terminal && !force && !first_discovery && (!config.snapshot_due || !config.snapshot_due())) return;
    if (first_discovery) discovery_snapshot_emitted = true;
    StageTimer timer(stats.diagnostic_ns);
    if (!diagnostic_source) diagnostic_source = std::make_shared<const DensePolynomial>(input);
    PiecewiseRootSnapshot saved;
    saved.source = diagnostic_source;
    saved.source_id = config.diagnostic_source_id;
    saved.multiplicity = config.diagnostic_multiplicity;
    saved.bits = o.bits;
    saved.terminal = terminal;
    saved.statistics = stats;
    saved.statistics.discovery_count = ownership ? zeros + ownership->count() : 0;
    saved.statistics.final_accurate_count = covered;
    if (ownership) {
        saved.singleton_ownership = ownership->singletons();
        // Reconstruct diagnostic membership from published proofs. Do not
        // enable the finisher or modify its scheduling ledger to observe it.
        std::vector<bool> accurate(owner_candidate.size(), false);
        for (const auto &root : result.roots) {
            const auto id = ownership->identifying(root);
            if (id && *id < accurate.size()) accurate[*id] = true;
        }
        for (size_t id = 0; id < owner_candidate.size(); ++id) {
            if (!ownership->active(id) || accurate[id]) continue;
            PiecewiseUnfinishedOwner item;
            item.id = id;
            item.proof = ownership->proof(id);
            item.candidate_id = owner_candidate[id];
            if (item.candidate_id < carried.size()) item.candidate = carried[item.candidate_id];
            if (item.candidate_id < candidate_history.size()) {
                const auto &[attempts, event] = candidate_history[item.candidate_id];
                item.attempts = attempts;
                if (attempts) item.last_attempt = event;
            }
            saved.unfinished.push_back(std::move(item));
        }
    }
    config.observe_snapshot(std::move(saved));
}
void PiecewiseRootState::Impl::configure(const cps_options &o) {
    // C1: no global finishing or deflated discovery (both work from the global count), and no regime exit: it
    // weighs piecewise's work against a global method's for all roots, and would stop an axis solve whose
    // passes need only more model bits. Leaving the axis mode is Auto's routing decision.
    // History: docs/code-history-20261021.md, piecewise_roots.cpp 12.
    if (config.axis_only) {
        config.source_finisher = false;
        config.deflated_discovery = false;
        config.regime_exit = false;
    }
    const bool changed = configured && (o.bits != settings.bits || o.goal != settings.goal ||
        !fmpq_equal(o.absolute_tolerance.significand.value(), settings.absolute_tolerance.significand.value()) ||
        !fmpz_equal(o.absolute_tolerance.exponent.value(), settings.absolute_tolerance.exponent.value()) ||
        o.absolute_tolerance.radix != settings.absolute_tolerance.radix);
    // Discovery stopping depends on the goal. (The solve's iteration
    // budget, cps_options::max_iterations, does not count this method's
    // local refinement, so it is not part of the test.)
    if (changed) {
        searched_model_bits = searched_precision = 0;
        progress_model_bits = progress_precision = 0;
        finishing.reset();
        pilot_work = marginal_pilot_work = 0;
        failed_finisher_precision = 0;
    }
    settings = o;
    settings.cancel = nullptr;
    settings.cancel_data = nullptr;
    if (!configured) precision = piecewise_start_precision(o.bits, model_bits);
    configured = true;
    if (!changed) return;
    // An axis-mode state answers one query: real_only_route.cpp builds a fresh one for each, and only the global
    // mode's states are advanced again with a new goal (PreparedPiecewiseState's rounds), for which the filter
    // below is written.
    if (config.axis_only)
        throw std::logic_error("an axis-mode piecewise state answers one query");
    pending.reset();
    models.reset();
    index = RootIndex{};
    guards = GuardIndex{};
    seed_guards = GuardIndex{};
    std::vector<Root> accurate;
    covered = 0;
    for (auto &root : result.roots) {
        if (meets_goal(root, o, std::max<slong>(precision, resources::precision_sum(o.bits, 32)))) {
            covered += root.count;
            index.insert(root, accurate.size());
            guards.insert(root, root);
            accurate.push_back(std::move(root));
        } else {
            // A cluster's center is a seed, never a claim that all of its
            // members coincide. Local counts/reseeding can split it later.
            if (!supplied) carried.push_back(root.box);
            retained.push_back(std::move(root));
        }
    }
    result.roots = std::move(accurate);
    if (config.source_finisher && ownership) {
        accurate_owner.assign(owner_candidate.size(), false);
        unfinished_owners.clear();
        for (size_t id = 0; id < owner_candidate.size(); ++id)
            if (ownership->active(id) && ownership->proof(id).count == 1) unfinished_owners.insert(id);
        for (const auto &root : result.roots) mark_finished_owner(root);
    }
    force_reseed = !discovered();
    precision = std::max(precision, resources::precision_sum(resources::precision_sum(resources::precision_product(model_bits, 2), o.bits), 64));
}
void PiecewiseRootState::Impl::observe_population(bool terminal) {
    if (!config.observe_population) return;
    StageTimer timer(stats.diagnostic_ns);
    PiecewisePopulationSnapshot saved;
    saved.source = config.diagnostic_source_id;
    saved.round = stats.rounds;
    saved.model_bits = static_cast<uint64_t>(model_bits);
    saved.arithmetic_bits = static_cast<uint64_t>(precision);
    saved.discovered = ownership ? zeros + ownership->count() : 0;
    saved.accurate = covered;
    saved.terminal = terminal;
    saved.points.reserve(carried.size());
    for (size_t i = 0; i < carried.size(); ++i) {
        PiecewisePopulationPoint point;
        point.point = carried[i];
        if (i < round_points.size()) point.indexed_point = round_points[i];
        point.owner = owner_of(i);
        point.inside_final_guard = guards.contains(carried[i].value());
        if (ownership) if (auto owner = ownership->containing(carried[i].value()))
            point.containing_owner = *owner;
        saved.points.push_back(std::move(point));
    }
    config.observe_population(std::move(saved));
}
void PiecewiseRootState::Impl::observe_selection(const RootModel &model, size_t i) {
    if (!config.observe_selection) return;
    StageTimer timer(stats.diagnostic_ns);
    PiecewiseSelectionEvent event;
    event.source = config.diagnostic_source_id;
    event.candidate = i;
    event.model = stats.models;
    event.model_bits = static_cast<uint64_t>(model_bits);
    event.arithmetic_bits = static_cast<uint64_t>(precision);
    event.moved = !acb_equal(round_points[i].value(), carried[i].value());
    ComplexBall coordinate;
    RealBall norm, one;
    arb_one(one.value());
    auto inspect = [&](const ComplexBall &point, bool &inside, double &approximate) {
        acb_sub(coordinate.value(), point.value(), model.center.value(), precision);
        acb_div_arb(coordinate.value(), coordinate.value(), model.radius.value(), precision);
        acb_abs(norm.value(), coordinate.value(), precision);
        inside = arb_lt(norm.value(), one.value());
        approximate = arf_get_d(arb_midref(norm.value()), ARF_RND_NEAR);
    };
    inspect(round_points[i], event.indexed_in_disk, event.indexed_norm);
    inspect(carried[i], event.current_in_disk, event.current_norm);
    arb_log_base_ui(norm.value(), model.radius.value(), 2, 64);
    event.radius_log2 = arf_get_d(arb_midref(norm.value()), ARF_RND_NEAR);
    config.observe_selection(event);
}
bool PiecewiseRootState::Impl::held_copy(const RootModel &model, const acb_t local) {
    ComplexBall point;
    acb_mul_arb(point.value(), local, model.radius.value(), precision);
    acb_add(point.value(), point.value(), model.center.value(), precision);
    acb_get_mid(point.value(), point.value());
    if (!acb_is_finite(point.value())) return false;
    if (guards.contains(point.value())) {
        ++stats.retired_points;
        ++stats.unpolished_copies;
        return true;
    }
    if (ownership)
        if (const auto other = ownership->containing(point.value()); other && ownership->proof(*other).count == 1) {
            ++stats.owned_proposals_retired;
            ++stats.unpolished_copies;
            return true;
        }
    return false;
}
bool PiecewiseRootState::Impl::certify_guard(Root &root, const RootModel &model, const acb_t point, const cps_options &o,
    ChartFailure *failure) {
    const auto query = guard_options(model, o);
    ChartFailure reason=ChartFailure::Candidate;
    const bool proved = certify_root_model(root, model, point, query, precision, &reason,
                                           nullptr, nullptr, nullptr, 0, config.lazy_guard_jet, &stats.work.certificate,
                                           &candidate_jet);
    if (failure) *failure=reason;
    if (config.observe_certificate && !o.cancelled()) {
        StageTimer timer(stats.diagnostic_ns);
        ComplexBall copy;acb_set(copy.value(),point);
        config.observe_certificate(model,copy,query,precision,proved,proved?-1:static_cast<int>(reason),false);
    }
    if (config.observe_guard) {
        ComplexBall copy;
        acb_set(copy.value(), point);
        config.observe_guard(model, copy, query, precision, proved);
    }
    return proved;
}
void PiecewiseRootState::Impl::publish_reflection() {
    if (!reflect || !strictly_nonreal(result.roots.back())) return;
    auto reflected = conjugate_proof(result.roots.back());
    if (index.overlaps(reflected, result.roots)) return;
    if (reflected.count > static_cast<uint64_t>(result.degree) - covered)
        throw std::logic_error("piecewise conjugates exceed input degree");
    covered += reflected.count;
    stats.reflected_roots += reflected.count;
    index.insert(reflected, result.roots.size());
    guards.insert(reflected, reflected);
    result.roots.push_back(std::move(reflected));
    mark_finished_owner(result.roots.back());
}
bool PiecewiseRootState::Impl::regime_exit_due() const {
    if (!config.regime_exit || config.regime_exit_units == 0 || finished())
        return false;
    // The counted work against the counts at their measured prices: ratios, in binary64.
    const double units = static_cast<double>(config.regime_exit_units);
    if (config.regime_exit_estimate > 0)
        return counted_work(stats.work) >
               binary64_operation_work *
                   (estimate_operation_overhead * static_cast<double>(config.regime_exit_estimate) + units);
    if (config.caller_method)
        return counted_work(stats.work) > binary64_operation_work * units;
    if (!stats.construction_units)
        return false;
    return counted_work(stats.work) * static_cast<double>(stats.construction_units) >
           construction_work(stats.work) * units;
}
cps_result PiecewiseRootState::Impl::advance(const cps_options &o) {
    // Reject before configure can filter proofs or discard continuation
    // state. An unsupported query must not inherit prior success flags.
    if (length && !valid(o)) {
        cps_result unsupported;
        unsupported.degree = result.degree;
        unsupported.goal = engine_goal(o);
        unsupported.termination = CPS_UNSUPPORTED_QUERY;
        unsupported.diagnostics = "piecewise solving requires an exact input and a supported complex root goal";
        return unsupported;
    }
    configure(o);
    result.goal = engine_goal(o);
    result.progress.work_units_this_call = 0;
    if (!length) { result.termination = CPS_ZERO_POLYNOMIAL; return result; }
    while (!finished() && !refused && !stats.regime_exits && !o.cancelled() && !axis_spent() &&
           result.progress.work_units_this_call < o.work_units &&
           static_cast<uint64_t>(precision) <= o.max_work_bits) {
        // Before every unit, so a solve in the wrong regime leaves within
        // a round, not only after one.
        if (regime_exit_due()) {
            ++stats.regime_exits;
            break;
        }
        snapshot(false, o);
        ++result.progress.work_units_this_call;
        ++result.statistics_v2.work_units;
        if (finishing) { advance_finisher(o); continue; }
        if (!models) { start_round(o); continue; }
        if (pending) { process_pending(models->model(), o); continue; }
        // An axis-mode round that has failed builds nothing more: it ends now (axis_settle).
        const auto next = config.axis_only && axis->failed() ? PiecewiseModelStep::Done : models->advance(o);
        if (next == PiecewiseModelStep::BatchBoundary) {
            begin_finishing(models->next_batch_work(), o);
        } else if (next == PiecewiseModelStep::Model) {
            prepare_model(models->model(), o);
            if (pending) process_pending(models->model(), o);
        } else if (next == PiecewiseModelStep::Refused) {
            refused = true;
        } else if (next == PiecewiseModelStep::Done) {
            if (finish_round(o) == RoundStep::Stop) break;
        }
    }
    // A cancelled candidate-generation pass may have produced only part
    // of its proposals. Do not record it as an exhausted discovery pass.
    if (force_reseed && o.cancelled()) recovery_interrupted = true;
    result.complete = finished();
    stats.discovery_count = ownership ? zeros + ownership->count() : 0;
    stats.final_accurate_count = covered;
    if (config.observe) config.observe(stats, "advance-complete");
    observe_population(true);
    snapshot(true, o);
    result.goal_met = result.complete;
    // What the caller receives at a stop short of the goal is a complete
    // cover (cover_at_stop); the state keeps its certified regions.
    std::vector<Root> returned = result.roots;
    bool returned_complete = result.complete;
    // An axis-mode stop is an attempt's end: its caller answers with another route, and a complex cover means
    // nothing to a real-only answer and can cost most of a short attempt.
    // History: docs/code-history-20261021.md, piecewise_roots.cpp 15.
    if (!result.complete && !o.cancelled() && !config.axis_only) {
        StageTimer timer(stats.certification_ns);
        returned = cover_at_stop(o);
        uint64_t counted = 0;
        for (const auto &root : returned) counted += root.count;
        returned_complete = length && counted == length - 1;
        if (!returned_complete) returned = result.roots;
    }
    // A range refusal is the input's (its exponents exceed the envelope's
    // working range), not a budget's.
    const bool range_refused = refused && stats.range_refusals;
    // In axis mode only the solve's spent work budget is a work limit; a refusal or the attempt's last round
    // is a stop of the attempt itself, whose caller takes another route.
    const bool budget_spent = result.progress.work_units_this_call >= o.work_units;
    result.termination = result.complete ? CPS_SUCCESS : o.cancelled() ? CPS_CANCELLED :
        range_refused ? CPS_UNSUPPORTED_QUERY :
        config.axis_only ? (!refused && !axis_spent() && budget_spent ? CPS_WORK_LIMIT : CPS_PRECISION_LIMIT) :
        refused || budget_spent ? CPS_WORK_LIMIT : CPS_PRECISION_LIMIT;
    result.progress.stop_reason = result.complete ? CPS_STOP_NONE : o.cancelled() ? CPS_STOP_CANCELLED :
        range_refused ? CPS_STOP_INPUT :
        result.termination == CPS_WORK_LIMIT ? CPS_STOP_WORK_BUDGET : CPS_STOP_PRECISION;
    result.progress.pending_components = result.complete ? 0 : 1;
    result.components = {{0, result.work_bits, 0, covered,
        result.complete ? CPS_COMPONENT_FINISHED : result.termination == CPS_PRECISION_LIMIT ?
        CPS_COMPONENT_PRECISION_LIMIT : CPS_COMPONENT_ACTIVE}};
    result.diagnostics = result.complete ? "complete persistent piecewise cover" :
        returned_complete ? "piecewise stop with a complete cover; candidates and certificates retained" :
        "partial persistent piecewise cover; candidates and certificates retained";
    if (result.complete) return result;
    cps_result out = result;
    out.roots = std::move(returned);
    out.complete = returned_complete;
    return out;
}
PiecewiseRootState::PiecewiseRootState(const DensePolynomial &p, PiecewiseRootConfig c,
    const std::vector<ComplexBall> *points) : impl_(std::make_unique<Impl>(p, c, points)) {}
PiecewiseRootState::~PiecewiseRootState() = default;
cps_result PiecewiseRootState::advance(const cps_options &o) { return impl_->advance(o); }
const PiecewiseRootStatistics &PiecewiseRootState::statistics() const { return impl_->stats; }
size_t PiecewiseRootState::certified_regions() const { return impl_->result.roots.size(); }
cps_result solve_piecewise_roots(const DensePolynomial &p, const cps_options &o, PiecewiseRootStatistics &s,
                                PiecewiseRootConfig config) {
    PiecewiseRootState state(p, config);
    auto result = state.advance(o);
    s = state.statistics();
    return result;
}
cps_result certify_piecewise_points(const DensePolynomial &p, const cps_options &o, PiecewiseRootStatistics &s,
                                   const std::vector<ComplexBall> &points, PiecewiseRootConfig config) {
    PiecewiseRootState state(p, config, &points);
    auto result = state.advance(o);
    s = state.statistics();
    return result;
}
} // namespace certpolysolve::detail
