// The query-level global driver: dispatch input routes, prove the degree,
// schedule factors, then publish the final stop, progress and diagnostics.
// Factor algorithms and cover assembly live in solver.cpp. Every phase uses
// the same latched cancellation callback, whose lifetime spans advance().
#include "solve/solve_state_impl.hpp"
#include "solve/factor_pass.hpp"
#include "solve/factor_state.hpp"
#include "core/telemetry.hpp"
#include "solve/program_route_state.hpp"
#include "certify/factor_cover.hpp"
#include "certify/separation.hpp"
#include "numeric/exact/secular_polynomial.hpp"
#include <algorithm>
#include <locale>
#include <optional>
#include <sstream>
#include <utility>
#include <vector>

namespace certpolysolve::detail {

bool SolveState::Impl::advance_input_route(const cps_options &o) {
    if (program_route) {
        program_route->solve(o, policy, result, initial_precision);
        return true;
    }
    if (secular_input && policy.backend && !secular_converted && !native_fallback && !o.real_only)
        if (native_advance(o))
            return true;
    if (secular_input && !secular_converted) {
        if (result.progress.work_units_this_call >= o.work_units || o.cancelled()) {
            result.termination = o.cancelled() ? CPS_CANCELLED : CPS_WORK_LIMIT;
            result.progress.stop_reason = o.cancelled() ? CPS_STOP_CANCELLED : CPS_STOP_WORK_BUDGET;
            return true;
        }
        ++result.statistics_v2.work_units;
        ++result.progress.work_units_this_call;
        StageTimer timer(result.statistics.preprocessing_ns);
        if (!secular_polynomial(input, *secular_input, o, result.diagnostics)) {
            result.termination = o.cancelled() ? CPS_CANCELLED : CPS_UNSUPPORTED_QUERY;
            result.progress.preprocessing_limited = !o.cancelled();
            result.progress.stop_reason = o.cancelled() ? CPS_STOP_CANCELLED : CPS_STOP_INPUT;
            return true;
        }
        secular_converted = true;
    }
    return false;
}

void SolveState::Impl::prepare_coefficient_query(const cps_options &o) {
    if (result.complete)
        result.goal_met = goal(result.roots, std::max<uint64_t>(64, result.progress.max_outer_bits));
    if (result.degree < 0) {
        result.termination = CPS_ZERO_POLYNOMIAL;
        result.progress.stop_reason = CPS_STOP_INPUT; // the input's verdict (I4; the program route has it too)
        terminal = true;
    }
    if (o.real_only && o.goal == CPS_CLUSTER_RADIUS) {
        result.termination = CPS_UNSUPPORTED_QUERY;
        result.progress.stop_reason = CPS_STOP_INPUT;
        terminal = true;
    }
    // A real-only query on input that is not exact and real can never be
    // answered, whatever the budget: say so before the first work unit.
    // Left to preprocessing, which runs on the second unit, a one-unit
    // budget would report a work limit, as if more work could help.
    // The stop code is the input's on every route, as the explicit
    // methods' rejection already has it. History:
    // docs/code-history-20261021.md, solver.cpp 24.
    if (o.real_only && !terminal && !preprocessed &&
        !std::all_of(input.coefficients.begin(), input.coefficients.end(),
                     [](const Coefficient &c) { return c.exact() && c.is_real(); })) {
        result.termination = CPS_UNSUPPORTED_QUERY;
        result.progress.stop_reason = CPS_STOP_INPUT;
        preprocessing_status = "real_only requires exact real input";
        terminal = true;
    }
    if (policy.prepared_subdivision && !terminal &&
        (o.real_only || !source_structure || !source_structure->exact)) {
        result.termination = CPS_UNSUPPORTED_QUERY;
        result.progress.stop_reason = CPS_STOP_INPUT;
        preprocessing_status = "standalone subdivision requires exact coefficient input and a complex query";
        terminal = true;
    }
}

SolveState::Impl::DriverStep SolveState::Impl::prove_input_degree(const cps_options &o) {
    ++result.statistics_v2.work_units;
    ++result.progress.work_units_this_call;
    ComplexPolynomial p;
    enclose_input(p.value(), static_cast<slong>(initial_precision));
    const auto leading = p.value()->coeffs + result.degree;
    if (acb_is_finite(leading) && !acb_contains_zero(leading)) {
        degree_proved = true;
        if (result.degree == 0) {
            result.complete = result.goal_met = true;
        } else if (!o.real_only) {
            Root global;
            if (cauchy_cluster(global, p.value(), static_cast<slong>(initial_precision))) {
                if (!result.complete)
                    result.roots = {std::move(global)};
                result.complete = true;
                result.goal_met = goal(result.roots, initial_precision);
            }
        }
    } else if (initial_precision == o.max_work_bits) {
        result.termination = CPS_INDETERMINATE_DEGREE;
        result.progress.stop_reason = CPS_STOP_INPUT;
        // A solve under a higher cap may prove it.
        return DriverStep::Stop;
    } else {
        initial_precision = std::min(o.max_work_bits, initial_precision * 2);
        ++result.statistics_v2.precision_increases;
    }
    return DriverStep::Continue;
}

SolveState::Impl::DriverStep SolveState::Impl::finish_factor_pass(const cps_options &o, FactorPass &pass) {
    // The pass counted finished factors and factors that a
    // budget blocked alike. If a budget blocked one (the
    // iteration budget, the cells, the precision, the
    // input), nothing the overlap branch below starts can
    // finish it: raising another factor's goal, or its own,
    // gives it no cells, precision or iterations. The call
    // ends with that budget's reason (FactorPass::stop_reason), as one
    // full pass of finished or blocked factors does: the
    // branch would otherwise raise a separation goal again
    // on every pass, at no work, up to the ceiling, and the
    // solve would stop "at the precision cap". History:
    // docs/code-history-20261021.md, solver.cpp 25.
    if (pass.has_blockage()) {
        result.progress.stop_reason = pass.stop_reason();
        return DriverStep::Stop;
    }
    // Individually accurate factors can still have overlapping
    // exports. Tighten only the conflicting factor covers.
    std::vector<const std::vector<Root> *> covers;
    for (const auto &g : factors)
        covers.push_back(&g->output_cover());
    // The extracted point is part of the final cover too. In the standalone
    // route it must keep a nearby factor active until their regions separate.
    // Its zero radius makes the overlapping nonzero region the wider one;
    // two exact coincident points would contradict coprime preparation.
    std::vector<Root> zero_cover;
    if ((policy.prepared_subdivision || heuristic_budget) && zero_order) {
        zero_cover.push_back(exact_zero_point(zero_order));
        covers.push_back(&zero_cover);
    }
    const auto chosen = factor_to_separate(covers);
    if (chosen && *chosen >= factors.size())
        throw std::logic_error("a prepared nonzero factor contains the exact zero root");
    FactorState *affected = chosen ? factors[*chosen].get() : nullptr;
    if (!affected) {
        // The pass also counts factors blocked by a budget
        // (the stop branch below). If one of those has not
        // met its goal, its budget is why the solve stops.
        if (!std::all_of(factors.begin(), factors.end(),
                         [&](const auto &g) { return g->ready && factor_goal(*g); }))
            result.progress.stop_reason = pass.stop_reason();
        return DriverStep::Stop;
    }
    // The tighter goal for the affected factor: its accuracy (the query's bits plus its separation
    // bits) doubles, by at least 16 bits. A factor whose working precision is already at the cap is
    // tightened all the same: the precision bounds the accuracy it can certify, not the accuracy
    // its cover has, and its own step reports the cap if it cannot meet the tighter goal there.
    // A goal set from the working precision instead would sit above that precision every round,
    // which would then have to grow past it: the goal would chase the precision to the cap.
    // History: docs/code-history-20261021.md, solver.cpp 26.
    if (affected->separation_bits >= resources::precision_max - settings.bits) {
        result.progress.stop_reason = CPS_STOP_PRECISION;
        return DriverStep::Stop;
    }
    affected->separation_bits = next_separation(settings.bits, affected->separation_bits);
    // A power composite's finished test then asks its reduced cover for the tighter goal
    // (factor_goal), so the factor is stepped rather than counted as finished at once. Where the
    // reduced cover meets it already (an exact root: w - 2 for z^2 - 2), no step follows, so
    // the lifted cover is recomputed here, at the precision the tighter goal sets (aggregate lifts
    // at the query's bits plus the separation bits plus 32), and the next pass sees the
    // tighter regions.
    affected->overlap_tightened = true;
    if (affected->power > 1)
        aggregate();
    const auto tighter = factor_options(*affected, o);
    if (affected->dovetail_search) {
        affected->dovetail_search->reopen(tighter);
        affected->dovetail_stop = CPS_STOP_NONE;
    }
    if (affected->search)
        affected->search->reopen(tighter);
    else if (!o.real_only)
        affected->local_refinement = true;
    else
        affected->precision = std::min(
            o.max_work_bits, std::max(affected->precision * 2, tighter.bits + 32));
    pass.restart();
    return DriverStep::Continue;
}

SolveState::Impl::DriverStep SolveState::Impl::advance_factor(const cps_options &o, FactorPass &pass) {
    // A synchronous local solve returns its useful prefix when the shared
    // allowance ends. Only its owning outer query starts complete service.
    if (settings.heuristic_budget && !heuristic_budget && settings.heuristic_budget->exhausted()) {
        result.progress.stop_reason = CPS_STOP_WORK_BUDGET;
        return DriverStep::Stop;
    }
    auto &f = *factors[cursor++ % factors.size()];
    if (f.ready && factor_goal(f)) {
        if (pass.visit_finished(factors.size()))
            return finish_factor_pass(o, pass);
        return DriverStep::Continue;
    }
    if (f.iteration_wait) {
        if (pass.visit_blocked(CPS_STOP_ITERATION_BUDGET, factors.size())) {
            result.progress.stop_reason = pass.stop_reason();
            return DriverStep::Stop;
        }
        return DriverStep::Continue;
    }
    ++result.statistics_v2.work_units;
    ++result.progress.work_units_this_call;
    step_kept_covers = false;
    auto reason = factor_step(f, o);
    if (settings.heuristic_budget && !f.complete_fallback) {
        auto &budget = *settings.heuristic_budget;
        budget.observe_precision(f.precision);
        if (f.search) budget.observe_precision(f.search->max_bits());
        if (f.dovetail_search) budget.observe_precision(f.dovetail_search->max_bits());
        if (heuristic_budget && reason == CPS_STOP_PRECISION &&
            o.max_work_bits > budget.ceiling()) {
            // The owner's numerical and retained terminal paths have no
            // further step at Q. A borrowed child's cap is only its local
            // stop: it cannot expire the whole query's remaining allowance.
            // The complete worker restores the caller's cap on the next step.
            budget.finish_at_precision_limit();
            f.blocked = false;
            reason = CPS_STOP_NONE;
        }
    }
    if (policy.snapshot && f.ready && factor_goal(f)) {
        NumericalSnapshot event;
        event.stage = "completed_factor";
        event.factor = f.id;
        event.source_precision = f.precision;
        for (const auto &root : f.output_cover()) {
            event.points.emplace_back();
            acb_get_mid(event.points.back().value(), root.box.value());
            event.counts.push_back(root.count);
        }
        policy.snapshot(std::move(event));
    }
    if (reason == CPS_STOP_CANCELLED) {
        result.progress.stop_reason = reason;
        return DriverStep::Stop;
    }
    // A bounded step may improve a factor's certified cover and
    // then report a precision/cell limit. Publish that improvement
    // before handling the stop; a resource limit is not rollback.
    // A dovetailed Subdivision step that adopted no cover changed
    // none of aggregate's inputs, and aggregate compares every
    // pair of regions: skip it, or a share spent in many small
    // steps costs far more than its metered work.
    if (!step_kept_covers) {
        StageTimer timer(result.statistics.coverage_validation_ns);
        aggregate();
    }
    if (FactorPass::blocks(reason)) {
        // A step whose cover met the goal (aggregate) ends the solve: its budget or precision report
        // says how the step ended, not why the solve stops, so no stall is recorded and the call
        // returns CPS_SUCCESS with CPS_STOP_NONE. History:
        // docs/code-history-20261021.md, solver.cpp 27.
        if (result.goal_met) return DriverStep::Stop;
        if (reason == CPS_STOP_ITERATION_BUDGET)
            f.iteration_wait = true;
        if (pass.visit_blocked(reason, factors.size())) {
            result.progress.stop_reason = pass.stop_reason();
            return DriverStep::Stop;
        }
        return DriverStep::Continue;
    }
    pass.restart();
    return DriverStep::Continue;
}

void SolveState::Impl::finish_coefficient_query(const cps_options &o) {
    // A stop short of the goal that is not a cancellation (a budget, the precision cap, a stall) completes every
    // factor that has had no step yet with its stop cover (aggregate), as piecewise's factor loop completes an
    // unfinished factor: one stop rule for both (the S8 map's rule 3,
    // question Q3; 23 October 2026). Before the preparation there are no factors (the degree proof's Cauchy
    // square stays), a verdict on the input completes nothing, and neither does a cancellation (a time limit is one).
    if (!terminal && !result.goal_met && preprocessed && !settings.real_only && !o.cancelled() &&
        result.progress.stop_reason != CPS_STOP_CANCELLED && result.termination != CPS_INDETERMINATE_DEGREE &&
        std::any_of(factors.begin(), factors.end(), [](const auto &f) { return !f->ready; })) {
        StageTimer timer(result.statistics.coverage_validation_ns);
        aggregate(&o);
    }
    // A met goal is a success whatever the last checkpoint polls, and a verdict on the input (a zero polynomial,
    // an unsupported query, an indeterminate degree) is not turned into a cancellation (docs/api-spec-0.5.md
    // section 6.2, invariants I2 and I4; change C25). History: docs/code-history-20261021.md, solver.cpp 28.
    const bool input_verdict = terminal || result.termination == CPS_INDETERMINATE_DEGREE;
    if (!result.goal_met && !input_verdict && o.cancelled())
        result.progress.stop_reason = CPS_STOP_CANCELLED;
    if (result.goal_met) {
        result.termination = CPS_SUCCESS;
        result.progress.stop_reason = CPS_STOP_NONE;
    } else if (result.progress.stop_reason == CPS_STOP_CANCELLED)
        result.termination = CPS_CANCELLED;
    else if (result.progress.stop_reason == CPS_STOP_WORK_BUDGET ||
             result.progress.stop_reason == CPS_STOP_CELL_BUDGET ||
             result.progress.stop_reason == CPS_STOP_ITERATION_BUDGET)
        result.termination = CPS_WORK_LIMIT;
    else if (!terminal && degree_proved)
        result.termination = CPS_PRECISION_LIMIT;
}

void SolveState::Impl::summarize_factor_progress() {
    result.progress.pending_components = 0;
    result.progress.live_cells = live_cells();
    result.progress.min_pending_bits = 0;
    result.progress.max_pending_bits = 0;
    result.components.clear();
    // Whether each factor is finished, tested once for both loops below (factor_goal is const, and they change no
    // factor).
    std::vector<char> finished(factors.size());
    for (size_t i = 0; i < factors.size(); ++i)
        finished[i] = factors[i]->ready && factor_goal(*factors[i]);
    for (size_t i = 0; i < factors.size(); ++i) {
        const auto &f = *factors[i];
        if (f.search)
            f.search->progress(result.components, i);
        else {
            uint64_t count = 0;
            for (const auto &root : f.output_cover())
                count += root.count;
            result.components.push_back(
                {i, f.precision, 0, count,
                 finished[i] ? CPS_COMPONENT_FINISHED
                             : (f.blocked ? CPS_COMPONENT_PRECISION_LIMIT : CPS_COMPONENT_ACTIVE)});
        }
    }
    for (size_t i = 0; i < factors.size(); ++i)
        if (const auto &f = factors[i]; !finished[i]) {
            result.progress.pending_components += f->search ? f->search->pending() : 1;
            const auto lo = f->search ? f->search->min_bits() : f->precision;
            const auto hi = f->search ? f->search->max_bits() : f->precision;
            if (!result.progress.min_pending_bits || lo < result.progress.min_pending_bits)
                result.progress.min_pending_bits = lo;
            result.progress.max_pending_bits = std::max(result.progress.max_pending_bits, hi);
        }
}

void SolveState::Impl::describe_coefficient_query() {
    std::ostringstream trace;
    trace.imbue(std::locale::classic()); // diagnostics never follow the host's locale
    if (secular_input)
        trace << "exact secular numerator; reconstructed identity verified; ";
    trace << preprocessing_status << "; persistent factor/component search; "
          << result.statistics_v2.work_units << " work units; " << result.progress.pending_components
          << " pending components; stop=" << static_cast<int>(result.progress.stop_reason);
    for (const auto &f : factors)
        if (f->real_line_seeds.computed)
            trace << "; real-line seeds: " << f->real_line_seeds.seeds.size() << " of "
                  << fmpz_poly_degree((f->real_line_seeds.centred ? f->centered : f->integer).value()) << " ("
                  << real_line_seed_stop_name(f->real_line_seeds.statistics.stop) << ")";
    if (result.progress.input_uncertainty)
        trace << "; coefficient uncertainty retained uniformly";
    if (result.progress.preprocessing_limited)
        trace << (policy.prepared_subdivision ? "; exact preparation refused for standalone subdivision"
                                             : "; exact preprocessing limited, unsupported multiplicities remain unknown");
    if (result.termination == CPS_ZERO_POLYNOMIAL)
        trace << "; zero polynomial has no finite all-roots enumeration";
    if (result.progress.stop_reason == CPS_STOP_PRECISION) {
        // The bound of the square-free part, the product of the distinct factors (roots of different
        // factors can be close too), when every factor has integer coefficients; Gaussian-rational and
        // inexact inputs give none here. Computed only at this stop.
        std::optional<SeparationBound> bound;
        if (!factors.empty() &&
            std::all_of(factors.begin(), factors.end(), [](const auto &f) { return f->kind == FactorState::Real; })) {
            IntegerPolynomial product;
            fmpz_poly_one(product.value());
            for (const auto &f : factors)
                fmpz_poly_mul(product.value(), product.value(), f->integer.value());
            bound = separation_bound(product.value());
        }
        trace << "; " << precision_cap_note << settings.max_work_bits << " bits";
        if (bound) {
            trace << "; by Mahler's bound the roots are more than 2^-" << bound->separation_bits
                  << " apart and less than 2^" << bound->size_bits << " in size, so a cap of "
                  << 2 * (bound->separation_bits + bound->size_bits + settings.bits)
                  << " bits suffices to isolate them (cps_options_set_max_work_bits)";
            // The same cap as a field (C14), saturating where the text's arithmetic would wrap.
            result.sufficient_work_bits = resources::multiply(
                resources::add(resources::add(bound->separation_bits, bound->size_bits), settings.bits), 2);
        }
    }
    if (heuristic_budget) {
        trace << "; heuristic-packets=" << heuristic_budget->ticks() << "/" << heuristic_budget->degree()
              << "; heuristic-Q=" << heuristic_budget->ceiling()
              << "; complete-service=" << std::count_if(factors.begin(), factors.end(),
                    [](const auto &f) { return f->complete_fallback; });
    }
    result.diagnostics = trace.str();
}

cps_result SolveState::Impl::advance(const cps_options *proposed) {
    configure(proposed); // Validate before changing any mathematical state.
    cancel_observed = false;
    struct Cancellation {
        Impl *self;
        CancelCallback callback;
        void *data;
        static int poll(void *ptr) {
            auto &c = *static_cast<Cancellation *>(ptr);
            if (!c.self->cancel_observed && c.callback && c.callback(c.data))
                c.self->cancel_observed = true;
            return c.self->cancel_observed ? 1 : 0;
        }
    } cancellation{this, proposed ? proposed->cancel : nullptr,
                   proposed ? proposed->cancel_data : nullptr};
    cps_options o = settings;
    o.cancel = Cancellation::poll;
    o.cancel_data = &cancellation;
    result.progress.work_units_this_call = 0;
    result.progress.stop_reason = CPS_STOP_NONE;
    if (advance_input_route(o))
        return result;
    prepare_coefficient_query(o);
    // One full pass without a productive step ends the call. The pass owns
    // stop precedence; factors retain their individual iteration waits.
    FactorPass pass;
    for (auto &g : factors)
        g->iteration_wait = false;
    while (!terminal && !result.goal_met) {
        if (o.cancelled()) {
            result.progress.stop_reason = CPS_STOP_CANCELLED;
            break;
        }
        if (result.progress.work_units_this_call >= o.work_units) {
            result.progress.stop_reason = CPS_STOP_WORK_BUDGET;
            break;
        }
        if (!degree_proved) {
            if (prove_input_degree(o) == DriverStep::Stop)
                break;
        } else if (!preprocessed) {
            TraceScope trace(&policy, "preprocessing", initial_precision,
                             static_cast<uint64_t>(result.degree), &result.statistics_v3,
                             &result.statistics_v2);
            ++result.statistics_v2.work_units;
            ++result.progress.work_units_this_call;
            preprocess(o);
        } else if (advance_factor(o, pass) == DriverStep::Stop) {
            break;
        }
    }
    finish_coefficient_query(o);
    summarize_factor_progress();
    describe_coefficient_query();
    sort_roots(result.roots);
    return result;
}

#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
ResearchServiceAction SolveState::Impl::ResearchControl::service(void *pointer,
                                                                const ResearchLocalPacketProgress &progress) {
    auto &control = *static_cast<ResearchControl *>(pointer);
    if (!control.advancing || control.servicing)
        throw std::logic_error("local packet service requires one active nonreentrant numerical turn");
    struct Guard {
        bool &flag;
        explicit Guard(bool &value) : flag(value) { flag = true; }
        ~Guard() { flag = false; }
    } guard(control.servicing);
    auto current = progress;
    current.work_units = control.owner->research_inflight_work_units();
    control.engine_peak_bits = std::max(control.engine_peak_bits, current.precision);
    const auto action = control.local_service ? control.local_service(control.local_service_data, current)
                                             : ResearchServiceAction::Continue;
    if (action == ResearchServiceAction::StopNumerical)
        control.local_service_stopped = true;
    return action;
}

void SolveState::research_start(const cps_options &options) {
    auto &s = *impl_;
    if (s.research)
        throw std::logic_error("private research driver is already started");
    if (s.program_route || s.secular_input || options.real_only)
        throw std::invalid_argument("private stepping requires a complex coefficient query");
    s.configure(&options);
    s.cancel_observed = false;
    s.research = std::make_unique<Impl::ResearchControl>();
    auto &r = *s.research;
    r.owner = &s;
    r.options = s.settings;
    r.callback = options.cancel;
    r.data = options.cancel_data;
    r.statistics_units_origin = s.result.statistics_v2.work_units;
    r.options.cancel = Impl::ResearchControl::poll;
    r.options.cancel_data = &r;
    if (s.policy.research_local_packet_meter) {
        s.policy.research_local_packet_service = Impl::ResearchControl::service;
        s.policy.research_local_packet_data = &r;
    }
    s.result.progress.work_units_this_call = 0;
    s.result.progress.stop_reason = CPS_STOP_NONE;
    s.prepare_coefficient_query(r.options);
}

void SolveState::research_set_local_packet_service(ResearchLocalPacketService callback, void *data) {
    auto &s = *impl_;
    if (!s.research || s.research->advancing || s.research->finished || s.preprocessed ||
        s.result.progress.work_units_this_call || !s.policy.research_local_packet_meter)
        throw std::logic_error("local service must be installed before an enabled private driver's first turn");
    auto &r = *s.research;
    r.local_service = callback;
    r.local_service_data = data;
    s.policy.research_local_packet_service = Impl::ResearchControl::service;
    s.policy.research_local_packet_data = &r;
}

bool SolveState::research_step(uint64_t remaining_work, uint64_t cell_limit, double subdivision_allowance) {
    auto &s = *impl_;
    if (!s.research) throw std::logic_error("research_start must precede a turn");
    auto &r = *s.research;
    if (r.advancing)
        throw std::logic_error("private numerical state cannot be reentered by its service callback");
    if (r.finished) return true;
    struct Guard {
        bool &flag;
        explicit Guard(bool &value) : flag(value) { flag = true; }
        ~Guard() { flag = false; }
    } guard(r.advancing);
    if (!(subdivision_allowance >= 0))
        throw std::invalid_argument("invalid subdivision work allowance");
    r.subdivision_allowance = subdivision_allowance;
    r.subdivision_turn = {};
    auto &o = r.options;
    // Cumulative absolute limit; a nested action may overrun it just as in the
    // one-shot driver. Never reset a counter or give either engine a new budget.
    o.work_units = std::min(s.settings.work_units,
        resources::add(s.result.progress.work_units_this_call, remaining_work));
    o.live_cells = std::min(s.settings.live_cells, cell_limit);
    bool stopped = s.terminal || s.result.goal_met;
    if (!stopped && o.cancelled()) {
        s.result.progress.stop_reason = CPS_STOP_CANCELLED;
        stopped = true;
    }
    if (!stopped && s.result.progress.work_units_this_call >= o.work_units) {
        s.result.progress.stop_reason = CPS_STOP_WORK_BUDGET;
        stopped = true;
    }
    if (!stopped) {
        if (!s.degree_proved) {
            stopped = s.prove_input_degree(o) == Impl::DriverStep::Stop;
        } else if (!s.preprocessed) {
            TraceScope trace(&s.policy, "preprocessing", s.initial_precision,
                static_cast<uint64_t>(s.result.degree), &s.result.statistics_v3, &s.result.statistics_v2);
            ++s.result.statistics_v2.work_units;
            ++s.result.progress.work_units_this_call;
            s.preprocess(o);
        } else {
            stopped = s.advance_factor(o, r.pass) == Impl::DriverStep::Stop;
        }
    }
    if (stopped || s.terminal || s.result.goal_met)
        research_finish(s.result.progress.stop_reason);
    return r.finished;
}

const cps_result &SolveState::research_result() const { return impl_->result; }

const cps_result &SolveState::research_finish(cps_stop_reason reason) {
    auto &s = *impl_;
    if (!s.research) throw std::logic_error("research_start must precede finish");
    if (s.research->servicing)
        throw std::logic_error("service callback cannot finalize its numerical state");
    if (!s.research->finished) {
        if (s.result.progress.stop_reason == CPS_STOP_NONE)
            s.result.progress.stop_reason = reason;
        s.finish_coefficient_query(s.research->options);
        // The private latch only provides the existing safe unwind path. It
        // is not an external cancellation or deadline; keep that distinction
        // in the public-shaped research result and the separate accessor.
        if (s.research->local_service_stopped && !s.cancel_observed &&
            !s.research->options.interrupted() && !s.result.goal_met) {
            s.result.termination = CPS_WORK_LIMIT;
            s.result.progress.stop_reason = CPS_STOP_WORK_BUDGET;
        }
        s.summarize_factor_progress();
        s.describe_coefficient_query();
        sort_roots(s.result.roots);
        s.research->finished = true;
    }
    return s.result;
}

uint64_t SolveState::research_live_cells() const { return impl_->live_cells(); }

uint64_t SolveState::research_inflight_work_units() const { return impl_->research_inflight_work_units(); }

bool SolveState::research_local_service_stopped() const {
    return impl_->research && impl_->research->local_service_stopped;
}

double SolveState::research_native_work() const {
    double work = 0;
    for (const auto &f : impl_->factors) {
        if (f->structured) work += f->structured->work();
        if (f->search) work += f->search->work();
        if (f->dovetail_search) work += f->dovetail_search->work();
    }
    return work;
}

bool SolveState::research_handoff_requested() const {
    return impl_->research && impl_->research->handoff_requested;
}

uint64_t SolveState::research_engine_bits() const {
    uint64_t bits = impl_->research ? impl_->research->engine_peak_bits : 0;
    for (const auto &f : impl_->factors)
        bits = std::max(bits, f->search ? f->search->max_bits() : f->precision);
    return bits;
}

SolveState::ResearchCheckpoint SolveState::research_checkpoint() const {
    const auto &s = *impl_;
    ResearchCheckpoint out;
    out.prepared = s.preprocessed;
    out.factors = s.factors.size();
    out.warm_eligible = s.research && s.preprocessed && !s.terminal && !s.result.goal_met &&
        !s.factors.empty() && !s.cancel_observed;
    for (const auto &f : s.factors) {
        out.ready_factors += f->ready;
        out.warm_eligible = out.warm_eligible && f->ready && f->kind != FactorState::Ball && f->power == 1;
    }
    return out;
}

SolveState::ResearchSubdivisionTurn SolveState::research_subdivision_turn() const {
    return impl_->research ? impl_->research->subdivision_turn : ResearchSubdivisionTurn{};
}

std::unique_ptr<SolveState> SolveState::research_warm_fork(uint64_t fallback_bits,
    void (*trace)(void *, const TraceEvent &), void *trace_data) const {
    const auto &s = *impl_;
    if (s.research && s.research->advancing)
        throw std::logic_error("cannot fork an in-flight numerical state");
    // The old research fork owns an independent persistent state. M4 instead
    // keeps its search in this query; its synchronous budget borrow must not
    // escape into a child that can outlive the owning query or renew the prefix.
    if (s.settings.heuristic_budget)
        throw std::logic_error("bounded heuristic queries retain their owned subdivision search");
    if (!research_checkpoint().warm_eligible || s.research->options.cancelled()) return nullptr;
    if (fallback_bits < 64 || fallback_bits > resources::precision_max)
        throw std::invalid_argument("invalid warm fallback precision cap");
    // One initial cell per certified region is an upper bound: the
    // constructor may merge overlapping support squares. Account for the
    // still-live parent before copying any factor/search storage.
    uint64_t required_cells = s.live_cells();
    for (const auto &f : s.factors)
        required_cells = resources::add(required_cells, f->cover.size());
    if (required_cells > s.research->options.live_cells) return nullptr;
    auto policy = s.policy;
    policy.backend = policy.dovetail_subdivision = policy.research_no_subdivision = false;
    policy.prepared_subdivision = false; // continue the requested goal from the cover
    policy.research_warm_subdivision = true;
    policy.research_integrated_dovetail = false;
    policy.research_scheduling_metrics = nullptr;
    policy.research_local_packet_meter = false;
    policy.research_local_packet_service = nullptr;
    policy.research_local_packet_data = nullptr;
    policy.trace = trace;
    policy.trace_data = trace_data;
    policy.snapshot = {};
    policy.study_statistics = nullptr;
    auto child = std::make_unique<SolveState>(s.input, std::move(policy));
    auto options = s.settings;
    options.max_work_bits = fallback_bits;
    // Borrow the external latch, never the parent's ResearchControl.
    options.cancel = s.research->callback;
    options.cancel_data = s.research->data;
    child->research_start(options);
    auto &c = *child->impl_;
    c.result = s.result;
    c.result.statistics = {};
    c.result.statistics_v2 = {};
    c.result.statistics_v3 = {};
    c.result.progress.work_units_this_call = 0;
    c.result.progress.stop_reason = CPS_STOP_NONE;
    c.result.termination = CPS_PRECISION_LIMIT;
    c.zero_order = s.zero_order;
    c.degree_proved = c.preprocessed = true;
    c.preprocessing_status = "copied exact prepared factors and certified checkpoint";
    for (const auto &f : s.factors) {
        auto g = std::make_unique<FactorState>();
        g->id = f->id;
        g->kind = f->kind;
        g->integer = f->integer;
        g->gaussian = f->gaussian;
        g->squarefree_proved = f->squarefree_proved;
        g->multiplicity = f->multiplicity;
        g->real_count = f->real_count;
        g->separation_bits = f->separation_bits;
        g->overlap_tightened = f->overlap_tightened;
        g->cover = f->cover;
        g->ready = true;
        g->precision = c.subdivision_start(options);
        g->search = c.make_research_warm_search(*g, options);
        c.factors.push_back(std::move(g));
    }
    return child;
}

bool SolveState::research_warm_handoff(uint64_t fallback_bits) {
    auto &s = *impl_;
    if (s.research && s.research->advancing)
        throw std::logic_error("cannot hand off an in-flight numerical state");
    if (!s.research || !s.preprocessed || s.terminal || s.result.goal_met ||
        s.research->options.cancelled() || s.factors.empty())
        return false;
    // The checkpoint's certified factor covers are in the original factor
    // coordinates. Power maps and refused exact preparation are deliberately
    // ineligible; never reinterpret numerical points as root-count proofs.
    for (const auto &f : s.factors)
        if (!f->ready || f->kind == FactorState::Ball || f->power != 1)
            return false;
    uint64_t cells = 0;
    for (const auto &f : s.factors) cells = resources::add(cells, f->cover.size());
    if (cells > s.research->options.live_cells) return false;
    if (fallback_bits) {
        if (fallback_bits < 64 || fallback_bits > resources::precision_max)
            throw std::invalid_argument("invalid warm fallback precision cap");
        s.settings.max_work_bits = s.research->options.max_work_bits = fallback_bits;
    }
    s.policy.research_no_subdivision = false;
    s.policy.research_warm_subdivision = true;
    s.policy.research_local_packet_meter = false;
    s.policy.research_local_packet_service = nullptr;
    s.policy.research_local_packet_data = nullptr;
    s.research->local_service = nullptr;
    s.research->local_service_data = nullptr;
    s.research->engine_peak_bits = 0;
    s.policy.dovetail_subdivision = false;
    for (auto &f : s.factors) {
        f->search = s.make_research_warm_search(*f, s.settings);
        f->structured.reset();
        f->dovetail_search.reset();
        f->candidates = CandidateState{};
        f->blocked = f->input_limited = f->iteration_wait = f->local_refinement = false;
        f->precision = s.subdivision_start(s.settings);
    }
    s.research->finished = false;
    s.research->pass.restart();
    s.result.progress.stop_reason = CPS_STOP_NONE;
    s.result.termination = CPS_PRECISION_LIMIT;
    return true;
}
#endif

} // namespace certpolysolve::detail
