// SolveState::Impl, the global driver's state for one advance (W8's S6: the class out of solver.cpp, its member
// functions defined in the files of their parts). Private to them: nothing else includes it. Its members, by
// part: the native secular route (native_advance) in
// solver_native.cpp, query-level phases in solver_advance.cpp, and factor preparation, steps and assembly in
// solver.cpp; goal and observe_precision, which several parts call, stay defined here. What the
// class does: solver.cpp's header.
#pragma once

#include "core/cover.hpp"
#include "core/engine_policy.hpp"
#include "core/query.hpp"
#include "input/coefficients/polynomial_input.hpp"
#include "numeric/exact/input_structure.hpp"
#include "solve/real_line_routing.hpp"
#include "solve/solve_state.hpp"
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <string>
#include <vector>
#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
#include "solve/factor_pass.hpp"
#endif

namespace certpolysolve::detail {

// Owned through pointers; their complete states belong in the implementation
// files that use them. Impl's constructors and destructor are out of line.
struct FactorState;
class Subdivision;
class SecularRefinement;
class ProgramRouteState;
class FactorPass;

struct SolveState::Impl {
    DensePolynomial input;
    std::unique_ptr<cps_secular> secular_input;
    // The all-real seed pass's outcome on the input before routing (hint_real_line_seeds): taken by the factor equal
    // to its polynomial at preprocessing.
    std::optional<RealLineSeedHint> real_line_hint;
    bool secular_converted = false;
    bool native_fallback = false;
    std::unique_ptr<SecularRefinement> native_secular;
    std::unique_ptr<ProgramRouteState> program_route;
    EnginePolicy policy;
    std::optional<InputStructure> source_structure;
    cps_options settings;
    cps_result result;
    std::vector<std::unique_ptr<FactorState>> factors;
    // The exact zero roots of coefficient input that preprocess sets apart (S8's B1): the point 0, counted
    // zero_order times, goes into every assembly beside the factors (aggregate). It stays 0 for a real-only query,
    // whose engine counts x among its factors, and for native secular input.
    uint64_t zero_order = 0;
    bool configured = false, degree_proved = false, preprocessed = false, terminal = false;
    bool cancel_observed = false;
    std::optional<HeuristicBudget> heuristic_budget;
    std::unique_ptr<Subdivision> make_research_warm_search(const FactorState &, const cps_options &);
#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
    struct ResearchControl {
        Impl *owner;
        cps_options options;
        CancelCallback callback = nullptr;
        void *data = nullptr;
        FactorPass pass;
        bool finished = false, handoff_requested = false;
        bool advancing = false, servicing = false, local_service_stopped = false;
        ResearchLocalPacketService local_service = nullptr;
        void *local_service_data = nullptr;
        uint64_t statistics_units_origin = 0;
        uint64_t engine_peak_bits = 0;
        double subdivision_allowance = std::numeric_limits<double>::infinity();
        SolveState::ResearchSubdivisionTurn subdivision_turn;
        static int poll(void *pointer) {
            auto &control = *static_cast<ResearchControl *>(pointer);
            if (!control.owner->cancel_observed && control.callback && control.callback(control.data))
                control.owner->cancel_observed = true;
            return control.owner->cancel_observed || control.local_service_stopped ? 1 : 0;
        }
        static ResearchServiceAction service(void *, const ResearchLocalPacketProgress &);
    };
    std::unique_ptr<ResearchControl> research;
    uint64_t research_inflight_work_units() const {
        const auto origin = research ? research->statistics_units_origin : 0;
        const auto charged = result.statistics_v2.work_units;
        return std::max(result.progress.work_units_this_call, charged >= origin ? charged - origin : 0);
    }
    cps_stop_reason research_handoff(cps_stop_reason reason) {
        if (research) research->handoff_requested = true;
        return reason;
    }
#endif
    uint64_t initial_precision = 64;
    size_t cursor = 0;
    std::string preprocessing_status;
    Impl(DensePolynomial p, EnginePolicy e);
    Impl(const cps_secular &p, EnginePolicy e);
    Impl(const cps_program &p, EnginePolicy e);
    ~Impl();
    // The solve's settings, fixed for its life: a solve state is advanced once (one-shot solves only, and the
    // real-only route's refinements are fresh solves: docs/session-state-design-note-20261017.md). History:
    // docs/code-history-20261021.md, solver.cpp 2.
    void configure(const cps_options *proposed);
    bool goal(const std::vector<Root> &cover, uint64_t precision) const {
        return std::all_of(cover.begin(), cover.end(), [&](const Root &r) {
            return meets_goal(r, settings, static_cast<slong>(precision));
        });
    }
    cps_options factor_options(const FactorState &f, const cps_options &base) const;
    bool factor_goal(const FactorState &f) const;
    // The input's or a factor's polynomial enclosed at `precision`; `observe` records the precision as the solve's
    // working precision (observe_precision), which a stop cover's enclosure does not (stop_cover).
    void enclose_input(acb_poly_t p, slong precision, bool observe = true);
    void observe_precision(slong precision) {
#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
        if (research)
            research->engine_peak_bits = std::max(research->engine_peak_bits, static_cast<uint64_t>(precision));
#endif
        result.work_bits = static_cast<uint64_t>(precision);
        result.progress.max_outer_bits = std::max(result.progress.max_outer_bits, result.work_bits);
        ++result.statistics.precision_attempts;
    }
    void enclose_factor(FactorState &f, acb_poly_t p, slong precision, bool observe = true);
    void preprocess(const cps_options &o);
    // Subdivision starts where every method starts: the working precision of
    // the factor's current goal, never below the precision that proved the
    // degree. Each component raises only as its own tests need; the precision
    // another method reached measures that method's needs, not counting's.
    uint64_t subdivision_start(const cps_options &o) const;
    // The candidate iteration has spent the precision budget: the complete
    // method continues from the certified cover. It starts at the goal's
    // working precision like every handoff (subdivision_start), not at the
    // budget the candidates reached: each component raises as its own tests
    // need, and most need far less.
    void candidate_handoff(FactorState &f, const cps_options &o);
    uint64_t live_cells() const;
    // Set by a factor step that provably left every factor's cover as it was.
    bool step_kept_covers = false;
    // The factors' covers assembled and published (solver.cpp). A factor that has no cover yet stops the assembly,
    // except at a stop short of the goal (`stop`, the call's options; advance passes them after its loop when the
    // stop is not a cancellation): the factor then contributes its stop cover. A cancellation during that completion
    // publishes nothing.
    void aggregate(const cps_options *stop = nullptr);
    // A factor's stop cover: complete_cover_at_stop on the factor's own polynomial in x (a power composite's reduced
    // polynomial spread back to x^power, so no lift of a region about the origin is needed), at the factor's
    // precision and at least 64 bits; its counts and multiplicities times the factor's multiplicity. False where the
    // leading coefficient's enclosure contains zero (no root bound) or the completion was cancelled.
    bool stop_cover(FactorState &f, const cps_options &o, std::vector<Root> &cover);
    cps_stop_reason factor_step(FactorState &f, const cps_options &base);
    // factor_step's steps (solver.cpp), in the order it tries them after a power's goal and a blocked factor: the
    // real-only query (RealRefinement); the research build's Portfolio share of Subdivision; after a pending local
    // refinement, a factor Subdivision owns; a structured factor's dovetailed turn (DovetailSchedule: Subdivision on
    // its metered share; no value when the turn falls through) and else its iteration (SecularRefinement::step, and at
    // its precision or input stop the handoff to Subdivision); and else the candidate stage (the structured entry on a
    // PolynomialSource, or the classical pass: exact centring, isolate_complex, the cover, the precision doubled).
    cps_stop_reason real_only_step(FactorState &f, const cps_options &o);
#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
    cps_stop_reason portfolio_subdivision_step(FactorState &f, const cps_options &o);
#endif
    cps_stop_reason subdivision_step(FactorState &f, const cps_options &o);
    cps_stop_reason start_retained_subdivision(FactorState &f, const cps_options &o);
    cps_stop_reason start_prepared_subdivision(FactorState &f, const cps_options &o);
    cps_stop_reason start_complete_fallback(FactorState &f, const cps_options &o);
    std::optional<cps_stop_reason> dovetail_turn(FactorState &f, const cps_options &o);
    cps_stop_reason structured_step(FactorState &f, const cps_options &o);
    cps_stop_reason candidate_step(FactorState &f, const cps_options &o);
    bool native_advance(const cps_options &o);
    // Query-level phases (solver_advance.cpp). Stop ends scheduling, not
    // publication: finish_coefficient_query still resolves proofs and limits.
    enum class DriverStep { Continue, Stop };
    // True when a special input route has produced the call's final result.
    bool advance_input_route(const cps_options &o);
    void prepare_coefficient_query(const cps_options &o);
    DriverStep prove_input_degree(const cps_options &o);
    DriverStep advance_factor(const cps_options &o, FactorPass &pass);
    DriverStep finish_factor_pass(const cps_options &o, FactorPass &pass);
    void finish_coefficient_query(const cps_options &o);
    void summarize_factor_progress();
    void describe_coefficient_query();
    cps_result advance(const cps_options *proposed);
};

} // namespace certpolysolve::detail
