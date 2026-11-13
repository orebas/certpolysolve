#include "solve/local_model_solver.hpp"
// The program route (ProgramRouteState and its private controller): the degree, leading
// coefficient and a root bound proved (prepare_program); a product of distinct factors solved by nested solves and a
// power c q^k through q; seeds placed by recurrence_seeds or program_seeds; the binary64 stage and ball packets of
// Aberth sweeps until the points converge or stall; then SeparatedRootsPass, and SecularRefinement on a ProgramSource
// where that pass fails. The finish labels the published regions: non-real off the axis, and real by conjugation
// (prove_real_by_reflection) when every constant of the program is real. Design: docs/slp-input-spec-20261015.md,
// docs/slp-second-version-note-20261017.md, docs/slp-structural-seeding-note-20261018.md.
#include "solve/program_route_state.hpp"
#include "numeric/exact/program_structure.hpp"
#include "core/cover.hpp"
#include "solve/solve_state.hpp"
#include "input/secular/secular.hpp"
#include "core/statistics_accumulation.hpp"
#include "core/telemetry.hpp"
#include <algorithm>
#include <cmath>
#include <locale>
#include <sstream>

namespace certpolysolve::detail {
class ProgramRouteState::Controller {
    ProgramRouteState &route;
    const cps_options &o;
    const EnginePolicy &policy;
    cps_result &result;
    const uint64_t initial_precision;
    const Program &program;
    const bool aberth;

    struct FinalStatus {
        cps_stop_reason stop;
        cps_termination termination;
        bool precision_stop;
        FinalStatus(cps_stop_reason reason, cps_termination ending, bool limited)
            : stop(reason), termination(ending), precision_stop(limited) {}
    };
    // A preliminary phase always ends this work unit, even when it converges.
    enum class Step { NextUnit, Stop, IterationLimit };
    // A failed disk pass hands over to refinement within the same charged unit.
    enum class DiskStep { NextUnit, RefineNow };

    const Program &solved() const { return route.base ? *route.base : program; }
    bool goal(const std::vector<Root> &cover, uint64_t precision) const {
        return std::all_of(cover.begin(), cover.end(),
                           [&](const Root &r) { return meets_goal(r, o, static_cast<slong>(precision)); });
    }
    void observe_precision(slong precision) {
        result.work_bits = static_cast<uint64_t>(precision);
        result.progress.max_outer_bits = std::max(result.progress.max_outer_bits, result.work_bits);
        ++result.statistics.precision_attempts;
    }

    void describe(std::ostringstream &) const;
    void with_zero(std::vector<Root> &) const;
    bool program_goal(const std::vector<Root> &, uint64_t precision) const;
    void finish(FinalStatus);
    Root global_cover() const;
    std::optional<FinalStatus> prepare();
    std::optional<FinalStatus> product_lift();
    void power_lift();
    bool seed_source();
    bool derivatives_unreliable(bool unsettled_only, slong bits) const;
    void raise_precision(const char *phase);
    bool raise_until_readable(const char *phase);
    std::optional<Step> machine_step();
    Step ball_step();
    Step preliminary_step();
    DiskStep disk_step();
    cps_stop_reason refinement_step();
    void update_cover();
    FinalStatus iterate();

  public:
    Controller(ProgramRouteState &state, const cps_options &options, const EnginePolicy &engine, cps_result &answer,
               uint64_t precision)
        : route(state), o(options), policy(engine), result(answer), initial_precision(precision),
          program(state.input.program), aberth(engine.refinement == RefinementMethod::Aberth) {}
    void run();
};

int64_t ProgramRouteState::degree_bound() const {
    const auto &program = input.program;
    if (!program.output || *program.output >= program.nodes.size())
        throw std::invalid_argument("program has no output");
    // The upper bound until preparation proves the degree.
    return static_cast<int64_t>(
        std::min<uint64_t>(program.degree_bounds()[*program.output], std::numeric_limits<int64_t>::max()));
}

void ProgramRouteState::solve(const cps_options &options, const EnginePolicy &policy, cps_result &result,
                              uint64_t initial_precision) {
    Controller(*this, options, policy, result, initial_precision).run();
}

void ProgramRouteState::Controller::run() {
    auto status = prepare();
    if (!status)
        status = iterate();
    // Exactly one publication applies the zero root and power lift together.
    finish(*status);
}

void ProgramRouteState::Controller::describe(std::ostringstream &trace) const {
    trace << "program input: " << program.nodes.size() << " operations, " << (aberth ? "Aberth" : "secular")
          << " mode; ";
    if (!route.preparation.ready()) {
        trace << route.preparation.diagnostic;
        return;
    }
    if (!route.product_note.empty())
        trace << route.product_note << "; ";
    if (route.lift.exponent > 1)
        trace << "the output is a power c q^" << route.lift.exponent << " of node " << route.lift.base << ", degree "
              << route.lift.degree << ": q solved (its cone, " << route.lift.nodes
              << " operations) and its cover lifted, counts times " << route.lift.exponent << "; q's ";
    trace << "degree " << route.preparation.degree << " and leading coefficient proved by the top "
          << "coefficients at width " << route.preparation.top_width;
    if (route.preparation.degree && route.preparation.monomial)
        trace << "; modular image: " << route.preparation.modular.note
              << "; every root zero (a monomial, proved by the program's majorant or its lowest coefficients)";
    else if (route.preparation.degree)
        trace << "; modular image: " << route.preparation.modular.note << "; root bound below 2^"
              << std::ceil(mag_get_d_log2_approx(route.preparation.root_bound.value()));
    if (route.preparation.zero_order)
        trace << "; a root of order " << route.preparation.zero_order
              << " at zero (the program's lowest coefficients), the rest solved as q = p / x^"
              << route.preparation.zero_order;
    if (const auto &seeds = route.preliminary.seeds; seeds.circles || seeds.evaluations)
        trace << "; seeds: " << seeds.placed << " of " << route.preparation.degree - route.preparation.zero_order
              << " at the roots' radii from " << seeds.circles << " circle counts (" << seeds.binary64_circles
              << " in binary64, " << seeds.unsharp << " unsharp, " << seeds.evaluations << " evaluations, "
              << seeds.annuli.size() << " annuli, " << seeds.raises << " precision raises, " << seeds.sample_doublings
              << " sample doublings, " << seeds.rotations << " rotations"
              << (seeds.budget_spent ? ", budget spent" : "") << "), " << seeds.precision << " bits";
    if (const auto &r = route.preliminary.recurrence; r.levels)
        trace << "; structural seeds: " << r.levels << " levels (" << r.compositions
              << " by a composition's preimages), from the roots of a base of degree " << r.base_degree << ", "
              << r.seeded << " seeded, " << r.leftover << " left over, " << r.sub_sweeps << " sweeps below"
              << (route.preparation.zero_order
                      ? ", the " + std::to_string(route.preparation.zero_order) + " nearest zero dropped"
                      : std::string());
    if (const auto &m = route.preliminary; m.machine_refused)
        trace << "; binary64 stage refused";
    else if (o.max_work_bits < static_cast<uint64_t>(machine::binary64_bits))
        trace << "; binary64 stage not run under a cap below 53 bits";
    else if (m.machine_sweeps || m.machine_done)
        trace << "; binary64 stage: " << m.machine_sweeps << " sweeps" << (m.machine_done ? "" : " so far") << ", "
              << m.machine_converged << " converged, " << m.machine_floor << " at the floor, " << m.machine_held
              << " out of range" << (m.machine_stalled ? ", stalled" : "")
              << (m.scaled_evaluations ? ", " + std::to_string(m.scaled_evaluations) + " scaled evaluations" : "");
    if (route.preliminary.sweeps || route.preliminary.machine_sweeps)
        trace << "; preliminary: " << route.preliminary.sweeps << " sweeps" << (route.preliminary.done ? "" : " so far")
              << ", " << route.preliminary.raises << " precision raises, " << route.preliminary.precision << " bits";
    if (const auto *pass = route.disks.pass.get()) {
        const auto state = pass->state();
        trace << "; Newton disks: "
              << (state == SeparatedRootsPass::State::Certified ? "a complete disjoint cover"
                  : state == SeparatedRootsPass::State::Failed  ? "no cover (" + pass->failure() + "), the refinement"
                                                                : std::string("running"))
              << " after " << pass->rounds() << " rounds (" << pass->evaluations() << " evaluations, " << pass->steps()
              << " Newton steps, " << pass->rehomed() << " re-homed, " << pass->raises() << " raises, "
              << pass->overlapping() << " overlapping, " << pass->wide() << " wide), " << pass->precision() << " bits";
    }
}

// A cover of q with the root of order v at zero, as the finish publishes it: the point 0, count and multiplicity
// v (proved by the program's lowest coefficients), or, when one of q's regions holds zero, that region counting
// the v roots too, its multiplicity unknown.
void ProgramRouteState::Controller::with_zero(std::vector<Root> &cover) const {
    const uint64_t v = route.zero_order;
    for (auto &root : cover)
        if (acb_contains_zero(root.box.value())) {
            root.count += v;
            root.multiplicity = 0;
            root.certificate += "; with the root of order " + std::to_string(v) + " at zero";
            return;
        }
    Root zero; // the exact point 0
    zero.count = zero.multiplicity = v;
    zero.certificate = "program-zero-root";
    zero.publish(static_cast<slong>(std::max<uint64_t>(64, result.work_bits)));
    cover.push_back(std::move(zero));
}

// The goal of the program route's covers. While a root at zero waits for its point (route.zero_order), a cover of q
// (its counts summing to q's degree) is judged as it will be published (with_zero), and the global cover, of the
// degree n, which holds zero, is no answer. Refusing a cover of q with a region holding zero outright would
// deadlock: the refinement, its own goal met, would have nothing left to do, and under the cluster goal the
// solve would sweep to its budget. History: docs/code-history-20261021.md, solver.cpp 9.
bool ProgramRouteState::Controller::program_goal(const std::vector<Root> &cover, uint64_t precision) const {
    if (!route.zero_order)
        return goal(cover, precision);
    uint64_t total = 0;
    for (const auto &root : cover)
        total += root.count;
    if (total != static_cast<uint64_t>(result.degree))
        return false;
    auto lifted = cover;
    with_zero(lifted);
    return goal(lifted, precision);
}

void ProgramRouteState::Controller::finish(FinalStatus status) {
    const auto [stop, termination, precision_stop] = status;
    // A root of order v at zero, before the power's lift (so the lift multiplies it too): q's cover, of q's
    // degree n - v, published with the point 0 (with_zero). A cover of the degree n (the global cover, which
    // holds zero) stays as it is.
    if (route.zero_order) {
        const uint64_t v = route.zero_order, q_degree = static_cast<uint64_t>(result.degree);
        uint64_t total = 0;
        for (const auto &root : result.roots)
            total += root.count;
        if (!result.roots.empty() && total == q_degree)
            with_zero(result.roots);
        // goal_met stays as the loop left it: the exact point meets every goal, so a cover that met it still
        // does, and one that did not (a stop) must not turn into a success here.
        result.degree = static_cast<int64_t>(q_degree + v);
    }
    // A power output's cover: q's regions, each holding its roots k times as often. A multiplicity is
    // multiplied only where q's region proves one (an unknown multiplicity stays unknown); the goal is
    // judged again on the lifted regions.
    if (route.lift.exponent > 1) {
        const auto k = route.lift.exponent;
        for (auto &root : result.roots) {
            root.count *= k;
            root.multiplicity *= k;
            root.certificate += "; lifted from the base q to the output c q^" + std::to_string(k);
        }
        result.degree = static_cast<int64_t>(route.lift.degree);
        // Lifting keeps a goal met or unmet (counts and multiplicities scale together, the regions stay), so
        // this only confirms a met one: an answer at a stop never turns into a success here.
        if (result.complete && result.goal_met)
            result.goal_met = goal(result.roots, std::max<uint64_t>(64, result.work_bits));
    }
    // The program contract (certpolysolve.h): a region off the real axis is reported non-real, whichever path made
    // it. The refinement's clusters counted on its model (secular-model-count) came without the label (the
    // program fuzz family's program_label_unresolved).
    for (auto &root : result.roots)
        if (root.realness == CPS_REAL_UNKNOWN && !arb_contains_zero(acb_imagref(root.box.value())))
            root.realness = CPS_NONREAL;
    // A region meeting the axis is reported real by conjugation (since 23 October 2026). A program whose
    // constants are all real has real coefficients (sums, differences and products of real numbers), so its
    // non-real roots come in conjugate pairs. On the complete disjoint cover as published (the point 0 and the
    // power's lift applied), a region holding one distinct root (count 1, or a proved multiplicity equal to its
    // count) whose reflection in the real axis meets no other region holds that root's conjugate too: the root
    // is real (prove_real_by_reflection, the coefficient route's proof). Such a region becomes its real segment,
    // labelled real, its certificate marked "+conjugation": its imaginary ends become 0 and its real ends stay as
    // published (keep_real_ends), so the regions the loop judged only shrink and a met goal stays met, and a stop
    // stays a stop. Any other region meeting the axis is left unknown, and so is every region of a program with a
    // non-real constant. The proof runs in the outermost solve only: a factor that the product lift solves in a
    // nested state (o.program_lift_factor) keeps its regions as they were, since the lift unites, merges and
    // judges the factors' covers (a factor's region shrunk to its real segment would merge less at a stop).
    if (!o.program_lift_factor && result.complete && !result.roots.empty() &&
        std::all_of(program.constants.begin(), program.constants.end(),
                    [](const Coefficient &c) { return c.is_real(); })) {
        uint64_t total = 0;
        for (const auto &root : result.roots)
            total += root.count;
        if (total == static_cast<uint64_t>(result.degree) && disjoint(result.roots))
            prove_real_by_reflection(result.roots, static_cast<slong>(std::max<uint64_t>(64, result.work_bits)), true);
    }
    result.progress.stop_reason = stop;
    result.termination = termination;
    result.progress.pending_components = result.goal_met ? 0 : route.refinement ? route.refinement->active() : 1;
    result.progress.min_pending_bits = result.progress.max_pending_bits = result.goal_met ? 0 : result.work_bits;
    result.components.clear();
    for (size_t i = 0; i < result.roots.size(); ++i)
        result.components.push_back(
            {0, result.work_bits, 0, result.roots[i].count,
             meets_goal(result.roots[i], o, static_cast<slong>(std::max<uint64_t>(64, result.work_bits)))
                 ? CPS_COMPONENT_FINISHED
                 : (precision_stop ? CPS_COMPONENT_PRECISION_LIMIT : CPS_COMPONENT_ACTIVE)});
    std::ostringstream trace;
    trace.imbue(std::locale::classic()); // diagnostics never follow the host's locale
    describe(trace);
    trace << "; " << result.statistics_v2.work_units << " work units; " << result.progress.pending_components
          << " pending components; stop=" << static_cast<int>(stop);
    if (precision_stop) {
        trace << "; " << precision_cap_note << o.max_work_bits << " bits; "
              << (route.preparation.modular.squarefree
                      ? "distinct roots, unresolved at the cap"
                      : "possibly a multiple root (the modular image did not prove the polynomial square-free)");
    }
    result.diagnostics = trace.str();
    sort_roots(result.roots);
}

// The global cover: the square about the origin of half-side R holds every root (the majorant bound of
// section 7 is a proof), so its count is the degree. The call returns it until a finer certified cover
// exists, as the coefficient path returns Cauchy's. A monomial's (program_monomial) is the point 0, the
// root of multiplicity the degree: it meets every goal, and the solve ends with it.
Root ProgramRouteState::Controller::global_cover() const {
    Root global;
    RealBall radius;
    arf_set_mag(arb_midref(radius.value()), route.preparation.root_bound.value());
    arb_add_error(acb_realref(global.box.value()), radius.value());
    arb_add_error(acb_imagref(global.box.value()), radius.value());
    global.count = route.preparation.degree;
    global.multiplicity = route.preparation.monomial      ? route.preparation.degree
                          : route.preparation.degree == 1 ? 1u
                                                          : 0u;
    global.certificate = route.preparation.monomial ? "program-monomial" : "program-majorant-bound";
    global.publish(static_cast<slong>(initial_precision));
    return global;
}

// Preparation and structural solves share the first charged unit.
std::optional<ProgramRouteState::Controller::FinalStatus> ProgramRouteState::Controller::prepare() {
    if (result.progress.work_units_this_call == o.work_units || o.cancelled()) {
        return FinalStatus(o.cancelled() ? CPS_STOP_CANCELLED : CPS_STOP_WORK_BUDGET,
                           o.cancelled() ? CPS_CANCELLED : CPS_WORK_LIMIT, false);
    }
    ++result.statistics_v2.work_units;
    ++result.progress.work_units_this_call;
    {
        StageTimer timer(result.statistics.preprocessing_ns);
        TraceScope trace(&policy, "program_preparation", initial_precision,
                         static_cast<uint64_t>(std::max<int64_t>(result.degree, 0)), &result.statistics_v3,
                         &result.statistics_v2);
        route.preparation = prepare_program(program, o);
    }
    // The program's own modular image (C14): proved square-free, or nothing known. Read before a power lift
    // replaces the preparation with its base's, whose verdict is not the program's.
    result.squarefree = route.preparation.modular.squarefree ? CPS_SQUAREFREE_PROVED : CPS_SQUAREFREE_UNKNOWN;
    using Refusal = ProgramPreparation::Refusal;
    if (!route.preparation.ready()) {
        switch (route.preparation.refusal) {
        case Refusal::ZeroPolynomial:
            result.degree = -1;
            return FinalStatus(CPS_STOP_INPUT, CPS_ZERO_POLYNOMIAL, false);
        case Refusal::LeadingUndetermined:
            return FinalStatus(CPS_STOP_INPUT, CPS_INDETERMINATE_DEGREE, false);
        case Refusal::NoOutput:
        case Refusal::Length:
        case Refusal::ConstantStorage:
        case Refusal::None:
            result.progress.preprocessing_limited = 1;
            return FinalStatus(CPS_STOP_INPUT, CPS_UNSUPPORTED_QUERY, false);
        }
    }
    result.degree = static_cast<int64_t>(route.preparation.degree);
    if (result.degree == 0) {
        // A nonzero constant: no roots, the answer a constant coefficient input gets.
        result.complete = result.goal_met = true;
        return FinalStatus(CPS_STOP_NONE, CPS_SUCCESS, false);
    }
    if (auto status = product_lift())
        return status;
    power_lift();
    // The storage rule (section 3.1) at the solve's cap, before anything of the degree's size is allocated.
    if (!resources::storage_fits(program_solve_storage(route.preparation.degree, o.max_work_bits), o.temporary_bytes)) {
        result.progress.preprocessing_limited = 1;
        route.preparation.diagnostic = "the solve's storage for degree " + std::to_string(route.preparation.degree) +
                                       " exceeds the temporary storage allowed";
        route.preparation.refusal = Refusal::Length;
        return FinalStatus(CPS_STOP_INPUT, CPS_UNSUPPORTED_QUERY, false);
    }
    // The global cover (above), of q's preparation when the output is a power.
    result.roots = {global_cover()};
    result.complete = true;
    observe_precision(static_cast<slong>(initial_precision));
    // A root of order v at zero: the solve goes on with q = p / x^v, of degree n - v (route.zero_order). The global
    // cover, of degree n, holds the point 0; it stays until q's cover replaces it, and is no answer meanwhile.
    if (route.preparation.zero_order) {
        route.zero_order = route.preparation.zero_order;
        result.degree -= static_cast<int64_t>(route.zero_order);
    }
    result.goal_met = program_goal(result.roots, initial_precision);
    return std::nullopt;
}

// A product of distinct factors c q_1^k_1 ... q_m^k_m (output_factors): each factor's cone is solved alone,
// by a nested solve state advanced once with this call's remaining budgets (its statistics added to this
// solve's), and its cover lifted as a power's (counts times k_i, a multiplicity times k_i where the region
// proves one); the covers are united. When regions of different factors overlap (a root shared by two
// factors, which the pass's counted clusters answer), or the proved degrees do not add up, or a factor is
// refused or its cover incomplete, the output is solved whole, like any other program; the factors' costs
// stay charged. A factor stopped by a budget, a time limit or a cancellation ends the solve with its stop.
// Neither lift is needed for a monomial: its cover is final. History:
// docs/code-history-20261021.md, solver.cpp 10.
std::optional<ProgramRouteState::Controller::FinalStatus> ProgramRouteState::Controller::product_lift() {
    if (!heuristics::program_product_lift || route.preparation.monomial)
        return std::nullopt;
    const auto factors = output_factors(program);
    if (!factors)
        return std::nullopt;
    std::vector<Root> united;
    std::string notes;
    // The output's roots placed so far: those of the factors solved with a complete cover. Only such a
    // cover's regions are lifted, since a region of one factor counts the output's roots only when
    // every other factor's roots are placed too.
    uint64_t degrees = 0, placed = 0, solved = 0;
    bool complete = true, refused = false;
    cps_termination termination = CPS_SUCCESS;
    cps_stop_reason stop = CPS_STOP_NONE;
    for (const auto &factor : *factors) {
        cps_options part_options = o;
        part_options.work_units = o.work_units > result.progress.work_units_this_call
                                      ? o.work_units - result.progress.work_units_this_call
                                      : 0;
        part_options.work_units_default = false;
        part_options.program_lift_factor = true; // the outer finish proves realness
        const uint64_t sweeps_left = o.sweeps_left(result.statistics_v3);
        if (o.max_iterations)
            part_options.max_iterations = sweeps_left;
        if (o.cancelled()) {
            termination = CPS_CANCELLED, stop = CPS_STOP_CANCELLED, complete = false;
            break;
        }
        if (!part_options.work_units || (o.max_iterations && !sweeps_left)) {
            termination = CPS_WORK_LIMIT, complete = false;
            stop = part_options.work_units ? CPS_STOP_ITERATION_BUDGET : CPS_STOP_WORK_BUDGET;
            break;
        }
        cps_program part;
        part.program = cone(program, factor.base);
        SolveState nested(part, policy);
        const cps_result answer = nested.advance(&part_options);
        // Nested solves use ordinary unsigned sums, including their factor counts. Unlike an
        // engine handoff, this is aggregation of the factors of the answering preparation.
        statistics::add_ordinary(result.statistics, answer.statistics, statistics::carried_stage_fields);
        result.statistics.exact_factors += answer.statistics.exact_factors;
        statistics::add_ordinary(result.statistics_v2, answer.statistics_v2, statistics::search_fields);
        statistics::add_ordinary(result.statistics_v3, answer.statistics_v3, statistics::candidate_fields);
        result.progress.work_units_this_call += answer.progress.work_units_this_call;
        result.progress.max_outer_bits = std::max(result.progress.max_outer_bits, answer.progress.max_outer_bits);
        result.work_bits = std::max(result.work_bits, answer.work_bits);
        ++solved;
        notes += "; factor " + std::to_string(solved) + " (node " + std::to_string(factor.base) + ", power " +
                 std::to_string(factor.exponent) + "): " + answer.diagnostics;
        // A budget, a time limit or a cancellation stops the lift whatever an earlier factor ended with
        // (an earlier factor's precision limit is the lift's only when every factor was reached), and
        // a factor it stopped before its preparation proved a degree is stopped, not refused. Keeping
        // the first termination instead would, after a factor at the precision cap, let a later
        // factor's stop leave the factors after it unsolved and the output solved whole. History:
        // docs/code-history-20261021.md, solver.cpp 11.
        const bool factor_stopped = answer.termination == CPS_CANCELLED || answer.termination == CPS_WORK_LIMIT;
        if (!factor_stopped &&
            (answer.termination == CPS_ZERO_POLYNOMIAL || answer.termination == CPS_INDETERMINATE_DEGREE ||
             answer.termination == CPS_UNSUPPORTED_QUERY || answer.degree <= 0)) {
            refused = true;
            break;
        }
        if (answer.degree > 0) {
            degrees += factor.exponent * static_cast<uint64_t>(answer.degree);
            if (answer.complete) {
                placed += factor.exponent * static_cast<uint64_t>(answer.degree);
                for (auto root : answer.roots) {
                    root.count *= factor.exponent;
                    root.multiplicity *= factor.exponent;
                    root.certificate += "; lifted from factor " + std::to_string(solved) + " of the output's " +
                                        std::to_string(factors->size()) + ", to its power " +
                                        std::to_string(factor.exponent);
                    united.push_back(std::move(root));
                }
            }
        }
        complete = complete && answer.complete;
        if (factor_stopped) {
            termination = answer.termination, stop = answer.progress.stop_reason;
            break;
        }
        if (answer.termination != CPS_SUCCESS && termination == CPS_SUCCESS)
            termination = answer.termination, stop = answer.progress.stop_reason;
    }
    const bool stopped = termination == CPS_CANCELLED || termination == CPS_WORK_LIMIT;
    const auto header = "the output is a product of " + std::to_string(factors->size()) + " distinct factors";
    if (!refused && stopped) {
        // A budget, a time limit or a cancellation. A factor not reached, or the roots an incomplete
        // cover leaves out, may lie in any solved factor's region (a root shared with it), and those
        // regions published with the solved factors' counts alone would undercount such a shared
        // root. So the answer is the factors' covers, overlaps merged, only when they place every root
        // of the output; otherwise it is the global cover, the answer of a stop before any factor.
        // History: docs/code-history-20261021.md, solver.cpp 12.
        const auto prec = static_cast<slong>(std::max<uint64_t>(64, result.work_bits));
        const bool all_placed = placed == route.preparation.degree;
        if (all_placed) {
            merge_overlapping(united, prec);
            result.roots = std::move(united);
        } else {
            result.roots = {global_cover()};
        }
        result.complete = true;
        result.goal_met = goal(result.roots, static_cast<uint64_t>(prec));
        route.product_note =
            header + ", " + std::to_string(solved) + " of them solved before the stop" +
            (all_placed ? ", their covers united" : ", so the answer is the global cover (not every root placed)") +
            notes;
        return FinalStatus(stop, termination, false);
    }
    if (!refused && complete && degrees == route.preparation.degree && disjoint(united)) {
        result.roots = std::move(united);
        result.complete = complete;
        result.goal_met = complete && goal(result.roots, std::max<uint64_t>(64, result.work_bits));
        route.product_note = header + ", each solved alone and its cover lifted" + notes;
        return FinalStatus(stop, termination, termination == CPS_PRECISION_LIMIT);
    }
    // A program's answer is complete from its global cover on, so only a stop before a factor's
    // preparation leaves its cover incomplete, and every factor's stop is the lift's (above): what
    // reaches this point is a refusal, overlapping covers (a shared root) or, for a wrong cone, a
    // degree mismatch. The output is then solved whole.
    route.product_note = header +
                         (refused                               ? ": a factor refused"
                          : !complete                           ? ": a factor's cover is incomplete"
                          : degrees != route.preparation.degree ? ": the factors' degrees do not add up"
                                                                : ": their covers overlap (a shared root)") +
                         ", so the output is solved whole";
    return std::nullopt;
}

// A power output c q^k: the route solves q, the node's cone, and lifts its cover when it publishes it. The
// roots of c q^k are q's, each k times as often, so the lift is exact; it is taken only when q's own
// preparation proves the degree deg p / k (a mismatch would mean a wrong cone: the input is solved).
void ProgramRouteState::Controller::power_lift() {
    if (heuristics::program_power_lift && !route.preparation.monomial) {
        if (const auto power = output_power(program)) {
            auto base = cone(program, power->base);
            ProgramPreparation base_preparation;
            {
                StageTimer timer(result.statistics.preprocessing_ns);
                base_preparation = prepare_program(base, o);
            }
            if (base_preparation.ready() && base_preparation.degree > 0 &&
                route.preparation.degree % power->exponent == 0 &&
                base_preparation.degree == route.preparation.degree / power->exponent) {
                route.lift = {power->base, power->exponent, route.preparation.degree,
                              static_cast<uint64_t>(base.nodes.size())};
                route.base = std::move(base);
                route.preparation = std::move(base_preparation);
                result.degree = static_cast<int64_t>(route.preparation.degree);
            }
        }
    }
}

bool ProgramRouteState::Controller::seed_source() {
    const auto n = static_cast<size_t>(result.degree);
    route.source =
        std::make_unique<ProgramSource>(solved(), result.degree, route.preparation.leading, route.zero_order);
    // Seeds at the roots' radii, from root counts on circles inside the root bound (program_seeds).
    // Section 6's start, every seed on the bound's circle, would make the preliminary phase crawl
    // inward by (n - 1) / (n + 1) a sweep, about (n / 2) ln(R / r) sweeps. The placement is
    // uncertified, as the polygon's is for coefficients. The phase starts at the initial precision
    // whatever the counts needed: the counts read balls, and a recurrence's balls at its roots are
    // pessimistic (at |z| = 1 a three-term recurrence's radii follow it with |2z| + 1 while its
    // values stay O(1), so its midpoints converge at a precision far below the one its balls need to
    // read, where the same sweeps cost less), while Laguerre's raises are real cancellation at its
    // large roots, which the phase meets as noise. The counts' precision is where the phase's first
    // raise goes (below: the balls are read at the first window's end, not only at a stall, while the
    // phase is below it). History: docs/code-history-20261021.md, solver.cpp 13.
    TraceScope trace(&policy, "program_seeds", initial_precision, n, &result.statistics_v3, &result.statistics_v2);
    RootVector seeds(result.degree);
    // Structural seeds first (program_recurrence_seeds.hpp): a squaring recurrence's output seeded from
    // its base's roots, which a rough solve below finds (the same rule again, then the binary64 stage).
    // Proposals, as the counts' placement is. History: docs/code-history-20261021.md, solver.cpp 14.
    bool structural = false;
    // With a root of order v at zero the structure is p's: p's n seeds, less the v nearest zero, which stand
    // for the root q does not have. Coincident seeds are spread apart (spread_duplicates), so at most one
    // lies at exactly zero, where q cannot be evaluated, and it goes first. The Mandelbrot iteration
    // p_(k+1) = p_k^2 + x is such a program: every p_k vanishes at zero. History:
    // docs/code-history-20261021.md, solver.cpp 15.
    if (heuristics::program_recurrence_seeds && o.max_work_bits >= static_cast<uint64_t>(machine::binary64_bits)) {
        const uint64_t v = route.zero_order, p_degree = static_cast<uint64_t>(result.degree) + v;
        RootVector p_seeds(static_cast<slong>(p_degree));
        structural = recurrence_seeds(p_seeds.data(), solved(), p_degree, route.preparation.root_bound.value(),
                                      static_cast<slong>(initial_precision), o, result.statistics_v3,
                                      route.preliminary.recurrence, heuristics::program_recurrence_seed_degree);
        if (structural) {
            std::vector<std::pair<double, size_t>> by_size(p_degree);
            for (size_t i = 0; i < p_degree; ++i) {
                const acb_srcptr z = p_seeds.data() + static_cast<slong>(i);
                by_size[i] = {std::hypot(arf_get_d(arb_midref(acb_realref(z)), ARF_RND_NEAR),
                                         arf_get_d(arb_midref(acb_imagref(z)), ARF_RND_NEAR)),
                              i};
            }
            std::sort(by_size.begin(), by_size.end());
            std::vector<uint8_t> dropped(p_degree, 0);
            for (size_t i = 0; i < v; ++i)
                dropped[by_size[i].second] = 1;
            slong kept = 0;
            for (size_t i = 0; i < p_degree; ++i)
                if (!dropped[i])
                    acb_set(seeds.data() + kept++, p_seeds.data() + static_cast<slong>(i));
        }
    }
    if (!structural) {
        // The counts read the program in binary64 where they can (program_machine.hpp), balls otherwise;
        // with a root at zero, q's counts, as the source's (the evaluator's zero order).
        std::unique_ptr<ProgramMachine<double>> counting;
        if (heuristics::program_machine_stage && o.max_work_bits >= static_cast<uint64_t>(machine::binary64_bits)) {
            counting = std::make_unique<ProgramMachine<double>>(route.source->program(), route.source->zero_order());
            if (!counting->valid())
                counting.reset();
        }
        if (!program_seeds(seeds.data(), *route.source, result.degree, route.preparation.root_bound.value(),
                           static_cast<slong>(initial_precision), o, result.statistics_v3, &policy,
                           &route.preliminary.seeds, counting.get()))
            return false; // cancelled
    }
    route.preliminary.precision = static_cast<slong>(initial_precision);
    if (route.preliminary.seeds.raises) {
        result.statistics_v2.precision_increases += route.preliminary.seeds.raises;
        observe_precision(route.preliminary.seeds.precision);
    }
    route.preliminary.points.resize(n);
    for (size_t i = 0; i < n; ++i)
        acb_set(route.preliminary.points[i].value(), seeds.data() + static_cast<slong>(i));
    return true;
}

// Any derivative ball without correct bits among the points (the unsettled ones only, when the
// mask is set) at the evaluation the packet makes at `bits` (16 guard bits, within the cap): the
// evaluation cancels beyond that precision there.
bool ProgramRouteState::Controller::derivatives_unreliable(bool unsettled_only, slong bits) const {
    const auto n = static_cast<size_t>(result.degree);
    const auto cap = static_cast<slong>(o.max_work_bits);
    const auto &preliminary = route.preliminary;
    const slong evaluation_bits = std::min<slong>(resources::precision_sum(bits, 16), cap);
    ComplexBall value, derivative;
    for (size_t i = 0; i < n; ++i) {
        if (unsettled_only && preliminary.settled.size() == n && preliminary.settled[i])
            continue;
        if (!route.source->evaluate(value.value(), derivative.value(), preliminary.points[i].value(),
                                    evaluation_bits) ||
            acb_contains_zero(derivative.value()))
            return true;
    }
    return false;
}

// The raise: doubled, or straight to the precision the root counts needed when that is more.
void ProgramRouteState::Controller::raise_precision(const char *phase) {
    const auto n = static_cast<size_t>(result.degree);
    const auto cap = static_cast<slong>(o.max_work_bits);
    auto &preliminary = route.preliminary;
    preliminary.precision = std::min<slong>(
        cap, std::max(resources::precision_product(preliminary.precision, 2), preliminary.seeds.precision));
    ++preliminary.raises;
    ++result.statistics_v2.precision_increases;
    preliminary.previous_window = std::numeric_limits<double>::infinity();
    observe_precision(preliminary.precision);
    if (policy.trace) {
        TraceEvent event{
            "program_preliminary_raise", phase, static_cast<uint64_t>(preliminary.precision), n, &result.statistics_v3,
            &result.statistics_v2};
        event.reason = "derivative-without-correct-bits";
        policy.trace(policy.trace_data, event);
    }
}

// Raised level by level until the derivative balls at the unsettled points have correct bits (or
// the cap): a packet steps on midpoints, and a midpoint step from a derivative without correct
// bits is noise that can throw a point to an enormous radius, from which it crawls back by
// (n - 1) / (n + 1) a sweep and halves its movement every window. Returns whether the balls read
// at the end. History: docs/code-history-20261021.md, solver.cpp 19.
bool ProgramRouteState::Controller::raise_until_readable(const char *phase) {
    const auto cap = static_cast<slong>(o.max_work_bits);
    auto &preliminary = route.preliminary;
    while (derivatives_unreliable(true, preliminary.precision)) {
        if (preliminary.precision >= cap)
            return false;
        raise_precision(phase);
    }
    return true;
}

// The binary64 stage first (program_machine.hpp): one quantum of sweeps a work unit, the stage
// kept across work units; then the handover, and the ball packets take the floor and held points.
// The binary64 stage needs a cap of at least binary64's 53 bits: under a lower cap it does not run.
// With a root at zero (every odd Chebyshev, Hermite and Legendre program has one) it iterates q's
// roots, as the source does (the evaluator's zero order). History:
// docs/code-history-20261021.md, solver.cpp 20.
std::optional<ProgramRouteState::Controller::Step> ProgramRouteState::Controller::machine_step() {
    const auto n = static_cast<size_t>(result.degree);
    auto &preliminary = route.preliminary;
    if (heuristics::program_machine_stage && !preliminary.machine_done &&
        o.max_work_bits >= static_cast<uint64_t>(machine::binary64_bits)) {
        if (!preliminary.machine) {
            {
                RootVector seeds(result.degree);
                for (size_t i = 0; i < n; ++i)
                    acb_set(seeds.data() + static_cast<slong>(i), preliminary.points[i].value());
                preliminary.machine = std::make_unique<ProgramMachineStage>(
                    route.source->program(), seeds.data(), n, heuristics::program_machine_chunk, o.temporary_bytes,
                    nullptr, route.source->zero_order());
                if (!preliminary.machine->valid()) {
                    preliminary.machine.reset();
                    preliminary.machine_refused = true;
                }
            }
            if (preliminary.machine_refused)
                preliminary.machine_done = true;
        }
        if (preliminary.machine) {
            const unsigned machine_sweeps =
                o.packet_sweeps(heuristics::program_machine_packet_sweeps, result.statistics_v3);
            if (!machine_sweeps) {
                return Step::IterationLimit;
            }
            TraceScope trace(&policy, "program_machine", static_cast<uint64_t>(machine::binary64_bits), n,
                             &result.statistics_v3, &result.statistics_v2);
            const auto step = preliminary.machine->advance(machine_sweeps, o, result.statistics_v3);
            preliminary.machine_sweeps += step.sweeps;
            preliminary.scaled_evaluations = step.scaled_evaluations;
            preliminary.machine_converged = step.converged;
            preliminary.machine_floor = step.floor;
            preliminary.machine_held = step.held;
            preliminary.machine_stalled = step.stalled;
            if (step.cancelled)
                return Step::Stop;
            if (!step.finished)
                return Step::NextUnit;
            // The handover: the points, the settled mask, a fresh window, and the read of the
            // derivative balls at the floor and held points.
            RootVector exported(result.degree);
            preliminary.machine->export_points(exported.data());
            for (size_t i = 0; i < n; ++i)
                if (acb_is_finite(exported.data() + static_cast<slong>(i)))
                    acb_set(preliminary.points[i].value(), exported.data() + static_cast<slong>(i));
            const std::vector<ProgramMachineStage::Stop> stops(preliminary.machine->stops().begin(),
                                                                       preliminary.machine->stops().end());
            preliminary.settled.assign(n, 0);
            for (size_t i = 0; i < n; ++i)
                preliminary.settled[i] = stops[i] == ProgramMachineStage::Stop::Converged;
            preliminary.machine.reset();
            preliminary.machine_done = true;
            preliminary.previous_window = std::numeric_limits<double>::infinity();
            preliminary.current_window = -std::numeric_limits<double>::infinity();
            preliminary.window_sweeps = 0;
            if (policy.trace) {
                TraceEvent event{
                    "program_machine",     "handover",           static_cast<uint64_t>(machine::binary64_bits), n,
                    &result.statistics_v3, &result.statistics_v2};
                event.model_settled = step.converged;
                event.reason = step.stalled ? "stalled" : "finished";
                policy.trace(policy.trace_data, event);
            }
            if (step.converged == n) {
                // Every point converged in binary64 (to about 2^-45): one short ball packet with every
                // point active polishes them for the refinement (cubic convergence takes the step from
                // 2^-45 to below the packet's 2^-(prec-12) in a sweep or two), which spares it a second
                // regeneration from seeds of binary64 accuracy.
                preliminary.settled.clear();
                preliminary.polish = true;
                return Step::NextUnit;
            }
            raise_until_readable("binary64-floor");
            return Step::NextUnit;
        }
    }
    return std::nullopt; // disabled, already handed over, or refused: a ball packet in this unit
}

ProgramRouteState::Controller::Step ProgramRouteState::Controller::ball_step() {
    const auto n = static_cast<size_t>(result.degree);
    auto &preliminary = route.preliminary;
    const slong precision = preliminary.precision;
    const uint64_t window = (n + 2) / 2; // ceil((n + 1) / 2)
    constexpr double halving_bits = 1;   // a window's movement must fall by a factor two
    const auto cap = static_cast<slong>(o.max_work_bits);
    const unsigned sweeps = o.packet_sweeps(preliminary.polish ? heuristics::program_machine_polish_sweeps
                                                               : heuristics::secular_packet_sweeps,
                                            result.statistics_v3);
    preliminary.polish = false;
    if (!sweeps) {
        return Step::IterationLimit;
    }
    TraceScope trace(&policy, "program_preliminary", initial_precision, n, &result.statistics_v3,
                     &result.statistics_v2);
    std::vector<ComplexBall> no_nodes(n), no_residues(n); // the direct packet reads no model
    // The points the binary64 stage settled are inactive: in the sums as fixed roots, not moved.
    std::vector<uint8_t> active = preliminary.settled.size() == n ? preliminary.settled : std::vector<uint8_t>(n, 0);
    std::vector<slong> bits(n, precision);
    RootVector initial(result.degree), output(result.degree);
    for (size_t i = 0; i < n; ++i)
        acb_set(initial.data() + static_cast<slong>(i), preliminary.points[i].value());
    const uint64_t before = result.statistics_v3.candidate_sweeps;
    const auto progress =
        secular_candidate_packet(output.data(), initial.data(), no_nodes, no_residues, active, bits, sweeps, o,
                                 result.statistics_v3, route.source.get(), true, true, nullptr, 0, &policy);
    if (progress.cancelled)
        return Step::Stop;
    const uint64_t swept = result.statistics_v3.candidate_sweeps - before;
    preliminary.sweeps += swept;
    if (!progress.finite) {
        preliminary.done = true; // the points before this packet start the refinement
        return Step::NextUnit;
    }
    double movement = -std::numeric_limits<double>::infinity();
    ComplexBall difference;
    Magnitude size, resolution;
    for (size_t i = 0; i < n; ++i) {
        acb_sub(difference.value(), output.data() + static_cast<slong>(i), preliminary.points[i].value(), precision);
        acb_get_mag(size.value(), difference.value());
        // A point that did not move counts as -inf, as the window rule below expects. FLINT's log2 of a
        // zero magnitude is the finite COEFF_MIN (about -2^62): once every point stands still the test
        // would compare -2^62 with -2^62 - 1, one double, so no window would be a stall again and the
        // phase would sweep until the work budget. So does a movement below the point's resolution at
        // this precision, |z| 2^-precision: a step the phase cannot use. A coordinate converging to
        // zero, a real root's imaginary part beside held points of a multiple root, would otherwise
        // move less every window without a floor and be read as progress. History:
        // docs/code-history-20261021.md, solver.cpp 21.
        acb_get_mag(resolution.value(), output.data() + static_cast<slong>(i));
        mag_mul_2exp_si(resolution.value(), resolution.value(), -precision);
        if (!mag_is_zero(size.value()) && mag_cmp(size.value(), resolution.value()) > 0)
            movement = std::max(movement, mag_get_d_log2_approx(size.value()));
        acb_set(preliminary.points[i].value(), output.data() + static_cast<slong>(i));
    }
    if (progress.converged) {
        preliminary.done = true;
        return Step::NextUnit;
    }
    preliminary.current_window = std::max(preliminary.current_window, movement);
    preliminary.window_sweeps += swept;
    if (preliminary.window_sweeps >= window) {
        // The first window at this precision (previous_window is infinite until a window ends).
        const bool first = !std::isfinite(preliminary.previous_window);
        // A window in which no point moved (every active point held: an evaluation that is not
        // finite at this precision, or none active) is a stall as well: -inf never exceeds a bound.
        const bool moved = preliminary.current_window > -std::numeric_limits<double>::infinity();
        const bool stalled = !moved || preliminary.current_window > preliminary.previous_window - halving_bits;
        if (policy.trace) {
            TraceEvent event{
                "program_preliminary_window", stalled ? "stall" : "progress", static_cast<uint64_t>(precision), n,
                &result.statistics_v3,        &result.statistics_v2};
            policy.trace(policy.trace_data, event);
        }
        preliminary.previous_window = std::min(preliminary.previous_window, preliminary.current_window);
        preliminary.current_window = -std::numeric_limits<double>::infinity();
        preliminary.window_sweeps = 0;
        const bool foretold = first && precision < preliminary.seeds.precision;
        if (stalled || foretold) {
            const bool unreliable = derivatives_unreliable(true, precision);
            if (unreliable && precision < cap)
                raise_until_readable(stalled ? "stall" : "first-window");
            else if (stalled)
                preliminary.done = true;
        }
    }
    return Step::NextUnit;
}

// The preliminary phase: packets of Aberth sweeps on the program's own evaluations (the direct
// packet, the interaction tree from heuristics::interaction_tree_degree) at the initial precision,
// from the seeds, until the points converge or stop making progress. SecularRefinement::preliminary
// does this for coefficient sources only. Without it the seeds, on a circle at the root bound
// (about the degree for a Mandelbrot polynomial whose roots lie within 2, section 7), reach the
// model far from any root: from a symmetric circle Aberth's step is z (n - 1) / (n + 1), a crawl
// of about (n / 2) ln(R / r) sweeps whose relative correction is constant, which the refinement
// reads as stagnation: it raises the precision every fourth round and, after each raise, restarts
// the one unresolved component on its own circle, the big one, and the crawl begins again
// until the cap. History: docs/code-history-20261021.md, solver.cpp 16.
// Progress is measured over windows of sweeps: the largest movement of a point over a packet,
// taken over the window, must halve against the smallest window so far (not the last one: noisy
// jumps halve a last window by chance, which delays the stall). The crawl contracts the points by
// ((n - 1) / (n + 1))^m over m sweeps, e^-1 at m = (n + 1) / 2 (the window), so a crawl always
// passes. A window that does not halve is a stall: the points are at a floor this precision sets.
// History: docs/code-history-20261021.md, solver.cpp 17.
// Two floors look alike in the movement and differ in the balls. A cluster the precision cannot
// resolve (a double root: its pair settles 2^-(prec/2) from the root) leaves the derivative with
// correct bits, since p' vanishes only at the root itself; the phase ends there and the
// refinement's model, growth and certificates take over, as they would have from the seeds.
// A point whose evaluation cancels beyond the precision (a Laguerre polynomial at its large roots)
// has a derivative ball containing zero, its Aberth step is noise and it jumps about the
// configuration; handed over like that, the refinement would certify nothing and its growth would
// re-seed every node on the big circle at each raise. So a stall with any derivative without
// correct bits raises the phase's precision and goes on, up to the cap: doubled, or straight to
// the precision the root counts needed (program_seeds) when that is more. While the phase is
// below that precision the balls are read at its first window's end as well, since noise at some
// points is then likely and the stall rule, which compares windows, can pass two windows of it
// before it says so; where the counts' raises were the balls' pessimism (Chebyshev at |z| = 1)
// the phase converges before its first window ends and reads nothing. History:
// docs/code-history-20261021.md, solver.cpp 18.
ProgramRouteState::Controller::Step ProgramRouteState::Controller::preliminary_step() {
    if (auto step = machine_step())
        return *step;
    return ball_step();
}

// The separated-roots pass (phases 2a and 2b): one round a work unit, from the preliminary phase's points,
// until it certifies (its cover is the answer, with no regeneration) or fails (the refinement answers,
// from the preliminary phase's points).
ProgramRouteState::Controller::DiskStep ProgramRouteState::Controller::disk_step() {
    const auto n = static_cast<size_t>(result.degree);
    const slong precision = route.preliminary.precision;
    if (heuristics::program_disk_pass) {
        auto &disks = route.disks;
        if (!disks.pass) {
            disks.source =
                std::make_unique<ProgramSource>(solved(), result.degree, route.preparation.leading, route.zero_order);
            if (const auto order = static_cast<int64_t>(route.zero_order))
                disks.undivided =
                    std::make_unique<ProgramSource>(solved(), result.degree + order, route.preparation.leading);
            RootVector start(result.degree);
            for (size_t i = 0; i < n; ++i)
                acb_set(start.data() + static_cast<slong>(i), route.preliminary.points[i].value());
            disks.pass = std::make_unique<SeparatedRootsPass>(start.data(), n, precision);
        }
        if (disks.pass->state() == SeparatedRootsPass::State::Running) {
            TraceScope trace(&policy, "program_disks", static_cast<uint64_t>(disks.pass->precision()), n,
                             &result.statistics_v3, &result.statistics_v2);
            ++result.statistics_v3.validation_calls;
            const auto state = disks.pass->advance(*disks.source, o, result.statistics_v3, disks.undivided.get());
            const auto bits = static_cast<uint64_t>(disks.pass->precision());
            observe_precision(disks.pass->precision());
            if (state == SeparatedRootsPass::State::Certified) {
                result.roots = disks.pass->cover();
                result.complete = true;
                result.work_bits = std::max(result.work_bits, bits);
                result.goal_met = program_goal(result.roots, std::max<uint64_t>(64, bits));
                return DiskStep::NextUnit;
            }
            if (state == SeparatedRootsPass::State::Running)
                return DiskStep::NextUnit; // the next round in the next work unit
        }
    }
    return DiskStep::RefineNow;
}

cps_stop_reason ProgramRouteState::Controller::refinement_step() {
    const auto n = static_cast<size_t>(result.degree);
    const slong precision = route.preliminary.precision;
    if (!route.refinement) {
        RootVector seeds(result.degree);
        for (size_t i = 0; i < n; ++i)
            acb_set(seeds.data() + static_cast<slong>(i), route.preliminary.points[i].value());
        // From the preliminary phase's points, not the failed pass's: a pass that fails may leave two points
        // at one root, which the refinement then spends its restarts undoing. So a failed pass costs its own
        // rounds and nothing more. History: docs/code-history-20261021.md, solver.cpp 22.
        route.refinement =
            std::make_unique<SecularRefinement>(std::move(route.source), seeds.data(), 1, precision, policy,
                                               solve_local_model);
    }
    // Under the accuracy goal, with a root at zero, the refinement works to the goal without its absolute
    // tolerance. A region of q that holds zero, published with the root at zero (with_zero), holds two
    // distinct roots and never meets the accuracy goal, however narrow, while q's own goal could accept it
    // by the tolerance and leave the refinement nothing to do: the deadlock that judging q's covers with
    // the root at zero (program_goal) avoids under the cluster goal. The relative test excludes zero (q(0)
    // is not zero), and the loop still stops at the caller's goal, judged at each update. History:
    // docs/code-history-20261021.md, solver.cpp 23.
    cps_options refinement_options = o;
    if (route.zero_order && o.goal == CPS_ROOT_ACCURACY)
        refinement_options.absolute_tolerance = Number();
    const auto reason = route.refinement->step(refinement_options, result.statistics_v2, result.statistics_v3);
    update_cover();
    return reason;
}

void ProgramRouteState::Controller::update_cover() {
    if (!route.refinement)
        return;
    result.work_bits = route.refinement->precision();
    result.progress.max_outer_bits = std::max(result.progress.max_outer_bits, result.work_bits);
    std::vector<Root> cover;
    if (route.refinement->cover(cover)) {
        const auto prec = std::max<uint64_t>(64, result.work_bits);
        if (program_goal(cover, prec) || !result.complete ||
            prefer_cover(cover, result.roots, static_cast<slong>(prec)))
            result.roots = std::move(cover);
        result.complete = true;
        result.goal_met = program_goal(result.roots, prec);
    }
}

ProgramRouteState::Controller::FinalStatus ProgramRouteState::Controller::iterate() {
    if (result.complete)
        result.goal_met = program_goal(result.roots, std::max<uint64_t>(64, result.work_bits));
    bool iteration_stop = false, precision_stop = false;
    while (!result.goal_met && !o.cancelled() && result.progress.work_units_this_call < o.work_units) {
        if (!o.sweeps_left(result.statistics_v3)) {
            iteration_stop = true;
            break;
        }
        ++result.statistics_v2.work_units;
        ++result.progress.work_units_this_call;
        if (!route.refinement && route.preliminary.points.empty() && !seed_source())
            break;
        if (!route.preliminary.done) {
            const auto step = preliminary_step();
            if (step == Step::IterationLimit)
                iteration_stop = true;
            if (step != Step::NextUnit)
                break;
            continue;
        }
        if (disk_step() == DiskStep::NextUnit)
            continue;
        const auto reason = refinement_step();
        if (reason == CPS_STOP_CANCELLED)
            break;
        if (reason == CPS_STOP_ITERATION_BUDGET) {
            iteration_stop = !result.goal_met;
            break;
        }
        if (reason == CPS_STOP_PRECISION && !result.goal_met) {
            precision_stop = true;
            break;
        }
        if (reason == CPS_STOP_INPUT)
            // Only an inexact source's residues are input-limited (regenerate_residues); a program is exact.
            throw std::logic_error("the program route reported an input-limited stop on exact input");
    }
    // A met goal is a success whatever the last checkpoint polls (I2; C25).
    const bool cancelled = !result.goal_met && o.cancelled();
    return FinalStatus(cancelled         ? CPS_STOP_CANCELLED
                       : result.goal_met ? CPS_STOP_NONE
                       : precision_stop  ? CPS_STOP_PRECISION
                       : iteration_stop  ? CPS_STOP_ITERATION_BUDGET
                                         : CPS_STOP_WORK_BUDGET,
                       cancelled         ? CPS_CANCELLED
                       : result.goal_met ? CPS_SUCCESS
                       : precision_stop  ? CPS_PRECISION_LIMIT
                                         : CPS_WORK_LIMIT,
                       precision_stop);
}

} // namespace certpolysolve::detail
