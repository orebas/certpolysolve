// Stage 3b-1's C additions (docs/api-spec-0.5.md, the order of work): the 0.5 statistics, engine counters,
// components, stop reason and preparation flag hold what the engines counted; C14's facts on known answers; the size
// protocol; a getter for every option, the options' copy and the split setters; the resource limits; and
// certpolysolve_testing.h's poll counter, which must stop a solve exactly where a callback that counts its calls stops
// it. Through the public C API, as an application would use it, except for three internal reads and writes
// (internal.hpp; no internal code runs): the published blocks are compared field by field with the result's internal
// counters, so a field the C API maps wrongly shows; the budget the removed combined setter set is set on the
// options' fields; and the counting callback, which no public call sets since stage 3b-3, is too.
#include <certpolysolve/certpolysolve.h>
#include <certpolysolve/certpolysolve_testing.h>
#include "internal.hpp"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
void require(bool condition, const std::string &message) {
    if (!condition) throw std::runtime_error(message);
}
void ok(cps_error e, const std::string &what) { require(e == CPS_OK, what + ": " + cps_last_error()); }

struct Polynomial {
    cps_polynomial *p = nullptr;
    explicit Polynomial(const std::vector<const char *> &ascending) {
        ok(cps_polynomial_create(&p), "create");
        for (size_t k = 0; k < ascending.size(); ++k) ok(cps_polynomial_set_coefficient(p, k, ascending[k], "0"), "set");
    }
    ~Polynomial() { cps_polynomial_destroy(p); }
};
struct Options {
    cps_options *o = nullptr;
    Options() { ok(cps_options_create(&o), "options"); }
    ~Options() { cps_options_destroy(o); }
};
struct Result {
    cps_result *r = nullptr;
    ~Result() { cps_result_destroy(r); }
};

// Every field the result reports but its times: what two solves on one path must agree in.
std::string describe(const cps_result *r) {
    std::string s = "termination " + std::to_string(cps_result_termination(r)) + " stop " +
                    std::to_string(cps_result_stop_reason(r)) + " complete " + std::to_string(cps_result_complete(r)) +
                    " met " + std::to_string(cps_result_goal_met(r)) + " bits " + std::to_string(cps_result_work_bits(r)) +
                    "\n";
    for (size_t i = 0; i < cps_result_size(r); ++i) {
        cps_root_view v{};
        ok(cps_result_get_root(r, i, &v), "root");
        for (const auto &e : {v.real_lower, v.real_upper, v.imag_lower, v.imag_upper})
            s += std::string(e.mantissa) + "@" + e.exponent + " ";
        s += std::to_string(v.count) + " " + std::to_string(v.multiplicity) + " " + std::to_string(v.realness) + " " +
             v.certificate + "\n";
    }
    cps_statistics t{};
    t.size = sizeof t;
    ok(cps_result_get_statistics(r, &t), "statistics");
    cps_engine_counters c{};
    c.size = sizeof c;
    ok(cps_result_get_engine_counters(r, &c), "engine counters");
    s += "work " + std::to_string(t.work_units) + " sweeps " + std::to_string(t.sweeps) + " evaluations " +
         std::to_string(c.evaluations) + " interactions " + std::to_string(c.root_interactions) + "\n";
    return s + cps_result_diagnostics(r);
}

// The 0.5 blocks against the engines' counters, field by field, on one result. The answering engine's stage
// counters are its own; after a hand-off (piecewise to secular) the result carries the earlier engine's apart, and the
// 0.5 totals add both: exactly their sum, the preparation time then strictly more than the answering engine's (both
// engines prepared the input), and nothing carried without a hand-off.
void blocks_agree(const cps_result *r, const std::string &name, bool handed_off = false) {
    const auto &a = r->statistics;
    const auto &k = r->carried_v1;
    const auto &b = r->statistics_v2;
    const auto &c = r->statistics_v3;
    const auto &p = r->progress;
    cps_statistics s{};
    s.size = sizeof s;
    ok(cps_result_get_statistics(r, &s), "statistics");
    cps_engine_counters e{};
    e.size = sizeof e;
    ok(cps_result_get_engine_counters(r, &e), "engine counters");
    uint64_t roots = 0, resolved = 0;
    for (size_t i = 0; i < cps_result_size(r); ++i) {
        cps_root_view v{};
        ok(cps_result_get_root(r, i, &v), "root");
        roots += v.count;
        resolved += v.multiplicity != 0;
    }
    const auto same = [&](uint64_t x, uint64_t y, const char *field) {
        require(x == y, name + ": " + field + " " + std::to_string(x) + " is not the engine's " + std::to_string(y));
    };
    // A stage counter: the answering engine's plus the earlier engines', and those only after a hand-off.
    const auto carried = [&](uint64_t x, uint64_t own, uint64_t earlier, const char *field) {
        same(x, own + earlier, field);
        require(handed_off || earlier == 0, name + ": " + field + " carried without a hand-off");
    };
    require(handed_off == (std::string(cps_result_diagnostics(r)).find(" to secular") != std::string::npos),
            name + ": the hand-off is not as expected");
    same(s.size, sizeof s, "size");
    same(s.work_units, b.work_units, "work_units");
    same(s.sweeps, c.candidate_sweeps, "sweeps");
    same(s.precision_increases, b.precision_increases, "precision_increases");
    same(s.final_work_bits, cps_result_work_bits(r), "final_work_bits");
    require(s.highest_work_bits >= s.final_work_bits && s.highest_work_bits >= p.max_outer_bits,
            name + ": highest_work_bits below a precision the solve reached");
    same(s.regions, cps_result_size(r), "regions");
    same(s.roots, roots, "roots");
    same(s.resolved_regions, resolved, "resolved_regions");
    same(s.exact_factors, a.exact_factors, "exact_factors");
    carried(s.preparation_ns, a.preprocessing_ns, k.preprocessing_ns, "preparation_ns");
    if (handed_off)
        require(s.preparation_ns > a.preprocessing_ns, name + ": the preparation before the hand-off is not counted");
    require(s.total_ns > 0 && s.total_ns >= s.preparation_ns, name + ": total_ns below a stage of the call");
    carried(e.input_enclosure_ns, a.input_enclosure_ns, k.input_enclosure_ns, "input_enclosure_ns");
    carried(e.isolation_ns, a.isolation_ns, k.isolation_ns, "isolation_ns");
    carried(e.coverage_validation_ns, a.coverage_validation_ns, k.coverage_validation_ns, "coverage_validation_ns");
    same(e.backend_ns, b.backend_ns, "backend_ns");
    same(e.taylor_shift_ns, b.taylor_shift_ns, "taylor_shift_ns");
    same(e.counting_ns, b.counting_ns, "counting_ns");
    same(e.refinement_ns, b.refinement_ns, "refinement_ns");
    same(e.search_ns, b.search_ns, "search_ns");
    same(e.candidate_ns, c.candidate_ns, "candidate_ns");
    same(e.evaluation_ns, c.evaluation_ns, "evaluation_ns");
    same(e.interaction_ns, c.interaction_ns, "interaction_ns");
    same(e.validation_ns, c.validation_ns, "validation_ns");
    same(e.regeneration_ns, c.regeneration_ns, "regeneration_ns");
    carried(e.precision_attempts, a.precision_attempts, k.precision_attempts, "precision_attempts");
    carried(e.isolation_passes, a.isolation_passes, k.isolation_passes, "isolation_passes");
    same(e.disk_tests, b.disk_tests, "disk_tests");
    same(e.graeffe_steps, b.graeffe_steps, "graeffe_steps");
    same(e.taylor_shifts, b.taylor_shifts, "taylor_shifts");
    same(e.taylor_cache_hits, b.taylor_cache_hits, "taylor_cache_hits");
    same(e.boxes_split, b.boxes_split, "boxes_split");
    same(e.boxes_excluded, b.boxes_excluded, "boxes_excluded");
    same(e.newton_attempts, b.newton_attempts, "newton_attempts");
    same(e.newton_successes, b.newton_successes, "newton_successes");
    same(e.backend_calls, b.backend_calls, "backend_calls");
    same(e.evaluations, c.evaluations, "evaluations");
    same(e.evaluation_terms, c.evaluation_terms, "evaluation_terms");
    same(e.root_interactions, c.root_interactions, "root_interactions");
    same(e.validation_calls, c.validation_calls, "validation_calls");
    same(e.secular_regenerations, c.secular_regenerations, "secular_regenerations");
    same(e.double_seed_calls, c.double_seed_calls, "double_seed_calls");
    same(e.candidate_stagnations, c.candidate_stagnations, "candidate_stagnations");
    same(e.real_refinement_steps, c.real_refinement_steps, "real_refinement_steps");
    same(e.real_bisections, c.real_bisections, "real_bisections");
    same(e.pending_components, p.pending_components, "pending_components");
    same(e.live_cells, p.live_cells, "live_cells");
    same(e.min_pending_bits, p.min_pending_bits, "min_pending_bits");
    same(e.max_pending_bits, p.max_pending_bits, "max_pending_bits");
    same(e.input_uncertainty, p.input_uncertainty, "input_uncertainty");
    same(static_cast<uint64_t>(cps_result_stop_reason(r)), static_cast<uint64_t>(p.stop_reason), "stop reason");
    same(static_cast<uint64_t>(cps_result_preparation_limited(r)), p.preprocessing_limited != 0, "preparation flag");
    same(cps_result_component_count(r), r->components.size(), "component count");
    for (size_t i = 0; i < cps_result_component_count(r); ++i) {
        const auto &old = r->components[i];
        cps_component_progress now{};
        ok(cps_result_get_component(r, i, &now), "component");
        require(old.factor_index == now.factor_index && old.work_bits == now.work_bits &&
                    old.live_cells == now.live_cells && old.certified_count == now.certified_count &&
                    old.state == now.state,
                name + ": a component differs from the engine's");
    }
    cps_component_progress past{};
    require(cps_result_get_component(r, cps_result_component_count(r), &past) == CPS_OUT_OF_RANGE,
            name + ": a component index past the count was accepted");
}

void statistics_protocol(const cps_result *r) {
    // Too small: refused, untouched.
    cps_statistics small{};
    small.size = sizeof small - 1;
    small.work_units = 12345;
    require(cps_result_get_statistics(r, &small) == CPS_INVALID_ARGUMENT && small.work_units == 12345 &&
                small.size == sizeof small - 1,
            "a size below the layout was accepted or changed the output");
    cps_engine_counters few{};
    few.size = 0;
    require(cps_result_get_engine_counters(r, &few) == CPS_INVALID_ARGUMENT, "engine counters: size 0 accepted");
    // Larger: the layout filled, the bytes beyond it untouched, the filled size reported.
    struct Bigger {
        cps_statistics s;
        unsigned char tail[32];
    } big;
    std::memset(&big, 0xA5, sizeof big);
    big.s.size = sizeof big;
    ok(cps_result_get_statistics(r, &big.s), "a larger size");
    require(big.s.size == sizeof(cps_statistics), "the filled size is not the layout's");
    for (unsigned char byte : big.tail) require(byte == 0xA5, "bytes beyond the layout were written");
    require(cps_result_get_statistics(nullptr, &big.s) == CPS_INVALID_ARGUMENT &&
                cps_result_get_statistics(r, nullptr) == CPS_INVALID_ARGUMENT &&
                cps_result_get_engine_counters(nullptr, &few) == CPS_INVALID_ARGUMENT,
            "NULL statistics arguments accepted");
    require(cps_result_stop_reason(nullptr) == CPS_STOP_INPUT && cps_result_preparation_limited(nullptr) == 0 &&
                cps_testing_result_polls(nullptr) == 0,
            "the NULL fallbacks");
}

void option_getters() {
    Options fresh;
    require(cps_options_method(fresh.o) == cps_options_method(nullptr) && cps_options_method(nullptr) == CPS_METHOD_AUTO,
            "method default");
    require(cps_options_precision(fresh.o) == CPS_DEFAULT_BITS && cps_options_precision(nullptr) == CPS_DEFAULT_BITS,
            "precision default");
    require(cps_options_max_work_bits(fresh.o) == CPS_DEFAULT_MAX_WORK_BITS &&
                cps_options_max_work_bits(nullptr) == CPS_DEFAULT_MAX_WORK_BITS,
            "cap default");
    require(cps_options_goal(fresh.o) == CPS_GOAL_DEFAULT && cps_options_goal(nullptr) == CPS_GOAL_DEFAULT,
            "goal default");
    require(std::strcmp(cps_options_cluster_radius(fresh.o), "0") == 0 &&
                std::strcmp(cps_options_cluster_radius(nullptr), "0") == 0 &&
                std::strcmp(cps_options_absolute_tolerance(fresh.o), "0") == 0 &&
                std::strcmp(cps_options_absolute_tolerance(nullptr), "0") == 0,
            "string defaults");
    require(cps_options_real_only(fresh.o) == 0 && cps_options_work_budget(fresh.o) == 0 &&
                cps_options_work_budget(nullptr) == 0 && cps_options_live_cells(fresh.o) == CPS_DEFAULT_LIVE_CELLS &&
                cps_options_max_iterations(fresh.o) == 0 && cps_options_time_limit(fresh.o) == 0,
            "the other defaults");
    const char *ends[4] = {"x", "x", "x", "x"};
    ok(cps_options_get_region(fresh.o, &ends[0], &ends[1], &ends[2], &ends[3]), "region getter");
    for (const char *end : ends) require(end == nullptr, "an unset region side is not NULL");
    require(cps_options_get_region(nullptr, &ends[0], &ends[1], &ends[2], &ends[3]) == CPS_INVALID_ARGUMENT &&
                cps_options_get_region(fresh.o, nullptr, &ends[1], &ends[2], &ends[3]) == CPS_INVALID_ARGUMENT,
            "region getter NULLs");

    Options o;
    ok(cps_options_set_method(o.o, CPS_METHOD_SECULAR), "method");
    ok(cps_options_set_precision(o.o, 100), "precision");
    ok(cps_options_set_max_work_bits(o.o, 4096), "cap");
    ok(cps_options_set_absolute_tolerance(o.o, "  1e-9 "), "tolerance");
    ok(cps_options_set_cluster_radius(o.o, " 1/1024"), "radius");
    ok(cps_options_set_real_only(o.o, 7), "real only");
    ok(cps_options_set_max_iterations(o.o, 99), "iterations");
    ok(cps_options_set_time_limit(o.o, 2.5), "time limit");
    ok(cps_options_set_region(o.o, " -1/3", "2", "", nullptr), "region");
    require(cps_options_method(o.o) == CPS_METHOD_SECULAR && cps_options_precision(o.o) == 100 &&
                cps_options_max_work_bits(o.o) == 4096 && cps_options_goal(o.o) == CPS_CLUSTER_RADIUS &&
                std::strcmp(cps_options_cluster_radius(o.o), "1/1024") == 0 &&
                std::strcmp(cps_options_absolute_tolerance(o.o), "1e-9") == 0 && cps_options_real_only(o.o) == 1 &&
                cps_options_max_iterations(o.o) == 99 && cps_options_time_limit(o.o) == 2.5,
            "a getter does not return what was set");
    ok(cps_options_get_region(o.o, &ends[0], &ends[1], &ends[2], &ends[3]), "region getter");
    require(ends[0] && std::strcmp(ends[0], "-1/3") == 0 && ends[1] && std::strcmp(ends[1], "2") == 0 && !ends[2] &&
                !ends[3],
            "the region getter");
    // Another goal clears the radius; the radius 0 clears it too.
    ok(cps_options_set_goal(o.o, CPS_ROOT_APPROXIMATION), "goal");
    require(cps_options_goal(o.o) == CPS_ROOT_APPROXIMATION && std::strcmp(cps_options_cluster_radius(o.o), "0") == 0,
            "a goal left the radius");
    ok(cps_options_set_cluster_radius(o.o, "1e-3"), "radius again");
    ok(cps_options_set_cluster_radius(o.o, "0"), "radius 0");
    require(std::strcmp(cps_options_cluster_radius(o.o), "0") == 0, "the radius 0 is not reported as 0");
    require(cps_options_goal(o.o) == CPS_GOAL_DEFAULT, "the radius 0 did not restore the default goal");
    ok(cps_options_set_goal(o.o, CPS_ROOT_ISOLATION), "isolation");
    ok(cps_options_set_goal(o.o, CPS_GOAL_DEFAULT), "the default goal set");
    require(cps_options_goal(o.o) == CPS_GOAL_DEFAULT, "the default goal is not reported");
    // Four unbounded sides clear the region.
    ok(cps_options_set_region(o.o, "", nullptr, "", ""), "clear the region");
    ok(cps_options_get_region(o.o, &ends[0], &ends[1], &ends[2], &ends[3]), "region getter");
    for (const char *end : ends) require(end == nullptr, "a cleared region side is not NULL");

    // The split setters.
    require(cps_options_set_max_work_bits(o.o, CPS_MIN_MAX_WORK_BITS - 1) == CPS_INVALID_ARGUMENT &&
                cps_options_max_work_bits(o.o) == 4096,
            "a cap below the minimum was accepted or changed the options");
    require(cps_options_set_max_work_bits(o.o, UINT64_MAX) == CPS_OUT_OF_RANGE, "an unrepresentable cap accepted");
    ok(cps_options_set_max_work_bits(o.o, CPS_MIN_MAX_WORK_BITS), "the minimum cap");
    require(cps_options_max_work_bits(o.o) == CPS_MIN_MAX_WORK_BITS && cps_options_precision(o.o) == 100,
            "the cap setter changed the precision");
    ok(cps_options_set_work_budget(o.o, 1234), "work budget");
    ok(cps_options_set_live_cells(o.o, 99), "live cells");
    require(cps_options_work_budget(o.o) == 1234 && cps_options_live_cells(o.o) == 99, "budget getters");
    ok(cps_options_set_work_budget(o.o, 0), "automatic work budget");
    ok(cps_options_set_live_cells(o.o, 0), "default live cells");
    require(cps_options_work_budget(o.o) == 0 && cps_options_live_cells(o.o) == CPS_DEFAULT_LIVE_CELLS,
            "0 does not select the default");

    // The resource limits (until stage 3b-3 also under _v1 names).
    cps_resource_limits limits;
    cps_resource_limits_default(&limits);
    require(limits.exact_materialization_bits == UINT64_MAX && limits.exact_operations == UINT64_MAX &&
                limits.temporary_bytes == UINT64_MAX && limits.flags == 0,
            "the default limits are not every allowance unspecified and no flag");
    cps_resource_limits_default(nullptr);
    limits.exact_operations = 77;
    limits.flags = CPS_RESOURCE_ALLOW_EXPONENT_EXPANSION;
    ok(cps_options_set_resource_limits(o.o, &limits), "limits");
    cps_resource_limits back{};
    ok(cps_options_get_resource_limits(o.o, &back), "limits back");
    require(std::memcmp(&limits, &back, sizeof limits) == 0, "the limits did not come back");
    limits.flags = 2;
    require(cps_options_set_resource_limits(o.o, &limits) == CPS_INVALID_ARGUMENT, "an unknown flag accepted");

    // The copy: every setting, independent of its source, and the same token.
    cps_options *copy = nullptr;
    ok(cps_options_copy(o.o, &copy), "copy");
    require(cps_options_method(copy) == cps_options_method(o.o) && cps_options_precision(copy) == cps_options_precision(o.o) &&
                cps_options_max_work_bits(copy) == cps_options_max_work_bits(o.o) &&
                cps_options_goal(copy) == cps_options_goal(o.o) &&
                std::strcmp(cps_options_absolute_tolerance(copy), cps_options_absolute_tolerance(o.o)) == 0 &&
                cps_options_real_only(copy) == cps_options_real_only(o.o) &&
                cps_options_max_iterations(copy) == cps_options_max_iterations(o.o) &&
                cps_options_time_limit(copy) == cps_options_time_limit(o.o),
            "the copy differs from its source");
    ok(cps_options_set_precision(copy, 64), "change the copy");
    ok(cps_options_set_max_work_bits(copy, 1024), "change the copy's cap");
    require(cps_options_precision(o.o) == 100, "changing the copy changed its source");
    cps_options_destroy(copy);
    require(cps_options_copy(nullptr, &copy) == CPS_INVALID_ARGUMENT && copy == nullptr &&
                cps_options_copy(o.o, nullptr) == CPS_INVALID_ARGUMENT,
            "copy NULLs");

    cps_cancel_token *token = nullptr;
    ok(cps_cancel_token_create(&token), "token");
    Options with_token;
    ok(cps_options_set_cancel_token(with_token.o, token), "set token");
    ok(cps_options_copy(with_token.o, &copy), "copy with a token");
    cps_cancel_token_request(token);
    Polynomial cubic({"-6", "11", "-6", "1"});
    Result stopped;
    ok(cps_solve(cubic.p, copy, &stopped.r), "solve with the copy");
    require(cps_result_termination(stopped.r) == CPS_CANCELLED && cps_result_stop_reason(stopped.r) == CPS_STOP_CANCELLED,
            "the copy does not name its source's token");
    cps_options_destroy(copy);
    cps_cancel_token_destroy(token);
}

// The split setters answer as the combined one did (the budget it set, on the options' fields): the same answer and
// counters.
void split_setters_answer_alike() {
    Polynomial p({"3", "-1", "0", "7", "0", "0", "2", "1"});
    Options combined, split, automatic, defaults;
    combined.o->work_units = 40;
    combined.o->work_units_default = false;
    combined.o->live_cells = 77;
    ok(cps_options_set_work_budget(split.o, 40), "work budget");
    ok(cps_options_set_live_cells(split.o, 77), "live cells");
    ok(cps_options_set_work_budget(automatic.o, 0), "automatic budget");
    Result a, b, c, d;
    ok(cps_solve(p.p, combined.o, &a.r), "solve");
    ok(cps_solve(p.p, split.o, &b.r), "solve");
    ok(cps_solve(p.p, automatic.o, &c.r), "solve");
    ok(cps_solve(p.p, defaults.o, &d.r), "solve");
    require(describe(a.r) == describe(b.r), "the split budget setters answer differently from the combined one");
    require(describe(c.r) == describe(d.r), "the automatic work budget answers differently from the default");
}

// The testing poll counter against a callback that counts its calls, on one solve and every stop of it: the same
// answer, stop and count at every k.
struct Countdown {
    uint64_t calls = 0, stop_at = 0;
};
int count_down(void *data) {
    auto &c = *static_cast<Countdown *>(data);
    ++c.calls;
    return c.stop_at && c.calls >= c.stop_at ? 1 : 0;
}
void polls_match_the_callback(const cps_polynomial *p, const std::string &name) {
    Options counting;
    ok(cps_testing_options_stop_at_poll(counting.o, UINT64_MAX), "count only");
    Result whole;
    ok(cps_solve(p, counting.o, &whole.r), "count");
    const uint64_t polls = cps_testing_result_polls(whole.r);
    Countdown all;
    Options with_callback;
    with_callback.o->cancel = count_down;
    with_callback.o->cancel_data = &all;
    Result by_callback;
    ok(cps_solve(p, with_callback.o, &by_callback.r), "callback count");
    require(polls > 0 && polls == all.calls, name + ": the counter's " + std::to_string(polls) +
                                                 " polls are not the callback's " + std::to_string(all.calls));
    require(describe(whole.r) == describe(by_callback.r), name + ": counting changed the answer");
    Options off;
    Result unarmed;
    ok(cps_solve(p, off.o, &unarmed.r), "off");
    require(cps_testing_result_polls(unarmed.r) == 0 && describe(unarmed.r) == describe(whole.r),
            name + ": an unarmed solve counted, or answered otherwise");
    // Every stop for the first polls, then a spread to the end.
    std::vector<uint64_t> stops;
    for (uint64_t k = 1; k <= std::min<uint64_t>(polls, 40); ++k) stops.push_back(k);
    for (uint64_t k = 41; k <= polls; k += std::max<uint64_t>(1, polls / 25)) stops.push_back(k);
    stops.push_back(polls);
    std::vector<std::string> answers;
    for (const uint64_t k : stops) {
        Options hook;
        ok(cps_testing_options_stop_at_poll(hook.o, k), "stop at");
        Countdown countdown{0, k};
        Options callback;
        callback.o->cancel = count_down;
        callback.o->cancel_data = &countdown;
        Result x, y;
        ok(cps_solve(p, hook.o, &x.r), "hook stop");
        ok(cps_solve(p, callback.o, &y.r), "callback stop");
        require(describe(x.r) == describe(y.r),
                name + ": the counter's stop at poll " + std::to_string(k) + " differs from the callback's:\n" +
                    describe(x.r) + "\n---\n" + describe(y.r));
        require(cps_testing_result_polls(x.r) == countdown.calls,
                name + ": at stop " + std::to_string(k) + " the polls differ from the callback's calls");
        require(cps_result_termination(x.r) == CPS_CANCELLED || k == polls,
                name + ": a stop before the last poll did not cancel");
        if (std::find(answers.begin(), answers.end(), describe(x.r)) == answers.end()) answers.push_back(describe(x.r));
    }
    // The stops land at different places: otherwise the comparison above would hold trivially.
    require(answers.size() > 2, name + ": the stops gave fewer than three different answers");
    std::cout << name << ": " << polls << " polls, " << stops.size() << " stops compared, " << answers.size()
              << " different answers\n";
}
// C14's facts on inputs whose answers are known: square-free verdicts on every route that prepares (and none where
// nothing is proved), the sufficient cap at a precision stop (the number the diagnostics state), the hand-off.
cps_result *solved(const cps_polynomial *p, int method = CPS_METHOD_AUTO, uint64_t cap = 0, bool real = false) {
    Options o;
    ok(cps_options_set_method(o.o, method), "method");
    if (cap) ok(cps_options_set_max_work_bits(o.o, cap), "cap");
    if (real) ok(cps_options_set_real_only(o.o, 1), "real only");
    cps_result *r = nullptr;
    ok(cps_solve(p, o.o, &r), "solve");
    return r;
}
void verdict(const cps_polynomial *p, cps_squarefree expected, const std::string &name, int method = CPS_METHOD_AUTO,
             bool real = false) {
    Result r;
    r.r = solved(p, method, 0, real);
    require(cps_result_squarefree(r.r) == expected,
            name + ": the square-free verdict is " + std::to_string(cps_result_squarefree(r.r)) + ", expected " +
                std::to_string(expected) + " (" + cps_result_diagnostics(r.r) + ")");
    // Forced piecewise leaves an input the selector would not send to it once it has spent a global solve's minimum
    // budget (its regime exit, since 21 October 2026): on these
    // small inputs that is a hand-off to secular, which answers, not a stop. Every other route answers alone.
    const bool regime_exit = method == CPS_METHOD_PIECEWISE && cps_result_handoff(r.r) == CPS_HANDOFF_REGIME_EXIT;
    require(cps_result_sufficient_work_bits(r.r) == 0 && (cps_result_handoff(r.r) == CPS_HANDOFF_NONE || regime_exit),
            name + ": a cap or a hand-off without a stop");
}
void facts() {
    Polynomial quadruple({"-1", "4", "-6", "4", "-1"}), trinomial({"1", "-3", "0", "0", "0", "0", "0", "0", "0", "0", "0", "1"}),
        composed({"-3", "0", "0", "0", "1"}), zero_double({"0", "0", "-1", "0", "0", "0", "1"}), double_root({"-2", "5", "-4", "1"}),
        zero_twice({"0", "0", "-3", "1"}), zero_once({"0", "-3", "1"});
    for (int method : {CPS_METHOD_AUTO, CPS_METHOD_CLASSICAL, CPS_METHOD_SECULAR, CPS_METHOD_PIECEWISE}) {
        const std::string m = " (method " + std::to_string(method) + ")";
        verdict(trinomial.p, CPS_SQUAREFREE_PROVED, "x^11 - 3x + 1" + m, method);
        verdict(double_root.p, CPS_SQUAREFREE_DISPROVED, "(x - 1)^2 (x - 2)" + m, method);
        verdict(zero_double.p, CPS_SQUAREFREE_DISPROVED, "x^6 - x^2 = x^2 (x^4 - 1), in x^2" + m, method);
        verdict(zero_twice.p, CPS_SQUAREFREE_DISPROVED, "x^2 (x - 3)" + m, method);
        verdict(zero_once.p, CPS_SQUAREFREE_PROVED, "x (x - 3)" + m, method);
    }
    verdict(quadruple.p, CPS_SQUAREFREE_DISPROVED, "(x - 1)^4");
    verdict(composed.p, CPS_SQUAREFREE_PROVED, "x^4 - 3, in x^4");
    verdict(trinomial.p, CPS_SQUAREFREE_PROVED, "x^11 - 3x + 1, real roots", CPS_METHOD_AUTO, true);
    verdict(double_root.p, CPS_SQUAREFREE_DISPROVED, "(x - 1)^2 (x - 2), real roots", CPS_METHOD_AUTO, true);
    cps_polynomial *gaussian = nullptr; // (x - i)^2
    ok(cps_polynomial_create(&gaussian), "create");
    ok(cps_polynomial_set_coefficient(gaussian, 0, "-1", "0"), "set");
    ok(cps_polynomial_set_coefficient(gaussian, 1, "0", "-2"), "set");
    ok(cps_polynomial_set_coefficient(gaussian, 2, "1", "0"), "set");
    verdict(gaussian, CPS_SQUAREFREE_DISPROVED, "(x - i)^2");
    cps_polynomial_destroy(gaussian);
    cps_polynomial *ball = nullptr; // x^2 - 3x + (2 +- 1e-9): nothing exact to prepare
    ok(cps_polynomial_create(&ball), "create");
    ok(cps_polynomial_set_ball_coefficient(ball, 0, "2", "0", "1e-9", "0"), "set");
    ok(cps_polynomial_set_coefficient(ball, 1, "-3", "0"), "set");
    ok(cps_polynomial_set_coefficient(ball, 2, "1", "0"), "set");
    Result inexact;
    inexact.r = solved(ball);
    require(cps_result_squarefree(inexact.r) == CPS_SQUAREFREE_UNKNOWN, "a ball input got a square-free verdict");
    cps_polynomial_destroy(ball);

    // Two roots 2^-200 apart, 1 and 1 + 2^-200 (2^200 x^2 - (2^201 + 1) x + 2^200 + 1), under a 64-bit cap: the
    // cap the diagnostics name, as a field; none without the stop.
    const std::string two200 = "1606938044258990275541962092341162602522202993782792835301376";
    Polynomial close({"1606938044258990275541962092341162602522202993782792835301377",
                      "-3213876088517980551083924184682325205044405987565585670602753", two200.c_str()});
    Result capped, free;
    capped.r = solved(close.p, CPS_METHOD_AUTO, 64);
    free.r = solved(close.p);
    require(cps_result_termination(capped.r) == CPS_PRECISION_LIMIT && cps_result_stop_reason(capped.r) == CPS_STOP_PRECISION,
            "the close pair did not stop at the cap");
    const std::string text = cps_result_diagnostics(capped.r);
    const auto at = text.find("so a cap of ");
    require(at != std::string::npos, "the precision stop names no sufficient cap: " + text);
    const uint64_t named = std::stoull(text.substr(at + std::strlen("so a cap of ")));
    require(named > 64 && cps_result_sufficient_work_bits(capped.r) == named,
            "the sufficient cap " + std::to_string(cps_result_sufficient_work_bits(capped.r)) +
                " is not the one the diagnostics name, " + std::to_string(named));
    require(cps_result_squarefree(capped.r) == CPS_SQUAREFREE_PROVED, "the close pair is square-free");
    require(cps_result_termination(free.r) == CPS_SUCCESS && cps_result_sufficient_work_bits(free.r) == 0,
            "a solve that did not stop at its cap has a sufficient cap");

    // The hand-off: forced piecewise under a 64-bit cap on x^60 + 5x^30 - 3x + 1 reaches its precision budget.
    std::vector<const char *> sixty(61, "0");
    sixty[0] = "1", sixty[1] = "-3", sixty[30] = "5", sixty[60] = "1";
    Polynomial handed(sixty);
    Result h;
    h.r = solved(handed.p, CPS_METHOD_PIECEWISE, 64);
    require(cps_result_handoff(h.r) == CPS_HANDOFF_PRECISION &&
                std::string(cps_result_diagnostics(h.r)).find("piecewise precision budget to secular") != std::string::npos,
            "the hand-off is not reported as the precision budget's");

    // Programs: the modular image proves square-freeness or says nothing.
    const auto program_of = [](bool cubed) {
        cps_program *program = nullptr;
        ok(cps_program_create(&program), "program");
        uint64_t x = 0, two = 0, one = 0, square = 0, shifted = 0, power = 0, linear = 0, product = 0;
        ok(cps_program_variable(program, &x), "x");
        ok(cps_program_constant(program, "2", "0", &two), "2");
        ok(cps_program_constant(program, "1", "0", &one), "1");
        ok(cps_program_power(program, x, 2, &square), "x^2");
        ok(cps_program_subtract(program, square, two, &shifted), "x^2 - 2");
        if (cubed) ok(cps_program_power(program, shifted, 3, &power), "power");
        else power = shifted;
        ok(cps_program_subtract(program, x, one, &linear), "x - 1");
        ok(cps_program_multiply(program, power, linear, &product), "product");
        ok(cps_program_set_output(program, product), "output");
        return program;
    };
    for (const bool cubed : {false, true}) {
        cps_program *program = program_of(cubed);
        Result r;
        ok(cps_solve_program(program, nullptr, &r.r), "program solve");
        require(cps_result_squarefree(r.r) == (cubed ? CPS_SQUAREFREE_UNKNOWN : CPS_SQUAREFREE_PROVED),
                std::string(cubed ? "(x^2 - 2)^3 (x - 1)" : "(x^2 - 2)(x - 1)") + ": the program's verdict is " +
                    std::to_string(cps_result_squarefree(r.r)));
        cps_program_destroy(program);
    }
    require(cps_result_sufficient_work_bits(nullptr) == 0 && cps_result_squarefree(nullptr) == CPS_SQUAREFREE_UNKNOWN &&
                cps_result_handoff(nullptr) == CPS_HANDOFF_NONE,
            "the facts' NULL fallbacks");
}
} // namespace

int main() {
    try {
        option_getters();
        split_setters_answer_alike();
        facts();
        // Results of several kinds: a double root and a cluster, a solve stopped by its work budget, real roots,
        // a secular equation and a program.
        // -(x - 1)^4 (x + 2): degree 5, since Auto's closed forms answer degree 4 and below in a few polls.
        Polynomial repeated({"-2", "7", "-8", "2", "2", "-1"});
        Polynomial mixed({"1", "-3", "0", "0", "0", "0", "0", "0", "0", "0", "0", "1"});
        Options plain, tight, real;
        ok(cps_options_set_work_budget(tight.o, 1), "tight budget");
        ok(cps_options_set_real_only(real.o, 1), "real only");
        Result r1, r2, r3;
        ok(cps_solve(repeated.p, plain.o, &r1.r), "repeated");
        ok(cps_solve(mixed.p, tight.o, &r2.r), "tight");
        ok(cps_solve(mixed.p, real.o, &r3.r), "real");
        require(cps_result_stop_reason(r2.r) == CPS_STOP_WORK_BUDGET, "the tight budget did not stop the solve");
        blocks_agree(r1.r, "(x - 1)^4 (x + 2)");
        blocks_agree(r2.r, "a work-budget stop");
        blocks_agree(r3.r, "real roots");
        statistics_protocol(r1.r);
        cps_secular_term_view terms[2] = {{"2/3", "0", "0", "0"}, {"-2/3", "0", "3", "0"}};
        cps_secular *secular = nullptr;
        ok(cps_secular_create(terms, 2, nullptr, &secular), "secular");
        Result r4;
        ok(cps_solve_secular(secular, nullptr, &r4.r), "secular solve");
        blocks_agree(r4.r, "secular");
        cps_secular_destroy(secular);
        // (x^2 - 2)^3 (x - 1), built node by node.
        cps_program *program = nullptr;
        ok(cps_program_create(&program), "program");
        uint64_t x = 0, two = 0, one = 0, square = 0, shifted = 0, cube = 0, linear = 0, product = 0;
        ok(cps_program_variable(program, &x), "x");
        ok(cps_program_constant(program, "2", "0", &two), "2");
        ok(cps_program_constant(program, "1", "0", &one), "1");
        ok(cps_program_power(program, x, 2, &square), "x^2");
        ok(cps_program_subtract(program, square, two, &shifted), "x^2 - 2");
        ok(cps_program_power(program, shifted, 3, &cube), "cube");
        ok(cps_program_subtract(program, x, one, &linear), "x - 1");
        ok(cps_program_multiply(program, cube, linear, &product), "product");
        ok(cps_program_set_output(program, product), "output");
        Result r5;
        ok(cps_solve_program(program, nullptr, &r5.r), "program solve");
        blocks_agree(r5.r, "program");
        cps_program_destroy(program);
        // Forced piecewise under a 64-bit cap hands x^60 + 5x^30 - 3x + 1 to secular, which then answers.
        std::vector<const char *> sixty(61, "0");
        sixty[0] = "1", sixty[1] = "-3", sixty[30] = "5", sixty[60] = "1";
        Polynomial handed(sixty);
        Options forced;
        ok(cps_options_set_method(forced.o, CPS_METHOD_PIECEWISE), "piecewise");
        ok(cps_options_set_max_work_bits(forced.o, 64), "cap 64");
        Result r6;
        ok(cps_solve(handed.p, forced.o, &r6.r), "hand-off solve");
        blocks_agree(r6.r, "a piecewise hand-off to secular", true);
        polls_match_the_callback(mixed.p, "x^11 - 3x + 1");
        polls_match_the_callback(repeated.p, "(x - 1)^4 (x + 2)");
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
    std::cout << "the 0.5 calls hold what the engines counted, and the poll counter stops where a callback does\n";
    return 0;
}
