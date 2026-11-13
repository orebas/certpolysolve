#include "solve/prepared_piecewise.hpp"
// Counter contracts at factor aggregation and engine handoff. In particular,
// ordinary unsigned sums, saturating handoff sums, maxima and local diagnostics
// must retain their different meanings when their field descriptions are shared.
#include "numeric/arithmetic/polynomial_arithmetic.hpp"
#include "methods/piecewise/piecewise_statistics.hpp"
#include "core/session_statistics.hpp"
#include "methods/piecewise/piecewise_roots.hpp"
#include "solve/auto_selection.hpp"
#include <certpolysolve/certpolysolve.hpp>
#include <array>
#include <bit>
#include <cstdint>
#include <iostream>
#include <stdexcept>

namespace d = certpolysolve::detail;
namespace {
void require(bool condition, const char *message) {
    if (!condition) throw std::runtime_error(message);
}
template <class T> using Counters = std::array<uint64_t, sizeof(T) / sizeof(uint64_t)>;
template <class T> T filled(uint64_t value) {
    Counters<T> fields;
    fields.fill(value);
    return std::bit_cast<T>(fields);
}
template <class T> Counters<T> counters(const T &value) {
    return std::bit_cast<Counters<T>>(value);
}
// Together with the type's size check, uniqueness prevents a duplicated field
// from hiding an omitted one, including overlap between operation categories.
template <class T, size_t... N>
consteval bool unique_fields(uint64_t T::*const (&...groups)[N]) {
    std::array<uint64_t T::*, (N + ...)> fields{};
    size_t next = 0;
    const auto append = [&](const auto &group) {
        for (const auto field : group) fields[next++] = field;
    };
    (append(groups), ...);
    for (size_t i = 0; i < fields.size(); ++i)
        for (size_t j = i + 1; j < fields.size(); ++j)
            if (fields[i] == fields[j]) return false;
    return true;
}
static_assert(unique_fields(d::statistics::carried_stage_fields));
static_assert(unique_fields(d::statistics::search_fields));
static_assert(unique_fields(d::statistics::candidate_fields));
static_assert(unique_fields(d::statistics::model_error_sum_fields));
static_assert(unique_fields(d::statistics::chart_sum_fields, d::statistics::chart_local_fields));
static_assert(unique_fields(d::statistics::piecewise_sum_fields, d::statistics::piecewise_max_fields,
                            d::statistics::piecewise_local_fields));

void handoff() {
    cps_result before;
    before.statistics = filled<d::StageStatistics>(1);
    before.statistics.isolation_ns = UINT64_MAX - 2;
    before.statistics.exact_factors = 99;
    before.carried_v1 = filled<d::StageStatistics>(3);
    before.statistics_v2 = filled<d::SearchStatistics>(5);
    before.statistics_v2.work_units = UINT64_MAX - 3;
    before.statistics_v3 = filled<d::CandidateStatistics>(7);
    before.statistics_v3.evaluations = UINT64_MAX - 5;
    before.progress.max_outer_bits = 256;
    before.work_bits = 128;

    d::CarriedStatistics history;
    history.carry(before);
    history.carry(before);
    cps_result after;
    after.statistics = filled<d::StageStatistics>(5);
    after.statistics.exact_factors = 7;
    const auto current_stage = counters(after.statistics);
    after.carried_v1 = filled<d::StageStatistics>(1);
    after.carried_v1.exact_factors = 11;
    after.statistics_v2 = filled<d::SearchStatistics>(1);
    after.statistics_v3 = filled<d::CandidateStatistics>(1);
    after.progress.max_outer_bits = 64;
    after.work_bits = 32;
    history.add_to(after);
    require(after.progress.max_outer_bits == 256 && after.work_bits == 32,
            "handoff lost the earlier highest precision or changed the final precision");

    auto stage = filled<d::StageStatistics>(9);
    stage.isolation_ns = UINT64_MAX;
    stage.exact_factors = 11;
    require(counters(after.statistics) == current_stage, "handoff changed the current engine's stage block");
    require(counters(after.carried_v1) == counters(stage), "handoff lost stage history or carried exact_factors");
    require(history.v1.exact_factors == 0, "history carried an earlier preparation's factor count");
    auto search = filled<d::SearchStatistics>(11);
    search.work_units = UINT64_MAX;
    require(counters(after.statistics_v2) == counters(search), "search counters did not saturate across handoffs");
    auto candidate = filled<d::CandidateStatistics>(15);
    candidate.evaluations = UINT64_MAX;
    require(counters(after.statistics_v3) == counters(candidate), "candidate counters did not saturate across handoffs");

    // A subsequent handoff must carry both the current and already-carried
    // stage work, while still excluding the answering preparation's factors.
    d::CarriedStatistics nested;
    nested.carry(after);
    stage = filled<d::StageStatistics>(14);
    stage.isolation_ns = UINT64_MAX;
    stage.exact_factors = 0;
    require(counters(nested.v1) == counters(stage), "a nested handoff lost or double-counted stage work");
    require(counters(nested.v2) == counters(search) && counters(nested.v3) == counters(candidate),
            "a nested handoff changed the cumulative engine counters");
    cps_result last;
    last.work_bits = 512;
    nested.add_to(last);
    require(last.progress.max_outer_bits == 512, "handoff lost the final engine's higher precision");

    cps_engine_counters published{};
    published.size = sizeof(published);
    require(cps_result_get_engine_counters(&after, &published) == CPS_OK, "counter publication failed");
    require(published.input_enclosure_ns == 14 && published.isolation_ns == UINT64_MAX &&
                published.evaluations == UINT64_MAX,
            "public engine counters lost carried work or saturation");
    cps_statistics summary{};
    summary.size = sizeof(summary);
    require(cps_result_get_statistics(&after, &summary) == CPS_OK, "statistics publication failed");
    require(summary.exact_factors == 7 && summary.preparation_ns == 14 && summary.work_units == UINT64_MAX,
            "public statistics mixed the current factor count with historical preparations");
}

void public_handoff_precision() {
    // Piecewise uses 149 bits here before its regime exit; the successful
    // secular handoff uses 85. A public maximum must include both engines.
    certpolysolve::Polynomial polynomial;
    polynomial.set(0, 1); polynomial.set(1, -3); polynomial.set(10, 5); polynomial.set(20, 1);
    cps_options options;
    options.method = CPS_METHOD_PIECEWISE;
    options.goal = CPS_ROOT_APPROXIMATION;
    options.bits = 53;
    options.max_work_bits = 160;
    require(cps_options_set_work_budget(&options, 64) == CPS_OK, "explicit work budget refused");
    auto input = d::dense(*polynomial.get());
    auto config = d::maintained_piecewise_config();
    config.caller_method = true;
    config.certificate_attempt = false;
    config.certificate_attempt_fallback = 0;
    auto automatic = options;
    automatic.method = CPS_METHOD_AUTO;
    const auto selection = d::select_auto_route(input, automatic);
    if (selection.route == d::AutoRoute::Piecewise) config.regime_exit_estimate = selection.piecewise_operations;
    d::PreparedPiecewiseState piecewise(input, config);
    const auto before = piecewise.advance(options);
    require(piecewise.statistics().regime_exits && before.progress.max_outer_bits > 0,
            "precision fixture did not execute the piecewise regime exit");
    cps_result *raw = nullptr;
    const auto error = cps_solve(polynomial.get(), &options, &raw);
    const auto result = certpolysolve::detail::make_result(error, raw);
    require(result.goal_met() && cps_result_handoff(result.get()) == CPS_HANDOFF_REGIME_EXIT,
            "precision fixture did not complete through the handoff");
    const auto statistics = result.statistics();
    require(statistics.final_work_bits < before.progress.max_outer_bits,
            "precision fixture no longer finishes at a lower precision");
    require(statistics.highest_work_bits >= before.progress.max_outer_bits,
            "public statistics discarded the earlier engine's highest precision");
}

void ordinary_attempt_counters() {
    auto total = filled<d::CandidateStatistics>(UINT64_MAX - 2);
    const auto attempt = filled<d::CandidateStatistics>(5);
    d::statistics::add_ordinary(total, attempt, d::statistics::candidate_fields);
    require(counters(total) == counters(filled<d::CandidateStatistics>(2)),
            "attempt aggregation saturated, omitted or duplicated a counter");
    d::statistics::add_ordinary(total, total, d::statistics::candidate_fields);
    require(counters(total) == counters(filled<d::CandidateStatistics>(4)),
            "attempt self-aggregation changed its unsigned sums");
}

void piecewise_aggregation() {
    d::PiecewiseRootStatistics total, other;
    total.certificate = filled<d::ChartCertificateStatistics>(UINT64_MAX - 2);
    other.certificate = filled<d::ChartCertificateStatistics>(5);
    total.work.enclosures = 1.5;
    other.work.enclosures = 2.25;
    total.finishing_work.finisher = 3.5;
    other.finishing_work.finisher = 4.25;
    total.proposals = UINT64_MAX - 2;
    other.proposals = 5;
    total.rounds = 3;
    other.rounds = 4;
    total.discovery_count = 11;
    other.discovery_count = 13;
    total.final_accurate_count = 5;
    other.final_accurate_count = 7;
    total.maximum_degree = 17;
    other.maximum_degree = 9;
    total.model_bits = 8;
    other.model_bits = 16;
    total.arithmetic_bits = 100;
    other.arithmetic_bits = 64;
    // These diagnostics were not carried by factor aggregation before this
    // refactor. Preserve that contract; changing their reporting is separate.
    total.axis_gap_failures = total.axis_undecided = total.axis_unrefined = total.axis_sign_conflicts = 2;
    other.axis_gap_failures = other.axis_undecided = other.axis_unrefined = other.axis_sign_conflicts = 9;
    total.certificate.round_model_shortfall_bits = total.predicted_rungs_skipped = 3;
    other.certificate.round_model_shortfall_bits = other.predicted_rungs_skipped = 8;
    total.certificate.round_shortfall_log2.fill(4);
    other.certificate.round_shortfall_log2.fill(10);
    total += other;
    auto expected_certificate = filled<d::ChartCertificateStatistics>(2);
    expected_certificate.round_model_shortfall_bits = 3;
    expected_certificate.round_shortfall_log2.fill(4);
    require(counters(total.certificate) == counters(expected_certificate),
            "factor aggregation lost shared certificate counters or carried local round state");
    require(total.work.enclosures == 3.75 && total.finishing_work.finisher == 7.75,
            "factor aggregation changed counted work");
    require(total.proposals == 2 && total.rounds == 7 && total.discovery_count == 24 && total.final_accurate_count == 12,
            "factor aggregation changed ordinary unsigned addition");
    require(total.maximum_degree == 17 && total.model_bits == 16 && total.arithmetic_bits == 100,
            "factor aggregation summed a maximum");
    require(total.axis_gap_failures == 2 && total.axis_undecided == 2 && total.axis_unrefined == 2 &&
                total.axis_sign_conflicts == 2 && total.certificate.round_model_shortfall_bits == 3 && total.predicted_rungs_skipped == 3,
            "factor aggregation changed local diagnostics");
    for (const auto count : total.certificate.round_shortfall_log2)
        require(count == 4, "factor aggregation changed the current round's shortfall histogram");
    total += total;
    expected_certificate = filled<d::ChartCertificateStatistics>(4);
    expected_certificate.round_model_shortfall_bits = 3;
    expected_certificate.round_shortfall_log2.fill(4);
    require(counters(total.certificate) == counters(expected_certificate),
            "certificate self-aggregation changed sums or local round state");
    require(total.proposals == 4 && total.rounds == 14 && total.maximum_degree == 17 &&
                total.model_bits == 16 && total.axis_gap_failures == 2 && total.work.enclosures == 7.5,
            "self-aggregation changed its existing semantics");
}
} // namespace

int main() {
    certpolysolve::ThreadScope lifetime;
    try {
        public_handoff_precision();
        handoff();
        ordinary_attempt_counters();
        piecewise_aggregation();
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
