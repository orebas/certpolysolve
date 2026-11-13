#pragma once
// A session's statistics are cumulative over its calls (certpolysolve.h), across a hand-over to another engine too:
// the engine that takes over counts from zero, so the session adds what the engines before it counted. Used by
// the session's piecewise-to-secular hand-off (session.cpp) and by the real-only route's attempts
// (real_only_route.cpp). History: docs/code-history-20261021.md, session_statistics.hpp 1.
#include "core/query.hpp"
#include "core/statistics_accumulation.hpp"
#include <algorithm>

namespace certpolysolve::detail {

struct CarriedStatistics {
    // The v1 block is carried apart (into cps_result::carried_v1): the old getter keeps returning the answering
    // engine's own v1, while the 0.5 totals (cps_statistics.preparation_ns and the engine counters) add the engines
    // before it (C10). exact_factors is not carried: it counts the factors of the answering preparation.
    StageStatistics v1{};
    SearchStatistics v2{};
    CandidateStatistics v3{};
    uint64_t highest_work_bits = 0;
    // What the engine that hands over had counted.
    void carry(const cps_result &before) {
        highest_work_bits = std::max({highest_work_bits, before.progress.max_outer_bits, before.work_bits});
        for (auto field : statistics::carried_stage_fields)
            v1.*field = resources::add(v1.*field, resources::add(before.statistics.*field, before.carried_v1.*field));
        statistics::add_saturating(v2, before.statistics_v2, statistics::search_fields);
        statistics::add_saturating(v3, before.statistics_v3, statistics::candidate_fields);
    }
    // The session's count: the engines before, plus the current one's.
    void add_to(cps_result &result) const {
        result.progress.max_outer_bits = std::max({result.progress.max_outer_bits, result.work_bits, highest_work_bits});
        statistics::add_saturating(result.carried_v1, v1, statistics::carried_stage_fields);
        statistics::add_saturating(result.statistics_v2, v2, statistics::search_fields);
        statistics::add_saturating(result.statistics_v3, v3, statistics::candidate_fields);
    }
};

} // namespace certpolysolve::detail
