#include "solve/query_execution.hpp"

namespace certpolysolve::detail {
QueryExecution::QueryExecution(const QueryInput &source, const cps_options *proposed)
    : next(proposed ? *proposed : cps_options{}) {
    if (!proposed) { next.cancel = nullptr; next.cancel_data = nullptr; }
    // The default goal for the input kind (with_resolved_goal), resolved first.
    next = with_resolved_goal(next, source.program_kind);
    degree = source.input_degree();
    zero_input = !source.program && source.input_length() == 0;
    // A work budget the caller left at its default grows with the degree
    // (call_work_units).
    // History: docs/code-history-20261021.md, session.cpp 7.
    next.work_units = call_work_units(next, degree);
    if (!source.supported(next)) {
        cps_result unsupported;
        // A program's degree bound can exceed int64 (it saturates); an input with no nonzero coefficient is -1.
        unsupported.degree =
            zero_input ? -1 : static_cast<int64_t>(std::min<uint64_t>(degree, std::numeric_limits<int64_t>::max()));
        unsupported.goal = next.goal;
        unsupported.termination = CPS_UNSUPPORTED_QUERY;
        unsupported.progress.stop_reason = CPS_STOP_INPUT;
        unsupported.diagnostics =
            next.goal == CPS_ROOT_ISOLATION && !next.real_only
                ? "complex isolation requires Auto and exact coefficient input"
            : next.goal == CPS_ROOT_ISOLATION || next.region.set()
                ? "the isolation goal and a region are answered for real-only Auto queries on exact real "
                  "polynomial input only"
            : source.program_kind ? "program input answers all-complex queries by Auto, classical and secular; "
                                    "real-only queries, the isolation goal, regions and the piecewise method are "
                                    "not available for programs"
                                  : "explicit methods require exact coefficients and a complex accuracy goal; "
                                    "native secular input supports Auto/Secular";
        refusal = std::move(unsupported);
        return;
    }
    // A region's ends as exact rationals, before any work (real_region.hpp): an end the exact storage policy
    // refuses makes the query unsupported.
    if (next.region.set()) {
        region = source.input ? real_query::region(next, *source.input) : real_query::region(next, *source.terms);
        if (region->refusal()) {
            cps_result unsupported;
            unsupported.degree = zero_input ? -1 : static_cast<int64_t>(degree);
            unsupported.goal = next.goal;
            unsupported.termination = CPS_UNSUPPORTED_QUERY;
            unsupported.progress.stop_reason = CPS_STOP_INPUT;
            unsupported.diagnostics = region->refusal();
            refusal = std::move(unsupported);
            return;
        }
    }
    cancellation_.callback = next.cancel;
    cancellation_.data = next.cancel_data;
    next.cancel = Cancellation::poll;
    next.cancel_data = &cancellation_;
}
} // namespace certpolysolve::detail
