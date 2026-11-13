#include "solve/local_model_solver.hpp"
#include "solve/solve_state.hpp"

namespace certpolysolve::detail {
cps_result solve_local_model(DensePolynomial &&input, const cps_options &options) {
    auto policy = maintained_method_policy(CPS_METHOD_SECULAR);
    policy.geometric_local_restart = false;
    SolveState state(std::move(input), std::move(policy));
    return state.advance(&options);
}
#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
cps_result solve_local_model_numerical(DensePolynomial &&input, const cps_options &options) {
    auto policy = maintained_method_policy(CPS_METHOD_SECULAR);
    policy.geometric_local_restart = false;
    policy.dovetail_subdivision = false;
    policy.research_no_subdivision = true;
    SolveState state(std::move(input), std::move(policy));
    return state.advance(&options);
}

#endif
cps_result solve_local_model_priced(DensePolynomial &&input, const cps_options &options) {
    ResearchSchedulingMetrics metrics;
    auto policy = maintained_method_policy(CPS_METHOD_SECULAR);
    policy.geometric_local_restart = false;
    policy.research_service_meter = true;
    policy.research_identity_shift_meter = true;
    policy.research_local_packet_meter = true;
    policy.research_scheduling_metrics = &metrics;
    // The state's lifetime ends before its stack-owned observer. Its settings
    // retain the caller's bounded work/sweep/cancellation options unchanged.
    SolveState state(std::move(input), std::move(policy));
    auto result = state.advance(&options);
    result.research_priced_work_available = true;
    result.research_scheduling_metrics = metrics;
    return result;
}
} // namespace certpolysolve::detail
