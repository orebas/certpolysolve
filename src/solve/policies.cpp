#include "solve/policies.hpp"
#include <stdexcept>

namespace certpolysolve::detail {
EnginePolicy complete_global_policy(cps_method method) {
    auto policy = maintained_method_policy(method);
    policy.research_integrated_dovetail = true;
    policy.research_dovetail_start = EnginePolicy::DovetailStart::PrecisionCluster;
    policy.research_dovetail_share_log2 = 3; // fixed 1/8 of priced practical service
    policy.research_dovetail_growth = 1;
    policy.research_service_meter = true;
    policy.research_identity_shift_meter = true;
    policy.research_local_packet_meter = true;
    policy.complete_fallback = true;
    return policy;
}

EnginePolicy standalone_subdivision_policy() {
    EnginePolicy policy;
    policy.backend = false;
    policy.prepared_subdivision = true;
    return policy;
}

EnginePolicy maintained_method_policy(cps_method method) {
    auto policy = production_policy();
    policy.complex_centering = policy.centering_refinement =
        method == CPS_METHOD_CLASSICAL || method == CPS_METHOD_SECULAR;
    policy.scale_invariant_centering = method == CPS_METHOD_CLASSICAL || method == CPS_METHOD_SECULAR;
    policy.centred_origin_basis = method == CPS_METHOD_CLASSICAL || method == CPS_METHOD_SECULAR;
    policy.pellet_real_shift = method == CPS_METHOD_CLASSICAL || method == CPS_METHOD_SECULAR;
    policy.centering_keeps_origin_scales = method == CPS_METHOD_SECULAR;
    if (method == CPS_METHOD_CLASSICAL) {
        policy.study = StudyMethod::Aberth;
        policy.refinement = RefinementMethod::Aberth;
        policy.classical_error_stopping = true;
        policy.defer_classical_validation = true;
        policy.newton_point_refinement = policy.candidate_allow_underflow = true;
        policy.coefficient_newton_bound = true;
    } else if (method == CPS_METHOD_SECULAR) {
        policy.study = StudyMethod::Secular;
        policy.refinement = RefinementMethod::Secular;
        // First try classical's machine/Newton certificate on every source; a
        // failed attempt whose machine run converged seeds the preliminary
        // phase ("variant C", the trade-off recorded in
        // docs/method-composition-design-20260928.md).
        // History: docs/code-history-20261021.md, session.cpp 1.
        policy.newton_point_refinement = policy.candidate_allow_underflow = true;
        policy.coefficient_newton_bound = policy.newton_point_handover = true;
        policy.model_stopping = SecularModelPolicy::Continuation;
        policy.regeneration.value_only = policy.regeneration.reuse_products = true;
        policy.regeneration.remember_precision = true;
        policy.regeneration.integer_node_products = true;
        policy.regeneration.limb_horner = true;
        policy.model_driven_refresh = policy.secular_preliminary = true;
        // A group of nodes that crawls onto a cluster whose root count the
        // model certifies exactly is re-seeded from its local model at the
        // roots' mean (EnginePolicy::geometric_local_restart).
        // History: docs/code-history-20261021.md, session.cpp 2.
        policy.geometric_local_restart = true;
        // A round whose active nodes iterate at different tiers sweeps the
        // ones within the initial stage in binary64 and the others in the
        // multiprecision engine, and a binary64 packet leaves the points it
        // does not move as they were (EnginePolicy::mixed_packets and
        // machine_keeps_unmoved).
        // History: docs/code-history-20261021.md, session.cpp 3.
        policy.mixed_packets = policy.machine_keeps_unmoved = true;
        // The binary64 packet's sums in lanes (EnginePolicy::packet_lane_sums), faster than the
        // one-accumulator loop. They change the sums' rounding, and with it which inputs' binary64
        // iteration ends with a node on the wrong cluster; they are on because the restart regroups such a
        // component (SecularRefinement::regroup_restart).
        // History: docs/code-history-20261021.md, session.cpp 4.
        policy.packet_lane_sums = true;
        // At a fresh model's floor the binary64 packet's sub-floor steps are
        // taken node-relative, so the validation there meets nodes 53 bits
        // closer (EnginePolicy::machine_node_offsets).
        // History: docs/code-history-20261021.md, session.cpp 5.
        policy.machine_node_offsets = true;
    } else if (method != CPS_METHOD_AUTO) {
        throw std::invalid_argument("method does not use the classical/secular driver");
    }
    return policy;
}

} // namespace certpolysolve::detail
