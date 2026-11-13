#pragma once
// Maintained controller recipes. EnginePolicy in core contains only the shared
// configuration data; lower modules never choose or call a controller recipe.
#include "core/engine_policy.hpp"

namespace certpolysolve::detail {
inline EnginePolicy production_policy() {
    EnginePolicy policy;
    // Qualified machine-to-adaptive candidate work uses the established source
    // and refinement choices. Explicit experiment policies keep their controls.
    policy.study_machine = true;
    policy.prepared_source_refinement = true;
    policy.early_power_composition = true;
    // No secular stall outlasts its budgets: once a secular factor has a
    // certified cover, the complete method (Subdivision) advances beside it
    // on a fixed share of the work. Inert where no factor runs secular.
    policy.dovetail_subdivision = true;
    return policy;
}
EnginePolicy maintained_method_policy(cps_method);
// Selected global candidate; the dispatcher enables it only in the separate
// Auto integration commit. Local solves keep their maintained inner recipe.
EnginePolicy complete_global_policy(cps_method);
// Exact coefficient isolation and the private standalone refinement driver.
// No candidate iteration, power reduction or numerical/fallback scheduler.
EnginePolicy standalone_subdivision_policy();
} // namespace certpolysolve::detail
