#pragma once
#include "numeric/candidates/candidates.hpp"
#include "core/engine_config.hpp"
#include <memory>

namespace certpolysolve::detail {
#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
// What an engine is and holds, for the research build's Portfolio trace (research/solvers/portfolio.cpp).
struct EngineCapabilities {
    bool coefficients = true, second_derivative = false, numerical_deflation = false;
    bool linear_storage = true, resumable = true;
    uint64_t preparation_bytes = 0, preparation_work = 0;
};
#endif
struct EngineAdvance {
    CandidateProgress candidates;
    bool arithmetic_failure = false, arithmetic_limited = false, exhausted = false;
    uint64_t active = 0, work = 0;
};
// A persistent numerical engine owns immutable prepared coefficients and all
// iteration state. It exports proposals, never Root/certificate objects.
class NumericalEngine {
  public:
    virtual ~NumericalEngine() = default;
#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
    virtual EngineCapabilities capabilities() const = 0;
#endif
    virtual EngineAdvance advance(unsigned sweeps, const cps_options &, certpolysolve::detail::CandidateStatistics &) = 0;
    virtual void export_candidates(acb_ptr output) const = 0;
};
std::unique_ptr<NumericalEngine> prepare_machine_engine(IterationMethod, ArithmeticTier,
                                                        const acb_poly_t polynomial, ConstRootView initial,
                                                        bool compensated = false, bool jit = false,
                                                        bool allow_underflow = false,
                                                        uint64_t storage_limit = UINT64_MAX,
                                                        bool nondyadic_scaling = false,
                                                        bool interaction_tree = false);
#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
// The research build's companion engines (research/solvers/structured_companion.cpp, dense_companion.cpp).
std::unique_ptr<NumericalEngine> prepare_structured_engine(IterationMethod, ArithmeticTier,
    const acb_poly_t, ConstRootView, uint64_t storage_limit = UINT64_MAX);
bool dense_eigen_available();
std::unique_ptr<NumericalEngine> prepare_dense_engine(const acb_poly_t, ConstRootView,
                                                      uint64_t storage_limit, uint64_t work_limit);
// The research build's machine engine for every request prepare_machine_engine does not serve (the other methods, tiers
// and evaluators, and the companions' dispatch): the engine as it was before W8's S10
// (research/solvers/machine_methods.cpp).
std::unique_ptr<NumericalEngine> prepare_research_machine_engine(IterationMethod, ArithmeticTier,
                                                                 const acb_poly_t polynomial, ConstRootView initial,
                                                                 bool compensated, bool jit, bool allow_underflow,
                                                                 uint64_t storage_limit, bool nondyadic_scaling,
                                                                 bool interaction_tree);
#endif
} // namespace certpolysolve::detail
