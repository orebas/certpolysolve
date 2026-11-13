#pragma once
// A factor of the input as the global driver (SolveState, solver.cpp) owns it: its exact polynomial in the
// input's and in centred coordinates, the states of the methods that work on it (candidates, subdivision, the
// real engines, the secular refinement), its centring (refine_center, prepare_secular_center), the separation
// bound a precision stop names, and the enclosers and shared polynomials the structured sources are built from.
// History: docs/code-history-20261021.md, factor_state.hpp 1.
#include "methods/classical/classical_isolation.hpp"
#include "core/dovetail.hpp"
#include "numeric/exact/exact_complex.hpp"
#include "core/cover.hpp"
#include "input/coefficients/polynomial_input.hpp"
#include "core/query.hpp"
#include "core/engine_policy.hpp"
#include "numeric/exact/preparation_diagnostics.hpp"
#include "core/raii.hpp"
#include "numeric/candidates/real_line_seeds.hpp"
#include "methods/real/real_refinement.hpp"
#include "methods/secular/secular_refinement.hpp"
#include "numeric/evaluate/structured_source.hpp"
#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
#include "sturm.hpp"
#endif
#include "methods/subdivision/subdivision.hpp"
#include <chrono>
#include <cstdint>
#include <memory>
#include <optional>
#include <vector>

namespace certpolysolve::detail {
struct FactorState {
    bool degree_policy_reported = false;
    uint64_t id = 0;
    IntegerPolynomial integer;
    GaussianPolynomial gaussian;
    enum Kind { Ball, Real, Gaussian } kind = Ball;
    // A Ball factor whose coefficients the modular test proved square-free as written, where the compact-exponent
    // rule refused exact preparation (SolveState::Impl::preprocess). Real and Gaussian factors are square-free by
    // construction.
    bool squarefree_proved = false;
    uint64_t multiplicity = 1, precision = 64, real_count = 0, separation_bits = 0;
    // The solve raised separation_bits because this factor's exported cover overlapped another factor's (the overlap
    // branch of solver.cpp): a power composite then counts as finished only once its reduced cover meets the raised
    // goal too (factor_goal). Never cleared.
    bool overlap_tightened = false;
    CandidateState candidates;
    bool ready = false, blocked = false, local_refinement = false;
    // The standalone route first isolates at bounded starting precision,
    // then reopens that same search for the caller's requested accuracy.
    bool subdivision_isolated = false;
    bool complete_fallback = false;
    uint64_t exact_height = 0;
    std::vector<Root> cover;
    ulong power = 1;
    std::vector<Root> lifted_cover;
    const std::vector<Root> &output_cover() const { return power == 1 ? cover : lifted_cover; }
    std::unique_ptr<Subdivision> search;
    // The persistent Subdivision run beside the factor's engine: the dovetail's (below), and the research Portfolio
    // engine's (CandidateEngine::Portfolio, every subdivision_period of its portfolio_steps).
    std::unique_ptr<Subdivision> dovetail_search;
    std::unique_ptr<RealRefinement> real_search;
    bool research_numerical_precision_increased = false;
#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
    // The research build's Portfolio engine: its Sturm step and its subdivision turns.
    std::unique_ptr<SturmRefinement> sturm;
    uint64_t portfolio_steps=0, real_stagnant=0;
    bool sturm_unavailable=false;
#endif
    // Dovetailing (policy.dovetail_subdivision): the persistent Subdivision is
    // dovetail_search. dovetail_stop names the budget every one of its
    // components awaits, if any.
    DovetailSchedule schedule;
    bool dovetail_engaged = false;
    cps_stop_reason dovetail_stop = CPS_STOP_NONE;
    // Blocked by the input: every active node's residue is the input's
    // (regenerate_residues, input_limited), which no precision changes.
    // History: docs/code-history-20261021.md, factor_state.hpp 2.
    bool input_limited = false;
    // The factor's next step needs sweeps and the solve's iteration budget
    // (cps_options::max_iterations) is spent: it is skipped until the solve ends.
    bool iteration_wait = false;
    IntegerPolynomial centered;
    GaussianPolynomial centered_gaussian;
    GaussianRational center;
    uint64_t center_precision = 0;
    bool have_center = false, center_rejected = false;
    std::unique_ptr<SecularRefinement> structured;
    // The all-real seed pass's outcome for this factor (cps_options::real_line_seeds), computed once, on the
    // polynomial the start seeds from (centred or not); its seeds, increasing binary64 values, replace the circles'
    // first points.
    struct RealLineSeedState {
        bool computed = false, centred = false;
        std::vector<double> seeds;
        RealLineSeedStatistics statistics;
    } real_line_seeds;
};
// Shared by explicit secular studies and admitted production factors. The
// transform is exact; source ownership and RootCoordinates keep certificates
// in the caller's coordinates. Candidate points are only iteration hints.
void observe_preparation(const EnginePolicy &policy, const char *stage,
                         const PreparationDiagnostic &diagnostic, uint64_t budget,
                         std::chrono::steady_clock::time_point started, uint64_t active = 0);
// What a precision stop tells the caller (the default cap is 65,536 bits; a caller who needs more raises
// it). History: docs/code-history-20261021.md, factor_state.hpp 4. For a square-free factor
// with integer coefficients, degree n >= 2, every two roots are more than sqrt(3) n^(-(n+2)/2)
// ||f||_2^(-(n-1)) apart (Mahler 1964) and at most 1 + max |a_i / a_n| in size (Cauchy). From bit
// lengths, rounded up: ||f||_2 <= sqrt(n + 1) 2^h with h the largest coefficient's bits, and the sqrt(3)
// dropped. A precision that resolves the roots' size and their separation plus the goal resolves every
// root; the ladder that raises the precision at most doubles it, so the cap named is twice that sum
// (benchmarks/precision-cap/separation_bound.py, benchmark_cap, the same rule for the benchmarks).
struct SeparationBound {
    uint64_t separation_bits = 0, size_bits = 0;
};
std::optional<SeparationBound> separation_bound(const fmpz_poly_t p);
// Moves the factor's centre to its residual root centroid -q_{n-1}/(n q_n), rounded like the first centre, while
// that at least halves the geometric mean of the root distances, about a bit a root (the rationale is with the
// definition).
void refine_center(FactorState &factor, const cps_options &options, certpolysolve::detail::SearchStatistics &statistics,
                   const EnginePolicy &policy);
// The maintained secular recipe's centring of a real (or, under complex_centering, Gaussian) factor before its
// secular stage: an exact Taylor shift to a rounded centre, accepted by the height precheck.
void prepare_secular_center(FactorState &factor, const acb_poly_t polynomial,
                            const cps_options &options, certpolysolve::detail::SearchStatistics &statistics,
                            const EnginePolicy &policy);
// Enclosers of the factor's own polynomial (owned: they outlive the solve's temporaries).
PolynomialEncloser owned_encloser(const FactorState &factor, const DensePolynomial &input);
ExactTaylorEncloser owned_exact_taylor(const FactorState &factor, uint64_t budget, bool allow_expansion, uint64_t operations
                                     , std::function<double()> *work_reader = nullptr
                                     );
// The exact polynomial a structured source evaluates, shared by its
// coefficient encloser and the exact value at a node that residue
// regeneration takes at the precision budget (one copy, as the encloser
// alone made before): the centred factor under a centred study, in whose
// coordinates its nodes live.
struct SharedFactorPolynomial {
    std::shared_ptr<const IntegerPolynomial> integer;
    std::shared_ptr<const GaussianPolynomial> gaussian;
    PolynomialEncloser encloser() const {
        if (integer)
            return [polynomial = integer](acb_poly_t out, slong bits) {
                acb_poly_set_fmpz_poly(out, polynomial->value(), bits);
            };
        return [polynomial = gaussian](acb_poly_t out, slong bits) { enclose_gaussian(out, *polynomial, bits); };
    }
};
std::optional<SharedFactorPolynomial> shared_factor_polynomial(const FactorState &factor, bool centered);
// The uncentred factor and its centre for a centred source (PolynomialSourceOptions::origin), with both
// polynomials' log2 coefficient sizes. Nothing when the centre is not exactly representable (centres are
// low-denominator dyadics, round_center).
void attach_origin(PolynomialSourceOptions &options, const FactorState &factor);
} // namespace certpolysolve::detail
