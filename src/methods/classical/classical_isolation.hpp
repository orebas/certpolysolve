#pragma once
// The classical route's isolation of one factor's roots from its candidates: the candidate state a factor
// keeps across its precision passes, isolate_complex (the candidates of one pass into certified roots:
// inclusion disks, the linear and binomial special cases; classified and published in the caller's
// coordinates). History: docs/code-history-20261021.md, classical_isolation.hpp 1.
#include "core/cover.hpp"
#include "core/query.hpp"
#include "core/engine_policy.hpp"
#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
#include "portfolio.hpp"
#endif
#include "core/raii.hpp"
#include "numeric/evaluate/structured_source.hpp"
#include "numeric/candidates/study_candidates.hpp"
#include <cstdint>
#include <limits>
#include <memory>
#include <vector>

namespace certpolysolve::detail {
struct CandidateState {
    std::unique_ptr<StudyCandidates> study;
#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
    std::unique_ptr<CandidatePortfolio> portfolio; // the research build's Portfolio engine
#endif
    std::unique_ptr<RootVector> centers;
    slong precision = 0;
    // The precision of the certificate attempt's certificates when it
    // succeeded (it may raise its own above the pass's).
    slong cover_bits = 0;
    unsigned packets = 0;
    // The sweeps those packets were allowed at this precision, heuristics::classical_packet_sweeps each (a packet
    // that the solve's iteration budget, cps_options::max_iterations, cuts short is not one of them).
    uint64_t charged = 0;
    // The iteration budget was spent before a packet could finish: no
    // verdict on this precision, and the caller stops the solve for that budget.
    bool budget_stop = false;
    double correction_log2 = std::numeric_limits<double>::infinity();
    bool repeat = false;
#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
    bool seeded = false; // the research build's DoubleSeed engine ran binary64_seeds
#endif
    bool newton_attempted = false;
    ComplexBall unresolved_center;
    slong unresolved_count = 0;
    std::vector<uint8_t> inactive;
    uint64_t last_goal_bits = 0;
    unsigned candidate_attempts = 0;
};

// The roots of `polynomial` (a factor, `multiplicity` times) from the candidates in `state` at `precision`, into
// `result`; true when every root is certified, false when this precision does not settle them (the caller raises
// it or changes method). With `allow_clusters` a cluster may be published as one counted region. The exact
// encloser, when given, refines the certificates' evaluations.
bool isolate_complex(std::vector<Root> &result, const acb_poly_t polynomial, uint64_t multiplicity,
                     CandidateState &state, const cps_options &options, slong precision, bool allow_clusters,
                     certpolysolve::detail::SearchStatistics &statistics, certpolysolve::detail::CandidateStatistics &candidate_statistics,
                     CandidateEngine engine, bool analytic_binomial, const EnginePolicy *policy,
                     const PolynomialEncloser *exact = nullptr);
} // namespace certpolysolve::detail
