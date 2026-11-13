#pragma once
#include <limits>
#include "numeric/evaluate/evaluator.hpp"
#include "numeric/candidates/model_stopping.hpp"
#include <numbers>

namespace certpolysolve::detail {
class MonicSource;
enum class CandidateContinuation {
    Ordinary, NextPacketAtRequestedPrecision, PrecisionDecision
};
struct CandidateProgress {
    double correction_log2 = 0;
    bool finite = false;
    bool converged = false;
    bool cancelled = false;
    CandidateContinuation continuation = CandidateContinuation::Ordinary;
    bool model_exhausted = false, model_immediate = false;
    double service_evaluation = 0, service_interaction = 0;
    uint64_t service_scaled_binary64 = 0;
};
#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
// A bounded hardware-binary64 initializer, the research build's DoubleSeed test
// engine's start. Range loss, coincident proposals, cancellation, or nonfinite
// arithmetic returns false for ordinary fallback.
bool binary64_seeds(acb_ptr output, const acb_poly_t polynomial, const cps_options &options,
                    certpolysolve::detail::CandidateStatistics &statistics);
#endif
// binomial_seeds: a sparse polynomial's start, each Newton-polygon edge's seeds at its binomial's roots
// (heuristics::binomial_seed_spacing_fraction); otherwise equally spaced (placement_angle).
bool initial_candidates(acb_ptr output, const acb_poly_t polynomial, slong precision, bool binomial_seeds = false);
// The angle of the j-th of `count` points placed on the `circle`-th circle of
// a start: equally spaced, then rotated by an irrational offset so that the
// placement shares no symmetry with the source. A Jacobi sweep of a real
// source preserves conjugate symmetry exactly, so a placement that is
// conjugation-invariant about a real centre keeps its number of real-axis
// points for ever, whatever the number of real roots: the m-th unit roots
// (m even) put a conjugate pair on the perpendicular bisector of two real
// roots, and there it stays. The half-step alone, pi/count, is still
// conjugation-invariant; the offsets are not a multiple of pi/count for any
// count, and differ between circles.
// History: docs/code-history-20261021.md, candidates.hpp 1.
inline double placement_angle(slong j, slong count, slong circle) {
    return 2 * std::numbers::pi * (static_cast<double>(j) + 0.5) / static_cast<double>(count) +
           0.2718281828459045 + static_cast<double>(circle) * 0.6180339887498948;
}
// `count` points on the circle of `radius` about `center` at placement_angle
// (circle 1), as exact midpoints.
void circle_candidates(acb_ptr output, slong count, const acb_t center, const arb_t radius, slong precision);
// The value and derivative of the midpoint polynomial at z, as the midpoint
// Aberth of small packets computes them: a dense Horner, or, when fewer than
// a quarter of the coefficients are nonzero, powers over the gaps between the
// nonzero ones. Both are candidate arithmetic (floating point at prec, the
// balls' midpoints only), never a certificate. Exposed for its test.
void midpoint_polynomial_evaluation(acb_t value, acb_t derivative, acb_srcptr coefficients, slong length,
                                    const acb_t z, slong precision);
// The same placement for a model whose coefficients the working precision may
// not resolve. A coefficient enclosure with no correct bit says only how large
// the coefficient can be, so it enters the Newton polygon with its upper
// bound: the polygon is then an upper envelope, and its radii are the largest
// the roots can be on this evidence. determined reports that every vertex of
// the polygon is a resolved coefficient; only then do the roots of the
// midpoint polynomial say more than these scales do. For a determined model
// the points equal initial_candidates' (an unresolved coefficient under the
// polygon stays under it with its bound).
// `reach`, when given, receives the largest of the polygon's circles that an unresolved coefficient bounds (an
// edge with such a vertex), as log2 of its radius, and the roots on it: where the model's noise, not its roots,
// sets a scale. An edge between two resolved vertices is at its roots' own scale, which no precision moves.
struct UndeterminedReach {
    double radius_log2 = -std::numeric_limits<double>::infinity();
    slong roots = 0;
};
bool scale_candidates(acb_ptr output, const acb_poly_t polynomial, slong precision, bool &determined,
                      UndeterminedReach *reach = nullptr);
CandidateProgress
secular_candidate_packet(acb_ptr output, acb_srcptr initial, std::span<const ComplexBall> nodes,
                         std::span<const ComplexBall> residues, std::span<const uint8_t> inactive,
                         std::span<const slong> root_precision, unsigned sweeps, const cps_options &options,
                         certpolysolve::detail::CandidateStatistics &statistics, MonicSource *direct_source = nullptr,
                         bool tree = false, bool jacobi = true,
                         class SecularCandidateWorkspace *workspace = nullptr, uint64_t revision = 0,
                         const EnginePolicy *policy = nullptr, ModelStoppingState *model = nullptr);
// The owning refinement increments revision whenever interpolation nodes or
// residues change. Precision/source changes invalidate the midpoint models;
// proposals are copied from the caller on every packet, never cached as proofs.
class SecularCandidateWorkspace {
    struct Impl;
    std::unique_ptr<Impl> impl_;
  public:
    SecularCandidateWorkspace();
    ~SecularCandidateWorkspace();
    friend CandidateProgress secular_candidate_packet(acb_ptr, acb_srcptr,
        std::span<const ComplexBall>, std::span<const ComplexBall>, std::span<const uint8_t>,
        std::span<const slong>, unsigned, const cps_options &, certpolysolve::detail::CandidateStatistics &, MonicSource *,
        bool, bool, SecularCandidateWorkspace *, uint64_t, const EnginePolicy *, ModelStoppingState *);
};
// A prepared candidate workspace owns scaled coefficient midpoints, sparse
// support, proposals and regenerated secular nodes. Its input is immutable at
// a fixed precision; callers replace it on representation/precision changes.
class CandidateWorkspace {
    struct Impl;
    std::unique_ptr<Impl> impl_;
  public:
    explicit CandidateWorkspace(const EnginePolicy *policy = nullptr);
    ~CandidateWorkspace();
    // Only the local-model caller installs this synchronous borrowed hook.
    // A default workspace otherwise has the transient packet's exact policy.
    void research_set_local_packet_service(ResearchLocalPacketService, void *);
    friend CandidateProgress candidate_packet(acb_ptr,const acb_poly_t,acb_srcptr,slong,unsigned,
        CandidateEngine,const cps_options&,certpolysolve::detail::CandidateStatistics&,std::span<const uint8_t>,bool,CandidateWorkspace*);
};
// Heuristic proposals ONLY. Output points must be checked against the original
// coefficient enclosures by the root inclusion/coverage layer.
CandidateProgress candidate_packet(acb_ptr output, const acb_poly_t polynomial, acb_srcptr initial,
                                   slong precision, unsigned sweeps, CandidateEngine engine,
                                   const cps_options &options, certpolysolve::detail::CandidateStatistics &statistics,
                                   std::span<const uint8_t> inactive = {}, bool tree = false, CandidateWorkspace *workspace = nullptr);

// Interpolation identity: for leading coefficient c and distinct nodes b_i,
// a_i = -p(b_i)/(c*product_{j!=i}(b_i-b_j)) gives
// p(x)/c = product(x-b_i) * (1-sum(a_i/(x-b_i))).
// Outward balls enclose these exact residues; zero-containing differences fail.
bool regenerate_secular(std::vector<ComplexBall> &residues, const acb_poly_t polynomial,
                        std::span<const ComplexBall> nodes, slong precision, const cps_options &options,
                        certpolysolve::detail::CandidateStatistics &statistics);
} // namespace certpolysolve::detail
