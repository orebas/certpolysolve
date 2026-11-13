#pragma once
// The real-only route of a session (session.cpp hands it every real-only query under Auto's selector on exact real
// polynomial input). A real-only query is answered by a solve on the complex route, within the same
// budgets: the real roots are its individual roots proved real by reflection. The complex route finishes more of
// the real inputs than the old default route's real engine (exact Descartes isolation), in less time. History:
// docs/code-history-20261021.md, real_only_route.hpp 1. A real-only query stops where the complex route stops at
// its precision cap, and the stop names a cap that suffices (solver.cpp, precision_cap_note): exact isolation past
// the cap by the real engine, tried as a second attempt, finished none of the corpus inputs that reached it.
// History: docs/code-history-20261021.md, real_only_route.hpp 2. The real engine answers no public real-only query
// on a polynomial: such a query reaches it only on input it refuses (inexact or not real). It answers a native
// secular equation's real-only query (on the exact numerator), and D5 routing (its tests) still reaches it.
//
// Before the complex route, two attempts, each dropped after its stop, its work and statistics carried:
// - on an input with few terms, the sparse recursion (B6, docs/b6-sparse-real-roots-note-20261015.md,
//   sparse_real_attempt.hpp), where its predicted cost is below the estimate of the route the selector would take
//   (heuristics::sparse_recursion_multiplications): certified signs along the chain of fractional derivatives, in
//   time counted in terms; it stops honestly where consecutive derivatives may share a root, at the cap, or over
//   its work ceiling. History: docs/code-history-20261021.md, real_only_route.hpp 3.
// - the axis charts alone (piecewise's axis mode, at most heuristics::axis_attempt_rounds rounds), under an ample
//   work budget and where the selector sends the complex query to the piecewise method: it answers the real-only
//   inputs whose axis its first models decide (mostly of high degree with few real roots, in less time than the
//   complex route) and otherwise stops early, leaving the query to the complex route within what it left of the
//   call's budget. History: docs/code-history-20261021.md, real_only_route.hpp 4.
//
// The contracts of D1 (docs/real-interval-contract-note-20261015.md) live here too: the isolation goal (no width
// asked; the nested query at the library's default goal) and a closed region (real_region.hpp: the roots in it,
// each decided exactly against an end it reaches, the route solved again at doubled accuracy until every root lies
// on one side; each such solve is a fresh one). History: docs/code-history-20261021.md, real_only_route.hpp 5.
#include "input/coefficients/polynomial_input.hpp"
#include "core/query.hpp"
#include "certify/real_region.hpp"
#include "core/session_statistics.hpp"
#include "methods/real/sparse_real_attempt.hpp"
#include <memory>
#include <string>

namespace certpolysolve::detail {

class RealOnlyRoute {
  public:
    // Decides the attempts for the query and keeps `input`, from which every solve of a route is built afresh.
    // `degree` and `zero_input` are the input's, as the session recorded them.
    RealOnlyRoute(DensePolynomial input, const cps_options &first, uint64_t degree, bool zero_input,
                  std::unique_ptr<RealRegionFilter> region);
    // The same on the input object's terms, borrowed for the call (SessionState's borrow, C40's second step): the
    // sparse recursion reads the terms, and the coefficient list is formed only when the axis charts or the complex
    // route are reached (or the selector needs it to decide the attempts).
    RealOnlyRoute(const cps_polynomial &borrowed, const cps_options &first, uint64_t degree, bool zero_input,
                  std::unique_ptr<RealRegionFilter> region);
    ~RealOnlyRoute();
    RealOnlyRoute(const RealOnlyRoute &) = delete;
    RealOnlyRoute &operator=(const RealOnlyRoute &) = delete;

    // One call of the session: the attempts in order, then the complex route, each answer filtered by the region.
    cps_result advance(const cps_options &next);

  private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace certpolysolve::detail
