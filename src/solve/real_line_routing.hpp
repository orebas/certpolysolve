#pragma once
// Where the all-real seed pass (real_line_seeds.hpp) meets the routes: its limits under a call's options, shared by
// a factor's start (solver.cpp) and Auto's check before routing (session.cpp), and that check
// (docs/allroots-seed-routing-note-20261020.md: W1, step (b)).
#include "core/query.hpp"
#include "solve/solve_state.hpp"
#include "numeric/candidates/real_line_seeds.hpp"
#include <optional>
#include <span>

namespace certpolysolve::detail {
// The pass's limits under a call's options: its precision cap, its storage, the cell share, the call's cancellation
// (the limits borrow `o`, which must outlive them).
RealLineSeedLimits real_line_seed_limits(const cps_options &o);

// The pass's outcome on an input's primitive integer polynomial, computed before routing. The factor whose integer
// polynomial equals `polynomial` takes these seeds instead of running the pass again (SolveState::
// hint_real_line_seeds); any other factor runs it on its own polynomial.
struct RealLineSeedHint {
    IntegerPolynomial polynomial;
    RealLineSeeds pass;
};

// The pass on the primitive integer polynomial of exact real coefficients, whatever it found: nothing only when the
// import is refused under the options' exact budgets or the degree is below the pass's floor
// (heuristics::real_line_seed_min_degree). Auto's routing (session.cpp) reads it before choosing between piecewise
// and the secular route, and keeps it for piecewise's hand-off to secular, whose start would otherwise run the same
// pass on the same polynomial again (dev/docs/slowdowns-20261029.md).
std::optional<RealLineSeedHint> real_line_check(std::span<const Coefficient> coefficients, const cps_options &o);

// Whether the pass seeded every root of its polynomial: Auto's condition for the secular route.
bool seeds_every_root(const RealLineSeedHint &hint);

// real_line_check for Auto's query for all roots that the selector would send to piecewise: a hint only when the pass
// seeds every root, which no input whose roots are not all real did in its screening (the note's tables). A wrong
// hint would cost time only: seeds are starting points, and every answer keeps its certificate.
std::optional<RealLineSeedHint> every_root_seeded(std::span<const Coefficient> coefficients, const cps_options &o);
} // namespace certpolysolve::detail
