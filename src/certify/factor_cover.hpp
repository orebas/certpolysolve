#pragma once
// The cover assembly both coefficient pipelines share (W8's S8, the map's section 4.2 and commit A1): the exact zero
// point, the ledger with the merge of meeting factor covers, and the publication rule. Since C1 to C4 and B1 the
// global driver (SolveState) and piecewise (PreparedPiecewiseState) assemble by one rule set, so the assembler takes
// no rules: each caller passes its own precision.
//
// State: none. Entry points: exact_zero_point (SolveState::Impl::aggregate in solver.cpp; piecewise's factor loop and
// its certificate attempt in prepared_piecewise.cpp), assemble_factor_covers (aggregate, the factor loop),
// publish_assembly (aggregate, the factor loop, native_advance in solver_native.cpp).
// Not here: the lifts of power factors (lift_power_cover, power_composition.cpp), the merge itself and the ranking of
// covers (merge_overlapping and prefer_cover, certification.cpp), the tightening of a factor's goal (separation.hpp).
#include "core/cover.hpp"
#include "core/query.hpp"
#include <vector>

namespace certpolysolve::detail {
// The exact zero roots set apart before the factors (B1): the point 0, counted `order` times with multiplicity
// `order`, real, published at 64 bits.
Root exact_zero_point(uint64_t order);
// One complete assembly of coprime factors' covers, the zero point among them: the counts must sum to `expected`
// (the degree, or a real-only query's count of real roots), regions that meet are merged into counted clusters
// (merge_overlapping at `precision`), and the result must be disjoint. A failure is the caller's logic error
// (std::logic_error).
void assemble_factor_covers(std::vector<Root> &cover, uint64_t expected, slong precision);
// Publication: the assembly replaces the answer when it meets the goal (`attained`), when no complete answer stands,
// or when prefer_cover ranks it ahead at `precision`; otherwise the better answer stays. Returns whether it replaced
// the answer; the result is complete either way.
bool publish_assembly(cps_result &result, std::vector<Root> &assembly, bool attained, slong precision);
} // namespace certpolysolve::detail
