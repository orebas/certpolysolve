// The cover assembly both coefficient pipelines share (factor_cover.hpp): the exact zero point, the ledger and merge,
// and the publication rule.
//
// State: none. Entry points: factor_cover.hpp's three functions. Not here: the callers' loops (solver.cpp's aggregate,
// prepared_piecewise.cpp's factor loop, solver_native.cpp's native_advance). Design notes: docs/mathematics.md
// ("Individual root certificates").
#include "certify/factor_cover.hpp"
#include <stdexcept>

namespace certpolysolve::detail {
Root exact_zero_point(uint64_t order) {
    Root zero;
    zero.count = zero.multiplicity = order;
    zero.realness = CPS_REAL;
    zero.certificate = "exact zero monomial factor before square-free preparation";
    zero.publish(64);
    return zero;
}

void assemble_factor_covers(std::vector<Root> &cover, uint64_t expected, slong precision) {
    uint64_t count = 0;
    for (const auto &root : cover) count += root.count;
    if (count != expected) throw std::logic_error("factor coverage count ledger disagrees with degree");
    // Square-free factors are coprime: merging meeting regions until the cover is disjoint sums their known roots.
    if (!merge_overlapping(cover, precision) && !disjoint(cover))
        throw std::logic_error("outward cover merging failed");
}

bool publish_assembly(cps_result &result, std::vector<Root> &assembly, bool attained, slong precision) {
    const bool published = attained || !result.complete || prefer_cover(assembly, result.roots, precision);
    if (published) result.roots = std::move(assembly);
    result.complete = true;
    return published;
}
} // namespace certpolysolve::detail
