#pragma once
#include "core/cover.hpp"
#include "core/query.hpp"

namespace certpolysolve::detail {
// Native Descartes isolation followed by interval Newton/bisection. Each step
// handles one retained node or one isolated interval; exact arithmetic calls
// remain atomic. cover() is false until every isolation node is accounted for.
class RealRefinement {
    struct Impl;
    std::unique_ptr<Impl> impl_;

  public:
    RealRefinement(const fmpz_poly_t polynomial, uint64_t multiplicity, slong precision,
                   certpolysolve::detail::SearchStatistics &statistics, certpolysolve::detail::CandidateStatistics &extended);
    ~RealRefinement();
    cps_stop_reason step(const cps_options &options);
    bool cover(std::vector<Root> &roots) const;
    uint64_t count() const;
    uint64_t precision() const;
    uint64_t cells() const;
};
} // namespace certpolysolve::detail
