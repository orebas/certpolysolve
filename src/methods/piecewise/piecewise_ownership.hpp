#pragma once
#include "core/cover.hpp"
#include <limits>
#include <memory>
#include <optional>

namespace certpolysolve::detail {
struct RootModel;
// Original-polynomial proofs only. Candidates cannot enter this ledger.
// IDs are never recycled, including when a counted cluster absorbs owners.
class PiecewiseOwnership {
    struct Impl;
    std::unique_ptr<Impl> impl_;
  public:
    using Id = size_t;
    static constexpr Id none = std::numeric_limits<Id>::max();
    explicit PiecewiseOwnership(uint64_t degree);
    ~PiecewiseOwnership();
    std::optional<Id> admit(const Root &proof);
    bool active(Id) const;
    const Root &proof(Id) const;
    bool contains(Id, const acb_t point) const;
    std::optional<Id> containing(const acb_t point) const;
    bool identifies(Id, const Root &refinement) const;
    std::optional<Id> identifying(const Root &refinement) const;
    uint64_t count() const;
    bool singletons() const;
    // Counts only full owned rectangles contained in this chart's disk.
    // The caller must separately establish the chart's original-input count.
    uint64_t contained_in(const RootModel &, slong precision) const;
};
} // namespace certpolysolve::detail
