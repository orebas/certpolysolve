#pragma once
#include "core/cover.hpp"
#include <optional>

namespace certpolysolve::detail {
// Local instrumentation for the focused growth tests, not solver work counters.
struct CoverSweepWork {
    uint64_t comparisons = 0, marked = 0;
    // Ranking's comparison charge counts the shorter mantissa, at least one
    // bit. Exponent encodings and fast early exits are not measured here.
    uint64_t rank_comparisons = 0, rank_mantissa_bits = 0, endpoint_mantissa_bits = 0;
};

// Owns order ranks of the actual Arb boxes' exact endpoints, not Root's
// exported endpoints. Equal values have equal ranks, preserving closed contact.
// Preparation refuses exponent expansion beyond a bound derived from the
// generating precision. Published native certificates fit that bound; arbitrary
// compact-exponent boxes may need the existing nonexpanding predicate instead.
class CoverSweep {
    struct Rectangle {
        std::array<size_t, 2> x;
        std::array<std::array<size_t, 2>, 2> y; // original and reflected
    };
    std::vector<Rectangle> rectangles_;
    explicit CoverSweep(size_t size) : rectangles_(size) {}
    struct Event;
    struct IntervalLess;
    std::vector<Event> events(bool reflected, CoverSweepWork *work) const;

  public:
    static std::optional<CoverSweep> prepare(const std::vector<Root> &, slong generating_precision,
                                             CoverSweepWork *work = nullptr);
    bool disjoint(CoverSweepWork *work = nullptr) const;
    // Precondition: disjoint() succeeded. Each active family's y intervals
    // are then disjoint. Self/reflection intersections do not mark a record.
    std::vector<bool> reflection_met(CoverSweepWork *work = nullptr) const;
};
} // namespace certpolysolve::detail
