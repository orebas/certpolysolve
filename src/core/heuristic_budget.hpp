#pragma once
#include "core/resource_limits.hpp"
#include <algorithm>
#include <cstdint>

namespace certpolysolve::detail {
// One nonrenewing allowance owned by the outer coefficient query. Synchronous
// local solves borrow it through their internal options. It is independent of
// caller cancellation and of the calibrated prices used for practical racing.
// A tick is 1/n of a full-size model packet: a subdivision action costs at
// least one tick, an outer numerical action or full-source local attempt n.
// Finish the admitted operation, then inspect the boundary. No rollback or
// interruption inside arithmetic is required. K is fixed in the model bound;
// frontier scans can contribute a K² constant (see the guarantee document).
class HeuristicBudget {
  public:
    enum class Stop { None, Packets, Precision };
    // Common allowance for both accuracy contracts. Calibration of the first
    // bounded candidate exposed useful source/model work up to about 85D,
    // D=max(64,n+tau+b), and packet cutoffs that interrupted the warm search.
    // 128D leaves precision headroom; 4096 is four times the initial packet
    // allowance. These fixed constants preserve the modeled exponents, but
    // must still qualify against caller limits and complete-service latency.
    static constexpr uint64_t default_packets = 4096;
    static constexpr uint64_t default_precision_multiplier = 128;

    HeuristicBudget(uint64_t degree, uint64_t height, uint64_t bits,
                    uint64_t packets = default_packets,
                    uint64_t precision_multiplier = default_precision_multiplier)
      : degree_(std::max<uint64_t>(1, degree)),
        limit_(resources::multiply(packets, degree_)),
        ceiling_(std::max<uint64_t>(64, std::min<uint64_t>(resources::precision_max,
            resources::multiply(precision_multiplier,
                std::max<uint64_t>(64, resources::add(resources::add(degree, height), bits)))))),
        effective_bits_(resources::add(ceiling_, resources::add(resources::add(degree, height), 1))) {}

    void quadratic() { charge(degree_); }
    void linear() { charge(1); }
    void observe_precision(uint64_t bits) {
        peak_ = std::max(peak_, bits);
    }
    // Q limits each admitted operation. An optional child may touch that
    // ceiling and return useful partial work without ending its parent's
    // allowance. Only the owning controller's terminal precision stop ends
    // the phase; the nonrenewing packet limit still bounds all work at Q.
    void finish_at_precision_limit() {
        if (stop_ == Stop::None) stop_ = Stop::Precision;
    }
    bool exhausted() const { return stop_ != Stop::None || used_ >= limit_; }
    uint64_t ceiling() const { return ceiling_; }
    uint64_t degree() const { return degree_; }
    uint64_t ticks() const { return used_; }
    uint64_t limit() const { return limit_; }
    uint64_t peak() const { return peak_; }
    // Existing exact-preparation admission uses bit/operation estimates. One
    // optional translation may consume a quadratic packet, including the
    // degree-dependent output growth; it cannot materialize an unlimited exact
    // representation merely because floating arithmetic is capped at Q.
    uint64_t exact_limit() const {
        return resources::multiply(resources::multiply(degree_, degree_), effective_bits_);
    }
    Stop stop() const { return stop_ == Stop::None && used_ >= limit_ ? Stop::Packets : stop_; }
  private:
    void charge(uint64_t ticks) { used_ = resources::add(used_, ticks); }
    uint64_t degree_, limit_, ceiling_, effective_bits_, used_ = 0, peak_ = 0;
    Stop stop_ = Stop::None;
};
} // namespace certpolysolve::detail
