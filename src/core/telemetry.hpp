#pragma once
#include "core/heuristics.hpp"
#include "core/engine_policy.hpp"
#include <chrono>

namespace certpolysolve::detail {
// Internal observer only. The observer must not throw, re-enter or mutate the
// solve. Entry events survive external process termination during atomic calls.
struct TraceScope {
    const EnginePolicy* policy;
    TraceEvent event;
    TraceScope(const EnginePolicy* p, const char* stage, uint64_t precision, uint64_t active,
               const certpolysolve::detail::CandidateStatistics* candidates = nullptr, const certpolysolve::detail::SearchStatistics* search = nullptr)
        : policy(p), event{stage, "begin", precision, active, candidates, search} {
        if (policy) {
            event.factor = policy->snapshot_factor;
            event.factor_power = policy->snapshot_power;
            event.centered = !fmpq_is_zero(policy->snapshot_offset.real.value()) ||
                             !fmpq_is_zero(policy->snapshot_offset.imag.value());
        }
        if (policy && policy->trace)
            policy->trace(policy->trace_data, event);
    }
    ~TraceScope() {
        if (policy && policy->trace) {
            event.phase = "end";
            policy->trace(policy->trace_data, event);
        }
    }
    void progress(uint64_t precision, uint64_t active) {
        event.precision = precision;
        event.active = active;
    }
};
// Eligibility diagnostics are deliberately labelled separately from execution.
// A disabled optional helper does not imply that the caller attempted it.
inline void trace_degree_policy(const EnginePolicy &policy, slong degree) {
    if (!policy.trace) return;
    const auto emit = [&](const char *name, slong cutoff, const char *fallback) {
        TraceEvent event{"algorithm_degree_policy", "eligibility", 0,
                         resources::nonnegative(degree), nullptr, nullptr};
        event.engine = name;
        event.algorithm_cutoff = resources::nonnegative(cutoff);
        event.reason = degree > cutoff ? fallback : "within optional helper crossover";
        policy.trace(policy.trace_data, event);
    };
    emit("automatic_secular_reconstruction", heuristics::secular_reconstruction_degree,
         "degree crossover: ordinary coefficient/direct route remains available");
    emit("secular_block_seeds", heuristics::secular_block_seed_degree,
         "degree crossover: ordinary secular initialization remains available");
#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
    // The research build's DoubleSeed test engine's start.
    emit("legacy_binary64_seeds", heuristics::binary64_seed_degree,
         "degree crossover: ordinary candidates remain available; modern machine factory unaffected");
#endif
    emit("candidate_cluster_recovery", heuristics::cluster_recovery_degree,
         "degree crossover: retained count cover and subdivision remain available");
}
struct StageTimer {
    uint64_t &elapsed;
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    explicit StageTimer(uint64_t &target) : elapsed(target) {}
    ~StageTimer() {
        elapsed += static_cast<uint64_t>(
            std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - start)
                .count());
    }
};
// Fine-grained private measurements can be disabled without reading the clock.
// Keep this separate from the always-on aggregate stage counters.
struct OptionalStageTimer {
    uint64_t *elapsed;
    std::chrono::steady_clock::time_point start;
    OptionalStageTimer(uint64_t &target, bool enabled)
        : elapsed(enabled ? &target : nullptr) {
        if (elapsed) start = std::chrono::steady_clock::now();
    }
    ~OptionalStageTimer() {
        if (elapsed)
            *elapsed += static_cast<uint64_t>(std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::steady_clock::now() - start).count());
    }
};
struct ExclusiveTimer {
    uint64_t &elapsed;
    const uint64_t *nested[4];
    uint64_t before;
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    uint64_t total() const {
        uint64_t sum = 0;
        for (auto p : nested)
            if (p)
                sum += *p;
        return sum;
    }
    ExclusiveTimer(uint64_t &target, const uint64_t *a = nullptr, const uint64_t *b = nullptr,
                   const uint64_t *c = nullptr, const uint64_t *d = nullptr)
        : elapsed(target), nested{a, b, c, d}, before(total()) {}
    ~ExclusiveTimer() {
        elapsed += static_cast<uint64_t>(std::chrono::duration_cast<std::chrono::nanoseconds>(
                                             std::chrono::steady_clock::now() - start)
                                             .count()) -
                   (total() - before);
    }
};
} // namespace certpolysolve::detail
