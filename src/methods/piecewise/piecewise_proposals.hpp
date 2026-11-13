#pragma once
// The piecewise method's proposal engine, kept apart from piecewise_roots.cpp for readability (History:
// docs/code-history-20261021.md, piecewise_proposals.hpp 1): the binary64 iteration of a chart's local model that
// proposes root points (proposals), their polish at the working precision (polish), the chart-exclusion test and
// the seed count behind them, the guard index that retires candidate work around published roots, and the options
// under which a provisional ownership guard is refined (guard_options). The small rectangle and realness helpers
// use cover.hpp for exact rectangle geometry; the rounding/realness helpers and guard index remain inline here.
#include "core/cover.hpp"
#include "numeric/arithmetic/modulus_bound.hpp"
#include "core/raii.hpp"
#include "methods/piecewise/root_index.hpp"
#include <algorithm>
#include <array>
#include <cstdint>
#include <limits>
#include <map>
#include <span>
#include <utility>
#include <vector>

namespace certpolysolve::detail {
struct CandidateStatistics;
struct PiecewiseRootConfig;
struct PiecewiseRootStatistics;
struct RootModel;
inline bool strictly_nonreal(const Root &root) {
    return arb_is_positive(acb_imagref(root.box.value())) || arb_is_negative(acb_imagref(root.box.value()));
}
inline Root conjugate_proof(Root root) {
    acb_conj(root.box.value(), root.box.value());
    std::swap(root.endpoints[2], root.endpoints[3]);
    for (size_t k : {size_t(2), size_t(3)})
        root.endpoints[k].negate();
    root.certificate += "; conjugation of exact real coefficients";
    return root;
}

// Larger, independently counted rectangles retire candidate work around an
// already published root. Index each guard at its own spatial resolution;
// point queries visit only the occupied levels, even across large exponents.
class GuardIndex {
    std::vector<Root> guards_;
    std::vector<std::array<RealBall, 4>> endpoints_;
    std::vector<size_t> owners_;
    std::map<slong, RootIndex> levels_;
    static bool contained(const std::array<RealBall, 4> &bounds, const acb_t box) {
        for (size_t k = 0; k < 4; ++k) {
            const auto component = k < 2 ? acb_realref(box) : acb_imagref(box);
            if (k % 2 == 0 ? !arb_le(bounds[k].value(), component) : !arb_le(component, bounds[k].value()))
                return false;
        }
        return true;
    }

  public:
    bool contains(const acb_t point, size_t before_owner = std::numeric_limits<size_t>::max()) const {
        if (!acb_is_finite(point))
            return false;
        Root query;
        acb_get_mid(query.box.value(), point);
        for (const auto &[level, index] : levels_) {
            (void)level;
            if (index.overlaps(query, guards_, [&](size_t i) {
                    return owners_[i] < before_owner && contained(endpoints_[i], query.box.value());
                }))
                return true;
        }
        return false;
    }
  private:
    bool store(Root guard, std::array<RealBall, 4> bounds, size_t owner) {
        Float size, norm;
        RealBall width, other;
        arb_get_rad_arb(width.value(), acb_realref(guard.box.value()));
        arb_get_rad_arb(other.value(), acb_imagref(guard.box.value()));
        arb_add(width.value(), width.value(), other.value(), 64);
        arb_get_ubound_arf(size.value(), width.value(), 64);
        acb_get_abs_ubound_arf(norm.value(), guard.box.value(), 64);
        const slong resolution = arf_is_zero(size.value()) ? 48 : std::clamp<slong>(
            arf_abs_bound_lt_2exp_si(norm.value()) - arf_abs_bound_lt_2exp_si(size.value()), 4, 48);
        auto [where, inserted] = levels_.try_emplace(resolution, resolution);
        (void)inserted;
        where->second.insert(guard, guards_.size());
        guards_.push_back(std::move(guard));
        endpoints_.push_back(std::move(bounds));
        owners_.push_back(owner);
        return true;
    }
  public:
    bool insert(Root guard, const Root &published, size_t owner = 0) {
        if (guard.count != published.count || !acb_is_finite(guard.box.value()) ||
            !acb_is_finite(published.box.value())) return false;
        auto bounds = exact_endpoints(guard);
        if (!contains_rectangle(bounds, exact_endpoints(published))) return false;
        return store(std::move(guard), std::move(bounds), owner);
    }
    bool insert_point(Root guard, const ComplexBall &point, size_t owner) {
        // An unpublished seed is a point, not a counted rectangle. A singleton
        // guard can retire later guesses in its region after containing this
        // earlier seed. No default Root endpoints participate in that test.
        if (guard.count != 1 || !acb_is_finite(guard.box.value()) || !acb_is_finite(point.value())) return false;
        auto bounds = exact_endpoints(guard);
        if (!contained(bounds, point.value())) return false;
        return store(std::move(guard), std::move(bounds), owner);
    }
};

// The options of a provisional ownership guard's refinement: the accuracy goal (a guard contains one root).
cps_options guard_options(const RootModel &model, const cps_options &options);
// Root proposals from the model's binary64 iteration into `points`, the count returned; with `full_recovery`
// the interior proposals are skipped (the rationale is with the definition).
slong proposals(RootVector &points, const RootModel &model, const cps_options &options,
                certpolysolve::detail::CandidateStatistics &stats, slong precision, PiecewiseRootStatistics &observed,
                const PiecewiseRootConfig &config, bool full_recovery = false);
// The proposals polished at `precision`; `stable` marks the points that settled, `owners` and `owner_counts`
// are the ownership diagnostic (the rationale is with the definition).
void polish(RootVector &points, slong count, const RootModel &model, const cps_options &options,
            slong precision, PiecewiseRootStatistics &stats, std::vector<bool> *stable, const PiecewiseRootConfig &config,
            const std::vector<size_t> *owners = nullptr, std::span<const uint64_t> owner_counts = {});
} // namespace certpolysolve::detail
