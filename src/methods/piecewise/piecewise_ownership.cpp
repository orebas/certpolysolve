// Stable proof identities and geometric ownership, independent of the piecewise controller.
#include "methods/piecewise/piecewise_ownership.hpp"
#include "core/cover.hpp"
#include "numeric/arithmetic/modulus_bound.hpp"
#include "methods/piecewise/root_index.hpp"
#include "numeric/models/root_model.hpp"
#include <algorithm>
#include <array>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

namespace certpolysolve::detail {

struct PiecewiseOwnership::Impl {
    using Bounds = std::array<RealBall, 4>;
    struct Bucket {
        RootIndex index;
        std::vector<size_t> members;
        explicit Bucket(slong resolution) : index(resolution) {}
    };
    uint64_t degree, total = 0, clusters = 0;
    std::vector<Root> proofs;
    std::vector<Bounds> bounds;
    std::vector<bool> live;
    std::map<slong, Bucket> buckets;
    explicit Impl(uint64_t n) : degree(n) {}
    static Bounds endpoints(const Root &root) {
        return exact_endpoints(root);
    }
    static bool contains(const Bounds &bounds, const acb_t box) {
        if (!acb_is_finite(box)) return false;
        for (size_t k = 0; k < 4; ++k) {
            const auto component = k < 2 ? acb_realref(box) : acb_imagref(box);
            if (k % 2 == 0 ? !arb_le(bounds[k].value(), component) : !arb_le(component, bounds[k].value()))
                return false;
        }
        return true;
    }
    std::vector<size_t> overlapping(const Root &query) const {
        std::vector<size_t> result;
        for (const auto &[resolution, bucket] : buckets) {
            (void)resolution;
            bucket.index.overlaps(query, proofs, [&](size_t i) {
                if (live[i]) result.push_back(i);
                return false; // Enumerate; a grid cell can mention an ID twice.
            }, &bucket.members);
        }
        std::sort(result.begin(), result.end());
        result.erase(std::unique(result.begin(), result.end()), result.end());
        return result;
    }
    void index(const Root &root, size_t id) {
        Float norm, radius;
        RealBall width, other;
        arb_get_rad_arb(width.value(), acb_realref(root.box.value()));
        arb_get_rad_arb(other.value(), acb_imagref(root.box.value()));
        arb_add(width.value(), width.value(), other.value(), 64);
        arb_get_ubound_arf(radius.value(), width.value(), 64);
        acb_get_abs_ubound_arf(norm.value(), root.box.value(), 64);
        const slong resolution = arf_is_zero(radius.value()) ? 48 : std::clamp<slong>(
            arf_abs_bound_lt_2exp_si(norm.value()) - arf_abs_bound_lt_2exp_si(radius.value()), 4, 48);
        auto [where, inserted] = buckets.try_emplace(resolution, resolution);
        (void)inserted;
        where->second.index.insert(root, id);
        where->second.members.push_back(id);
    }
};
PiecewiseOwnership::PiecewiseOwnership(uint64_t degree) : impl_(std::make_unique<Impl>(degree)) {}
PiecewiseOwnership::~PiecewiseOwnership() = default;
bool PiecewiseOwnership::active(Id id) const { return id < impl_->live.size() && impl_->live[id]; }
const Root &PiecewiseOwnership::proof(Id id) const {
    if (!active(id)) throw std::logic_error("inactive piecewise owner");
    return impl_->proofs[id];
}
bool PiecewiseOwnership::contains(Id id, const acb_t point) const {
    return active(id) && Impl::contains(impl_->bounds[id], point);
}
std::optional<PiecewiseOwnership::Id> PiecewiseOwnership::containing(const acb_t point) const {
    if (!acb_is_finite(point)) return {};
    Root query;
    acb_set(query.box.value(), point);
    for (auto id : impl_->overlapping(query)) if (contains(id, point)) return id;
    return {};
}
bool PiecewiseOwnership::identifies(Id id, const Root &refinement) const {
    if (!active(id) || refinement.count != impl_->proofs[id].count ||
        !acb_is_finite(refinement.box.value())) return false;
    const auto bounds = Impl::endpoints(refinement);
    return contains_rectangle(impl_->bounds[id], bounds) ||
           contains_rectangle(bounds, impl_->bounds[id]);
}
std::optional<PiecewiseOwnership::Id> PiecewiseOwnership::identifying(const Root &refinement) const {
    if (!acb_is_finite(refinement.box.value())) return {};
    ComplexBall midpoint;
    acb_get_mid(midpoint.value(), refinement.box.value());
    const auto nearby = containing(midpoint.value());
    if (nearby && identifies(*nearby, refinement)) return nearby;
    // A valid larger proof can contain its owner while its midpoint lies
    // outside that owner's narrower rectangle. Check proof nesting itself.
    for (const auto id : impl_->overlapping(refinement))
        if (identifies(id, refinement)) return id;
    return {};
}
std::optional<PiecewiseOwnership::Id> PiecewiseOwnership::admit(const Root &proof) {
    if (!proof.count || !acb_is_finite(proof.box.value())) return {};
    auto bounds = Impl::endpoints(proof);
    const auto overlaps = impl_->overlapping(proof);
    if (overlaps.size() == 1 && identifies(overlaps[0], proof)) return overlaps[0];
    uint64_t replaced = 0;
    for (auto id : overlaps) {
        // A cluster may absorb owners only by enclosing every old proof.
        // Partial geometric overlap makes identity inconclusive; mutate nothing.
        if (!contains_rectangle(bounds, impl_->bounds[id])) return {};
        replaced += impl_->proofs[id].count;
    }
    if (proof.count < replaced) return {};
    if (proof.count - replaced > impl_->degree - impl_->total)
        throw std::logic_error("piecewise ownership exceeds factor degree");
    const uint64_t next_total = impl_->total - replaced + proof.count;
    const Id id = impl_->proofs.size();
    impl_->proofs.push_back(proof);
    impl_->bounds.push_back(std::move(bounds));
    impl_->live.push_back(true);
    impl_->index(proof, id);
    // No old proof is destroyed. Commit the active count cover atomically after
    // all containment/count checks and storage operations have succeeded.
    for (auto old : overlaps) {
        impl_->live[old] = false;
        impl_->clusters -= impl_->proofs[old].count > 1;
    }
    impl_->clusters += proof.count > 1;
    impl_->total = next_total;
    return id;
}
uint64_t PiecewiseOwnership::count() const { return impl_->total; }
bool PiecewiseOwnership::singletons() const { return impl_->clusters == 0; }
uint64_t PiecewiseOwnership::contained_in(const RootModel &model, slong precision) const {
    Root query;
    query.box = model.center;
    arb_add_error(acb_realref(query.box.value()), model.radius.value());
    arb_add_error(acb_imagref(query.box.value()), model.radius.value());
    uint64_t total = 0;
    for (auto id : impl_->overlapping(query)) {
        ComplexBall difference;
        RealBall distance;
        acb_sub(difference.value(), impl_->proofs[id].box.value(), model.center.value(), precision);
        upper(distance.value(), difference.value());
        if (arb_lt(distance.value(), model.radius.value())) total += impl_->proofs[id].count;
    }
    return total;
}

} // namespace certpolysolve::detail
