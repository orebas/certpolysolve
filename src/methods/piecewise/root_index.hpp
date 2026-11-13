#pragma once
#include <functional>
#include "core/cover.hpp"
#include "core/wide_int.hpp"
#include <map>
#include <tuple>

namespace certpolysolve::detail {
// An exponent-normalized grid indexes complete rectangle bounds, not only
// centers. Narrow nonzero boxes can overlap only in neighboring exponent bins.
// Wide boxes use a conservative fallback; the scaling contract uses >=53 bits.
class RootIndex {
    using Key = std::tuple<slong, slong, slong>;
    std::map<Key, std::vector<size_t>> cells_;
    std::vector<size_t> wide_;
    slong resolution_;
    bool grid(const Root &r, slong exponent, std::array<slong, 4> &bounds, bool *representable = nullptr) const {
        if (representable) *representable = false;
        slong scale;
        if (subtract_overflows(resolution_, exponent, &scale)) return false;
        Float lo, hi;
        for (int axis = 0; axis < 2; ++axis) {
            arb_get_interval_arf(lo.value(), hi.value(),
                                 axis ? acb_imagref(r.box.value()) : acb_realref(r.box.value()), 64);
            arf_mul_2exp_si(lo.value(), lo.value(), scale);
            arf_mul_2exp_si(hi.value(), hi.value(), scale);
            bounds[static_cast<size_t>(2 * axis)] = arf_get_si(lo.value(), ARF_RND_FLOOR);
            bounds[static_cast<size_t>(2 * axis + 1)] = arf_get_si(hi.value(), ARF_RND_FLOOR);
        }
        if (representable) *representable = true;
        return bounds[1] - bounds[0] <= 2 && bounds[3] - bounds[2] <= 2;
    }
    static bool bin(const Root &r, slong &exponent) {
        Float lo, hi;
        acb_get_abs_lbound_arf(lo.value(), r.box.value(), 64);
        acb_get_abs_ubound_arf(hi.value(), r.box.value(), 64);
        if (!arf_is_finite(hi.value()) || arf_sgn(lo.value()) <= 0)
            return false;
        exponent = arf_abs_bound_lt_2exp_si(hi.value());
        // The neighbor loop visits e-1 through e+1 and then increments once.
        if (exponent <= WORD_MIN + 1 || exponent >= WORD_MAX - 1) return false;
        arf_mul_2exp_si(hi.value(), hi.value(), -1);
        return arf_cmp(lo.value(), hi.value()) > 0;
    }

  public:
    explicit RootIndex(slong resolution = 24) : resolution_(resolution) {
        if (resolution < 0 || resolution > FLINT_BITS - 3)
            throw std::invalid_argument("spatial index resolution exceeds signed grid coordinates");
    }
    bool overlaps(const Root &r, const std::vector<Root> &roots,
                  const std::function<bool(size_t)> &accept = {},
                  const std::vector<size_t> *members = nullptr) const {
        auto hit = [&](size_t i) {
            return acb_overlaps(r.box.value(), roots[i].box.value()) && (!accept || accept(i));
        };
        auto scan = [&] {
            if (members) {
                for (auto i : *members) if (hit(i)) return true;
            } else {
                for (size_t i = 0; i < roots.size(); ++i) if (hit(i)) return true;
            }
            return false;
        };
        slong e;
        if (!bin(r, e)) return scan();
        std::array<slong, 4> finest;
        bool representable;
        grid(r, e - 1, finest, &representable);
        if (!representable) return scan();
        // floor(floor(x)/2) = floor(x/2), including negative coordinates.
        // The outward endpoints are the same at all three scales. Extract
        // them once at the finest level, then coarsen integer coordinates.
        auto coarsen = [](std::array<slong, 4> bounds) {
            for (auto &x : bounds) x = x / 2 - (x < 0 && x % 2 != 0);
            return bounds;
        };
        const auto base = coarsen(finest);
        if (base[1] - base[0] > 2 || base[3] - base[2] > 2) return scan();
        for (auto i : wide_)
            if (hit(i))
                return true;
        auto b = finest;
        for (slong level = e - 1; level <= e + 1; ++level) {
            for (slong x = b[0]; x <= b[1]; ++x)
                for (slong y = b[2]; y <= b[3]; ++y) {
                    const auto found = cells_.find({level, x, y});
                    if (found != cells_.end())
                        for (auto i : found->second)
                            if (hit(i))
                                return true;
                }
            b = coarsen(b);
        }
        return false;
    }
    void insert(const Root &r, size_t index) {
        slong e;
        std::array<slong, 4> b;
        if (!bin(r, e) || !grid(r, e, b)) {
            wide_.push_back(index);
            return;
        }
        for (slong x = b[0]; x <= b[1]; ++x)
            for (slong y = b[2]; y <= b[3]; ++y)
                cells_[{e, x, y}].push_back(index);
    }
};

} // namespace certpolysolve::detail
