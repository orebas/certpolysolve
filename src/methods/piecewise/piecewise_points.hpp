#pragma once
#include "numeric/models/root_model.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace certpolysolve::detail {
// Candidate-only spatial index. Outward chart boxes prune a static kd tree;
// actual original-input certificates never depend on floating-point pruning.
class PiecewisePointIndex {
    struct Node {
        size_t point = 0;
        slong left = -1, right = -1;
        std::array<double, 2> lo{}, hi{};
    };
    std::vector<Node> nodes_;
    std::vector<std::array<double, 2>> coordinates_;
    slong root_ = -1, scale_ = 0;
    slong build(std::vector<size_t> &order, size_t lo, size_t hi, unsigned depth) {
        if (lo == hi)
            return -1;
        const size_t mid = lo + (hi - lo) / 2, axis = depth % 2;
        std::nth_element(order.begin() + static_cast<ptrdiff_t>(lo),
                         order.begin() + static_cast<ptrdiff_t>(mid),
                         order.begin() + static_cast<ptrdiff_t>(hi),
                         [&](size_t a, size_t b) { return coordinates_[a][axis] < coordinates_[b][axis]; });
        Node &node = nodes_[mid];
        node.point = order[mid];
        node.lo = node.hi = coordinates_[node.point];
        node.left = build(order, lo, mid, depth + 1);
        node.right = build(order, mid + 1, hi, depth + 1);
        for (auto child : {node.left, node.right})
            if (child >= 0)
                for (size_t k = 0; k < 2; ++k) {
                    node.lo[k] = std::min(node.lo[k], nodes_[static_cast<size_t>(child)].lo[k]);
                    node.hi[k] = std::max(node.hi[k], nodes_[static_cast<size_t>(child)].hi[k]);
                }
        return static_cast<slong>(mid);
    }
    void query_node(slong i, const std::array<double, 4> &box, std::vector<size_t> &out) const {
        if (i < 0)
            return;
        const auto &node = nodes_[static_cast<size_t>(i)];
        if (node.hi[0] < box[0] || node.lo[0] > box[1] || node.hi[1] < box[2] || node.lo[1] > box[3])
            return;
        const auto &p = coordinates_[node.point];
        if (p[0] >= box[0] && p[0] <= box[1] && p[1] >= box[2] && p[1] <= box[3])
            out.push_back(node.point);
        query_node(node.left, box, out);
        query_node(node.right, box, out);
    }

  public:
    explicit PiecewisePointIndex(const std::vector<ComplexBall> &points)
        : nodes_(points.size()), coordinates_(points.size()) {
        Float bound, component;
        for (const auto &p : points) {
            if (!acb_is_finite(p.value()))
                throw std::invalid_argument("nonfinite supplied point");
            acb_get_abs_ubound_arf(bound.value(), p.value(), 64);
            if (!arf_is_zero(bound.value()))
                scale_ = std::max(scale_, arf_abs_bound_lt_2exp_si(bound.value()));
        }
        for (size_t i = 0; i < points.size(); ++i)
            for (size_t k = 0; k < 2; ++k) {
                arf_mul_2exp_si(component.value(),
                                arb_midref(k ? acb_imagref(points[i].value()) : acb_realref(points[i].value())),
                                -scale_);
                coordinates_[i][k] = arf_get_d(component.value(), ARF_RND_NEAR);
            }
        std::vector<size_t> order(points.size());
        std::iota(order.begin(), order.end(), 0);
        root_ = build(order, 0, order.size(), 0);
    }
    std::vector<size_t> query(const RootModel &model) const {
        std::array<double, 4> box;
        RealBall component;
        Float lo, hi;
        for (size_t k = 0; k < 2; ++k) {
            arb_set(component.value(), k ? acb_imagref(model.center.value()) : acb_realref(model.center.value()));
            arb_add_error(component.value(), model.radius.value());
            arb_mul_2exp_si(component.value(), component.value(), -scale_);
            arb_get_interval_arf(lo.value(), hi.value(), component.value(), 64);
            box[2 * k] =
                std::nextafter(arf_get_d(lo.value(), ARF_RND_FLOOR), -std::numeric_limits<double>::infinity());
            box[2 * k + 1] =
                std::nextafter(arf_get_d(hi.value(), ARF_RND_CEIL), std::numeric_limits<double>::infinity());
        }
        std::vector<size_t> result;
        query_node(root_, box, result);
        return result;
    }
};
} // namespace certpolysolve::detail
