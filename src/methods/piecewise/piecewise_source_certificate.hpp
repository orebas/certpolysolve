#pragma once
#include "core/cover.hpp"
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace certpolysolve::detail {
class PiecewiseFinisherSource;

// Success of certification against one source's exact coefficient vector.
// Only that source can create it. The controller still proves owner identity,
// disjointness and available degree before admission; this type alone is not
// a complete cover. Call-scoped: the originating source must remain alive.
class PiecewiseSourceCertificate {
    friend class PiecewiseFinisherSource;
    const PiecewiseFinisherSource *source_;
    Root root_;

    PiecewiseSourceCertificate(const PiecewiseFinisherSource &source, Root &&root)
        : source_(&source), root_(std::move(root)) {}

  public:
    PiecewiseSourceCertificate(const PiecewiseSourceCertificate &) = delete;
    PiecewiseSourceCertificate &operator=(const PiecewiseSourceCertificate &) = delete;
    PiecewiseSourceCertificate(PiecewiseSourceCertificate &&other) noexcept(std::is_nothrow_move_constructible_v<Root>)
        : source_(std::exchange(other.source_, nullptr)), root_(std::move(other.root_)) {}
    PiecewiseSourceCertificate &operator=(PiecewiseSourceCertificate &&) = delete;

    bool belongs_to(const PiecewiseFinisherSource &source) const { return source_ == &source; }
    const Root &proof() const & {
        if (!source_) throw std::logic_error("consumed piecewise source certificate");
        return root_;
    }
    const Root &proof() const && = delete;
    Root release() && {
        if (!source_) throw std::logic_error("consumed piecewise source certificate");
        source_ = nullptr;
        return std::move(root_);
    }
};
} // namespace certpolysolve::detail
