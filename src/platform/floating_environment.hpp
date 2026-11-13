#pragma once
#include <cfenv>

namespace certpolysolve::detail {
// The floating-point environment of a solve call (docs/api-spec-0.5.md section 8.2, change C28): the default one,
// whatever the calling thread has set (round-to-nearest, every exception masked and its flag clear, flush-to-zero
// and denormals-are-zero off), so the solve's path and answer do not depend on the caller's rounding mode or MXCSR.
// The caller's environment, its raised flags included, is saved at entry and given back when the call returns (the
// thread-safety audit's C3 and E4). A cancellation callback runs in the default environment too. If the
// environment cannot be saved, it is left as it is.
class DefaultFloatingPoint {
    std::fenv_t caller_{};
    bool saved_;

  public:
    DefaultFloatingPoint() noexcept : saved_(std::fegetenv(&caller_) == 0) {
        if (saved_) std::fesetenv(FE_DFL_ENV);
    }
    ~DefaultFloatingPoint() {
        if (saved_) std::fesetenv(&caller_);
    }
    DefaultFloatingPoint(const DefaultFloatingPoint &) = delete;
    DefaultFloatingPoint &operator=(const DefaultFloatingPoint &) = delete;
};
} // namespace certpolysolve::detail
