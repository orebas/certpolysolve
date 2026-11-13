#pragma once
// Caps the address space for the enclosing scope, so that a regression in a
// memory bound fails with an allocation error instead of exhausting the
// machine. (Not under a sanitizer, whose shadow memory needs the space.)
#include <algorithm>
#if defined(__linux__)
#include <sys/resource.h>
#endif

// AddressSanitizer: GCC defines __SANITIZE_ADDRESS__, Clang answers __has_feature(address_sanitizer) only (clang 21
// does not define the macro; checked 20 October 2026, when piecewise_geometry failed under Clang's ASan with the
// runtime's own mmap refused under the cap).
#if defined(__SANITIZE_ADDRESS__)
#define CERTPOLYSOLVE_TEST_ADDRESS_SANITIZER 1
#elif defined(__has_feature)
#if __has_feature(address_sanitizer)
#define CERTPOLYSOLVE_TEST_ADDRESS_SANITIZER 1
#endif
#endif

class AddressSpaceCap {
#if defined(__linux__) && !defined(CERTPOLYSOLVE_TEST_ADDRESS_SANITIZER)
    rlimit saved_{};
    bool active_ = false;

  public:
    explicit AddressSpaceCap(rlim_t bytes) {
        if (getrlimit(RLIMIT_AS, &saved_) != 0) return;
        rlimit capped = saved_;
        capped.rlim_cur = saved_.rlim_max == RLIM_INFINITY ? bytes : std::min(bytes, saved_.rlim_max);
        active_ = setrlimit(RLIMIT_AS, &capped) == 0;
    }
    ~AddressSpaceCap() {
        if (active_) setrlimit(RLIMIT_AS, &saved_);
    }
#else
  public:
    explicit AddressSpaceCap(unsigned long long) {}
#endif
    AddressSpaceCap(const AddressSpaceCap &) = delete;
    AddressSpaceCap &operator=(const AddressSpaceCap &) = delete;
};
