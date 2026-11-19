// The choice of a kernel's variant (kernel_variant.hpp): the CPU is asked once, at the first call.
#include "platform/kernel_variant.hpp"
#include <atomic>

namespace certpolysolve::detail {
namespace {
constexpr int unset = -1;
// The variant in force, or unset before the first question. Relaxed order: every value a thread can read is a
// variant this CPU runs, and all of them return the same bits.
std::atomic<int> chosen{unset};
KernelVariant best() noexcept {
    return kernel_variant_available(KernelVariant::Avx2Fma) ? KernelVariant::Avx2Fma : KernelVariant::Reference;
}
} // namespace

bool kernel_variant_available(KernelVariant variant) noexcept {
    if (variant == KernelVariant::Reference)
        return true;
#ifdef CERTPOLYSOLVE_HAVE_AVX2_FMA_KERNELS
    // The compiler's runtime asks CPUID, and the operating system whether it saves the wide registers.
    return __builtin_cpu_supports("avx2") && __builtin_cpu_supports("fma");
#else
    return false;
#endif
}

KernelVariant kernel_variant() noexcept {
    int now = chosen.load(std::memory_order_relaxed);
    if (now == unset) {
        now = static_cast<int>(best());
        chosen.store(now, std::memory_order_relaxed);
    }
    return static_cast<KernelVariant>(now);
}

KernelVariant set_kernel_variant(KernelVariant variant) noexcept {
    const KernelVariant before = kernel_variant();
    if (kernel_variant_available(variant))
        chosen.store(static_cast<int>(variant), std::memory_order_relaxed);
    return before;
}

namespace {
// Whether the selector's prices are pinned to the reference's: a hook of the tests (set_reference_prices). The
// library reads no environment, so nothing outside the process pins them.
std::atomic<bool> pinned{false};
} // namespace

KernelVariant price_variant() noexcept {
    return pinned.load(std::memory_order_relaxed) ? KernelVariant::Reference : kernel_variant();
}

bool set_reference_prices(bool pin) noexcept { return pinned.exchange(pin, std::memory_order_relaxed); }
} // namespace certpolysolve::detail
