#pragma once
#include <chrono>
#include <type_traits>

namespace certpolysolve::detail {
// A nonnegative limit in seconds, on the solve's monotonic clock. Clamp before
// converting to integer ticks and before adding to the clock: either operation
// can overflow even though the public limit is a valid finite double.
inline std::chrono::steady_clock::time_point deadline_after(
    double seconds, std::chrono::steady_clock::time_point now) {
    using Clock = std::chrono::steady_clock;
    using Duration = Clock::duration;
    static_assert(std::is_integral_v<Duration::rep>);
    const double ticks = std::chrono::duration<double, Duration::period>(
        std::chrono::duration<double>(seconds)).count();
    // The integer maximum may round up when converted to double; strict <
    // still admits only representable integers. Infinity also takes this arm.
    if (!(ticks < static_cast<double>(Duration::max().count())))
        return Clock::time_point::max();
    const Duration delay{static_cast<Duration::rep>(ticks)};
    const auto latest = Clock::time_point::max();
    return now > latest - delay ? latest : now + delay;
}
} // namespace certpolysolve::detail
