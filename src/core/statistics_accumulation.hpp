#pragma once
// A caller chooses both a field list and its arithmetic. Ordinary unsigned
// accumulation retains wraparound; handoff/publication accounting saturates.
// Neither operation applies to maxima, current gauges or local diagnostics.
#include "core/resource_limits.hpp"
#include <cstddef>
#include <cstdint>

namespace certpolysolve::detail::statistics {
template <class Statistics, size_t N>
constexpr void add_ordinary(Statistics &total, const Statistics &other,
                            uint64_t Statistics::*const (&fields)[N]) {
    for (const auto field : fields) total.*field += other.*field;
}
template <class Statistics, size_t N>
constexpr void add_saturating(Statistics &total, const Statistics &other,
                              uint64_t Statistics::*const (&fields)[N]) {
    for (const auto field : fields) total.*field = resources::add(total.*field, other.*field);
}
} // namespace certpolysolve::detail::statistics
