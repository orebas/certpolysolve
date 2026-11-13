#pragma once
// The separation of factor covers that meet (W8's S8, rule 2), one rule for the global driver (SolveState) and
// piecewise (PreparedPiecewiseState): a factor whose cover meets another factor's is solved to a tighter goal, the
// caller's goal with separation bits added, and only the factor whose region is the wider of the first meeting pair
// is tightened, its goal's accuracy doubled by at least 16 bits.
//
// State: none. Entry points: separated_options (SolveState::Impl::factor_options, piecewise's factor loop),
// next_separation and factor_to_separate (SolveState::Impl::advance's overlap branch, piecewise's factor loop).
// Not here: the covers' merging (merge_overlapping, certification.cpp) and the precision a lifted cover needs.
#include "core/query.hpp"
#include <cstddef>
#include <optional>
#include <span>
#include <vector>

namespace certpolysolve::detail {
// The caller's goal with `separation` more bits: the absolute tolerance dropped, and a cluster radius scaled by
// 2^-separation with it. The goal itself when `separation` is 0.
cps_options separated_options(const cps_options &base, uint64_t separation);
// A factor's next separation bits when its cover still meets another's: the goal's accuracy (`bits` plus the
// separation) doubled, by at least 16 bits, at most precision_max - bits.
uint64_t next_separation(uint64_t bits, uint64_t separation);
// The factor to tighten when two factors' covers meet: of the first meeting pair in the order (i, j) with j < i (and
// each cover's regions in order), the factor whose region is the wider (the sum of its box's real and imaginary
// radii; i on a tie). None when no two covers meet.
std::optional<size_t> factor_to_separate(std::span<const std::vector<Root> *const> covers);
} // namespace certpolysolve::detail
