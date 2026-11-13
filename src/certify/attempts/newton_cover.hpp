#pragma once
// The certificate attempt (machine_newton_cover): a binary64 packet on the whole input, then a Newton certificate at
// every point (certify_newton_point) against the original source, repair rounds on the refused points, and what a
// refused attempt keeps (PartialNewtonCover). Used by
// the classical and secular routes, the selector's attempt before piecewise, and piecewise. Defined in
// attempts/newton_cover.cpp. Proof predicates live in the parent certify module.
#include <certpolysolve/certpolysolve.h>
#include "core/cover.hpp"
#include "core/number.hpp"
#include <algorithm>
#include <cstdint>
#include <memory>
#include <span>
#include <vector>

namespace certpolysolve::detail {
class PolynomialSource;
struct CandidateStatistics;
struct EnginePolicy;
// What a refused certificate attempt had certified when it ended
// (machine_newton_cover's `partial`): every point, at its certified centre or
// where the packet and the rounds left it, and for each certified point its
// region, which holds at least one root of the source. The certified regions
// are pairwise disjoint. A refusal throws none of them away: a cluster below
// the certificate's resolution (two roots 2^-384 apart at a 333-bit goal)
// refuses the attempt while every other point holds its certificate, and a
// caller may keep those and iterate only the refused points.
struct PartialNewtonCover {
    std::unique_ptr<RootVector> points;
    std::vector<Root> roots; // indexed like the points; meaningful where certified
    std::vector<uint8_t> certified;
    slong certificate_bits = 0;
    // Every point was tried (the pass ended within the refusals' allowance)
    // and the rounds iterated the refused ones: the attempt ended on points
    // its own multiprecision iteration did not certify.
    bool pass_complete = false;
    size_t held() const { return static_cast<size_t>(std::count(certified.begin(), certified.end(), uint8_t(1))); }
};
// Optional owned continuation hints have no proof authority. A failed shortcut
// leaves the supplied cover unchanged; finite hints may still be returned.
// Without scaled_packet, a polynomial outside binary64's range is refused at
// once (no packet in ScaledBinary64): for a caller whose alternative is
// cheaper than that packet (the attempt before piecewise). With partial, a
// refusal after the certificates began returns what they hold (above).
//
// settled (indexed like the points, with partial and the coefficient bound): a
// retry from an iteration's points. No binary64 packet runs (it can only round
// points an iteration refined), the settled points are neither moved nor
// certified, and what the others certify is returned through partial; the
// attempt then never returns a cover itself.
//
// fallback_operations: what the caller's route costs if the attempt fails, in the selector's complex coefficient
// operations (auto_selection.cpp): the attempt as a whole (its binary64 packet's bound, its certificates and
// refusals, its repair rounds) spends no more, so a failed attempt at most doubles that route's cost. Zero means
// the global recipes' own minimum budget, initial_sweep_allowance sweeps of degree x (degree + support).
bool machine_newton_cover(std::vector<Root> &, PolynomialSource &, acb_srcptr initial,
                          const cps_options &, slong precision, bool allow_underflow,
                          certpolysolve::detail::CandidateStatistics &, bool coefficient_bound = false, const EnginePolicy *diagnostic = nullptr,
                          std::unique_ptr<RootVector> *retained_proposals = nullptr,
                          bool *machine_converged = nullptr, slong *certificate_bits = nullptr,
                          bool scaled_packet = true, PartialNewtonCover *partial = nullptr,
                          std::span<const uint8_t> settled = {}, uint64_t fallback_operations = 0);
} // namespace certpolysolve::detail
