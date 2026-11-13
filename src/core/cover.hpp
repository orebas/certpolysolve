#pragma once
// A certified cover and its regions: Endpoint (an exported dyadic end), Root (one region: its ball, its four ends,
// count, multiplicity, realness and certificate), and the predicates every route applies to a cover (disjointness,
// merging overlapping regions, realness by reflection, the completion of a partial cover at a stop, the goal tests,
// the cover order dovetailing records, the Cauchy disk). Defined in certification.cpp. Part of internal.hpp.
#include <certpolysolve/certpolysolve.h>
#include "core/number.hpp"
#include <array>
#include <cstdint>
#include <string>
#include <tuple>
#include <vector>

namespace certpolysolve::detail {
// An end's text as the API returns it: the mantissa and the exponent of mantissa@exponent, in decimal.
struct EndpointText {
    std::string mantissa = "0", exponent = "0";
    cps_dyadic_view view() const { return {mantissa.c_str(), exponent.c_str()}; }
};
// A published region's end: an exact dyadic number, held as an arf. It stays binary for the whole solve. The API's
// text is made once, when a result is published (publish_text): text exists only at input and output (since 26
// October 2026; until then every publish wrote the ends as decimal strings and the certificates, the midpoint order
// and the real regions parsed them back). Its text is the mantissa and exponent arf_get_fmpz_2exp gives (an odd
// mantissa, or 0 with exponent 0). Since 27 October 2026 the end is the arf itself, not that mantissa and exponent as
// two fmpz: FLINT keeps an arf's mantissa in place up to ARF_NOPTR_LIMBS limbs (128 bits), so making or copying a
// 53- to 101-bit end allocates nothing, where an fmpz above 62 bits takes an mpz (and text() makes the API's text
// without one).
class Endpoint {
    arf_t value_;

  public:
    Endpoint() { arf_init(value_); }
    explicit Endpoint(const arf_t x) : Endpoint() { arf_set(value_, x); }
    ~Endpoint() { arf_clear(value_); }
    Endpoint(const Endpoint &x) : Endpoint() { arf_set(value_, x.value_); }
    Endpoint(Endpoint &&x) noexcept : Endpoint() { arf_swap(value_, x.value_); }
    Endpoint &operator=(Endpoint x) noexcept {
        arf_swap(value_, x.value_);
        return *this;
    }
    arf_srcptr value() const & noexcept { return value_; }
    arf_srcptr value() const && = delete;
    void negate() noexcept { arf_neg(value_, value_); }
    // The end as mantissa * 2^exponent (arf_get_fmpz_2exp), and the API's two decimal texts of them (text(), equal
    // to fmpz_get_str of that mantissa and exponent).
    void get(fmpz_t mantissa, fmpz_t exponent) const { arf_get_fmpz_2exp(mantissa, exponent, value_); }
    EndpointText text() const;
    std::string mantissa_text() const;
    std::string exponent_text() const;
};
struct Root {
    ComplexBall box;
    std::array<Endpoint, 4> endpoints;
    uint64_t count = 1, multiplicity = 0;
    cps_realness realness = CPS_REAL_UNKNOWN;
    std::string certificate;
    // Outward-rounded exact dyadic endpoints. box is replaced with a containing
    // ball so later disjointness and accuracy checks also cover the exported box.
    void publish(slong precision);
};
// Exact rectangle geometry shared by ownership, guards and cover checks.
inline std::array<RealBall, 4> exact_endpoints(const Root &root) {
    std::array<RealBall, 4> bounds;
    for (size_t k = 0; k < bounds.size(); ++k)
        arb_set_arf(bounds[k].value(), root.endpoints[k].value());
    return bounds;
}
inline bool contains_rectangle(const std::array<RealBall, 4> &a, const std::array<RealBall, 4> &b) {
    return arb_le(a[0].value(), b[0].value()) && arb_le(b[1].value(), a[1].value()) &&
           arb_le(a[2].value(), b[2].value()) && arb_le(b[3].value(), a[3].value());
}
inline bool contains_rectangle(const Root &outer, const Root &inner) {
    return contains_rectangle(exact_endpoints(outer), exact_endpoints(inner));
}

// A generating precision enables the bounded exact-endpoint sweep. Zero keeps
// the compact-exponent pair path; neither choice changes the exact predicate.
bool disjoint(const std::vector<Root> &roots, slong generating_precision = 0);
// By real part, then imaginary part, of the midpoints.
void sort_roots(std::vector<Root> &roots);
// Merge overlapping regions of covers of coprime factors, or of a complete
// algebraic root enumeration, until disjoint. For the latter, each root copy
// must be assigned exactly once to a containing branch box, with multiplicity
// proved separately; the whole root multiset must be present. A count ledger
// on arbitrary numerical candidates is not enough. Each resulting disjoint
// union holds exactly the sum of its assigned counts. Merging clears known
// multiplicity and retains realness only when both regions were proved real.
// Returns true when the cover was disjoint already (nothing merged): the caller
// need not test it again.
bool merge_overlapping(std::vector<Root> &cover, slong precision);
// For finite regions: whether another region meets the reflection of region i
// in the real axis (acb_overlaps of the conjugate box and the other box).
std::vector<bool> reflection_met(const std::vector<Root> &roots, slong generating_precision = 0);
// For a complete disjoint cover of the roots of a polynomial with real
// coefficients: every region of one distinct root whose reflection in the
// real axis meets no other region holds a real root; it is replaced by its
// real segment, marked real and published at `precision`. With
// `keep_real_ends`, for a cover already published and judged, only the
// imaginary ends become 0 and the real ends stay as they were (publishing
// again would round them outward once more).
void prove_real_by_reflection(std::vector<Root> &cover, slong precision, bool keep_real_ends = false);
// Complete a partial cover at a stop short of the goal. The certified regions
// stay as they are; the roots not yet certified are enclosed by one square
// about the origin whose count is proved on the polynomial
// (certify_cluster_box: Pellet's test on the disks sandwiching it). The square
// is the smallest of a geometric ladder of half-sides, from Cauchy's bound
// downward, that meets every certified region either wholly inside or wholly
// outside it and whose count, with the regions outside, makes the degree: the
// regions inside are then among its roots and are dropped, so the cover is
// disjoint and its counts sum to the degree. When no half-side of the ladder
// qualifies, the cover is Cauchy's global square alone. Honest, not a rescue:
// the coarse region meets no goal. Returns false, leaving the cover as it was,
// only when the leading coefficient contains zero (no root bound exists). With
// the caller's options, its cancellation is polled before every half-side and,
// inside the disk tests, before every Taylor shift and Graeffe step; a
// cancelled ladder keeps the best cover it reached. Null never stops it.
// History: docs/code-history-20261021.md, internal.hpp 2.
bool complete_cover_at_stop(std::vector<Root> &cover, const acb_poly_t polynomial, slong precision,
                            const cps_options *cancellation = nullptr);
// Mark every rectangle intersecting another; outward sweep bounds only prune
// impossible overlaps. Used when a disjoint partial cover must be preserved.
std::vector<bool> overlapping_roots(const std::vector<Root> &roots);
// The relative test of accuracy_met and approximation_met (coordinate_accuracy_met, certification.cpp) in its two
// forms, exported for their test: with E the box's radii summed and L the larger absolute lower bound of its parts,
// E <= 2^-bits L in ball arithmetic at `precision`, and a binary64 form, tried first, that holds only where the ball
// form does.
bool relative_accuracy_ball(const acb_t box, uint64_t bits, slong precision);
bool relative_accuracy_binary64(const acb_t box, uint64_t bits, slong precision);
// acb_contains(outer, inner), decided first by two rounded sums a part: inner lies in outer where, in each part,
// |m_outer - m_inner| + r_inner <= r_outer, and that holds where the left side rounded up at 64 bits is at most
// r_outer; acb_contains decides the rest, so the answer is always acb_contains's. FLINT's own first test rounds the
// ends to 30 bits, too coarse for a box 2^-60 of its centre wide: the compensated certificate's check reached FLINT's
// exact four-term sum at every root (dev/docs/small-degree-design-pass-20261026.md).
bool box_contains(const acb_t outer, const acb_t inner);
bool approximation_met(const Root &root, uint64_t bits, const Number &absolute_tolerance, slong precision);
bool accuracy_met(const Root &root, uint64_t bits, const Number &absolute_tolerance, slong precision);
bool cluster_radius_met(const Root &root, const Number &radius, slong precision);
bool meets_goal(const Root &root, const cps_options &options, slong precision);
// A certified cover's structure: (roots in unresolved clusters, the largest
// such cluster, minus the number of regions). Smaller is better.
std::tuple<uint64_t, uint64_t, int64_t> cover_structure(const std::vector<Root> &roots);
// log2 of the largest enclosure radius among unresolved clusters (-inf if none).
double unresolved_radius_log2(const std::vector<Root> &roots);
// The order dovetailing records progress in: a better cover structure, then
// a smaller largest unresolved radius. True when cover a is ahead of cover b.
bool cover_ahead(const std::vector<Root> &a, const std::vector<Root> &b);
bool prefer_cover(const std::vector<Root> &candidate, const std::vector<Root> &previous, slong precision);
bool cauchy_cluster(Root &root, const acb_poly_t polynomial, slong precision);
// The ball that Number(mantissa@exponent).enclose(out, bits) gives for an end: the mantissa over 1 rounded at bits
// (arb_set_fmpq, as Number::enclose rounds it) times 2^exponent. History: until 26 October 2026 it built and parsed
// the literal, then (A6) parsed the ends' decimal text; the ends are binary since.
// The end rounded to `bits` (ARB_RND, an ulp of radius where inexact). Until 27 October 2026 the mantissa was rounded
// by arb_set_fmpq and then scaled by 2^exponent; arb_set_round of the exact end gives the same ball: the same rounded
// mantissa, and the radius arf_mag_set_ulp gave, 2^(exponent of the rounded end - bits), as arf_mag_add_ulp adds it
// to a zero radius (FLINT 3.6, arb/div.c and arf.h).
inline void endpoint_ball(arb_t out, const Endpoint &e, slong bits) {
    arb_set_arf(out, e.value());
    arb_set_round(out, out, bits);
}

} // namespace certpolysolve::detail
