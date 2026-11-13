#pragma once
#include <optional>
#include <vector>
#include "numeric/exact/exact_complex.hpp"
#include "core/number.hpp"
#include "numeric/exact/preparation_diagnostics.hpp"
namespace certpolysolve::detail {
bool center_materialization_fits(const acb_t, uint64_t bits);
uint64_t center_owned_bits(const fmpz_poly_t);
uint64_t center_owned_bits(const GaussianPolynomial &);
// Exact squared coefficient norms compare monic heights without a floating
// logarithm or a dependence on the polynomial's nonzero scalar factor.
bool gaussian_center_improves_height(const GaussianPolynomial &, const GaussianPolynomial &);
// Necessary condition for the existing height-improvement admission: a monic
// coefficient height cannot fall below zero (the leading coefficient is one).
bool center_height_allows_gain(const fmpz_poly_t);
// A lower bound on the monic coefficient height (largest coefficient's bit length minus the leading
// one's) of the exact translate of p to `center`, as center_integer builds it, from the constant term
// alone: that term is p(center) up to the common factor of the translate, and one ball evaluation
// gives its binary exponent. No value when the evaluation cannot exclude p(center) = 0. A cost
// predictor for the height rules that accept or refuse a centering; never a certificate.
std::optional<slong> translated_height_lower_bound(const fmpz_poly_t p, const fmpq_t center);
bool center_height_allows_gain(const GaussianPolynomial &);
// Scale-invariant centering heuristics (EnginePolicy::scale_invariant_centering).
// Plain monic height is not scale-invariant: when every root lies well inside
// the unit disk, every coefficient is below the leading one, so a cloud offset
// by about its own radius was never centered. These compare heights at the
// polynomial's own root scale instead. Cost heuristics, never proofs.
//
// A profile estimates log2|c_k| for each nonzero coefficient from bit lengths
// (within 1.5 bits; no exact products, which cost more than the shift itself
// on large rational coefficients).
struct HeightProfile {
    std::vector<slong> bits;
    std::vector<bool> present;
};
HeightProfile height_profile(const fmpz_poly_t);
HeightProfile height_profile(const GaussianPolynomial &);
// Monic height at scale t (x -> 2^t x): max_k (bits_k + t k) - (bits_n + t n).
// It never increases with t.
double height_at(const HeightProfile &, double t);
// log2 of the geometric-mean modulus of the nonzero roots, the scale that
// balances the lowest-order and leading terms. It is not rounded to an
// integer: rounding can put the scale up to half a bit inside a circle of
// roots, where the height reads up to n/2 bits that no translation removes.
double root_scale(const HeightProfile &);
// The largest excess of `after`'s height over `before`'s at any scale. Each
// height is a maximum of lines in t minus a line, so the excess is largest at
// a breakpoint of either upper hull of (k, bits_k), or in the limits: it
// tends to 0 as t -> +inf, and as t -> -inf to a constant, or to +inf when
// `after` vanishes at zero to a lower order.
double worst_height_excess(const HeightProfile &after, const HeightProfile &before);
// The roots a translation by c carries out of binary64's reach with nothing for a cluster restart to hold.
// A root x with |x| < |c| 2^-reach needs `reach_bits` more bits than its own in the translated coordinates:
// the secular iteration's binary64 stage cannot tell it from its neighbours there, where at the origin it
// could. p's Newton polygon (the profile's upper hull) says how many roots lie below that scale. One certified
// group that holds them all (Pellet's test at the origin, certified_group_scale, at a radius below |c|/2) is a cluster
// about the image of the origin, which a restart places whole (fresh clustered-719's 180 roots within 2^-85).
// Without such a gap they shade into the other roots scale by scale, the translated iteration peels them one
// scale at a time in multiprecision, and the untranslated one finds them in binary64 (roots q^k beside a
// factor away from the origin: the spiral beside a ring, little q-Legendre times a shifted Chebyshev).
// Returns their number when there are two or more and no such group, zero otherwise. A cost estimate for the
// admission rule, never a proof. History: docs/code-history-20261021.md, centering.hpp 2.
slong unseparated_origin_roots(const HeightProfile &, double center_log2, double reach_bits);
// A proposed centre computed in ball arithmetic: each part is rounded to `bits`
// significant bits (a low-denominator dyadic keeps the exact translation
// cheap), and a part whose enclosure contains zero is set to zero. Such a part
// is not known to differ from zero, and its midpoint is rounding noise: the
// centroid -c_{n-1}/(n c_n) of a complex multiple of a real polynomial is
// real, and as a 64-bit ball quotient its imaginary part can be about 2^-60
// of the real one. Rounded to 8 bits of its own it would give the centre an
// 80-bit denominator, and the exact shift would grow every coefficient by
// degree x 80 bits. A part that is known to be nonzero is rounded as before,
// however small it is beside the other: a real polynomial moved by 3i/8 has
// the centroid's imaginary part -3/8 exactly, and the translation by it is
// what makes the centred polynomial real again. Only midpoints are rounded; a
// non-finite centre is left as it is.
// History: docs/code-history-20261021.md, centering.hpp 3.
void round_center(acb_t center, slong bits);
// Exact change of variable: output is a nonzero constant multiple of p(x+c).
// Preflight limits bound coefficient growth before materializing denominator powers.
bool center_integer(fmpz_poly_t out, const fmpz_poly_t p, const fmpq_t center, uint64_t budget,
                    PreparationDiagnostic * = nullptr
                    , double *research_work = nullptr
                    );
bool center_gaussian(GaussianPolynomial &out, const GaussianPolynomial &p, const GaussianRational &center,
                     uint64_t budget, PreparationDiagnostic * = nullptr
                     , double *research_work = nullptr
                     );

// Owned, single-source exact translation cache. A failed preparation is cached
// as unavailable and can never expose the polynomial from an earlier centre.
class ExactTaylorCache {
    IntegerPolynomial integer_, shifted_integer_;
    GaussianPolynomial gaussian_, shifted_gaussian_;
    ComplexBall center_;
    uint64_t budget_, operations_;
    bool allow_expansion_;
    bool attempted_ = false, available_ = false, use_integer_ = false;
    bool research_service_meter_ = false;
    double research_work_ = 0;

  public:
    ExactTaylorCache(const IntegerPolynomial &polynomial, uint64_t budget, bool allow_expansion = false, uint64_t operations = UINT64_MAX)
        : integer_(polynomial), budget_(budget), operations_(operations), allow_expansion_(allow_expansion) {}
    ExactTaylorCache(const GaussianPolynomial &polynomial, uint64_t budget, bool allow_expansion = false, uint64_t operations = UINT64_MAX)
        : gaussian_(polynomial), budget_(budget), operations_(operations), allow_expansion_(allow_expansion) {}
    // A nonzero scalar multiple of p(center+x); the source normalizes to monic.
    bool enclose(acb_poly_t out, const acb_t center, slong precision);
    // Cumulative raw word-volume estimate. A reader shares this cache's
    // ownership; it neither prices the work nor changes exact admission.
    void research_service_meter(bool enabled) { research_service_meter_ = enabled; }
    double research_work() const { return research_work_; }
};
} // namespace certpolysolve::detail
