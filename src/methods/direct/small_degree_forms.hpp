#pragma once
// The closed forms of degree 3 for Auto's direct route (small_degree.cpp; the design and its measurements:
// dev/docs/small-degree-closed-forms-20261026/). The exact part is division-free: the cubic's coefficients as Gaussian
// integers, and its homogeneous invariants D0 = b^2 - 3ac, D1 = 2b^3 - 9abc + 27a^2 d and E = D1^2 - 4 D0^3, which is
// -27 a^2 times the discriminant. E != 0 proves the cubic square-free; E = 0 gives its multiple roots exactly (a
// triple root -b/(3a) where D0 = 0, else the double root (9ad - bc)/(2 D0) and the simple root (4abc - 9a^2 d - b^3) /
// (a D0)). The roots of a square-free cubic are boxes in ball arithmetic at one precision:
// - real coefficients, E < 0 (three real roots): Viete's form, (2 sqrt(D0) cos((theta + 2 pi k)/3) - b)/(3a) with
//   theta = atan2(sqrt(-E), -D1), in real balls, so each box is exactly real;
// - real coefficients, E > 0 (one real root and a conjugate pair): the real Cardano form, C = ((|D1| + sqrt E)/2)^(1/3),
//   K = C + D0/C and L = C - D0/C each in the form without cancellation for the exact signs of D0, the real root
//   (-b - sign(D1) K)/(3a) and the pair (-b + sign(D1) K/2)/(3a) +- i sqrt(3) L/(6a);
// - complex coefficients: Cardano's form, C a cube root of the larger of (D1 +- sqrt E)/2, the square and cube roots
//   rotated away from the cut, the roots -(b + xi^k C + D0/(xi^k C))/(3a), xi = exp(2 pi i/3).
// The quartic a x^4 + b x^3 + c x^2 + d x + e (ExactQuartic) has D0 = c^2 - 3bd + 12ae,
// D1 = 2c^3 - 9bcd + 27b^2 e + 27ad^2 - 72ace and E = D1^2 - 4 D0^3 = -27 disc, and in y = x + b/(4a) the depressed
// y^4 + p y^2 + q y + r with p = P8/(8a^2), q = Q8/(8a^3), r = R256/(256a^4) (P8 = 8ac - 3b^2, Q8 = b^3 - 4abc + 8a^2 d,
// R256 = -3b^4 + 16ab^2 c - 64a^2 bd + 256a^3 e). E = 0: for real coefficients the exact square-free factorization
// (factors of degree at most 2, each with its multiplicity); complex coefficients with a multiple root are declined.
// Q8 = 0: the biquadratic, w = y^2 a root of w^2 + p w + r and y = +- sqrt(w). Otherwise Descartes' factors
// (y^2 + s y + t)(y^2 - s y + v) with z = s^2 a root of the resolvent, z = (-2p + (Q + D0/Q)/a)/3 for a cube root Q of
// (D1 + sqrt E)/2, t, v = (p + z -+ q/s)/2: for real coefficients the largest real root (positive since q != 0),
// 2 sqrt(D0) cos(theta/3) with theta = atan2(sqrt(-E), D1) for Q + D0/Q where E < 0 and the real Cardano form's
// sign(D1) K where E > 0; for complex coefficients the root of largest modulus. Each factor's discriminant, decided in
// balls (the precision raised where its ball holds 0), gives two exactly real roots or a conjugate pair proved not real.
// Real coefficients are normalized to a positive leading coefficient first (the real forms assume it). Every value is
// a ball, so each box holds its exact root; the caller publishes the boxes and tests the goal and disjointness.
//
// State: the exact coefficients and invariants (ExactCubic, ExactQuartic). Entry point: direct_small_degree
// (small_degree.cpp); tests/test_small_degree_cubic.cpp, tests/test_small_degree_quartic.cpp.
#include "core/cover.hpp"
#include "core/raii.hpp"
#include <optional>
#include <span>
#include <vector>

namespace certpolysolve::detail {
struct GaussianInteger {
    Integer re, im;
    bool zero() const { return fmpz_is_zero(re.value()) && fmpz_is_zero(im.value()); }
    bool real() const { return fmpz_is_zero(im.value()); }
};

class ExactCubic {
  public:
    // The coefficients lowest first (d, c, b, a), a nonzero.
    explicit ExactCubic(std::span<const GaussianInteger, 4> coefficients);
    bool square_free() const { return !E_.zero(); }
    bool real() const { return real_; }
    // The roots' boxes at `bits`: three simple roots, or the exact multiple roots as one box each with its count and
    // multiplicity; nullopt where a choice of the forms is not decided at this precision (a cube root's argument
    // whose ball holds 0).
    std::optional<std::vector<Root>> roots(slong bits) const;

  private:
    GaussianInteger a_, b_, c_, d_, D0_, D1_, E_;
    bool real_ = false;
};

class ExactQuartic {
  public:
    // The coefficients lowest first (e, d, c, b, a), a nonzero.
    explicit ExactQuartic(std::span<const GaussianInteger, 5> coefficients);
    bool square_free() const { return !E_.zero(); }
    // A multiple root of a quartic with complex coefficients (its factorization would need Gaussian gcds): the
    // general route answers it.
    bool declined() const { return declined_; }
    // As ExactCubic::roots; nullopt also where a sign the forms need (the resolvent root's, a factor's discriminant's)
    // is not decided at this precision.
    std::optional<std::vector<Root>> roots(slong bits) const;

  private:
    struct Factor {
        std::vector<GaussianInteger> coefficients;  // lowest first, degree 1 or 2, real
        uint64_t multiplicity;
    };
    GaussianInteger a_, b_, c_, d_, e_, D0_, D1_, E_, P8_, Q8_, R256_;
    std::vector<Factor> factors_;
    bool real_ = false, declined_ = false;
};
} // namespace certpolysolve::detail
