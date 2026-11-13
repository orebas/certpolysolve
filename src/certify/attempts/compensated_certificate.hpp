#pragma once
// The certificate attempt's first tier (proposal B of dev/docs/easy-inputs-profile-20261025.md; the design and the
// proof are in dev/docs/small-input-fast-path-20261026.md): one root in a goal-sized rectangle, proved in binary64
// arithmetic with error-free transformations, where the ball certificate (certificate_attempt.cpp) would spend about a
// hundred 85-bit ball operations a root.
//
// At the binary64 point x nearest a candidate, p(x) is evaluated by the compensated Horner scheme (its error-free
// transformations as in Cameron and Graillat, "Accurate Horner methods in real and complex floating-point arithmetic",
// BIT 64:17, 2024, Theorems 2.4 and 2.5): the rounded value s and a correction c with p(x) = s + E(x) exactly, E the
// polynomial of the rounding errors, and c approximating E(x) to within a computed bound. p'(x) is the derivative
// recurrence on the rounded values, with its own bound, and |p''| is bounded on a disk about x by the coefficient
// majorant. One Newton step gives the centre z1 = x - c_hat as an exact two-word number, and Rouché's theorem on the
// circles of radius r and 2r about z1 (r the goal's half-width) proves exactly one root in the rectangle z1 +- r(1 +
// i): the disk of radius r inside it holds one, the disk of radius 2r around it no other. Where the majorant's term
// fails those circles (at a coarse goal), the same test runs on circles of a smaller radius from the same evaluation
// (heuristics::compensated_circle_curvature_share), whose rectangle meets the goal as well. Every bound is rigorous for
// every polynomial in the coefficient balls: rounding is bounded by the standard model (each operation within a
// relative 2^-53, products possibly off by an absolute 2^-1075 below the normal range, sums exact there), so the proof
// needs neither interval arithmetic nor the floating-point environment's flags.
//
// State: the coefficients, prepared once per attempt (CompensatedCertificate). Entry point: machine_newton_cover's
// certificate (certificate_attempt.cpp), before the ball certificate; tests/test_compensated_certificate.cpp.
// Not here: the ball certificate and the attempt's rounds (certificate_attempt.cpp), the goal tests (certification.cpp).
#include "core/cover.hpp"
#include "core/query.hpp"
#include <flint/acb_poly.h>
#include <cstdint>
#include <vector>

namespace certpolysolve::detail {
// What the compensated evaluation at a binary64 point gives (CompensatedCertificate::evaluate): for every polynomial g
// in the coefficient balls (scaled by 2^-scale), |g(x) - (value + correction)| <= value_error and
// |g'(x) - slope| <= slope_error.
// Where prove's test fails at a binary64 point and Newton's iteration from it is in its quadratic range, the binary64
// point nearest the Newton step (CompensatedCertificate::prove; certify evaluates there again).
struct CompensatedRetry {
    double re = 0, im = 0;
    bool available = false;
};
struct CompensatedEvaluation {
    double x_re = 0, x_im = 0;
    double value_re = 0, value_im = 0, correction_re = 0, correction_im = 0;
    double slope_re = 0, slope_im = 0;
    double value_error = 0, slope_error = 0;
    // An upper bound on |x|, and one on max(1, |x|)^degree: the factor the underflow allowances use.
    double modulus = 0;
    double amplification = 1;
    bool ok = false;
};

class CompensatedCertificate {
  public:
    // From the coefficient enclosure of the attempt's source (the polynomial itself; its roots do not depend on the
    // scale). Unavailable (available() false) where the degree is below 1, a coefficient is not finite, or the
    // leading coefficient's ball may hold 0.
    explicit CompensatedCertificate(const acb_poly_t polynomial);
    bool available() const { return !terms_.empty(); }
    int64_t degree() const { return static_cast<int64_t>(terms_.size()) - 1; }
    // The coefficients are held as a_k 2^-scale(): the evaluation's values and bounds are in those units.
    int64_t scale() const { return scale_; }
    // One root in a rectangle about the Newton step from the binary64 point nearest `point`, of the goal's half-width
    // (goal_half_width) or a smaller one (prove), published at `precision` and meeting the goal; false, with `root`
    // untouched, where an
    // admission or a test fails, so that the caller's ball certificate runs. It asks may_certify first, then prove;
    // where prove's test fails in Newton's quadratic range, both again from the binary64 point nearest the Newton step
    // (heuristics::compensated_newton_retries). Where `seen` is given it receives the evaluation at that first binary64
    // point, or one whose ok is false where none was made there.
    bool certify(Root &root, const acb_t point, const cps_options &options, slong precision,
                 CompensatedEvaluation *seen = nullptr) const;
    // Whether Rouché's test can pass at (x_re, x_im) at all, from the coefficients' errors alone, before the
    // evaluation: false only where sum_k error_k |x|^k, a part of the value's bound, already reaches 2 R U (U a bound
    // on the computed slope's modulus, R one on the half-width of every rectangle the test could accept there), so
    // that the test must fail. Where binary64 cannot hold the coefficients (a span beyond its range: what the words
    // leave out goes into the errors) the evaluation would be spent for nothing at every point. True where no
    // coefficient carries an error, or where an absolute tolerance sets the half-width.
    bool may_certify(double x_re, double x_im, const cps_options &options) const;
    // Whether some coefficient carries an error (may_certify's test runs only then).
    bool carries_errors() const { return errors_; }
    // The evaluation at (x_re, x_im) and Rouché's test, without may_certify (certify's second step; the tests check
    // that may_certify never refuses a point this proves). Where the test fails and `retry` is given, it is set where
    // the evaluation shows Newton's iteration from x in its quadratic range; where `seen` is given it receives the
    // evaluation.
    bool prove(Root &root, double x_re, double x_im, const cps_options &options, slong precision,
               CompensatedRetry *retry = nullptr, CompensatedEvaluation *seen = nullptr) const;
    // The evaluation at (x_re, x_im), its bounds rigorous (for the tests and certify).
    CompensatedEvaluation evaluate(double x_re, double x_im) const;
    // An upper bound on |g''| over the disk of radius `reach` about a point of modulus at most `modulus`, for every g
    // in the scaled coefficient balls (infinite where the majorant leaves binary64's range).
    double second_derivative_bound(double modulus, double reach) const;
    // Proposal 1 (dev/docs/cluster-double-double-design-20261027.md): simultaneous Aberth steps in double-double for
    // the points `moving` (indices into the n `points`), every other point held where it is; p and p' by Horner in
    // double-double from the two-word coefficients. A proposal only, nothing proved: certify proves. Where binary64
    // leaves a root of a small cluster at about u / s^(k-1) (s the cluster's size, k its roots), double-double
    // reaches about u^2 / s^(k-1), close enough for certify's binary64 point. Each point is first turned by its own
    // tiny angle (heuristics::double_double_nudge), so that a pair of close real roots that binary64 placed as a
    // conjugate pair can be reached: a real polynomial's iteration keeps a conjugate pair conjugate. A point stops when
    // its correction is below 2^-(goal_bits + margin) of it, or when the reach test finds its cluster too tight for
    // that (heuristics::double_double_reach_margin_bits); the repair stops when every point has, or after max_steps.
    // Returns whether a point
    // moved; a point whose step is not finite stays where it was. Adds to `evaluated` the steps it took at a point
    // (a converged point takes no more), the unit its caller charges.
    bool refine(acb_ptr points, int64_t n, const std::vector<int64_t> &moving, int64_t goal_bits, uint64_t &evaluated,
                unsigned max_steps) const;

  private:
    // a_k 2^-scale_ = (re_hi + re_lo) + i (im_hi + im_lo) + delta, |delta| <= error, for every a_k in its ball.
    struct Term {
        double re_hi = 0, im_hi = 0, re_lo = 0, im_lo = 0, error = 0;
    };
    std::vector<Term> terms_;
    // k times an upper bound on |a_k 2^-scale_|: the coefficients of the majorant of p' (may_certify).
    std::vector<double> first_;
    // k (k - 1) times an upper bound on |a_k 2^-scale_|: the coefficients of the majorant of p''.
    std::vector<double> second_;
    // Whether some coefficient carries an error (its ball's radius, or what its two words leave out): only then can
    // may_certify refuse a point.
    bool errors_ = false;
    int64_t scale_ = 0;
};
} // namespace certpolysolve::detail
