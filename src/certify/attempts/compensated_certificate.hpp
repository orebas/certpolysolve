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
#include "certify/attempts/compensated_kernel.hpp"
#include "core/cover.hpp"
#include "core/query.hpp"
#include <flint/acb_poly.h>
#include <cstdint>
#include <span>
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
// What prove derives from an evaluation before its tests (CompensatedCertificate::geometry): Newton's step q with
// |q| from above and the bound on its rounding, |D| from below, the centre z1 = x - q as two words a part and as a
// ball, the goal's half-width r about z1 from below, and the reach |q| + 2r of the larger circle from x. ok is false
// where prove declines before its tests.
struct CompensatedGeometry {
    bool ok = false;
    double step = 0, step_error = 0, slope_lower = 0;
    double z_re_hi = 0, z_re_lo = 0, z_im_hi = 0, z_im_lo = 0;
    ComplexBall center;
    double r = 0, reach_out = 0;
};
// What CompensatedCertificate::evaluate_ahead found at a point, for certify: nothing (certify decides alone), that
// may_certify refuses the binary64 point (x_re, x_im) nearest it, or the evaluation there with what prove derives from
// it at `precision` before its tests: the geometry and, where that is ok, the bound on |g''| over the larger circle's
// disk.
struct CompensatedAhead {
    enum class State : uint8_t { None, Refused, Evaluated };
    State state = State::None;
    double x_re = 0, x_im = 0;
    slong precision = 0;
    CompensatedEvaluation evaluation;
    CompensatedGeometry geometry;
    double second = 0;
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
    // point, or one whose ok is false where none was made there. Where `ahead` is given and is for that first binary64
    // point and this precision (evaluate_ahead), its refusal stands for may_certify's answer, and its evaluation,
    // geometry and bound on |g''| for what prove computes first: the same values, made before.
    bool certify(Root &root, const acb_t point, const cps_options &options, slong precision,
                 CompensatedEvaluation *seen = nullptr, const CompensatedAhead *ahead = nullptr) const;
    // certify's tries one at a time, for a caller that makes the later tries of many points together (the attempt's
    // first pass): certify is certify_first, then certify_again from each `retry` while one is available, at most
    // heuristics::compensated_newton_retries times. certify_first is certify's first try at the binary64 point
    // nearest `point`; certify_again its try at `from`, the point the try before asked for. Each sets `retry` to the
    // point of the next try, or to none. `ahead`, where given and made for the try's own point and this precision,
    // stands for what the try computes first (evaluate_ahead, evaluate_ahead_at).
    bool certify_first(Root &root, const acb_t point, const cps_options &options, slong precision,
                       CompensatedEvaluation *seen, const CompensatedAhead *ahead, CompensatedRetry &retry) const;
    bool certify_again(Root &root, const CompensatedRetry &from, const cps_options &options, slong precision,
                       const CompensatedAhead *ahead, CompensatedRetry &retry) const;
    // What certify would compute first at each of points[indices[m]], made ahead for all of them and together, into
    // ahead[m]: the binary64 point nearest the point, may_certify there, and where it admits the point the evaluation,
    // prove's geometry at `precision` and its bound on |g''|. Each is the value certify and prove compute one point
    // at a time; made together, the evaluations run four points at a time and the majorants' recurrences sixteen
    // (compensated_kernel.hpp, platform/kernel_variant.hpp). A point that is not finite or is zero, and every point
    // where the tier is unavailable or the goal beyond it, gets State::None. The attempt's first pass asks every
    // point in turn; this changes none of its results.
    void evaluate_ahead(acb_srcptr points, std::span<const slong> indices, const cps_options &options, slong precision,
                        std::span<CompensatedAhead> ahead) const;
    // The same at the binary64 points (x_re[m], x_im[m]) themselves, into ahead[m]: what a later try computes first
    // at the point the try before asked for (certify_again).
    void evaluate_ahead_at(std::span<const double> x_re, std::span<const double> x_im, const cps_options &options,
                           slong precision, std::span<CompensatedAhead> ahead) const;
    // Whether Rouché's test can pass at (x_re, x_im) at all, from the coefficients' errors alone, before the
    // evaluation: false only where sum_k error_k |x|^k, a part of the value's bound, already reaches 2 R U (U a bound
    // on the computed slope's modulus, R one on the half-width of every rectangle the test could accept there), so
    // that the test must fail. Where binary64 cannot hold the coefficients (a span beyond its range: what the words
    // leave out goes into the errors) the evaluation would be spent for nothing at every point. True where no
    // coefficient carries an error, or where an absolute tolerance sets the half-width.
    bool may_certify(double x_re, double x_im, const cps_options &options) const;
    // Whether some coefficient carries an error (may_certify's test runs only then).
    bool carries_errors() const { return errors_; }
    // Whether every polynomial in the coefficient balls has real coefficients: every ball's imaginary part is exactly
    // zero. The roots of each then come in conjugate pairs with their multiplicities, so the reflection in the real
    // axis of a rectangle that holds exactly one root holds exactly one root.
    bool real_coefficients() const { return real_; }
    // The evaluation at (x_re, x_im) and Rouché's test, without may_certify (certify's second step; the tests check
    // that may_certify never refuses a point this proves). Where the test fails and `retry` is given, it is set where
    // the evaluation shows Newton's iteration from x in its quadratic range; where `seen` is given it receives the
    // evaluation. `ahead`, where not null, holds evaluate(x_re, x_im), the geometry at `precision` and the bound on
    // |g''|, made before (evaluate_ahead): they are used in place of computing them.
    bool prove(Root &root, double x_re, double x_im, const cps_options &options, slong precision,
               CompensatedRetry *retry = nullptr, CompensatedEvaluation *seen = nullptr,
               const CompensatedAhead *ahead = nullptr) const;
    // The evaluation at (x_re, x_im), its bounds rigorous (for the tests and certify).
    CompensatedEvaluation evaluate(double x_re, double x_im) const;
    // The evaluations at the points (x_re[m], x_im[m]) into out[m], each what evaluate returns there: the checks and
    // the bounds are made a point, the loop between them (compensated_kernel.hpp) four points at a time where the CPU
    // has the instructions, a last group filled with copies of its last point.
    void evaluate_many(std::span<const double> x_re, std::span<const double> x_im,
                       std::span<CompensatedEvaluation> out) const;
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
    using Term = compensated_kernel::Term;
    // Before the loop at out's point: the checks, |x| from above and the amplification. False where the evaluation
    // declines (out.ok stays false).
    bool begin_evaluation(CompensatedEvaluation &out) const;
    // After the loop: the values and their bounds from the loop's sums.
    void finish_evaluation(CompensatedEvaluation &out, const compensated_kernel::Sums &sums) const;
    // prove's steps. geometry: from the evaluation to the larger circle's reach, before any test. finish_proof: the
    // tests on the circles, the retry's point where they fail in Newton's quadratic range, and the rectangle's
    // publication, from the evaluation, its geometry and second = second_derivative_bound(e.modulus, reach_out).
    CompensatedGeometry geometry(const CompensatedEvaluation &e, const cps_options &options, slong precision) const;
    bool finish_proof(Root &root, const CompensatedEvaluation &e, const CompensatedGeometry &g, double second,
                      const cps_options &options, slong precision, CompensatedRetry *retry) const;
    // second_derivative_bound's steps about its majorant's sum (degree at least 2). begin: Y >= |z| on the disk and
    // the amplification, false where the bound is infinite without the sum. finish: the bound from the sum at Y.
    bool second_derivative_begin(double modulus, double reach, double &Y, double &amplification) const;
    double second_derivative_finish(double sum, double amplification) const;
    // may_certify's steps about its sums. begin: X >= |x| and the amplification, false where the answer is true
    // without the sums. finish: the answer from S = sum_k error_k X^k, slope = sum_k k A_k X^(k - 1) and
    // K0 = second_derivative_bound(X, 0).
    bool may_certify_begin(double x_re, double x_im, const cps_options &options, double &X, double &amplification) const;
    bool may_certify_finish(double X, double amplification, double S, double slope, double K0,
                            const cps_options &options) const;
    // may_certify at the points (x_re[m], x_im[m]) into admitted[m], and a majorant's recurrence at the points x[m]
    // into out[m]: the values of one point at a time, the recurrences made compensated_kernel::majorant_group points
    // at a time where the CPU has the instructions.
    void may_certify_many(std::span<const double> x_re, std::span<const double> x_im, const cps_options &options,
                          std::span<uint8_t> admitted) const;
    static void majorants(const double *coefficients, uint64_t top, uint64_t bottom, std::span<const double> x,
                          std::span<double> out);
    // evaluate_ahead's steps after the points are known, for at most heuristics::compensated_ahead_points of them:
    // ahead[asked[k]] holds the binary64 point (x_re[k], x_im[k]) and the precision; may_certify, the evaluations,
    // the geometries and the bounds on |g''| are made for all of them together. x_re and x_im are overwritten.
    void ahead_group(const size_t *asked, double *x_re, double *x_im, size_t asking, const cps_options &options,
                     slong precision, std::span<CompensatedAhead> ahead) const;
    std::vector<Term> terms_;
    // terms_[k].error, in an array of their own: the coefficients of may_certify's first sum.
    std::vector<double> error_;
    // k times an upper bound on |a_k 2^-scale_|: the coefficients of the majorant of p' (may_certify).
    std::vector<double> first_;
    // k (k - 1) times an upper bound on |a_k 2^-scale_|: the coefficients of the majorant of p''.
    std::vector<double> second_;
    // Whether some coefficient carries an error (its ball's radius, or what its two words leave out): only then can
    // may_certify refuse a point.
    bool errors_ = false;
    // Whether every coefficient ball's imaginary part is exactly zero (real_coefficients).
    bool real_ = false;
    int64_t scale_ = 0;
};
} // namespace certpolysolve::detail
