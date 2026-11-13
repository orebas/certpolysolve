#pragma once
// The all-real seed pass (docs/real-line-seed-pass-design-20261020.md): starting points on the real line for an exact
// real polynomial whose roots are all real, from one interpolant of f = p 2^-phi on the interval of its roots
// (prototype v4). The seeds only start an iteration: the certified route decides every answer, so a wrong seed costs
// time and never an answer, and every failure of the pass is "no seeds", never an error. History:
// docs/code-history-20261021.md, real_line_seeds.hpp 1.
//
// The same bits on every CPU: every value that decides anything (a sign, a count, a precision, a node count, a seed)
// comes from FLINT ball arithmetic, from integer arithmetic, or from a fixed sequence of IEEE additions,
// subtractions, multiplications, divisions and comparisons (the library is built with -ffp-contract=off). No libm
// transcendental is called in binary64: cosines, exponentials, logarithms and the transforms' twiddles come from arb
// and are rounded to binary64 once.
#include "core/raii.hpp"
#include <cstdint>
#include <functional>
#include <vector>

namespace certpolysolve::detail {

// Why the pass ended. Seeds is the only outcome with seeds; every other one leaves the caller's own start.
enum class RealLineSeedStop : int {
    Seeds,
    Degree,             // degree below 2, or a zero constant term (the caller passes a factor without the root 0)
    NewtonInequality,   // the coefficients violate Newton's inequalities: p has a non-real root
    LaguerreInequality, // (n-1) p'^2 < n p p'' certainly at an evaluated point: p has a non-real root
    LaguerreStalled,    // Laguerre's iteration toward an extreme root did not converge within its steps
    FewCells,           // the certified sign changes at the nodes are too few against the degree
    NoisyInterpolant,   // the interpolant changes sign more often than p can (more than its degree): it is noise
    WideRange,          // log2 |p| ranges too widely over the interval for the nodes (before the change of basis)
    Memory,             // the predicted storage of the pass exceeds the limit
    Precision,          // a precision the pass needs exceeds the limit
    Cancelled,          // the caller's poll said stop (cancellation or the time limit)
    Numerical,          // a value that cannot be used (not finite, an empty interval)
};
const char *real_line_seed_stop_name(RealLineSeedStop stop);

struct RealLineSeedLimits {
    // No evaluation, change of basis or transform above this many bits (the call's work-bits cap).
    slong max_precision = 1 << 20;
    // The pass's own storage (balls of the nodes and of the change of basis), predicted before it is allocated.
    uint64_t max_bytes = uint64_t(4) << 30;
    // Polled between stages and inside every loop over nodes, seeds or steps. May be empty.
    std::function<bool()> cancelled;
    // The certified sign changes at the nodes must number at least this share of the degree, or the pass ends
    // there (FewCells): the check that bounds what a polynomial with many non-real roots costs. Zero disables it
    // (the measurement of its threshold).
    double min_cell_fraction = 0;
    // The per-seed tests of seeds that share a cell (two ball evaluations each). Off: such seeds are kept untested,
    // up to the degree (the measurement of the tests' worth, design note "Fewer seeds than roots").
    bool seed_tests = true;
};

struct RealLineSeedStatistics {
    RealLineSeedStop stop = RealLineSeedStop::Seeds;
    slong degree = 0;
    slong newton_failure = -1;        // the first k where Newton's inequality fails, or -1
    uint64_t laguerre_steps = 0;      // both ends
    uint64_t ball_evaluations = 0;    // of p, p' or p'' at one point
    slong largest_evaluation_precision = 0; // the most bits one evaluation needed
    double lower_root = 0, upper_root = 0;   // the extreme roots Laguerre's method reached
    double range_bits = 0;            // log2 |p| over the interval, max - min of the samples
    double cancellation_bits = 0;     // what the change of basis loses to the monomial basis (its bound)
    slong conversion_precision = 0;   // P_conv
    slong transform_precision = 0;    // P_dct
    slong nodes = 0;                  // N (the last one, after any doubling)
    uint64_t predicted_bytes = 0;
    slong envelope_terms = 0;
    double envelope_raise_bits = 0;
    double tail_bits = 0;             // log2 of the interpolant's tail against its largest coefficient
    slong cells_with_sign_change = 0;
    slong interpolant_seeds = 0, touching_seeds = 0, node_seeds = 0;
    slong alone = 0, tested = 0, certified_by_test = 0, unverified = 0, dropped = 0;
    uint64_t interval_ns = 0, range_ns = 0, conversion_ns = 0, transform_ns = 0, envelope_ns = 0, scale_ns = 0,
             roots_ns = 0, cells_ns = 0, total_ns = 0;
};

struct RealLineSeeds {
    std::vector<double> seeds; // increasing, at most the degree; empty unless statistics.stop == Seeds
    RealLineSeedStatistics statistics;
};

// The pass on p (integer coefficients, degree at least 2). A root of p at zero is divided out exactly and seeded
// by the point 0 itself.
RealLineSeeds real_line_seeds(const fmpz_poly_struct *p, const RealLineSeedLimits &limits);

// The parts below are exposed for the tests.

// Newton's inequalities for a polynomial with only real roots, with the binomials divided out:
// a_k^2 k (n - k) >= a_{k-1} a_{k+1} (k + 1) (n - k + 1) for 0 < k < n, and a zero a_k only between coefficients of
// opposite signs. Returns -1 when every one holds, else the first k that fails. Exact.
slong newton_inequality_failure(const fmpz_poly_struct *p);

// A radix-2 discrete Fourier transform in binary64, in a fixed order, with twiddles from arb: built once for the
// largest length a pass needs and applied to any power of two up to it.
class Binary64Transform {
  public:
    explicit Binary64Transform(slong max_length); // a power of two, at least 2
    slong max_length() const { return max_length_; }
    // In place: (re, im) of length L (a power of two <= max_length) becomes X_k = sum_j x_j exp(-2 pi i j k / L).
    void forward(std::vector<double> &re, std::vector<double> &im) const;
    // DCT-I: y_j = x_0 + (-1)^j x_M + 2 sum_{k=1}^{M-1} x_k cos(pi j k / M), j = 0..M, M = x.size() - 1 a power of
    // two with 2 M <= max_length (the transform of the even extension of x, length 2 M).
    std::vector<double> dct1(const std::vector<double> &x) const;
    // cos(pi k / m) as held in the twiddle table, for m a power of two with 2 m <= max_length and 0 <= k <= m.
    double cos_pi(slong k, slong m) const;

  private:
    slong max_length_;
    std::vector<double> cos_, sin_; // cos(2 pi k / max_length), sin(...), 0 <= k < max_length / 2
};

// The Chebyshev coefficients of the interpolant through values at t_j = cos(pi j / M), j = 0..M, and back.
std::vector<double> chebyshev_coefficients(const Binary64Transform &t, const std::vector<double> &values);
std::vector<double> chebyshev_values(const Binary64Transform &t, const std::vector<double> &coefficients, slong m);
// sum_k c_k T_k(x) by Clenshaw's recurrence, in a fixed order.
double chebyshev_evaluate(const std::vector<double> &c, slong length, double x);

} // namespace certpolysolve::detail
