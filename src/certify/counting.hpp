#pragma once
#include "core/cover.hpp"
#include "core/query.hpp"
#include <bit>
#include <cmath>
#include <limits>
#include <optional>
#include <utility>
#include <vector>

namespace certpolysolve::detail {
struct EnginePolicy;

struct DiskCount {
    slong count = -1; // -1 is inconclusive, never an assertion of zero roots.
    unsigned graeffe_steps = 0;
    bool needs_precision = false; // Scheduling information, never a root count.
    // The direct selected-coefficient or zero-exclusion comparison, before
    // Graeffe, was undecided because of arithmetic width. This permits a
    // rotated-shift retry; it is distinct from the all-index precision request.
    bool width_limited = false;
    // The test stopped before an operation its context's allowance did not
    // cover (CountContext::allowance): nothing was decided. Calling again with
    // the same arguments resumes it where it stopped.
    bool suspended = false;
    // The direct zero-exclusion comparison was width-obstructed. Subdivision
    // uses this narrower finite-cap stop heuristic: uncertainty about a
    // positive count need not prevent a useful split into smaller cells.
    bool exclusion_limited = false;
};

struct CountContext {
    certpolysolve::detail::SearchStatistics &statistics;
    bool adaptive = false;
    bool real_shift = false; // EnginePolicy::pellet_real_shift (pellet_count)
    bool research_identity_shift_meter = false;
    // Borrowed only from the owning SolveState, installed before each step.
    const EnginePolicy *research_service_policy = nullptr;
    // Operations this context's Taylor shifts and Graeffe steps cost by
    // FLINT's block structure (dovetail.hpp), unweighted; the dovetail meter
    // charges them at the precision's weight.
    double metered = 0;
    // Dovetailing gives a test a work allowance (dovetail.hpp): no Taylor shift
    // or Graeffe step starts that would take `metered` past it. Each such
    // operation's cost is known before it runs, from its operand's block
    // structure, so a test that would need one it is not allowed stops
    // (DiskCount::suspended) with what it has computed. This bounds the
    // transforms charged here; other work in the caller has its own costs.
    // A single disk test is not: one test can need many times the
    // iteration's whole work under its share. Infinite, the default, never
    // stops a test. History: docs/code-history-20261021.md, counting.hpp 1.
    double allowance = std::numeric_limits<double>::infinity();
    // Whether the last pellet_count call suspended, and the cost (unweighted)
    // of the operation it waits for.
    bool suspended = false;
    double pending = 0;
    // The caller's cancellation, polled where the allowance is checked: before
    // every Taylor shift and Graeffe step. A cancelled test suspends as one
    // past its allowance does, so a step of the complete method returns within
    // one such operation of the caller's stop. Without it one test runs to its
    // end, however far past the caller's deadline. Null (the default) never
    // stops a test. History: docs/code-history-20261021.md, counting.hpp 2.
    const cps_options *cancellation = nullptr;
    bool bounded() const { return std::isfinite(allowance); }
    bool covers(double cost) const {
        return metered + cost <= allowance && !(cancellation && cancellation->cancelled());
    }
    // A test's arguments, the key of a suspended or completed test.
    struct TestKey {
        ComplexPolynomial source;
        ComplexBall center;
        RealBall radius; // exact: the test's radius as its midpoint
        slong precision = 0, expected = -1;
        unsigned max_graeffe = 0;
    };
    // The suspended test: its stage (0: on the direct shift; 1: on the rotated
    // real shift's retry, the direct result kept) and, once the stage's shift
    // ran, the scaled shifted polynomial after `step` Graeffe steps.
    struct Suspension {
        TestKey key;
        int stage = 0;
        bool shifted_ready = false, width_limited = false, exclusion_limited = false;
        unsigned step = 0;
        ComplexPolynomial shifted;
        DiskCount direct;
    };
    std::optional<Suspension> suspension;
    // Tests completed under an allowance, so that a step resumed after a later
    // test of the same certificate suspended does not compute them again. A
    // step runs at most three tests (a component's count, then the inner and
    // outer tests of its certificate), so three entries cover any one step.
    static constexpr size_t completed_tests = 3;
    std::vector<std::pair<TestKey, DiskCount>> completed;
    explicit CountContext(certpolysolve::detail::SearchStatistics& s, bool adapt = false) : statistics(s), adaptive(adapt) {}
};

// S(x) = p(u (x + |c|)) = p(c + u x) with u = c/|c|: p shifted to the centre
// and rotated, with every coefficient modulus of p(c + y) and every disk
// count about c unchanged. The rotation costs one rounding per coefficient
// and the shift is by a real number, so Arb's rectangular balls are not
// re-boxed at every step as they are by an off-axis centre.
void rotated_real_shift(acb_poly_t out, const acb_poly_t polynomial, const acb_t center, slong precision);

// Uniform count in |z-center| < radius, including multiplicity. The strict
// dominance certificate also excludes roots on the boundary. Center and radius
// must be finite exact dyadics, radius positive, and the degree must be fixed.
// `expected` restricts the count being tested; -1 asks for any provable count.
// A disk test is one Taylor shift and at most this many Graeffe iterations.
inline constexpr unsigned default_graeffe_steps = 4;
// The subdivision's disk tests allow more. Becker, Sagraloff, Sharma and Yap's test (arXiv:1509.06231, Lemma 1;
// reference-library/text/complex_isolation2015.txt) applies N = ceil(log2(1 + log2 n)) + 5 root-squarings, after
// which it succeeds on every disk whose k roots lie within 2 sqrt(2) / 3 of its radius and the others beyond 4/3:
// a separation ratio of about 1.41, constant in n. Keep the existing allowance
// of nine below degree 32,769; thereafter use the degree-dependent bound.
// The predicate stops at the first decisive iterate.
inline constexpr unsigned subdivision_graeffe_steps(slong degree) {
    constexpr unsigned established_allowance = 9;
    if (degree <= 1)
        return established_allowance;
    // ell = ceil(log2 n). N changes only when log2 n crosses an integer
    // of the form 2^k - 1, so rounding that inner logarithm up is exact here.
    const unsigned ell = static_cast<unsigned>(std::bit_width(static_cast<ulong>(degree - 1)));
    const unsigned bound = 5 + static_cast<unsigned>(std::bit_width(ell));
    return bound > established_allowance ? bound : established_allowance;
}
DiskCount pellet_count(const acb_poly_t polynomial, const acb_t center, const arf_t radius, slong precision,
                       unsigned max_graeffe = default_graeffe_steps, slong expected = -1,
                       CountContext *context = nullptr);

// Certify both an inner and an outer disk with the same count, sandwiching the
// EXPORTED rectangle. A disk count by itself does not establish a square count.
// `factor_multiplicity` is known from exact squarefree preprocessing, or 1.
// The output is changed only on success; count > 1 does not prove multiplicity.
bool certify_cluster_box(Root &out, const acb_poly_t polynomial, const acb_t center, const arf_t radius,
                         slong count, uint64_t factor_multiplicity, slong precision,
                         unsigned max_graeffe = default_graeffe_steps, CountContext *context = nullptr,
                         bool *needs_precision = nullptr);

// Bounded, candidate-guided recovery, not a subdivision/completeness algorithm.
// Hints can be wrong or duplicated: only count predicates can certify output.
// Appends a full disjoint cover of this factor on success; otherwise no change.
// If further precision is available, propagate an inconclusive arithmetic
// comparison so the source owner can re-enclose before trying more geometry.
bool recover_clusters(std::vector<Root> &out, const acb_poly_t polynomial, acb_srcptr hints, slong hint_count,
                      uint64_t factor_multiplicity, const cps_options &options, slong precision,
                      CountContext *context = nullptr, bool *needs_precision = nullptr);

} // namespace certpolysolve::detail
