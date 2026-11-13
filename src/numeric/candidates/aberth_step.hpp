#pragma once
// The operations of Aberth-Ehrlich's step that its kernels share (W8's S10): Aberth's correction N / (1 - N S), with N
// = p/p' Newton's correction and S = sum_{j != i} 1/(z_i - z_j) the interaction sum over the other points, and the
// interaction sum's terms. Each function is the sequence of operations its callers wrote at their sites before, in the
// same order, so every caller gets the bits it got before: tests/test_aberth_step.cpp runs each old expression against
// its function. What differs between the kernels by design stays at the kernels: Newton's evaluator, the order of the
// sweep (Gauss-Seidel or Jacobi), the guards and the stops.
//
// The kernels that take Aberth's step, and what each keeps at its site. Every one moves a point z_i by N / (1 - N S);
// the five binary64 kernels and the two multiprecision packets take it from the functions below, and the others
// differ in an operation, kept for the reason given. A note such as machine_candidates.cpp 4 is that entry of
// docs/code-history-20261021.md; the S10 map (section 2) has each
// kernel's lines.
// - MachineEngine<R>::advance (machine_candidates.cpp; binary64 and ScaledDouble): Gauss-Seidel, each point moved
//   before the next is summed. N from PreparedNewton in binary64 (the reversed polynomial above |z| = 1), from
//   PreparedJets in ScaledDouble. S by interaction_term, one accumulator, j ascending; or, in binary64 from
//   heuristics::machine_interaction_tree_degree points, a MachineInteractionTree's far field and its near terms.
//   aberth_correction, or logarithmic_aberth_correction where p' is zero or the exponent fields put N or N S beyond
//   binary64's range (without it a five-term input lost its packet: machine_candidates.cpp 4). A point stops at a
//   relative step below 2^-45 or at |p| below Horner's noise; in binary64 an N that is not finite ends the packet
//   (exhausted); a coincident point, a step that is not finite, or an invalid, divide-by-zero, overflow or (binary64,
//   unless allowed) underflow flag fails it.
// - SecularPoints<R>::advance (study_candidates.cpp; binary64 and ScaledDouble): Jacobi. N from the secular model (on
//   a node, its removable limit). S, in binary64 on the secular recipe while every part is within 2^+-400, in four
//   fixed-order lanes, conj(d) times 1 / max(|d|^2, DBL_MIN) (12.7 -> 2.8 ns a pair: study_candidates.cpp 2);
//   otherwise by interaction_term in one accumulator. aberth_denominator with a zero test, or aberth_pole_limit at a
//   pole of the model. A zero denominator, a zero S at a pole, a coincident point, a step that is not finite or a
//   flag refuses the packet; the model's stops (ModelStoppingState) or a largest relative step below 2^-45 end it.
// - NodeOffsets::refine (study_candidates.cpp; binary64): one Jacobi step, node-relative, at a fresh model's floor.
//   N the node's limit from two-part node differences; S by interaction_term; aberth_correction. Only steps between
//   2^(2 - stage) and 2^(2 - 53) of the node are taken; a step that is not finite or a flag refuses them all.
// - ProgramMachineStage::sweep (program_candidates.cpp; binary64): Jacobi, every point evaluated at the sweep's
//   start positions, then every correction. N = value / derivative from the batched program evaluator, or from its
//   range tier (BlockDouble) converted back. S transposed, conj(d) / |d|^2 by two divisions, each sum its own
//   accumulator in j order; or, from heuristics::machine_interaction_tree_degree points with
//   heuristics::program_machine_tree_active active, the program route's tree. aberth_correction. No flag window: the
//   stops Converged, Floor and Held, a patience, and the window rule (program_machine.hpp).
// - small_roots (program_recurrence_seeds.cpp; binary64): Gauss-Seidel on F(r, .) of a small degree. N by Horner,
//   skipped where p' is zero; S by interaction_term, coincident points skipped; aberth_correction, a step that is not
//   finite skipped. It stops at a squared relative step below 2^-100, or after 200 sweeps.
// - candidate_packet and secular_candidate_packet (candidates.cpp; acf): Gauss-Seidel, and Jacobi or Gauss-Seidel by
//   a flag with a precision per point. N from a dense or sparse Horner, the regenerated secular model, or a
//   MonicSource. S by add_interactions, after an InteractionTree's far field from heuristics::interaction_tree_degree
//   points (the secular packet's in Jacobi order only). aberth_delta; -1/S at a secular pole and the research build's
//   Durand-Kerner quotient stay at their sites.
// - The chart proposals (piecewise_proposals.cpp; std::complex values): not on these functions. Their quotient is
//   the same algorithm as Complex's (a private copy), and S adds every j, j = i as zero, with conj(d) times
//   1 / |d|^2; but 1.0 - N S in std::complex keeps -x.im as the imaginary part where Complex's 1 + (-x) gives +0 when
//   x.im is +0, so the sign of a zero differs. Newton's step where 1 - N S is zero; the step's size at most 1/2.
//   Their Smith division and squared sizes replaced __divdc3 and hypot (19.4 -> 11.4 s on hyperbolic-24576-1:
//   piecewise_proposals.cpp 1).
// - Deflated discovery (piecewise_rounds.cpp) and the separated-roots pass (program_certificate.cpp): not on these
//   functions. N and the step are in balls, S in balls or in binary64 std::complex (the C runtime's __divdc3).
// - binary64_seeds (candidates.cpp, the research build only): std::complex, __divdc3; not on these functions.
//
// Include this header from the library's .cpp files only (and from its test, built with the library's floating-point
// flags), never from a header that other code includes: a test built with other flags must not compile these functions.
#include "numeric/arithmetic/machine_arithmetic.hpp"
#include <cstddef>
#include <cstdint>
#include <flint/acf.h>
#include <span>

// Binary64 and scaled binary64, in machine::Complex<R> with R double or machine::ScaledDouble, whose operators fix
// every operation and its order (machine_arithmetic.hpp: a - b is a + (-b), a product is two sums of two products, a
// quotient is Smith's division), under the library's -ffp-contract=off: MachineEngine::advance
// (machine_candidates.cpp), SecularPoints::advance and NodeOffsets::refine (study_candidates.cpp),
// ProgramMachineStage::sweep (program_candidates.cpp) and small_roots (program_recurrence_seeds.cpp).
namespace certpolysolve::detail::machine {
// 1 - N S: the product N S, then 1 + (-(N S)). For a caller that refuses a zero denominator before it divides.
template <class R> inline Complex<R> aberth_denominator(Complex<R> newton, Complex<R> sum) {
    return Complex<R>(R(1)) - newton * sum;
}
// Aberth's correction N / (1 - N S). No test of its own: a zero denominator gives a result that is not finite (and the
// division-by-zero flag), which every caller rejects.
template <class R> inline Complex<R> aberth_correction(Complex<R> newton, Complex<R> sum) {
    return newton / aberth_denominator(newton, sum);
}
// The same correction as 1 / (L - S), L = p'/p the logarithmic derivative: for MachineEngine where p' is zero or N or
// N S is out of binary64's range.
template <class R> inline Complex<R> logarithmic_aberth_correction(Complex<R> logarithmic_derivative, Complex<R> sum) {
    return Complex<R>(R(1)) / (logarithmic_derivative - sum);
}
// The correction's limit -1/S where N is infinite: for SecularPoints at a pole of the secular model.
template <class R> inline Complex<R> aberth_pole_limit(Complex<R> sum) { return -Complex<R>(R(1)) / sum; }
// A term of the interaction sum, 1/d for d = z_i - z_j by Smith's division. The loop over j, its order and what it does
// at a coincident point stay with the callers.
template <class R> inline Complex<R> interaction_term(Complex<R> difference) { return Complex<R>(R(1)) / difference; }
} // namespace certpolysolve::detail::machine

// Multiprecision, in FLINT's acf (arf parts, no radius) rounded to nearest at the point's precision: the packets
// candidate_packet and secular_candidate_packet (candidates.cpp).
namespace certpolysolve::detail {
// Adds to `sum` the interaction sum of z, the terms 1/(z - z_j) over the neighbours j != i in the neighbours'
// order: each difference z - z_j rounded to nearest at prec, FLINT's approximate inverse, then the sum. False at a
// point equal to z, the sum then holding the terms before it. `points` is the caller's vector of points, each an acf
// `value`; `difference` and `inverse` are scratch; `interactions` counts the terms added.
template <class Points>
inline bool add_interactions(acf_t sum, const acf_t z, const Points &points, std::span<const size_t> neighbours,
                             size_t i, slong prec, acf_t difference, acf_t inverse, uint64_t &interactions) {
    for (const size_t j : neighbours)
        if (j != i) {
            acf_sub(difference, z, points[j].value, prec, ARF_RND_NEAR);
            if (arf_is_zero(acf_realref(difference)) && arf_is_zero(acf_imagref(difference)))
                return false;
            acf_approx_inv(inverse, difference, prec, ARF_RND_NEAR);
            acf_add(sum, sum, inverse, prec, ARF_RND_NEAR);
            ++interactions;
        }
    return true;
}
// Aberth's correction delta = N / (1 - N S), rounded to nearest at prec: the product N S, its negation, 1 added to its
// real part, then FLINT's approximate division. False where 1 - N S is zero, delta then untouched. `denominator` is
// scratch, and neither `newton` nor `sum` may be it.
inline bool aberth_delta(acf_t delta, const acf_t newton, const acf_t sum, acf_t denominator, slong prec) {
    acf_mul(denominator, newton, sum, prec, ARF_RND_NEAR);
    acf_neg(denominator, denominator);
    arf_add_ui(acf_realref(denominator), acf_realref(denominator), 1, prec, ARF_RND_NEAR);
    if (arf_is_zero(acf_realref(denominator)) && arf_is_zero(acf_imagref(denominator)))
        return false;
    acf_approx_div(delta, newton, denominator, prec, ARF_RND_NEAR);
    return true;
}
} // namespace certpolysolve::detail
