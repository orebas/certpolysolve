#pragma once
#include <optional>
#include "input/coefficients/polynomial_input.hpp"
#include "core/query.hpp"
#include "numeric/exact/input_structure.hpp"

namespace certpolysolve::detail {
// A single experimental admission rule. These constants schedule work; they
// do not weaken certificates or constrain the public solver's input domain.
// The caps are integers, compared exactly with the estimate's counts (docs/same-routes-every-platform-note-20261021.md);
// an absent cap is no cap.
struct PiecewiseAdmissionPolicy {
    std::optional<uint64_t> maximum_envelope_bend = 16; // initial local-model accuracy, in bits
    std::optional<uint64_t> maximum_work_ratio = 128;   // coarse development cost calibration, a unit of degree x support
    uint64_t minimum_degree = 0;
    bool require_attempt_fits_call = false;
    // Refuse when the whole estimate reaches this many sweeps of the global
    // method, each priced at degree x (degree + support) coefficient
    // operations: work >= sweeps x degree x (degree + support). The
    // estimate stops walking its rings as soon as that holds, so its cost is
    // bounded by the bound and not by the radial width of the root window.
    std::optional<uint64_t> maximum_global_sweeps;
    // The last condition, a modular proof that the input is square-free (a gcd at the input's degree, which can be
    // most of the selector's time). Off, the admission is its estimate's verdict alone: the axis charts' attempt
    // asks only that (session.cpp), and its own exact preparation proves square-freeness or stops it.
    // History: docs/code-history-20261021.md, piecewise_admission.hpp 1.
    bool prove_squarefree = true;
};
struct PiecewiseAdmission {
    bool selected = false, squarefree_proved = false;
    const char *reason = "unsupported_query";
    uint64_t degree = 0, support = 0, support_stride = 0, zero_order = 0;
    uint64_t model_bits = 0, estimated_models = 0, elapsed_ns = 0;
    uint64_t estimated_attempt_work = 0;
    // The estimate: model construction and local solves over the rungs, in operations (exact, saturating).
    uint64_t work = 0;
    // Statistics in binary64: the envelope's bend in bits, and work / (degree x support).
    double envelope_bend = 0, work_ratio = 0;
};
PiecewiseAdmission admit_piecewise(const DensePolynomial &, const cps_options &, PiecewiseAdmissionPolicy = {});
// The same on the input object's terms (C40's second step), without the square-free proof: the policy's
// prove_squarefree must be false. The cancellation is polled once a term, not once a coefficient.
PiecewiseAdmission admit_piecewise(const cps_polynomial &, const cps_options &, PiecewiseAdmissionPolicy);
// The admission decided from the input's structure alone, or nothing when it needs the coefficients: the decision,
// the reason and the fields admit_piecewise reaches for an input of that structure, its cancellation poll included
// (C40's second step: the selector on the input object's terms).
std::optional<PiecewiseAdmission> admission_by_structure(const InputStructure &, const cps_options &,
                                                         PiecewiseAdmissionPolicy = {});
// The D7/D8 selected-or-fallback wrapper (PiecewiseSelectionState,
// solve_piecewise_selected) was removed with the private
// piecewise_selection route.
// History: docs/code-history-20261021.md, piecewise_admission.hpp 2.
} // namespace certpolysolve::detail
