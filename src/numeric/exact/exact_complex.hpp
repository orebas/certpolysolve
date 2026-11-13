#pragma once
#include "core/query.hpp"
#include "core/engine_policy.hpp"
#include <span>

namespace certpolysolve::detail {
using GaussianPolynomial = std::vector<GaussianRational>;
struct GaussianFactor { GaussianPolynomial polynomial; uint64_t multiplicity = 1; };
struct GaussianDecomposition {
    std::vector<GaussianFactor> factors;
    std::string status;
    bool cancelled = false;
    // No factors, but the split prime's test proved the coefficients as written square-free: the compact-exponent
    // rule refused their write-out (exact_complex_factors). The caller keeps the ball path with that proof.
    bool squarefree_proved = false;
    bool available() const { return !factors.empty(); }
};

// Exact, transactional squarefree decomposition over Q(i). A degree/work/storage
// guard or cancellation returns no factors; callers retain the ball path.
// Normal inputs first use a one-sided modular squarefreeness certificate: the
// split prime's (modular_squarefree.hpp), then the extension field's modulo
// 2^31 - 1, both admitted by modular_proof_admitted and neither charged beyond
// that admission (a per-stage estimate, docs/resource-policy.md). Below the
// coefficient count the exact_operations allowance admits no stage at all and
// the decomposition is refused before any runs.
GaussianDecomposition exact_complex_factors(std::span<const Coefficient> coefficients,
                                            const cps_options& options, const EnginePolicy *trace = nullptr);
// Words an element of the extension-field proof takes (its admission): an
// fq_nmod element, a polynomial of two residues with its header, in each of
// its three polynomials.
extern const uint64_t extension_proof_words;
void enclose_gaussian(acb_poly_t out, const GaussianPolynomial& polynomial, slong precision);
}
