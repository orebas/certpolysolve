#pragma once
#include "input/coefficients/polynomial_input.hpp"
#include "core/query.hpp"

namespace certpolysolve::detail {
// One-sided proof over F_p, p = 998244353, i mapped to a square root of -1.
// Includes zero roots: callers may pass an already deflated nonzero part, but
// this function never silently removes a repeated zero factor. Degree loss,
// bad denominators, uncertainty or cancellation are inconclusive: nothing
// about the input follows. A resource refusal (modular_proof_admitted) is
// Refused: the allowances did not admit the proof, which is no evidence
// either, and no costlier proof may stand in for it. Compressed
// binary/decimal exponents need no expansion.
// History: docs/code-history-20261021.md, modular_squarefree.hpp 1.
enum class ModularSquarefreeResult { Inconclusive, Proved, Cancelled, Refused };
ModularSquarefreeResult squarefree_modular(std::span<const Coefficient>, const cps_options &);
// One-sided proof for the exact nonzero part after removing its known power of x.
// Repeated zero roots are deliberately excluded; false proves nothing about that part.
bool nonzero_part_squarefree_modular(const DensePolynomial &, const cps_options &);
// The admission of a modular square-free proof over `size` coefficients, one
// rule for both residue fields (docs/resource-policy.md): exact_operations
// admits one unit an imported element, and temporary_bytes the proof's
// residue polynomials at `words` machine words an element.
bool modular_proof_admitted(size_t size, uint64_t words, const cps_options &);
// Words an element of the split-prime proof takes: one residue in each of its
// three polynomials (the input, its derivative, their gcd) and one of slack.
inline constexpr uint64_t split_proof_words = 4;
} // namespace certpolysolve::detail
