// A one-sided modular proof that an exact polynomial over Q(i) is square-free (squarefree_modular), and the resource
// admission that both modular square-free proofs share (modular_proof_admitted).
//
// p = 998244353 is 1 mod 4, so it splits in Z[i]: i maps to imaginary_unit = 3^((p-1)/4), a square root of -1 in F_p
// (a static_assert checks it). Each coefficient is reduced from its compressed form significand x radix^exponent
// (reduce_number), the exponent modulo p - 1. If no denominator is divisible by p, the degree survives and the gcd of
// f and f' over F_p (FLINT's nmod_poly_gcd) is constant, f is square-free over Q(i): a repeated factor would reduce to
// a nonconstant common divisor. Every other outcome, Inconclusive, Cancelled or Refused (the admission declined), says
// nothing about f. modular_proof_admitted admits n coefficients when the exact effort allowance covers n and
// temporary_bytes covers n x words machine words (split_proof_words for this proof).
//
// State: none.
// Entry points: squarefree_modular (SolveState::Impl::preprocess in solver.cpp, after the real import or on the
// coefficients as written where it is refused; nonzero_part_squarefree_modular here;
// exact_complex_factors in exact_complex.cpp, as its first proof, or alone where the write-out is refused);
// modular_proof_admitted and split_proof_words (exact_complex.cpp). nonzero_part_squarefree_modular removes
// only the known zero factor before the same proof (exact_preprocessing.cpp and piecewise admission tests).
// Not here: the second proof, over F_p[i]/(i^2 + 1) with p = 2^31 - 1 (ModularSquarefreeTest, exact_complex.cpp).
// Design notes: docs/piecewise-selection.md ("Mathematical capability check"), docs/resource-policy.md (the
// admission), docs/preparation-selection-findings-20260924.md (why this proof runs first).
#include "input/coefficients/polynomial_input.hpp"
#include "numeric/exact/modular_squarefree.hpp"
#include <flint/nmod_poly.h>

namespace certpolysolve::detail {
namespace {
std::pair<size_t, size_t> nonzero_range(const DensePolynomial &input) {
    size_t length = input.coefficients.size(), start = 0;
    while (length && input.coefficients[length - 1].zero())
        --length;
    while (start + 1 < length && input.coefficients[start].zero())
        ++start;
    return {start, length};
}
} // namespace

bool nonzero_part_squarefree_modular(const DensePolynomial &input, const cps_options &options) {
    const auto [start, length] = nonzero_range(input);
    if (length <= start + 1)
        return false;
    return squarefree_modular(std::span(input.coefficients).subspan(start, length-start), options) ==
           ModularSquarefreeResult::Proved;
}


namespace {
constexpr uint64_t prime = 998244353;
constexpr uint64_t power_mod(uint64_t a, uint64_t exponent) {
    uint64_t result = 1;
    for (; exponent; exponent >>= 1, a = a * a % prime)
        if (exponent & 1) result = result * a % prime;
    return result;
}
constexpr uint64_t imaginary_unit = power_mod(3, (prime - 1) / 4);
static_assert(imaginary_unit * imaginary_unit % prime == prime - 1);

class ModularPolynomial {
    nmod_poly_t value_;
  public:
    ModularPolynomial() { nmod_poly_init(value_, static_cast<ulong>(prime)); }
    ~ModularPolynomial() { nmod_poly_clear(value_); }
    ModularPolynomial(const ModularPolynomial &) = delete;
    ModularPolynomial &operator=(const ModularPolynomial &) = delete;
    nmod_poly_struct *value() & { return value_; }
    const nmod_poly_struct *value() const & { return value_; }
    nmod_poly_struct *value() && = delete;
    const nmod_poly_struct *value() const && = delete;
};
bool reduce_number(uint64_t &out, const Number &number) {
    const uint64_t denominator = fmpz_fdiv_ui(fmpq_denref(number.significand.value()), prime);
    if (!denominator) return false;
    const uint64_t numerator = fmpz_fdiv_ui(fmpq_numref(number.significand.value()), prime);
    // Each power below is about sixty modular products (sixty divisions): it is skipped where it multiplies zero (a
    // zero part), and where it is 1 (an integer significand; an exponent that is 0 mod p - 1, by Fermat's little
    // theorem as the radix is prime to p). The residue is the same.
    if (!numerator) {
        out = 0;
        return true;
    }
    out = denominator == 1 ? numerator : numerator * power_mod(denominator, prime - 2) % prime;
    const uint64_t exponent = fmpz_fdiv_ui(number.exponent.value(), prime - 1);
    if (exponent)
        out = out * power_mod(number.radix, exponent) % prime;
    return true;
}
} // namespace

bool modular_proof_admitted(size_t size, uint64_t words, const cps_options &options) {
    return options.exact_effort_fits(size) &&
           resources::storage_fits(resources::multiply(resources::multiply(size, words), sizeof(ulong)),
                                   options.temporary_bytes);
}

ModularSquarefreeResult squarefree_modular(std::span<const Coefficient> input, const cps_options &options) {
    if (input.size() < 2) return ModularSquarefreeResult::Inconclusive;
    if (!modular_proof_admitted(input.size(), split_proof_words, options)) return ModularSquarefreeResult::Refused;
    ModularPolynomial polynomial, derivative, gcd;
    for (size_t j = 0; j < input.size(); ++j) {
        if (options.cancelled()) return ModularSquarefreeResult::Cancelled;
        if (!input[j].exact()) return ModularSquarefreeResult::Inconclusive;
        uint64_t real, imag;
        if (!reduce_number(real, input[j].real) || !reduce_number(imag, input[j].imag))
            return ModularSquarefreeResult::Inconclusive;
        const auto residue = static_cast<ulong>((real + imaginary_unit * imag % prime) % prime);
        nmod_poly_set_coeff_ui(polynomial.value(), static_cast<slong>(j), residue);
    }
    if (nmod_poly_degree(polynomial.value()) != static_cast<slong>(input.size()-1))
        return ModularSquarefreeResult::Inconclusive;
    nmod_poly_derivative(derivative.value(), polynomial.value());
    nmod_poly_gcd(gcd.value(), polynomial.value(), derivative.value());
    if (options.cancelled()) return ModularSquarefreeResult::Cancelled;
    return nmod_poly_degree(gcd.value()) == 0 ? ModularSquarefreeResult::Proved :
                                              ModularSquarefreeResult::Inconclusive;
}
} // namespace certpolysolve::detail
