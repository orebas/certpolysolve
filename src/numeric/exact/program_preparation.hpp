#pragma once
// Program input, milestone 2 (docs/slp-input-spec-20261015.md, sections 4, 5 and 7): what a solve needs before it
// starts, without expanding the program over the integers or the rationals. The proved degree and leading
// coefficient, the modular image's verdict on square-freeness (under a cost rule), and a bound of the roots'
// moduli; or a refusal with its diagnostic.
#include "core/query.hpp"
#include "input/program/program.hpp"
#include <string>

namespace certpolysolve::detail {

// The image modulo a word prime, priced before it runs (section 5's cost rule). Word operations: an addition or a
// negation its length, a product its length times the length's bit count (an FFT's), a power twice its last
// squaring's, and the square-free test a half-gcd of the output and its derivative (a product's operations times
// the bit count again). Words: the most held at once under the interpreter's liveness, plus the test's
// derivative and gcd. Saturates at UINT64_MAX.
struct ModularCost {
    uint64_t operations = 0, words = 0;
};
ModularCost modular_cost(const Program &program);

// The modular image's verdict (section 5): Proved when, modulo a prime that keeps the degree, the image and its
// derivative are coprime; Unknown otherwise (over the cost rule, two primes not square-free, or no usable prime).
// Never a proof of a repeated root.
struct ModularVerdict {
    bool squarefree = false;
    bool cancelled = false;
    std::string note; // for the diagnostics and the tests
};
ModularVerdict modular_squarefree_image(const Program &program, uint64_t degree, const cps_options &options);

// A bound of every root's modulus from the program's majorant M (section 7), for a polynomial of degree n >= 1 and
// leading coefficient a_n: the smaller of Cauchy's bound M(1) / |a_n| and the least power of two R found with
// M(R) < 2 |a_n| R^n. Both hold for any program, also one whose leading terms cancel (where the second may not
// exist): M(r) bounds the sum of |a_i| r^i of the expanded polynomial at every r.
Magnitude program_root_bound(const Program &program, uint64_t degree, const GaussianRational &leading);

// Whether every root is zero, the program computing a_n x^n: proved when the majorant test holds at R = 2^-(D+L+1),
// with D an upper bound of the bits of a common denominator of the coefficients (from the program's constants, up a
// product's and a sum's operands and a power's exponent) and |a_n| < 2^L. The test gives |a_i| < |a_n| R < 2^-D for
// every i < n, where a nonzero a_i has |a_i| >= 2^-D. False whenever that is not shown: the bound D too large to
// hold, or a majorant too loose (a program computing x^2 as (x+1)^2 - 2x - 1); never a proof of a nonzero root.
bool program_monomial(const Program &program, uint64_t degree, const GaussianRational &leading);

// The order v >= 1 of the program's root at zero (p = x^v q, q(0) != 0), proved; 0 when there is none or it is not
// decided. The lowest coefficients of p come from the program run over truncated power series at zero with complex
// ball coefficients, at precisions doubling from 64 bits. A coefficient is decided zero when its ball lies within
// 2^-D of zero, D an upper bound of the bits of a common denominator of the coefficients (a nonzero one is at least
// 2^-D), and nonzero when its ball excludes zero. v is the first nonzero coefficient's index with every one below it
// zero. Undecided coefficients at the precision ceiling, or an order past the width limit, give 0: no claim. A program
// whose value at zero is clearly nonzero, the common case, costs one evaluation at 64 bits.
uint64_t program_zero_order(const Program &program, uint64_t degree);

struct ProgramPreparation {
    // The refusals (section 4), and the resource rule of section 3.1 (the length against exact_operations, a
    // constant against the exact storage limit).
    enum class Refusal { None, NoOutput, Length, ConstantStorage, ZeroPolynomial, LeadingUndetermined };
    Refusal refusal = Refusal::None;
    std::string diagnostic; // the refusal's text; empty when the program is ready
    uint64_t degree = 0;
    size_t top_width = 0; // the width of the top series that decided the degree
    GaussianRational leading;
    // Degree 1 or more only: a constant has no roots to bound and no square-freeness to test.
    ModularVerdict modular;
    // Every root zero (program_monomial): the root bound is then zero, and the answer is the point 0, its count and
    // its multiplicity the degree.
    bool monomial = false;
    // The order of a root at zero (program_zero_order), when not a monomial: the route solves q = p / x^v (its
    // ProgramSource deflated) and publishes the point 0, count and multiplicity v, beside q's regions.
    uint64_t zero_order = 0;
    Magnitude root_bound;
    bool ready() const { return refusal == Refusal::None; }
};
// Validates and prepares a program: an output, a length within exact_operations, the proved degree and leading
// coefficient (section 4; refusals: identically zero, the leading coefficient not determined, a constant over the
// exact storage budget), and for degree 1 or more the modular verdict and the root bound. A program proved a
// nonzero constant is ready with degree 0 (no roots), as a constant coefficient input is.
ProgramPreparation prepare_program(const Program &program, const cps_options &options);

// The resident storage of a solve on the structured route, in bytes, for the rule of section 3.1 (temporary_bytes
// bounds the storage; the operations are bounded by prepare_program): refused before anything of the degree's
// size is allocated. An estimate of what is held for the solve's life at the working precision, not a bound on
// peak memory (certpolysolve.h says so of every resource limit). The route prices it at the precision cap
// (max_work_bits, 65,536 by default), the most the solve may hold, although a solve of well-separated roots
// stays near the initial precision: with the default limits nothing is refused, and a caller who sets
// temporary_bytes sets it against the cap. Saturates at UINT64_MAX.
// History: docs/slp-input-spec-20261015.md, section 13.
uint64_t program_solve_storage(uint64_t degree, uint64_t precision_bits);

} // namespace certpolysolve::detail
