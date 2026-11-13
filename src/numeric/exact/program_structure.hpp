#pragma once
// Structure read from a program's graph, for the structural seeds and lifts of
// docs/slp-structural-seeding-note-20261018.md: the sub-programs a structural rule solves first, and the rules'
// detection. Everything here is exact and syntactic: it reads the nodes the caller wrote and never evaluates them.
#include "input/program/program.hpp"
#include <cstdint>
#include <optional>

namespace certpolysolve::detail {

// The program computing `node`'s value: the nodes `node` depends on, in their original order and renumbered, with
// `node` as the output. The constants are kept whole, so the cone's Constant nodes index them as before. Throws
// std::invalid_argument when `node` is not a node of the program.
Program cone(const Program &program, uint64_t node);

// The output as c q^k for an earlier node q, a nonzero constant c and k >= 2: descending from the output through
// Power(q, e) and Multiply(q, q) (the exponents multiplying), Negate, and products with a constant operand (a node
// of degree bound 0, which a program whose degree is proved positive cannot have zero). The roots of the output are
// then those of q, each k times as often. Empty when the output is no such power, or when the exponent would pass
// UINT64_MAX.
struct ProgramPower {
    uint64_t base = 0, exponent = 1;
};
std::optional<ProgramPower> output_power(const Program &program);

// The output as c q_1^k_1 ... q_m^k_m with m >= 2 distinct factors of positive degree bound: descending from the
// output through Multiply (both operands of positive degree bound are factors to descend into), Power(q, e) and
// Multiply(q, q) (the exponents multiplying), Negate, and products with a constant operand (degree bound 0). The
// factors are the nodes reached that are none of those, with their exponents (a node reached twice adds them), in
// increasing node order. The roots of the output are the factors' roots, each k_i times as often. Empty when fewer
// than two distinct factors of positive degree bound remain (one is output_power's case), or when an exponent
// would pass UINT64_MAX.
std::optional<std::vector<ProgramPower>> output_factors(const Program &program);

// The output as F(x, q) for an earlier node q of degree 2 or more, where F has degree d in q, 2 <= d <=
// `max_degree_in_base`, and its coefficients have degree at most `max_x_degree` in x: a squaring recurrence such as
// Mandelbrot's p_(k+1) = x p_k^2 + 1 (q = p_k, d = 2, x degree 1). Read by degree bounds with q as a second
// variable: a node that does not depend on q keeps its degree bound in x; q is (0, 1); sums take the larger of each,
// products add them, powers multiply them. Any node of high degree in x that reaches the output around q makes the
// x degree too large, so such a q is no base. Among the nodes of degree bound at most deg p / 2, the `candidates` of
// largest degree bound are tried, largest first (the later node first among equal degrees), and the first that
// passes is returned. `leftover` is the output's degree bound minus d times q's, by the bounds: the proved degrees
// can be lower (the recorder's x 0 + 1 has bound 1).
struct ProgramRecurrence {
    uint64_t base = 0, degree_in_base = 0, degree_in_x = 0, leftover = 0;
};
std::optional<ProgramRecurrence> squaring_recurrence(const Program &program, uint64_t max_x_degree = 4,
                                                     uint64_t max_degree_in_base = 8, size_t candidates = 64);

// The program above a cut: when every path from the variable to the output passes through `base` (the output is
// f(q), squaring_recurrence with no x outside q), the program computing f, with q as its variable. Nodes above the
// base are copied in order, the base becomes the variable, and the nodes they read that do not depend on the base
// are copied with their cones, which must be constants (degree bound 0). Empty otherwise.
std::optional<Program> above_cut(const Program &program, uint64_t base);

} // namespace certpolysolve::detail
