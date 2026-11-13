#pragma once
// Exact straight-line input: owned nodes and constants, checked construction and text parsing.
// Numerical interpretation and arithmetic belong to numeric/evaluate.
#include "core/number.hpp"
#include <cstdint>
#include <optional>
#include <string_view>
#include <vector>

namespace certpolysolve::detail {

enum class ProgramOp : uint8_t { Variable, Constant, Add, Subtract, Multiply, Negate, Power };

struct ProgramNode {
    ProgramOp op = ProgramOp::Variable;
    uint64_t a = 0, b = 0; // operands: earlier nodes
    uint64_t value = 0;    // Constant: index into `constants`; Power: the exponent (at least 2)
};

struct Program {
    std::vector<ProgramNode> nodes;
    std::vector<Coefficient> constants; // exact: zero radii
    std::optional<uint64_t> output;

    // Appends a node after checking its operands and exponent; returns its index.
    uint64_t append(const ProgramNode &node);
    uint64_t append_constant(Coefficient constant);
    // For each node, the index of the last node reading it, or nodes.size() for the output and a node
    // nothing reads after it is computed (freed at once). A value can be released after its last reader.
    std::vector<uint64_t> last_uses() const;
    // An upper bound for each node's degree, saturating at UINT64_MAX.
    std::vector<uint64_t> degree_bounds() const;
};

// The text form (PSPROGRAM 1; one operation a line: x, c REAL [IMAG], add A B, sub A B, mul A B, neg A,
// pow A E, out N; '#' starts a comment). Throws std::invalid_argument naming the line.
Program parse_program(std::string_view text);

} // namespace certpolysolve::detail

struct cps_program {
    certpolysolve::detail::Program program;
};
