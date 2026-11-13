#pragma once
// Reusable value slots for a program in node order. Metadata only: no borrowed owner.
#include "input/program/program.hpp"
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace certpolysolve::detail {

struct ProgramSlots {
    std::vector<uint64_t> last;
    std::vector<size_t> slot;
    size_t count = 0;

    explicit ProgramSlots(const Program &program)
        : last(program.last_uses()), slot(program.nodes.size()) {
        if (!program.output || *program.output >= program.nodes.size())
            throw std::invalid_argument("program has no output");
        std::vector<size_t> free;
        for (size_t i = 0; i < program.nodes.size(); ++i) {
            // The output takes its own slot before operands are released.
            // An arithmetic operation may therefore read both operands while writing it.
            if (free.empty()) slot[i] = count++;
            else {
                slot[i] = free.back();
                free.pop_back();
            }
            const auto &node = program.nodes[i];
            if (node.op != ProgramOp::Variable && node.op != ProgramOp::Constant) {
                if (last[node.a] == i && node.a != *program.output)
                    free.push_back(slot[node.a]);
                if (node.op != ProgramOp::Negate && node.op != ProgramOp::Power && node.b != node.a &&
                    last[node.b] == i && node.b != *program.output)
                    free.push_back(slot[node.b]);
            }
            // An unused node still executes, but its value need not survive it.
            if (last[i] == program.nodes.size() && i != *program.output)
                free.push_back(slot[i]);
        }
    }
};

} // namespace certpolysolve::detail
