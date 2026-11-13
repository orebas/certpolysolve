#pragma once
// Interpret a program with an arithmetic supplied by the caller. Release each value at its last use.
#include "input/program/program.hpp"
#include "numeric/evaluate/program_slots.hpp"
#include <algorithm>
#include <stdexcept>
#include <utility>

namespace certpolysolve::detail {
// The interpreter. Arithmetic provides a movable Value and variable(), constant(const Coefficient&) or
// constant_at(index) (an arithmetic that prepared the program's constants once takes the index),
// add(a, b), subtract(a, b), multiply(a, b), negate(a), power(a, exponent). Each value is released after
// its last reader, so a program as long as a three-term recurrence of length n keeps a few values live,
// not n. `peak_live`, when given, receives the most values held at once.
template <class Arithmetic>
typename Arithmetic::Value run_program(const Program &program, Arithmetic &arithmetic, size_t *peak_live = nullptr) {
    if (!program.output || *program.output >= program.nodes.size())
        throw std::invalid_argument("program has no output");
    const auto last = program.last_uses();
    std::vector<std::optional<typename Arithmetic::Value>> values(program.nodes.size());
    size_t live = 0, peak = 0;
    for (size_t i = 0; i < program.nodes.size(); ++i) {
        const auto &node = program.nodes[i];
        const auto &a = values[node.a];
        const auto &b = values[node.b];
        switch (node.op) {
        case ProgramOp::Variable: values[i].emplace(arithmetic.variable()); break;
        case ProgramOp::Constant:
            if constexpr (requires { arithmetic.constant_at(node.value); })
                values[i].emplace(arithmetic.constant_at(node.value));
            else
                values[i].emplace(arithmetic.constant(program.constants[node.value]));
            break;
        case ProgramOp::Add: values[i].emplace(arithmetic.add(*a, *b)); break;
        case ProgramOp::Subtract: values[i].emplace(arithmetic.subtract(*a, *b)); break;
        case ProgramOp::Multiply: values[i].emplace(arithmetic.multiply(*a, *b)); break;
        case ProgramOp::Negate: values[i].emplace(arithmetic.negate(*a)); break;
        case ProgramOp::Power: values[i].emplace(arithmetic.power(*a, node.value)); break;
        }
        peak = std::max(peak, ++live);
        // Release the operands this node read for the last time, and a value nothing reads.
        const auto release = [&](uint64_t j) {
            if (values[j] && last[j] == i && j != *program.output) {
                values[j].reset();
                --live;
            }
        };
        if (node.op != ProgramOp::Variable && node.op != ProgramOp::Constant) {
            release(node.a);
            if (node.op != ProgramOp::Negate && node.op != ProgramOp::Power && node.b != node.a)
                release(node.b);
        }
        if (last[i] == program.nodes.size() && i != *program.output) {
            values[i].reset();
            --live;
        }
    }
    if (peak_live)
        *peak_live = peak;
    return std::move(*values[*program.output]);
}

// Borrows an immutable program that must outlive it. Prepare last uses and reusable
// slots once; run_prepared keeps only the peak live value storage between evaluations.
struct PreparedProgram {
    const Program *program;
    ProgramSlots slots;
    explicit PreparedProgram(const Program &p) : program(&p), slots(p) {}
    PreparedProgram(Program &&) = delete;
    PreparedProgram(const Program &&) = delete;
};
template <class Arithmetic>
typename Arithmetic::Value run_prepared(const PreparedProgram &prepared, Arithmetic &arithmetic,
                                        std::vector<std::optional<typename Arithmetic::Value>> &values) {
    const auto &program = *prepared.program;
    const auto &last = prepared.slots.last;
    const auto &slot = prepared.slots.slot;
    values.resize(prepared.slots.count);
    // A failed arithmetic call must release its live operands too. Keep capacity
    // for reuse, but never retain a half-evaluated program after an exception.
    struct Cleanup {
        std::vector<std::optional<typename Arithmetic::Value>> &values;
        bool needed = true;
        ~Cleanup() { if (needed) for (auto &value : values) value.reset(); }
    } cleanup{values};
    const uint64_t output = *program.output;
    for (size_t i = 0; i < program.nodes.size(); ++i) {
        const auto &node = program.nodes[i];
        auto &target = values[slot[i]];
        switch (node.op) {
        case ProgramOp::Variable: target.emplace(arithmetic.variable()); break;
        case ProgramOp::Constant:
            if constexpr (requires { arithmetic.constant_at(node.value); })
                target.emplace(arithmetic.constant_at(node.value));
            else
                target.emplace(arithmetic.constant(program.constants[node.value]));
            break;
        case ProgramOp::Add: target.emplace(arithmetic.add(*values[slot[node.a]], *values[slot[node.b]])); break;
        case ProgramOp::Subtract: target.emplace(arithmetic.subtract(*values[slot[node.a]], *values[slot[node.b]])); break;
        case ProgramOp::Multiply: target.emplace(arithmetic.multiply(*values[slot[node.a]], *values[slot[node.b]])); break;
        case ProgramOp::Negate: target.emplace(arithmetic.negate(*values[slot[node.a]])); break;
        case ProgramOp::Power: target.emplace(arithmetic.power(*values[slot[node.a]], node.value)); break;
        }
        if (node.op != ProgramOp::Variable && node.op != ProgramOp::Constant) {
            if (last[node.a] == i && node.a != output) values[slot[node.a]].reset();
            if (node.op != ProgramOp::Negate && node.op != ProgramOp::Power && node.b != node.a && last[node.b] == i &&
                node.b != output)
                values[slot[node.b]].reset();
        }
        if (last[i] == program.nodes.size() && i != output) target.reset();
    }
    auto result = std::move(*values[slot[output]]);
    values[slot[output]].reset();
    cleanup.needed = false;
    return result;
}

} // namespace certpolysolve::detail
