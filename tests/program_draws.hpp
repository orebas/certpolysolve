#pragma once
// The programs the program tests draw, shared by tests/test_program.cpp and tests/test_aberth_step.cpp (which runs the
// program kernels' correction against its shared function on machine_bits' draws): a node and a Gaussian constant
// appended by hand, a random program of a bounded degree, and the spec's Mandelbrot iteration with the internal
// program it records. Moved from test_program.cpp unchanged but for `inline` (W8's S10, 23 October 2026).
#include "input/program/program.hpp"
#include <cstdint>
#include <certpolysolve/certpolysolve.hpp>
#include <random>

namespace certpolysolve::detail {

inline Coefficient gaussian(const char *real, const char *imag = "0") {
    Coefficient c;
    c.real = Number(real);
    c.imag = Number(imag);
    return c;
}
inline uint64_t node(Program &p, ProgramOp op, uint64_t a = 0, uint64_t b = 0, uint64_t value = 0) {
    ProgramNode n;
    n.op = op;
    n.a = a;
    n.b = b;
    n.value = value;
    return p.append(n);
}

// A random program of `length` operations whose degree stays at most `max_degree`.
inline Program random_program(std::mt19937_64 &rng, size_t length, uint64_t max_degree) {
    static const char *reals[] = {"1", "-1", "3/2", "0", "2", "5@-3", "-7/3", "1@4"};
    static const char *imags[] = {"0", "0", "0", "1", "-1/2", "2"};
    Program p;
    node(p, ProgramOp::Variable);
    while (p.nodes.size() < length) {
        const auto degrees = p.degree_bounds();
        // Half the operands from the last three nodes, so that products build on each other.
        const auto pick = [&] {
            const uint64_t size = p.nodes.size();
            if (size > 3 && rng() % 2) return size - 1 - rng() % 3;
            return std::uniform_int_distribution<uint64_t>(0, size - 1)(rng);
        };
        switch (std::uniform_int_distribution<int>(0, 6)(rng)) {
        case 0: node(p, ProgramOp::Variable); break;
        case 1: p.append_constant(gaussian(reals[rng() % 8], imags[rng() % 6])); break;
        case 2: node(p, ProgramOp::Add, pick(), pick()); break;
        case 3: node(p, ProgramOp::Subtract, pick(), pick()); break;
        case 4: {
            const auto a = pick(), b = pick();
            if (degrees[a] + degrees[b] <= max_degree) node(p, ProgramOp::Multiply, a, b);
            break;
        }
        case 5: node(p, ProgramOp::Negate, pick()); break;
        case 6: {
            const auto a = pick();
            const uint64_t e = 2 + rng() % 2;
            if (degrees[a] * e <= max_degree) node(p, ProgramOp::Power, a, 0, e);
            break;
        }
        }
    }
    // The output: the node of the highest degree bound (the latest of equals).
    const auto degrees = p.degree_bounds();
    uint64_t output = 0;
    for (uint64_t i = 0; i < degrees.size(); ++i)
        if (degrees[i] >= degrees[output]) output = i;
    p.output = output;
    return p;
}

// The spec's Mandelbrot iteration: p_0 = 1, p_(k+1) = x p_k^2 + 1 (degree 2^k - 1; umand31 is k = 5).
inline certpolysolve::Program mandelbrot(int k) {
    return certpolysolve::Program::from([k](auto x) {
        auto p = x * 0 + 1;
        for (int j = 0; j < k; ++j) p = x * p * p + 1;
        return p;
    });
}
inline const Program &detail_of(const certpolysolve::Program &p) { return p.get()->program; }

} // namespace certpolysolve::detail
