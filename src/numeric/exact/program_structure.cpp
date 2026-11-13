// Structure read from a program's graph (declared in program_structure.hpp), for the program route's lifts and
// structural seeds: the sub-program a node computes (cone), and the output's forms, a power c q^k (output_power), a
// product of distinct factors (output_factors), a squaring recurrence or composition F(x, q) (squaring_recurrence), and
// the outer function f of a composition f(q) through a cut node (above_cut).
//
// Everything here is exact and syntactic: the functions read the nodes and their degree bounds (Program::degree_bounds)
// and never evaluate the program. cone marks the nodes a node reads in one backward pass (operands are earlier nodes)
// and copies them renumbered, the constants whole. output_power descends from the output through Power, Multiply(q, q),
// Negate and products with one constant operand, multiplying the exponents. output_factors descends the same way and
// also into each operand of positive degree bound of any other product, adding the exponents of a factor reached twice.
// squaring_recurrence tries the nodes of the output's cone of degree bound at least 2 and at most half the output's,
// largest first and the later node among equal degrees, at most `candidates` of them: for each q it reads every node's
// degree bounds in x and in q, with q as a second variable (sums take the larger, products add, powers multiply), and
// returns the first whose output has degree d in q between 2 and max_degree_in_base, degree in x at most max_x_degree,
// and d times q's degree bound at most the output's. above_cut copies the nodes above the base with the base as the
// variable, provided every node they read that does not depend on the base is a constant (degree bound 0), which it
// copies with its cone. Degree arithmetic saturates; an exponent past UINT64_MAX gives no form.
//
// State: none.
// Entry points: cone, from solver_program.cpp (program_advance: a factor of the product lift, a power's base) and
// program_recurrence_seeds.cpp (the sub-programs of q and of g); output_factors and output_power, from
// solver_program.cpp (program_advance: the product lift and the power lift); squaring_recurrence and above_cut, from
// program_recurrence_seeds.cpp.
// Not here: solving the parts and lifting their roots (program_advance in solver_program.cpp), the seeds built on these
// forms (program_recurrence_seeds.cpp).
// Design notes: docs/slp-structural-seeding-note-20261018.md (section 2, the structures read from the graph;
// section 10, the products).
#include "numeric/exact/program_structure.hpp"
#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>

namespace certpolysolve::detail {

Program cone(const Program &program, uint64_t node) {
    if (node >= program.nodes.size())
        throw std::invalid_argument("cone of node " + std::to_string(node) + ", which is not in the program");
    // Operands are earlier nodes, so one backward pass marks everything `node` reads.
    std::vector<uint8_t> needed(node + 1, 0);
    needed[node] = 1;
    for (uint64_t i = node + 1; i-- > 0;) {
        if (!needed[i])
            continue;
        const auto &n = program.nodes[i];
        switch (n.op) {
        case ProgramOp::Variable:
        case ProgramOp::Constant: break;
        case ProgramOp::Negate:
        case ProgramOp::Power: needed[n.a] = 1; break;
        case ProgramOp::Add:
        case ProgramOp::Subtract:
        case ProgramOp::Multiply: needed[n.a] = needed[n.b] = 1; break;
        }
    }
    Program out;
    out.constants = program.constants;
    std::vector<uint64_t> renumbered(node + 1, 0);
    for (uint64_t i = 0; i <= node; ++i) {
        if (!needed[i])
            continue;
        ProgramNode copy = program.nodes[i];
        switch (copy.op) {
        case ProgramOp::Variable:
        case ProgramOp::Constant: break;
        case ProgramOp::Negate:
        case ProgramOp::Power: copy.a = renumbered[copy.a]; break;
        case ProgramOp::Add:
        case ProgramOp::Subtract:
        case ProgramOp::Multiply:
            copy.a = renumbered[copy.a];
            copy.b = renumbered[copy.b];
            break;
        }
        renumbered[i] = out.append(copy);
    }
    out.output = renumbered[node];
    return out;
}

std::optional<ProgramPower> output_power(const Program &program) {
    if (!program.output || *program.output >= program.nodes.size())
        return std::nullopt;
    const auto degrees = program.degree_bounds();
    ProgramPower power{*program.output, 1};
    for (;;) {
        const auto &n = program.nodes[power.base];
        uint64_t exponent = 1;
        if (n.op == ProgramOp::Power) {
            exponent = n.value;
        } else if (n.op == ProgramOp::Multiply && n.a == n.b) {
            exponent = 2;
        } else if (n.op == ProgramOp::Negate) {
            power.base = n.a;
            continue;
        } else if (n.op == ProgramOp::Multiply && (degrees[n.a] == 0) != (degrees[n.b] == 0)) {
            power.base = degrees[n.a] == 0 ? n.b : n.a; // a constant factor
            continue;
        } else {
            break;
        }
        if (power.exponent > UINT64_MAX / exponent)
            return std::nullopt;
        power.exponent *= exponent;
        power.base = n.a;
    }
    if (power.exponent < 2)
        return std::nullopt;
    return power;
}

std::optional<std::vector<ProgramPower>> output_factors(const Program &program) {
    if (!program.output || *program.output >= program.nodes.size())
        return std::nullopt;
    const auto degrees = program.degree_bounds();
    std::vector<uint64_t> exponent(program.nodes.size(), 0);
    std::vector<ProgramPower> stack{{*program.output, 1}};
    while (!stack.empty()) {
        const auto [node, k] = stack.back();
        stack.pop_back();
        const auto &n = program.nodes[node];
        const auto times = [&](uint64_t e) -> std::optional<uint64_t> {
            if (k > UINT64_MAX / e)
                return std::nullopt;
            return k * e;
        };
        if (n.op == ProgramOp::Power || (n.op == ProgramOp::Multiply && n.a == n.b)) {
            const auto e = times(n.op == ProgramOp::Power ? n.value : 2);
            if (!e)
                return std::nullopt;
            stack.push_back({n.a, *e});
        } else if (n.op == ProgramOp::Negate) {
            stack.push_back({n.a, k});
        } else if (n.op == ProgramOp::Multiply) {
            if (degrees[n.a] != 0)
                stack.push_back({n.a, k});
            if (degrees[n.b] != 0)
                stack.push_back({n.b, k});
        } else if (degrees[node] != 0) {
            if (exponent[node] > UINT64_MAX - k)
                return std::nullopt;
            exponent[node] += k;
        }
    }
    std::vector<ProgramPower> factors;
    for (uint64_t i = 0; i < exponent.size(); ++i)
        if (exponent[i])
            factors.push_back({i, exponent[i]});
    if (factors.size() < 2)
        return std::nullopt;
    return factors;
}

std::optional<ProgramRecurrence> squaring_recurrence(const Program &program, uint64_t max_x_degree,
                                                     uint64_t max_degree_in_base, size_t candidates) {
    if (!program.output || *program.output >= program.nodes.size())
        return std::nullopt;
    const uint64_t output = *program.output;
    const auto degrees = program.degree_bounds();
    const uint64_t top = degrees[output];
    if (top < 3 || top == UINT64_MAX)
        return std::nullopt;
    // The output's cone, and the candidates in it: degree at least 1 and at most half the output's.
    std::vector<uint8_t> in_cone(output + 1, 0);
    in_cone[output] = 1;
    for (uint64_t i = output + 1; i-- > 0;) {
        if (!in_cone[i])
            continue;
        const auto &n = program.nodes[i];
        if (n.op == ProgramOp::Negate || n.op == ProgramOp::Power)
            in_cone[n.a] = 1;
        else if (n.op == ProgramOp::Add || n.op == ProgramOp::Subtract || n.op == ProgramOp::Multiply)
            in_cone[n.a] = in_cone[n.b] = 1;
    }
    // A base of degree 1 (x itself, or x - 3) makes any program F(x, q) trivially, so the base has degree 2 or more.
    std::vector<uint64_t> order;
    for (uint64_t i = 0; i < output; ++i)
        if (in_cone[i] && degrees[i] >= 2 && degrees[i] <= top / 2)
            order.push_back(i);
    // Largest degree first; among equal degrees the later node, nearer the output (x (x^2 + 1)^3 - 2 is also
    // F(x, x^2), but x^2's root is double, a poor base for seeds).
    std::sort(order.begin(), order.end(),
              [&](uint64_t a, uint64_t b) { return degrees[a] != degrees[b] ? degrees[a] > degrees[b] : a > b; });
    if (order.size() > candidates)
        order.resize(candidates);
    // Saturating arithmetic on the degree pairs.
    const auto add = [](uint64_t a, uint64_t b) { return a > UINT64_MAX - b ? UINT64_MAX : a + b; };
    const auto times = [](uint64_t a, uint64_t e) { return a && e > UINT64_MAX / a ? UINT64_MAX : a * e; };
    std::vector<uint64_t> in_x(output + 1), in_q(output + 1);
    std::vector<uint8_t> above(output + 1);
    for (const uint64_t q : order) {
        std::fill(above.begin(), above.end(), uint8_t{0});
        for (uint64_t i = 0; i <= output; ++i) {
            if (!in_cone[i])
                continue;
            if (i == q) {
                in_x[i] = 0, in_q[i] = 1, above[i] = 1;
                continue;
            }
            const auto &n = program.nodes[i];
            const bool reads_q = (n.op != ProgramOp::Variable && n.op != ProgramOp::Constant && above[n.a]) ||
                                 ((n.op == ProgramOp::Add || n.op == ProgramOp::Subtract || n.op == ProgramOp::Multiply) &&
                                  above[n.b]);
            if (!reads_q) {
                in_x[i] = degrees[i], in_q[i] = 0;
                continue;
            }
            above[i] = 1;
            switch (n.op) {
            case ProgramOp::Variable:
            case ProgramOp::Constant: break;
            case ProgramOp::Negate: in_x[i] = in_x[n.a], in_q[i] = in_q[n.a]; break;
            case ProgramOp::Power:
                in_x[i] = times(in_x[n.a], n.value), in_q[i] = times(in_q[n.a], n.value);
                break;
            case ProgramOp::Add:
            case ProgramOp::Subtract:
                in_x[i] = std::max(in_x[n.a], in_x[n.b]), in_q[i] = std::max(in_q[n.a], in_q[n.b]);
                break;
            case ProgramOp::Multiply:
                in_x[i] = add(in_x[n.a], in_x[n.b]), in_q[i] = add(in_q[n.a], in_q[n.b]);
                break;
            }
        }
        if (!above[output] || in_q[output] < 2 || in_q[output] > max_degree_in_base || in_x[output] > max_x_degree)
            continue;
        const uint64_t covered = times(in_q[output], degrees[q]);
        if (covered > top)
            continue;
        return ProgramRecurrence{q, in_q[output], in_x[output], top - covered};
    }
    return std::nullopt;
}

std::optional<Program> above_cut(const Program &program, uint64_t base) {
    if (!program.output || *program.output >= program.nodes.size() || base >= *program.output)
        return std::nullopt;
    const uint64_t output = *program.output;
    const auto degrees = program.degree_bounds();
    std::vector<uint8_t> above(output + 1, 0);
    above[base] = 1;
    for (uint64_t i = base + 1; i <= output; ++i) {
        const auto &n = program.nodes[i];
        if (n.op == ProgramOp::Negate || n.op == ProgramOp::Power)
            above[i] = above[n.a];
        else if (n.op == ProgramOp::Add || n.op == ProgramOp::Subtract || n.op == ProgramOp::Multiply)
            above[i] = above[n.a] || above[n.b];
    }
    if (!above[output])
        return std::nullopt;
    Program out;
    out.constants = program.constants;
    std::vector<uint64_t> renumbered(output + 1, UINT64_MAX);
    ProgramNode variable;
    variable.op = ProgramOp::Variable;
    renumbered[base] = out.append(variable);
    // A node read from above that does not depend on the base: a constant, copied with its cone.
    bool constant = true;
    const auto copy = [&](auto &self, uint64_t j) -> uint64_t {
        if (renumbered[j] != UINT64_MAX)
            return renumbered[j];
        if (degrees[j] != 0) {
            constant = false;
            return 0;
        }
        ProgramNode n = program.nodes[j];
        switch (n.op) {
        case ProgramOp::Variable: constant = false; return 0;
        case ProgramOp::Constant: break;
        case ProgramOp::Negate:
        case ProgramOp::Power: n.a = self(self, n.a); break;
        case ProgramOp::Add:
        case ProgramOp::Subtract:
        case ProgramOp::Multiply:
            n.a = self(self, n.a);
            n.b = self(self, n.b);
            break;
        }
        if (!constant)
            return 0;
        return renumbered[j] = out.append(n);
    };
    for (uint64_t i = base + 1; i <= output && constant; ++i) {
        if (!above[i])
            continue;
        ProgramNode n = program.nodes[i];
        switch (n.op) {
        case ProgramOp::Variable:
        case ProgramOp::Constant: break;
        case ProgramOp::Negate:
        case ProgramOp::Power: n.a = copy(copy, n.a); break;
        case ProgramOp::Add:
        case ProgramOp::Subtract:
        case ProgramOp::Multiply:
            n.a = copy(copy, n.a);
            n.b = copy(copy, n.b);
            break;
        }
        if (!constant)
            break;
        renumbered[i] = out.append(n);
    }
    if (!constant)
        return std::nullopt;
    out.output = renumbered[output];
    return out;
}

} // namespace certpolysolve::detail
