#pragma once
// Owned input, or the C API's immutable term borrow for one synchronous call.
// Numerical operations and route selection live outside this storage object.
#include "numeric/arithmetic/polynomial_arithmetic.hpp"
#include "input/program/program.hpp"
#include "input/secular/secular.hpp"
#include "solve/real_query.hpp"
#include <algorithm>
#include <memory>
#include <stdexcept>

namespace certpolysolve::detail {
struct QueryInput {
    std::unique_ptr<DensePolynomial> input;
    // Or the input object's terms, borrowed for the one advance (cps_solve's call; C40's second step): the questions
    // asked before a route is built are answered from them, and `input` is formed only for a route that reads every
    // coefficient (dense_input). The real-only route takes the borrow on.
    const cps_polynomial *terms = nullptr;
    bool polynomial() const { return input || terms; }
    DensePolynomial &dense_input() {
        if (!input) input = std::make_unique<DensePolynomial>(dense(*terms));
        return *input;
    }
    std::unique_ptr<cps_secular> secular;
    // Program input (docs/slp-input-spec-20261015.md): the structured route on the program, in Aberth mode
    // under the classical recipe's policy (Auto and classical) or in secular mode under the secular recipe's.
    std::unique_ptr<cps_program> program;
    // The input is a program: kept for the session's life (the program itself is handed to the solve state at
    // the first call), since every call resolves the default goal for the input kind (with_resolved_goal).
    bool program_kind = false;
    explicit QueryInput(const DensePolynomial &p) : input(std::make_unique<DensePolynomial>(p)) {}
    explicit QueryInput(DensePolynomial &&p) : input(std::make_unique<DensePolynomial>(std::move(p))) {}
    explicit QueryInput(const cps_polynomial &borrowed) : terms(&borrowed) {}
    explicit QueryInput(const cps_secular &p) : secular(std::make_unique<cps_secular>(p)) {}
    explicit QueryInput(const cps_program &p) : program(std::make_unique<cps_program>(p)), program_kind(true) {
        if (!program->program.output || *program->program.output >= program->program.nodes.size())
            throw std::invalid_argument("program has no output");
    }
    // Engine orchestration borrows this stable storage object. Moving the owned
    // polynomial into an engine is explicit; moving the storage itself is forbidden.
    QueryInput(const QueryInput &) = delete;
    QueryInput &operator=(const QueryInput &) = delete;
    QueryInput(QueryInput &&) = delete;
    QueryInput &operator=(QueryInput &&) = delete;

    // The input's degree, for the work budget and an unsupported query's result: the coefficients' length, the
    // secular terms' count, or the program's degree bound (its exact degree is proved at the first advance, with
    // the options' budgets; the budget is a soft limit, and the bound is at least the degree).
    uint64_t input_degree() const {
        if (program) {
            const auto &p = program->program;
            return p.degree_bounds()[*p.output];
        }
        const size_t length = input_length();
        return length ? static_cast<uint64_t>(length - 1) : 0;
    }
    // Coefficient or secular input: the coefficients' length without trailing zeros, or the secular terms' count
    // and one (0: an input with no nonzero coefficient).
    size_t input_length() const {
        if (!input && terms) return terms->terms.empty() ? 0 : static_cast<size_t>(terms->terms.rbegin()->first) + 1;
        size_t length = input ? input->coefficients.size() : secular->terms.size() + 1;
        if (input) while (length && input->coefficients[length - 1].zero()) --length;
        return length;
    }
    // real_query::applies on whichever form the input has.
    bool real_route_applies(const cps_options &o) const {
        return input ? real_query::applies(o, *input) : real_query::applies(o, *terms);
    }
    bool exact_coefficients() const {
        if (input)
            return std::all_of(input->coefficients.begin(), input->coefficients.end(),
                               [](const Coefficient &c) { return c.exact(); });
        return terms && std::all_of(terms->terms.begin(), terms->terms.end(),
                                    [](const auto &t) { return t.second.exact(); });
    }
    bool supported(const cps_options &o) const {
        // Complex isolation uses the standalone exact-coefficient path. A region
        // still belongs to a real-only query; secular/program isolation is unavailable.
        if (o.goal == CPS_ROOT_ISOLATION && !o.real_only)
            return polynomial() && o.method == CPS_METHOD_AUTO && !o.region.set() && exact_coefficients();
        if (o.goal == CPS_ROOT_ISOLATION || o.region.set())
            return polynomial() && real_route_applies(o) &&
                   o.goal != CPS_CLUSTER_RADIUS;
        // A program: all complex roots, every other goal, by Auto (Aberth mode), classical (the same) or secular
        // (secular mode); never real-only, never piecewise (docs/slp-input-spec-20261015.md, section 2).
        if (program) return !o.real_only && o.method != CPS_METHOD_PIECEWISE;
        if (o.method == CPS_METHOD_AUTO) return true;
        if (o.real_only || (o.goal != CPS_ROOT_ACCURACY && o.goal != CPS_ROOT_APPROXIMATION)) return false;
        if (secular) return o.method == CPS_METHOD_SECULAR;
        return exact_coefficients();
    }
};
} // namespace certpolysolve::detail
