// C40's second step: a solve that borrows the input object's terms (SessionState(const cps_polynomial &), what
// cps_solve does) answers exactly as one given the expanded coefficient list (SessionState(DensePolynomial)): the same
// regions, counts, multiplicities, realness and certificates, the same diagnostics, progress, components and
// counters; only timings may differ. Real-only queries on few terms (the sparse recursion, deciding its attempt from
// the terms), the same with a region and under the isolation goal, a budget too small for the attempts, a dense real
// input, all-complex queries, an inexact input and a zero polynomial. A miniature of the identity run that is this
// change's gate.
#include "internal.hpp"
#include "numeric/arithmetic/polynomial_arithmetic.hpp"
#include <certpolysolve/certpolysolve.h>
#include <flint/fmpz.h>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace d = certpolysolve::detail;

namespace {
void require(bool condition, const std::string &message) {
    if (!condition) throw std::runtime_error(message);
}
// A polynomial by its terms (exponent, real text, imaginary text, real radius).
struct Term {
    uint64_t exponent;
    const char *real, *imag = "0", *radius = "0";
};
cps_polynomial polynomial(const std::vector<Term> &terms) {
    cps_polynomial *p = nullptr;
    require(cps_polynomial_create(&p) == CPS_OK, "create");
    for (const auto &t : terms)
        require(cps_polynomial_set_ball_coefficient(p, t.exponent, t.real, t.imag, t.radius, "0") == CPS_OK, "set");
    cps_polynomial copy = *p;
    cps_polynomial_destroy(p);
    return copy;
}
cps_polynomial from_roots(int count) {
    std::vector<std::string> texts;
    for (int k = 1; k <= count; ++k) texts.push_back(std::to_string(k));
    std::vector<const char *> roots;
    for (const auto &t : texts) roots.push_back(t.c_str());
    cps_polynomial *p = nullptr;
    require(cps_polynomial_from_roots(roots.data(), nullptr, nullptr, roots.size(), &p) == CPS_OK, "from_roots");
    cps_polynomial copy = *p;
    cps_polynomial_destroy(p);
    return copy;
}
// (3x - 1)(5x - 7)(7x + 17) prod_{i=1..126} (x^2 + 9^i), test_real_interval's input for the axis charts: the selector
// sends its complex query to piecewise, so the real-only route attempts the axis charts and forms the list for them.
cps_polynomial axis_input() {
    std::vector<std::string> real{"1/3", "7/5", "-17/7"}, imag{"0", "0", "0"};
    fmpz_t power;
    fmpz_init(power);
    for (ulong i = 1; i <= 126; ++i) {
        fmpz_ui_pow_ui(power, 3, i);
        char *text = fmpz_get_str(nullptr, 10, power);
        real.push_back("0"), imag.push_back(text);
        real.push_back("0"), imag.push_back(std::string("-") + text);
        flint_free(text);
    }
    fmpz_clear(power);
    std::vector<const char *> re, im;
    for (size_t k = 0; k < real.size(); ++k) re.push_back(real[k].c_str()), im.push_back(imag[k].c_str());
    cps_polynomial *p = nullptr;
    require(cps_polynomial_from_roots(re.data(), im.data(), nullptr, re.size(), &p) == CPS_OK, "from_roots");
    cps_polynomial copy = *p;
    cps_polynomial_destroy(p);
    return copy;
}
// Every field of the result but its timings.
std::string describe(const cps_result &r) {
    std::string s = "termination " + std::to_string(r.termination) + " goal " + std::to_string(r.goal) + " met " +
                    std::to_string(r.goal_met) + " complete " + std::to_string(r.complete) + " degree " +
                    std::to_string(r.degree) + " bits " + std::to_string(r.work_bits) + "\n";
    for (const auto &root : r.roots) {
        for (const auto &e : root.endpoints) s += e.mantissa_text() + "@" + e.exponent_text() + " ";
        s += "count " + std::to_string(root.count) + " multiplicity " + std::to_string(root.multiplicity) +
             " realness " + std::to_string(root.realness) + " " + root.certificate + "\n";
    }
    const auto &p = r.progress;
    s += "progress " + std::to_string(p.work_units_this_call) + " " + std::to_string(p.pending_components) + " " +
         std::to_string(p.live_cells) + " " + std::to_string(p.min_pending_bits) + " " +
         std::to_string(p.max_pending_bits) + " " + std::to_string(p.max_outer_bits) + " " +
         std::to_string(p.preprocessing_limited) + " " + std::to_string(p.input_uncertainty) + " " +
         std::to_string(p.stop_reason) + "\n";
    const auto &v = r.statistics_v3;
    s += "counters " + std::to_string(v.candidate_sweeps) + " " + std::to_string(v.evaluations) + " " +
         std::to_string(v.evaluation_terms) + " " + std::to_string(v.root_interactions) + " " +
         std::to_string(v.validation_calls) + " " + std::to_string(v.secular_regenerations) + " " +
         std::to_string(v.double_seed_calls) + " " + std::to_string(v.candidate_stagnations) + " " +
         std::to_string(v.real_refinement_steps) + " " + std::to_string(v.real_bisections) + "\n";
    for (const auto &c : r.components)
        s += "component " + std::to_string(c.factor_index) + " " + std::to_string(c.work_bits) + " " +
             std::to_string(c.live_cells) + " " + std::to_string(c.certified_count) + " " + std::to_string(c.state) + "\n";
    return s + r.diagnostics;
}
// The two sessions on one query must describe the same result, answered by the route named in `by` if one is.
void same(const std::string &name, const cps_polynomial &input, const cps_options &options, const char *by = nullptr) {
    d::SessionState dense_session(d::dense(input));
    const cps_result from_list = dense_session.advance(&options);
    d::SessionState borrowing(input);
    const cps_result from_terms = borrowing.advance(&options);
    const auto a = describe(from_list), b = describe(from_terms);
    if (a != b) {
        std::cerr << name << ": the coefficient list gives\n" << a << "\nthe terms give\n" << b << '\n';
        throw std::runtime_error(name + ": the borrowed terms answer differently");
    }
    if (by && from_terms.diagnostics.find(by) == std::string::npos)
        throw std::runtime_error(name + ": not answered by " + by + ": " + from_terms.diagnostics);
}
cps_options real_query() {
    cps_options o;
    o.real_only = true;
    return o;
}
} // namespace

int main() {
    try {
        const auto trinomial = polynomial({{100000, "1"}, {7, "-3"}, {0, "1"}});
        same("real roots of x^100000 - 3x^7 + 1", trinomial, real_query(), "by the sparse recursion");
        auto isolation = real_query();
        isolation.goal = CPS_ROOT_ISOLATION;
        same("isolation on x^100000 - 3x^7 + 1", trinomial, isolation, "by the sparse recursion");
        auto region = real_query();
        region.region.real_lower = d::Number("0");
        region.region.real_upper = d::Number("2");
        const auto small = polynomial({{300, "1"}, {7, "-3"}, {0, "1"}});
        same("x^300 - 3x^7 + 1 in [0, 2]", small, region);
        auto tight = real_query();
        tight.work_units = 10;
        tight.work_units_default = false;
        same("a budget too small for the attempts", small, tight, "by the complex route");
        same("real roots of (x - 1)...(x - 20)", from_roots(20), real_query());
        const auto mixed = polynomial({{60, "2/3"}, {41, "-5"}, {40, "1@-30"}, {3, "7"}, {0, "-1"}});
        same("real roots of a five-term polynomial", mixed, real_query());
        same("all roots of x^300 - 3x^7 + 1", small, cps_options{});
        same("all roots of a five-term polynomial", mixed, cps_options{});
        same("all roots of (x - 1)...(x - 20)", from_roots(20), cps_options{});
        const auto inexact = polynomial({{2, "1"}, {1, "-3", "0", "1e-9"}, {0, "2", "0", "1e-9"}});
        cps_options clusters;
        clusters.goal = CPS_CLUSTER_RADIUS;
        clusters.cluster_radius = d::Number("1e-6");
        same("an inexact quadratic", inexact, clusters);
        const auto axis = axis_input();
        same("real roots by the axis charts", axis, real_query(), "by the axis charts");
        same("isolation by the axis charts", axis, isolation, "by the axis charts");
        same("the zero polynomial", polynomial({}), cps_options{});
        same("the zero polynomial, real roots", polynomial({}), real_query());
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
    std::cout << "the borrowed terms answer as the coefficient list does\n";
    return 0;
}
