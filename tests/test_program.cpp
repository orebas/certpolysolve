// Program input (docs/slp-input-spec-20261015.md), milestone 1: the representation, the text form, the
// interpreter's liveness and every arithmetic against the exact expansion; milestone 2: the preparation (refusals,
// the modular verdict, the root bound).
#include "certify/winding_count.hpp"
#include "numeric/evaluate/program_arithmetic.hpp"
#include "numeric/evaluate/program_interpreter.hpp"
#include "internal.hpp"
#include "api_test_support.hpp"
#include "input/program/program.hpp"
#include "program_reference_arithmetics.hpp"
#include "numeric/exact/program_preparation.hpp"
#include "certify/attempts/separated_roots.hpp"
#include "program_draws.hpp"
#include "numeric/evaluate/program_machine.hpp"
#include "numeric/candidates/program_candidates.hpp"
#include "numeric/candidates/program_seeds.hpp"
#include "numeric/evaluate/program_source.hpp"
#include "numeric/exact/program_structure.hpp"
#include <algorithm>
#include <cctype>
#include <flint/acb_poly.h>
#include <flint/arb_poly.h>
#include <flint/nmod.h>
#include <flint/ulong_extras.h>
#include <chrono>
#include <cmath>
#include <numbers>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <memory>
#include <limits>
#include <optional>
#include <certpolysolve/certpolysolve.hpp>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>

using namespace certpolysolve::detail;
namespace {
void require(bool value, const std::string &message) {
    if (!value) throw std::runtime_error(message);
}
constexpr uint64_t unlimited = UINT64_MAX;
// gaussian, node, random_program, mandelbrot and detail_of: tests/program_draws.hpp.

// The expansion's coefficient j modulo the prime, with i mapped.
ulong reduce_coefficient(const GaussianExpansion &e, slong j, const GaussianPrime &prime) {
    const auto residue = [&](const fmpq_poly_struct *part) -> ulong {
        Rational c;
        fmpq_poly_get_coeff_fmpq(c.value(), part, j);
        const ulong numerator = fmpz_fdiv_ui(fmpq_numref(c.value()), prime.prime);
        const ulong denominator = fmpz_fdiv_ui(fmpq_denref(c.value()), prime.prime);
        require(denominator != 0, "test constants have denominators coprime to the prime");
        nmod_t mod;
        nmod_init(&mod, prime.prime);
        return nmod_mul(numerator, nmod_inv(denominator, mod), mod);
    };
    nmod_t mod;
    nmod_init(&mod, prime.prime);
    return nmod_add(residue(e.real()), nmod_mul(residue(e.imag()), prime.imaginary_unit, mod), mod);
}

void check_against_expansion(const Program &p, const std::string &label) {
    ExpansionArithmetic exact(unlimited);
    const auto e = run_program(p, exact);
    const slong degree = e.degree();

    // Degree and leading coefficient from the top coefficients.
    const auto top = top_of_program(p, 64, unlimited);
    if (degree < 0) {
        require(top.verdict == ProgramTop::Verdict::ZeroPolynomial, label + ": zero polynomial not recognised");
    } else {
        require(top.verdict == ProgramTop::Verdict::Determined, label + ": degree undetermined");
        require(top.degree == static_cast<uint64_t>(degree), label + ": degree " + std::to_string(top.degree) +
                                                                 " against " + std::to_string(degree));
        Rational re, im;
        fmpq_poly_get_coeff_fmpq(re.value(), e.real(), degree);
        fmpq_poly_get_coeff_fmpq(im.value(), e.imag(), degree);
        require(fmpq_equal(re.value(), top.leading.real.value()) && fmpq_equal(im.value(), top.leading.imag.value()),
                label + ": leading coefficient differs");
    }

    // The modular image equals the expansion reduced.
    const auto prime = gaussian_prime(0);
    ModularArithmetic modular(prime);
    const auto image = run_program(p, modular);
    require(!modular.failed(), label + ": modular image failed");
    for (slong j = 0; j <= std::max<slong>(degree, nmod_poly_degree(image.value())); ++j) {
        const ulong expected = degree >= 0 && j <= degree ? reduce_coefficient(e, j, prime) : 0;
        require(nmod_poly_get_coeff_ui(image.value(), j) == expected, label + ": modular image differs at " +
                                                                           std::to_string(j));
    }

    // Value and derivative in balls contain the expansion's.
    constexpr slong precision = 256;
    ComplexPolynomial expansion;
    {
        arb_poly_t re, im;
        arb_poly_init(re);
        arb_poly_init(im);
        arb_poly_set_fmpq_poly(re, e.real(), precision);
        arb_poly_set_fmpq_poly(im, e.imag(), precision);
        acb_poly_set2_arb_poly(expansion.value(), re, im);
        arb_poly_clear(re);
        arb_poly_clear(im);
    }
    const char *points[][2] = {{"3/8", "-5/16"}, {"-2", "1/4"}, {"1/1024", "0"}, {"7/4", "7/4"}};
    for (const auto &z : points) {
        ComplexBall point, value, derivative;
        gaussian(z[0], z[1]).enclose(point.value(), precision);
        DualBallArithmetic dual(point.value(), precision);
        const auto v = run_program(p, dual);
        acb_poly_evaluate2(value.value(), derivative.value(), expansion.value(), point.value(), precision);
        require(acb_overlaps(v.value.value(), value.value()), label + ": value disagrees");
        require(acb_overlaps(v.derivative.value(), derivative.value()), label + ": derivative disagrees");
    }

    // The majorant bounds sum |a_i| r^i.
    Magnitude r;
    mag_set_ui_2exp_si(r.value(), 3, -1);
    MajorantArithmetic majorant(r.value());
    const auto bound = run_program(p, majorant);
    RealBall sum, term, power, rr;
    arb_set_ui(rr.value(), 3);
    arb_mul_2exp_si(rr.value(), rr.value(), -1);
    arb_one(power.value());
    for (slong j = 0; j <= degree; ++j) {
        ComplexBall c;
        acb_poly_get_coeff_acb(c.value(), expansion.value(), j);
        acb_abs(term.value(), c.value(), precision);
        arb_mul(term.value(), term.value(), power.value(), precision);
        arb_add(sum.value(), sum.value(), term.value(), precision);
        arb_mul(power.value(), power.value(), rr.value(), precision);
    }
    RealBall limit;
    arf_set_mag(arb_midref(limit.value()), bound.value());
    require(!arb_gt(sum.value(), limit.value()), label + ": majorant below the sum of moduli");
}

void random_programs() {
    std::mt19937_64 rng(20261015);
    int zero = 0, high = 0, complex_valued = 0;
    for (int t = 0; t < 300; ++t) {
        const auto p = random_program(rng, 10 + t % 30, 40);
        check_against_expansion(p, "random program " + std::to_string(t));
        ExpansionArithmetic exact(unlimited);
        const auto e = run_program(p, exact);
        zero += e.degree() < 0;
        high += e.degree() >= 10;
        complex_valued += fmpq_poly_degree(e.imag()) >= 0;
    }
    // The draw must reach the cases that matter, not only constants.
    require(high >= 50 && complex_valued >= 50, "random programs too simple: " + std::to_string(high) +
                                                    " of degree >= 10, " + std::to_string(complex_valued) +
                                                    " with complex coefficients");
    std::cout << "random programs: 300, " << high << " of degree >= 10, " << complex_valued
              << " with complex coefficients, " << zero << " identically zero\n";
}

// The structure read from a program's graph (program_structure.hpp): the cone of a node computes exactly that
// node's polynomial, with the node's degree bound, and the output is read as c q^k only where it is one.
void structure() {
    std::mt19937_64 rng(20261018);
    size_t cones = 0;
    for (int t = 0; t < 100; ++t) {
        const auto p = random_program(rng, 10 + t % 30, 40);
        const auto bounds = p.degree_bounds();
        for (uint64_t node = 0; node < p.nodes.size(); node += 1 + p.nodes.size() / 8) {
            Program same = p;
            same.output = node;
            const Program c = cone(p, node);
            require(c.output && c.nodes.size() <= node + 1, "cone: at most the node's prefix");
            ExpansionArithmetic exact(unlimited), exact_cone(unlimited);
            auto whole = run_program(same, exact);
            auto part = run_program(c, exact_cone);
            require(fmpq_poly_equal(whole.real(), part.real()) && fmpq_poly_equal(whole.imag(), part.imag()),
                    "cone: node " + std::to_string(node) + " of random program " + std::to_string(t) +
                        " computes the node's polynomial");
            require(c.degree_bounds()[*c.output] == bounds[node], "cone: the node's degree bound");
            ++cones;
        }
    }
    const auto power_of = [](const char *text) { return output_power(parse_program(text)); };
    const auto cube = power_of("PSPROGRAM 1\nx\npow 0 3\nout 1\n");
    require(cube && cube->base == 0 && cube->exponent == 3, "output power: x^3");
    const auto square = power_of("PSPROGRAM 1\nx\nc 2\nsub 0 1\nmul 2 2\nout 3\n");
    require(square && square->base == 2 && square->exponent == 2, "output power: (x - 2)^2");
    const auto nested = power_of("PSPROGRAM 1\nx\nc 1\nadd 0 1\nmul 2 2\npow 3 3\nout 4\n");
    require(nested && nested->base == 2 && nested->exponent == 6, "output power: ((x + 1)^2)^3");
    const auto scaled = power_of("PSPROGRAM 1\nx\nc 1\nadd 0 1\npow 2 2\nc 3\nmul 4 3\nneg 5\nout 6\n");
    require(scaled && scaled->base == 2 && scaled->exponent == 2, "output power: -(3 (x + 1)^2)");
    require(!power_of("PSPROGRAM 1\nx\nc 1\nadd 0 1\nc 3\nmul 2 3\nout 4\n"), "output power: 3 (x + 1) is none");
    // A product of distinct factors is not a power (the lift waits for a decision on products), nor is
    // Mandelbrot's output, a sum.
    require(!power_of("PSPROGRAM 1\nx\nc 1\nsub 0 1\nc 1\nsub 0 3\nmul 2 4\nc 2\nadd 0 6\nmul 5 7\nout 8\n"),
            "output power: (x - 1)(x - 1)(x + 2), two nodes for x - 1, is none");
    std::cout << "structure: " << cones << " cones of random programs computed their nodes' polynomials\n";
}

void cancellations() {
    // (x + 1)^2 - x^2 - 2x = 1: two leading coefficients cancel.
    {
        Program p;
        const auto x = node(p, ProgramOp::Variable);
        const auto one = p.append_constant(gaussian("1"));
        const auto two = p.append_constant(gaussian("2"));
        const auto s = node(p, ProgramOp::Add, x, one);
        const auto sq = node(p, ProgramOp::Power, s, 0, 2);
        const auto x2 = node(p, ProgramOp::Multiply, x, x);
        const auto d = node(p, ProgramOp::Subtract, sq, x2);
        const auto tx = node(p, ProgramOp::Multiply, two, x);
        p.output = node(p, ProgramOp::Subtract, d, tx);
        check_against_expansion(p, "(x+1)^2 - x^2 - 2x");
        const auto top = top_of_program(p, 64, unlimited);
        require(top.degree == 0 && top.width > 1, "cancellation needs a wider run");
    }
    // x^20 - x^20 = 0 and (x+1)^3 - (x^3 + 3x^2 + 3x + 1) = 0.
    {
        Program p;
        const auto x = node(p, ProgramOp::Variable);
        const auto a = node(p, ProgramOp::Power, x, 0, 20);
        p.output = node(p, ProgramOp::Subtract, a, a);
        require(top_of_program(p, 64, unlimited).verdict == ProgramTop::Verdict::ZeroPolynomial, "x^20 - x^20");
    }
    {
        const auto p = parse_program("PSPROGRAM 1\nx\nc 1\nadd 0 1\npow 2 3\npow 0 3\nc 3\npow 0 2\nmul 5 6\n"
                                     "mul 5 0\nadd 4 7\nadd 9 8\nadd 10 1\nsub 3 11\nout 12\n");
        check_against_expansion(p, "(x+1)^3 - expansion");
        require(top_of_program(p, 64, unlimited).verdict == ProgramTop::Verdict::ZeroPolynomial, "(x+1)^3 cancels");
    }
    // A width too small to decide is reported, not guessed.
    {
        const auto p = parse_program("PSPROGRAM 1\nx\nc 1\nadd 0 1\npow 2 9\npow 0 9\nsub 3 4\nout 5\n");
        require(top_of_program(p, 1, unlimited).verdict == ProgramTop::Verdict::Undetermined, "width 1 undetermined");
        require(top_of_program(p, 2, unlimited).degree == 8, "width 2 decides degree 8");
    }
}

void liveness() {
    // Chebyshev T_n by T_(k+1) = 2 x T_k - T_(k-1): about 3n nodes, a few live at once.
    Program p;
    const auto x = node(p, ProgramOp::Variable);
    const auto two = p.append_constant(gaussian("2"));
    uint64_t previous = p.append_constant(gaussian("1")), current = x;
    const int n = 300;
    for (int k = 1; k < n; ++k) {
        const auto twice = node(p, ProgramOp::Multiply, two, x);
        const auto product = node(p, ProgramOp::Multiply, twice, current);
        const auto next = node(p, ProgramOp::Subtract, product, previous);
        previous = current;
        current = next;
    }
    p.output = current;
    ComplexBall point;
    acb_set_d_d(point.value(), 0.3, 0.2);
    DualBallArithmetic dual(point.value(), 128);
    size_t peak = 0;
    run_program(p, dual, &peak);
    require(peak <= 8, "liveness: " + std::to_string(peak) + " values live at once for " +
                           std::to_string(p.nodes.size()) + " nodes");
    require(top_of_program(p, 4, unlimited).degree == static_cast<uint64_t>(n), "Chebyshev degree");
}

void text_form() {
    const auto p = parse_program("# Mandelbrot, three steps\nPSPROGRAM 1\nx\nc 1\nmul 0 0\nadd 2 1   # x^2+1\n"
                                 "mul 3 3\nmul 4 0\nadd 5 1\nout 6\n");
    require(p.nodes.size() == 7 && p.output == 6u, "parsed size and output");
    const char *bad[] = {
        "x\nout 0\n",                          // no header
        "PSPROGRAM 1\nx\nadd 0 1\nout 2\n",     // forward reference
        "PSPROGRAM 1\nx\npow 0 1\nout 1\n",     // exponent below 2
        "PSPROGRAM 1\nx\nfoo 0\nout 0\n",       // unknown operation
        "PSPROGRAM 1\nx\nc 1 2 3\nout 0\n",     // too many fields
        "PSPROGRAM 1\nx\n",                     // no output
        "PSPROGRAM 1\nx\nc 1/0\nout 1\n",       // zero denominator
    };
    for (const auto *text : bad) {
        bool refused = false;
        try {
            parse_program(text);
        } catch (const std::invalid_argument &) {
            refused = true;
        }
        require(refused, std::string("accepted: ") + text);
    }
}

void prime() {
    for (unsigned j = 0; j < 3; ++j) {
        const auto g = gaussian_prime(j);
        require(g.prime % 4 == 1 && n_is_prime(g.prime) && g.prime > (UWORD(1) << 62), "Gaussian prime");
        nmod_t mod;
        nmod_init(&mod, g.prime);
        require(nmod_mul(g.imaginary_unit, g.imaginary_unit, mod) == g.prime - 1, "square root of -1");
        if (j) require(g.prime > gaussian_prime(j - 1).prime, "primes increase");
    }
}

void c_api() {
    cps_program *p = nullptr;
    require(cps_program_create(&p) == CPS_OK && p, "create");
    uint64_t x = 0, half = 0, product = 0, cube = 0, bad = 0;
    require(cps_program_variable(p, &x) == CPS_OK && x == 0, "variable");
    require(cps_program_constant(p, "1/2", "3", &half) == CPS_OK && half == 1, "constant");
    require(cps_program_multiply(p, x, half, &product) == CPS_OK && product == 2, "multiply");
    require(cps_program_power(p, product, 3, &cube) == CPS_OK && cube == 3, "power");
    require(cps_program_add(p, x, 9, &bad) == CPS_INVALID_ARGUMENT && cps_program_size(p) == 4, "forward operand refused");
    require(cps_program_power(p, x, 1, &bad) == CPS_INVALID_ARGUMENT && cps_program_size(p) == 4, "exponent 1 refused");
    require(cps_program_constant(p, "1/0", nullptr, &bad) == CPS_INVALID_ARGUMENT && cps_program_size(p) == 4,
            "bad constant refused");
    require(cps_program_output(p) == CPS_NO_NODE, "no output yet");
    require(cps_program_set_output(p, 7) == CPS_INVALID_ARGUMENT, "output out of range refused");
    require(cps_program_set_output(p, cube) == CPS_OK && cps_program_output(p) == 3, "output");
    cps_program_destroy(p);
    cps_program *q = nullptr;
    require(cps_program_parse("PSPROGRAM 1\nx\nc 1\nadd 0 1\nout 2\n", &q) == CPS_OK && cps_program_size(q) == 3, "parse");
    cps_program_destroy(q);
    require(cps_program_parse("nonsense", &q) == CPS_INVALID_ARGUMENT && !q, "parse refusal");
}

// Whether a type can be compared or tested: what the recording type must not allow (a constraint, so that an
// invalid expression is false rather than an error).
template <class T>
constexpr bool recorder_comparisons = requires(const T &a) { a < a; } || requires(const T &a) { a == a; } ||
                                      requires(const T &a) { a > 0; } || requires(const T &a) { static_cast<bool>(a); };
static_assert(recorder_comparisons<int>, "the constraint recognises a comparable type");

// The C++ recorder: a generic lambda recorded once equals the hand-written program.
void recorder() {
    const auto mandelbrot = certpolysolve::Program::from([](auto x) {
        auto p = x * 0 + 1;
        for (int k = 0; k < 4; ++k) p = x * p * p + 1;
        return p;
    });
    ExpansionArithmetic exact(unlimited);
    const auto recorded = run_program(mandelbrot.get()->program, exact);
    // p_(k+1) = x p_k^2 + 1 from p_0 = 1, expanded directly.
    GaussianExpansion expected;
    fmpq_poly_set_si(expected.real(), 1);
    for (int k = 0; k < 4; ++k) {
        GaussianExpansion x;
        fmpq_poly_set_coeff_si(x.real(), 1, 1);
        ExpansionArithmetic e(unlimited);
        auto square = e.multiply(expected, expected);
        auto next = e.multiply(x, square);
        GaussianExpansion one;
        fmpq_poly_set_si(one.real(), 1);
        expected = e.add(next, one);
    }
    require(fmpq_poly_equal(recorded.real(), expected.real()) && fmpq_poly_is_zero(recorded.imag()),
            "recorded Mandelbrot differs");
    require(recorded.degree() == 15, "Mandelbrot p_4 has degree 15");
    // A double enters at its exact value: 0.1 is 3602879701896397 / 2^55.
    const auto tenth = certpolysolve::Program::from(
        [](auto x) { return x * 0.1 + certpolysolve::Number::complex(certpolysolve::Number("1/3"), certpolysolve::Number("-2")); });
    const auto top = top_of_program(tenth.get()->program, 4, unlimited);
    Rational expected_lead;
    fmpq_set_si(expected_lead.value(), 3602879701896397, 1);
    fmpq_div_2exp(expected_lead.value(), expected_lead.value(), 55);
    require(top.verdict == ProgramTop::Verdict::Determined && top.degree == 1 &&
                fmpq_equal(top.leading.real.value(), expected_lead.value()), "0.1 recorded exactly");
    // Powers, negation, compound assignment, a constant result.
    const auto mixed = certpolysolve::Program::from([](auto x) {
        auto p = pow(x - 2, 3);
        p *= -x;
        p += 7;
        return p;
    });
    ExpansionArithmetic e(unlimited);
    const auto m = run_program(mixed.get()->program, e);
    Rational c;
    fmpq_poly_get_coeff_fmpq(c.value(), m.real(), 0);
    require(m.degree() == 4 && fmpq_equal_si(c.value(), 7), "-(x-2)^3 x + 7");
    const auto constant = certpolysolve::Program::from([](auto) { return 5; });
    require(constant.ok() && constant.get()->program.nodes.size() == 2 && constant.get()->program.output == 1u,
            "a constant result is a constant node");
    // The recording type has no comparisons, so a program cannot branch on the variable's value (spec, section
    // 3.2): the negative-compilation test as constraints, checked when this file compiles (recorder_comparisons).
    static_assert(!recorder_comparisons<certpolysolve::ProgramValue>, "ProgramValue must not compare or convert to bool");
}

// Binary64 value and derivative agree with the balls on a short program.
void binary64() {
    const auto p = parse_program("PSPROGRAM 1\nx\nc 1\nmul 0 0\nadd 2 1\npow 3 3\nc 1/3 2\nmul 4 5\nout 6\n");
    const std::complex<double> z(0.25, -0.5);
    ComplexDoubleArithmetic d(p, z);
    const auto v = run_program(p, d);
    ComplexBall point;
    acb_set_d_d(point.value(), z.real(), z.imag());
    DualBallArithmetic dual(point.value(), 128);
    const auto w = run_program(p, dual);
    const auto close = [](std::complex<double> a, const acb_t b) {
        const std::complex<double> c(arf_get_d(arb_midref(acb_realref(b)), ARF_RND_NEAR),
                                     arf_get_d(arb_midref(acb_imagref(b)), ARF_RND_NEAR));
        return std::abs(a - c) <= 1e-13 * std::max(1.0, std::abs(c));
    };
    require(close(v.value, w.value.value()) && close(v.derivative, w.derivative.value()), "binary64 arithmetic");
}
} // namespace

// ------------------------------------------------------------------- milestone 2: the preparation

// The largest root modulus from below, by the expansion: |a_(n-k) / a_n| <= C(n, k) max|r|^k (Vieta), so
// max|r| >= (|a_(n-k)| / (C(n, k) |a_n|))^(1/k) for every k. A root bound under it is wrong.
double largest_root_lower_bound(const GaussianExpansion &e) {
    const slong n = e.degree();
    const auto modulus = [&](slong j) {
        Rational re, im;
        fmpq_poly_get_coeff_fmpq(re.value(), e.real(), j);
        fmpq_poly_get_coeff_fmpq(im.value(), e.imag(), j);
        return std::hypot(fmpq_get_d(re.value()), fmpq_get_d(im.value()));
    };
    const double lead = modulus(n);
    double best = 0, binomial = 1;
    for (slong k = 1; k <= n; ++k) {
        binomial = binomial * static_cast<double>(n - k + 1) / static_cast<double>(k);
        if (const double c = modulus(n - k); c > 0)
            best = std::max(best, std::pow(c / (binomial * lead), 1.0 / static_cast<double>(k)));
    }
    return best;
}
// The root bound of a prepared program is at least `largest` (a root's modulus, or a lower bound of the largest);
// returns it.
double check_root_bound(const ProgramPreparation &prepared, double largest, const std::string &label) {
    require(prepared.ready(), label + ": refused (" + prepared.diagnostic + ")");
    const double bound = mag_get_d(prepared.root_bound.value()); // rounded up
    require(bound >= largest * (1 - 1e-9), label + ": root bound " + std::to_string(bound) + " below " +
                                               std::to_string(largest));
    return bound;
}
// prod (x - r_j) over Gaussian rationals.
Program from_roots(const std::vector<std::pair<const char *, const char *>> &roots) {
    Program p;
    const auto x = node(p, ProgramOp::Variable);
    std::optional<uint64_t> product;
    for (const auto &[re, im] : roots) {
        const auto factor = node(p, ProgramOp::Subtract, x, p.append_constant(gaussian(re, im)));
        product = product ? node(p, ProgramOp::Multiply, *product, factor) : factor;
    }
    p.output = *product;
    return p;
}
// The same polynomial in Horner form, from its exact expanded coefficients: an output that is a sum, not a product,
// so the product lift (output_factors) leaves it to the numerical pipeline the tests exercise.
Program horner_from_roots(const std::vector<std::pair<const char *, const char *>> &roots) {
    ExpansionArithmetic exact(unlimited);
    auto e = run_program(from_roots(roots), exact);
    const auto coefficient = [&](slong j) {
        Rational re, im;
        fmpq_poly_get_coeff_fmpq(re.value(), e.real(), j);
        fmpq_poly_get_coeff_fmpq(im.value(), e.imag(), j);
        char *a = fmpq_get_str(nullptr, 10, re.value());
        char *b = fmpq_get_str(nullptr, 10, im.value());
        Coefficient c = gaussian(a, b);
        flint_free(a);
        flint_free(b);
        return c;
    };
    Program p;
    const auto x = node(p, ProgramOp::Variable);
    uint64_t acc = p.append_constant(coefficient(e.degree()));
    for (slong j = e.degree() - 1; j >= 0; --j)
        acc = node(p, ProgramOp::Add, node(p, ProgramOp::Multiply, acc, x), p.append_constant(coefficient(j)));
    p.output = acc;
    return p;
}

void refusals() {
    const cps_options options;
    using Refusal = ProgramPreparation::Refusal;
    const auto refused = [&](const char *text, Refusal refusal, const std::string &expected) {
        const auto prepared = prepare_program(parse_program(text), options);
        require(!prepared.ready() && prepared.refusal == refusal && prepared.diagnostic.find(expected) != std::string::npos,
                std::string("refusal '") + expected + "' expected, got '" + prepared.diagnostic + "'");
    };
    refused("PSPROGRAM 1\nx\nsub 0 0\nout 1\n", Refusal::ZeroPolynomial, "identically zero");
    // (x^100 + x) - x^100 = x: its leading terms cancel over 99 coefficients, past the widest top series.
    refused("PSPROGRAM 1\nx\npow 0 100\nadd 1 0\nsub 2 1\nout 3\n", Refusal::LeadingUndetermined,
            "leading coefficient not determined");
    // (x + 1) - x = 1: a nonzero constant is ready with degree 0 (no roots), not refused.
    const auto constant = prepare_program(parse_program("PSPROGRAM 1\nx\nc 1\nadd 0 1\nsub 2 0\nout 3\n"), options);
    require(constant.ready() && constant.degree == 0 && fmpq_is_one(constant.leading.real.value()), "a constant program");
    // (x^50 + x) - x^50 = x is within it.
    const auto within = prepare_program(parse_program("PSPROGRAM 1\nx\npow 0 50\nadd 1 0\nsub 2 1\nout 3\n"), options);
    require(within.ready() && within.degree == 1, "(x^50 + x) - x^50 is x");
    cps_options short_budget;
    short_budget.exact_operations = 3;
    const auto longer = prepare_program(parse_program("PSPROGRAM 1\nx\nc 1\nadd 0 1\nmul 2 2\nout 3\n"), short_budget);
    require(!longer.ready() && longer.refusal == Refusal::Length &&
                longer.diagnostic.find("exact operations") != std::string::npos,
            "length refusal");
}

void modular_verdicts() {
    const cps_options options;
    const auto simple = prepare_program(from_roots({{"1", "0"}, {"2", "0"}, {"3", "1/2"}}), options);
    require(simple.ready() && simple.modular.squarefree, "three distinct roots: " + simple.modular.note);
    const auto repeated = prepare_program(from_roots({{"1", "0"}, {"1", "0"}, {"-2", "0"}}), options);
    require(repeated.ready() && !repeated.modular.squarefree &&
                repeated.modular.note.find("not square-free modulo 2 primes") != std::string::npos,
            "a double root: " + repeated.modular.note);
    // A constant whose denominator the first prime divides, and a leading coefficient it divides: the second prime.
    const auto first = std::to_string(gaussian_prime(0).prime), second = std::to_string(gaussian_prime(1).prime);
    const auto denominator = prepare_program(from_roots({{("1/" + first).c_str(), "0"}, {"2", "0"}}), options);
    require(denominator.modular.squarefree && denominator.modular.note.find(second) != std::string::npos,
            "a denominator divisible by the first prime: " + denominator.modular.note);
    const auto leading = prepare_program(parse_program("PSPROGRAM 1\nx\nc " + first + "\nmul 0 1\nc 1\nadd 2 3\n"
                                                       "mul 4 0\nout 5\n"), options);
    require(leading.ready() && leading.degree == 2 && leading.modular.squarefree &&
                leading.modular.note.find(second) != std::string::npos,
            "a leading coefficient divisible by the first prime: " + leading.modular.note);
    cps_options priced;
    priced.exact_operations = 64;
    const auto skipped = prepare_program(from_roots({{"1", "0"}, {"2", "0"}, {"3", "0"}, {"4", "0"}}), priced);
    require(skipped.ready() && !skipped.modular.squarefree &&
                skipped.modular.note.find("skipped") != std::string::npos,
            "the cost rule: " + skipped.modular.note);
}

void root_bounds() {
    const cps_options options;
    // x - 1000: Cauchy's M(1) / |a_n| = 1001 against the root's 1000.
    const double far = check_root_bound(prepare_program(from_roots({{"1000", "0"}}), options), 1000, "x - 1000");
    require(far <= 1001 * (1 + 1e-6), "x - 1000: bound " + std::to_string(far)); // 30-bit magnitudes
    // Known roots of several sizes, complex ones included.
    check_root_bound(prepare_program(from_roots({{"1/3", "0"}, {"-5/2", "7"}, {"0", "-1/1024"}, {"9", "9"}}), options),
                     std::hypot(9.0, 9.0), "four roots");
    // Mandelbrot p_(k+1) = p_k^2 + c (roots in |c| <= 2). Its coefficients are nonnegative, so the majorant is the
    // polynomial itself and every bound from coefficient moduli is about the degree: M(R) < 2 R^n needs
    // (1 + 1/R)^(n/2) < 2, R > n / (2 ln 2), 369 at n = 512; the least power of two is 512.
    {
        Program p;
        const auto c = node(p, ProgramOp::Variable);
        uint64_t q = c;
        for (int k = 1; k < 10; ++k) q = node(p, ProgramOp::Add, node(p, ProgramOp::Multiply, q, q), c);
        p.output = q;
        ExpansionArithmetic exact(unlimited);
        const auto e = run_program(p, exact);
        const auto prepared = prepare_program(p, options);
        require(prepared.degree == 512 && prepared.modular.squarefree, "Mandelbrot p_10: " + prepared.modular.note);
        const double bound = check_root_bound(prepared, largest_root_lower_bound(e), "Mandelbrot p_10");
        require(bound <= 512 * (1 + 1e-6), "Mandelbrot p_10: bound " + std::to_string(bound));
        std::cout << "root bound: Mandelbrot p_10 (degree 512) " << bound << '\n';
    }
    // (x + 1)^n - x^n: its leading terms cancel, so the test at a power of two never holds and Cauchy's bound
    // stands; the roots 1 / (w^k - 1) reach 1 / (2 sin(pi / n)).
    for (const uint64_t n : {8u, 32u}) {
        Program p;
        const auto x = node(p, ProgramOp::Variable);
        const auto s = node(p, ProgramOp::Add, x, p.append_constant(gaussian("1")));
        p.output = node(p, ProgramOp::Subtract, node(p, ProgramOp::Power, s, 0, n), node(p, ProgramOp::Power, x, 0, n));
        const double largest = 1 / (2 * std::sin(std::acos(-1.0) / static_cast<double>(n)));
        const double bound = check_root_bound(prepare_program(p, options), largest, "(x+1)^n - x^n");
        std::cout << "root bound: (x+1)^" << n << " - x^" << n << " " << bound << " (largest root " << largest << ")\n";
    }
    // Random programs: never below the expansion's lower bound of the largest root.
    std::mt19937_64 rng(20261016);
    int checked = 0;
    for (int t = 0; t < 300; ++t) {
        const auto p = random_program(rng, 10 + t % 30, 40);
        const auto prepared = prepare_program(p, options);
        if (!prepared.ready()) continue;
        ExpansionArithmetic exact(unlimited);
        check_root_bound(prepared, largest_root_lower_bound(run_program(p, exact)), "random program " + std::to_string(t));
        ++checked;
    }
    // At least the 50 of degree 10 or more that random_programs requires (no refusal applies to them: a width of
    // 64 decides any degree below 64).
    require(checked >= 50, "random programs prepared: " + std::to_string(checked));
    std::cout << "root bounds: " << checked << " random programs\n";
}

// ------------------------------------------------------------------------- milestone 3: the solve

namespace {
// A program solved through the public entry (cps_solve_program), as a caller would: the C options the 0.5 header's
// solve sets from `options` (certpolysolve_test::COptions), the answer its Result. An error code fails the test, as the
// pre-0.5 header's throw did.
certpolysolve::Result solve_program(const Program &p, const certpolysolve::Options &options) {
    cps_program handle;
    handle.program = p;
    const certpolysolve_test::COptions settings(options);
    require(settings.ok(), "the options were refused: " + settings.status());
    cps_result *result = nullptr;
    const cps_error error = cps_solve_program(&handle, settings.get(), &result);
    require(error == CPS_OK, std::string("cps_solve_program failed: ") + cps_last_error());
    return certpolysolve_test::from_c(error, result);
}
certpolysolve::Result solve_program(const certpolysolve::Program &p, const certpolysolve::Options &options) {
    auto answer = certpolysolve::solve(p, options);
    require(answer.stop() != certpolysolve::Stop::error, "the solve failed: " + answer.message());
    return answer;
}
// The same with C options (api_test_support.hpp): the settings the 0.5 Options leaves out (a cancellation callback,
// certpolysolve_testing.h's poll counter, real-only and region queries a program refuses).
certpolysolve::Result solve_program(const certpolysolve::Program &p, const cps_options *options) {
    auto answer = certpolysolve_test::solve(p, options);
    require(answer.stop() != certpolysolve::Stop::error, "the solve failed: " + answer.message());
    return answer;
}
certpolysolve::Result solve_program(const Program &p, const cps_options *options) {
    cps_program handle;
    handle.program = p;
    cps_result *result = nullptr;
    const cps_error error = cps_solve_program(&handle, options, &result);
    require(error == CPS_OK, std::string("cps_solve_program failed: ") + cps_last_error());
    return certpolysolve_test::from_c(error, result);
}
// C options as the 0.5 header's solve sets them from `options`, with certpolysolve_testing.h's poll counter armed: a stop
// at the k-th cancellation poll (k >= 1), or counting only (k = UINT64_MAX). Until stage 3b a callback counting its
// calls made those stops; the counter polls where it was called.
certpolysolve_test::COptions stop_at_poll(const certpolysolve::Options &options, uint64_t k) {
    certpolysolve_test::COptions settings(options);
    require(settings.ok() && cps_testing_options_stop_at_poll(settings.get(), k) == CPS_OK, "the poll counter");
    return settings;
}
// The coefficient twin: the program expanded exactly (tests only) and given to the coefficient path.
certpolysolve::detail::DensePolynomial twin_of(const Program &p) {
    ExpansionArithmetic exact(unlimited);
    const auto e = run_program(p, exact);
    certpolysolve::detail::DensePolynomial twin;
    const auto text = [](const fmpq_poly_struct *part, slong j) {
        Rational c;
        fmpq_poly_get_coeff_fmpq(c.value(), part, j);
        char *s = fmpq_get_str(nullptr, 10, c.value());
        std::string out(s);
        flint_free(s);
        return out;
    };
    for (slong j = 0; j <= e.degree(); ++j) {
        Coefficient c;
        c.real = Number(text(e.real(), j));
        c.imag = Number(text(e.imag(), j));
        twin.coefficients.push_back(std::move(c));
    }
    return twin;
}
certpolysolve::Result solve_twin(const certpolysolve::detail::DensePolynomial &twin, const certpolysolve::Options &options) {
    cps_polynomial handle; // the public input object holds the twin's nonzero terms (C40)
    for (size_t k = 0; k < twin.coefficients.size(); ++k)
        if (!twin.coefficients[k].zero()) handle.terms.emplace(k, twin.coefficients[k]);
    const certpolysolve_test::COptions settings(options);
    require(settings.ok(), "the options were refused: " + settings.status());
    cps_result *result = nullptr;
    const cps_error error = cps_solve(&handle, settings.get(), &result);
    require(error == CPS_OK, std::string("cps_solve failed: ") + cps_last_error());
    return certpolysolve_test::from_c(error, result);
}
Rational endpoint(const certpolysolve::Dyadic &d) {
    Rational value;
    Integer mantissa;
    require(fmpz_set_str(mantissa.value(), d.mantissa().c_str(), 10) == 0, "endpoint mantissa");
    fmpq_set_fmpz(value.value(), mantissa.value());
    const long long e = std::stoll(d.exponent());
    if (e >= 0) fmpq_mul_2exp(value.value(), value.value(), static_cast<ulong>(e));
    else fmpq_div_2exp(value.value(), value.value(), static_cast<ulong>(-e));
    return value;
}
struct Box {
    Rational lo, hi, ilo, ihi;
    uint64_t count = 0, multiplicity = 0;
    cps_realness realness = CPS_REAL_UNKNOWN;
    explicit Box(const certpolysolve::Root &r)
        : lo(endpoint(r.box.real.lo)), hi(endpoint(r.box.real.hi)), ilo(endpoint(r.box.imag.lo)),
          ihi(endpoint(r.box.imag.hi)), count(r.count), multiplicity(r.multiplicity.value_or(0)),
          realness(r.realness == certpolysolve::Realness::real       ? CPS_REAL
                   : r.realness == certpolysolve::Realness::not_real ? CPS_NONREAL
                                                                 : CPS_REAL_UNKNOWN) {}
    bool contains(const Rational &x, const Rational &y) const {
        return fmpq_cmp(lo.value(), x.value()) <= 0 && fmpq_cmp(x.value(), hi.value()) <= 0 &&
               fmpq_cmp(ilo.value(), y.value()) <= 0 && fmpq_cmp(y.value(), ihi.value()) <= 0;
    }
    bool overlaps(const Box &o) const {
        return fmpq_cmp(hi.value(), o.lo.value()) >= 0 && fmpq_cmp(o.hi.value(), lo.value()) >= 0 &&
               fmpq_cmp(ihi.value(), o.ilo.value()) >= 0 && fmpq_cmp(o.ihi.value(), ilo.value()) >= 0;
    }
};
std::vector<Box> boxes_of(const certpolysolve::Result &result) {
    std::vector<Box> out;
    for (const auto &r : result.roots()) out.emplace_back(r);
    return out;
}
// The exported cover's own contract: the counts sum to the degree, the rectangles are disjoint, and under the
// accuracy goal every region is one root (multiplicity 1) whose half-width is within 2^-bits of its modulus
// (exact_reference.hpp's check_accuracy, in FLINT rationals).
void check_cover(const std::vector<Box> &boxes, uint64_t degree, const certpolysolve::Result &result, uint64_t bits,
                 const std::string &label) {
    uint64_t total = 0;
    for (size_t i = 0; i < boxes.size(); ++i) {
        total += boxes[i].count;
        for (size_t j = 0; j < i; ++j)
            require(!boxes[i].overlaps(boxes[j]), label + ": exported rectangles overlap");
        if (cps_result_goal(result.get()) == CPS_ROOT_ACCURACY) {
            const auto &b = boxes[i];
            require(b.count == 1 && b.multiplicity == 1, label + ": accuracy goal met with a counted cluster");
            Rational error, lower, limit, candidate;
            fmpq_sub(error.value(), b.hi.value(), b.lo.value());
            fmpq_sub(candidate.value(), b.ihi.value(), b.ilo.value());
            fmpq_add(error.value(), error.value(), candidate.value());
            fmpq_div_2exp(error.value(), error.value(), 1);
            // The lower bound of |Re z| over the box, then of |Im z|; the larger serves.
            const auto lower_abs = [](const Rational &a, const Rational &c) {
                Rational out;
                if (fmpq_sgn(a.value()) > 0) fmpq_set(out.value(), a.value());
                else if (fmpq_sgn(c.value()) < 0) fmpq_neg(out.value(), c.value());
                return out;
            };
            lower = lower_abs(b.lo, b.hi);
            candidate = lower_abs(b.ilo, b.ihi);
            if (fmpq_cmp(candidate.value(), lower.value()) > 0) lower = candidate;
            fmpq_div_2exp(limit.value(), lower.value(), static_cast<ulong>(bits));
            require(fmpq_cmp(error.value(), limit.value()) <= 0, label + ": exported accuracy contract failed");
        }
    }
    require(total == degree, label + ": counts sum to " + std::to_string(total) + ", not the degree " +
                                 std::to_string(degree));
}
// Two certified covers of the same roots at the same goal agree (the gate's cross-check, in exact arithmetic):
// in the graph whose edges are overlapping rectangles across the covers, every component has the same count on
// both sides (each root lies in one region of each cover, and those two overlap), and every region's midpoint,
// which approximates each of its roots to 2^-bits of the root's modulus, is within that distance of a region
// of the other cover, componentwise (the L1 bound of the rectangle's corners bounds the modulus).
void cross_check(const certpolysolve::Result &a, const certpolysolve::Result &b, uint64_t degree, uint64_t bits,
                 const std::string &label) {
    require(a.complete() && a.goal_met(), label + ": the program's answer is incomplete: " + a.message());
    require(b.complete() && b.goal_met(), label + ": the twin's answer is incomplete: " + b.message());
    require(cps_result_goal(a.get()) == cps_result_goal(b.get()), label + ": the twins were asked different goals");
    const auto boxes_a = boxes_of(a), boxes_b = boxes_of(b);
    check_cover(boxes_a, degree, a, bits, label + " (program)");
    check_cover(boxes_b, degree, b, bits, label + " (twin)");
    // Components by union-find over both covers' regions.
    std::vector<size_t> parent(boxes_a.size() + boxes_b.size());
    for (size_t i = 0; i < parent.size(); ++i) parent[i] = i;
    const auto find = [&](size_t i) {
        while (parent[i] != i) i = parent[i] = parent[parent[i]];
        return i;
    };
    for (size_t i = 0; i < boxes_a.size(); ++i)
        for (size_t j = 0; j < boxes_b.size(); ++j)
            if (boxes_a[i].overlaps(boxes_b[j])) parent[find(i)] = find(boxes_a.size() + j);
    std::vector<int64_t> balance(parent.size(), 0);
    std::vector<uint8_t> sides(parent.size(), 0);
    for (size_t i = 0; i < boxes_a.size(); ++i) {
        balance[find(i)] += static_cast<int64_t>(boxes_a[i].count);
        sides[find(i)] |= 1;
    }
    for (size_t j = 0; j < boxes_b.size(); ++j) {
        balance[find(boxes_a.size() + j)] -= static_cast<int64_t>(boxes_b[j].count);
        sides[find(boxes_a.size() + j)] |= 2;
    }
    for (size_t i = 0; i < parent.size(); ++i)
        if (find(i) == i) {
            require(sides[i] == 3, label + ": a region of one cover meets no region of the other");
            require(balance[i] == 0, label + ": a group of overlapping regions counts differently in the two covers");
        }
    const auto near = [&](const std::vector<Box> &from, const std::vector<Box> &to, const char *which) {
        for (const auto &box : from) {
            Rational x, y;
            fmpq_add(x.value(), box.lo.value(), box.hi.value());
            fmpq_div_2exp(x.value(), x.value(), 1);
            fmpq_add(y.value(), box.ilo.value(), box.ihi.value());
            fmpq_div_2exp(y.value(), y.value(), 1);
            bool found = false;
            for (const auto &other : to) {
                // The modulus bound: the larger |real corner| plus the larger |imaginary corner|.
                Rational bound, part, tolerance, dx, dy;
                fmpq_abs(bound.value(), other.lo.value());
                fmpq_abs(part.value(), other.hi.value());
                if (fmpq_cmp(part.value(), bound.value()) > 0) bound = part;
                fmpq_abs(part.value(), other.ilo.value());
                Rational imag_part;
                fmpq_abs(imag_part.value(), other.ihi.value());
                if (fmpq_cmp(imag_part.value(), part.value()) > 0) part = imag_part;
                fmpq_add(bound.value(), bound.value(), part.value());
                fmpq_div_2exp(tolerance.value(), bound.value(), static_cast<ulong>(bits));
                // The distance from the midpoint to the rectangle, in each coordinate.
                const auto distance = [](Rational &out, const Rational &v, const Rational &lo, const Rational &hi) {
                    if (fmpq_cmp(v.value(), lo.value()) < 0) fmpq_sub(out.value(), lo.value(), v.value());
                    else if (fmpq_cmp(v.value(), hi.value()) > 0) fmpq_sub(out.value(), v.value(), hi.value());
                    else fmpq_zero(out.value());
                };
                distance(dx, x, other.lo, other.hi);
                distance(dy, y, other.ilo, other.ihi);
                if (fmpq_cmp(dx.value(), tolerance.value()) <= 0 && fmpq_cmp(dy.value(), tolerance.value()) <= 0) {
                    found = true;
                    break;
                }
            }
            require(found, label + ": a " + which + " region's midpoint is not within the goal's distance of the other cover");
        }
    };
    near(boxes_a, boxes_b, "program");
    near(boxes_b, boxes_a, "twin");
}
// A known-root check (test_structured.cpp's `known`): every region holds exactly its count of the planted roots.
void check_known(const certpolysolve::Result &result, const std::vector<std::pair<const char *, const char *>> &roots,
                 uint64_t bits, const std::string &label) {
    require(result.complete() && result.goal_met(), label + ": incomplete: " + result.message());
    const auto boxes = boxes_of(result);
    check_cover(boxes, roots.size(), result, bits, label);
    uint64_t total = 0;
    for (const auto &box : boxes) {
        uint64_t inside = 0;
        for (const auto &[re, im] : roots) {
            Rational x, y;
            fmpq_set_str(x.value(), re, 10);
            fmpq_set_str(y.value(), im, 10);
            inside += box.contains(x, y);
        }
        require(inside == box.count, label + ": a region holds " + std::to_string(inside) + " planted roots, counts " +
                                         std::to_string(box.count));
        total += inside;
    }
    require(total == roots.size(), label + ": roots lost");
}
// Three-term recurrences (the real corpus's families), by the recorder.
certpolysolve::Program chebyshev(int n) {
    return certpolysolve::Program::from([n](auto x) {
        auto previous = x * 0 + 1, current = x;
        for (int k = 1; k < n; ++k) {
            auto next = 2 * x * current - previous;
            previous = current;
            current = next;
        }
        return current;
    });
}
certpolysolve::Program hermite(int n) {
    return certpolysolve::Program::from([n](auto x) {
        auto previous = x * 0 + 1, current = 2 * x;
        for (int k = 1; k < n; ++k) {
            auto next = 2 * x * current - (2 * k) * previous;
            previous = current;
            current = next;
        }
        return current;
    });
}
certpolysolve::Program laguerre(int n) {
    return certpolysolve::Program::from([n](auto x) {
        auto previous = x * 0 + 1, current = 1 - x;
        for (int k = 1; k < n; ++k) {
            auto next = ((2 * k + 1 - x) * current - k * previous) * certpolysolve::Number("1/" + std::to_string(k + 1));
            previous = current;
            current = next;
        }
        return current;
    });
}
// det(x I - T) for a symmetric tridiagonal T with diagonal a_k = (k mod 3) - 1 and off-diagonal b_k = 1 + (k mod 2) / 2:
// p_0 = 1, p_1 = x - a_1, p_k = (x - a_k) p_(k-1) - b_(k-1)^2 p_(k-2).
certpolysolve::Program tridiagonal(int n) {
    return certpolysolve::Program::from([n](auto x) {
        auto previous = x * 0 + 1, current = x - ((1 % 3) - 1);
        for (int k = 2; k <= n; ++k) {
            const int a = (k % 3) - 1;
            const certpolysolve::Number b_squared((k - 1) % 2 ? "9/4" : "1");
            auto next = (x - a) * current - b_squared * previous;
            previous = current;
            current = next;
        }
        return current;
    });
}

// The batched binary64 evaluator (program_machine.hpp): the cone and the slots, against the balls, its range,
// its power's exactness.
std::complex<double> midpoint_of(const acb_t b) {
    return {arf_get_d(arb_midref(acb_realref(b)), ARF_RND_NEAR), arf_get_d(arb_midref(acb_imagref(b)), ARF_RND_NEAR)};
}
std::complex<double> standard(machine::Complex<double> z) { return {z.re, z.im}; }
bool same_bits(double a, double b) {
    uint64_t x, y;
    std::memcpy(&x, &a, sizeof x);
    std::memcpy(&y, &b, sizeof y);
    return x == y;
}
bool same_bits(const ProgramMachine<double>::Evaluation &a, const ProgramMachine<double>::Evaluation &b) {
    return same_bits(a.value.re, b.value.re) && same_bits(a.value.im, b.value.im) &&
           same_bits(a.derivative.re, b.derivative.re) && same_bits(a.derivative.im, b.derivative.im) &&
           same_bits(a.noise, b.noise) && same_bits(a.noise_derivative, b.noise_derivative) && a.finite == b.finite;
}
void machine_evaluator() {
    // The short program of binary64(), against the balls at 128 bits.
    {
        const auto p = parse_program("PSPROGRAM 1\nx\nc 1\nmul 0 0\nadd 2 1\npow 3 3\nc 1/3 2\nmul 4 5\nout 6\n");
        ProgramMachine<double> m(p);
        require(m.valid() && m.terms() == 7, "machine evaluator: the short program's cone has 7 operations");
        const std::complex<double> z(0.25, -0.5);
        const auto e = m.evaluate_one({z.real(), z.imag()});
        ComplexBall point;
        acb_set_d_d(point.value(), z.real(), z.imag());
        DualBallArithmetic dual(point.value(), 128);
        const auto w = run_program(p, dual);
        const auto close = [](machine::Complex<double> a, const acb_t b) {
            const auto c = midpoint_of(b);
            return std::abs(standard(a) - c) <= 1e-13 * std::max(1.0, std::abs(c));
        };
        require(e.finite && close(e.value, w.value.value()) && close(e.derivative, w.derivative.value()),
                "machine evaluator: the short program against the balls");
        require(e.noise > 0 && e.noise < 1e-12 && e.noise_derivative > 0 && e.noise_derivative < 1e-11,
                "machine evaluator: a short program's noise is a few roundings");
    }
    // Range: x^2000 at 2 overflows binary64 and is finite in the scaled arithmetic, at the right magnitude.
    {
        const auto p = parse_program("PSPROGRAM 1\nx\npow 0 2000\nout 1\n");
        ProgramMachine<double> m(p);
        require(m.valid() && !m.evaluate_one({2.0, 0.0}).finite, "machine evaluator: x^2000 at 2 is not finite in binary64");
        ProgramMachine<machine::BlockDouble> s(p);
        const auto f = s.evaluate_one({machine::BlockDouble(2.0), machine::BlockDouble(0.0)});
        require(f.finite && std::abs(machine::logabs(f.value) - 2000.0) < 1e-9 &&
                    std::abs(machine::logabs(f.derivative) - (1999.0 + std::log2(2000.0))) < 1e-9,
                "machine evaluator: x^2000 at 2 in the scaled arithmetic");
    }
    // The power's binary powering never squares past the last exponent bit: x^3 at 2^300 fits binary64
    // (2^900), where a square past the last bit (x^4 = 2^1200) would not.
    {
        const auto p = parse_program("PSPROGRAM 1\nx\npow 0 3\nout 1\n");
        ProgramMachine<double> m(p);
        const auto e = m.evaluate_one({0x1p300, 0.0});
        require(e.finite && e.value.re == 0x1p900 && e.value.im == 0 && e.derivative.re == 3 * 0x1p600,
                "machine evaluator: x^3 at 2^300 is finite in binary64");
    }
    // An exact zero: x * x at 0.
    {
        const auto p = parse_program("PSPROGRAM 1\nx\nmul 0 0\nout 1\n");
        ProgramMachine<double> m(p);
        const auto e = m.evaluate_one({0.0, 0.0});
        require(e.finite && e.value == machine::Complex<double>{} && e.derivative == machine::Complex<double>{} &&
                    e.noise == 0 && e.noise_derivative == 0,
                "machine evaluator: x * x at 0 is an exact zero");
    }
    // The cone leaves dead nodes out, and the slot plan keeps a recurrence to a few arrays: Chebyshev T_300
    // (the liveness program) in at most 9 slots; a dead overflowing node does not poison the output.
    {
        ProgramMachine<double> m(detail_of(chebyshev(300)));
        require(m.valid() && m.slots() <= 9, "machine evaluator: Chebyshev T_300 in " + std::to_string(m.slots()) + " slots");
        const auto p = parse_program("PSPROGRAM 1\nx\npow 0 5000\nadd 0 0\nout 2\n"); // node 1 is dead
        ProgramMachine<double> d(p);
        require(d.valid() && d.terms() == 2 && d.evaluate_one({3.0, 0.0}).finite &&
                    d.evaluate_one({3.0, 0.0}).value.re == 6.0,
                "machine evaluator: a dead node is left out of the cone");
    }
    // A constant outside binary64's range refuses the evaluator.
    {
        const auto p = parse_program("PSPROGRAM 1\nx\nc 1@-1100\nadd 0 1\nout 2\n");
        require(!ProgramMachine<double>(p).valid(), "machine evaluator: a constant below binary64's range refuses");
        require(ProgramMachine<machine::BlockDouble>(p).valid(), "machine evaluator: the range tier takes it");
    }
}

// The range tier's arithmetic (block_double.hpp): where binary64 holds the operands and the result as normal
// numbers, every operation returns binary64's bits (its scalings are exact, and an operand it drops lies below
// half an ulp of the result); beyond binary64's range it keeps the exponent; the ordering, the square root and
// the conversions from and to balls.
void block_double() {
    using B = machine::BlockDouble;
    std::mt19937_64 rng(20261017);
    std::uniform_real_distribution<double> mantissa(0.5, 1.0), exponent(-500, 500);
    std::bernoulli_distribution sign(0.5);
    const auto draw = [&] {
        const double x = std::ldexp(mantissa(rng), static_cast<int>(exponent(rng)));
        return sign(rng) ? -x : x;
    };
    const auto normal = [](double x) { return std::isfinite(x) && (x == 0 || std::abs(x) >= 0x1p-1022); };
    size_t checked = 0;
    for (int t = 0; t < 200000; ++t) {
        const double x = draw(), y = t % 7 == 0 ? x * 0.75 : draw();
        const B a(x), b(y);
        require(same_bits(static_cast<double>(a), x), "block double: a double converts both ways exactly");
        for (const auto &[exact, block] : {std::pair{x + y, a + b}, std::pair{x - y, a - b}, std::pair{x * y, a * b},
                                           std::pair{x / y, a / b}}) {
            if (!normal(exact)) continue;
            require(same_bits(static_cast<double>(block), exact), "block double: an operation returns binary64's bits");
            ++checked;
        }
        if (x > 0 && normal(std::sqrt(x)))
            require(same_bits(static_cast<double>(machine::sqrt(a)), std::sqrt(x)), "block double: the square root's bits");
        require((a < b) == (x < y) && (b < a) == (y < x), "block double: the order");
    }
    require(checked > 700000, "block double: " + std::to_string(checked) + " operations compared");
    // Beyond binary64: 2^900 squared, a sum of two such, a quotient back into range, the order across exponents.
    const B big = B(0x1p900) * B(0x1p900);
    require(std::abs(machine::logabs(big) - 1800) < 1e-12 && machine::finite(big) &&
                std::isinf(static_cast<double>(big)),
            "block double: 2^1800 beyond binary64");
    require(std::abs(machine::logabs(big + big) - 1801) < 1e-12, "block double: a sum beyond binary64");
    require(same_bits(static_cast<double>(big / (B(0x1p900) * B(0x1p899))), 2.0), "block double: back into range");
    require(B(0x1p900) < big && -big < B(-0x1p900) && B() < big && -big < B(), "block double: the order across exponents");
    require(std::abs(machine::logabs(machine::sqrt(big)) - 900) < 1e-12, "block double: the square root beyond binary64");
    // Through balls: 2^-5000 and 3 * 2^7001 both ways, exactly.
    for (const auto &[mantissa_value, power] : {std::pair{1.0, -5000L}, std::pair{3.0, 7001L}}) {
        Float in, out;
        arf_set_d(in.value(), mantissa_value);
        arf_mul_2exp_si(in.value(), in.value(), power);
        B x;
        require(machine::from_arf(x, in.value()), "block double: from a ball's midpoint");
        machine::to_arf(out.value(), x);
        require(arf_equal(in.value(), out.value()), "block double: through a ball exactly");
    }
}

// The batched loops return the bits of the one-point path (the same loops with a chunk of one), on random
// programs and on Mandelbrot k = 8 at 64 points; where finite and above the floor, within 2^-40 relative of
// the balls' midpoints or within eight times the radius of the same evaluation in 53-bit balls (the rigorous
// forward bound of a 53-bit evaluation, which a cancellation above the output's last sum exceeds while the
// floor, a lower bound by design, does not see it).
void machine_bits() {
    std::mt19937_64 rng(20261016);
    size_t compared = 0, checked = 0;
    std::vector<machine::Complex<double>> zs;
    for (int j = 0; j < 64; ++j) zs.push_back({0.3 * std::cos(0.7 * j + 0.1), 0.3 * std::sin(0.7 * j + 0.1)});
    zs[1] = {0x1p40, 0.5};
    zs[2] = {0x1p-40, -0x1p-41};
    zs[3] = {-1.5, 2.25};
    const auto check = [&](const Program &p, const std::string &label) {
        ProgramMachine<double> m(p);
        if (!m.valid()) return;
        std::vector<ProgramMachine<double>::Evaluation> batched(zs.size()), small(zs.size());
        m.evaluate(zs, batched, 64);
        m.evaluate(zs, small, 7);
        for (size_t j = 0; j < zs.size(); ++j) {
            const auto one = m.evaluate_one(zs[j]);
            require(same_bits(one, batched[j]) && same_bits(one, small[j]), label + ": the chunks of 64 and 7 return the one-point bits");
            ++checked;
            if (!one.finite) continue;
            ComplexBall point;
            acb_set_d_d(point.value(), zs[j].re, zs[j].im);
            DualBallArithmetic dual(point.value(), 256);
            const auto w = run_program(p, dual);
            if (acb_contains_zero(w.value.value()) || acb_contains_zero(w.derivative.value())) continue;
            const auto cv = midpoint_of(w.value.value()), cd = midpoint_of(w.derivative.value());
            if (std::abs(cv) <= 16 * one.noise || std::abs(cd) <= 16 * one.noise_derivative) continue; // at the floor
            DualBallArithmetic dual53(point.value(), 53);
            const auto w53 = run_program(p, dual53);
            const auto radius = [](const acb_t b) {
                return mag_get_d(arb_radref(acb_realref(b))) + mag_get_d(arb_radref(acb_imagref(b)));
            };
            const double rv = radius(w53.value.value()), rd = radius(w53.derivative.value());
            require(std::abs(standard(one.value) - cv) <= std::max(0x1p-40 * std::abs(cv), 8 * rv),
                    label + ": the value within 2^-40 of the balls' or eight 53-bit radii");
            require(std::abs(standard(one.derivative) - cd) <= std::max(0x1p-40 * std::abs(cd), 8 * rd),
                    label + ": the derivative within 2^-40 of the balls' or eight 53-bit radii");
            ++compared;
        }
    };
    for (int t = 0; t < 300; ++t) {
        const auto p = random_program(rng, 10 + t % 30, 40);
        if (!p.output) continue;
        check(p, "random program " + std::to_string(t));
    }
    check(detail_of(mandelbrot(8)), "Mandelbrot k = 8");
    require(checked >= 10000 && compared >= 2000,
            "machine bits: " + std::to_string(checked) + " points checked, " + std::to_string(compared) + " compared");
}

// Grow and reuse the packed machine buffers across cache-line and chunk boundaries.
// Each call must retain the one-point bits, including both propagated noise fields.
void machine_storage() {
    const auto p = detail_of(chebyshev(20));
    const auto exercise = [&]<class R>() {
        using Machine = ProgramMachine<R>;
        Machine reference(p), reused(p);
        std::vector<typename Machine::C> points(1025);
        std::vector<typename Machine::Evaluation> expected(points.size()), actual(points.size());
        for (size_t i = 0; i < points.size(); ++i) {
            points[i] = {R(0.125 + 0.001 * static_cast<double>(i % 31)),
                         R(-0.125 + 0.001 * static_cast<double>(i % 19))};
            expected[i] = reference.evaluate_one(points[i]);
        }
        for (size_t count : {size_t{0}, size_t{1}, size_t{7}, size_t{8}, size_t{9}, size_t{31},
                             size_t{32}, size_t{33}, size_t{255}, size_t{256}, size_t{257},
                             size_t{510}, size_t{511}, size_t{512}, size_t{513}, size_t{1022},
                             size_t{1023}, size_t{1024}, size_t{1025}, size_t{7}, size_t{512}}) {
            reused.evaluate(std::span<const typename Machine::C>(points.data(), count),
                            std::span<typename Machine::Evaluation>(actual.data(), count), std::max(count, size_t{1}));
            for (size_t i = 0; i < count; ++i) {
                const auto &a = actual[i], &b = expected[i];
                const auto same = [](R x, R y) { return same_bits(static_cast<double>(x), static_cast<double>(y)); };
                require(a.finite == b.finite && same(a.value.re, b.value.re) && same(a.value.im, b.value.im) &&
                            same(a.derivative.re, b.derivative.re) && same(a.derivative.im, b.derivative.im) &&
                            same(a.noise, b.noise) && same(a.noise_derivative, b.noise_derivative),
                        "machine storage: growing/reusing buffers changes a point's bits");
            }
            const uint64_t fields = 6 * reused.slots();
            require(reused.storage_bytes(count) >= fields * count * sizeof(R) &&
                        reused.storage_bytes(count) <= fields * (count * sizeof(R) + 2 * 64),
                    "machine storage: padding is missing from the bounded storage estimate");
        }
        require(reused.storage_bytes(0) == 0 &&
                    reused.storage_bytes(std::numeric_limits<size_t>::max()) == UINT64_MAX,
                "machine storage: zero and unrepresentable sizes");
    };
    exercise.template operator()<double>();
    exercise.template operator()<machine::BlockDouble>();

}

// The evaluator's zero order (program_machine.hpp): with v > 0 it stands for q = p / x^v, scaled by z^v. Its Newton
// ratio and Cauchy term against q's own, from q's exact coefficients (p's expansion shifted down by v) in 256-bit
// balls, within the evaluation's noise, wherever it is above the floor; the chunks return the one-point bits; the
// range tier carries the order; a point at 0 is not finite. Then the counts and the stage on a planted product with
// a double root at zero: every circle read in binary64, q's three roots converged.
void deflated_evaluator() {
    struct Case {
        std::string label;
        Program program;
        uint64_t order;
    };
    std::vector<Case> cases;
    cases.push_back({"x^4 - x^2", parse_program("PSPROGRAM 1\nx\nc 1\nmul 1 0\nmul 2 0\nc -1\nadd 3 4\nmul 5 0\nmul 6 0\nout 7\n"), 2});
    cases.push_back({"Chebyshev T_25", detail_of(chebyshev(25)), 1});
    cases.push_back({"a sum added and taken away (seed 13957833835543924094)",
                     parse_program("PSPROGRAM 1\nx\nc 1 0\nsub 0 1\nc 1 0\nmul 3 0\nc -1/2251799813685248 -1/4503599627370496\n"
                                   "add 4 5\nmul 6 0\nc 0 0\nadd 7 8\nmul 2 9\nc -119/111 0\nmul 11 10\nc -1/6 1/9\nmul 13 0\n"
                                   "c -2 0\nadd 14 15\nadd 12 16\nsub 17 16\nout 18\n"),
                     1});
    std::vector<machine::Complex<double>> zs;
    for (int j = 0; j < 16; ++j) zs.push_back({0.3 * std::cos(0.7 * j + 0.1), 0.3 * std::sin(0.7 * j + 0.1)});
    zs.push_back({1.7, 0.4});
    zs.push_back({-0.6, -0.9});
    zs.push_back({3.0, 1.0});
    zs.push_back({0x1p-20, 0x1p-21});
    size_t compared = 0;
    for (const auto &c : cases) {
        require(prepare_program(c.program, cps_options{}).zero_order == c.order, c.label + ": the order at zero");
        ProgramMachine<double> m(c.program, c.order);
        require(m.valid(), c.label + ": the deflated evaluator");
        ProgramMachine<machine::BlockDouble> scaled(m);
        ExpansionArithmetic exact(unlimited);
        const auto e = run_program(c.program, exact);
        ComplexPolynomial q; // p's coefficients from v on
        for (slong j = static_cast<slong>(c.order); j <= e.degree(); ++j) {
            Rational re, im;
            fmpq_poly_get_coeff_fmpq(re.value(), e.real(), j);
            fmpq_poly_get_coeff_fmpq(im.value(), e.imag(), j);
            ComplexBall coefficient;
            arb_set_fmpq(acb_realref(coefficient.value()), re.value(), 256);
            arb_set_fmpq(acb_imagref(coefficient.value()), im.value(), 256);
            acb_poly_set_coeff_acb(q.value(), j - static_cast<slong>(c.order), coefficient.value());
        }
        std::vector<ProgramMachine<double>::Evaluation> batched(zs.size()), small(zs.size());
        m.evaluate(zs, batched, 64);
        m.evaluate(zs, small, 7);
        for (size_t j = 0; j < zs.size(); ++j) {
            const auto one = m.evaluate_one(zs[j]);
            require(same_bits(one, batched[j]) && same_bits(one, small[j]), c.label + ": the chunks return the one-point bits");
            require(one.finite, c.label + ": finite at a point away from zero");
            const double value = machine::norm(one.value), derivative = machine::norm(one.derivative);
            if (value <= 16 * one.noise || derivative <= 16 * one.noise_derivative) continue; // at the floor
            ComplexBall z, qv, qd;
            acb_set_d_d(z.value(), zs[j].re, zs[j].im);
            acb_poly_evaluate2(qv.value(), qd.value(), q.value(), z.value(), 256);
            const auto cv = midpoint_of(qv.value()), cd = midpoint_of(qd.value());
            const std::complex<double> point = standard(zs[j]);
            // The evaluation's relative error is about its noise over its size, the floor being a lower bound of it.
            const double tolerance = 0x1p-30 + 64 * (one.noise / value + one.noise_derivative / derivative);
            const auto ratio = standard(one.value) / standard(one.derivative), reference = cv / cd;
            require(std::abs(ratio - reference) <= tolerance * std::abs(reference), c.label + ": the Newton ratio is q / q'");
            const auto term = point * standard(one.derivative) / standard(one.value), cauchy = point * cd / cv;
            require(std::abs(term - cauchy) <= tolerance * std::max(1.0, std::abs(cauchy)), c.label + ": the Cauchy term is z q' / q");
            const auto s = scaled.evaluate_one({machine::BlockDouble(zs[j].re), machine::BlockDouble(zs[j].im)});
            const machine::Complex<machine::BlockDouble> scaled_ratio = s.value / s.derivative;
            const std::complex<double> r{static_cast<double>(scaled_ratio.re), static_cast<double>(scaled_ratio.im)};
            require(s.finite && std::abs(r - ratio) <= 0x1p-40 * std::abs(ratio), c.label + ": the range tier carries the order");
            ++compared;
        }
        require(!m.evaluate_one({0.0, 0.0}).finite, c.label + ": a point at 0 is not finite");
    }
    require(compared >= 50, "deflated evaluator: " + std::to_string(compared) + " points compared");
    // x^2 (x - 1) (x - 2i) (x + 3/2 + i): the counts of q's three roots in binary64, and the stage converges them.
    const cps_options options;
    certpolysolve::detail::CandidateStatistics statistics{};
    const Program p = from_roots({{"0", "0"}, {"0", "0"}, {"1", "0"}, {"0", "2"}, {"-3/2", "-1"}});
    const auto prepared = prepare_program(p, options);
    require(prepared.ready() && prepared.zero_order == 2 && prepared.degree == 5, "the planted product's double root at zero");
    ProgramSource source(p, 3, prepared.leading, 2);
    ProgramMachine<double> counting(p, 2);
    ProgramSeedReport report;
    RootVector seeds(3);
    require(program_seeds(seeds.data(), source, 3, prepared.root_bound.value(), 85, options, statistics, nullptr, &report,
                          &counting),
            "deflated evaluator: the seeds placed");
    require(report.circles > 0 && report.binary64_circles == report.circles,
            "deflated evaluator: every circle counted in binary64 (" + std::to_string(report.binary64_circles) + " of " +
                std::to_string(report.circles) + ")");
    ProgramMachineStage stage(p, seeds.data(), 3, 1024, UINT64_MAX, nullptr, 2);
    require(stage.valid(), "deflated evaluator: the stage from the seeds");
    const auto step = stage.advance(200, options, statistics);
    require(step.finished && step.converged == 3 && step.sweeps < 100,
            "deflated evaluator: q's roots converged in binary64 in " + std::to_string(step.sweeps) + " sweeps");
    RootVector out(3);
    stage.export_points(out.data());
    std::vector<std::complex<double>> planted{{1, 0}, {0, 2}, {-1.5, -1}};
    for (slong i = 0; i < 3; ++i) {
        const auto z = midpoint_of(out.data() + i);
        const auto nearest = std::min_element(planted.begin(), planted.end(),
                                              [&](auto a, auto b) { return std::abs(a - z) < std::abs(b - z); });
        require(std::abs(*nearest - z) <= 0x1p-40 * std::abs(*nearest), "deflated evaluator: a point within 2^-40 of a root of q");
        planted.erase(nearest);
    }
    std::cout << "deflated evaluator: " << cases.size() << " programs with a root at zero, " << compared
              << " points against q's coefficients; " << report.circles << " circles in binary64, q's roots in "
              << step.sweeps << " sweeps\n";
}

// The floor (program_machine.hpp) in the stage's terms: a point whose value is at the floor has converged when
// its correction is small, and is noise when it is not. At every root of Laguerre L_64, Chebyshev T_128 and
// Mandelbrot k = 8 the evaluation reads as converged (L_64's cancellation at its large roots is matched by its
// derivative, so the Newton correction is accurate; the reset at sums keeps the floor at a few roundings where
// a propagated bound would grow like the balls' radii); a difference of two large expressions near its root is
// at the floor (the value is roundoff of those, the distance to the root below it); a double root is at the
// derivative's floor.
void noise_floor() {
    struct Read {
        bool floor;            // |p| below 16 times its noise
        bool derivative_floor; // |p'| below 16 times its noise
        double step;           // the Newton correction |p / p'| relative to |z|
    };
    const auto read_at = [](ProgramMachine<double> &m, double x, double y) {
        const auto e = m.evaluate_one({x, y});
        require(e.finite, "noise floor: a finite evaluation");
        const double p = std::abs(standard(e.value)), d = std::abs(standard(e.derivative));
        const double z = std::max(std::hypot(x, y), 1e-300);
        return Read{p < 16 * e.noise, d < 16 * e.noise_derivative, p / std::max(d, 1e-300) / z};
    };
    const auto read = [&](ProgramMachine<double> &m, const Box &box) {
        return read_at(m, (fmpq_get_d(box.lo.value()) + fmpq_get_d(box.hi.value())) / 2,
                       (fmpq_get_d(box.ilo.value()) + fmpq_get_d(box.ihi.value())) / 2);
    };
    for (const auto &[program, label] : {std::pair{laguerre(64), "Laguerre L_64"}, std::pair{chebyshev(128), "Chebyshev T_128"},
                                         std::pair{mandelbrot(8), "Mandelbrot k = 8"}}) {
        const auto answer = solve_program(program, certpolysolve::Options{});
        require(answer.complete() && answer.goal_met(), std::string("noise floor: ") + label + " solved");
        ProgramMachine<double> m(detail_of(program));
        size_t converged = 0;
        const auto boxes = boxes_of(answer);
        for (const auto &box : boxes) {
            const auto r = read(m, box);
            converged += !r.derivative_floor && r.step < 0x1p-45;
        }
        require(converged == boxes.size(), std::string("noise floor: ") + label + ": " + std::to_string(converged) + " of " +
                                               std::to_string(boxes.size()) + " roots read as converged");
    }
    {
        // (x^20 + (x - 2)) - (x^10)^2 is x - 2, computed as the difference of two values near 2^20 reached by
        // different roundings (a power of 20, a power of 10 squared): at 2 + 2^-40 the computed value is their
        // roundoff, far above the distance to the root, while the derivative is 1 (its cancellation is matched),
        // so the Newton correction is roundoff, not the distance: noise, at the floor.
        const auto p = parse_program("PSPROGRAM 1\nx\npow 0 20\nc 2\nsub 0 2\nadd 1 3\npow 0 10\nmul 5 5\nsub 4 6\nout 7\n");
        ProgramMachine<double> m(p);
        // Within the floor's band the roundings are deterministic and, for a program this short, exact: the two
        // paths agree bit for bit and the computed value is 0 (a zero correction: doubles stop there, which the
        // stage reads as converged, 2^-40 from the root). Roundoff that moves a point is the long program's
        // case; the double root below is the floor's other reading.
        for (const double factor : {1.0, 1.3, 1.7, 2.1, 2.9})
            require(read_at(m, 2 + factor * 0x1p-40, 0).floor && !read_at(m, 2 + factor * 0x1p-40, 0).derivative_floor,
                    "noise floor: a cancelled difference near its root is at the floor");
        require(!read_at(m, 2.5, 0).floor, "noise floor: the same program away from its root is readable");
    }
    {
        // A double root: (x - 1)^2 (x + 2) at 1 + 2^-50 is at the derivative's floor (at 1 + 2^-40 the
        // derivative, 6 x 2^-40, is still readable and Newton's correction, 2^-41, still moves the point: the
        // linear approach to a double root). The simple root at -2 + 2^-50 reads as converged.
        const auto p = parse_program("PSPROGRAM 1\nx\nc 1\nsub 0 1\nmul 2 2\nc 2\nadd 0 4\nmul 3 5\nout 6\n");
        ProgramMachine<double> m(p);
        require(read_at(m, 1 + 0x1p-50, 0).derivative_floor, "noise floor: a double root is at the derivative's floor");
        require(!read_at(m, 1 + 0x1p-40, 0).derivative_floor && read_at(m, 1 + 0x1p-40, 0).step > 0x1p-45,
                "noise floor: farther from the double root the derivative is readable and the correction moves");
        require(!read_at(m, -2 + 0x1p-50, 0).derivative_floor && read_at(m, -2 + 0x1p-50, 0).step < 0x1p-45,
                "noise floor: the simple root of the same program reads as converged");
    }
}

// The binary64 stage (program_candidates.hpp) on the planted roots from the counts' placement: every point
// converged within 2^-40 of a planted root in a few sweeps; cancellation before a sweep; a square factor's
// pair at the floor and its third root converged; a seed out of binary64's range held while the others
// converge; a storage refusal.
void machine_stage() {
    const cps_options options;
    certpolysolve::detail::CandidateStatistics statistics{};
    struct Always {
        static int poll(void *) { return 1; }
    };
    const auto place = [&](const Program &q, slong degree, const ProgramPreparation &prepared, RootVector &seeds) {
        ProgramSource source(q, degree, prepared.leading);
        ProgramSeedReport report;
        require(program_seeds(seeds.data(), source, degree, prepared.root_bound.value(), 85, options, statistics, nullptr, &report),
                "the stage test's seeds placed");
    };
    {
        // One Jacobi sweep of p(z) = z^2 - 1 from (2, -3) gives (19/17, -13/11):
        // both interaction sums use the original pair. Updating the first point
        // before forming the second sum would instead give -73/71 there.
        const Program p = from_roots({{"1", "0"}, {"-1", "0"}});
        for (const bool reversed : {false, true}) {
            RootVector seeds(2), out(2);
            acb_set_si(seeds[0], reversed ? -3 : 2);
            acb_set_si(seeds[1], reversed ? 2 : -3);
            ProgramMachineStage stage(p, seeds.data(), 2, 2, UINT64_MAX);
            require(stage.valid(), "the two-point Jacobi stage");
            const auto step = stage.advance(1, options, statistics);
            require(step.sweeps == 1 && !step.cancelled, "exactly one Jacobi sweep");
            stage.export_points(out.data());
            // Comparison slack for evaluation, complex quotients and the update;
            // the sequential-update error is about 0.15, far above roundoff.
            constexpr double slack = 32 * std::numeric_limits<double>::epsilon();
            for (slong i = 0; i < 2; ++i) {
                const double expected = ((i == 0) != reversed) ? 19.0 / 17 : -13.0 / 11;
                const auto actual = midpoint_of(out[i]);
                require(actual.imag() == 0 && std::abs(actual.real() - expected) <= slack * std::abs(expected),
                        "every correction reads the sweep's initial points, in either seed order");
            }
        }
    }
    {
        const std::vector<std::pair<const char *, const char *>> roots{{"1", "0"}, {"0", "2"}, {"-3/2", "-1"}, {"7", "0"}};
        const Program p = from_roots(roots);
        const auto prepared = prepare_program(p, options);
        RootVector seeds(4);
        place(p, 4, prepared, seeds);
        ProgramMachineStage stage(p, seeds.data(), 4, 1024, UINT64_MAX);
        require(stage.valid(), "the stage from the seeds");
        const uint64_t before = statistics.candidate_sweeps;
        const auto step = stage.advance(200, options, statistics);
        require(step.finished && step.converged == 4 && !step.stalled && step.sweeps > 0 && step.sweeps < 100 &&
                    statistics.candidate_sweeps == before + step.sweeps,
                "every planted root converged in binary64 in " + std::to_string(step.sweeps) + " sweeps");
        RootVector out(4);
        stage.export_points(out.data());
        std::vector<std::complex<double>> planted{{1, 0}, {0, 2}, {-1.5, -1}, {7, 0}};
        for (slong i = 0; i < 4; ++i) {
            const auto z = midpoint_of(out.data() + i);
            const auto nearest = std::min_element(planted.begin(), planted.end(), [&](auto a, auto b) {
                return std::abs(a - z) < std::abs(b - z);
            });
            require(std::abs(*nearest - z) <= 0x1p-40 * std::abs(*nearest), "a converged point within 2^-40 of a planted root");
            planted.erase(nearest);
        }
        ProgramMachineStage fresh(p, seeds.data(), 4, 1024, UINT64_MAX);
        cps_options cancelling = options;
        cancelling.cancel = Always::poll;
        const auto cancelled = fresh.advance(10, cancelling, statistics);
        require(cancelled.cancelled && cancelled.sweeps == 0 && cancelled.active == 4, "a cancelled stage took no sweep");
        require(!ProgramMachineStage(p, seeds.data(), 4, 1024, 1).valid(), "a stage over the storage allowance is refused");
    }
    {
        const Program q = from_roots({{"1", "0"}, {"1", "0"}, {"-2", "0"}});
        const auto prepared = prepare_program(q, options);
        RootVector seeds(3);
        place(q, 3, prepared, seeds);
        ProgramMachineStage stage(q, seeds.data(), 3, 1024, UINT64_MAX);
        require(stage.valid(), "the square factor's stage");
        // The pair at the double root converges in doubles to the roots of the polynomial as computed, a pair
        // split by about the square root of the noise (Aberth's corrections there go to zero), or stops at the
        // floor: either way it finishes, the pair within 2^-20 of 1 and the third root within 2^-40 of -2, and
        // the ball phase and the refinement take the cluster from there (square_factor).
        const auto step = stage.advance(400, options, statistics);
        require(step.finished && step.converged + step.floor == 3,
                "the square factor: " + std::to_string(step.converged) + " converged, " + std::to_string(step.floor) +
                    " at the floor in " + std::to_string(step.sweeps) + " sweeps");
        RootVector out(3);
        stage.export_points(out.data());
        size_t near_one = 0, near_minus_two = 0;
        for (slong i = 0; i < 3; ++i) {
            const auto z = midpoint_of(out.data() + i);
            near_one += std::abs(z - std::complex<double>(1, 0)) <= 0x1p-20;
            near_minus_two += std::abs(z - std::complex<double>(-2, 0)) <= 0x1p-40 * 2;
        }
        require(near_one == 2 && near_minus_two == 1, "the square factor's pair near 1 and its third root at -2");
    }
    {
        const Program m = detail_of(mandelbrot(6));
        const auto prepared = prepare_program(m, options);
        RootVector seeds(63);
        place(m, 63, prepared, seeds);
        acb_set_d_d(seeds.data(), 1e6, 0); // p_6 there is about 10^378: out of binary64's range
        ProgramMachineStage stage(m, seeds.data(), 63, 1024, UINT64_MAX);
        require(stage.valid(), "Mandelbrot k = 6's stage with a far seed");
        const auto step = stage.advance(800, options, statistics);
        require(step.finished && step.held == 0 && step.scaled_evaluations > 0 && step.converged + step.floor == 63,
                "a seed out of binary64's range is read in the scaled arithmetic and converges with the others: " +
                    std::to_string(step.held) + " held, " + std::to_string(step.converged) + " converged, " +
                    std::to_string(step.floor) + " at the floor, " + std::to_string(step.scaled_evaluations) +
                    " scaled evaluations");
    }
}

// The second version's certificate, phase 2a (program_certificate.hpp): Newton disks from the program's
// evaluations in a disjoint cover. Separated planted roots certify, one region each; two points at one root, and a
// double root, never certify (the pass fails and the refinement would answer); Chebyshev T_128 from its analytic
// roots and Mandelbrot k = 8 from a solve's midpoints certify (the sort-and-sweep path, n >= 128); cancellation
// leaves the round to do again.
void separated_roots() {
    const cps_options options; // the accuracy goal at 53 bits: count and multiplicity 1, each root to 2^-53
    const auto run_with = [](MonicSource &source, const std::vector<std::complex<double>> &start, slong precision,
                             const cps_options &goal) {
        RootVector points(static_cast<slong>(start.size()));
        for (size_t i = 0; i < start.size(); ++i)
            acb_set_d_d(points.data() + static_cast<slong>(i), start[i].real(), start[i].imag());
        auto pass = std::make_unique<SeparatedRootsPass>(points.data(), start.size(), precision);
        certpolysolve::detail::CandidateStatistics statistics{};
        for (int call = 0; call < 200 && pass->state() == SeparatedRootsPass::State::Running; ++call)
            pass->advance(source, goal, statistics);
        return pass;
    };
    const auto run = [&](MonicSource &source, const std::vector<std::complex<double>> &start, slong precision) {
        return run_with(source, start, precision, options);
    };
    const auto contains = [](const Root &region, std::complex<double> z) {
        ComplexBall point;
        acb_set_d_d(point.value(), z.real(), z.imag());
        return acb_contains(region.box.value(), point.value());
    };
    {
        const std::vector<std::complex<double>> planted{{1, 0}, {0, 2}, {-1.5, -1}, {7, 0}};
        const Program p = from_roots({{"1", "0"}, {"0", "2"}, {"-3/2", "-1"}, {"7", "0"}});
        const auto prepared = prepare_program(p, options);
        ProgramSource source(p, 4, prepared.leading);
        std::vector<std::complex<double>> start;
        for (const auto &z : planted) start.push_back(z + std::complex<double>(1e-3, -7e-4));
        const auto pass = run(source, start, 85);
        require(pass->state() == SeparatedRootsPass::State::Certified,
                "separated roots: the planted roots certify (" + pass->failure() + ")");
        for (const auto &z : planted) {
            size_t holding = 0;
            for (const auto &region : pass->cover()) holding += contains(region, z);
            require(holding == 1, "separated roots: each planted root in exactly one region");
        }
        for (const auto &region : pass->cover())
            require(region.count == 1 && region.multiplicity == 1 &&
                        region.certificate == "original-source Newton inclusion disk in a complete disjoint cover",
                    "separated roots: count and multiplicity 1, the disk certificate");
        // Two points near the root 1 and none near -3/2 - i: Aberth's repulsion may part them and certify, or the
        // pass fails; a certificate it gives must hold every planted root in exactly one region.
        const auto duplicate = run(source, {{1 + 1e-6, 0}, {1 - 1e-6, 0}, {1e-3, 2}, {7 - 1e-3, 0}}, 85);
        if (duplicate->state() == SeparatedRootsPass::State::Certified)
            for (const auto &z : planted) {
                size_t holding = 0;
                for (const auto &region : duplicate->cover()) holding += contains(region, z);
                require(holding == 1, "separated roots: a certificate from two points at one root holds each root once");
            }
        else
            require(duplicate->state() == SeparatedRootsPass::State::Failed, "separated roots: two points at one root end");
        std::cout << "separated roots: two points at one root "
                  << (duplicate->state() == SeparatedRootsPass::State::Certified ? "parted and certified" : "failed")
                  << " after " << duplicate->rounds() << " rounds\n";
        // Two points exactly on the root 1 (their ball difference contains 0, so the sums skip it and Aberth cannot part
        // them) and none near -3/2 - i: the pass re-homes at least one of them (nudged off, then iterated with the
        // occupant in its sum) and certifies a cover that holds each planted root once.
        const auto coincident = run(source, {{1, 0}, {1, 0}, {1e-3, 2}, {7 - 1e-3, 0}}, 85);
        require(coincident->state() == SeparatedRootsPass::State::Certified && coincident->rehomed() >= 1,
                "separated roots: two points exactly on one root are re-homed and certified (" + coincident->failure() +
                    ", " + std::to_string(coincident->rehomed()) + " re-homed)");
        for (const auto &z : planted) {
            size_t holding = 0;
            for (const auto &region : coincident->cover()) holding += contains(region, z);
            require(holding == 1, "separated roots: the re-homed cover holds each planted root once");
        }
        std::cout << "separated roots: two points exactly on one root re-homed (" << coincident->rehomed()
                  << ") and certified after " << coincident->rounds() << " rounds\n";
        // Cancellation: the round is left to do again.
        cps_options cancelling = options;
        struct Always {
            static int poll(void *) { return 1; }
        };
        cancelling.cancel = Always::poll;
        RootVector points(4);
        for (slong i = 0; i < 4; ++i) acb_set_d_d(points.data() + i, start[static_cast<size_t>(i)].real(), start[static_cast<size_t>(i)].imag());
        SeparatedRootsPass cancelled(points.data(), 4, 85);
        certpolysolve::detail::CandidateStatistics statistics{};
        require(cancelled.advance(source, cancelling, statistics) == SeparatedRootsPass::State::Running &&
                    cancelled.evaluations() == 0 && cancelled.rounds() == 0,
                "separated roots: a cancelled round does nothing");
    }
    {
        // A double root: (x - 1)^2 (x + 2). Under the accuracy goal, which needs every root alone, never certified.
        const Program q = from_roots({{"1", "0"}, {"1", "0"}, {"-2", "0"}});
        const auto prepared = prepare_program(q, options);
        ProgramSource source(q, 3, prepared.leading);
        const std::vector<std::complex<double>> start{{1 + 1e-3, 1e-4}, {1 - 1e-3, -1e-4}, {-2 + 1e-3, 0}};
        const auto pass = run(source, start, 85);
        require(pass->state() == SeparatedRootsPass::State::Failed && pass->counts() == 0,
                "separated roots: a double root must not certify under the accuracy goal");
        // Under the approximation goal (phase 2b): the two points at the double root keep overlapping (Aberth's
        // rate there is linear), so the pass's cluster step counts them on a circle, shrinks it to the goal and
        // certifies a counted cluster of 2 with no multiplicity claim beside a Newton disk about -2.
        cps_options approximation;
        approximation.goal = CPS_ROOT_APPROXIMATION;
        const auto clustered = run_with(source, start, 85, approximation);
        require(clustered->state() == SeparatedRootsPass::State::Certified && clustered->clusters() == 1 &&
                    clustered->cover().size() == 2,
                "separated roots: a double root is a counted cluster under the approximation goal (" +
                    clustered->failure() + ", " + std::to_string(clustered->counts()) + " counts)");
        uint64_t total = 0;
        for (const auto &region : clustered->cover()) {
            const bool double_root = contains(region, {1, 0}), simple = contains(region, {-2, 0});
            require(double_root != simple, "separated roots: each region holds one of the two roots");
            if (double_root)
                require(region.count == 2 && region.multiplicity == 0 &&
                            region.certificate == "original-source argument-principle count in a complete disjoint cover",
                        "separated roots: the double root's region counts 2, with no multiplicity claim");
            else
                require(region.count == 1 && region.multiplicity == 1, "separated roots: -2's region is one simple root");
            require(meets_goal(region, approximation, clustered->precision()), "separated roots: the cluster meets the goal");
            total += region.count;
        }
        require(total == 3, "separated roots: the clustered cover counts the degree");
        // The step runs once every region meets the goal with overlaps left (the factored program reads the pair's
        // values exactly, so the points converge far below the precision and never stall): well before the rounds
        // allowed run out.
        require(clustered->rounds() < heuristics::program_disk_rounds / 2,
                "separated roots: the cluster step runs before the rounds run out (" +
                    std::to_string(clustered->rounds()) + " rounds)");
        std::cout << "separated roots: a double root counted as a cluster of 2 after " << clustered->rounds() << " rounds, "
                  << clustered->counts() << " counts (" << clustered->count_arcs() << " arcs), " << clustered->raises()
                  << " raises, " << clustered->precision() << " bits\n";
    }
    {
        // Chebyshev T_128 from its analytic roots cos((2k - 1) pi / 256), perturbed.
        const auto program = detail_of(chebyshev(128));
        const auto prepared = prepare_program(program, options);
        ProgramSource source(program, 128, prepared.leading);
        std::vector<std::complex<double>> roots, start;
        for (int k = 1; k <= 128; ++k) {
            roots.push_back({std::cos((2 * k - 1) * std::numbers::pi / 256), 0});
            start.push_back(roots.back() * (1 + 1e-9) + std::complex<double>(0, 1e-10));
        }
        const auto pass = run(source, start, 85);
        require(pass->state() == SeparatedRootsPass::State::Certified,
                "separated roots: Chebyshev T_128 certifies (" + pass->failure() + ", " + std::to_string(pass->rounds()) +
                    " rounds, " + std::to_string(pass->precision()) + " bits)");
        for (size_t i = 0; i < roots.size(); ++i) {
            const auto mid = midpoint_of(pass->cover()[i].box.value());
            require(std::abs(mid - roots[i]) < 1e-12, "separated roots: T_128's region at its root");
        }
        std::cout << "separated roots: Chebyshev T_128 certified in " << pass->rounds() << " rounds, " << pass->evaluations()
                  << " evaluations, " << pass->raises() << " raises, " << pass->precision() << " bits\n";
    }
    {
        // Mandelbrot k = 8 from a solve's midpoints: the pass certifies them without the refinement's cover, and
        // each certified region meets exactly one region of the solve.
        const auto program = mandelbrot(8);
        const auto answer = solve_program(program, certpolysolve::Options{});
        const auto boxes = boxes_of(answer);
        std::vector<std::complex<double>> start;
        for (const auto &box : boxes)
            start.push_back({(fmpq_get_d(box.lo.value()) + fmpq_get_d(box.hi.value())) / 2,
                             (fmpq_get_d(box.ilo.value()) + fmpq_get_d(box.ihi.value())) / 2});
        const auto prepared = prepare_program(detail_of(program), options);
        ProgramSource source(detail_of(program), 255, prepared.leading);
        const auto pass = run(source, start, 85);
        require(pass->state() == SeparatedRootsPass::State::Certified,
                "separated roots: Mandelbrot k = 8 certifies (" + pass->failure() + ")");
        for (size_t i = 0; i < start.size(); ++i)
            require(contains(pass->cover()[i], start[i]) || std::abs(midpoint_of(pass->cover()[i].box.value()) - start[i]) < 1e-12,
                    "separated roots: Mandelbrot's certified region at the solve's root");
        std::cout << "separated roots: Mandelbrot k = 8 certified in " << pass->rounds() << " rounds, " << pass->evaluations()
                  << " evaluations, " << pass->raises() << " raises\n";
    }
}

// Phase 2b's root counts on circles (certified_root_count): planted roots inside and outside circles, a double root
// counted twice, a circle through a root left without a count, cancellation; and the root bound by doubling, from
// counts alone, on Mandelbrot and on a three-term recurrence.
void winding_counts() {
    const cps_options options;
    certpolysolve::detail::CandidateStatistics statistics{};
    const auto count = [&](MonicSource &source, std::complex<double> centre, double radius, const cps_options &o) {
        ComplexBall c;
        acb_set_d_d(c.value(), centre.real(), centre.imag());
        Magnitude r;
        mag_set_d(r.value(), radius);
        return certified_root_count(source, c.value(), r.value(), 85, 4096, o, statistics);
    };
    const Program p = from_roots({{"1", "0"}, {"0", "2"}, {"-3/2", "-1"}, {"7", "0"}});
    ProgramSource source(p, 4, prepare_program(p, options).leading);
    struct Case {
        std::complex<double> centre;
        double radius;
        uint64_t roots;
    };
    uint64_t arcs = 0;
    for (const auto &[centre, radius, roots] : std::vector<Case>{{{0, 0}, 0.5, 0},
                                                              {{1, 0}, 0.5, 1},
                                                              {{0, 2}, 1e-3, 1},
                                                              {{-1.5, -1}, 0.25, 1},
                                                              {{4, 0}, 3.5, 2},
                                                              {{0, 0}, 3, 3},
                                                              {{0, 0}, 10, 4}}) {
        const auto counted = count(source, centre, radius, options);
        require(counted.roots && *counted.roots == roots,
                "winding counts: the circle about (" + std::to_string(centre.real()) + ", " +
                    std::to_string(centre.imag()) + ") of radius " + std::to_string(radius) + " holds " +
                    std::to_string(roots) + " roots, counted " +
                    (counted.roots ? std::to_string(*counted.roots) : std::string("none")));
        arcs += counted.arcs;
    }
    // A circle exactly through the root 1 (radius 1 set exactly: mag_set_d rounds up, and the circle of that
    // radius holds the root 1, which a count may then rightly include): no count, and not for cancellation.
    ComplexBall origin;
    Magnitude one;
    mag_one(one.value());
    const auto through = certified_root_count(source, origin.value(), one.value(), 85, 4096, options, statistics);
    require(!through.roots && !through.cancelled,
            "winding counts: a circle through a root gives no count (counted " +
                (through.roots ? std::to_string(*through.roots) : std::string("none")) + ")");
    // Just outside it, a 2^-29 farther: the root 1 inside, counted with the others within the circle.
    const auto outside = count(source, {0, 0}, 1, options);
    require(!outside.roots || *outside.roots == 1, "winding counts: a circle just outside the root 1 holds it alone");
    // Multiplicity: the double root of (x - 1)^2 (x + 2) counts twice.
    const Program q = from_roots({{"1", "0"}, {"1", "0"}, {"-2", "0"}});
    ProgramSource double_root(q, 3, prepare_program(q, options).leading);
    const auto twice = count(double_root, {1, 0}, 0.5, options);
    require(twice.roots && *twice.roots == 2, "winding counts: a double root counts twice");
    // Cancellation: no count, said so.
    cps_options cancelling = options;
    struct Always {
        static int poll(void *) { return 1; }
    };
    cancelling.cancel = Always::poll;
    const auto cancelled = count(source, {0, 0}, 3, cancelling);
    require(!cancelled.roots && cancelled.cancelled, "winding counts: a cancelled count gives no count");
    std::cout << "winding counts: 7 circles certified in " << arcs << " arcs; a circle through a root refused after "
              << through.arcs << " arcs; one just outside it "
              << (outside.roots ? "counted 1 in " + std::to_string(outside.arcs) + " arcs" : std::string("refused"))
              << '\n';
}

// The twins: the program's answer against the expanded polynomial's, at the same goal and precision.
void twins() {
    struct Case {
        std::string label;
        certpolysolve::Program program;
    };
    std::vector<Case> cases;
    // k <= 6 keeps the test to about a second (test_program --mandelbrot K prints the ladder beyond it, each rung
    // against its twin).
    for (int k = 2; k <= 6; ++k) cases.push_back({"Mandelbrot k=" + std::to_string(k), mandelbrot(k)});
    for (int n : {8, 32, 64}) cases.push_back({"Chebyshev T_" + std::to_string(n), chebyshev(n)});
    for (int n : {8, 24}) cases.push_back({"Hermite H_" + std::to_string(n), hermite(n)});
    for (int n : {8, 20}) cases.push_back({"Laguerre L_" + std::to_string(n), laguerre(n)});
    for (int n : {8, 32}) cases.push_back({"tridiagonal " + std::to_string(n), tridiagonal(n)});
    for (const auto &[label, program] : cases) {
        const auto twin = twin_of(detail_of(program));
        const uint64_t degree = twin.coefficients.size() - 1;
        // The default goal of a program is CPS_ROOT_APPROXIMATION; the twin is asked for it explicitly.
        certpolysolve::Options twin_options;
        twin_options.goal = certpolysolve::Goal::approximation();
        const auto answer = solve_program(program, certpolysolve::Options{});
        require(cps_result_goal(answer.get()) == CPS_ROOT_APPROXIMATION,
                label + ": a program's default goal is approximation");
        require(answer.message().find("Aberth mode") != std::string::npos,
                label + ": Aberth mode is the default: " + answer.message());
        const auto twin_answer = solve_twin(twin, twin_options);
        cross_check(answer, twin_answer, degree, 53, label + " (approximation, Auto)");
        // The start at the roots' radii (program_seeds) keeps the program's sweeps near the twin's: from the
        // bound's circle the preliminary phase crawled about (n / 2) ln(R / r) sweeps (Mandelbrot k = 6: about
        // 150 against the twin's 38; T_64: 94), and today's are 27 and 33 (16 October 2026).
        // A twin of degree 3 or below is answered by Auto's closed forms, without sweeps: no comparison then.
        const uint64_t program_sweeps = answer.statistics().sweeps;
        const uint64_t twin_sweeps = twin_answer.statistics().sweeps;
        const bool twin_direct = twin_answer.message().find("auto-route=direct") != std::string::npos;
        require(twin_direct || program_sweeps <= 2 * twin_sweeps + 8,
                label + ": the program's sweeps (" + std::to_string(program_sweeps) + ") within twice the twin's (" +
                    std::to_string(twin_sweeps) + ") and eight");
        // The separated-roots pass (program_certificate.hpp) certifies Mandelbrot's and Chebyshev's roots by
        // Newton disks, with no regeneration; the others may fall back to the refinement (a root at 0, which no
        // disk about a nonzero point meets a relative goal for, or a cluster).
        const bool by_disks = answer.message().find("Newton disks: a complete disjoint cover") != std::string::npos;
        if (label.rfind("Mandelbrot", 0) == 0 || label.rfind("Chebyshev", 0) == 0)
            require(by_disks && answer.engine_counters().secular_regenerations == 0,
                    label + ": certified by Newton disks without a regeneration: " + answer.message());
        std::cout << "twins: " << label << " degree " << degree << " approximation: sweeps " << program_sweeps
                  << " (twin " << twin_sweeps << "), regenerations " << answer.engine_counters().secular_regenerations
                  << ", work bits " << answer.work_bits() << (by_disks ? ", Newton disks" : ", the refinement") << '\n';
    }
    // The accuracy goal (every root in its own region, multiplicity 1), explicit classical (the Aberth mode) and
    // secular mode, and a tighter precision, on a few of them.
    {
        certpolysolve::Options accuracy;
        accuracy.goal = certpolysolve::Goal::accuracy();
        std::vector<Case> accurate;
        accurate.push_back({"Mandelbrot k=6", mandelbrot(6)});
        accurate.push_back({"Hermite H_16", hermite(16)});
        for (const auto &[label, program] : accurate) {
            const auto twin = twin_of(detail_of(program));
            cross_check(solve_program(program, accuracy), solve_twin(twin, accuracy), twin.coefficients.size() - 1, 53,
                        label + " (accuracy)");
        }
        certpolysolve::Options classical;
        classical.method = certpolysolve::Method::classical;
        const auto c = chebyshev(24);
        const auto classical_answer = solve_program(c, classical);
        require(classical_answer.message().find("Aberth mode") != std::string::npos, "classical is the Aberth mode");
        cross_check(classical_answer,
                    solve_twin(twin_of(detail_of(c)), certpolysolve::Options{.goal = certpolysolve::Goal::approximation()}), 24,
                    53, "Chebyshev T_24 (classical)");
        certpolysolve::Options secular;
        secular.method = certpolysolve::Method::secular;
        const auto secular_answer = solve_program(c, secular);
        require(secular_answer.message().find("secular mode") != std::string::npos, "secular is the secular mode");
        cross_check(secular_answer,
                    solve_twin(twin_of(detail_of(c)), certpolysolve::Options{.goal = certpolysolve::Goal::approximation()}), 24,
                    53, "Chebyshev T_24 (secular)");
        std::cout << "twins: Chebyshev T_24 secular mode: sweeps " << secular_answer.statistics().sweeps
                  << ", regenerations " << secular_answer.engine_counters().secular_regenerations << '\n';
        certpolysolve::Options tight;
        tight.bits = 200; tight.goal = certpolysolve::Goal::accuracy();
        const auto l = laguerre(12);
        cross_check(solve_program(l, tight), solve_twin(twin_of(detail_of(l)), tight), 12, 200, "Laguerre L_12 (200 bits)");
    }
}

// Planted roots: a non-monic program (leading coefficient 1/1000) with roots on and off the axis, under every goal.
void planted() {
    const std::vector<std::pair<const char *, const char *>> roots{{"1", "0"}, {"0", "2"}, {"-3/2", "-1"}, {"7", "0"}};
    // In Horner form: a product of the four linear factors would be solved factor by factor (the product lift).
    Program p = horner_from_roots(roots);
    p.output = node(p, ProgramOp::Multiply, *p.output, p.append_constant(gaussian("1/1000")));
    check_known(solve_program(p, certpolysolve::Options{}), roots, 53, "planted, approximation");
    check_known(solve_program(p, certpolysolve::Options{.goal = certpolysolve::Goal::accuracy()}), roots, 53,
                "planted, accuracy");
    check_known(solve_program(p, certpolysolve::Options{.bits = 300, .goal = certpolysolve::Goal::accuracy()}), roots, 300,
                "planted, accuracy at 300 bits");
    const auto clusters =
        solve_program(p, certpolysolve::Options{.goal = certpolysolve::Goal::clusters(certpolysolve::Number("1@-30"))});
    require(cps_result_goal(clusters.get()) == CPS_CLUSTER_RADIUS, "cluster goal kept");
    check_known(clusters, roots, 53, "planted, clusters");
    // These roots are not closed under conjugation, so the program has non-real constants and no root is reported
    // real (the reflection proof needs real coefficients: real_by_reflection). A region off the real axis is
    // non-real whatever the coefficients; the two roots on the axis are left unknown.
    unsigned nonreal = 0, unknown = 0;
    for (const auto &r : clusters.roots()) {
        require(r.realness != certpolysolve::Realness::real, "a program's root was reported real");
        nonreal += r.realness == certpolysolve::Realness::not_real;
        unknown += r.realness == certpolysolve::Realness::unknown;
    }
    require(nonreal == 2 && unknown == 2, "two regions off the axis are non-real, two on it unknown");
}

// The start at the roots' radii (program_seeds): the counts on circles and the placement, on planted roots.
void seed_placement() {
    const cps_options options;
    certpolysolve::detail::CandidateStatistics statistics{};
    // Roots of moduli 1, 2, sqrt(13)/2 (about 1.8) and 7: inside 2^-1, 2^0.5, 2^1.5 and 2^3 lie 0, 1, 3 and 4 of them.
    const std::vector<std::pair<const char *, const char *>> roots{{"1", "0"}, {"0", "2"}, {"-3/2", "-1"}, {"7", "0"}};
    const Program p = from_roots(roots);
    const auto prepared = prepare_program(p, options);
    require(prepared.ready() && prepared.degree == 4, "the planted program prepared");
    ProgramSource source(p, 4, prepared.leading);
    for (const auto &[log2_radius, inside] : std::vector<std::pair<double, double>>{{-1, 0}, {0.5, 1}, {1.5, 3}, {3, 4}}) {
        const auto count = circle_count(source, log2_radius, 64, 0, 85, statistics);
        const std::string where = " at 2^" + std::to_string(log2_radius);
        require(count.finite && !count.near_root && !count.cancellation, "a count taken" + where);
        require(std::abs(count.imag) < 0.01 && count.radius < 0.01, "a sharp count" + where);
        require(std::abs(count.real - inside) < 0.01, std::to_string(inside) + " roots inside" + where + ": " + std::to_string(count.real));
    }
    // A root at a sample point is reported; from the samples turned by half a spacing the count is one half for that
    // root (on the circle, its (xi / rho)^q is -1 against the turned samples) plus about 0 for the root outside.
    {
        ComplexBall sample;
        circle_sample(sample.value(), 1, 64, 0, 5);
        Coefficient at_sample;
        at_sample.real = Number::from_double(arf_get_d(arb_midref(acb_realref(sample.value())), ARF_RND_NEAR));
        at_sample.imag = Number::from_double(arf_get_d(arb_midref(acb_imagref(sample.value())), ARF_RND_NEAR));
        Program q;
        const auto x = node(q, ProgramOp::Variable);
        const auto first = node(q, ProgramOp::Subtract, x, q.append_constant(at_sample));
        const auto second = node(q, ProgramOp::Subtract, x, q.append_constant(gaussian("3")));
        q.output = node(q, ProgramOp::Multiply, first, second);
        const auto prepared_q = prepare_program(q, options);
        require(prepared_q.ready() && prepared_q.degree == 2, "the program with a root at a sample prepared");
        ProgramSource source_q(q, 2, prepared_q.leading);
        require(circle_count(source_q, 1, 64, 0, 85, statistics).near_root, "a root at a sample is reported");
        const auto turned = circle_count(source_q, 1, 64, 0.5, 85, statistics);
        require(turned.finite && !turned.near_root && std::abs(turned.real - 0.5) < 0.01,
                "a root on the circle counts one half from the turned samples: " + std::to_string(turned.real));
    }
    // The placement: four seeds from the counts (with four roots no annulus is bisected: the halving ladder's annuli
    // hold at most the capacity), each within a factor sqrt 2 of a planted root's modulus (its annulus's middle),
    // none outside the bound.
    RootVector seeds(4);
    ProgramSeedReport report;
    require(program_seeds(seeds.data(), source, 4, prepared.root_bound.value(), 85, options, statistics, nullptr, &report),
            "the seeds placed");
    require(report.placed == 4 && !report.budget_spent && !report.cancelled,
            "every seed from a count: " + std::to_string(report.placed) + " placed");
    uint64_t in_annuli = 0;
    for (const auto &annulus : report.annuli) in_annuli += annulus.count;
    require(in_annuli == 4, "the annuli hold the degree");
    const double top = mag_get_d_log2_approx(prepared.root_bound.value());
    std::vector<double> moduli{0, 1, std::log2(std::sqrt(13.0) / 2), std::log2(7.0)};
    for (slong i = 0; i < 4; ++i) {
        Magnitude size;
        acb_get_mag(size.value(), seeds.data() + i);
        const double log2_modulus = mag_get_d_log2_approx(size.value());
        require(log2_modulus <= top + 1e-6, "a seed inside the bound");
        const auto nearest = std::min_element(moduli.begin(), moduli.end(), [&](double a, double b) {
            return std::abs(a - log2_modulus) < std::abs(b - log2_modulus);
        });
        require(std::abs(*nearest - log2_modulus) <= 0.5 + 1e-6,
                "a seed within a factor sqrt 2 of a root's modulus: 2^" + std::to_string(log2_modulus));
        moduli.erase(nearest);
    }
    std::cout << "seeds: " << report.circles << " circle counts (" << report.unsharp << " unsharp), " << report.evaluations << " evaluations, "
              << report.annuli.size() << " annuli, " << report.sample_doublings << " sample doublings, " << report.rotations
              << " rotations, " << report.raises << " precision raises\n";
    // The same from the binary64 evaluator (program_machine.hpp): every count read in doubles with the same
    // integer, the same placement; a root at a sample reported from doubles; a circle where the program
    // overflows binary64 read in balls.
    {
        ProgramMachine<double> machine(p);
        require(machine.valid(), "the planted program in binary64");
        for (const auto &[log2_radius, inside] : std::vector<std::pair<double, double>>{{-1, 0}, {0.5, 1}, {1.5, 3}, {3, 4}}) {
            const auto count = circle_count(source, log2_radius, 64, 0, 85, statistics, &machine);
            require(count.binary64 && count.finite && std::abs(count.real - inside) < 0.01 && count.radius < 0.01,
                    "a count read in binary64 at 2^" + std::to_string(log2_radius) + ": " + std::to_string(count.real));
        }
        RootVector again(4);
        ProgramSeedReport in_doubles;
        require(program_seeds(again.data(), source, 4, prepared.root_bound.value(), 85, options, statistics, nullptr,
                              &in_doubles, &machine),
                "the seeds placed from binary64 counts");
        require(in_doubles.binary64_circles == in_doubles.circles && in_doubles.circles == report.circles &&
                    in_doubles.placed == 4,
                "every circle read in binary64: " + std::to_string(in_doubles.binary64_circles) + " of " +
                    std::to_string(in_doubles.circles));
        for (slong i = 0; i < 4; ++i)
            require(acb_equal(again.data() + i, seeds.data() + i), "the same placement from the binary64 counts");
        ComplexBall sample;
        circle_sample(sample.value(), 1, 64, 0, 5);
        Coefficient at_sample;
        at_sample.real = Number::from_double(arf_get_d(arb_midref(acb_realref(sample.value())), ARF_RND_NEAR));
        at_sample.imag = Number::from_double(arf_get_d(arb_midref(acb_imagref(sample.value())), ARF_RND_NEAR));
        Program q;
        const auto x = node(q, ProgramOp::Variable);
        q.output = node(q, ProgramOp::Multiply, node(q, ProgramOp::Subtract, x, q.append_constant(at_sample)),
                        node(q, ProgramOp::Subtract, x, q.append_constant(gaussian("3"))));
        const auto prepared_q = prepare_program(q, options);
        ProgramSource source_q(q, 2, prepared_q.leading);
        ProgramMachine<double> machine_q(q);
        require(circle_count(source_q, 1, 64, 0, 85, statistics, &machine_q).near_root,
                "a root at a sample is reported from binary64");
        const auto big = parse_program("PSPROGRAM 1\nx\npow 0 1100\nc 1\nadd 1 2\nout 3\n"); // x^1100 + 1
        const auto prepared_big = prepare_program(big, options);
        ProgramSource source_big(big, 1100, prepared_big.leading);
        ProgramMachine<double> machine_big(big);
        const auto overflowing = circle_count(source_big, 1, 64, 0, 85, statistics, &machine_big); // 2^1100 overflows
        require(!overflowing.binary64 && overflowing.finite && std::abs(overflowing.real - 1100) < 0.5,
                "a circle where binary64 overflows is read in balls without the range tier: " +
                    std::to_string(overflowing.real));
        ProgramMachine<machine::BlockDouble> block_big(machine_big);
        const auto in_range_tier = circle_count(source_big, 1, 64, 0, 85, statistics, &machine_big, &block_big);
        require(in_range_tier.binary64 && in_range_tier.finite && std::abs(in_range_tier.real - 1100) < 0.25 &&
                    in_range_tier.evaluations == 128,
                "a circle where binary64 overflows is read by the range tier: " + std::to_string(in_range_tier.real) +
                    ", " + std::to_string(in_range_tier.evaluations) + " evaluations");
    }
    // Cancelled: false, and the report says so.
    struct Always {
        static int poll(void *) { return 1; }
    };
    cps_options cancelling = options;
    cancelling.cancel = Always::poll;
    require(!program_seeds(seeds.data(), source, 4, prepared.root_bound.value(), 85, cancelling, statistics, nullptr, &report) &&
                report.cancelled,
            "a cancelled placement says so");
}

// (x - 1)^2 (x + 2): a square factor. ACCURACY ends at the cap with "possibly a multiple root" and the counted
// cluster; APPROXIMATION and the cluster goal answer with a counted cluster of 2.
void square_factor() {
    const auto p = from_roots({{"1", "0"}, {"1", "0"}, {"-2", "0"}});
    const auto approximation = solve_program(p, certpolysolve::Options{});
    check_known(approximation, {{"1", "0"}, {"1", "0"}, {"-2", "0"}}, 53, "square factor, approximation");
    bool pair = false;
    for (const auto &r : approximation.roots()) pair = pair || (r.count == 2 && r.multiplicity.value_or(0) == 0);
    require(pair, "the double root is a counted cluster of 2 with no multiplicity claim");
    // Since phase 2b the separated-roots pass answers it (its cluster step), with no regeneration.
    require(approximation.message().find("Newton disks: a complete disjoint cover") != std::string::npos &&
                approximation.engine_counters().secular_regenerations == 0,
            "the square factor's cluster comes from the pass's counts: " + approximation.message());
    check_known(solve_program(p, certpolysolve::Options{.goal = certpolysolve::Goal::clusters(certpolysolve::Number("1@-20"))}),
                {{"1", "0"}, {"1", "0"}, {"-2", "0"}}, 53, "square factor, clusters");
    certpolysolve::Options accuracy;
    accuracy.bits = 53; accuracy.max_bits = 256; accuracy.goal = certpolysolve::Goal::accuracy();
    const auto stopped = solve_program(p, accuracy);
    require(cps_result_termination(stopped.get()) == CPS_PRECISION_LIMIT &&
                cps_result_stop_reason(stopped.get()) == CPS_STOP_PRECISION && stopped.complete() && !stopped.goal_met(),
            "accuracy on a square factor stops at the cap: " + stopped.message());
    require(stopped.message().find(precision_cap_note) != std::string::npos &&
                stopped.message().find("possibly a multiple root") != std::string::npos,
            "the cap stop names the unproved square-freeness: " + stopped.message());
    // The counted cover at the stop still holds every root.
    const auto boxes = boxes_of(stopped);
    uint64_t total = 0;
    for (const auto &box : boxes) {
        Rational one, minus_two, zero;
        fmpq_one(one.value());
        fmpq_set_si(minus_two.value(), -2, 1);
        const uint64_t inside = (box.contains(one, zero) ? 2u : 0u) + (box.contains(minus_two, zero) ? 1u : 0u);
        require(inside == box.count, "the cap stop's cover counts its roots");
        total += inside;
    }
    require(total == 3, "the cap stop's cover is complete");
    std::cout << "square factor at the cap: " << stopped.message() << '\n';
}

// Phase 2b end to end: a squared composition times a linear factor, (f^6 + g^3 + 5)^2 (x + 3) with f = x^2 + x - 1
// and g = x^3 - 2x + 3, whose output is no power (so the power lift leaves it to the pass): twelve double roots and
// a simple one. The preliminary phase brings two points to each double root, the separated-roots pass counts each
// pair as a cluster of 2 and certifies the cover, and the twin (the expanded polynomial) agrees.
void squared_composition() {
    const auto recorded = certpolysolve::Program::from([](auto x) {
        const auto f = x * x + x - 1;
        const auto g = x * x * x - 2 * x + 3;
        const auto inner = pow(f, 6) + pow(g, 3) + 5;
        return inner * inner * (x + 3);
    });
    // A sum at the top (+ 0, the same polynomial): as a product, the product lift would solve inner and x + 3 alone.
    Program program = detail_of(recorded);
    program.output = node(program, ProgramOp::Add, *program.output, program.append_constant(gaussian("0")));
    const auto twin = twin_of(program);
    const uint64_t degree = twin.coefficients.size() - 1;
    require(degree == 25, "squared composition: degree 25");
    const auto answer = solve_program(program, certpolysolve::Options{});
    certpolysolve::Options twin_options;
    twin_options.goal = certpolysolve::Goal::approximation();
    const auto twin_answer = solve_twin(twin, twin_options);
    cross_check(answer, twin_answer, degree, 53, "squared composition");
    require(answer.message().find("Newton disks: a complete disjoint cover") != std::string::npos &&
                answer.message().find("the output is a power") == std::string::npos &&
                answer.engine_counters().secular_regenerations == 0,
            "squared composition: certified by the pass's counted clusters: " + answer.message());
    size_t clusters = 0, simple = 0;
    for (const auto &r : answer.roots()) {
        if (r.count == 2) {
            require(r.multiplicity.value_or(0) == 0 &&
                        r.certificate == "original-source argument-principle count in a complete disjoint cover",
                    "squared composition: a double root is a counted cluster of 2");
            ++clusters;
        } else {
            require(r.count == 1 && r.multiplicity.value_or(0) == 1, "squared composition: the simple root alone");
            ++simple;
        }
    }
    require(clusters == 12 && simple == 1, "squared composition: twelve clusters and one simple root");
    std::cout << "squared composition: 12 clusters of 2 and a simple root, sweeps " << answer.statistics().sweeps
              << " (twin " << twin_answer.statistics().sweeps << "), work bits " << answer.work_bits() << '\n';
}

// The power lift (program_structure.hpp; 18 October 2026): an output c q^k is solved as q, and q's
// cover is published with counts times k, a multiplicity times k where q's region proves one.
void power_lift() {
    const auto lifted = [](const certpolysolve::Result &answer) {
        return answer.message().find("the output is a power c q^") != std::string::npos;
    };
    // Each planted root in exactly one region, which counts it `times` times.
    const auto holds = [](const certpolysolve::Result &answer,
                          const std::vector<std::pair<const char *, const char *>> &roots, uint64_t times,
                          const std::string &label) {
        const auto boxes = boxes_of(answer);
        uint64_t total = 0;
        for (const auto &box : boxes) {
            uint64_t inside = 0;
            for (const auto &[re, im] : roots) {
                Rational x, y;
                fmpq_set_str(x.value(), re, 10);
                fmpq_set_str(y.value(), im, 10);
                inside += box.contains(x, y);
            }
            require(inside * times == box.count, label + ": a region counts " + std::to_string(box.count) + " for " +
                                                     std::to_string(inside) + " planted roots");
            total += box.count;
        }
        require(total == roots.size() * times, label + ": the counts sum to the degree");
        return boxes;
    };
    {
        // A square against its twin: every region two roots, multiplicity 2 proved.
        const auto program = certpolysolve::Program::from([](auto x) {
            const auto f = x * x + x - 1;
            const auto g = x * x * x - 2 * x + 3;
            const auto inner = pow(f, 6) + pow(g, 3) + 5;
            return inner * inner;
        });
        const auto answer = solve_program(program, certpolysolve::Options{});
        certpolysolve::Options twin_options;
        twin_options.goal = certpolysolve::Goal::approximation();
        const auto twin_answer = solve_twin(twin_of(detail_of(program)), twin_options);
        cross_check(answer, twin_answer, 24, 53, "power lift: a squared composition");
        require(lifted(answer), "power lift: the square is lifted: " + answer.message());
        for (const auto &r : answer.roots())
            require(r.count == 2 && r.multiplicity.value_or(0) == 2 &&
                        r.certificate.find("lifted from the base q") != std::string::npos,
                    "power lift: each region two roots of multiplicity 2");
        std::cout << "power lift: a squared composition of degree 24 in " << answer.statistics().sweeps
                  << " sweeps (twin " << twin_answer.statistics().sweeps << ")\n";
    }
    {
        // The accuracy goal needs each region to be one root of known multiplicity: met on a square, which the
        // pass alone could not (square_factor: a double root stops at the cap).
        // The base in Horner form: as a product of linear factors the product lift would take it.
        Program p = horner_from_roots({{"1", "0"}, {"0", "2"}, {"-3/2", "-1"}});
        p.output = node(p, ProgramOp::Multiply, *p.output, *p.output);
        certpolysolve::Options accuracy;
        accuracy.goal = certpolysolve::Goal::accuracy();
        const auto answer = solve_program(p, accuracy);
        require(answer.complete() && answer.goal_met() && lifted(answer),
                "power lift: the accuracy goal met on a square: " + answer.message());
        for (const auto &box : holds(answer, {{"1", "0"}, {"0", "2"}, {"-3/2", "-1"}}, 2, "power lift, accuracy"))
            require(box.multiplicity == 2, "power lift: each root of the square has multiplicity 2");
    }
    {
        // Nested powers and a cube: ((x + 1)^2)^3 is one root of multiplicity 6; (x^2 + 1)^3 two of multiplicity 3.
        const auto sixfold = solve_program(parse_program("PSPROGRAM 1\nx\nc 1\nadd 0 1\nmul 2 2\npow 3 3\nout 4\n"),
                                           certpolysolve::Options{});
        require(sixfold.complete() && sixfold.goal_met() && lifted(sixfold),
                "power lift: ((x + 1)^2)^3: " + sixfold.message());
        const auto one = holds(sixfold, {{"-1", "0"}}, 6, "power lift, ((x + 1)^2)^3");
        require(one.size() == 1 && one[0].multiplicity == 6, "power lift: -1 of multiplicity 6");
        const auto cube = solve_program(parse_program("PSPROGRAM 1\nx\nmul 0 0\nc 1\nadd 1 2\npow 3 3\nout 4\n"),
                                        certpolysolve::Options{});
        require(cube.complete() && cube.goal_met() && lifted(cube), "power lift: (x^2 + 1)^3: " + cube.message());
        for (const auto &box : holds(cube, {{"0", "1"}, {"0", "-1"}}, 3, "power lift, (x^2 + 1)^3"))
            require(box.multiplicity == 3, "power lift: +-i of multiplicity 3");
    }
    {
        // A power of a Mandelbrot node against its twin.
        const auto m5 = mandelbrot(5);
        Program p = detail_of(m5);
        p.output = node(p, ProgramOp::Power, *p.output, 0, 2);
        const auto answer = solve_program(p, certpolysolve::Options{});
        certpolysolve::Options twin_options;
        twin_options.goal = certpolysolve::Goal::approximation();
        cross_check(answer, solve_twin(twin_of(p), twin_options), 62, 53, "power lift: Mandelbrot k = 5, squared");
        require(lifted(answer), "power lift: Mandelbrot's square is lifted");
        for (const auto &r : answer.roots())
            require(r.count == 2 && r.multiplicity.value_or(0) == 2,
                    "power lift: Mandelbrot's roots of multiplicity 2");
    }
    {
        // A stop inside q: ((x - 1)^2 (x + 2))^2 under the accuracy goal with a cap of 256 bits. q's double root stops
        // q's solve at the cap (as square_factor's does), and the lifted answer says so: the precision limit, the cover
        // complete, the goal not met, the double root's region four roots with its multiplicity still unknown.
        Program p = from_roots({{"1", "0"}, {"1", "0"}, {"-2", "0"}});
        p.output = node(p, ProgramOp::Multiply, *p.output, *p.output);
        certpolysolve::Options accuracy;
        accuracy.bits = 53; accuracy.max_bits = 256; accuracy.goal = certpolysolve::Goal::accuracy();
        const auto stopped = solve_program(p, accuracy);
        require(lifted(stopped) && cps_result_termination(stopped.get()) == CPS_PRECISION_LIMIT && stopped.complete() &&
                    !stopped.goal_met(),
                "power lift: q's cap stop is the output's: " + stopped.message());
        Rational one, minus_two, zero;
        fmpq_one(one.value());
        fmpq_set_si(minus_two.value(), -2, 1);
        uint64_t total = 0;
        for (const auto &box : boxes_of(stopped)) {
            if (box.contains(one, zero))
                require(box.count == 4 && box.multiplicity == 0, "power lift: the double root's region, 4 roots, multiplicity unknown");
            else if (box.contains(minus_two, zero))
                require(box.count == 2 && (box.multiplicity == 0 || box.multiplicity == 2), "power lift: -2's region, 2 roots");
            else
                require(false, "power lift: a region holds no root");
            total += box.count;
        }
        require(total == 6, "power lift: the stopped cover counts the degree");
    }
    {
        // Zero deflation and a power lift must publish together, exactly once, at every stop.
        // q = x^2 (x - 1) (x + 2), in Horner form plus zero so the product lift leaves it whole;
        // p = q^3 has the planted roots 0 (six times), 1 and -2 (three times each).
        const auto p = parse_program("PSPROGRAM 1\nx\nc 1\nadd 0 1\nmul 2 0\nc 2\nsub 3 4\n"
                                     "mul 5 0\nmul 6 0\nc 0\nadd 7 8\npow 9 3\nout 10\n");
        const auto check = [&](const certpolysolve::Result &answer, const std::string &label) {
            require(answer.degree() == 12, label + ": the degree is restored after both lifts");
            if (answer.complete())
                holds(answer, {{"0", "0"}, {"0", "0"}, {"1", "0"}, {"-2", "0"}}, 3, label);
            else
                require(answer.roots().empty() && !answer.goal_met(), label + ": stopped before a cover");
        };
        const auto whole = solve_program(p, stop_at_poll(certpolysolve::Options{}, UINT64_MAX).get());
        require(whole.goal_met() && lifted(whole), "zero and power lift: the completed solve uses the lift");
        check(whole, "zero and power lift: complete");
        const auto units = whole.statistics().work_units;
        const auto polls = cps_testing_result_polls(whole.get());
        require(units > 1 && units <= 64 && polls > 1 && polls <= 4096,
                "zero and power lift: bounded work and cancellation sweep");
        for (uint64_t budget = 1; budget <= units; ++budget) {
            const auto answer = solve_program(p, certpolysolve::Options{.work_budget = budget});
            check(answer, "zero and power lift: work budget " + std::to_string(budget));
            require(answer.statistics().work_units <= budget, "zero and power lift: work budget respected");
            if (!answer.goal_met())
                require(cps_result_stop_reason(answer.get()) == CPS_STOP_WORK_BUDGET,
                        "zero and power lift: a work stop remains a stop after publication");
        }
        for (uint64_t poll = 1; poll <= polls + 1; ++poll) {
            const auto answer = solve_program(p, stop_at_poll(certpolysolve::Options{}, poll).get());
            check(answer, "zero and power lift: cancellation at " + std::to_string(poll));
            if (!answer.goal_met())
                require(cps_result_termination(answer.get()) == CPS_CANCELLED &&
                            cps_result_stop_reason(answer.get()) == CPS_STOP_CANCELLED,
                        "zero and power lift: cancellation remains a stop after publication");
        }
        std::cout << "zero and power lift: planted counts at " << units << " work budgets and " << polls + 1
                  << " cancellation checkpoints\n";
    }
    std::cout << "power lift: squares, a cube, nested powers, a Mandelbrot square and stops inside q\n";
}

// Squaring recurrences read from the graph (program_structure.hpp, the structural seeds' second step): Mandelbrot's
// output is F(x, p_(k-1)) = x p_(k-1)^2 + 1, with one root left over; x q^3 - 2 is one of degree 3 in q; a
// three-term recurrence and a sum of powers of different nodes are none.
void recurrence_structure() {
    const auto m6 = mandelbrot(6);
    const auto &m = detail_of(m6);
    const auto found = squaring_recurrence(m);
    // The recorder's p_0 = x 0 + 1 has degree bound 1, so the bounds run about twice the true degrees: the base's
    // cone is p_5, of degree 31, whatever its bound.
    require(found && found->degree_in_base == 2 && found->degree_in_x == 1,
            "recurrence: Mandelbrot k = 6 is F(x, q) with d = 2 and x degree 1");
    ExpansionArithmetic exact(unlimited);
    auto base = run_program(cone(m, found->base), exact);
    require(base.degree() == 31, "recurrence: Mandelbrot k = 6's base is p_5, of degree 31");
    const auto cubic = squaring_recurrence(parse_program("PSPROGRAM 1\nx\nmul 0 0\nc 1\nadd 1 2\npow 3 3\nmul 0 4\nc 2\nsub 5 6\nout 7\n"));
    require(cubic && cubic->base == 3 && cubic->degree_in_base == 3 && cubic->degree_in_x == 1 && cubic->leftover == 1,
            "recurrence: x (x^2 + 1)^3 - 2 is F(x, q) with d = 3");
    const auto t32 = chebyshev(32);
    require(!squaring_recurrence(detail_of(t32)), "recurrence: a three-term recurrence is none");
    const auto sum = certpolysolve::Program::from([](auto x) {
        const auto f = x * x + x - 1;
        const auto g = x * x * x - 2 * x + 3;
        return pow(f, 6) + pow(g, 3) + 5;
    });
    require(!squaring_recurrence(detail_of(sum)), "recurrence: f^6 + g^3 + 5 is none");
    std::cout << "recurrence structure: Mandelbrot and x q^3 - 2 read, a three-term recurrence and a sum of powers not\n";
}

// The recurrence seeds end to end (program_recurrence_seeds.hpp, the structural seeds' second step): Mandelbrot
// k = 8 and x q^3 - 2 for a Mandelbrot node q of degree 31 are seeded from their bases' roots and agree with their
// twins; Chebyshev T_128, a three-term recurrence, keeps the circle counts.
void recurrence_seeding() {
    const auto seeded = [](const certpolysolve::Result &answer) {
        return answer.message().find("; structural seeds: ") != std::string::npos;
    };
    certpolysolve::Options twin_options;
    twin_options.goal = certpolysolve::Goal::approximation();
    {
        const auto m8 = mandelbrot(8);
        const auto answer = solve_program(m8, certpolysolve::Options{});
        const auto twin_answer = solve_twin(twin_of(detail_of(m8)), twin_options);
        cross_check(answer, twin_answer, 255, 53, "recurrence seeds: Mandelbrot k = 8");
        require(seeded(answer), "recurrence seeds: Mandelbrot k = 8 is seeded from p_7's roots: " + answer.message());
        std::cout << "recurrence seeds: Mandelbrot k = 8 in " << answer.statistics().sweeps
                  << " sweeps, every level included (twin " << twin_answer.statistics().sweeps << ")\n";
    }
    {
        const auto m5 = mandelbrot(5);
        Program p = detail_of(m5);
        const auto q = *p.output;
        const auto cube = node(p, ProgramOp::Power, q, 0, 3);
        const auto x = node(p, ProgramOp::Variable);
        const auto product = node(p, ProgramOp::Multiply, x, cube);
        p.output = node(p, ProgramOp::Subtract, product, p.append_constant(gaussian("2")));
        const auto answer = solve_program(p, certpolysolve::Options{});
        cross_check(answer, solve_twin(twin_of(p), twin_options), 94, 53, "recurrence seeds: x q^3 - 2");
        require(seeded(answer), "recurrence seeds: x q^3 - 2 is seeded from q's roots: " + answer.message());
    }
    {
        const auto t128 = chebyshev(128);
        const auto answer = solve_program(t128, certpolysolve::Options{});
        require(answer.complete() && answer.goal_met() && !seeded(answer),
                "recurrence seeds: a three-term recurrence keeps the circle counts: " + answer.message());
    }
    {
        // A composition through a cut node (the third step): T_16(T_8(x)) = T_128(x), seeded by the exact preimages
        // (T_16's roots s, then T_8(x) = s for each).
        const auto nested = certpolysolve::Program::from([](auto x) {
            auto previous = x * 0 + 1, current = x;
            for (int k = 1; k < 8; ++k) {
                auto next = 2 * x * current - previous;
                previous = current;
                current = next;
            }
            const auto y = current;
            auto outer_previous = y * 0 + 1, outer = y;
            for (int k = 1; k < 16; ++k) {
                auto next = 2 * y * outer - outer_previous;
                outer_previous = outer;
                outer = next;
            }
            return outer;
        });
        const auto answer = solve_program(nested, certpolysolve::Options{});
        cross_check(answer, solve_twin(twin_of(detail_of(nested)), twin_options), 128, 53,
                    "composition seeds: T_16(T_8(x))");
        require(answer.message().find("(1 by a composition's preimages)") != std::string::npos,
                "composition seeds: T_16(T_8(x)) seeded by its preimages: " + answer.message());
        std::cout << "composition seeds: T_16(T_8(x)) in " << answer.statistics().sweeps << " sweeps\n";
    }
    {
        // Coincident seeds: (f^20 + g^13 + 5)^2 (x + 3) reads as F(x, q) = q^2 (x + 3), whose F(r, .) has the double
        // root 0, so the rule places two equal seeds at each root of q. They are spread on a small circle (equal seeds
        // would hold each other in Aberth's sums: this shape at degree 201 ran for minutes before, 18 October 2026),
        // and the pass's counted clusters certify the double roots.
        const auto doubled_recorded = certpolysolve::Program::from([](auto x) {
            const auto f = x * x + x - 1;
            const auto g = x * x * x - 2 * x + 3;
            const auto inner = pow(f, 20) + pow(g, 13) + 5;
            return inner * inner * (x + 3);
        });
        // A sum at the top (+ 0): as a product, the product lift would take it.
        Program doubled = detail_of(doubled_recorded);
        doubled.output = node(doubled, ProgramOp::Add, *doubled.output, doubled.append_constant(gaussian("0")));
        const auto answer = solve_program(doubled, certpolysolve::Options{});
        require(answer.complete() && answer.goal_met() && seeded(answer),
                "recurrence seeds: a squared node times x + 3, seeded and certified: " + answer.message());
        size_t clusters = 0, simple = 0;
        for (const auto &r : answer.roots()) {
            clusters += r.count == 2;
            simple += r.count == 1;
        }
        require(clusters == 40 && simple == 1, "recurrence seeds: 40 clusters of 2 and the simple root -3");
        // A composition whose f has a double root: g^2 (g + 1) with g = x^22 + x + 1, through the cut g; its seeds
        // for the double root s = 0 are spread the same way.
        const auto composed_recorded = certpolysolve::Program::from([](auto x) {
            const auto g = pow(x, 22) + x + 1;
            return g * g * (g + 1);
        });
        Program composed = detail_of(composed_recorded);
        composed.output = node(composed, ProgramOp::Add, *composed.output, composed.append_constant(gaussian("0")));
        const auto composed_answer = solve_program(composed, certpolysolve::Options{});
        cross_check(composed_answer, solve_twin(twin_of(composed), twin_options), 66, 53,
                    "composition seeds: g^2 (g + 1)");
        require(composed_answer.message().find("(1 by a composition's preimages)") != std::string::npos,
                "composition seeds: g^2 (g + 1) seeded by its preimages: " + composed_answer.message());
    }
    std::cout << "recurrence seeds: Mandelbrot k = 8 and x q^3 - 2 seeded, T_128 not\n";
}

// The product lift (output_factors; 18 October 2026): a product of distinct factors solved factor by factor,
// its covers lifted and united; a root shared by two factors sends the output whole to the pipeline; a factor's stop
// at the cap is the product's.
void product_lift() {
    const auto lifted = [](const certpolysolve::Result &answer) {
        return answer.message().find("distinct factors, each solved alone and its cover lifted") != std::string::npos;
    };
    certpolysolve::Options twin_options;
    twin_options.goal = certpolysolve::Goal::approximation();
    {
        // Mandelbrot k = 5 times Chebyshev T_16, against the twin.
        const auto program = certpolysolve::Program::from([](auto x) {
            auto p = x * 0 + 1;
            for (int j = 0; j < 5; ++j)
                p = x * p * p + 1;
            auto previous = x * 0 + 1, current = x;
            for (int k = 1; k < 16; ++k) {
                auto next = 2 * x * current - previous;
                previous = current;
                current = next;
            }
            return p * current;
        });
        const auto answer = solve_program(program, certpolysolve::Options{});
        const auto twin_answer = solve_twin(twin_of(detail_of(program)), twin_options);
        cross_check(answer, twin_answer, 47, 53, "product lift: Mandelbrot k = 5 times T_16");
        require(lifted(answer), "product lift: two factors solved alone: " + answer.message());
        std::cout << "product lift: Mandelbrot k = 5 times T_16 in " << answer.statistics().sweeps
                  << " sweeps (twin " << twin_answer.statistics().sweeps << ")\n";
    }
    {
        // (x^2 + 1)^2 (x^3 - 2) under the accuracy goal: the double roots +-i of multiplicity 2 proved.
        const auto answer = solve_program(
            parse_program("PSPROGRAM 1\nx\nmul 0 0\nc 1\nadd 1 2\npow 3 2\nmul 1 0\nc 2\nsub 5 6\nmul 4 7\nout 8\n"),
            certpolysolve::Options{.goal = certpolysolve::Goal::accuracy()});
        require(answer.complete() && answer.goal_met() && lifted(answer),
                "product lift: the accuracy goal met on (x^2 + 1)^2 (x^3 - 2): " + answer.message());
        size_t doubles = 0, simples = 0;
        for (const auto &r : answer.roots()) {
            if (r.count == 2) {
                require(r.multiplicity.value_or(0) == 2, "product lift: +-i of multiplicity 2");
                ++doubles;
            } else {
                require(r.count == 1 && r.multiplicity.value_or(0) == 1, "product lift: the cube roots of 2 simple");
                ++simples;
            }
        }
        require(doubles == 2 && simples == 3, "product lift: two double roots and three simple ones");
    }
    {
        // A root shared by two factors: (x - 1)(x - 1)(x + 2), two nodes for x - 1. Their covers overlap, so the
        // output is solved whole, square_factor's answer (a counted cluster of 2).
        const auto shared = solve_program(from_roots({{"1", "0"}, {"1", "0"}, {"-2", "0"}}), certpolysolve::Options{});
        require(shared.message().find("their covers overlap (a shared root), so the output is solved whole") !=
                    std::string::npos,
                "product lift: a shared root sends the output whole: " + shared.message());
        check_known(shared, {{"1", "0"}, {"1", "0"}, {"-2", "0"}}, 53, "product lift: a shared root, solved whole");
    }
    {
        // A factor that stops at the cap: (x^2 - 2x + 1)(x + 2) under the accuracy goal, a cap of 256 bits. The double
        // root's factor stops there, and the product's answer says so: the precision limit, the cover complete, the
        // goal not met, the double root's region two roots with its multiplicity unknown.
        const auto stopped = solve_program(
            parse_program("PSPROGRAM 1\nx\nmul 0 0\nc 2\nmul 2 0\nsub 1 3\nc 1\nadd 4 5\nc 2\nadd 0 7\nmul 6 8\nout 9\n"),
            certpolysolve::Options{.bits = 53, .max_bits = 256, .goal = certpolysolve::Goal::accuracy()});
        require(lifted(stopped) && cps_result_termination(stopped.get()) == CPS_PRECISION_LIMIT && stopped.complete() &&
                    !stopped.goal_met(),
                "product lift: a factor's cap stop is the product's: " + stopped.message());
        Rational one, minus_two, zero;
        fmpq_one(one.value());
        fmpq_set_si(minus_two.value(), -2, 1);
        uint64_t total = 0;
        for (const auto &box : boxes_of(stopped)) {
            if (box.contains(one, zero))
                require(box.count == 2 && box.multiplicity == 0, "product lift: the double root's region, multiplicity unknown");
            else
                require(box.contains(minus_two, zero) && box.count == 1, "product lift: -2's region");
            total += box.count;
        }
        require(total == 3, "product lift: the stopped cover counts the degree");
    }
    std::cout << "product lift: two structured factors, multiplicities proved, a shared root solved whole, a cap stop\n";
}

// An answer against the exact roots. The reference holds the expanded polynomial's distinct roots, each in a region
// far smaller than the answer's, with its proved multiplicity. The answer's regions are disjoint; each holds exactly
// its count of roots (every reference region that meets it lies inside it); a multiplicity it claims is that of the
// one distinct root inside; and a complete answer places every root, its counts summing to the degree. An answer
// may name an exact root as a point (the program route's -1 of Mandelbrot's p_5), which no reference region of
// positive width fits inside: a reference region that holds the point, where the expansion vanishes to the
// region's count, holds that root alone, and is replaced by the point.
bool box_inside(const Box &inner, const Box &outer) {
    return fmpq_cmp(outer.lo.value(), inner.lo.value()) <= 0 && fmpq_cmp(inner.hi.value(), outer.hi.value()) <= 0 &&
           fmpq_cmp(outer.ilo.value(), inner.ilo.value()) <= 0 && fmpq_cmp(inner.ihi.value(), outer.ihi.value()) <= 0;
}
// The order to which the exact expansion vanishes at x + iy (0 where it does not), by Horner's rule over Q(i) on
// it and its derivatives.
uint64_t vanishing_order(const GaussianExpansion &expanded, const Rational &x, const Rational &y) {
    GaussianExpansion e = expanded;
    for (uint64_t order = 0; e.degree() >= 0; ++order) {
        Rational re, im, c, t, u;
        for (slong j = e.degree(); j >= 0; --j) {
            fmpq_mul(t.value(), re.value(), x.value());
            fmpq_mul(u.value(), im.value(), y.value());
            fmpq_sub(t.value(), t.value(), u.value());
            fmpq_mul(u.value(), re.value(), y.value());
            fmpq_mul(c.value(), im.value(), x.value());
            fmpq_add(im.value(), u.value(), c.value());
            fmpq_poly_get_coeff_fmpq(c.value(), e.real(), j);
            fmpq_add(re.value(), t.value(), c.value());
            fmpq_poly_get_coeff_fmpq(c.value(), e.imag(), j);
            fmpq_add(im.value(), im.value(), c.value());
        }
        if (!fmpq_is_zero(re.value()) || !fmpq_is_zero(im.value()))
            return order;
        fmpq_poly_derivative(e.real(), e.real());
        fmpq_poly_derivative(e.imag(), e.imag());
    }
    return 0; // the zero polynomial: not a program's output here
}
std::string describe(const Box &b) {
    std::ostringstream out;
    out.precision(17);
    out << "[" << fmpq_get_d(b.lo.value()) << ", " << fmpq_get_d(b.hi.value()) << "] x [" << fmpq_get_d(b.ilo.value())
        << ", " << fmpq_get_d(b.ihi.value()) << "], count " << b.count << ", multiplicity " << b.multiplicity
        << (b.realness == CPS_REAL ? ", real" : b.realness == CPS_NONREAL ? ", not real" : "");
    return out.str();
}
void check_exact_counts(const certpolysolve::Result &answer, std::vector<Box> reference, const GaussianExpansion &expanded,
                        uint64_t degree, const std::string &label) {
    const auto boxes = boxes_of(answer);
    for (auto &r : reference)
        for (const auto &b : boxes)
            if (fmpq_equal(b.lo.value(), b.hi.value()) && fmpq_equal(b.ilo.value(), b.ihi.value()) &&
                r.contains(b.lo, b.ilo) && vanishing_order(expanded, b.lo, b.ilo) == r.count) {
                r.lo = r.hi = b.lo;
                r.ilo = r.ihi = b.ilo;
                break;
            }
    std::vector<uint8_t> placed(reference.size(), 0);
    uint64_t total = 0;
    for (size_t i = 0; i < boxes.size(); ++i) {
        for (size_t j = 0; j < i; ++j)
            require(!boxes[i].overlaps(boxes[j]), label + ": two regions overlap: " + answer.message());
        uint64_t held = 0, distinct = 0, multiplicity = 0;
        for (size_t r = 0; r < reference.size(); ++r) {
            if (!boxes[i].overlaps(reference[r]))
                continue;
            if (!box_inside(reference[r], boxes[i]))
                require(false, label + ": a root's reference region " + describe(reference[r]) +
                                   " crosses the edge of the region " + describe(boxes[i]));
            held += reference[r].count;
            ++distinct;
            multiplicity = reference[r].multiplicity;
            placed[r] = 1;
        }
        require(held == boxes[i].count, label + ": a region counts " + std::to_string(boxes[i].count) +
                                            " roots and holds " + std::to_string(held) + ": " + answer.message());
        if (boxes[i].multiplicity)
            require(distinct == 1 && multiplicity == boxes[i].multiplicity,
                    label + ": a region's multiplicity " + std::to_string(boxes[i].multiplicity) +
                        " is not its root's: " + answer.message());
        total += boxes[i].count;
    }
    if (answer.complete()) {
        require(total == degree, label + ": a complete answer counts " + std::to_string(total) + " roots of " +
                                     std::to_string(degree) + ": " + answer.message());
        for (size_t r = 0; r < reference.size(); ++r)
            require(placed[r], label + ": a complete answer leaves a root outside its regions: " + answer.message());
    }
}


// The product lift under every kind of stop. Until 19 October 2026, a budget that stopped it after its first factor
// published that factor's region at a root the next factor shares, counting 1 where the output has 2 (and the first
// factor's root-bound square, counting 1, where it has all 33). Each answer is checked exactly against the expanded
// polynomial's roots from the coefficient route (the accuracy goal at 256 bits): work budgets of 1 to 120 units, then
// iteration budgets, a cancellation at a checkpoint further on each time, and time limits. A factor's cover is
// incomplete only after a stop before its preparation (a program's answer is complete from its global cover on), so
// the precision cap's stop, which has its own test in product_lift, leaves every cover complete.
void product_lift_stops() {
    struct Case {
        const char *name;
        certpolysolve::Program program;
    };
    std::vector<Case> cases;
    // (x - 1) times (x - 1) r, r Mandelbrot's p_5, its second factor a difference (x r - r), not a product node: the
    // two factors share the root 1, a double root of the output (degree 33).
    cases.push_back({"a root shared by two factors", certpolysolve::Program::from([](auto x) {
                         auto first = x - 1;
                         auto r = x * 0 + 1;
                         for (int j = 0; j < 5; ++j)
                             r = x * r * r + 1;
                         return first * (x * r - r);
                     })});
    // Mandelbrot's p_5 times T_16: no root shared (degree 47).
    cases.push_back({"two factors, no root shared", certpolysolve::Program::from([](auto x) {
                         auto p = x * 0 + 1;
                         for (int j = 0; j < 5; ++j)
                             p = x * p * p + 1;
                         auto previous = x * 0 + 1, current = x;
                         for (int k = 1; k < 16; ++k) {
                             auto next = 2 * x * current - previous;
                             previous = current;
                             current = next;
                         }
                         return p * current;
                     })});
    // (x^2 + 1) T_8 (x^2 + 1) m, m Mandelbrot's p_3, the third factor a sum (x^2 m + m): three factors, the first and
    // the third sharing +-i (degree 19).
    cases.push_back({"three factors, the first and the third sharing +-i", certpolysolve::Program::from([](auto x) {
                         auto first = x * x + 1;
                         auto previous = x * 0 + 1, current = x;
                         for (int k = 1; k < 8; ++k) {
                             auto next = 2 * x * current - previous;
                             previous = current;
                             current = next;
                         }
                         auto m = x * 0 + 1;
                         for (int j = 0; j < 3; ++j)
                             m = x * m * m + 1;
                         return first * current * (x * x * m + m);
                     })});
    for (const auto &c : cases) {
        const std::string name = c.name;
        const auto reference_answer = solve_twin(twin_of(detail_of(c.program)),
                                                 certpolysolve::Options{.bits = 256, .goal = certpolysolve::Goal::accuracy()});
        require(reference_answer.complete() && reference_answer.goal_met(),
                name + ": the reference is incomplete: " + reference_answer.message());
        const auto reference = boxes_of(reference_answer);
        const auto degree = static_cast<uint64_t>(reference_answer.degree());
        ExpansionArithmetic exact(unlimited);
        const auto expanded = run_program(detail_of(c.program), exact);
        const auto check = [&](const certpolysolve::Result &answer, const std::string &label) {
            check_exact_counts(answer, reference, expanded, degree, name + ", " + label);
        };
        {
            const auto whole = solve_program(c.program, certpolysolve::Options{});
            require(whole.complete() && whole.goal_met(), name + ": the unlimited solve: " + whole.message());
            check(whole, "no limit");
        }
        // Every stop after a factor was solved, and how it was answered: the global cover (the fixed path) or the
        // factors' covers, united when they place every root.
        uint64_t stops = 0, global_covers = 0, united = 0;
        const auto tally = [&](const certpolysolve::Result &answer) {
            const auto diagnostics = answer.message();
            if (diagnostics.find("solved before the stop") == std::string::npos ||
                diagnostics.find(", 0 of them solved before the stop") != std::string::npos)
                return;
            ++stops;
            global_covers += diagnostics.find("the answer is the global cover") != std::string::npos;
            united += diagnostics.find("their covers united") != std::string::npos;
        };
        // Work budgets from 1 unit until the budget no longer binds (the solves are deterministic, so a larger
        // budget gives the same answer).
        uint64_t units = 1;
        for (; units <= 120; ++units) {
            const auto answer = solve_program(c.program, certpolysolve::Options{.work_budget = units});
            check(answer, std::to_string(units) + " work units");
            tally(answer);
            if (answer.statistics().work_units < units)
                break;
        }
        require(units <= 120, name + ": 120 work units did not finish the solve");
        // Iteration budgets, then cancellations, each a little later, until the solve no longer stops.
        uint64_t iteration_stops = 0, cancellations = 0;
        for (uint64_t sweeps = 1;; sweeps += 1 + sweeps / 4) {
            const auto answer = solve_program(c.program, certpolysolve::Options{.iteration_budget = sweeps});
            check(answer, std::to_string(sweeps) + " sweeps");
            tally(answer);
            if (cps_result_stop_reason(answer.get()) != CPS_STOP_ITERATION_BUDGET)
                break;
            ++iteration_stops;
        }
        for (uint64_t checkpoints = 0;; checkpoints += 1 + checkpoints / 4) {
            // A stop at the poll after `checkpoints` polls.
            const auto answer =
                solve_program(c.program, stop_at_poll(certpolysolve::Options{}, checkpoints + 1).get());
            check(answer, "cancelled after " + std::to_string(checkpoints) + " checkpoints");
            tally(answer);
            if (cps_result_termination(answer.get()) != CPS_CANCELLED)
                break;
            ++cancellations;
        }
        require(stops > 0 && global_covers > 0 && united > 0,
                name + ": the stops did not reach both answers after a factor (" + std::to_string(global_covers) +
                    " global covers, " + std::to_string(united) + " united)");
        // Time limits (how many stop depends on the machine; each answer is checked either way).
        uint64_t timed_out = 0;
        for (const double seconds : {1e-6, 1e-5, 1e-4, 1e-3, 1e-2}) {
            const auto answer =
                solve_program(c.program, certpolysolve::Options{.time_limit = std::chrono::duration<double>(seconds)});
            check(answer, "a time limit of " + std::to_string(seconds) + " s");
            tally(answer);
            timed_out += cps_result_stop_reason(answer.get()) == CPS_STOP_TIME_LIMIT;
        }
        // Precision caps, under the default goal and the accuracy goal: a factor stopped at its cap keeps a complete
        // cover (its global cover at least), so the lift unites the covers or solves the output whole, and never for
        // an incomplete cover.
        uint64_t capped = 0, capped_lifted = 0;
        for (const bool accuracy : {false, true})
            for (const uint64_t cap : {32, 48, 64, 85, 96, 128, 170, 256}) {
                certpolysolve::Options options;
                options.bits = 53; options.max_bits = cap;
                if (accuracy)
                    options.goal = certpolysolve::Goal::accuracy();
                const auto answer = solve_program(c.program, options);
                check(answer, std::string(accuracy ? "the accuracy goal" : "the default goal") + " at a cap of " +
                                  std::to_string(cap) + " bits");
                require(answer.message().find("a factor's cover is incomplete") == std::string::npos,
                        name + ": a factor's cover incomplete without a stop: " + answer.message());
                if (cps_result_termination(answer.get()) == CPS_PRECISION_LIMIT) {
                    ++capped;
                    capped_lifted += answer.message().find("each solved alone and its cover lifted") !=
                                     std::string::npos;
                }
            }
        require(capped > 0, name + ": no cap stopped the solve");
        // A cap with a work budget. The lift keeps the first termination it meets, so after a factor ends at the cap
        // a later factor's stop is not a stop of the lift: the covers are lifted when complete and disjoint, and
        // otherwise the output is solved whole. A budget always leaves a later factor its preparation (the loop
        // records the stop when no unit is left); a cancellation at that factor's first checkpoint does not, the one
        // way to an incomplete cover without a recorded stop. Stepping a cancellation through every checkpoint
        // under caps of 32 and 48 bits reached it on the second and third products, every count exact, but took
        // minutes (the record of 19 October 2026); it is left to a fuzz family.
        uint64_t combined = 0;
        for (const uint64_t cap : {64, 96})
            for (uint64_t budget = 1; budget <= 120; ++budget) {
                certpolysolve::Options options;
                options.bits = 53;
                options.max_bits = cap;
                options.goal = certpolysolve::Goal::accuracy();
                options.work_budget = budget;
                const auto answer = solve_program(c.program, options);
                check(answer, "a cap of " + std::to_string(cap) + " bits and " + std::to_string(budget) + " work units");
                tally(answer);
                ++combined;
                if (answer.statistics().work_units < budget)
                    break;
            }
        std::cout << "product lift stops (" << name << "): every count exact under work budgets 1-" << units << ", "
                  << iteration_stops + 1 << " iteration budgets, " << cancellations + 1 << " cancellation points, "
                  << timed_out << " of 5 time limits stopping, " << capped << " of 16 caps stopping (" << capped_lifted
                  << " of them lifted), " << combined << " caps with budgets; " << stops << " stops after a factor, "
                  << global_covers
                  << " answered by the global cover, " << united << " by the factors' covers\n";
    }
}

// Every root zero (program_monomial, src/program_preparation.cpp): the preparation proves p = a_n x^n from the majorant
// test at a radius below the constants' denominator bound, and the answer is the point 0, its count and multiplicity
// the degree, at preparation under every goal. Before, the root bound's search went 2^32 binades below Cauchy's bound,
// and under a loose goal the answer was that square, of half-width 2^-(2^32 - 2) (the program fuzz family's seed
// 8817851590760132761, a factor (x - 0)^2 under the cluster goal 2^-49). x^2 written as (x + 1)^2 - 2x - 1, whose
// majorant is too loose for that proof, is proved by its lowest coefficients (program_zero_order: every one below the
// degree zero). The other side: x^2 + 2^-200 x (a root 2^-200 beside zero) is no monomial, its counts exact.
void monomials() {
    const cps_options defaults;
    const auto text = [](const std::string &body) { return "PSPROGRAM 1\n" + body + "\n"; };
    std::string third_power;
    {
        Integer power;
        fmpz_ui_pow_ui(power.value(), 3, 1000);
        char *digits = fmpz_get_str(nullptr, 10, power.value());
        third_power = digits;
        flint_free(digits);
    }
    struct Case {
        std::string label, body;
        uint64_t degree;
    };
    const std::vector<Case> monomial_cases = {
        {"x", "x\nout 0", 1},
        {"x - 0", "x\nc 0\nsub 0 1\nout 2", 1},
        {"x^5, a power", "x\npow 0 5\nout 1", 5},
        {"3 x x x, a chain", "x\nc 3\nmul 1 0\nmul 2 0\nmul 3 0\nout 4", 3},
        {"(x^5 + x) - x^5", "x\npow 0 5\nadd 1 0\nsub 2 1\nout 3", 1},
        {"(-5/7 + 9i/11) x^4", "x\npow 0 4\nc -5/7 9/11\nmul 2 1\nout 3", 4},
        {"x^2 / 3^1000", "x\npow 0 2\nc 1/" + third_power + "\nmul 2 1\nout 3", 2},
        // Its majorant is too loose for the first proof; its lowest coefficients prove it (program_zero_order).
        {"(x + 1)^2 - 2x - 1", "x\nc 1\nadd 0 1\npow 2 2\nc 2\nmul 4 0\nsub 3 5\nsub 6 1\nout 7", 2},
    };
    for (const auto &c : monomial_cases) {
        const auto parsed = parse_program(text(c.body));
        // The point 0 of a program whose constants are all real is labelled real by conjugation (real_by_reflection);
        // the one case with a non-real constant keeps it unknown.
        const bool real_constants = std::all_of(parsed.constants.begin(), parsed.constants.end(),
                                                [](const Coefficient &k) { return k.is_real(); });
        const auto prepared = prepare_program(parsed, defaults);
        require(prepared.ready() && prepared.degree == c.degree && prepared.monomial, c.label + ": not proved a monomial");
        for (int goal = 0; goal < 3; ++goal) {
            certpolysolve::Options options;
            if (goal == 0) options.goal = certpolysolve::Goal::accuracy();
            else if (goal == 1) options.goal = certpolysolve::Goal::approximation();
            else options.goal = certpolysolve::Goal::clusters(certpolysolve::Number("1@-49"));
            const auto r = solve_program(certpolysolve::Program(text(c.body)), options);
            const auto label = c.label + " (goal " + std::to_string(goal) + ")";
            require(cps_result_termination(r.get()) == CPS_SUCCESS && r.goal_met() && r.complete(),
                    label + ": not met at preparation: " + r.message());
            require(r.statistics().work_units == 1, label + ": more work than the preparation's unit");
            require(r.roots().size() == 1, label + ": one region");
            const auto &e = r.roots()[0];
            const certpolysolve::Dyadic zero(0);
            require(e.box.real.lo == zero && e.box.real.hi == zero && e.box.imag.lo == zero && e.box.imag.hi == zero,
                    label + ": not the point 0");
            require(e.count == c.degree && e.multiplicity.value_or(0) == c.degree, label + ": count and multiplicity");
            require(e.certificate == (real_constants ? "program-monomial+conjugation" : "program-monomial"),
                    label + ": certificate " + e.certificate);
            require((e.realness == certpolysolve::Realness::real) == real_constants &&
                        (real_constants || e.realness == certpolysolve::Realness::unknown),
                    label + ": realness");
            require(r.message().find("every root zero") != std::string::npos, label + ": diagnostics");
        }
    }
    // The fuzz family's case: (x - 0)^2 (x - 1) under the cluster goal 2^-49, a product whose first factor is a
    // monomial's power: the point 0 with count and multiplicity 2, beside the root 1.
    {
        certpolysolve::Options options;
        options.goal = certpolysolve::Goal::clusters(certpolysolve::Number("1@-49"));
        const auto r =
            solve_program(certpolysolve::Program(text("x\nc 0\nsub 0 1\npow 2 2\nc -1\nadd 0 4\nmul 3 5\nout 6")), options);
        require(cps_result_termination(r.get()) == CPS_SUCCESS && r.goal_met() && r.roots().size() == 2,
                "(x - 0)^2 (x - 1): " + r.message());
        Rational zero, one;
        fmpq_one(one.value());
        bool point = false, unit = false;
        for (const auto &b : boxes_of(r)) {
            if (fmpq_is_zero(b.lo.value()) && fmpq_is_zero(b.hi.value()) && fmpq_is_zero(b.ilo.value()) &&
                fmpq_is_zero(b.ihi.value()))
                point = b.count == 2 && b.multiplicity == 2;
            else
                unit = b.count == 1 && b.contains(one, zero);
        }
        require(point && unit, "(x - 0)^2 (x - 1): the point 0 counted twice and the root 1: " + r.message());
    }
    // Not a monomial: the certificate does not claim it, and its counts are exact against the coefficient twin.
    for (const auto &[label, body] : std::vector<std::pair<std::string, std::string>>{
             {"x^2 + 2^-200 x", "x\npow 0 2\nc 1@-200\nmul 2 0\nadd 1 3\nout 4"}}) {
        const auto program = parse_program(text(body));
        require(!prepare_program(program, defaults).monomial, label + ": taken for a monomial");
        const auto reference_answer =
            solve_twin(twin_of(program), certpolysolve::Options{.bits = 256, .goal = certpolysolve::Goal::accuracy()});
        require(reference_answer.complete(), label + ": the reference is incomplete: " + reference_answer.message());
        ExpansionArithmetic exact(unlimited);
        const auto expanded = run_program(program, exact);
        const auto answer = solve_program(certpolysolve::Program(text(body)), certpolysolve::Options{});
        require(answer.complete(), label + ": incomplete: " + answer.message());
        check_exact_counts(answer, boxes_of(reference_answer), expanded, 2, label);
        require(answer.message().find("every root zero") == std::string::npos, label + ": diagnostics");
    }
    std::cout << "monomials: " << monomial_cases.size() << " proved under three goals, the point 0 at preparation; "
              << "a non-monomial left to the route, counts exact\n";
}

// Points that stand still end the preliminary phase's window (src/solver.cpp, the window rule): a window in which no
// point moved is a stall. FLINT's log2 of a zero magnitude is a finite -2^62, so before 19 October 2026 the rule
// compared -2^62 with -2^62 - 1, one double, and swept until the work budget (the program fuzz family: 7 of its
// 10,000 cases ran to 100,000 work units, 83 solves to the 10 s deadline). Two of them: a pair 2^-31 apart under a
// cap of 96 bits (799,984 sweeps; now certified in a few work units), and a product whose factors share a root, solved
// whole (at the stall its derivatives have no correct bits, so the stall raises the precision, which it never did).
void still_points() {
    const auto text = [](const std::string &body) { return "PSPROGRAM 1\n" + body + "\n"; };
    struct Case {
        std::string label, body;
        uint64_t cap;
    };
    const std::vector<Case> cases = {
        {"a pair 2^-31 apart (seed 6176160324952821556)",
         "x\nc -5/7 6/11\nmul 1 0\nc -11349786683/82678120448 307845136303/330712481792\nadd 2 3\nmul 4 0\n"
         "c 24425202127829/169324790677504 8449488329705/42331197669376\nadd 5 6\nout 7",
         96},
        {"a shared multiple root (seed 3617104339256384825)",
         "x\nc 1 0\nmul 1 0\nc -505/873 533/407\nadd 2 3\nc 505/873 -533/407\nsub 0 5\nc -1017/1024 -15397/32768\n"
         "sub 0 7\nmul 6 8\nc 15/16 1\nsub 0 10\nmul 9 11\npow 12 2\nmul 4 13\nout 14",
         256},
    };
    for (const auto &c : cases) {
        certpolysolve::Options options;
        options.bits = 53;
        options.max_bits = c.cap;
        options.goal = certpolysolve::Goal::approximation();
        options.work_budget = 2000;
        const auto r = solve_program(certpolysolve::Program(text(c.body)), options);
        const auto sweeps = r.statistics().sweeps;
        require(cps_result_termination(r.get()) == CPS_SUCCESS && r.goal_met(), c.label + ": " + r.message());
        require(r.statistics().work_units <= 200 && sweeps <= 2000,
                c.label + ": " + std::to_string(r.statistics().work_units) + " work units, " +
                    std::to_string(sweeps) + " sweeps");
        const auto program = parse_program(text(c.body));
        const auto reference_answer =
            solve_twin(twin_of(program), certpolysolve::Options{.bits = 256, .goal = certpolysolve::Goal::accuracy()});
        ExpansionArithmetic exact(unlimited);
        check_exact_counts(r, boxes_of(reference_answer), run_program(program, exact),
                           static_cast<uint64_t>(reference_answer.degree()), c.label);
        std::cout << "still points: " << c.label << ", " << r.statistics().work_units << " work units, " << sweeps
                  << " sweeps\n";
    }
    // A product sharing multiple roots, solved whole, under the accuracy goal at 100 bits and a cap of 512 (seed
    // 5177962170700613951): its multiple roots' points are held, and one point's imaginary part fell toward a real root
    // with no floor (2^-198145 by its thousandth window), read as progress until the 30 s deadline. A movement below
    // the point's resolution at the precision is none: the phase ends, and the solve stops at the cap, counts exact.
    {
        const std::string body =
            "x\nc -6825/16777216 0\nsub 0 1\nc 51/32768 1/256\nsub 0 3\nmul 2 4\nc 51/32768 -1/256\nsub 0 6\nmul 5 7\n"
            "c 941/524288 287/52352\nsub 0 9\nmul 8 10\nc 941/524288 -287/52352\nsub 0 12\nmul 11 13\nc -405/94 0\n"
            "mul 15 0\nc -2764125/1577058304 0\nadd 16 17\npow 18 2\nc -39725/131072 -14550021/16777216\nsub 0 20\n"
            "c -39725/131072 14550021/16777216\nsub 0 22\nmul 21 23\nc 447/889 0\nsub 0 25\nmul 24 26\npow 27 2\n"
            "mul 14 19\nmul 29 28\nc 1/16384 0\nmul 31 30\nout 32";
        certpolysolve::Options options;
        options.bits = 100;
        options.max_bits = 512;
        options.goal = certpolysolve::Goal::accuracy();
        options.work_budget = 2000;
        const auto r = solve_program(certpolysolve::Program(text(body)), options);
        require(cps_result_stop_reason(r.get()) == CPS_STOP_PRECISION && r.statistics().work_units < 2000,
                "a real root's imaginary part: " + r.message());
        const auto program = parse_program(text(body));
        // The reference finer than the answer's 512 bits, so that its regions fit inside the answer's.
        const auto reference_answer =
            solve_twin(twin_of(program), certpolysolve::Options{.bits = 1024, .goal = certpolysolve::Goal::accuracy()});
        ExpansionArithmetic exact(unlimited);
        check_exact_counts(r, boxes_of(reference_answer), run_program(program, exact),
                           static_cast<uint64_t>(reference_answer.degree()), "a real root's imaginary part");
        std::cout << "still points: a real root's imaginary part (seed 5177962170700613951), "
                  << r.statistics().work_units << " work units, stopped at the cap\n";
    }
    // The second at a cap of 256 bits raises the phase's precision at its stall ("preliminary: N sweeps, R precision
    // raises, B bits").
    certpolysolve::Options options;
    options.bits = 53;
    options.max_bits = 256;
    options.goal = certpolysolve::Goal::approximation();
    options.work_budget = 2000;
    const auto r = solve_program(certpolysolve::Program(text(cases[1].body)), options);
    const auto diagnostics = r.message();
    const auto clause = diagnostics.find("; preliminary: ");
    require(clause != std::string::npos, "the shared root: no preliminary phase: " + diagnostics);
    const auto comma = diagnostics.find(", ", clause);
    const auto raises = std::stoull(diagnostics.substr(comma + 2));
    require(raises > 0 && diagnostics.find(" precision raises, ", comma) == diagnostics.find(" ", comma + 2),
            "the shared root: the stall raises the phase's precision: " + diagnostics);
}

// A program's region off the real axis is labelled non-real, one meeting the axis is not (certpolysolve.h). The
// refinement's clusters counted on its model (secular-model-count) came unlabelled: here a cluster of four roots
// 2^-25 apart off the axis, in Horner form, under the approximation goal at 16 bits and a cap of 80 (the program
// fuzz family's seed 2528404422577375018, its planted roots).
void off_axis_labels() {
    const auto program = horner_from_roots(
        {{"324/73", "-7593/256"}, {"169869385/38273024", "-31100931/1048576"}, {"339738405/76546048", "-31100927/1048576"},
         {"339738551/76546048", "-31100929/1048576"}, {"3236/105", "0"}, {"1559/64", "-16"}, {"-679/64", "9248/331"},
         {"-27/2", "0"}, {"-45/4", "0"}, {"2711/1024", "8"}, {"-3227/128", "45/2"}, {"12384/221", "372765/65536"},
         {"423377/16384", "14385/512"}, {"-3585/256", "3984/119"}, {"-32", "-31"}, {"2463177/131072", "26767/1024"},
         {"-22816/597", "-256/17"}, {"5056/899", "0"}, {"16695953/524288", "0"}, {"-2896173/131072", "0"},
         {"-33/2", "2528/127"}, {"2685335/131072", "-17696/467"}, {"1580/57", "8197671/262144"}, {"24", "441/16"},
         {"123255/4096", "0"}});
    uint64_t off_axis = 0, model_counts = 0;
    for (const uint64_t cap : {64u, 80u, 96u, 128u}) {
        certpolysolve::Options options;
        options.bits = 16; options.max_bits = cap; options.goal = certpolysolve::Goal::approximation();
        const auto answer = solve_program(program, options);
        for (const auto &e : answer.roots()) {
            const Box b(e);
            const bool meets_axis = fmpq_sgn(b.ilo.value()) <= 0 && fmpq_sgn(b.ihi.value()) >= 0;
            require(meets_axis ? b.realness != CPS_NONREAL : b.realness == CPS_NONREAL,
                    "cap " + std::to_string(cap) + ": a region " + (meets_axis ? "meeting" : "off") +
                        " the axis labelled " + std::to_string(b.realness) + " (" + e.certificate + ")");
            off_axis += !meets_axis;
            model_counts += !meets_axis && e.certificate == "secular-model-count";
        }
    }
    // The path that came unlabelled must still be taken, or the test checks nothing new (cap 80: one model count).
    require(model_counts > 0, "off-axis labels: the refinement's model counts were not reached");
    std::cout << "off-axis labels: " << off_axis << " regions off the axis labelled non-real, " << model_counts
              << " of them the refinement's model counts\n";
}

// A program answer's labels against the program contract (certpolysolve.h): a region off the axis is non-real; a region
// labelled real belongs to a program whose constants are all real, is a real interval, holds one distinct root and
// carries the conjugation mark; and on a complete answer of such a program a region of one distinct root left unknown
// has a reflection in the real axis that meets another region. Returns the number of regions labelled real.
uint64_t check_labels(const certpolysolve::Result &answer, bool real_constants, const std::string &label) {
    const auto boxes = boxes_of(answer);
    const auto &roots = answer.roots();
    // Whether the reflection of b in the real axis meets c.
    const auto meets_reflection = [](const Box &b, const Box &c) {
        Rational lo, hi;
        fmpq_neg(lo.value(), b.ihi.value());
        fmpq_neg(hi.value(), b.ilo.value());
        return fmpq_cmp(b.hi.value(), c.lo.value()) >= 0 && fmpq_cmp(c.hi.value(), b.lo.value()) >= 0 &&
               fmpq_cmp(hi.value(), c.ilo.value()) >= 0 && fmpq_cmp(c.ihi.value(), lo.value()) >= 0;
    };
    uint64_t real = 0;
    for (size_t i = 0; i < boxes.size(); ++i) {
        const auto &b = boxes[i];
        const std::string where = label + ": region " + std::to_string(i) + " (" + roots[i].certificate + ")";
        const bool axis = fmpq_sgn(b.ilo.value()) <= 0 && fmpq_sgn(b.ihi.value()) >= 0;
        const bool one_root = b.count == 1 || b.multiplicity == b.count;
        if (!axis) {
            require(b.realness == CPS_NONREAL, where + ": off the axis but not labelled non-real");
            continue;
        }
        if (b.realness == CPS_REAL) {
            ++real;
            require(real_constants, where + ": labelled real on a program with a non-real constant");
            require(fmpq_is_zero(b.ilo.value()) && fmpq_is_zero(b.ihi.value()), where + ": labelled real off the line");
            require(one_root, where + ": labelled real but may hold distinct roots");
            require(roots[i].certificate.find("+conjugation") != std::string::npos, where + ": no conjugation mark");
            continue;
        }
        require(b.realness == CPS_REAL_UNKNOWN, where + ": meets the axis but labelled non-real");
        if (real_constants && answer.complete() && one_root) {
            bool met = false;
            for (size_t j = 0; j < boxes.size() && !met; ++j)
                met = j != i && meets_reflection(b, boxes[j]);
            require(met, where + ": one distinct root left unknown though its reflection meets no other region");
        }
    }
    return real;
}

// Realness by conjugation on program input (23 October 2026: step (1) of "Real roots on program input";
// src/solver_program.cpp, the finish). A program whose constants are all real has real coefficients, so on the
// complete disjoint cover the finish publishes, a region holding one distinct root whose reflection in the real axis
// meets no other region holds a real root: it is labelled real, its imaginary interval 0, its certificate marked
// "+conjugation". Every answer here is checked against that rule (check_labels), and the real labels are counted
// against the coefficient twin's, which the coefficient route proves the same way.
void real_by_reflection() {
    // Every root of T_12 is real: all twelve labelled real, under the default goal and under the accuracy goal.
    const auto t12 = chebyshev(12);
    for (const auto &goal : {certpolysolve::Goal::approximation(), certpolysolve::Goal::accuracy()}) {
        const auto answer = solve_program(t12, certpolysolve::Options{.goal = goal});
        require(cps_result_goal_met(answer.get()), "T_12 solved");
        require(check_labels(answer, true, "T_12") == 12, "T_12: every root labelled real");
    }
    // Real roots among pairs (Mandelbrot k = 4, degree 15) and a recurrence: as many real labels as the twin's.
    const auto twin_agrees = [](const char *name, const certpolysolve::Program &program) {
        const certpolysolve::Options accuracy{.goal = certpolysolve::Goal::accuracy()};
        const auto answer = solve_program(program, accuracy);
        const auto twin = solve_twin(twin_of(detail_of(program)), accuracy);
        require(cps_result_goal_met(answer.get()) && cps_result_goal_met(twin.get()), std::string(name) + " solved");
        uint64_t twin_real = 0;
        for (const auto &r : twin.roots()) twin_real += r.realness == certpolysolve::Realness::real;
        const uint64_t real = check_labels(answer, true, name);
        require(real == twin_real, std::string(name) + ": " + std::to_string(real) +
                                       " regions labelled real against the twin's " + std::to_string(twin_real));
        return real;
    };
    const uint64_t mandelbrot_real = twin_agrees("Mandelbrot k = 4", mandelbrot(4));
    const uint64_t laguerre_real = twin_agrees("Laguerre L_8", laguerre(8));
    require(mandelbrot_real > 0 && mandelbrot_real < 15 && laguerre_real == 8, "the twins' real roots");
    // The roots i, -i and 2 in Horner form have real constants: 2 is real and the pair non-real. As linear factors
    // the same polynomial has the constants i and -i, and nothing is labelled real: the product lift solves x - 2 on
    // its own, but only the outer finish proves realness, on the whole program's constants.
    const std::vector<std::pair<const char *, const char *>> mixed{{"0", "1"}, {"0", "-1"}, {"2", "0"}};
    const auto horner_mixed = solve_program(horner_from_roots(mixed), certpolysolve::Options{});
    require(check_labels(horner_mixed, true, "i, -i, 2 in Horner form") == 1,
            "i, -i, 2 in Horner form: the root 2 labelled real");
    const auto linear_mixed = solve_program(from_roots(mixed), certpolysolve::Options{});
    require(check_labels(linear_mixed, false, "i, -i, 2 as linear factors") == 0,
            "i, -i, 2 as linear factors: nothing labelled real");
    // (x^2 - 2)^2, a power: the lift's regions hold one distinct root each, multiplicity 2 proved; both real.
    Program square;
    {
        const auto x = node(square, ProgramOp::Variable);
        const auto q = node(square, ProgramOp::Add, node(square, ProgramOp::Multiply, x, x),
                            square.append_constant(gaussian("-2")));
        square.output = node(square, ProgramOp::Power, q, 0, 2);
    }
    const auto lifted = solve_program(square, certpolysolve::Options{.goal = certpolysolve::Goal::accuracy()});
    require(check_labels(lifted, true, "(x^2 - 2)^2") == 2, "(x^2 - 2)^2: both double roots labelled real");
    for (const auto &r : lifted.roots())
        require(r.count == 2 && r.multiplicity.value_or(0) == 2, "(x^2 - 2)^2: multiplicity 2 proved");
    // x^3 (x - 1) in Horner form: the exact point 0 (count and multiplicity 3) and the root 1, both real.
    const auto zero =
        solve_program(horner_from_roots({{"0", "0"}, {"0", "0"}, {"0", "0"}, {"1", "0"}}), certpolysolve::Options{});
    require(check_labels(zero, true, "x^3 (x - 1)") == 2, "x^3 (x - 1): the point 0 and the root 1 labelled real");
    // (x - 1)^2 (x - 3) in Horner form, not a power: the double root's region is a cluster of count 2 whose
    // multiplicity nothing proves (two distinct real roots would look the same), so it stays unknown, and 3 is real.
    // Under the accuracy goal the cluster cannot be resolved: the solve stops at its cap, and the cover it publishes
    // is labelled all the same.
    const auto double_root = horner_from_roots({{"1", "0"}, {"1", "0"}, {"3", "0"}});
    const auto clustered = solve_program(double_root, certpolysolve::Options{});
    require(check_labels(clustered, true, "(x - 1)^2 (x - 3)") == 1, "(x - 1)^2 (x - 3): the root 3 labelled real");
    uint64_t clusters = 0;
    for (const auto &r : clustered.roots()) clusters += r.count == 2 && !r.multiplicity;
    require(clusters == 1, "(x - 1)^2 (x - 3): one cluster of two roots, its multiplicity unknown");
    const auto capped = solve_program(
        double_root, certpolysolve::Options{.bits = 53, .max_bits = 128, .goal = certpolysolve::Goal::accuracy()});
    require(cps_result_termination(capped.get()) == CPS_PRECISION_LIMIT && capped.complete(),
            "(x - 1)^2 (x - 3), accuracy: stopped at the cap with a complete cover");
    require(check_labels(capped, true, "(x - 1)^2 (x - 3) at the cap") == 1,
            "(x - 1)^2 (x - 3) at the cap: the root 3 labelled real");
    // A product whose first factor stops at its cap: a pair 2^-100 apart that 64 bits cannot part, x^2 + 1 and
    // x^2 - x - 6, every constant real. The product lift solves each factor alone and lifts their covers (not the
    // whole output); the answer is the cap's stop, and the factors' covers, united, are labelled: 3 and -2 real, the
    // pair's cluster unknown, +-i non-real. (A lift stopped by a budget or a cancellation, "their covers united", is
    // the program fuzz family's: stops_after_the_cap and product_lift_stops check its counts.)
    const auto product = certpolysolve::Program::from([](auto x) {
        auto pair = x * x - certpolysolve::Number("2535301200456458802993406410753@-100") * x +
                    certpolysolve::Number("1267650600228229401496703205377@-100");
        return pair * (x * x + 1) * (x * x - x - 6);
    });
    const auto product_stop =
        solve_program(product, certpolysolve::Options{.bits = 53, .max_bits = 64, .goal = certpolysolve::Goal::accuracy()});
    require(cps_result_termination(product_stop.get()) == CPS_PRECISION_LIMIT && product_stop.complete() &&
                product_stop.message().find("each solved alone and its cover lifted") != std::string::npos,
            "the product at its cap, through the lift: " + product_stop.message());
    require(check_labels(product_stop, true, "the product at its cap") == 2, "the product at its cap: 3 and -2 real");
    // The outer finish proved them, after the lift (its note precedes the mark): a factor's nested solve proves no
    // region real, so the lift unites, merges and judges the factors' covers as before the proof existed.
    for (const auto &r : product_stop.roots())
        if (r.realness == certpolysolve::Realness::real) {
            const auto lifted = r.certificate.find("; lifted from factor"), mark = r.certificate.find("+conjugation");
            require(lifted != std::string::npos && mark != std::string::npos && mark > lifted,
                    "the product at its cap: proved real by the outer finish, after the lift: " + r.certificate);
        }
    // Every stop on the way: T_12 under the accuracy goal with 1 to 40 work units, the labels by the rule whatever
    // stopped the solve.
    uint64_t stops = 0, stops_with_real = 0;
    for (uint64_t units = 1; units <= 40; ++units) {
        const auto answer =
            solve_program(t12, certpolysolve::Options{.goal = certpolysolve::Goal::accuracy(), .work_budget = units});
        const uint64_t real = check_labels(answer, true, "T_12 with " + std::to_string(units) + " work units");
        if (cps_result_termination(answer.get()) != CPS_SUCCESS) {
            ++stops;
            stops_with_real += real > 0;
        }
    }
    std::cout << "real by reflection: T_12 all real; Mandelbrot k = 4 " << mandelbrot_real
              << " real of 15, as its twin; real labels at a cap stop (one cluster) and a product's cap stop; " << stops
              << " budget stops of T_12 by the rule (" << stops_with_real << " with real labels)\n";
}

// A stop after a factor at the precision cap is the lift's stop (src/solver.cpp, the product lift). Three factors, the
// first a pair 2^-100 apart (expanded, one factor) that ends at a cap of 64 bits under the accuracy goal; a work budget
// or a cancellation then stops a later one. Until 19 October 2026 the lift kept the first termination, the cap's: the
// later factors were left unsolved, and the output was solved whole under "the factors' degrees do not add up" or "a
// factor's cover is incomplete" (the program fuzz family: 331 and 11 solves). Every answer here is the lift's stop
// (the factors' covers united when they place every root, else the global cover), its stop honest, its counts exact.
void stops_after_the_cap() {
    const auto program = certpolysolve::Program::from([](auto x) {
        auto pair = x * x - certpolysolve::Number("2535301200456458802993406410753@-100") * x +
                    certpolysolve::Number("1267650600228229401496703205377@-100"); // (x - 1)(x - 1 - 2^-100)
        auto circle = x * x + 1;
        auto line = x * x - x - 6; // (x - 3)(x + 2)
        return pair * circle * line;
    });
    const auto &detail = detail_of(program);
    const auto reference_answer =
        solve_twin(twin_of(detail), certpolysolve::Options{.bits = 256, .goal = certpolysolve::Goal::accuracy()});
    ExpansionArithmetic exact(unlimited);
    const auto expanded = run_program(detail, exact);
    uint64_t stops = 0, after_the_cap = 0;
    const auto check = [&](const certpolysolve::Result &answer, const std::string &label) {
        const auto d = answer.message();
        for (const char *old : {"do not add up", "cover is incomplete", "a factor refused"})
            require(d.find(old) == std::string::npos, label + ": " + d);
        check_exact_counts(answer, boxes_of(reference_answer), expanded, 6, label);
        const auto stop = cps_result_stop_reason(answer.get());
        // A stop before the preparation (the first checkpoint) never reaches the lift.
        if ((cps_result_termination(answer.get()) == CPS_WORK_LIMIT ||
             cps_result_termination(answer.get()) == CPS_CANCELLED) &&
            d.find("distinct factors") != std::string::npos) {
            ++stops;
            require(d.find("solved before the stop") != std::string::npos, label + ": not the lift's stop: " + d);
            // The first factor ended at the cap before this stop: the stop is the later factor's, not the cap.
            if (d.find(", 0 of them solved before the stop") == std::string::npos &&
                d.find(", 1 of them solved before the stop") == std::string::npos)
                ++after_the_cap;
        }
        if (stop == CPS_STOP_WORK_BUDGET)
            require(answer.statistics().work_units >= 1, label + ": a budget stop without work");
    };
    for (uint64_t units = 1; units <= 60; ++units) {
        certpolysolve::Options options;
        options.bits = 53;
        options.max_bits = 64;
        options.goal = certpolysolve::Goal::accuracy();
        options.work_budget = units;
        const auto answer = solve_program(program, options);
        check(answer, "budget " + std::to_string(units));
        if (cps_result_stop_reason(answer.get()) == CPS_STOP_WORK_BUDGET)
            require(answer.statistics().work_units == units, "budget " + std::to_string(units) + ": unspent");
    }
    const certpolysolve::Options capped{.bits = 53, .max_bits = 64, .goal = certpolysolve::Goal::accuracy()};
    uint64_t checkpoints = 0; // the whole solve's polls: counted, never stopped
    {
        const auto whole = solve_program(program, stop_at_poll(capped, UINT64_MAX).get());
        require(cps_result_termination(whole.get()) == CPS_PRECISION_LIMIT,
                "the first factor ends at the cap: " + whole.message());
        checkpoints = cps_testing_result_polls(whole.get());
    }
    for (uint64_t at = 1; at <= checkpoints; ++at)
        check(solve_program(program, stop_at_poll(capped, at).get()),
              "cancelled at checkpoint " + std::to_string(at));
    require(stops > 0 && after_the_cap > 0, "no stop after the factor at the cap");
    std::cout << "stops after the cap: " << stops << " stops over 60 budgets and " << checkpoints
              << " cancellation points, " << after_the_cap << " of them after the first factor ended at the cap\n";
}

// A root of order v at exactly zero (program_zero_order, src/program_preparation.cpp): proved from the program's lowest
// coefficients, divided out of the source (ProgramSource's zero order), and published as the point 0 with count and
// multiplicity v beside the rest. Before 19 October 2026 a point converging to the root at zero moved less every window
// (no noise floor there), so the preliminary phase never stalled (U_25 under a cap of 32 bits: the program fuzz
// family's seed 7340420364318503455, to the work budget), and Newton disks certified it with exponents near -2.7e15
// (Legendre P_21; seed 6596401102943691335). The other side: a constant term 2^-400 is no root at zero.
void zero_roots() {
    const cps_options defaults;
    const auto recurrence = [](int n, bool legendre) {
        return certpolysolve::Program::from([n, legendre](auto x) {
            auto previous = x * 0 + 1, current = legendre ? x : 2 * x;
            for (int k = 1; k < n; ++k) {
                auto next = legendre ? certpolysolve::Number(std::to_string(2 * k + 1) + "/" + std::to_string(k + 1)) * x *
                                               current -
                                           certpolysolve::Number(std::to_string(k) + "/" + std::to_string(k + 1)) * previous
                                     : 2 * x * current - previous;
                previous = current;
                current = next;
            }
            return current;
        });
    };
    struct Case {
        std::string label;
        certpolysolve::Program program;
        uint64_t order, degree;
        bool structural = false; // seeded by the structure (a squaring recurrence)
    };
    std::vector<Case> cases;
    // The order decided past twice 64 bits (zero_order_precision_limit's comment): the constant terms cancel 317 bits
    // at zero. Until 19 October 2026 (evening) undecided, and the root at zero stopped the solve at the precision cap.
    cases.push_back({"(x + 3^100)^2 - 3^200",
                     certpolysolve::Program(std::string(
                         "PSPROGRAM 1\nx\nc 515377520732011331036461129765621272702107522001\nadd 0 1\npow 2 2\n"
                         "c 265613988875874769338781322035779626829233452653394495974574961739092490901302182994384699044001\n"
                         "sub 3 4\nout 5\n")),
                     1, 2});
    // The Mandelbrot iteration in its c form vanishes at zero: seeded from p_7's roots, less the seed at zero.
    cases.push_back({"Mandelbrot p_(k+1) = p_k^2 + x, k = 8", certpolysolve::Program::from([](auto x) {
                         auto p = x;
                         for (int j = 1; j < 8; ++j) p = p * p + x;
                         return p;
                     }),
                     1, 128, true});
    cases.push_back({"Legendre P_21", recurrence(21, true), 1, 21});
    cases.push_back({"Chebyshev U_25", recurrence(25, false), 1, 25});
    cases.push_back({"x^4 - x^2 (Horner)", certpolysolve::Program(std::string("PSPROGRAM 1\nx\nc 1\nmul 1 0\nmul 2 0\nc -1\nadd 3 4\nmul 5 0\nmul 6 0\nout 7\n")), 2, 4});
    cases.push_back({"(x + 1)^2 - 1", certpolysolve::Program(std::string("PSPROGRAM 1\nx\nc 1\nadd 0 1\npow 2 2\nsub 3 1\nout 4\n")), 1, 2});
    cases.push_back({"(x (x - 1) + x^0 - 1)^3", certpolysolve::Program(std::string("PSPROGRAM 1\nx\nc 1\nsub 0 1\nmul 0 2\nc 1\nadd 3 4\nsub 5 4\npow 6 3\nout 7\n")), 0, 6});
    for (const auto &c : cases) {
        const auto &detail = detail_of(c.program);
        const auto prepared = prepare_program(detail, defaults);
        // The power's order is its base's, read by the power lift's own preparation: the output's is not taken.
        if (c.order)
            require(prepared.ready() && prepared.zero_order == c.order && !prepared.monomial,
                    c.label + ": order " + std::to_string(prepared.zero_order));
        const auto reference_answer =
            solve_twin(twin_of(detail), certpolysolve::Options{.bits = 256, .goal = certpolysolve::Goal::accuracy()});
        ExpansionArithmetic exact(unlimited);
        const auto expanded = run_program(detail, exact);
        for (int goal = 0; goal < 3; ++goal) {
            certpolysolve::Options options;
            if (goal == 0) options.goal = certpolysolve::Goal::accuracy();
            else if (goal == 1) options.goal = certpolysolve::Goal::approximation();
            else options.goal = certpolysolve::Goal::clusters(certpolysolve::Number("1@-49"));
            const auto r = solve_program(c.program, options);
            const auto label = c.label + " (goal " + std::to_string(goal) + ")";
            require(cps_result_termination(r.get()) == CPS_SUCCESS && r.goal_met(), label + ": " + r.message());
            // q's roots go through the binary64 stage, as any program's do (off for a root at zero until 19 October
            // 2026, evening).
            require(r.message().find("; binary64 stage: ") != std::string::npos,
                    label + ": no binary64 stage: " + r.message());
            if (c.structural)
                require(r.message().find("; structural seeds: ") != std::string::npos &&
                            r.message().find(", the 1 nearest zero dropped") != std::string::npos,
                        label + ": not seeded by its structure: " + r.message());
            check_exact_counts(r, boxes_of(reference_answer), expanded, c.degree, label);
            uint64_t at_zero = 0;
            const certpolysolve::Dyadic zero(0);
            for (const auto &e : r.roots())
                if (e.box.real.lo == zero && e.box.real.hi == zero && e.box.imag.lo == zero && e.box.imag.hi == zero) {
                    at_zero += e.count;
                    require(e.multiplicity.value_or(0) == e.count, label + ": the point 0's multiplicity");
                }
            const uint64_t expected = c.order ? c.order : 3; // the power: its base's simple root at zero, three times
            require(at_zero == expected, label + ": the point 0 counts " + std::to_string(at_zero));
        }
    }
    // U_25 under a cap of 32 bits, where the fuzz family's solve spun: it ends well inside a budget.
    {
        certpolysolve::Options options;
        options.bits = 53;
        options.max_bits = 32;
        options.goal = certpolysolve::Goal::approximation();
        options.work_budget = 2000;
        const auto u25 = std::find_if(cases.begin(), cases.end(), [](const Case &c) { return c.label == "Chebyshev U_25"; });
        require(u25 != cases.end(), "U_25 among the cases");
        const auto r = solve_program(u25->program, options);
        require(cps_result_stop_reason(r.get()) != CPS_STOP_WORK_BUDGET && r.statistics().work_units < 2000,
                "U_25 at a cap of 32 bits: " + r.message());
    }
    // A root of q 2^-51 from the root at zero (the program fuzz family's seed 13957833835543924094):
    // -119/111 (x - 1) x (x - 2^-51 - 2^-52 i), with (-1/6 + i/9) x - 2 added and taken away, so that q's values near
    // zero carry an error of about 2^-precision against 2^-51. Under the cluster goal 2^-11 at a cap of 48 bits, q's
    // cover meets its own goal with a region holding zero; until 19 October 2026 (evening) that cover was refused, the
    // refinement had nothing left to do, and the solve swept to its budget of 1,000 units. That region, with the root
    // at zero merged into it, is a counted cluster of two that meets the goal. Under the accuracy goal with an absolute
    // tolerance of 2^-20, the same region meets q's goal by the tolerance at 85 bits, but merged it holds two distinct
    // roots: the refinement works without the tolerance, and the answer is the point 0 and a region of its own.
    {
        const std::string text = "PSPROGRAM 1\nx\nc 1 0\nsub 0 1\nc 1 0\nmul 3 0\nc -1/2251799813685248 -1/4503599627370496\n"
                                 "add 4 5\nmul 6 0\nc 0 0\nadd 7 8\nmul 2 9\nc -119/111 0\nmul 11 10\nc -1/6 1/9\nmul 13 0\n"
                                 "c -2 0\nadd 14 15\nadd 12 16\nsub 17 16\nout 18\n";
        const auto program = parse_program(text);
        const auto reference_answer =
            solve_twin(twin_of(program), certpolysolve::Options{.bits = 1024, .goal = certpolysolve::Goal::accuracy()});
        ExpansionArithmetic exact(unlimited);
        const auto expanded = run_program(program, exact);
        for (const bool accuracy : {false, true}) {
            const certpolysolve::Options options =
                accuracy ? certpolysolve::Options{.bits = 53, .max_bits = 512, .goal = certpolysolve::Goal::accuracy(),
                                              .work_budget = 1000, .absolute_tolerance = certpolysolve::Number("1@-20")}
                         : certpolysolve::Options{.bits = 53, .max_bits = 48,
                                              .goal = certpolysolve::Goal::clusters(certpolysolve::Number("1@-11")),
                                              .work_budget = 1000};
            const auto r = solve_program(certpolysolve::Program(text), options);
            const std::string label = std::string("a root 2^-51 from the root at zero, ") + (accuracy ? "accuracy" : "cluster") + " goal";
            require(cps_result_termination(r.get()) == CPS_SUCCESS && r.goal_met() && r.statistics().work_units <= 50,
                    label + ": " + std::to_string(r.statistics().work_units) + " work units: " + r.message());
            check_exact_counts(r, boxes_of(reference_answer), expanded, 3, label);
            const Rational origin;
            uint64_t holding_zero = 0;
            for (const auto &b : boxes_of(r))
                if (b.contains(origin, origin)) {
                    holding_zero += b.count;
                    if (accuracy)
                        require(fmpq_is_zero(b.lo.value()) && fmpq_is_zero(b.hi.value()) && fmpq_is_zero(b.ilo.value()) &&
                                    fmpq_is_zero(b.ihi.value()) && b.multiplicity == 1,
                                label + ": the root at zero is not the point 0 alone: " + r.message());
                }
            require(holding_zero == (accuracy ? 1u : 2u), label + ": the region at zero counts " + std::to_string(holding_zero));
        }
    }
    // The program fuzz family's seed 8002516022173716811 (harness fuzz-8, `--program --describe`): a product of a
    // degree-8 and a degree-5 polynomial in Horner form, with a degree-8 polynomial Q added and taken away
    // ((P + Q) - Q), degree 13, with a simple root at zero and a double root at 3995/262144 (about 2^-6). Under its
    // options (the cluster goal 2^-4, the classical method, 100 bits, a cap of 512, 27 work units) q = p / x is
    // solved. The cluster step's count on a circle about the double root ran out of arcs on q: q' = (p' - p / x) / x
    // over a ball carries p's value, whose ball loses the cancellation of the last sums (about 100 times as wide as
    // p'). The pass failed, and the refinement answered in 39 units (from 19 October 2026, evening; 8 before the root
    // at zero was divided out). The count reads p itself on a circle that keeps zero outside, where p has q's roots.
    uint64_t fuzz_units = 0;
    {
        const std::string text =
            "PSPROGRAM 1\nx\nc 1 0\nmul 1 0\nc 333727720417/46276280320 0\nadd 2 3\nmul 4 0\n"
            "c 718441532123877459052061/34263905924083508838400 0\nadd 5 6\nmul 7 0\n"
            "c 574401821246784398953559608497/17964154709125894681867059200 0\nadd 8 9\nmul 10 0\n"
            "c 1056085476735972985843215803093180169/37673562976552788283866866935398400 0\nadd 11 12\nmul 13 0\n"
            "c 4471633873627719302101255390246817452457/308621827903920441621437373934783692800 0\nadd 14 15\n"
            "mul 16 0\nc 130338183771828233345847642034730171443/30862182790392044162143737393478369280 0\nadd 17 18\n"
            "mul 19 0\nc 6842673919902988142536084240932780111001/12641150070944581288814074836368740057088 0\n"
            "add 20 21\nmul 22 0\nc 0 0\nadd 23 24\nc 1 0\nmul 26 0\nc 209917/1048576 2953803/16777216\nadd 27 28\n"
            "mul 29 0\nc 4235345612911/17592186044416 17899345481831/17592186044416\nadd 30 31\nmul 32 0\n"
            "c -752323210080645365/2305843009213693952 -779163113517710319/2305843009213693952\nadd 33 34\nmul 35 0\n"
            "c 11815583787879744228375/1208925819614629174706176 11590991184051204686295/1208925819614629174706176\n"
            "add 36 37\nmul 38 0\n"
            "c -366606129416140290013375/4951760157141521099596496896 "
            "-352802523019763595911425/4951760157141521099596496896\n"
            "add 39 40\nmul 25 41\nc -1 0\nmul 43 0\nc 5/6 0\nadd 44 45\nmul 46 0\nc -5/8 2/9\nadd 47 48\nmul 49 0\n"
            "c -2 0\nadd 50 51\nmul 52 0\nc 8/5 0\nadd 53 54\nmul 55 0\nc -5/7 -8\nadd 56 57\nmul 58 0\nc 1/2 -3/5\n"
            "add 59 60\nmul 61 0\nc -1 -3/2\nadd 62 63\nmul 64 0\nc 5/2 -7/3\nadd 65 66\n"
            "add 42 67\nsub 68 67\nout 69\n";
        const auto program = parse_program(text);
        const auto reference_answer =
            solve_twin(twin_of(program), certpolysolve::Options{.bits = 1024, .goal = certpolysolve::Goal::accuracy()});
        ExpansionArithmetic exact(unlimited);
        const auto expanded = run_program(program, exact);
        const auto r = solve_program(certpolysolve::Program(text),
                                     certpolysolve::Options{.bits = 100, .max_bits = 512,
                                                        .goal = certpolysolve::Goal::clusters(certpolysolve::Number("1@-4")),
                                                        .method = certpolysolve::Method::classical, .work_budget = 27,
                                                        .cell_budget = 65536});
        const std::string label = "the program fuzz family's seed 8002516022173716811";
        fuzz_units = r.statistics().work_units;
        require(cps_result_termination(r.get()) == CPS_SUCCESS && r.goal_met(),
                label + ": " + std::to_string(fuzz_units) + " work units: " + r.message());
        require(r.message().find("; Newton disks: a complete disjoint cover") != std::string::npos,
                label + ": not certified by the Newton disks' pass: " + r.message());
        check_exact_counts(r, boxes_of(reference_answer), expanded, 13, label);
    }
    // No root at zero: a constant term of 2^-399 (roots -2^-200 +- 2^-200 i), and x^2 - 2^-400 (roots +-2^-200).
    for (const auto &[label, text, order] : std::vector<std::tuple<std::string, std::string, uint64_t>>{
             {"x^2 + 2^-199 x + 2^-399", "PSPROGRAM 1\nx\npow 0 2\nc 1@-199\nmul 2 0\nadd 1 3\nc 1@-399\nadd 4 5\nout 6\n", 0},
             {"x^2 - 2^-400", "PSPROGRAM 1\nx\npow 0 2\nc 1@-400\nsub 1 2\nout 3\n", 0}}) {
        const auto prepared = prepare_program(parse_program(text), defaults);
        require(prepared.zero_order == order, label + ": order " + std::to_string(prepared.zero_order));
        const auto program = parse_program(text);
        const auto reference_answer =
            solve_twin(twin_of(program), certpolysolve::Options{.bits = 512, .goal = certpolysolve::Goal::accuracy()});
        ExpansionArithmetic exact(unlimited);
        const auto answer =
            solve_program(certpolysolve::Program(text),
                          certpolysolve::Options{.bits = 53, .max_bits = 4096, .goal = certpolysolve::Goal::accuracy()});
        require(answer.complete(), label + ": " + answer.message());
        check_exact_counts(answer, boxes_of(reference_answer), run_program(program, exact), 2, label);
    }
    std::cout << "zero roots: " << cases.size() << " programs with a root at zero under three goals, the point 0 with its "
              << "order; U_25 at 32 bits ends; a root 2^-51 from zero merged under the cluster goal, separated under the "
              << "accuracy goal's tolerance; the fuzz family's seed 8002516022173716811 certified in " << fuzz_units
              << " work units; two near-zero programs left alone\n";
}

// A close pair, 1 and 1 + 2^-100, under the accuracy goal: unresolved at a cap of 96 bits ("distinct roots" by
// the modular image), separated under a cap of 1,024. (Until 18 October 2026 one session raised its cap; the 0.5
// API has one-shot solves only.)
void close_pair_at_the_cap() {
    // In Horner form, (x - 1)(x - 1 - 2^-100) expanded: as a product its two linear factors would be solved one by one
    // (the product lift), which separates them at any cap.
    const Program p = horner_from_roots({{"1", "0"}, {"1267650600228229401496703205377/1267650600228229401496703205376", "0"}});
    const auto stopped =
        solve_program(p, certpolysolve::Options{.bits = 53, .max_bits = 96, .goal = certpolysolve::Goal::accuracy()});
    require(cps_result_termination(stopped.get()) == CPS_PRECISION_LIMIT && !stopped.goal_met() && stopped.complete(),
            "the close pair stops at a 96-bit cap: " + stopped.message());
    require(stopped.message().find("distinct roots, unresolved at the cap") != std::string::npos,
            "the cap stop names the proved distinctness: " + stopped.message());
    const auto finished =
        solve_program(p, certpolysolve::Options{.bits = 53, .max_bits = 1024, .goal = certpolysolve::Goal::accuracy()});
    require(finished.goal_met() && finished.roots().size() == 2,
            "a cap of 1,024 bits separates the pair: " + finished.message());
    Rational one, other, zero;
    fmpq_one(one.value());
    fmpq_one(other.value());
    Rational tiny;
    fmpq_one(tiny.value());
    fmpq_div_2exp(tiny.value(), tiny.value(), 100);
    fmpq_add(other.value(), other.value(), tiny.value());
    uint64_t found = 0;
    for (const auto &box : boxes_of(finished)) {
        require(box.count == 1 && box.multiplicity == 1, "each of the pair is one root");
        found += box.contains(one, zero) + box.contains(other, zero);
    }
    require(found == 2, "the pair's regions hold the roots");
    std::cout << "close pair: " << stopped.message() << '\n';
}

// The program's default goal at 53 and at 200 bits, and through the C entry with null options (the defaults).
// (Until 18 October 2026 one session tightened its goal and had to refuse a change of the query kind.)
void default_goal() {
    const auto program = hermite(10);
    certpolysolve::Options o;
    o.bits = 53;
    const auto first = solve_program(program, o);
    require(first.goal_met() && cps_result_goal(first.get()) == CPS_ROOT_APPROXIMATION, "the default goal at 53 bits");
    o.bits = 200;
    const auto refined = solve_program(program, o);
    require(refined.goal_met() && cps_result_goal(refined.get()) == CPS_ROOT_APPROXIMATION,
            "the default goal at 200 bits");
    cross_check(refined,
                solve_twin(twin_of(detail_of(program)),
                           certpolysolve::Options{.bits = 200, .goal = certpolysolve::Goal::approximation()}),
                10, 200, "Hermite H_10 at 200 bits");
    const cps_options *none = nullptr; // the C entry's defaults
    const auto defaults = solve_program(program, none);
    require(cps_result_goal(defaults.get()) == CPS_ROOT_APPROXIMATION && defaults.goal_met(),
            "null options: the program default");
    // The radius 0 restores CPS_GOAL_DEFAULT (C13), so a program then takes its default, approximation; until stage
    // 3b-3 the radius 0 set CPS_ROOT_ACCURACY, the coefficients' default, for every input kind.
    const auto radius_zero = certpolysolve_test::raw_options();
    require(cps_options_set_cluster_radius(radius_zero.get(), "1e-3") == CPS_OK &&
                cps_options_set_cluster_radius(radius_zero.get(), "0") == CPS_OK &&
                cps_options_goal(radius_zero.get()) == CPS_GOAL_DEFAULT,
            "the radius 0 did not restore the default goal");
    const auto after_zero = solve_program(program, radius_zero.get());
    require(cps_result_goal(after_zero.get()) == CPS_ROOT_APPROXIMATION && after_zero.goal_met(),
            "after the radius 0, a program does not take its default goal");
}

// The refusals of the contract and of the preparation, through the public entry.
void solve_refusals() {
    const auto program = chebyshev(6);
    // Real-only and piecewise: usage errors (C7), the solve's error and no result, whatever the method or goal.
    for (int variant = 0; variant < 4; ++variant) {
        const auto options = certpolysolve_test::raw_options();
        cps_options *o = options.get();
        const bool set = variant == 0   ? cps_options_set_real_only(o, 1) == CPS_OK
                         : variant == 1 ? cps_options_set_real_only(o, 1) == CPS_OK &&
                                              cps_options_set_method(o, CPS_METHOD_SECULAR) == CPS_OK
                         : variant == 2 ? cps_options_set_method(o, CPS_METHOD_PIECEWISE) == CPS_OK
                                        : cps_options_set_real_only(o, 1) == CPS_OK &&
                                              cps_options_set_goal(o, CPS_ROOT_ACCURACY) == CPS_OK;
        require(set, "the refused query's options");
        const auto refused = certpolysolve_test::solve(program, o); // not solve_program, which requires an answer
        require(refused.stop() == certpolysolve::Stop::error && !refused.get() && !refused.message().empty(),
                "real-only or piecewise was not a usage error: " + refused.message());
    }
    // The isolation goal and a region (the real-only route's, D1) are usage errors for programs, under every method,
    // with or without real_only.
    for (int variant = 0; variant < 4; ++variant) {
        const auto options = certpolysolve_test::raw_options();
        cps_options *o = options.get();
        const bool set = variant == 0   ? cps_options_set_goal(o, CPS_ROOT_ISOLATION) == CPS_OK
                         : variant == 1 ? cps_options_set_goal(o, CPS_ROOT_ISOLATION) == CPS_OK &&
                                              cps_options_set_method(o, CPS_METHOD_CLASSICAL) == CPS_OK
                         : variant == 2 ? cps_options_set_region(o, "-1", "1", "", "") == CPS_OK
                                        : cps_options_set_real_only(o, 1) == CPS_OK &&
                                              cps_options_set_region(o, "0", "2", "", "") == CPS_OK;
        require(set, "the refused query's options");
        const auto refused = certpolysolve_test::solve(program, o);
        require(refused.stop() == certpolysolve::Stop::error && !refused.get(),
                "the isolation goal or a region was not a usage error for a program: " + refused.message());
        require(refused.message().find("isolation goal or region") != std::string::npos,
                "the error names the isolation goal and regions: " + refused.message());
    }
    // The preparation's refusals as terminations.
    const auto zero = solve_program(parse_program("PSPROGRAM 1\nx\nsub 0 0\nout 1\n"), certpolysolve::Options{});
    require(cps_result_termination(zero.get()) == CPS_ZERO_POLYNOMIAL && zero.degree() == -1,
            "zero program: " + zero.message());
    const auto undetermined =
        solve_program(parse_program("PSPROGRAM 1\nx\npow 0 100\nadd 1 0\nsub 2 1\nout 3\n"), certpolysolve::Options{});
    require(cps_result_termination(undetermined.get()) == CPS_INDETERMINATE_DEGREE && undetermined.degree() == 100,
            "undetermined program: " + undetermined.message());
    const auto constant =
        solve_program(parse_program("PSPROGRAM 1\nx\nc 1\nadd 0 1\nsub 2 0\nout 3\n"), certpolysolve::Options{});
    require(cps_result_termination(constant.get()) == CPS_SUCCESS && constant.goal_met() && constant.complete() &&
                constant.degree() == 0 && constant.roots().size() == 0,
            "constant program: " + constant.message());
    certpolysolve::ResourceLimits limits;
    limits.exact_operations = 3;
    const auto too_long = solve_program(program, certpolysolve::Options{.resources = limits});
    require(cps_result_termination(too_long.get()) == CPS_UNSUPPORTED_QUERY &&
                cps_result_stop_reason(too_long.get()) == CPS_STOP_INPUT && too_long.preparation_limited(),
            "length over the operations allowed: " + too_long.message());
    certpolysolve::ResourceLimits tiny;
    tiny.temporary_bytes = 16;
    const auto storage = solve_program(program, certpolysolve::Options{.resources = tiny});
    require(cps_result_termination(storage.get()) == CPS_UNSUPPORTED_QUERY &&
                storage.message().find("storage") != std::string::npos,
            "storage over the temporary bytes allowed: " + storage.message());
    // The rule's boundary at the precision cap: the estimate itself is admitted and a byte less is refused. The
    // estimate holds at least the mantissas it names (since 26 October 2026 the region ends' binary mantissas): per
    // root three balls of two parts, and two covers of a ball of two parts and four ends each.
    {
        constexpr uint64_t degree = 6, bits_per_byte = 8, mantissas_per_root = 3 * 2 + 2 * (2 + 4);
        const uint64_t estimate = certpolysolve::detail::program_solve_storage(degree, CPS_DEFAULT_MAX_WORK_BITS);
        require(estimate >= degree * mantissas_per_root * (CPS_DEFAULT_MAX_WORK_BITS / bits_per_byte),
                "the storage estimate holds the mantissas it names");
        certpolysolve::ResourceLimits enough, short_by_one;
        enough.temporary_bytes = estimate;
        short_by_one.temporary_bytes = estimate - 1;
        const auto admitted = solve_program(program, certpolysolve::Options{.resources = enough});
        const auto refused = solve_program(program, certpolysolve::Options{.resources = short_by_one});
        require(admitted.degree() == static_cast<int64_t>(degree) &&
                    admitted.message().find("the solve's storage") == std::string::npos,
                "the storage estimate itself is admitted: " + admitted.message());
        require(cps_result_termination(refused.get()) == CPS_UNSUPPORTED_QUERY &&
                    refused.message().find("the solve's storage") != std::string::npos,
                "a byte below the storage estimate is refused: " + refused.message());
    }
    // No output: an invalid argument at the entry.
    Program incomplete;
    node(incomplete, ProgramOp::Variable);
    cps_program handle;
    handle.program = incomplete;
    cps_result *result = nullptr;
    require(cps_solve_program(&handle, nullptr, &result) == CPS_INVALID_ARGUMENT && !result, "no output refused");
}

// The budgets: the work budget, the iteration budget and cancellation, each stopping one solve; and the solve
// without them. (Until 18 October 2026 a session went on under each budget to the goal.)
void budgets() {
    const auto program = tridiagonal(12);
    certpolysolve::Options one_unit;
    one_unit.work_budget = 1;
    const auto first = solve_program(program, one_unit);
    require(cps_result_termination(first.get()) == CPS_WORK_LIMIT &&
                cps_result_stop_reason(first.get()) == CPS_STOP_WORK_BUDGET && first.complete() &&
                first.roots().size() == 1 && first.roots()[0].count == 12,
            "one work unit prepares and returns the global cover: " + first.message());
    require(first.roots()[0].certificate == "program-majorant-bound", "the global cover's certificate");
    const auto result = solve_program(program, certpolysolve::Options{});
    require(result.goal_met(), "the unlimited solve reaches the goal");
    cross_check(result,
                solve_twin(twin_of(detail_of(program)), certpolysolve::Options{.goal = certpolysolve::Goal::approximation()}),
                12, 53, "tridiagonal 12");
    // The iteration budget: one sweep.
    certpolysolve::Options one_sweep;
    one_sweep.iteration_budget = 1;
    const auto swept = solve_program(program, one_sweep);
    require(cps_result_termination(swept.get()) == CPS_WORK_LIMIT &&
                cps_result_stop_reason(swept.get()) == CPS_STOP_ITERATION_BUDGET && swept.statistics().sweeps == 1,
            "one sweep stops for the iteration budget: " + swept.message());
    // Cancellation: at once (the first poll), and after 40 polls (the 41st).
    const auto cancelled = solve_program(program, stop_at_poll(certpolysolve::Options{}, 1).get());
    require(cps_result_termination(cancelled.get()) == CPS_CANCELLED &&
                cps_result_stop_reason(cancelled.get()) == CPS_STOP_CANCELLED,
            "cancelled at once: " + cancelled.message());
    const auto cancelled_later = solve_program(program, stop_at_poll(certpolysolve::Options{}, 41).get());
    require(cps_result_termination(cancelled_later.get()) == CPS_CANCELLED &&
                cps_testing_result_polls(cancelled_later.get()) > 40,
            "cancelled after 40 polls: " + cancelled_later.message());
}

// The number that ends just before `end` in `text` (the diagnostics print counts without separators).
uint64_t number_ending_at(const std::string &text, size_t end) {
    size_t start = end;
    while (start > 0 && std::isdigit(static_cast<unsigned char>(text[start - 1]))) --start;
    require(start < end, "a number before position " + std::to_string(end) + " of: " + text);
    return std::stoull(text.substr(start, end - start));
}

// The preliminary phase's precision on Laguerre L_128, whose recurrence cancels about 100 bits at its large
// roots (program_advance, 16 October 2026): the root counts read most circles in binary64 and need no raise;
// the binary64 stage converges all but a straggler at the largest roots (its noise is above the floor
// estimate; the per-point stagnation rule stops it); the handover reads the derivative balls there, which
// have no correct bits at 85 bits, and raises once, to 170; the ball phase then finishes the straggler in a
// few sweeps with the others settled. The sweeps stay within twice the twin's and eight.
void precision_hint() {
    const auto program = laguerre(128);
    const auto answer = solve_program(program, certpolysolve::Options{});
    const std::string diagnostics = answer.message();
    const size_t seeds = diagnostics.find("; seeds: "), phase = diagnostics.find("; preliminary: ");
    require(seeds != std::string::npos && phase != std::string::npos && seeds < phase,
            "L_128's diagnostics name the seeds and the preliminary phase: " + diagnostics);
    const size_t counts_bits = diagnostics.find(" bits;", seeds), phase_bits = diagnostics.find(" bits;", phase);
    const size_t phase_raises = diagnostics.find(" precision raises, ", phase);
    const size_t phase_sweeps = diagnostics.find(" sweeps, ", phase);
    const size_t circles = diagnostics.find(" circle counts (", seeds), in_binary64 = diagnostics.find(" in binary64, ", seeds);
    require(counts_bits != std::string::npos && counts_bits < phase && phase_bits != std::string::npos &&
                phase_raises != std::string::npos && phase_sweeps != std::string::npos && phase_sweeps < phase_raises &&
                phase_raises < phase_bits && circles != std::string::npos && in_binary64 != std::string::npos &&
                circles < in_binary64 && in_binary64 < counts_bits,
            "L_128's diagnostics give the counts' bits and circles and the phase's sweeps, raises and bits: " + diagnostics);
    const uint64_t counts = number_ending_at(diagnostics, counts_bits);
    const uint64_t raises = number_ending_at(diagnostics, phase_raises), bits = number_ending_at(diagnostics, phase_bits);
    const uint64_t circle_count = number_ending_at(diagnostics, circles), binary64_count = number_ending_at(diagnostics, in_binary64);
    require(counts == 85 && 2 * binary64_count >= circle_count && circle_count > 0,
            "L_128's counts read most circles in binary64 (the outermost overflow it) at the initial precision: " +
                std::to_string(binary64_count) + " of " + std::to_string(circle_count) + ", " + std::to_string(counts) + " bits");
    require(diagnostics.find("; binary64 stage: ") != std::string::npos, "L_128's diagnostics name the binary64 stage: " + diagnostics);
    require(raises == 1 && bits == 170,
            "L_128's phase raised once, to 170 bits, from the handover read: " + std::to_string(raises) + " raises, " +
                std::to_string(bits) + " bits");
    // The twin from the circles: since 20 October 2026 Auto starts an exact real factor of degree 64 and up whose
    // roots are all real from the all-real seed pass's seeds (L_128's twin in 1 sweep), which program input does not
    // have; forced secular, the coefficient path Auto takes here, keeps the circles (83 sweeps that day).
    certpolysolve::Options twin_options;
    twin_options.goal = certpolysolve::Goal::approximation();
    twin_options.method = certpolysolve::Method::secular;
    const auto twin_answer = solve_twin(twin_of(detail_of(program)), twin_options);
    cross_check(answer, twin_answer, 128, 53, "Laguerre L_128 (the precision hint)");
    const uint64_t program_sweeps = answer.statistics().sweeps;
    const uint64_t twin_sweeps = twin_answer.statistics().sweeps;
    require(program_sweeps <= 2 * twin_sweeps + 8,
            "L_128's sweeps (" + std::to_string(program_sweeps) + ") within twice the twin's (" + std::to_string(twin_sweeps) +
                ") and eight");
    std::cout << "precision hint: Laguerre L_128: counts " << counts << " bits (" << binary64_count << " of " << circle_count
              << " circles in binary64), phase " << raises << " raise to " << bits << " bits, sweeps " << program_sweeps
              << " (twin " << twin_sweeps << ")\n";
}

// The Mandelbrot ladder, untimed: a probe for milestone 4 (test_program --mandelbrot K), each rung against its twin.
int mandelbrot_ladder(int top) {
    for (int k = 2; k <= top; ++k) {
        const auto program = mandelbrot(k);
        const auto start = std::chrono::steady_clock::now();
        const auto answer = solve_program(program, certpolysolve::Options{});
        const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
        const auto s = answer.engine_counters();
        std::cout << "mandelbrot k=" << k << " degree " << ((1 << k) - 1) << ": " << seconds << " s (untimed), goal "
                  << (answer.goal_met() ? "met" : "NOT met") << ", sweeps " << answer.statistics().sweeps
                  << ", regenerations " << s.secular_regenerations << ", regeneration " << s.regeneration_ns / 1e9
                  << " s, candidates " << s.candidate_ns / 1e9 << " s, validation " << s.validation_ns / 1e9
                  << " s, work bits " << answer.work_bits() << ", " << answer.roots().size() << " regions\n"
                  << "  " << answer.message() << '\n';
        if (k <= 11) {
            const auto twin = twin_of(detail_of(program));
            const auto twin_start = std::chrono::steady_clock::now();
            const auto reference = solve_twin(twin, certpolysolve::Options{.goal = certpolysolve::Goal::approximation()});
            std::cout << "  twin: "
                      << std::chrono::duration<double>(std::chrono::steady_clock::now() - twin_start).count()
                      << " s (untimed), sweeps " << reference.statistics().sweeps << ", regenerations "
                      << reference.engine_counters().secular_regenerations << ", work bits " << reference.work_bits()
                      << '\n';
            cross_check(answer, reference, twin.coefficients.size() - 1, 53, "Mandelbrot k=" + std::to_string(k));
        }
    }
    return 0;
}
} // namespace

int main(int argc, char **argv) {
    try {
        if (argc == 3 && std::string(argv[1]) == "--mandelbrot") return mandelbrot_ladder(std::stoi(argv[2]));
        if (argc == 4 && std::string(argv[1]) == "--budget") {
            // test_program --budget FILE UNITS: a program's diagnostics after UNITS work units (a probe).
            std::ifstream in(argv[2]);
            const std::string text((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
            certpolysolve::Options budget;
            budget.work_budget = static_cast<uint64_t>(std::stoull(argv[3]));
            const auto result = solve_program(certpolysolve::Program{text}, budget);
            std::cout << "termination " << cps_result_termination(result.get()) << ", stop "
                      << cps_result_stop_reason(result.get()) << ", goal met " << result.goal_met() << ", sweeps "
                      << result.statistics().sweeps << ", evaluations " << result.engine_counters().evaluations
                      << ", work units " << result.statistics().work_units << "\n"
                      << result.message() << "\n";
            return 0;
        }
        text_form();
        prime();
        random_programs();
        structure();
        cancellations();
        liveness();
        c_api();
        recorder();
        binary64();
        refusals();
        modular_verdicts();
        root_bounds();
        solve_refusals();
        planted();
        seed_placement();
        machine_evaluator();
        block_double();
        machine_bits();
        machine_storage();
        deflated_evaluator();
        noise_floor();
        machine_stage();
        separated_roots();
        winding_counts();
        square_factor();
        squared_composition();
        power_lift();
        product_lift();
        product_lift_stops();
        monomials();
        still_points();
        off_axis_labels();
        real_by_reflection();
        stops_after_the_cap();
        zero_roots();
        recurrence_structure();
        recurrence_seeding();
        close_pair_at_the_cap();
        default_goal();
        budgets();
        twins();
        precision_hint();
    } catch (const std::exception &error) {
        std::cerr << "FAIL " << error.what() << '\n';
        return 1;
    }
    std::cout << "program tests passed\n";
    return 0;
}
