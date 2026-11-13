// What a program solve needs before it starts, proved without expanding the program (prepare_program, declared in
// program_preparation.hpp with its parts): the exact degree and leading coefficient, the modular image's verdict on
// square-freeness, whether every root is zero or else the order of a root at zero, and a bound of the roots' moduli; or
// a refusal with its diagnostic. program_solve_storage prices a solve's resident storage for the resource rule.
//
// The stages of prepare_program, in order. The program must have an output and at most exact_operations nodes. The top
// coefficients decide the degree and the leading coefficient (top_of_program, at widths up to top_width_limit), or
// refuse: the zero polynomial, leading terms that cancel past the widest series, a constant over the exact storage
// budget; a nonzero constant is ready with degree 0. For degree 1 or more: the modular verdict
// (modular_squarefree_image: priced first by modular_cost against exact_operations and temporary_bytes, then at most
// squarefree_tests tests over the first modular_prime_limit primes from gaussian_prime, skipping a prime that divides a
// denominator or the leading coefficient; square-free when the image is coprime to its derivative, and never a proof of
// a repeated root); every root zero (program_monomial: the majorant test at R = 2^-(D+L+1), with D from
// DenominatorBits, an upper bound of the bits of a common denominator of the coefficients, and |a_n| < 2^L); otherwise
// the order of a root at zero (program_zero_order: the lowest coefficients from the program over truncated power series
// at zero, LowSeries, at precisions doubling from 64 bits, a coefficient zero when its ball lies within 2^-D of zero
// and nonzero when it excludes zero; an order equal to the degree is a monomial too); then, unless every root is zero,
// the root bound (program_root_bound: Cauchy's bound from the majorant at 1, then the least power of two R below it
// where the majorant shows M(R) < 2 |a_n| R^n, by gallop_steps doubling steps and bisection).
//
// State: none.
// Entry points: prepare_program, from solver_program.cpp (program_advance: the output, and a power output's base) and
// program_recurrence_seeds.cpp (the sub-programs a structural seed solves first); program_solve_storage, from
// solver_program.cpp (program_advance's storage admission). The other functions are reached only through
// prepare_program.
// Not here: the arithmetics run here (TopArithmetic, ModularArithmetic and MajorantArithmetic, in program.cpp), the
// solve and what the preparation's answers change in it (program_advance in solver_program.cpp).
// Design notes: docs/slp-input-spec-20261015.md (section 4, the degree and leading coefficient; section 5, the modular
// image and its cost rule; section 7, the root bound; section 3.1, the resource limits).
#include "numeric/evaluate/program_arithmetic.hpp"
#include "numeric/evaluate/program_interpreter.hpp"
#include "numeric/exact/program_preparation.hpp"
#include "core/cover.hpp"
#include "core/resource_limits.hpp"
#include <algorithm>
#include <cmath>
#include <flint/fmpz.h>
#include <flint/mag.h>
#include <flint/nmod_poly.h>
#include <bit>
#include <stdexcept>

namespace certpolysolve::detail {
namespace {
using resources::add;
using resources::multiply;
constexpr uint64_t saturated = resources::unlimited;

// Word operations of a product of the given length by FFT: the length times its bit count.
uint64_t product_operations(uint64_t length) {
    return multiply(length, static_cast<uint64_t>(std::bit_width(length)));
}

// The widest top series tried (section 4). A cost hypothesis, not a mathematical limit: a run of width m costs
// about m^2 / 2 Gaussian products a product node (2,048 at 64). A program whose leading terms still cancel over
// 63 coefficients at some node is refused ("leading coefficient not determined"), the case the spec expects to be
// rare.
constexpr size_t top_width_limit = 64;

// The modular image's primes (section 5): a fixed sequence of 63-bit primes just above 2^62, the next one where a
// constant's denominator or the leading coefficient vanishes modulo the last (a program has to be built for that), and
// a second square-freeness test before "unknown". No answer depends on the verdict in this version.
constexpr unsigned modular_prime_limit = 4;
constexpr unsigned squarefree_tests = 2;

// The root bound's search below Cauchy's (section 7): halving steps of 1, 2, 4, ... in the exponent, then
// bisection. A polynomial whose test holds 2^32 binades below Cauchy's bound (all its roots at or next to 0) keeps
// that bound: past any use as the solve's start circle. A program whose roots are all exactly zero is told apart
// before the search (program_monomial).
constexpr int gallop_steps = 32;

// An upper bound, in bits, of a common denominator of a node's coefficients: each coefficient is (u + i v) / Q with
// integers u and v and one Q of at most this many bits. A constant's Q is the product of its parts' denominators (a
// significand's, times radix^k for an exponent -k: ceil(log2 radix) bits a unit of k), a sum's, a difference's and
// a product's the product of the operands' Q, a power's the operand's Q to the exponent, the variable's 1.
// Saturates (no bound).
class DenominatorBits {
    static uint64_t bits(const Number &n) {
        if (n.zero()) return 0;
        const fmpz *denominator = fmpq_denref(n.significand.value());
        uint64_t result = fmpz_is_one(denominator) ? 0 : fmpz_bits(denominator);
        if (fmpz_sgn(n.exponent.value()) < 0) {
            Integer units;
            fmpz_neg(units.value(), n.exponent.value());
            if (!fmpz_abs_fits_ui(units.value())) return saturated;
            const auto per_unit = static_cast<uint64_t>(std::bit_width(n.radix - 1u));
            result = resources::add(result, resources::multiply(fmpz_get_ui(units.value()), per_unit));
        }
        return result;
    }

  public:
    using Value = uint64_t;
    Value variable() { return 0; }
    Value constant(const Coefficient &c) { return resources::add(bits(c.real), bits(c.imag)); }
    Value add(Value a, Value b) { return resources::add(a, b); }
    Value subtract(Value a, Value b) { return resources::add(a, b); }
    Value multiply(Value a, Value b) { return resources::add(a, b); }
    Value negate(Value a) { return a; }
    Value power(Value a, uint64_t exponent) { return resources::multiply(a, exponent); }
};

// A value's lowest `width` coefficients, as complex balls (index 0 the value at zero): the program over power series
// at zero, truncated (program_zero_order).
class LowSeries {
    size_t width_;
    slong precision_;

  public:
    using Value = std::vector<ComplexBall>;
    LowSeries(size_t width, slong precision) : width_(width), precision_(precision) {}
    Value variable() {
        Value v(width_);
        if (width_ > 1) acb_one(v[1].value());
        return v;
    }
    Value constant(const Coefficient &c) {
        Value v(width_);
        c.enclose(v[0].value(), precision_);
        return v;
    }
    Value add(const Value &a, const Value &b) {
        Value v(width_);
        for (size_t k = 0; k < width_; ++k) acb_add(v[k].value(), a[k].value(), b[k].value(), precision_);
        return v;
    }
    Value subtract(const Value &a, const Value &b) {
        Value v(width_);
        for (size_t k = 0; k < width_; ++k) acb_sub(v[k].value(), a[k].value(), b[k].value(), precision_);
        return v;
    }
    Value negate(const Value &a) {
        Value v(width_);
        for (size_t k = 0; k < width_; ++k) acb_neg(v[k].value(), a[k].value());
        return v;
    }
    Value multiply(const Value &a, const Value &b) {
        Value v(width_);
        for (size_t k = 0; k < width_; ++k)
            for (size_t i = 0; i <= k; ++i) acb_addmul(v[k].value(), a[i].value(), b[k - i].value(), precision_);
        return v;
    }
    Value power(const Value &a, uint64_t exponent) {
        Value result(width_), base = a;
        acb_one(result[0].value());
        for (; exponent; exponent >>= 1) {
            if (exponent & 1) result = multiply(result, base);
            if (exponent > 1) base = multiply(base, base);
        }
        return result;
    }
};
// The order of a root at zero is decided from the lowest coefficients, widened 1, 2, 4, ... up to the top series' limit
// (an order of 64 or more makes no claim), at precisions doubling from 64 bits to twice D + 64 (a coefficient bound
// D above 2^19 bits makes none either: the evaluations would cost more than the claim saves), or further, to cover
// the terms a coefficient cancels: its radius at 64 bits, 2^b, says its terms carry about b bits above the unit, so
// b + 64 + D bits, and a guard of 64, bring the radius below 2^-D; never past 2^20 bits, the most the first rule
// reaches. The second rule is for a coefficient that cancels terms far larger than itself at zero, which the first
// rule's precision leaves undecided: an undecided root at zero, which no iteration lands on exactly, would stop
// every solve at the precision cap.
// History: docs/slp-input-spec-20261015.md, section 13.
constexpr uint64_t zero_order_denominator_limit = uint64_t{1} << 19;
constexpr slong zero_order_precision_limit = slong{1} << 20;
} // namespace

ModularCost modular_cost(const Program &program) {
    if (!program.output || *program.output >= program.nodes.size()) return {saturated, saturated};
    const auto degrees = program.degree_bounds();
    for (const auto d : degrees)
        if (d == saturated) return {saturated, saturated};
    const auto last = program.last_uses();
    const auto length = [&](uint64_t j) { return degrees[j] + 1; };
    ModularCost cost;
    uint64_t live = 0; // words held, released as run_program releases values
    for (size_t i = 0; i < program.nodes.size(); ++i) {
        const auto &node = program.nodes[i];
        const uint64_t own = length(i);
        uint64_t operations = own; // the variable, a constant, a sum, a difference, a negation
        if (node.op == ProgramOp::Multiply) operations = product_operations(own);
        if (node.op == ProgramOp::Power) operations = multiply(2, product_operations(own));
        cost.operations = add(cost.operations, operations);
        live = add(live, own);
        if (live == saturated) return {saturated, saturated};
        cost.words = std::max(cost.words, live);
        const auto release = [&](uint64_t j) {
            if (last[j] == i && j != *program.output) live -= length(j);
        };
        if (node.op != ProgramOp::Variable && node.op != ProgramOp::Constant) {
            release(node.a);
            if (node.op != ProgramOp::Negate && node.op != ProgramOp::Power && node.b != node.a) release(node.b);
        }
        if (last[i] == program.nodes.size() && i != *program.output) live -= own;
    }
    const uint64_t n = length(*program.output);
    cost.operations = add(cost.operations, multiply(product_operations(n), static_cast<uint64_t>(std::bit_width(n))));
    cost.words = add(cost.words, multiply(2, n));
    return cost;
}

ModularVerdict modular_squarefree_image(const Program &program, uint64_t degree, const cps_options &options) {
    ModularVerdict out;
    const auto cost = modular_cost(program);
    if (!options.exact_effort_fits(cost.operations) ||
        !resources::storage_fits(multiply(cost.words, sizeof(ulong)), options.temporary_bytes)) {
        out.note = "modular image skipped: over the exact operations or the temporary storage allowed";
        return out;
    }
    unsigned tests = 0;
    for (unsigned index = 0; index < modular_prime_limit && tests < squarefree_tests; ++index) {
        if (options.cancelled()) {
            out.cancelled = true;
            out.note = "modular image cancelled";
            return out;
        }
        const auto prime = gaussian_prime(index);
        ModularArithmetic arithmetic(prime);
        const auto image = run_program(program, arithmetic);
        if (arithmetic.failed()) continue;                                         // a denominator vanishes mod p
        if (nmod_poly_degree(image.value()) != static_cast<slong>(degree)) continue; // p divides a_n
        ++tests;
        ModularPolynomial derivative(prime.prime), gcd(prime.prime);
        nmod_poly_derivative(derivative.value(), image.value());
        nmod_poly_gcd(gcd.value(), image.value(), derivative.value());
        if (nmod_poly_degree(gcd.value()) == 0) {
            out.squarefree = true;
            out.note = "square-free modulo the prime " + std::to_string(prime.prime);
            return out;
        }
    }
    out.note = tests ? "not square-free modulo " + std::to_string(tests) + " primes: unknown"
                     : "no prime kept the degree and the constants: unknown";
    return out;
}

Magnitude program_root_bound(const Program &program, uint64_t degree, const GaussianRational &leading) {
    if (!degree || leading.zero()) throw std::invalid_argument("a root bound needs a proved degree of at least 1");
    // |a_n| from below; a bound needs no more than a word's precision.
    constexpr slong bound_bits = 64;
    Magnitude lead;
    {
        ComplexBall a;
        arb_set_fmpq(acb_realref(a.value()), leading.real.value(), bound_bits);
        arb_set_fmpq(acb_imagref(a.value()), leading.imag.value(), bound_bits);
        acb_get_mag_lower(lead.value(), a.value());
    }
    const auto majorant = [&](const mag_t r) {
        MajorantArithmetic arithmetic(r);
        return run_program(program, arithmetic);
    };
    // Cauchy: every root has |z| < 1 + max_{i<n} |a_i| / |a_n|, at most M(1) / |a_n| since M(1) bounds the sum of
    // the |a_i|, a_n's included.
    Magnitude unit, cauchy;
    mag_one(unit.value());
    const auto at_one = majorant(unit.value());
    mag_div(cauchy.value(), at_one.value(), lead.value());
    // The test at R = 2^e: M(R) < 2 |a_n| R^n, the left side rounded up, the right side down (a lower bound times a
    // power of two is one). It shows sum_{i<n} |a_i| R^i < |a_n| R^n, which then holds at every larger radius:
    // no root of modulus R or more.
    Integer exponent;
    const auto holds = [&](slong e) {
        Magnitude r, right;
        mag_one(r.value());
        mag_mul_2exp_si(r.value(), r.value(), e);
        const auto left = majorant(r.value());
        fmpz_set_si(exponent.value(), e);
        fmpz_mul_ui(exponent.value(), exponent.value(), degree);
        fmpz_add_ui(exponent.value(), exponent.value(), 1);
        mag_mul_2exp_fmpz(right.value(), lead.value(), exponent.value());
        return mag_cmp(left.value(), right.value()) < 0;
    };
    // Cauchy's bound is below 2^high. Where the test fails there (the program's leading terms cancel, so M grows
    // faster than R^n), Cauchy's bound stands.
    if (mag_is_zero(cauchy.value()) || !fmpz_fits_si(MAG_EXPREF(cauchy.value()))) return cauchy;
    slong high = fmpz_get_si(MAG_EXPREF(cauchy.value()));
    if (!holds(high)) return cauchy;
    slong low = high, step = 1;
    for (int k = 0; k < gallop_steps; ++k, step *= 2) {
        low = high - step;
        if (!holds(low)) break;
        high = low;
    }
    if (low < high) {
        while (high - low > 1) {
            const slong middle = low + (high - low) / 2;
            if (holds(middle)) high = middle;
            else low = middle;
        }
    }
    Magnitude power;
    mag_one(power.value());
    mag_mul_2exp_si(power.value(), power.value(), high);
    if (mag_cmp(power.value(), cauchy.value()) < 0) return power;
    return cauchy;
}

bool program_monomial(const Program &program, uint64_t degree, const GaussianRational &leading) {
    if (!degree || leading.zero()) return false;
    DenominatorBits denominators;
    const uint64_t d = run_program(program, denominators);
    if (d == saturated) return false;
    // |a_n| from below for the test, from above for L: |a_n| <= upper < 2^L, L the upper bound's exponent.
    constexpr slong bound_bits = 64;
    Magnitude lower, upper;
    {
        ComplexBall a;
        arb_set_fmpq(acb_realref(a.value()), leading.real.value(), bound_bits);
        arb_set_fmpq(acb_imagref(a.value()), leading.imag.value(), bound_bits);
        acb_get_mag_lower(lower.value(), a.value());
        acb_get_mag(upper.value(), a.value());
    }
    // R = 2^e, e = -(D + L + 1). a_n's own denominator divides Q, so |a_n| >= 2^-D and L + D >= 1: e <= -2, R <= 1/4.
    Integer e;
    fmpz_set_ui(e.value(), d);
    fmpz_add(e.value(), e.value(), MAG_EXPREF(upper.value()));
    fmpz_add_ui(e.value(), e.value(), 1);
    fmpz_neg(e.value(), e.value());
    Magnitude r;
    mag_one(r.value());
    mag_mul_2exp_fmpz(r.value(), r.value(), e.value());
    MajorantArithmetic majorant(r.value());
    const auto left = run_program(program, majorant);
    // M(R) < 2 |a_n| R^n, the left side rounded up and the right side a lower bound (lower times a power of two):
    // then sum_{i<n} |a_i| R^i < |a_n| R^n, so |a_i| < |a_n| R^(n-i) <= |a_n| R < 2^-D, and a_i = 0.
    Integer exponent;
    fmpz_mul_ui(exponent.value(), e.value(), degree);
    fmpz_add_ui(exponent.value(), exponent.value(), 1);
    Magnitude right;
    mag_mul_2exp_fmpz(right.value(), lower.value(), exponent.value());
    return mag_cmp(left.value(), right.value()) < 0;
}

uint64_t program_zero_order(const Program &program, uint64_t degree) {
    if (!degree) return 0;
    DenominatorBits denominators;
    const uint64_t d = run_program(program, denominators);
    if (d == saturated || d > zero_order_denominator_limit) return 0;
    Magnitude smallest; // 2^-D: no nonzero coefficient is smaller
    mag_one(smallest.value());
    mag_mul_2exp_si(smallest.value(), smallest.value(), -static_cast<slong>(d));
    for (size_t width = 1;; width = static_cast<size_t>(std::min<uint64_t>(2 * width, degree))) {
        if (width > top_width_limit) return 0;
        auto ceiling = static_cast<slong>(2 * (d + 64));
        for (slong precision = 64;; precision = std::min<slong>(2 * precision, ceiling)) {
            LowSeries arithmetic(width, precision);
            const auto low = run_program(program, arithmetic);
            bool undecided = false;
            for (size_t k = 0; k < width; ++k) {
                if (!acb_contains_zero(low[k].value())) return k; // the first nonzero, every one below it zero
                Magnitude size;
                acb_get_mag(size.value(), low[k].value());
                if (mag_cmp(size.value(), smallest.value()) < 0) continue; // within 2^-D of zero: zero
                undecided = true;
                if (precision == 64) { // the terms this coefficient cancels (zero_order_precision_limit's comment)
                    const double needed = mag_get_d_log2_approx(size.value()) + 64 + static_cast<double>(d) + 64;
                    const double limit = static_cast<double>(zero_order_precision_limit);
                    ceiling = std::max(ceiling, static_cast<slong>(std::isfinite(needed) ? std::min(needed, limit) : limit));
                }
                break;
            }
            if (!undecided) break;              // `width` zero coefficients: wider, or all of them
            if (precision >= ceiling) return 0; // undecided at the ceiling: no claim
        }
        // Every coefficient below the degree zero: p = a_n x^n (the leading coefficient is proved nonzero).
        if (width >= degree) return degree;
    }
}

ProgramPreparation prepare_program(const Program &program, const cps_options &options) {
    ProgramPreparation out;
    const auto refuse = [&](ProgramPreparation::Refusal refusal, std::string diagnostic) {
        out.refusal = refusal;
        out.diagnostic = std::move(diagnostic);
        return out;
    };
    if (!program.output || *program.output >= program.nodes.size())
        return refuse(ProgramPreparation::Refusal::NoOutput, "the program has no output");
    if (!options.exact_effort_fits(program.nodes.size()))
        return refuse(ProgramPreparation::Refusal::Length,
                      "the program has more operations than the exact operations allowed");
    ProgramTop top;
    try {
        top = top_of_program(program, top_width_limit, options.exact_storage_limit());
    } catch (const std::out_of_range &error) {
        return refuse(ProgramPreparation::Refusal::ConstantStorage, error.what());
    }
    switch (top.verdict) {
    case ProgramTop::Verdict::ZeroPolynomial:
        return refuse(ProgramPreparation::Refusal::ZeroPolynomial, "the program's polynomial is identically zero");
    case ProgramTop::Verdict::Undetermined:
        return refuse(ProgramPreparation::Refusal::LeadingUndetermined,
                      "leading coefficient not determined (the program's leading terms cancel over more than " +
                          std::to_string(top_width_limit - 1) + " coefficients)");
    case ProgramTop::Verdict::Determined: break;
    }
    out.degree = top.degree;
    out.top_width = top.width;
    out.leading = std::move(top.leading);
    if (!out.degree) return out; // a nonzero constant: no roots, nothing to bound or test
    out.modular = modular_squarefree_image(program, out.degree, options);
    // Every root exactly zero: the root bound stays zero (no search). Otherwise the order of a root at zero, which
    // proves a monomial too when it is the degree (a majorant too loose for program_monomial: (x + 1) - 1).
    out.monomial = program_monomial(program, out.degree, out.leading);
    if (!out.monomial) {
        out.zero_order = program_zero_order(program, out.degree);
        if (out.zero_order == out.degree) {
            out.monomial = true;
            out.zero_order = 0;
        }
    }
    if (!out.monomial) out.root_bound = program_root_bound(program, out.degree, out.leading);
    return out;
}

uint64_t program_solve_storage(uint64_t degree, uint64_t precision_bits) {
    // Per root, for the solve's life: the refinement's proposal, node and residue (SecularRefinement), complex
    // balls of two parts whose mantissas hold the working precision; its certified cover and the result's copy of
    // it (result.roots, solver_program.cpp), each a Root (a complex ball and four published endpoints, each an
    // Integer mantissa of the precision's bits, an mpz past the small-integer range, and a small exponent). The
    // packets' point vectors are of the same order and transient; the interaction tree is smaller. History: until
    // 26 October 2026 the refinement also held a restart copy of its cover (three covers, of which the count named
    // the refinement's two), and the ends were decimal strings, priced at a digit for every log2(10) bits. The ends are binary since, and so is their
    // price: at the precision cap the same
    // temporary storage admits programs of about 1.6 times the degree.
    constexpr uint64_t balls_per_root = 3, parts_per_ball = 2, covers_per_root = 2, endpoints_per_root = 4;
    constexpr uint64_t bits_per_byte = 8;
    const uint64_t mantissa_bytes = add(precision_bits / bits_per_byte, 1);
    const uint64_t ball_bytes = multiply(parts_per_ball, add(sizeof(arb_struct), mantissa_bytes));
    const uint64_t endpoint_bytes = add(sizeof(Endpoint), add(sizeof(__mpz_struct), mantissa_bytes));
    const uint64_t root_bytes = add(sizeof(Root), add(ball_bytes, multiply(endpoints_per_root, endpoint_bytes)));
    const uint64_t per_root = add(multiply(balls_per_root, ball_bytes), multiply(covers_per_root, root_bytes));
    return multiply(degree, per_root);
}

} // namespace certpolysolve::detail
