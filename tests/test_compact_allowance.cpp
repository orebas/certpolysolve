// Exact preparation of m@e coefficients under the compact-exponent rule (docs/resource-policy.md, "Compact exponent
// rule"): the written-out numbers may exceed the input's own storage by one precision cap a number, and where the
// write-out is still refused, the modular test runs on the coefficients as written. The regressions come from the
// easy-input panel's finding (dev/docs/easy-inputs-profile-20261025.md, section 4): an exact double root among
// binary64 roots at degree 48, real and complex. The rule's room was the size of the number records alone, so it
// refused them, and the accuracy goal climbed to the precision cap and stopped.
#include "core/number.hpp"
#include "core/raii.hpp"
#include <certpolysolve/certpolysolve.hpp>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <random>
#include <set>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace d = certpolysolve::detail;
using namespace std::chrono_literals;

namespace {
void require(bool ok, const char *why) {
    if (!ok) throw std::runtime_error(why);
}

// The degree of the regression inputs: past the old rule's reach (degree 40 already failed), small enough to stay
// a quick test.
constexpr int regression_degree = 48;
// A guard against a regression that climbs to the cap again: the fixed solves take milliseconds.
constexpr auto solve_limit = 60s;

// A dyadic rational, written as the m@e literal the panel's inputs use (an odd mantissa times a power of two).
std::string literal(const fmpq *q) {
    if (fmpq_is_zero(q)) return "0";
    const slong denominator_bits = static_cast<slong>(fmpz_val2(fmpq_denref(q)));
    d::Integer power;
    fmpz_one(power.value());
    fmpz_mul_2exp(power.value(), power.value(), static_cast<ulong>(denominator_bits));
    require(fmpz_equal(power.value(), fmpq_denref(q)), "a coefficient is not dyadic");
    d::Integer mantissa;
    const slong shift = static_cast<slong>(fmpz_val2(fmpq_numref(q)));
    fmpz_tdiv_q_2exp(mantissa.value(), fmpq_numref(q), static_cast<ulong>(shift));
    char *digits = fmpz_get_str(nullptr, 10, mantissa.value());
    std::string out = std::string(digits) + "@" + std::to_string(shift - denominator_bits);
    flint_free(digits);
    return out;
}

// A binary64-like dyadic in (-1, 1): a 53-bit odd mantissa over 2^53 to 2^60.
d::Rational dyadic(std::mt19937_64 &rng) {
    d::Integer mantissa, denominator;
    fmpz_set_ui(mantissa.value(), (rng() >> 11) | (UINT64_C(1) << 52) | 1);
    if (rng() & 1) fmpz_neg(mantissa.value(), mantissa.value());
    fmpz_one(denominator.value());
    fmpz_mul_2exp(denominator.value(), denominator.value(), 53 + rng() % 8);
    d::Rational out;
    fmpq_set_fmpz_frac(out.value(), mantissa.value(), denominator.value());
    return out;
}

struct Gaussian {
    d::Rational real, imag;
};

// The coefficients of the product of (x - r) over `roots`, lowest first, exactly.
std::vector<Gaussian> product(const std::vector<Gaussian> &roots) {
    std::vector<Gaussian> c(1);
    fmpq_one(c[0].real.value());
    d::Rational t;
    for (const auto &r : roots) {
        std::vector<Gaussian> next(c.size() + 1);
        for (size_t i = 0; i < c.size(); ++i) {
            fmpq_add(next[i + 1].real.value(), next[i + 1].real.value(), c[i].real.value());
            fmpq_add(next[i + 1].imag.value(), next[i + 1].imag.value(), c[i].imag.value());
            // next[i] -= r c[i]
            fmpq_mul(t.value(), r.real.value(), c[i].real.value());
            fmpq_sub(next[i].real.value(), next[i].real.value(), t.value());
            fmpq_mul(t.value(), r.imag.value(), c[i].imag.value());
            fmpq_add(next[i].real.value(), next[i].real.value(), t.value());
            fmpq_mul(t.value(), r.real.value(), c[i].imag.value());
            fmpq_sub(next[i].imag.value(), next[i].imag.value(), t.value());
            fmpq_mul(t.value(), r.imag.value(), c[i].real.value());
            fmpq_sub(next[i].imag.value(), next[i].imag.value(), t.value());
        }
        c = std::move(next);
    }
    return c;
}

struct Input {
    certpolysolve::Polynomial polynomial;
    std::vector<d::Coefficient> coefficients;  // the same, as exact preparation reads them
};

Input written(const std::vector<Gaussian> &coefficients) {
    Input out;
    for (size_t i = 0; i < coefficients.size(); ++i) {
        const auto re = literal(coefficients[i].real.value()), im = literal(coefficients[i].imag.value());
        out.polynomial.set(i, certpolysolve::Number::complex(certpolysolve::Number(re), certpolysolve::Number(im)));
        d::Coefficient c;
        c.real = d::Number(re);
        c.imag = d::Number(im);
        out.coefficients.push_back(std::move(c));
    }
    return out;
}

// A double root and regression_degree - 2 simple roots, all distinct dyadics; real roots, or complex ones in no
// conjugate pairs (so that the coefficients are complex and the Q(i) path answers).
Input double_root_input(bool complex, uint64_t seed) {
    std::mt19937_64 rng(seed);
    std::vector<Gaussian> roots;
    std::set<std::string> seen;
    while (roots.size() + 1 < regression_degree) {
        Gaussian r;
        r.real = dyadic(rng);
        if (complex) r.imag = dyadic(rng);
        if (seen.insert(literal(r.real.value()) + " " + literal(r.imag.value())).second) roots.push_back(std::move(r));
    }
    roots.push_back(roots.front());  // the double root
    return written(product(roots));
}

void rule() {
    // x^2 + 2^e: the expanded constant grows by about e bits. The twelve number records own 3,098 bits plus the
    // literals' own bits (32 bytes a record), so 2^10000 outgrows them and 2^1000 does not.
    const auto list = [](const char *constant) {
        std::vector<d::Coefficient> c(3);
        c[0].real = d::Number(constant);
        c[2].real = d::Number("1");
        return c;
    };
    const auto small = list("1@10000"), large = list("1@1000000");
    require(d::compact_expansion_fits(list("1@1000"), 0), "an expansion within the owned records refused");
    require(!d::compact_expansion_fits(small, 0), "growth admitted with no allowance");
    require(d::compact_expansion_fits(small, CPS_DEFAULT_MAX_WORK_BITS), "growth within the allowance refused");
    require(!d::compact_expansion_fits(large, CPS_DEFAULT_MAX_WORK_BITS), "growth past the allowance admitted");
    // x^2 + 2^1000000 grows by 1,000,004 - 3,106 = 996,898 bits, and only its two nonzero numbers carry the
    // allowance: it fits from 498,449 bits a number. Counting the ten zero records as well would admit it at 400,000.
    require(!d::compact_expansion_fits(large, 498448), "the allowance boundary moved");
    require(d::compact_expansion_fits(large, 498449), "the allowance of the two written-out numbers refused");
    require(!d::compact_expansion_fits(large, 400000), "zero numbers counted in the allowance");
    // An exponent beyond 64 bits has no expanded size: never admitted.
    const auto huge = list("1@100000000000000000000000");
    require(!d::compact_expansion_fits(huge, UINT64_MAX), "an exponent beyond 64 bits admitted");
    // The boundary is inclusive, and the allowance saturates.
    require(d::expansion_fits(10, 16, 2, 3) && !d::expansion_fits(10, 17, 2, 3), "allowance boundary moved");
    require(d::expansion_fits(10, UINT64_MAX - 1, 2, UINT64_MAX), "saturated allowance refused");
    require(!d::expansion_fits(10, UINT64_MAX, 2, UINT64_MAX), "an unbounded expansion admitted");
}

void double_root(bool complex) {
    const auto why = complex ? "complex" : "real";
    const auto input = double_root_input(complex, complex ? 20261026 : 20261025);
    // The regime of the finding: the old rule refused this input, the allowance admits it.
    require(!d::compact_expansion_fits(input.coefficients, 0), "the input is not one the old rule refused");
    require(d::compact_expansion_fits(input.coefficients, CPS_DEFAULT_MAX_WORK_BITS), "the allowance refuses it");
    certpolysolve::Options options;
    options.goal = certpolysolve::Goal::accuracy();
    options.time_limit = solve_limit;
    const auto result = certpolysolve::solve(input.polynomial, options);
    if (!result.goal_met()) std::cerr << why << ": " << result.message() << '\n';
    require(result.goal_met(), "the accuracy goal is not met on an exact double root");
    require(result.squarefree() == certpolysolve::Squarefree::disproved, "the double root is not proved");
    uint64_t doubles = 0, total = 0;
    for (const auto &root : result.roots()) {
        require(root.multiplicity.has_value(), "a multiplicity is unknown");
        require(*root.multiplicity == root.count, "a region holds more than one distinct root");
        doubles += *root.multiplicity == 2;
        require(*root.multiplicity <= 2, "a multiplicity above two");
        total += root.count;
    }
    require(doubles == 1, "not exactly one double root");
    require(total == regression_degree, "roots lost");
    if (complex) return;
    // The real-only queries failed alike (a double real root cannot be proved real without its multiplicity): every
    // root here is real.
    for (const auto &goal : {certpolysolve::Goal::accuracy(), certpolysolve::Goal::isolation()}) {
        options.goal = goal;
        const auto real = certpolysolve::real_roots(input.polynomial, options);
        if (!real.goal_met()) std::cerr << "real-only: " << real.message() << '\n';
        require(real.goal_met(), "a real-only goal is not met on an exact double root");
        uint64_t real_doubles = 0, real_total = 0;
        for (const auto &root : real.roots()) {
            require(root.multiplicity.has_value() && *root.multiplicity == root.count, "a real multiplicity is unknown");
            real_doubles += *root.multiplicity == 2;
            real_total += root.count;
        }
        require(real_doubles == 1 && real_total == regression_degree, "the real roots are not all found once");
    }
}

void beyond_allowance() {
    // x^2 - 2^-200000: square-free, its write-out past the allowance (about 200,000 bits of growth against 2 x 65,536).
    certpolysolve::Polynomial square_free;
    square_free.set(0, certpolysolve::Number("-1@-200000")).set(2, certpolysolve::Number("1"));
    std::vector<d::Coefficient> c(3);
    c[0].real = d::Number("-1@-200000");
    c[2].real = d::Number("1");
    require(!d::compact_expansion_fits(c, CPS_DEFAULT_MAX_WORK_BITS), "the input is within the allowance");
    certpolysolve::Options options;
    options.goal = certpolysolve::Goal::accuracy();
    options.time_limit = solve_limit;
    const auto proved = certpolysolve::solve(square_free, options);
    require(proved.goal_met(), "the square-free input is not answered");
    require(proved.squarefree() == certpolysolve::Squarefree::proved, "the test on the coefficients as written is lost");
    require(!proved.preparation_limited(), "a proved input reported as limited");
    require(proved.message().find("squarefreeness proved modulo 998244353 on the coefficients as written") !=
                std::string::npos,
            "the proof is not in the diagnostics");
    // (x - 2^-100000)^2: refused, and the test fails, as it must on a double root. Nothing is claimed: the
    // approximation goal answers with a counted cluster and the multiplicity stays unknown.
    certpolysolve::Polynomial doubled;
    doubled.set(0, certpolysolve::Number("1@-200000")).set(1, certpolysolve::Number("-1@-99999")).set(2, certpolysolve::Number("1"));
    options.goal = certpolysolve::Goal::approximation();
    const auto cluster = certpolysolve::solve(doubled, options);
    if (!cluster.goal_met()) std::cerr << "counted compact double root: " << cluster.message() << '\n';
    require(cluster.goal_met(), "the double root's cluster is not answered");
    require(cluster.squarefree() == certpolysolve::Squarefree::unknown, "a verdict claimed without exact preparation");
    require(cluster.roots().size() == 1 && cluster.roots().front().count == 2, "the cluster's count is lost");
}
} // namespace

int main() {
    try {
        rule();
        double_root(false);
        double_root(true);
        beyond_allowance();
        std::cout << "compact-exponent allowance: the rule, a double root at degree " << regression_degree
                  << " (real and complex), and the test on the coefficients as written\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << "test_compact_allowance: " << e.what() << '\n';
        return 1;
    }
}
