#include "solve/policies.hpp"
// The all-real seed pass (src/real_line_seeds.cpp): its gate never refuses a polynomial whose roots are all real and
// refuses the plain cases with non-real roots; its seeds sit next to roots, one a root, on inputs whose roots are not
// dyadic; the same bits every run; every stop is "no seeds", including a cancellation at each of its polls; and which
// queries it seeds through the public API (Auto on exact real coefficients, not a forced method).
#include "numeric/arithmetic/polynomial_arithmetic.hpp"
#include "core/heuristics.hpp"
#include "internal.hpp"
#include "solve/real_line_routing.hpp"
#include "numeric/candidates/real_line_seeds.hpp"
#include <certpolysolve/certpolysolve.h>
#include <flint/arb.h>
#include <flint/arb_fmpz_poly.h>
#include <flint/fmpq.h>
#include <flint/fmpq_poly.h>
#include <flint/fmpz_poly.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

using namespace certpolysolve::detail;

namespace {
int failures = 0;
void check(bool ok, const char *what) {
    if (!ok) {
        std::fprintf(stderr, "FAIL: %s\n", what);
        ++failures;
    }
}

// prod (d x - r_i) over integer pairs: roots r_i / d, not dyadic when d is odd and r_i not a multiple of d.
void product_of_linears(fmpz_poly_t out, const std::vector<long> &numerators, long denominator) {
    fmpz_poly_one(out);
    fmpz_poly_t factor;
    fmpz_poly_init(factor);
    for (const long r : numerators) {
        fmpz_poly_zero(factor);
        fmpz_poly_set_coeff_si(factor, 1, denominator);
        fmpz_poly_set_coeff_si(factor, 0, -r);
        fmpz_poly_mul(out, out, factor);
    }
    fmpz_poly_clear(factor);
}

// The sign of p at x (a binary64), certain or 0.
int sign_at(const fmpz_poly_t p, double x) {
    arb_t point, value;
    arb_init(point);
    arb_init(value);
    arb_set_d(point, x);
    int s = 0;
    for (slong prec = 128; prec <= 1 << 16 && s == 0; prec *= 2) {
        arb_fmpz_poly_evaluate_arb(value, p, point, prec);
        s = arb_is_positive(value) ? 1 : (arb_is_negative(value) ? -1 : 0);
    }
    arb_clear(point);
    arb_clear(value);
    return s;
}

// Every seed has a sign change of p within s (1 +- 2^-10), and the seeds are increasing and at most deg p.
bool seeds_bracket_roots(const fmpz_poly_t p, const std::vector<double> &seeds) {
    if (static_cast<slong>(seeds.size()) > fmpz_poly_degree(p)) return false;
    for (size_t i = 0; i < seeds.size(); ++i) {
        if (i > 0 && !(seeds[i] > seeds[i - 1])) return false;
        const double s = seeds[i], d = std::ldexp(std::fabs(s), -10);
        if (sign_at(p, s - d) * sign_at(p, s + d) >= 0) return false;
    }
    return true;
}

void test_newton() {
    fmpz_poly_t p;
    fmpz_poly_init(p);
    fmpz_poly_chebyshev_t(p, 40);
    // T_40 has zero odd coefficients: the zero rule applies at every odd k, and holds (its neighbours alternate).
    check(newton_inequality_failure(p) == -1, "Newton: T_40 passes");
    product_of_linears(p, {-7, -2, 1, 3, 4, 11, 13, 20}, 3);
    check(newton_inequality_failure(p) == -1, "Newton: a product of linear factors passes");
    fmpz_poly_zero(p);
    fmpz_poly_set_coeff_si(p, 2, 1);
    fmpz_poly_set_coeff_si(p, 0, 1); // x^2 + 1
    check(newton_inequality_failure(p) == 1, "Newton: x^2 + 1 fails at k = 1");
    fmpz_poly_zero(p);
    fmpz_poly_set_coeff_si(p, 3, 1);
    fmpz_poly_set_coeff_si(p, 0, 1); // x^3 + 1: Newton's inequality holds as 0 >= 0, the zero rule refuses it
    check(newton_inequality_failure(p) == 1, "Newton: x^3 + 1 fails by the zero rule");
    // A polynomial with random coefficients of random signs (fixed seed): fails, as every one in the corpus did.
    fmpz_poly_zero(p);
    unsigned state = 12345;
    for (int k = 0; k <= 60; ++k) {
        state = state * 1103515245u + 12345u;
        const long c = static_cast<long>((state >> 8) % 2001) - 1000;
        fmpz_poly_set_coeff_si(p, k, c == 0 ? 1 : c);
    }
    check(newton_inequality_failure(p) >= 0, "Newton: a random polynomial fails");
    fmpz_poly_clear(p);
}

// The check as it was before its leading-bits estimate (exact products at every k): the estimate must decide as it does.
slong newton_reference(const fmpz_poly_struct *p) {
    const slong n = fmpz_poly_degree(p);
    fmpz_t left, right;
    fmpz_init(left);
    fmpz_init(right);
    slong failure = -1;
    for (slong k = 1; k < n && failure < 0; ++k) {
        const fmpz *before = p->coeffs + k - 1, *at = p->coeffs + k, *after = p->coeffs + k + 1;
        if (fmpz_is_zero(at)) {
            if (fmpz_sgn(before) * fmpz_sgn(after) >= 0) failure = k;
            continue;
        }
        fmpz_mul(left, at, at);
        fmpz_mul_si(left, left, k * (n - k));
        fmpz_mul(right, before, after);
        fmpz_mul_si(right, right, (k + 1) * (n - k + 1));
        if (fmpz_cmp(left, right) < 0) failure = k;
    }
    fmpz_clear(left);
    fmpz_clear(right);
    return failure;
}

void test_newton_estimate() {
    flint_rand_t state;
    flint_rand_init(state);
    fmpz_poly_t p;
    fmpz_poly_init(p);
    fmpz_t t, a, scale;
    fmpz_init(t);
    fmpz_init(a);
    fmpz_init(scale);
    int disagreements = 0, ties = 0, failures_seen = 0;
    for (int trial = 0; trial < 3000; ++trial) {
        const slong n = 3 + static_cast<slong>(n_randint(state, 60));
        const flint_bitcnt_t bits = 1 + n_randint(state, trial % 3 == 0 ? 3000 : 80);
        fmpz_poly_zero(p);
        // A real-rooted product of linear factors (every k holds), a random polynomial, or either with one
        // coefficient set to an exact tie of Newton's inequality, nudged by -1, 0 or +1.
        if (trial % 2 == 0) {
            fmpz_poly_one(p);
            fmpz_poly_t linear;
            fmpz_poly_init(linear);
            for (slong i = 0; i < n; ++i) {
                fmpz_poly_zero(linear);
                fmpz_randtest_not_zero(a, state, 1 + bits / static_cast<flint_bitcnt_t>(n));
                fmpz_poly_set_coeff_fmpz(linear, 0, a);
                fmpz_poly_set_coeff_si(linear, 1, 1 + static_cast<slong>(n_randint(state, 3)));
                fmpz_poly_mul(p, p, linear);
            }
            fmpz_poly_clear(linear);
        } else {
            for (slong i = 0; i <= n; ++i) {
                fmpz_randtest_not_zero(a, state, bits);
                fmpz_poly_set_coeff_fmpz(p, i, a);
            }
        }
        const slong m = fmpz_poly_degree(p);
        if (trial % 3 == 1 && m >= 2) {
            // a_(k-1) = k (m - k) t^2 and a_(k+1) = (k + 1) (m - k + 1), so a_k = (k + 1) (m - k + 1) t is a tie.
            const slong k = 1 + static_cast<slong>(n_randint(state, static_cast<ulong>(m - 1)));
            fmpz_randtest_not_zero(t, state, bits);
            fmpz_mul(a, t, t);
            fmpz_mul_si(a, a, k * (m - k));
            fmpz_poly_set_coeff_fmpz(p, k - 1, a);
            fmpz_poly_set_coeff_si(p, k + 1, (k + 1) * (m - k + 1));
            fmpz_mul_si(a, t, (k + 1) * (m - k + 1));
            fmpz_add_si(a, a, static_cast<slong>(n_randint(state, 3)) - 1);
            if (fmpz_is_zero(a)) fmpz_one(a);
            fmpz_poly_set_coeff_fmpz(p, k, a);
            ++ties;
        }
        const slong reference = newton_reference(p);
        if (newton_inequality_failure(p) != reference) ++disagreements;
        if (reference >= 0) ++failures_seen;
    }
    check(disagreements == 0, "Newton: the leading-bits estimate decides as the exact products do");
    check(ties > 500 && failures_seen > 500, "Newton: the trials include ties and failures");
    fmpz_clear(t);
    fmpz_clear(a);
    fmpz_clear(scale);
    fmpz_poly_clear(p);
    flint_rand_clear(state);
}

void test_transform() {
    // dct1 against the direct sum in arb, M = 64.
    const slong m = 64;
    Binary64Transform t(2 * m);
    std::vector<double> x(m + 1);
    for (slong k = 0; k <= m; ++k) x[static_cast<size_t>(k)] = std::ldexp(static_cast<double>((k * 37) % 23) - 11.0, -3);
    const auto y = t.dct1(x);
    double worst = 0;
    arb_t sum, term, c;
    arb_init(sum);
    arb_init(term);
    arb_init(c);
    fmpq_t angle;
    fmpq_init(angle);
    for (slong j = 0; j <= m; ++j) {
        arb_set_d(sum, x[0] + ((j % 2) ? -x[static_cast<size_t>(m)] : x[static_cast<size_t>(m)]));
        for (slong k = 1; k < m; ++k) {
            fmpq_set_si(angle, j * k, static_cast<ulong>(m));
            arb_cos_pi_fmpq(c, angle, 128);
            arb_set_d(term, 2 * x[static_cast<size_t>(k)]);
            arb_mul(term, term, c, 128);
            arb_add(sum, sum, term, 128);
        }
        worst = std::max(worst, std::fabs(arf_get_d(arb_midref(sum), ARF_RND_NEAR) - y[static_cast<size_t>(j)]));
    }
    arb_clear(sum);
    arb_clear(term);
    arb_clear(c);
    fmpq_clear(angle);
    check(worst < 1e-12, "transform: dct1 equals the direct sum");
    // Coefficients and values round trip; Clenshaw equals the values at the nodes.
    const auto coefficients = chebyshev_coefficients(t, x);
    const auto back = chebyshev_values(t, coefficients, m);
    double trip = 0, clenshaw = 0;
    for (slong j = 0; j <= m; ++j) {
        trip = std::max(trip, std::fabs(back[static_cast<size_t>(j)] - x[static_cast<size_t>(j)]));
        clenshaw = std::max(clenshaw, std::fabs(chebyshev_evaluate(coefficients, m + 1, t.cos_pi(j, m)) - x[static_cast<size_t>(j)]));
    }
    check(trip < 1e-12, "transform: coefficients and values round trip");
    check(clenshaw < 1e-12, "transform: Clenshaw at the nodes");
}

void test_seeds() {
    fmpz_poly_t p;
    fmpz_poly_init(p);
    RealLineSeedLimits limits;
    // Chebyshev T_n in x^2 (the solver's power composition hands the pass q(y), T_n(x) = q(x^2)): n/2 roots.
    fmpz_poly_chebyshev_t(p, 96);
    {
        fmpz_poly_t q;
        fmpz_poly_init(q);
        for (slong k = 0; k <= 48; ++k) fmpz_poly_set_coeff_fmpz(q, k, p->coeffs + 2 * k);
        const auto r = real_line_seeds(q, limits);
        check(r.statistics.stop == RealLineSeedStop::Seeds, "seeds: T_96 in x^2 is seeded");
        check(static_cast<slong>(r.seeds.size()) == 48, "seeds: T_96 in x^2 gets 48 seeds");
        check(seeds_bracket_roots(q, r.seeds), "seeds: T_96 in x^2's seeds bracket roots");
        fmpz_poly_clear(q);
    }
    // Roots r/7 (not dyadic), including a close pair 1/7 apart among wider gaps.
    std::vector<long> numerators;
    for (long r = -59; r <= 58; r += 3) numerators.push_back(r); // none zero: the pass needs p(0) != 0
    numerators.push_back(60); // 2/7 from 58/7: the close pair
    product_of_linears(p, numerators, 7);
    const auto r = real_line_seeds(p, limits);
    check(r.statistics.stop == RealLineSeedStop::Seeds, "seeds: the product of linear factors is seeded");
    check(static_cast<slong>(r.seeds.size()) == fmpz_poly_degree(p), "seeds: one seed a root");
    check(seeds_bracket_roots(p, r.seeds), "seeds: the product's seeds bracket roots");
    // A root at zero (as a factor centred on one of its roots has): divided out, and 0 itself a seed.
    {
        fmpz_poly_t z;
        fmpz_poly_init(z);
        fmpz_poly_shift_left(z, p, 1);
        const auto rz = real_line_seeds(z, limits);
        check(rz.statistics.stop == RealLineSeedStop::Seeds, "seeds: a root at zero is seeded");
        check(static_cast<slong>(rz.seeds.size()) == fmpz_poly_degree(z), "seeds: a root at zero gets every seed");
        check(std::binary_search(rz.seeds.begin(), rz.seeds.end(), 0.0), "seeds: 0 is a seed");
        fmpz_poly_clear(z);
    }
    // The same bits every run.
    const auto again = real_line_seeds(p, limits);
    check(again.seeds.size() == r.seeds.size() &&
              std::memcmp(again.seeds.data(), r.seeds.data(), r.seeds.size() * sizeof(double)) == 0,
          "seeds: the same bits every run");
    fmpz_poly_clear(p);
}

void test_refusals() {
    fmpz_poly_t p, f;
    fmpz_poly_init(p);
    fmpz_poly_init(f);
    RealLineSeedLimits limits;
    // (x^2 + 1) times real linear factors: a complex pair. Whatever stops it, it is not seeded as if all real.
    product_of_linears(p, {-5, -1, 2, 3, 8}, 3);
    fmpz_poly_zero(f);
    fmpz_poly_set_coeff_si(f, 2, 1);
    fmpz_poly_set_coeff_si(f, 0, 1);
    fmpz_poly_mul(p, p, f);
    {
        const auto r = real_line_seeds(p, limits);
        check(r.statistics.stop != RealLineSeedStop::Seeds || static_cast<slong>(r.seeds.size()) < fmpz_poly_degree(p),
              "refusal: a complex pair is not seeded as all real");
        check(seeds_bracket_roots(p, r.seeds), "refusal: any seeds it gives bracket roots");
    }
    // f^2 + 1 with f real-rooted: no real root at all.
    product_of_linears(f, {-3, 1, 4, 9}, 5);
    fmpz_poly_mul(p, f, f);
    fmpz_poly_set_coeff_si(p, 0, fmpz_get_si(p->coeffs) + 1);
    {
        const auto r = real_line_seeds(p, limits);
        check(r.seeds.empty(), "refusal: f^2 + 1 gets no seeds");
    }
    // Degree 1 and a zero constant term.
    fmpz_poly_zero(p);
    fmpz_poly_set_coeff_si(p, 1, 3);
    fmpz_poly_set_coeff_si(p, 0, -1);
    check(real_line_seeds(p, limits).statistics.stop == RealLineSeedStop::Degree, "refusal: degree 1");
    fmpz_poly_shift_left(p, p, 1); // x (3x - 1): divided by x, a linear factor is left
    check(real_line_seeds(p, limits).statistics.stop == RealLineSeedStop::Degree, "refusal: degree 1 after the root at zero");
    fmpz_poly_clear(p);
    fmpz_poly_clear(f);
}

void test_limits() {
    fmpz_poly_t p;
    fmpz_poly_init(p);
    std::vector<long> numerators;
    for (long r = -41; r <= 39; r += 4) numerators.push_back(r);
    product_of_linears(p, numerators, 5);
    // A cancellation at every poll: each run stops there with no seeds; past the last poll the pass completes.
    RealLineSeedLimits counting;
    long polls = 0;
    counting.cancelled = [&] { ++polls; return false; };
    const auto full = real_line_seeds(p, counting);
    check(full.statistics.stop == RealLineSeedStop::Seeds && polls > 3, "limits: the pass polls");
    for (long k = 1; k <= polls; ++k) {
        RealLineSeedLimits stopping;
        long seen = 0;
        stopping.cancelled = [&] { return ++seen >= k; };
        const auto r = real_line_seeds(p, stopping);
        if (r.statistics.stop != RealLineSeedStop::Cancelled || !r.seeds.empty()) {
            std::fprintf(stderr, "cancellation at poll %ld of %ld: stop %s, %zu seeds\n", k, polls,
                         real_line_seed_stop_name(r.statistics.stop), r.seeds.size());
            check(false, "limits: a cancellation at every poll stops the pass with no seeds");
            break;
        }
    }
    RealLineSeedLimits narrow;
    narrow.max_precision = 64;
    check(real_line_seeds(p, narrow).statistics.stop == RealLineSeedStop::Precision, "limits: the precision cap");
    RealLineSeedLimits small;
    small.max_bytes = 1024;
    check(real_line_seeds(p, small).statistics.stop == RealLineSeedStop::Memory, "limits: the memory cap");
    RealLineSeedLimits fraction;
    fraction.min_cell_fraction = 2.0; // more sign changes than the degree can give
    check(real_line_seeds(p, fraction).statistics.stop == RealLineSeedStop::FewCells, "limits: the cell fraction");
    fmpz_poly_clear(p);
}
} // namespace

// Which queries the pass seeds, through the public API: Auto on exact real coefficients, for all roots (the default
// and the approximation goal) and for the real roots; a forced method keeps the circles. The input: 64! L_64, the
// Laguerre polynomial (irreducible, its roots real, positive and not dyadic), whose one factor has the pass's floor
// degree, 64.
std::string solve_diagnostics(const cps_polynomial *p, int method, int goal, int real_only, size_t *regions,
                              bool *complete) {
    cps_options *o = nullptr;
    cps_result *r = nullptr;
    std::string diagnostics;
    if (cps_options_create(&o) == CPS_OK && cps_options_set_method(o, method) == CPS_OK &&
        cps_options_set_goal(o, goal) == CPS_OK && cps_options_set_real_only(o, real_only) == CPS_OK &&
        cps_solve(p, o, &r) == CPS_OK) {
        diagnostics = cps_result_diagnostics(r);
        *regions = cps_result_size(r);
        *complete = cps_result_complete(r) && cps_result_goal_met(r);
    }
    cps_result_destroy(r);
    cps_options_destroy(o);
    return diagnostics;
}
// The Laguerre polynomial L_n with its denominators cleared: irreducible, its roots real, positive, not dyadic.
void laguerre(fmpz_poly_t out, slong n) {
    fmpq_poly_t l;
    fmpq_poly_init(l);
    fmpq_poly_laguerre_l(l, static_cast<ulong>(n));
    fmpq_poly_get_numerator(out, l);
    fmpq_poly_clear(l);
}
// p as the public input object (the caller destroys it).
cps_polynomial *public_polynomial(const fmpz_poly_t p) {
    cps_polynomial *out = nullptr;
    check(cps_polynomial_create(&out) == CPS_OK, "the public polynomial");
    for (slong k = 0; k <= fmpz_poly_degree(p); ++k) {
        char *text = fmpz_get_str(nullptr, 10, p->coeffs + k);
        check(cps_polynomial_set_coefficient(out, static_cast<uint64_t>(k), text, "0") == CPS_OK, "a coefficient");
        flint_free(text);
    }
    return out;
}
void test_routes() {
    constexpr slong degree = 64;
    static_assert(degree >= heuristics::real_line_seed_min_degree, "the input must reach the pass's floor");
    fmpz_poly_t integer;
    fmpz_poly_init(integer);
    laguerre(integer, degree);
    cps_polynomial *p = public_polynomial(integer);
    const std::string seeded = "real-line seeds: 64 of 64 (seeds)";
    size_t regions = 0;
    bool complete = false;
    for (const int goal : {CPS_GOAL_DEFAULT, CPS_ROOT_APPROXIMATION, CPS_ROOT_ACCURACY}) {
        regions = 0;
        complete = false;
        const auto d = solve_diagnostics(p, CPS_METHOD_AUTO, goal, 0, &regions, &complete);
        check(d.find(seeded) != std::string::npos, "routes: Auto's query for all roots is seeded");
        check(complete && regions == degree, "routes: Auto's query for all roots finds the 64 roots");
    }
    regions = 0;
    complete = false;
    const auto real = solve_diagnostics(p, CPS_METHOD_AUTO, CPS_GOAL_DEFAULT, 1, &regions, &complete);
    check(real.find(seeded) != std::string::npos, "routes: Auto's real-only query is seeded");
    check(complete && regions == degree, "routes: Auto's real-only query finds the 64 roots");
    regions = 0;
    complete = false;
    const auto forced = solve_diagnostics(p, CPS_METHOD_SECULAR, CPS_ROOT_ACCURACY, 0, &regions, &complete);
    check(forced.find("real-line seeds") == std::string::npos, "routes: forced secular starts from the circles");
    check(complete && regions == degree, "routes: forced secular finds the 64 roots");
    cps_polynomial_destroy(p);
    fmpz_poly_clear(integer);
}

// The check before Auto routes a query for all roots (real_line_routing.hpp): a hint on L_64 (its primitive
// polynomial and its 64 seeds); none on (x^2 + 1) L_62, two of whose roots are not real, nor below the pass's floor.
// And the hand-over: a solve given a hint takes its seeds (here only the first 10 of L_64's, so that the note shows
// whose they are) instead of running the pass, and still finds every root.
void test_hint() {
    constexpr slong degree = 64;
    fmpz_poly_t l64, product, l32;
    fmpz_poly_init(l64);
    fmpz_poly_init(product);
    fmpz_poly_init(l32);
    laguerre(l64, degree);
    laguerre(product, degree - 2);
    {
        fmpz_poly_t square_plus_one;
        fmpz_poly_init(square_plus_one);
        fmpz_poly_set_coeff_si(square_plus_one, 2, 1);
        fmpz_poly_set_coeff_si(square_plus_one, 0, 1);
        fmpz_poly_mul(product, product, square_plus_one);
        fmpz_poly_clear(square_plus_one);
    }
    laguerre(l32, 32);
    cps_options o;
    o.bits = 53;
    cps_polynomial *p = public_polynomial(l64);
    const DensePolynomial input = dense(*p);
    auto hint = every_root_seeded(input.coefficients, o);
    fmpz_poly_t primitive;
    fmpz_poly_init(primitive);
    fmpz_poly_primitive_part(primitive, l64);
    check(hint && static_cast<slong>(hint->pass.seeds.size()) == degree, "hint: L_64 has every root seeded");
    check(hint && fmpz_poly_equal(hint->polynomial.value(), primitive), "hint: on L_64's primitive polynomial");
    for (const fmpz_poly_struct *other : {&product[0], &l32[0]}) {
        cps_polynomial *q = public_polynomial(other);
        check(!every_root_seeded(dense(*q).coefficients, o), "hint: none where a root is not real or below the floor");
        cps_polynomial_destroy(q);
    }
    if (hint) {
        constexpr size_t kept = 10;
        hint->pass.seeds.resize(kept);
        SolveState state(input, maintained_method_policy(CPS_METHOD_SECULAR));
        state.hint_real_line_seeds(std::move(*hint));
        cps_options solve = o;
        solve.goal = CPS_ROOT_ACCURACY;
        solve.real_line_seeds = true;
        const auto r = state.advance(&solve);
        check(r.diagnostics.find("real-line seeds: 10 of 64 (seeds)") != std::string::npos,
              "hint: the solve takes the hint's seeds instead of running the pass");
        check(r.complete && r.goal_met && r.roots.size() == static_cast<size_t>(degree),
              "hint: the hinted solve finds the 64 roots");
    }
    cps_polynomial_destroy(p);
    fmpz_poly_clear(primitive);
    fmpz_poly_clear(l32);
    fmpz_poly_clear(product);
    fmpz_poly_clear(l64);
}

int main() {
    test_newton();
    test_newton_estimate();
    test_transform();
    test_seeds();
    test_refusals();
    test_limits();
    test_routes();
    test_hint();
    if (failures) {
        std::fprintf(stderr, "%d failures\n", failures);
        return 1;
    }
    std::printf("real_line_seeds: all checks passed\n");
    return 0;
}
