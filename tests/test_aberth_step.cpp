// Aberth's step: the shared operations of src/aberth_step.hpp against the code they replaced at the kernels, bit for
// bit (W8's S10). Each reference function below is a kernel's code as it was at e454afc, verbatim but for its
// variables passed as arguments, and named by its file and lines there. The test runs a reference and the shared
// function on the same inputs and requires the same results: in binary64 and ScaledDouble the same bits (NaN payloads
// included), the same refusal and the same floating-point flags; in FLINT's acf equal values (arf_equal on both
// parts, which compares every bit of the mantissa and the exponent), the same returned flag and count, and the same
// scratch. Built with the library's floating-point flags (-fno-fast-math -ffp-contract=off; /fp:strict), so that no
// reference is contracted.
#include "numeric/candidates/aberth_step.hpp"
#include "numeric/arithmetic/block_double.hpp"
#include "numeric/arithmetic/machine_arithmetic.hpp"
#include "program_draws.hpp"
#include "numeric/evaluate/program_machine.hpp"
#include <algorithm>
#include <bit>
#include <cfenv>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <flint/acf.h>
#include <flint/arb.h>
#include <flint/arf.h>
#include <flint/flint.h>
#include <limits>
#include <memory>
#include <numeric>
#include <random>
#include <type_traits>
#include <utility>
#include <vector>

namespace {
namespace d = certpolysolve::detail;

int failures = 0;
size_t compared = 0;
// The refusals met: a point equal to z in the interaction loops, 1 - N S zero in the corrections.
size_t coincident = 0, zero_denominators = 0;
void require(bool ok, const char *what, size_t trial) {
    if (!ok && failures++ < 20)
        std::fprintf(stderr, "FAIL: %s (case %zu)\n", what, trial);
}

// A fixed sequence (SplitMix64), the same on every run.
uint64_t state = 0x9e3779b97f4a7c15ULL;
uint64_t next() {
    uint64_t z = (state += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

// ---------------------------------------------------------------------------------------------------------------------
// Binary64 and scaled binary64 (W8's S10.2): MachineEngine::advance (machine_candidates.cpp), SecularPoints::advance
// and NodeOffsets::refine (study_candidates.cpp), ProgramMachineStage::sweep (program_candidates.cpp) and
// small_roots (program_recurrence_seeds.cpp).
//
// Each computation runs between a clear and a read of the floating-point flags, its inputs read and its result written
// through volatile objects, so that the compiler cannot move it out from between them. The flags so read are those of
// this file's code, which the compiler may schedule differently from the library's (the library has no FENV_ACCESS
// either): equal flags show that the shared function and the old expression perform the same operations. The values'
// bits are the proof for the expressions; the library's own bytes are shown by cmp, or by the gate's identity run.

namespace machine = certpolysolve::detail::machine;
template <class R> using Complex = machine::Complex<R>;
using machine::ScaledDouble;

double launder(double x) {
    volatile double v = x;
    return v;
}
int64_t launder(int64_t x) {
    volatile int64_t v = x;
    return v;
}
ScaledDouble launder(ScaledDouble x) {
    ScaledDouble y;
    y.mantissa = launder(x.mantissa);
    y.exponent = launder(x.exponent);
    return y;
}
template <class R> Complex<R> launder(Complex<R> z) { return {launder(z.re), launder(z.im)}; }
volatile double kept_double;
volatile int64_t kept_integer;
void keep(double x) { kept_double = x; }
void keep(ScaledDouble x) {
    kept_double = x.mantissa;
    kept_integer = x.exponent;
}
template <class R> void keep(Complex<R> z) {
    keep(z.re);
    keep(z.im);
}

// Two NaNs count as equal. Which NaN an operation returns when both operands are NaN depends on their order, which C++
// leaves free for + and *, and MSVC orders them differently in the shared function and in the old forms (GCC and Clang
// do not). Every caller rejects a step that is not finite, so a NaN's sign or payload never reaches an answer.
bool same_bits(double a, double b) {
    return (std::isnan(a) && std::isnan(b)) || std::bit_cast<uint64_t>(a) == std::bit_cast<uint64_t>(b);
}
bool same_bits(ScaledDouble a, ScaledDouble b) { return same_bits(a.mantissa, b.mantissa) && a.exponent == b.exponent; }
template <class R> bool same_bits(Complex<R> a, Complex<R> b) { return same_bits(a.re, b.re) && same_bits(a.im, b.im); }

// A step, whether the kernel refused it (SecularPoints at a zero denominator), and the floating-point flags raised.
template <class R> struct Outcome {
    Complex<R> value;
    bool refused = false;
    int flags = 0;
};
template <class R, class F> Outcome<R> observe(F f) {
    std::feclearexcept(FE_ALL_EXCEPT);
    Outcome<R> outcome = f();
    keep(outcome.value);
    outcome.flags = std::fetestexcept(FE_ALL_EXCEPT);
    return outcome;
}
size_t binary64_compared = 0, scaled_compared = 0, secular_refusals = 0;
template <class R> void check(const Outcome<R> &shared, const Outcome<R> &reference, const char *what, size_t trial) {
    require(shared.refused == reference.refused && same_bits(shared.value, reference.value) &&
                shared.flags == reference.flags,
            what, trial);
    ++(std::is_same_v<R, double> ? binary64_compared : scaled_compared);
}

// The old expressions, at e454afc.
// MachineEngine::advance, the first form (machine_candidates.cpp:446-457, 507-512): correction = N and denominator = 1,
// then denominator -= correction * sum and delta = correction / denominator.
template <class R> Complex<R> reference_engine_first_form(Complex<R> newton, Complex<R> sum) {
    using C = Complex<R>;
    C correction, denominator;
    correction = newton;
    denominator = C(R(1));
    denominator -= correction * sum;
    const C delta = correction / denominator;
    return delta;
}
// The second form (machine_candidates.cpp:443-445, or 503-505 after the sums, then 507-508 and 512): correction = 1 and
// denominator = p'/p, then denominator -= sum and the same quotient.
template <class R> Complex<R> reference_engine_second_form(Complex<R> logarithmic_derivative, Complex<R> sum) {
    using C = Complex<R>;
    C correction, denominator;
    correction = C(R(1));
    denominator = logarithmic_derivative;
    denominator -= sum;
    const C delta = correction / denominator;
    return delta;
}
// The direct sums' term: MachineEngine's (machine_candidates.cpp:478, 497) and SecularPoints' general loop's
// (study_candidates.cpp:555).
template <class R> Complex<R> reference_engine_term(Complex<R> difference) {
    using C = Complex<R>;
    return C(R(1)) / difference;
}
// SecularPoints::advance at a pole (study_candidates.cpp:570).
template <class R> Complex<R> reference_secular_pole(Complex<R> sum) {
    using C = Complex<R>;
    return -C(R(1)) / sum;
}
// SecularPoints::advance elsewhere (study_candidates.cpp:572-574): a zero denominator refuses the packet.
template <class R> Outcome<R> reference_secular_step(Complex<R> q, Complex<R> sum) {
    using C = Complex<R>;
    const C denominator = C(R(1)) - q * sum;
    if (denominator == C{})
        return {C{}, true};
    return {q / denominator, false};
}
// NodeOffsets::refine (study_candidates.cpp:696, 699).
Complex<double> reference_offsets_term(Complex<double> between) {
    using C = Complex<double>;
    return C(1.0) / between;
}
Complex<double> reference_offsets_step(Complex<double> q, Complex<double> interaction) {
    using C = Complex<double>;
    const C step = q / (C(1.0) - q * interaction);
    return step;
}
// ProgramMachineStage::sweep (program_candidates.cpp).
Complex<double> reference_stage_delta(Complex<double> newton, Complex<double> sum) {
    using R = double;
    using C = Complex<R>;
    const C denominator = C{R(1), R{}} - newton * sum;
    const C delta = newton / denominator;
    return delta;
}
// small_roots (program_recurrence_seeds.cpp:114, 115); the term is taken of y[i] - y[j] there.
Complex<double> reference_small_term(Complex<double> difference) {
    using C = Complex<double>;
    return C{1, 0} / difference;
}
Complex<double> reference_small_step(Complex<double> ratio, Complex<double> sum) {
    using C = Complex<double>;
    const C step = ratio / (C{1, 0} - ratio * sum);
    return step;
}
// SecularPoints' step on the shared denominator, as its site now reads.
template <class R> Outcome<R> shared_secular_step(Complex<R> q, Complex<R> sum) {
    using C = Complex<R>;
    const C denominator = machine::aberth_denominator(q, sum);
    if (denominator == C{})
        return {C{}, true};
    return {q / denominator, false};
}

// Every comparison for the pair (a, b): a as N, p'/p, q, ratio or a difference, b as S.
template <class R> void compare(Complex<R> a, Complex<R> b, size_t trial) {
    using O = Outcome<R>;
    const O correction = observe<R>([&] { return O{machine::aberth_correction(launder(a), launder(b))}; });
    check(correction, observe<R>([&] { return O{reference_engine_first_form(launder(a), launder(b))}; }),
          "aberth_correction against MachineEngine's first form", trial);
    const O logarithmic = observe<R>([&] { return O{machine::logarithmic_aberth_correction(launder(a), launder(b))}; });
    check(logarithmic, observe<R>([&] { return O{reference_engine_second_form(launder(a), launder(b))}; }),
          "logarithmic_aberth_correction against MachineEngine's second form", trial);
    const O term = observe<R>([&] { return O{machine::interaction_term(launder(a))}; });
    check(term, observe<R>([&] { return O{reference_engine_term(launder(a))}; }),
          "interaction_term against MachineEngine's and SecularPoints' term", trial);
    check(observe<R>([&] { return O{machine::aberth_pole_limit(launder(b))}; }),
          observe<R>([&] { return O{reference_secular_pole(launder(b))}; }),
          "aberth_pole_limit against SecularPoints' pole", trial);
    const O secular = observe<R>([&] { return reference_secular_step(launder(a), launder(b)); });
    check(observe<R>([&] { return shared_secular_step(launder(a), launder(b)); }), secular,
          "aberth_denominator against SecularPoints' step", trial);
    secular_refusals += secular.refused;
    if constexpr (std::is_same_v<R, double>) {
        check(correction, observe<R>([&] { return O{reference_offsets_step(launder(a), launder(b))}; }),
              "aberth_correction against NodeOffsets' step", trial);
        check(correction, observe<R>([&] { return O{reference_stage_delta(launder(a), launder(b))}; }),
              "aberth_correction against the program stage's", trial);
        check(correction, observe<R>([&] { return O{reference_small_step(launder(a), launder(b))}; }),
              "aberth_correction against small_roots' step", trial);
        check(term, observe<R>([&] { return O{reference_offsets_term(launder(a))}; }),
              "interaction_term against NodeOffsets' term", trial);
        check(term, observe<R>([&] { return O{reference_small_term(launder(a))}; }),
              "interaction_term against small_roots' term", trial);
    }
}

// The values a component draws from: zeros of both signs, the subnormal and normal extremes, the infinities, NaNs (of
// both signs, one with a payload), one and its neighbours, and powers of two at the limits where a square or a
// product leaves the range.
std::vector<double> special_values() {
    using limits = std::numeric_limits<double>;
    const double infinity = limits::infinity(), nan = limits::quiet_NaN();
    const double payload = std::bit_cast<double>(std::bit_cast<uint64_t>(nan) | 0x5a5a5ULL);
    return {0.0, -0.0, limits::denorm_min(), -limits::denorm_min(), std::nextafter(limits::min(), 0.0),
            limits::min(), -limits::min(), limits::max(), -limits::max(), infinity, -infinity, nan, -nan, payload,
            1.0, -1.0, 0.5, 2.0, -2.0, std::nextafter(1.0, 2.0), std::nextafter(1.0, 0.0), 3.0, 0x1p-537, 0x1p511,
            -0x1p512, 0x1p1023};
}
// A double from every bit pattern alike (every exponent, both signs, the infinities and NaNs among them), or a normal
// one of exponent within +-80 (where sums cancel and quotients stay in range).
double any_double() { return std::bit_cast<double>(next()); }
double moderate_double() {
    const double mantissa = 1.0 + static_cast<double>(next() >> 12) * 0x1p-52;
    const double x = std::ldexp(mantissa, static_cast<int>(next() % 161) - 80);
    return next() % 2 ? -x : x;
}
Complex<double> draw(bool moderate) {
    Complex<double> z{moderate ? moderate_double() : any_double(), moderate ? moderate_double() : any_double()};
    if (next() % 8 == 0)
        z.im = next() % 2 ? z.re : -z.re; // |re| = |im|: Smith's division at its branch
    return z;
}

void binary64_inputs() {
    size_t trial = 0;
    const auto specials = special_values();
    std::vector<Complex<double>> grid;
    for (const double re : specials)
        for (const double im : specials)
            grid.push_back({re, im});
    for (const auto a : grid)
        for (const auto b : grid)
            compare(a, b, trial++);
    for (int t = 0; t < 400000; ++t, ++trial)
        compare(draw(t % 2 == 0), draw(t % 2 == 0), trial);
    // 1 - N S exactly zero: N = 2^k and S = 2^-k; N = 2^k (1 + i) and S = 2^-(k+1) (1 - i). Every k from the least
    // normal exponent to the largest: 2^-k and 2^-(k+1) are then normal or subnormal, both exact.
    for (int k = std::numeric_limits<double>::min_exponent - 1; k < std::numeric_limits<double>::max_exponent; ++k) {
        compare(Complex<double>{std::ldexp(1.0, k), 0.0}, Complex<double>{std::ldexp(1.0, -k), 0.0}, trial++);
        compare(Complex<double>{std::ldexp(1.0, k), std::ldexp(1.0, k)},
                Complex<double>{std::ldexp(1.0, -k - 1), -std::ldexp(1.0, -k - 1)}, trial++);
    }
}

// ScaledDouble's values: a binary64 mantissa normalized by frexp, and exponents over a wide range or within +-60.
ScaledDouble draw_scaled(int t) {
    if (next() % 32 == 0)
        return ScaledDouble(0.0);
    const double mantissa = 0.5 + static_cast<double>(next() >> 12) * 0x1p-53;
    const int64_t span = t % 3 == 0 ? 60 : t % 3 == 1 ? 3000 : int64_t(1) << 40;
    const int64_t exponent = static_cast<int64_t>(next() % static_cast<uint64_t>(2 * span + 1)) - span;
    return ScaledDouble(next() % 2 ? -mantissa : mantissa, exponent);
}
void scaled_inputs() {
    size_t trial = 0;
    for (int t = 0; t < 200000; ++t, ++trial) {
        Complex<ScaledDouble> a{draw_scaled(t), draw_scaled(t)}, b{draw_scaled(t), draw_scaled(t)};
        if (next() % 8 == 0)
            a.im = next() % 2 ? a.re : -a.re;
        compare(a, b, trial);
    }
    // 1 - N S exactly zero: N = 2^k and S = 2^-k, beyond binary64's range too.
    for (int64_t k = -5000; k <= 5000; k += 7)
        compare(Complex<ScaledDouble>{ScaledDouble(1.0, k)}, Complex<ScaledDouble>{ScaledDouble(1.0, -k)}, trial++);
}

// The program kernels on their own inputs, as machine_bits draws them (tests/test_program.cpp): its 300 random
// programs and Mandelbrot k = 8, at its 64 points and at 64 points beyond binary64's reach for them. N from the
// binary64 evaluator (program_machine.cpp:452) and, where the values leave binary64, from the range tier
// (machine::BlockDouble, its quotient converted to binary64: 481-482); S the stage's transposed sums over the points
// (conj(d) / |d|^2, j ascending, j != i: 556-570) and small_roots' (Smith's term, j ascending). Then small_roots'
// iteration itself on polynomials whose coefficients are those programs' values.
size_t range_tier_inputs = 0, trajectory_steps = 0;

// small_roots' unit vector (program_recurrence_seeds.cpp:64-70).
Complex<double> unit(uint64_t j, uint64_t m) {
    arb_t t, s, c;
    arb_init(t);
    arb_init(s);
    arb_init(c);
    arb_set_ui(t, 2 * j + 1);
    arb_div_ui(t, t, m, 64);
    arb_sin_cos_pi(s, c, t, 64);
    const Complex<double> u{arf_get_d(arb_midref(c), ARF_RND_NEAR), arf_get_d(arb_midref(s), ARF_RND_NEAR)};
    arb_clear(t);
    arb_clear(s);
    arb_clear(c);
    return u;
}
// small_roots' Aberth iteration (program_recurrence_seeds.cpp:91-124, verbatim but for its helpers' names), every term
// and step compared with the shared functions on the iteration's own values.
void small_roots_iteration(std::vector<Complex<double>> c, size_t &trial) {
    using C = Complex<double>;
    const auto norm = [](C a) { return a.re * a.re + a.im * a.im; };
    const auto finite = [](C a) { return std::isfinite(a.re) && std::isfinite(a.im); };
    while (c.size() > 1 && c.back().re == 0 && c.back().im == 0)
        c.pop_back();
    const size_t d = c.size() - 1;
    if (d < 3)
        return;
    double radius = 0;
    for (size_t j = 0; j < d; ++j)
        radius = std::max(radius, std::sqrt(norm(c[j] / c[d])));
    radius += 1;
    std::vector<C> y(d);
    for (size_t j = 0; j < d; ++j)
        y[j] = unit(j, d) * C{radius, 0};
    for (int sweep = 0; sweep < 200; ++sweep) {
        double largest = 0;
        for (size_t i = 0; i < d; ++i) {
            C value = c[d], derivative{0, 0};
            for (size_t j = d; j-- > 0;) {
                derivative = derivative * y[i] + value;
                value = value * y[i] + c[j];
            }
            if (norm(derivative) == 0)
                continue;
            const C ratio = value / derivative;
            C sum{0, 0};
            for (size_t j = 0; j < d; ++j)
                if (j != i && norm(y[i] - y[j]) > 0) {
                    using O = Outcome<double>;
                    check(observe<double>([&] { return O{machine::interaction_term(launder(y[i]) - launder(y[j]))}; }),
                          observe<double>([&] { return O{reference_small_term(launder(y[i]) - launder(y[j]))}; }),
                          "interaction_term against small_roots' term, on its iteration", trial);
                    sum = sum + C{1, 0} / (y[i] - y[j]);
                }
            using O = Outcome<double>;
            check(observe<double>([&] { return O{machine::aberth_correction(launder(ratio), launder(sum))}; }),
                  observe<double>([&] { return O{reference_small_step(launder(ratio), launder(sum))}; }),
                  "aberth_correction against small_roots' step, on its iteration", trial);
            ++trial;
            ++trajectory_steps;
            const C step = ratio / (C{1, 0} - ratio * sum);
            if (!finite(step))
                continue;
            y[i] = y[i] - step;
            largest = std::max(largest, norm(step) / std::max(norm(y[i]), 1e-300));
        }
        if (largest < 0x1p-100)
            return;
    }
}

// A point set with the stage's sums over it (program_machine.cpp:556-570, the transposed direct sum) and small_roots'
// (its Smith term, j ascending).
struct ProgramPoints {
    std::vector<Complex<double>> z, stage_sums, smith_sums;
    explicit ProgramPoints(std::vector<Complex<double>> points)
        : z(std::move(points)), stage_sums(z.size()), smith_sums(z.size()) {
        const size_t n = z.size();
        for (size_t j = 0; j < n; ++j)
            for (size_t i = 0; i < n; ++i) {
                if (i == j)
                    continue;
                const Complex<double> difference = z[i] - z[j];
                const double q = difference.re * difference.re + difference.im * difference.im;
                if (q == 0)
                    continue;
                stage_sums[i] += Complex<double>{difference.re / q, -difference.im / q};
            }
        for (size_t i = 0; i < n; ++i)
            for (size_t j = 0; j < n; ++j)
                if (j != i)
                    smith_sums[i] = smith_sums[i] + reference_small_term(z[i] - z[j]);
    }
};

void program_inputs() {
    using S = machine::BlockDouble;
    using SC = Complex<S>;
    std::mt19937_64 rng(20261016); // machine_bits' seed: the same programs in the same order
    std::vector<Complex<double>> near_points, far_points;
    for (int j = 0; j < 64; ++j)
        near_points.push_back({0.3 * std::cos(0.7 * j + 0.1), 0.3 * std::sin(0.7 * j + 0.1)});
    near_points[1] = {0x1p40, 0.5};
    near_points[2] = {0x1p-40, -0x1p-41};
    near_points[3] = {-1.5, 2.25};
    // And 64 points beyond binary64's reach for these programs, where the stage takes N from its range tier: a
    // degree-40 program's values overflow from about |z| = 2^(1024/40) = 2^25.6, so |z| runs from 2^26 to 2^60.
    constexpr int far_low_log2 = 26, far_high_log2 = 60;
    for (int j = 0; j < 64; ++j) {
        const double radius = std::ldexp(1.0, far_low_log2 + j % (far_high_log2 - far_low_log2 + 1));
        far_points.push_back({radius * std::cos(0.7 * j + 0.1), radius * std::sin(0.7 * j + 0.1)});
    }
    const ProgramPoints sets[] = {ProgramPoints(near_points), ProgramPoints(far_points)};
    size_t trial = 0;
    const auto check_program = [&](const d::Program &p, size_t draw) {
        d::ProgramMachine<double> m(p);
        if (!m.valid())
            return;
        std::unique_ptr<d::ProgramMachine<S>> range_tier;
        std::vector<d::ProgramMachine<double>::Evaluation> near_values;
        for (const ProgramPoints &set : sets) {
            const size_t n = set.z.size();
            std::vector<d::ProgramMachine<double>::Evaluation> e(n);
            m.evaluate(set.z, e, n);
            for (size_t k = 0; k < n; ++k) {
                Complex<double> newton;
                if (e[k].finite) {
                    if (e[k].derivative == Complex<double>{})
                        continue;
                    newton = e[k].value / e[k].derivative;
                } else {
                    if (!range_tier)
                        range_tier = std::make_unique<d::ProgramMachine<S>>(m);
                    if (!range_tier->valid())
                        continue;
                    const auto f = range_tier->evaluate_one({S(set.z[k].re), S(set.z[k].im)});
                    if (!f.finite || f.derivative == SC{})
                        continue;
                    const SC quotient = f.value / f.derivative;
                    newton = {static_cast<double>(quotient.re), static_cast<double>(quotient.im)};
                    ++range_tier_inputs;
                }
                compare(newton, set.stage_sums[k], trial++);
                compare(newton, set.smith_sums[k], trial++);
            }
            if (near_values.empty())
                near_values = std::move(e);
        }
        // small_roots on c_0 + ... + c_d y^d, the c_j the program's values at machine_bits' first d + 1 points.
        const size_t degree = 3 + draw % 6;
        std::vector<Complex<double>> coefficients;
        for (size_t j = 0; j <= degree; ++j)
            if (near_values[j].finite)
                coefficients.push_back(near_values[j].value);
        small_roots_iteration(coefficients, trial);
    };
    for (int t = 0; t < 300; ++t) {
        const auto p = d::random_program(rng, 10 + t % 30, 40);
        if (!p.output)
            continue;
        check_program(p, static_cast<size_t>(t));
    }
    check_program(d::detail_of(d::mandelbrot(8)), 300);
    // Mandelbrot k = 8 (degree 255) alone overflows binary64 at every far point, so the range tier serves at least
    // those.
    require(range_tier_inputs >= far_points.size(), "the range tier's corrections are among the program inputs",
            range_tier_inputs);
    require(trajectory_steps > 1000, "small_roots' iteration took its steps", trajectory_steps);
}

// ---------------------------------------------------------------------------------------------------------------------
// Multiprecision: candidate_packet and secular_candidate_packet (candidates.cpp), W8's S10.3.

// candidates.cpp's Point and its zero and zero test.
struct Point {
    acf_t value;
    Point() { acf_init(value); }
    ~Point() { acf_clear(value); }
    Point(const Point &x) : Point() { acf_set(value, x.value); }
    Point &operator=(const Point &x) {
        acf_set(value, x.value);
        return *this;
    }
};
constexpr auto rounding = ARF_RND_NEAR;
void zero(acf_t x) {
    arf_zero(acf_realref(x));
    arf_zero(acf_imagref(x));
}
bool is_zero(const acf_t x) { return arf_is_zero(acf_realref(x)) && arf_is_zero(acf_imagref(x)); }

// candidate_packet's interaction loop, its Aberth arm (candidates.cpp:1189-1203; z is points[i].value there).
bool reference_packet_interactions(acf_t interaction, acf_srcptr z, const std::vector<Point> &points,
                                   const std::vector<size_t> &neighbours, size_t i, slong prec, acf_t difference,
                                   acf_t inverse, uint64_t &root_interactions) {
    bool valid = true;
    for (size_t j : neighbours)
        if (j != i) {
            acf_sub(difference, z, points[j].value, prec, rounding);
            if (is_zero(difference)) {
                valid = false;
                break;
            }
            acf_approx_inv(inverse, difference, prec, rounding);
            acf_add(interaction, interaction, inverse, prec, rounding);
            ++root_interactions;
        }
    return valid;
}
// secular_candidate_packet's interaction loop (candidates.cpp:676-686).
bool reference_secular_interactions(acf_t interaction, const std::vector<Point> &points,
                                    const std::vector<size_t> &neighbours, size_t i, slong prec, acf_t difference,
                                    acf_t inverse, uint64_t &root_interactions) {
    bool valid = true;
    for (size_t j : neighbours)
        if (j != i) {
            acf_sub(difference, points[i].value, points[j].value, prec, rounding);
            if (is_zero(difference)) {
                valid = false;
                break;
            }
            acf_approx_inv(inverse, difference, prec, rounding);
            acf_add(interaction, interaction, inverse, prec, rounding);
            ++root_interactions;
        }
    return valid;
}
// candidate_packet's correction (candidates.cpp:1212-1220): false where the packet skips the point.
bool reference_packet_delta(acf_t delta, const acf_t correction, const acf_t interaction, acf_t denominator,
                            slong prec) {
    acf_mul(denominator, correction, interaction, prec, rounding);
    acf_neg(denominator, denominator);
    arf_add_ui(acf_realref(denominator), acf_realref(denominator), 1, prec, rounding);
    if (is_zero(denominator))
        return false;
    acf_approx_div(delta, correction, denominator, prec, rounding);
    return true;
}
// secular_candidate_packet's correction away from a pole (candidates.cpp:702-707).
bool reference_secular_delta(acf_t delta, const acf_t correction, const acf_t interaction, acf_t denominator,
                             slong prec) {
    acf_mul(denominator, correction, interaction, prec, rounding);
    acf_neg(denominator, denominator);
    arf_add_ui(acf_realref(denominator), acf_realref(denominator), 1, prec, rounding);
    const bool valid = !is_zero(denominator);
    if (valid)
        acf_approx_div(delta, correction, denominator, prec, rounding);
    return valid;
}

// The precisions of the multiprecision packets, from binary64's neighbourhood to a long ladder's top.
constexpr slong precisions[] = {53, 64, 65, 128, 200, 256, 333, 512, 1000, 2048, 4096};

// A random acf value: FLINT's randtest mantissas (long runs of ones and zeros among them) of up to prec + 8 bits, with
// exponents up to 2^magnitude_bits in size.
void random_value(acf_t x, flint_rand_t random, slong prec, slong magnitude_bits) {
    acf_randtest(x, random, prec + 8, magnitude_bits);
}

void interactions(flint_rand_t random) {
    Point sum[2], difference[2], inverse[2];
    size_t trial = 0;
    for (const slong prec : precisions)
        for (int t = 0; t < 400; ++t, ++trial) {
            const size_t n = 1 + next() % 12;
            const slong magnitude_bits = t % 3 == 0 ? 2 : t % 3 == 1 ? 8 : 20;
            std::vector<Point> points(n);
            for (auto &p : points)
                random_value(p.value, random, prec, magnitude_bits);
            const size_t i = next() % n;
            // A point on z (a coincident pair), or one a unit in the last place of z's real part away.
            if (n > 1 && next() % 6 == 0) {
                const size_t j = (i + 1 + next() % (n - 1)) % n;
                points[j] = points[i];
                if (next() % 2) {
                    const arf_struct *re = acf_realref(points[i].value);
                    Point unit;
                    arf_one(acf_realref(unit.value));
                    arf_mul_2exp_si(acf_realref(unit.value), acf_realref(unit.value),
                                    (arf_is_zero(re) ? 0 : arf_abs_bound_lt_2exp_si(re)) - prec);
                    arf_add(acf_realref(points[j].value), re, acf_realref(unit.value), ARF_PREC_EXACT, ARF_RND_DOWN);
                }
            }
            // The neighbours: every point, or a tree's near list (a subset in its own order, which may hold i).
            std::vector<size_t> neighbours(n);
            std::iota(neighbours.begin(), neighbours.end(), size_t{0});
            if (next() % 2) {
                for (size_t k = n; k > 1; --k)
                    std::swap(neighbours[k - 1], neighbours[next() % k]);
                neighbours.resize(next() % (n + 1));
            }
            // The sum starts at zero or at a far field.
            Point start;
            if (next() % 2)
                random_value(start.value, random, prec, magnitude_bits);
            for (int form = 0; form < 2; ++form) {
                uint64_t counts[2] = {7, 7};
                bool results[2];
                for (int k = 0; k < 2; ++k) {
                    acf_set(sum[k].value, start.value);
                    zero(difference[k].value);
                    zero(inverse[k].value);
                }
                results[0] = form == 0 ? reference_packet_interactions(sum[0].value, points[i].value, points,
                                                                       neighbours, i, prec, difference[0].value,
                                                                       inverse[0].value, counts[0])
                                       : reference_secular_interactions(sum[0].value, points, neighbours, i, prec,
                                                                        difference[0].value, inverse[0].value,
                                                                        counts[0]);
                results[1] = d::add_interactions(sum[1].value, points[i].value, points, neighbours, i, prec,
                                                 difference[1].value, inverse[1].value, counts[1]);
                require(results[0] == results[1] && counts[0] == counts[1] && acf_equal(sum[0].value, sum[1].value) &&
                            acf_equal(difference[0].value, difference[1].value) &&
                            acf_equal(inverse[0].value, inverse[1].value),
                        form == 0 ? "add_interactions against candidate_packet's loop"
                                  : "add_interactions against secular_candidate_packet's loop",
                        trial);
                ++compared;
                coincident += !results[0];
            }
        }
}

void corrections(flint_rand_t random) {
    Point newton, sum, delta[2], denominator[2];
    size_t trial = 0;
    for (const slong prec : precisions)
        for (int t = 0; t < 600; ++t, ++trial) {
            const slong magnitude_bits = t % 3 == 0 ? 2 : t % 3 == 1 ? 8 : 20;
            random_value(newton.value, random, prec, magnitude_bits);
            random_value(sum.value, random, prec, magnitude_bits);
            // 1 - N S exactly zero: N = 2^k and S = 2^-k, or N = 1 + i and S = (1 - i) / 2.
            if (t % 10 == 0) {
                const slong k = static_cast<slong>(next() % 2001) - 1000;
                zero(newton.value);
                arf_one(acf_realref(newton.value));
                arf_mul_2exp_si(acf_realref(newton.value), acf_realref(newton.value), k);
                zero(sum.value);
                arf_one(acf_realref(sum.value));
                arf_mul_2exp_si(acf_realref(sum.value), acf_realref(sum.value), -k);
            } else if (t % 10 == 1) {
                arf_one(acf_realref(newton.value));
                arf_one(acf_imagref(newton.value));
                arf_set_d(acf_realref(sum.value), 0.5);
                arf_set_d(acf_imagref(sum.value), -0.5);
            } else if (t % 10 == 2)
                zero(newton.value);
            else if (t % 10 == 3)
                zero(sum.value);
            for (int form = 0; form < 2; ++form) {
                bool results[2];
                for (int k = 0; k < 2; ++k) {
                    // A value the correction must leave alone where it declines.
                    arf_set_si(acf_realref(delta[k].value), 12345);
                    arf_set_si(acf_imagref(delta[k].value), -678);
                    zero(denominator[k].value);
                }
                results[0] = form == 0 ? reference_packet_delta(delta[0].value, newton.value, sum.value,
                                                                denominator[0].value, prec)
                                       : reference_secular_delta(delta[0].value, newton.value, sum.value,
                                                                 denominator[0].value, prec);
                results[1] = d::aberth_delta(delta[1].value, newton.value, sum.value, denominator[1].value, prec);
                require(results[0] == results[1] && acf_equal(delta[0].value, delta[1].value) &&
                            acf_equal(denominator[0].value, denominator[1].value),
                        form == 0 ? "aberth_delta against candidate_packet's correction"
                                  : "aberth_delta against secular_candidate_packet's correction",
                        trial);
                ++compared;
                zero_denominators += !results[0];
            }
        }
}
} // namespace

int main() {
    flint_rand_t random;
    flint_rand_init(random);
    flint_rand_set_seed(random, 20261023, 10);
    interactions(random);
    corrections(random);
    binary64_inputs();
    scaled_inputs();
    const size_t before_programs = binary64_compared;
    program_inputs();
    const size_t program_compared = binary64_compared - before_programs;
    flint_rand_clear(random);
    flint_cleanup();
    require(coincident > 0 && zero_denominators > 0, "the refusals are among the multiprecision cases", 0);
    require(secular_refusals > 0, "SecularPoints' refusals are among the binary64 cases", 0);
    if (failures) {
        std::fprintf(stderr, "%d of %zu comparisons differ\n", failures,
                     binary64_compared + scaled_compared + compared);
        return 1;
    }
    std::printf("aberth_step: %zu binary64 comparisons (%zu SecularPoints refusals; %zu on the program kernels' own "
                "values: %zu with N from the range tier, %zu small_roots steps), %zu scaled binary64, %zu "
                "multiprecision (%zu at a coincident point, %zu at a zero denominator); every one the same\n",
                binary64_compared, secular_refusals, program_compared, range_tier_inputs, trajectory_steps,
                scaled_compared, compared, coincident, zero_denominators);
    return 0;
}
