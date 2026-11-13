// portable::log2 and portable::exp2 (src/core/portable_math.hpp) against the C library's: within a stated number of units
// in the last place over sweeps of binary64 values (for log2 every binade, the integers that decisions take it of, the
// neighbourhoods of 1 and of the branch point sqrt(1/2); for exp2 the meter's range densely, the neighbourhoods of the
// half-integers where the reduction's integer changes, and the whole domain), exact at the powers of two and the
// integers, NaN outside the domain. The C library is the reference here because glibc's log2 and exp2 are accurate
// to within one unit; the functions themselves never call them.
#include "core/portable_math.hpp"
#include <bit>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <numbers>

namespace {
int failures = 0;
void require(bool ok, const char *what) {
    if (!ok) {
        std::fprintf(stderr, "FAIL: %s\n", what);
        ++failures;
    }
}
// The distance between a and b in units in the last place of b (b finite).
double ulps(double a, double b) {
    if (a == b) return 0;
    const double unit = std::nextafter(std::abs(b), std::numeric_limits<double>::infinity()) - std::abs(b);
    return std::abs(a - b) / unit;
}
// A fixed sequence (SplitMix64), the same on every run.
uint64_t state = 0x9e3779b97f4a7c15ULL;
uint64_t next() {
    uint64_t z = (state += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}
} // namespace

int main() {
    namespace p = certpolysolve::detail::portable;
    // The bound the header states: a few units; 8 leaves room over the analysis's 6.
    constexpr double bound = 8;
    double worst = 0, worst_at = 0;
    const auto check = [&](double x) {
        const double error = ulps(p::log2(x), std::log2(x));
        if (error > worst) {
            worst = error;
            worst_at = x;
        }
    };
    // The series' coefficients are the correctly rounded quotients.
    for (size_t j = 0; j < p::atanh_series.size(); ++j)
        require(p::atanh_series[j] == 1.0 / static_cast<double>(2 * j + 1), "a series coefficient");
    // Exact at every power of two, normal and subnormal.
    for (int k = std::numeric_limits<double>::min_exponent - std::numeric_limits<double>::digits;
         k < std::numeric_limits<double>::max_exponent; ++k)
        require(p::log2(std::ldexp(1.0, k)) == static_cast<double>(k), "exact at a power of two");
    // The integers 1..2^20 (degrees, supports, their quotients' arguments) and their quotients by small supports.
    for (uint64_t n = 1; n <= (uint64_t{1} << 20); ++n) {
        check(static_cast<double>(n));
        check(static_cast<double>(n) / 3);
        check(static_cast<double>(n) / 127);
    }
    // Around 1 and the branch point.
    for (int k = 1; k <= 60; ++k) {
        check(1 + std::ldexp(1.0, -k));
        check(1 - std::ldexp(1.0, -k - 1));
    }
    const double branch = std::sqrt(0.5);
    double below = branch, above = branch;
    for (int k = 0; k < 1000; ++k) {
        below = std::nextafter(below, 0.0);
        above = std::nextafter(above, 1.0);
        check(below);
        check(above);
        check(2 * below);
        check(2 * above);
    }
    // Random significands in every binade, normal and subnormal.
    constexpr int fraction_bits = std::numeric_limits<double>::digits - 1;
    constexpr uint64_t infinite_exponent = 0x7ff;
    for (int k = 0; k < 2000000; ++k) {
        const uint64_t bits = next() >> 1; // the sign bit clear
        if ((bits >> fraction_bits) == infinite_exponent) continue;
        const double x = std::bit_cast<double>(bits);
        if (x > 0) check(x);
    }
    check(std::numeric_limits<double>::max());
    check(std::numeric_limits<double>::min());
    check(std::numeric_limits<double>::denorm_min());
    std::printf("portable::log2 against the C library: at most %.2f units in the last place (at %a)\n", worst,
                worst_at);
    require(worst <= bound, "within the stated bound");
    // Outside the domain.
    for (const double x : {0.0, -0.0, -1.0, std::numeric_limits<double>::infinity(),
                           -std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()})
        require(std::isnan(p::log2(x)), "NaN outside the domain");

    // portable::exp2. Its coefficients: the first two exact (1, and ln 2 correctly rounded), every next one ln 2 / j
    // times the one before, to within the roundings of the two values and of this check (six at most, of 2^-53).
    require(p::exp2_series[0] == 1 && p::exp2_series[1] == std::numbers::ln2, "an exp2 coefficient");
    for (size_t j = 2; j < p::exp2_series.size(); ++j) {
        const double ratio = p::exp2_series[j] * static_cast<double>(j) / (p::exp2_series[j - 1] * std::numbers::ln2);
        require(std::abs(ratio - 1) <= 0x1p-50, "an exp2 coefficient's ratio to the one before");
    }
    // Exact at every integer whose power is a binary64 value, normal and subnormal.
    for (int k = std::numeric_limits<double>::min_exponent - std::numeric_limits<double>::digits;
         k < std::numeric_limits<double>::max_exponent; ++k)
        require(p::exp2(static_cast<double>(k)) == std::ldexp(1.0, k), "exp2 exact at an integer");
    // The bound the header states: a unit or two; 4 leaves room.
    constexpr double exp2_bound = 4;
    double worst_exp2 = 0, worst_exp2_at = 0;
    const auto check_exp2 = [&](double x) {
        const double error = ulps(p::exp2(x), std::exp2(x));
        if (error > worst_exp2) {
            worst_exp2 = error;
            worst_exp2_at = x;
        }
    };
    const auto uniform = [] { return static_cast<double>(next() >> 11) * 0x1p-53; }; // in [0, 1)
    // The meter's range: the base-2 logarithms of weights, counts and passes and their interpolations.
    constexpr double meter_range = 64;
    for (int k = 0; k < 2000000; ++k)
        check_exp2(meter_range * (2 * uniform() - 1));
    // Around the half-integers, where round(x) changes.
    for (int j = -70; j <= 70; ++j) {
        double below = j + 0.5, above = j + 0.5;
        check_exp2(below);
        for (int s = 0; s < 100; ++s) {
            below = std::nextafter(below, -std::numeric_limits<double>::infinity());
            above = std::nextafter(above, std::numeric_limits<double>::infinity());
            check_exp2(below);
            check_exp2(above);
        }
    }
    // The whole domain, subnormal results and the last finite binade included.
    const double lowest = std::numeric_limits<double>::min_exponent - std::numeric_limits<double>::digits - 1,
                 highest = std::numeric_limits<double>::max_exponent;
    for (int k = 0; k < 1000000; ++k)
        check_exp2(lowest + (highest - lowest) * uniform());
    check_exp2(std::nextafter(highest, 0.0));
    std::printf("portable::exp2 against the C library: at most %.2f units in the last place (at %a)\n", worst_exp2,
                worst_exp2_at);
    require(worst_exp2 <= exp2_bound, "exp2 within the stated bound");
    // The ends of the domain.
    require(std::isnan(p::exp2(std::numeric_limits<double>::quiet_NaN())), "exp2 of NaN");
    require(p::exp2(highest) == std::numeric_limits<double>::infinity() &&
                p::exp2(std::numeric_limits<double>::infinity()) == std::numeric_limits<double>::infinity(),
            "exp2 above the finite range");
    require(p::exp2(lowest) == 0 && p::exp2(-std::numeric_limits<double>::infinity()) == 0, "exp2 below the subnormals");
    if (failures) return EXIT_FAILURE;
    std::puts("portable_math: all checks pass");
    return EXIT_SUCCESS;
}
