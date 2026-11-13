// arf_to_binary64, next_up and next_down (arf_binary64.hpp) against what they replace: arf_get_d in every rounding
// mode, and std::nextafter toward both infinities, bit for bit. Random arf numbers of one to four limbs with exponents
// across binary64's range and beyond it, both signs; then the edges: mantissas of exactly 53 bits, ties to even
// both ways, a carry into the next binade, the exponents where the fast path hands over, zero, the infinities and NaN.
#include "numeric/arithmetic/arf_binary64.hpp"
#include "core/raii.hpp"
#include <flint/fmpz.h>
#include <bit>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>

namespace d = certpolysolve::detail;

namespace {
void require(bool ok, const std::string &why) {
    if (!ok) throw std::runtime_error(why);
}
constexpr arf_rnd_t modes[] = {ARF_RND_DOWN, ARF_RND_UP, ARF_RND_FLOOR, ARF_RND_CEIL, ARF_RND_NEAR};
bool same(double a, double b) {
    return (std::isnan(a) && std::isnan(b)) || std::bit_cast<uint64_t>(a) == std::bit_cast<uint64_t>(b);
}
uint64_t compare_modes(const arf_t x, const std::string &what) {
    for (const arf_rnd_t rnd : modes)
        require(same(d::arf_to_binary64(x, rnd), arf_get_d(x, rnd)),
                what + ": arf_to_binary64 differs from arf_get_d in mode " + std::to_string(rnd));
    return std::size(modes);
}
// A decimal integer scaled by a power of two so that 2^(e-1) <= |x| < 2^e: the arf exponent e.
void set(arf_t x, const char *mantissa, slong e) {
    d::Integer m;
    fmpz_set_str(m.value(), mantissa, 10);
    arf_set_fmpz(x, m.value());
    arf_mul_2exp_si(x, x, e - arf_abs_bound_lt_2exp_si(x));
}
} // namespace

int main() {
    try {
        std::mt19937_64 rng(20261026);
        uint64_t checked = 0;
        d::Float x;
        d::Integer m;
        for (int trial = 0; trial < 200000; ++trial) {
            const int limbs = 1 + trial % 4;
            fmpz_zero(m.value());
            for (int k = 0; k < limbs; ++k) {
                fmpz_mul_2exp(m.value(), m.value(), FLINT_BITS);
                fmpz_add_ui(m.value(), m.value(), static_cast<ulong>(rng()));
            }
            if (trial % 5 == 0)  // trailing zeros: short mantissas, exact values
                fmpz_tdiv_q_2exp(m.value(), m.value(), static_cast<flint_bitcnt_t>(rng() % (FLINT_BITS * limbs)));
            if (fmpz_is_zero(m.value()))
                fmpz_one(m.value());
            if (rng() & 1)
                fmpz_neg(m.value(), m.value());
            arf_set_fmpz(x.value(), m.value());
            const slong e = static_cast<slong>(rng() % 2400) - 1200 - static_cast<slong>(fmpz_bits(m.value()));
            arf_mul_2exp_si(x.value(), x.value(), e);
            checked += compare_modes(x.value(), "random arf");
        }
        // The edges: 53 bits exactly, 54-bit ties (to even, down and up), a carry from an all-ones mantissa, and the
        // exponents about the fast path's range (-1020, 1020].
        for (const slong e : {slong(-1080), slong(-1075), slong(-1060), slong(-1021), slong(-1020), slong(-1019),
                              slong(-60), slong(0), slong(60), slong(1019), slong(1020), slong(1021), slong(1030)})
            for (const char *mantissa : {"1", "9007199254740991", "9007199254740993", "9007199254740995",
                                         "18014398509481983", "18446744073709551615", "-9007199254740993",
                                         "-18446744073709551615", "340282366920938463463374607431768211455"}) {
                set(x.value(), mantissa, e);
                checked += compare_modes(x.value(), std::string("edge ") + mantissa + " at exponent " + std::to_string(e));
            }
        arf_zero(x.value());
        checked += compare_modes(x.value(), "zero");
        arf_pos_inf(x.value());
        checked += compare_modes(x.value(), "+inf");
        arf_neg_inf(x.value());
        checked += compare_modes(x.value(), "-inf");
        arf_nan(x.value());
        checked += compare_modes(x.value(), "nan");
        // A bignum exponent (beyond a small fmpz).
        set(x.value(), "3", 0);
        const slong large = WORD(1) << (FLINT_BITS - 2);
        arf_mul_2exp_si(x.value(), x.value(), large);
        arf_mul_2exp_si(x.value(), x.value(), large);
        checked += compare_modes(x.value(), "bignum exponent");
        // next_up and next_down against std::nextafter: random bit patterns, then the edges.
        constexpr double inf = std::numeric_limits<double>::infinity();
        uint64_t neighbours = 0;
        auto neighbour_check = [&](double v) {
            require(same(d::next_up(v), std::nextafter(v, inf)), "next_up differs from nextafter");
            require(same(d::next_down(v), std::nextafter(v, -inf)), "next_down differs from nextafter");
            ++neighbours;
        };
        for (int trial = 0; trial < 1000000; ++trial)
            neighbour_check(std::bit_cast<double>(rng()));
        for (const double v : {0.0, -0.0, std::numeric_limits<double>::denorm_min(),
                               -std::numeric_limits<double>::denorm_min(), std::numeric_limits<double>::min(),
                               -std::numeric_limits<double>::min(), std::numeric_limits<double>::max(),
                               -std::numeric_limits<double>::max(), inf, -inf, 1.0, -1.0,
                               std::numeric_limits<double>::quiet_NaN()})
            neighbour_check(v);
        std::cout << "arf to binary64: " << checked << " conversions as arf_get_d gives them; " << neighbours
                  << " neighbours as std::nextafter gives them\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << "test_arf_binary64: " << e.what() << '\n';
        return 1;
    }
}
