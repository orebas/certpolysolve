// Number's two conversions on their fast paths since 26 October 2026 (number.cpp). rational_with_budget writes an
// integer mantissa times 2^exponent out by shifts: the rational must be the canonical one the power and the gcd gave,
// for every mantissa (odd, even with many factors of two, negative, one) and every exponent of either sign, with the
// same budget charged; decimal and fractional numbers keep the general path and are checked the same way. enclose
// sets an integer significand's ball with arb_set_round_fmpz: the ball must be arb_set_fmpq's, midpoint and radius,
// at every precision, also for integers longer than the precision.
#include "core/number.hpp"
#include "core/raii.hpp"
#include <flint/arb.h>
#include <flint/fmpq.h>
#include <cstdint>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>

namespace d = certpolysolve::detail;

namespace {
void require(bool ok, const std::string &why) {
    if (!ok) throw std::runtime_error(why);
}
// The reference: the significand times radix^exponent, made with a power and canonicalised.
void reference(fmpq_t out, const d::Number &x) {
    fmpq_set(out, x.significand.value());
    d::Integer magnitude, power;
    fmpz_abs(magnitude.value(), x.exponent.value());
    fmpz_ui_pow_ui(power.value(), x.radix, fmpz_get_ui(magnitude.value()));
    if (fmpz_sgn(x.exponent.value()) < 0)
        fmpz_mul(fmpq_denref(out), fmpq_denref(out), power.value());
    else
        fmpz_mul(fmpq_numref(out), fmpq_numref(out), power.value());
    fmpq_canonicalise(out);
}
void check(const std::string &text) {
    const d::Number x(text);
    d::Rational got, expected;
    uint64_t budget = UINT64_MAX / 4;
    require(x.rational_with_budget(got.value(), budget, true), text + ": refused");
    reference(expected.value(), x);
    require(fmpq_equal(got.value(), expected.value()) != 0, text + ": not the canonical rational");
    require(fmpq_is_canonical(got.value()) != 0, text + ": not canonical");
}
// The ball Number::enclose gave before its integer path: arb_set_fmpq, then the power of the radix.
void enclose_reference(arb_t out, const d::Number &x, slong precision) {
    arb_set_fmpq(out, x.significand.value(), precision);
    if (x.zero() || fmpz_is_zero(x.exponent.value()))
        return;
    if (x.radix == 2) {
        arb_mul_2exp_fmpz(out, out, x.exponent.value());
    } else {
        d::RealBall power, base;
        arb_set_ui(base.value(), 10);
        arb_pow_fmpz(power.value(), base.value(), x.exponent.value(), precision);
        arb_mul(out, out, power.value(), precision);
    }
}
uint64_t check_enclose(const std::string &text) {
    const d::Number x(text);
    uint64_t checked = 0;
    for (const slong precision : {slong(2), slong(10), slong(53), slong(64), slong(65), slong(128), slong(200),
                                  slong(333)}) {
        d::RealBall got, expected;
        x.enclose(got.value(), precision);
        enclose_reference(expected.value(), x, precision);
        require(arb_equal(got.value(), expected.value()) != 0,
                text + ": enclose differs from arb_set_fmpq's ball at precision " + std::to_string(precision));
        ++checked;
    }
    return checked;
}
} // namespace

int main() {
    try {
        std::mt19937_64 rng(20261026);
        std::uniform_int_distribution<int64_t> mantissa(-(int64_t(1) << 62), int64_t(1) << 62);
        std::uniform_int_distribution<int> exponent(-3000, 3000), zeros(0, 70);
        uint64_t checked = 0;
        for (int trial = 0; trial < 20000; ++trial) {
            int64_t m = mantissa(rng);
            if (trial % 4 == 0) m = (m >> 20) << (zeros(rng) % 40); // many factors of two
            if (trial % 97 == 0) m = trial % 2 ? 1 : -1;
            if (m == 0) m = 3;
            check(std::to_string(m) + "@" + std::to_string(exponent(rng)));
            ++checked;
        }
        for (const char *text : {"1@0", "-1@-1", "6@-1", "6@-2", "1024@-10", "1024@-11", "3@4000", "-5@-4000",
                                 "0.125", "-2.5e-3", "7/12", "1e300", "123456789012345678901234567890@-77"}) {
            check(text);
            ++checked;
        }
        // The budget charged is the one before the shifts: the significand's bits, and one bit for every power of
        // two (four for every power of ten).
        struct Charge {
            const char *text;
            uint64_t exponent_magnitude, bits_per_power;
        };
        for (const Charge &c : {Charge{"12345@-300", 300, 1}, Charge{"-12288@41", 41, 1}, Charge{"-2.5e-3", 4, 4}}) {
            const d::Number x(c.text);
            const uint64_t start = uint64_t(1) << 20;
            uint64_t left = start;
            d::Rational q;
            require(x.rational_with_budget(q.value(), left, true), std::string(c.text) + ": refused");
            const uint64_t charged = fmpz_bits(fmpq_numref(x.significand.value())) +
                                     fmpz_bits(fmpq_denref(x.significand.value())) +
                                     c.exponent_magnitude * c.bits_per_power;
            require(start - left == charged,
                    std::string(c.text) + ": the budget charged is not the significand's bits and the exponent's");
            ++checked;
        }
        // enclose: integer significands of 1 to 448 bits (longer than most of the precisions), binary and decimal
        // exponents, and fractions and decimals on the general path.
        uint64_t enclosed = 0;
        d::Integer big;
        for (int trial = 0; trial < 2000; ++trial) {
            fmpz_zero(big.value());
            const int words = 1 + trial % 7;
            for (int k = 0; k < words; ++k) {
                fmpz_mul_2exp(big.value(), big.value(), FLINT_BITS);
                fmpz_add_ui(big.value(), big.value(), static_cast<ulong>(rng()));
            }
            fmpz_tdiv_q_2exp(big.value(), big.value(), static_cast<flint_bitcnt_t>(rng() % (FLINT_BITS * words)));
            if (trial % 2)
                fmpz_neg(big.value(), big.value());
            char *digits = fmpz_get_str(nullptr, 10, big.value());
            const std::string m(digits);
            flint_free(digits);
            const int e = exponent(rng) / 10;
            enclosed += check_enclose(m + "@" + std::to_string(e));
            enclosed += check_enclose(m + "e" + std::to_string(e / 10));
        }
        for (const char *text : {"0", "1", "-1", "7/12", "-22/7", "0.1", "1e-300", "3@-1100"})
            enclosed += check_enclose(text);
        std::cout << "number rational: " << checked << " numbers written out, each the canonical rational; "
                  << enclosed << " balls enclosed as arb_set_fmpq encloses them\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << "test_number_rational: " << e.what() << '\n';
        return 1;
    }
}
