// Endpoint::text (cover.hpp, certification.cpp) against what it replaces: arf_get_fmpz_2exp's mantissa and exponent
// written by fmpz_get_str, character for character. Random arf numbers of one to four limbs (the fast path takes one
// and two, FLINT's path the rest), both signs, exponents from tiny to beyond a word (an mpz exponent takes FLINT's
// path); then the edges: zero, powers of two, all-ones limbs, a mantissa with its lowest bit at each limb boundary,
// the largest 128-bit mantissa, the decimal lengths and chunk boundaries (10^k +- 1 and the like, to 10^39), and the
// exponents at COEFF_MAX. Also the copy, move and negation of an Endpoint.
#include "core/cover.hpp"
#include "core/raii.hpp"
#include <flint/fmpz.h>
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
std::string fmpz_text(const fmpz_t x) {
    char *s = fmpz_get_str(nullptr, 10, x);
    std::string out(s);
    flint_free(s);
    return out;
}
// The texts the API wrote before 27 October 2026.
d::EndpointText reference(const arf_t x) {
    d::Integer m, e;
    arf_get_fmpz_2exp(m.value(), e.value(), x);
    return {fmpz_text(m.value()), fmpz_text(e.value())};
}
uint64_t check(const arf_t x, const std::string &what) {
    const d::Endpoint end(x);
    const auto want = reference(x), got = end.text();
    require(got.mantissa == want.mantissa && got.exponent == want.exponent,
            what + ": " + got.mantissa + "@" + got.exponent + " where FLINT gives " + want.mantissa + "@" + want.exponent);
    require(end.mantissa_text() == want.mantissa && end.exponent_text() == want.exponent, what + ": the single texts");
    return 1;
}
// x = (a mantissa of `limbs` random limbs, its top bit set and its low bits cleared below `low`) * 2^e, sign s.
void random_arf(arf_t x, std::mt19937_64 &rng, int limbs, int low, const fmpz_t e, bool negative) {
    d::Integer m;
    fmpz_zero(m.value());
    for (int i = 0; i < limbs; ++i) {
        fmpz_mul_2exp(m.value(), m.value(), FLINT_BITS);
        fmpz_add_ui(m.value(), m.value(), rng());
    }
    fmpz_setbit(m.value(), static_cast<ulong>(limbs) * FLINT_BITS - 1);
    fmpz_fdiv_q_2exp(m.value(), m.value(), static_cast<ulong>(low));
    fmpz_mul_2exp(m.value(), m.value(), static_cast<ulong>(low));
    if (fmpz_is_zero(m.value())) fmpz_one(m.value());
    if (negative) fmpz_neg(m.value(), m.value());
    arf_set_fmpz_2exp(x, m.value(), e);
}
} // namespace

int main() {
    try {
        uint64_t checks = 0;
        d::Float x;
        std::mt19937_64 rng(20261027);
        d::Integer e;
        for (int round = 0; round < 20000; ++round) {
            const int limbs = 1 + static_cast<int>(rng() % 4);
            const int low = static_cast<int>(rng() % (static_cast<unsigned>(limbs) * FLINT_BITS));
            switch (rng() % 4) {
            case 0: fmpz_set_si(e.value(), static_cast<slong>(rng() % 200) - 100); break;
            case 1: fmpz_set_si(e.value(), static_cast<slong>(rng() % 2000000) - 1000000); break;
            case 2: fmpz_set_si(e.value(), static_cast<slong>(rng() >> 2) * ((rng() & 1) ? 1 : -1) / 4); break;
            default: // an exponent beyond a word: FLINT's path
                fmpz_one(e.value());
                fmpz_mul_2exp(e.value(), e.value(), 64 + rng() % 64);
                if (rng() & 1) fmpz_neg(e.value(), e.value());
            }
            random_arf(x.value(), rng, limbs, low, e.value(), rng() & 1);
            checks += check(x.value(), "random arf, " + std::to_string(limbs) + " limbs, round " + std::to_string(round));
        }
        // Edges.
        arf_zero(x.value());
        checks += check(x.value(), "zero");
        for (slong k : {slong(-200), slong(-1), slong(0), slong(1), slong(63), slong(64), slong(127), slong(128)}) {
            arf_one(x.value());
            arf_mul_2exp_si(x.value(), x.value(), k);
            checks += check(x.value(), "2^" + std::to_string(k));
            arf_neg(x.value(), x.value());
            checks += check(x.value(), "-2^" + std::to_string(k));
        }
        for (int bits : {1, 2, 52, 53, 63, 64, 65, 101, 127, 128, 129, 192}) {
            d::Integer m;
            fmpz_one(m.value());
            fmpz_mul_2exp(m.value(), m.value(), static_cast<ulong>(bits));
            fmpz_sub_ui(m.value(), m.value(), 1); // all ones: 2^bits - 1, the largest of its size at 128
            arf_set_fmpz(x.value(), m.value());
            checks += check(x.value(), std::to_string(bits) + " ones");
            for (int lowest : {0, 1, 63, 64, 65}) {
                if (lowest >= bits) continue;
                fmpz_one(m.value());
                fmpz_mul_2exp(m.value(), m.value(), static_cast<ulong>(bits - 1));
                fmpz_setbit(m.value(), static_cast<ulong>(lowest));
                arf_set_fmpz(x.value(), m.value());
                checks += check(x.value(), std::to_string(bits) + " bits, lowest bit " + std::to_string(lowest));
            }
        }
        // Decimal lengths and the 19-digit chunks: 10^k - 1, 10^k + 1, 10^k + 10^(k/2) + 1 (zeros inside a chunk)
        // and 10^19 c + 1 (a chunk of zeros and a one), each of both signs, k up to 39 (past 128 bits).
        for (ulong k = 0; k <= 39; ++k) {
            d::Integer power, m;
            fmpz_ui_pow_ui(power.value(), 10, k);
            for (int variant = 0; variant < 4; ++variant) {
                fmpz_set(m.value(), power.value());
                if (variant == 0) fmpz_sub_ui(m.value(), m.value(), 1);
                else if (variant == 1) fmpz_add_ui(m.value(), m.value(), 1);
                else if (variant == 2) {
                    d::Integer half;
                    fmpz_ui_pow_ui(half.value(), 10, k / 2);
                    fmpz_add(m.value(), m.value(), half.value());
                    fmpz_add_ui(m.value(), m.value(), 1);
                } else {
                    fmpz_ui_pow_ui(m.value(), 10, 19);
                    fmpz_mul_ui(m.value(), m.value(), k + 1);
                    fmpz_add_ui(m.value(), m.value(), 1);
                }
                if (fmpz_is_zero(m.value())) continue;
                for (bool negative : {false, true}) {
                    arf_set_fmpz(x.value(), m.value());
                    if (negative) arf_neg(x.value(), x.value());
                    checks += check(x.value(), "decimal edge " + std::to_string(k) + "/" + std::to_string(variant) +
                                                   (negative ? ", negative" : ""));
                }
            }
        }
        for (slong edge : {COEFF_MAX, COEFF_MIN, COEFF_MAX - 200, COEFF_MIN + 200}) {
            fmpz_set_si(e.value(), edge);
            arf_one(x.value());
            arf_mul_2exp_fmpz(x.value(), x.value(), e.value());
            checks += check(x.value(), "exponent " + std::to_string(edge));
        }
        // Copy, move, negation, assignment.
        arf_set_si(x.value(), -12345);
        arf_mul_2exp_si(x.value(), x.value(), 70);
        d::Endpoint a(x.value());
        d::Endpoint b = a;
        require(arf_equal(a.value(), b.value()), "a copy differs");
        b.negate();
        arf_neg(x.value(), x.value());
        require(arf_equal(b.value(), x.value()) && b.text().mantissa == "12345" && b.text().exponent == "70",
                "negation: " + b.text().mantissa + "@" + b.text().exponent);
        d::Endpoint c = std::move(b);
        require(c.text().mantissa == "12345" && d::Endpoint{}.text().mantissa == "0", "move or default");
        a = c;
        require(arf_equal(a.value(), c.value()), "assignment");
        checks += 4;
        std::cout << "endpoint text: " << checks << " checks ok\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "test_endpoint_text: " << error.what() << '\n';
        return 1;
    }
}
