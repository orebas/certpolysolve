// Exact numbers: Number (number.hpp), a rational significand times 2 or 10 to an integer exponent kept unexpanded,
// read from the literal grammar (Number::Number) or a double (Number::from_double); its ball at a precision
// (Number::enclose, and Coefficient::enclose, which adds the radii as errors); its expansion into one rational.
//
// The parser reads "p/q" (q != 0), "m@e" (m 2^e) or a decimal, its fractional digits moved into the base-10 exponent,
// and throws std::invalid_argument on anything else. enclose multiplies the significand's ball by 2^e exactly, or by a
// ball of 10^e. expanded_storage_bits adds |e| bits for radix 2 and 4|e| for radix 10 (4 >= log2 10) to the
// significand's; compact_expansion_fits admits a coefficient list whose summed expanded size exceeds the storage it owns
// (owned_storage_bits) by at most the caller's allowance for each nonzero number (expansion_fits); Number::rational_with_budget expands a value only when that size fits the caller's running
// balance, which it charges, and, with allow_exponent_expansion false, its owned storage too; else it returns false.
//
// State: none.
// Entry points: Number::Number from c_api.cpp, program.cpp and piecewise_proposals.*; from_double from c_api.cpp;
// enclose wherever a ball is formed (solver.cpp, piecewise_roots.cpp, among others); compact_expansion_fits from
// solver.cpp, exact_preprocessing.cpp, exact_complex.cpp; rational_with_budget from those, polynomial_arithmetic.cpp,
// program.cpp, exact_text.cpp, real_region.hpp and c_api.cpp, the last two also reading the storage counts.
// Design notes: docs/api.md (Numbers; Resource limits), docs/resource-policy.md (Compact exponent rule).
#include "core/number.hpp"
#include <algorithm>
#include <cmath>
#include <climits>
#include <limits>
#include <stdexcept>
#include <string>

namespace certpolysolve::detail {
namespace {
bool digit(char c) { return c >= '0' && c <= '9'; }
std::string_view trim(std::string_view s) {
    constexpr std::string_view spaces = " \t\r\n";
    const auto begin = s.find_first_not_of(spaces);
    if (begin == s.npos) return {};
    return s.substr(begin, s.find_last_not_of(spaces) - begin + 1);
}
void parse_integer(fmpz_t out, std::string_view text) {
    if (text.empty()) throw std::invalid_argument("expected an integer");
    const auto offset = (text.front() == '+' || text.front() == '-') ? 1u : 0u;
    if (offset == text.size() || !std::all_of(text.begin() + offset, text.end(), digit))
        throw std::invalid_argument("invalid integer: " + std::string(text.substr(0, 80)));
    if (text.front() == '+') text.remove_prefix(1);
    if (fmpz_set_str(out, std::string(text).c_str(), 10) != 0)
        throw std::invalid_argument("invalid integer");
}
}

Number::Number(std::string_view source) {
    source = trim(source);
    if (source.empty()) throw std::invalid_argument("coefficient must be nonempty");
    if (source.size() > static_cast<size_t>(WORD_MAX) || source.size() > std::string().max_size())
        throw std::out_of_range("coefficient literal exceeds addressable storage");
    const auto slash = source.find('/');
    const auto at = source.find('@');
    if (slash != source.npos) {
        parse_integer(fmpq_numref(significand.value()), source.substr(0, slash));
        parse_integer(fmpq_denref(significand.value()), source.substr(slash + 1));
        if (fmpz_is_zero(fmpq_denref(significand.value()))) throw std::invalid_argument("zero denominator");
        fmpq_canonicalise(significand.value());
    } else if (at != source.npos) {
        parse_integer(fmpq_numref(significand.value()), source.substr(0, at));
        parse_integer(exponent.value(), source.substr(at + 1));
    } else {
        radix = 10;
        const auto e = source.find_first_of("eE");
        auto mantissa = source.substr(0, e);
        if (e != source.npos) parse_integer(exponent.value(), source.substr(e + 1));
        const auto dot = mantissa.find('.');
        std::string digits(mantissa);
        if (dot != mantissa.npos) {
            const auto fractional_digits = mantissa.size() - dot - 1;
            digits.erase(dot, 1);
            fmpz_sub_ui(exponent.value(), exponent.value(), static_cast<ulong>(fractional_digits));
        }
        parse_integer(fmpq_numref(significand.value()), digits);
    }
    if (zero()) fmpz_zero(exponent.value());
}

Number Number::from_double(double x) {
    if (!std::isfinite(x)) throw std::invalid_argument("NaN and infinity are not polynomial coefficients");
    Number result;
    Float value;
    arf_set_d(value.value(), x);
    arf_get_fmpz_2exp(fmpq_numref(result.significand.value()), result.exponent.value(), value.value());
    return result;
}

void Number::enclose(arb_t result, slong precision) const {
    // An integer significand (every m@e coefficient): arb_set_round_fmpz, the ball arb_set_fmpq gives it without its
    // division. arb_set_fmpq is arb_fmpz_div_fmpz, which rounds x / 1 correctly toward zero (ARB_RND), as
    // arf_set_round_fmpz rounds x, and sets the same one-ulp radius where the rounding is inexact; its Newton branch
    // needs a denominator of more than half the precision's bits (FLINT 3.6, arb/div.c).
    if (fmpz_is_one(fmpq_denref(significand.value())))
        arb_set_round_fmpz(result, fmpq_numref(significand.value()), precision);
    else
        arb_set_fmpq(result, significand.value(), precision);
    if (zero() || fmpz_is_zero(exponent.value())) return;
    if (radix == 2) {
        arb_mul_2exp_fmpz(result, result, exponent.value());
    } else {
        RealBall power, base;
        arb_set_ui(base.value(), 10);
        arb_pow_fmpz(power.value(), base.value(), exponent.value(), precision);
        arb_mul(result, result, power.value(), precision);
    }
}

uint64_t Number::owned_storage_bits() const {
    return resources::add(sizeof(Number) * CHAR_BIT,
        resources::add(fmpz_bits(fmpq_numref(significand.value())),
          resources::add(fmpz_bits(fmpq_denref(significand.value())), fmpz_bits(exponent.value()))));
}
uint64_t Number::expanded_storage_bits() const {
    if (zero()) return 0;
    Integer magnitude;
    fmpz_abs(magnitude.value(), exponent.value());
    if (!fmpz_abs_fits_ui(magnitude.value())) return UINT64_MAX;
    const auto base = resources::add(fmpz_bits(fmpq_numref(significand.value())),
                                    fmpz_bits(fmpq_denref(significand.value())));
    return resources::add(base, resources::multiply(fmpz_get_ui(magnitude.value()), radix == 2 ? 1 : 4));
}
bool expansion_fits(uint64_t owned, uint64_t expanded, uint64_t numbers, uint64_t allowance) {
    return expanded != UINT64_MAX && expanded <= resources::add(owned, resources::multiply(numbers, allowance));
}
bool compact_expansion_fits(std::span<const Coefficient> coefficients, uint64_t allowance) {
    uint64_t owned = 0, expanded = 0, numbers = 0;
    for (const auto &c : coefficients)
        for (const auto *n : {&c.real, &c.imag, &c.real_radius, &c.imag_radius}) {
            owned = resources::add(owned, n->owned_storage_bits());
            expanded = resources::add(expanded, n->expanded_storage_bits());
            numbers += !n->zero();
        }
    return expansion_fits(owned, expanded, numbers, allowance);
}
bool Number::rational_with_budget(fmpq_t result, uint64_t& remaining_bits,
                                  bool allow_exponent_expansion) const {
    if (zero()) { fmpq_zero(result); return true; }
    const auto expanded = expanded_storage_bits();
    if (expanded == UINT64_MAX || expanded > remaining_bits ||
        (!allow_exponent_expansion && expanded > owned_storage_bits()) ||
        expanded / CHAR_BIT > static_cast<uint64_t>(WORD_MAX)) return false;
    Integer magnitude;
    fmpz_abs(magnitude.value(), exponent.value());
    if (!fmpz_abs_fits_ui(magnitude.value())) return false;
    const auto exp = static_cast<uint64_t>(fmpz_get_ui(magnitude.value()));
    const auto base_bits = static_cast<uint64_t>(fmpz_bits(fmpq_numref(significand.value())))
                         + static_cast<uint64_t>(fmpz_bits(fmpq_denref(significand.value())));
    const uint64_t factor = radix == 2 ? 1 : 4; // ceil(log2(10)), conservative.
    if (base_bits > remaining_bits || exp > (remaining_bits - base_bits) / factor) return false;
    remaining_bits -= base_bits + exp * factor;
    fmpq_set(result, significand.value());
    if (exp && radix == 2 && fmpz_is_one(fmpq_denref(significand.value()))) {
        // An integer mantissa times 2^exponent: fmpq_mul_2exp and fmpq_div_2exp move the powers of two between the
        // numerator and the denominator and leave a canonical rational canonical, so this is the rational the power
        // and the gcd below give, without either. (The design pass of 26 October 2026 found the power's allocation
        // two thirds of a degree-8 solve's exact preparation.)
        if (fmpz_sgn(exponent.value()) < 0)
            fmpq_div_2exp(result, result, static_cast<flint_bitcnt_t>(exp));
        else
            fmpq_mul_2exp(result, result, static_cast<flint_bitcnt_t>(exp));
        return true;
    }
    if (exp) {
        Integer power;
        fmpz_ui_pow_ui(power.value(), radix, static_cast<ulong>(exp));
        if (fmpz_sgn(exponent.value()) < 0)
            fmpz_mul(fmpq_denref(result), fmpq_denref(result), power.value());
        else
            fmpz_mul(fmpq_numref(result), fmpq_numref(result), power.value());
        fmpq_canonicalise(result);
    }
    return true;
}

void Coefficient::enclose(acb_t result, slong precision) const {
    real.enclose(acb_realref(result), precision);
    imag.enclose(acb_imagref(result), precision);
    RealBall radius;
    real_radius.enclose(radius.value(), precision);
    arb_add_error(acb_realref(result), radius.value());
    imag_radius.enclose(radius.value(), precision);
    arb_add_error(acb_imagref(result), radius.value());
}
} // namespace certpolysolve::detail
