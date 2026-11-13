// Exact values as text for the C API: parse_dyadic reads a root end m 2^e from its two decimal integers (an optional
// '-' and digits each); dyadic_format writes a dyadic to at most n significant digits, dyadic_rational exactly as
// "p/q" or "p", dyadic_text as "m@e" with m odd; number_text writes a coefficient part exactly, in the literal grammar.
//
// dyadic_format sets x exactly in an MPFR number with as many bits as m and writes it as C's %g does (g_format): the
// significant digits from mpfr_get_str, correctly rounded to nearest with ties to even, laid out here, so that the
// host's locale cannot change the decimal point (MPFR's printf takes it from localeconv()); an exponent beyond MPFR's
// range, or a bit count or exponent that a 32-bit long cannot hold (mpfr_precision, mpfr_exponent), throws
// std::out_of_range. dyadic_rational's m/2^-e is in lowest terms
// because m is odd. number_text writes a binary Number with an integer significand as "p/q" or "p" while |e| <= 64
// (readable_exponent_bits, m odd) and as "m@e" beyond (dyadic_readable), one with exponent 0 as its significand, and
// any other as its expanded rational (Number::rational_with_budget): each text parses back to the same value.
//
// State: none; MpfrFlagsKept leaves MPFR's exception flags, the calling thread's, as dyadic_format found them.
// Entry points: all five, called only from c_api.cpp: parse_dyadic by the six exact-end calls (cps_dyadic_compare,
// cps_dyadic_to_double, cps_dyadic_format, cps_dyadic_rational, cps_dyadic_midpoint, cps_dyadic_width), dyadic_format and
// dyadic_rational by their namesakes, dyadic_text by the last two, number_text by cps_polynomial_get_coefficient.
// Not here: parsing a literal (number.cpp). Design notes: docs/api.md (Exact ends; Reading a polynomial).
#include "numeric/arithmetic/exact_text.hpp"
#include <flint/arb.h>
#include <flint/fmpq.h>
#include <mpfr.h>
#include <algorithm>
#include <charconv>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string_view>

namespace certpolysolve::detail {
namespace {
// An optional '-' and one or more decimal digits.
bool integer_text(const char *s) {
    if (!s) return false;
    if (*s == '-') ++s;
    if (!*s) return false;
    for (; *s; ++s)
        if (*s < '0' || *s > '9') return false;
    return true;
}
std::string flint_text(char *s) {
    std::string t(s);
    flint_free(s);
    return t;
}
std::string fmpz_decimal(const fmpz_t z) { return flint_text(fmpz_get_str(nullptr, 10, z)); }
std::string fmpq_decimal(const fmpq_t q) { return flint_text(fmpq_get_str(nullptr, 10, q)); }

// GMP and MPFR values that free themselves, for the decimal conversion below.
class GmpInteger {
    mpz_t z_;

  public:
    explicit GmpInteger(const fmpz_t z) {
        mpz_init(z_);
        fmpz_get_mpz(z_, z);
    }
    ~GmpInteger() { mpz_clear(z_); }
    GmpInteger(const GmpInteger &) = delete;
    GmpInteger &operator=(const GmpInteger &) = delete;
    mpz_srcptr get() const { return z_; }
};
class MpfrNumber {
    mpfr_t v_;

  public:
    explicit MpfrNumber(mpfr_prec_t bits) { mpfr_init2(v_, bits); }
    ~MpfrNumber() { mpfr_clear(v_); }
    MpfrNumber(const MpfrNumber &) = delete;
    MpfrNumber &operator=(const MpfrNumber &) = delete;
    mpfr_ptr get() { return v_; }
};
// MPFR's exception flags are the calling thread's state: a conversion leaves them as it found them.
class MpfrFlagsKept {
    mpfr_flags_t saved_ = mpfr_flags_save();

  public:
    MpfrFlagsKept() = default;
    ~MpfrFlagsKept() { mpfr_flags_restore(saved_, MPFR_FLAGS_ALL); }
    MpfrFlagsKept(const MpfrFlagsKept &) = delete;
    MpfrFlagsKept &operator=(const MpfrFlagsKept &) = delete;
};
} // namespace

void parse_dyadic(arf_t out, cps_dyadic_view end) {
    if (!integer_text(end.mantissa) || !integer_text(end.exponent))
        throw std::invalid_argument("a root end needs a decimal integer mantissa and exponent");
    Integer m, e;
    fmpz_set_str(m.value(), end.mantissa, 10);
    fmpz_set_str(e.value(), end.exponent, 10);
    arf_set_fmpz_2exp(out, m.value(), e.value());
}
namespace {
// mpfr_prec_t and mpfr_exp_t are a long: 64 bits here, 32 on Windows (LLP64), where a count of bits or an exponent
// beyond them is refused rather than wrapped. Where the long has 64 bits the checks compile to nothing.
mpfr_prec_t mpfr_precision(ulong bits) {
    if constexpr (sizeof(mpfr_prec_t) < sizeof(ulong))
        if (bits > static_cast<ulong>(MPFR_PREC_MAX)) throw std::out_of_range("the value has more bits than MPFR's precision allows");
    return static_cast<mpfr_prec_t>(bits);
}
mpfr_exp_t mpfr_exponent(slong exponent) {
    if constexpr (sizeof(mpfr_exp_t) < sizeof(slong))
        if (exponent > std::numeric_limits<mpfr_exp_t>::max() || exponent < std::numeric_limits<mpfr_exp_t>::min())
            throw std::out_of_range("the value's binary exponent is beyond MPFR's range");
    return static_cast<mpfr_exp_t>(exponent);
}
// C's "%.*g" with precision `digits` (at least 1) for a nonzero value, without the C library's formatting.
// mpfr_get_str gives the first `digits` significant digits, correctly rounded to nearest with ties to even, and the
// point's place: value = 0.d1 d2 ... times 10^point. The exponent %e would write is x = point - 1. %g writes fixed
// notation when -4 <= x < digits, else one digit, the rest after the point and e+xx (two exponent digits at least);
// trailing zeros are dropped, and the point with them when nothing follows it.
std::string g_format(mpfr_srcptr value, unsigned digits) {
    mpfr_exp_t point = 0;
    char *raw = mpfr_get_str(nullptr, &point, 10, digits, value, MPFR_RNDN);
    if (!raw) throw std::bad_alloc();
    const std::unique_ptr<char, decltype(&mpfr_free_str)> owned(raw, mpfr_free_str);
    std::string_view significand(owned.get());
    std::string out;
    if (significand.front() == '-') {
        out += '-';
        significand.remove_prefix(1);
    }
    const auto trimmed = [](std::string_view s) {
        while (!s.empty() && s.back() == '0') s.remove_suffix(1);
        return s;
    };
    const auto fraction = [&out](std::string_view s) {
        if (s.empty()) return;
        out += '.';
        out.append(s);
    };
    const long long x = static_cast<long long>(point) - 1;
    if (x >= -4 && x < static_cast<long long>(digits)) {
        if (x >= 0) {
            const auto whole = static_cast<size_t>(x) + 1;
            out.append(significand.substr(0, whole));
            fraction(trimmed(significand.substr(whole)));
        } else {
            out += "0.";
            out.append(static_cast<size_t>(-x - 1), '0');
            out.append(trimmed(significand));
        }
        return out;
    }
    out += significand.front();
    fraction(trimmed(significand.substr(1)));
    out += x < 0 ? "e-" : "e+";
    char exponent[24];
    const char *written = std::to_chars(exponent, exponent + sizeof exponent, x < 0 ? -x : x).ptr;
    if (written - exponent < 2) out += '0';
    out.append(exponent, static_cast<size_t>(written - exponent));
    return out;
}
} // namespace
std::string dyadic_format(const arf_t x, unsigned digits) {
    if (arf_is_zero(x)) return "0";
    if (digits > static_cast<unsigned>(std::numeric_limits<int>::max()))
        throw std::out_of_range("more significant digits than a format precision can ask for");
    Integer m, e;
    arf_get_fmpz_2exp(m.value(), e.value(), x); // x = m 2^e
    if (!fmpz_fits_si(e.value())) throw std::out_of_range("the value's binary exponent is beyond MPFR's range");
    // x exactly in MPFR (as many bits as its mantissa), then C's %g, correctly rounded to nearest with ties to even,
    // which is also what C++'s std::format("{:.N}") does with a double.
    const MpfrFlagsKept flags;
    MpfrNumber value(std::max<mpfr_prec_t>(MPFR_PREC_MIN, mpfr_precision(fmpz_bits(m.value()))));
    // exact unless the exponent leaves MPFR's range (an overflow or an underflow then rounds it)
    const int rounded = mpfr_set_z_2exp(value.get(), GmpInteger(m.value()).get(), mpfr_exponent(fmpz_get_si(e.value())), MPFR_RNDN);
    if (rounded != 0 || !mpfr_regular_p(value.get()))
        throw std::out_of_range("the value's binary exponent is beyond MPFR's range");
    return g_format(value.get(), digits);
}
std::string dyadic_rational(const arf_t x) {
    if (arf_is_zero(x)) return "0";
    Integer m, e;
    arf_get_fmpz_2exp(m.value(), e.value(), x); // x = m 2^e, m odd
    Integer shift;
    fmpz_abs(shift.value(), e.value());
    if (!fmpz_abs_fits_ui(shift.value())) throw std::out_of_range("the exact value is too large to write");
    if (fmpz_sgn(e.value()) >= 0) {
        fmpz_mul_2exp(m.value(), m.value(), fmpz_get_ui(shift.value()));
        return fmpz_decimal(m.value());
    }
    Integer q;
    fmpz_one(q.value());
    fmpz_mul_2exp(q.value(), q.value(), fmpz_get_ui(shift.value()));
    return fmpz_decimal(m.value()) + "/" + fmpz_decimal(q.value());
}
std::string dyadic_text(const arf_t x) {
    if (arf_is_zero(x)) return "0";
    Integer m, e;
    arf_get_fmpz_2exp(m.value(), e.value(), x);
    return fmpz_decimal(m.value()) + "@" + fmpz_decimal(e.value());
}
// A dyadic m 2^e (m odd) as the fraction or integer it is while |e| is at most this many bits (a denominator of up to
// 2^64, twenty digits, still reads well), and as "m@e" past that, where the fraction would only grow.
constexpr slong readable_exponent_bits = 64;
std::string dyadic_readable(const arf_t x) {
    Integer m, e;
    arf_get_fmpz_2exp(m.value(), e.value(), x);
    if (fmpz_cmp_si(e.value(), -readable_exponent_bits) >= 0 && fmpz_cmp_si(e.value(), readable_exponent_bits) <= 0)
        return dyadic_rational(x);
    return fmpz_decimal(m.value()) + "@" + fmpz_decimal(e.value());
}
std::string number_text(const Number &n) {
    if (n.zero()) return "0";
    const fmpq *s = n.significand.value();
    if (n.radix == 2 && fmpz_is_one(fmpq_denref(s))) {
        Float x;
        arf_set_fmpz_2exp(x.value(), fmpq_numref(s), n.exponent.value());
        return dyadic_readable(x.value());
    }
    if (fmpz_is_zero(n.exponent.value())) return fmpq_decimal(s);
    // A decimal exponent, or a fraction with a binary one: the exact value, expanded.
    Rational q;
    uint64_t budget = UINT64_MAX;
    if (!n.rational_with_budget(q.value(), budget, true))
        throw std::out_of_range("the exact value is too large to write");
    return fmpq_decimal(q.value());
}
} // namespace certpolysolve::detail
