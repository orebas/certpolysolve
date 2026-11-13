// Stage 3a's C additions for exact values (docs/api-spec-0.5.md; the C++ layer is built on them): root ends m 2^e
// compared, converted to binary64 with each rounding, formatted, written as fractions, their midpoints and widths;
// the names of terminations and stop reasons; reading a polynomial's degree, coefficients (exact text that sets the
// same value again) and terms. Public API only.
#include <certpolysolve/certpolysolve.h>
#include <cfloat>
#include <clocale>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace {
void require(bool condition, const std::string &message) {
    if (!condition) throw std::runtime_error(message);
}
cps_dyadic_view end(const char *m, const char *e) { return {m, e}; }
int compare(cps_dyadic_view a, cps_dyadic_view b) {
    int c = 7;
    require(cps_dyadic_compare(a, b, &c) == CPS_OK, std::string("compare failed: ") + cps_last_error());
    return c;
}
double to_double(cps_dyadic_view x, cps_rounding r) {
    double d = 0;
    require(cps_dyadic_to_double(x, r, &d) == CPS_OK, std::string("to_double failed: ") + cps_last_error());
    return d;
}
template <class Call> std::string text(Call call) {
    size_t length = 0;
    require(call(nullptr, 0, &length) == CPS_OK, std::string("length query failed: ") + cps_last_error());
    std::string s(length + 1, '\0');
    size_t again = 0;
    require(call(s.data(), s.size(), &again) == CPS_OK && again == length, "text call failed");
    s.resize(length);
    return s;
}

void dyadics() {
    require(compare(end("1", "0"), end("3", "-1")) == -1, "1 < 3/2");
    require(compare(end("2", "0"), end("1", "1")) == 0, "2 = 1 2^1");
    require(compare(end("-5", "-3"), end("-1", "0")) == 1, "-5/8 > -1");
    require(compare(end("1", "1000000000"), end("1", "999999999")) == 1, "huge exponents");
    // binary64 with each rounding: 3/2 exactly; 2^-1075, half the smallest subnormal, a tie that goes to even (0);
    // 2^2000 beyond the range
    require(to_double(end("3", "-1"), CPS_ROUND_NEAREST) == 1.5, "3/2");
    const double tiny = std::numeric_limits<double>::denorm_min();
    require(to_double(end("1", "-1075"), CPS_ROUND_DOWN) == 0 && to_double(end("1", "-1075"), CPS_ROUND_UP) == tiny &&
                to_double(end("1", "-1075"), CPS_ROUND_NEAREST) == 0,
            "2^-1075 rounded");
    require(std::isinf(to_double(end("1", "2000"), CPS_ROUND_NEAREST)) &&
                to_double(end("1", "2000"), CPS_ROUND_DOWN) == DBL_MAX && std::isinf(to_double(end("1", "2000"), CPS_ROUND_UP)),
            "2^2000 rounded");
    require(to_double(end("-1", "2000"), CPS_ROUND_UP) == -DBL_MAX, "-2^2000 rounded up");
    // 17 significant digits read back as the nearest binary64; a shorter format of 1/3's nearest dyadic
    for (const auto &[m, e] : std::vector<std::pair<const char *, const char *>>{{"1", "-1"}, {"6004799503160661", "-54"},
                                                                                {"-123456789", "-40"}, {"7", "300"}}) {
        const auto s = text([&](char *b, size_t n, size_t *l) { return cps_dyadic_format(end(m, e), 17, b, n, l); });
        require(std::strtod(s.c_str(), nullptr) == to_double(end(m, e), CPS_ROUND_NEAREST), "17 digits: " + s);
    }
    const auto third = text([&](char *b, size_t n, size_t *l) { return cps_dyadic_format(end("6004799503160661", "-54"), 5, b, n, l); });
    require(third == "0.33333", "5 digits of 1/3: " + third);
    // the layout and the rounding of C's "%.*g" for a double (and of C++'s std::format("{:.N}")): exact ties to even
    // (0.125 to 2 digits is 0.12, 2.5 to 1 digit is 2), trailing zeros dropped, the exponent form below 10^-4 and
    // from 10^digits, every precision from 1 to 17 on doubles of many sizes
    std::vector<double> values{0.125, 0.375, 2.5, 3.5, 0.5, 1.5, 9.5, 99.5, 0.95, 1048576, 0.0009765625, 6.103515625e-05,
                               -2, 123456789, 1e-30, 7e-30, 0.1, 1.0 / 3, 5480471259.847853, -1e300, 1e-310};
    for (int k = 1; k <= 200; ++k) values.push_back(std::ldexp(std::sin(k * 0.7) + 1.3, (k * 37) % 400 - 200));
    for (const double v : values) {
        int exponent = 0;
        const double fraction = std::frexp(v, &exponent); // v = f 2^e exactly, with f 2^53 an integer
        const auto mantissa = std::to_string(static_cast<long long>(std::ldexp(fraction, 53)));
        const auto power = std::to_string(exponent - 53);
        for (unsigned digits = 1; digits <= 17; ++digits) {
            const auto ours = text([&](char *b, size_t n, size_t *l) {
                return cps_dyadic_format(end(mantissa.c_str(), power.c_str()), digits, b, n, l);
            });
            char expected[64];
            std::snprintf(expected, sizeof expected, "%.*g", static_cast<int>(digits), v);
            require(ours == expected, "format of " + std::string(expected) + " to " + std::to_string(digits) +
                                          " digits: " + ours);
        }
    }
    // 2^-1 to 100 digits is exact; a value beyond MPFR's exponent range is refused, not written wrongly
    require(text([&](char *b, size_t n, size_t *l) { return cps_dyadic_format(end("1", "-1"), 100, b, n, l); }) == "0.5",
            "0.5 to 100 digits");
    // more digits than a double has, both layouts (the texts MPFR's "%.*RNg" writes in the C locale)
    const char *wide = "123456789012345678901234567890123456789012345678901";
    const auto fixed = text([&](char *b, size_t n, size_t *l) { return cps_dyadic_format(end(wide, "-60"), 43, b, n, l); });
    require(fixed == "107081695084215790682679598959835.4881438168", "43 digits of a 167-bit value: " + fixed);
    const std::string negative_wide = std::string("-") + wide;
    const auto scientific = text([&](char *b, size_t n, size_t *l) {
        return cps_dyadic_format(end(negative_wide.c_str(), "-300"), 30, b, n, l);
    });
    require(scientific == "-6.06060916187146340863956649125e-41", "30 digits of a small value: " + scientific);
    require(cps_dyadic_format(end("1", "4000000000"), 17, nullptr, 0, nullptr) == CPS_OUT_OF_RANGE &&
                cps_dyadic_format(end("1", "-4000000000"), 17, nullptr, 0, nullptr) == CPS_OUT_OF_RANGE,
            "a value beyond MPFR's exponent range");
    // exact fractions, lowest terms
    auto rational = [&](const char *m, const char *e) {
        return text([&](char *b, size_t n, size_t *l) { return cps_dyadic_rational(end(m, e), b, n, l); });
    };
    require(rational("3", "-2") == "3/4" && rational("6", "-2") == "3/2" && rational("5", "3") == "40" &&
                rational("0", "7") == "0" && rational("-1", "-1") == "-1/2",
            "fractions");
    // midpoint and width, exact
    require(text([&](char *b, size_t n, size_t *l) { return cps_dyadic_midpoint(end("1", "0"), end("2", "0"), b, n, l); }) == "3@-1",
            "midpoint of 1 and 2");
    require(text([&](char *b, size_t n, size_t *l) { return cps_dyadic_width(end("1", "-1"), end("3", "0"), b, n, l); }) == "5@-1",
            "width from 1/2 to 3");
    // truncation: at most size bytes with the NUL, the full length reported
    char small[4] = {'x', 'x', 'x', 'x'};
    size_t length = 0;
    require(cps_dyadic_rational(end("12345", "0"), small, sizeof small, &length) == CPS_OK && length == 5 &&
                std::string(small) == "123",
            "truncation");
    // errors
    int c = 0;
    require(cps_dyadic_compare(end("1.5", "0"), end("1", "0"), &c) == CPS_INVALID_ARGUMENT, "a decimal point accepted");
    require(cps_dyadic_compare(end(nullptr, "0"), end("1", "0"), &c) == CPS_INVALID_ARGUMENT, "a null mantissa accepted");
    require(cps_dyadic_compare(end("1", ""), end("1", "0"), &c) == CPS_INVALID_ARGUMENT, "an empty exponent accepted");
    // An out-of-range rounding is refused: tested from C (test_c_contract.c, enum_arguments), where passing 9 is
    // legal; in C++ the cast to the enum type is itself undefined for a value outside its range.
    require(cps_dyadic_format(end("1", "0"), 0, small, sizeof small, &length) == CPS_INVALID_ARGUMENT, "0 digits accepted");
    require(cps_dyadic_rational(end("1", "0"), nullptr, 5, &length) == CPS_INVALID_ARGUMENT, "a null buffer accepted");
}

// The host's locale does not change the text. cps_dyadic_format once wrote with MPFR's printf, which takes the
// decimal point from localeconv(): "1,5" under a German locale. Checked under the first comma-decimal locale the
// machine has (a locale built by localedef can be named with LOCPATH); with none installed the check says it is
// skipped.
void locale_free() {
    const std::vector<std::tuple<const char *, const char *, unsigned>> cases{
        {"3", "-1", 5}, {"-123456789", "-40", 4}, {"7", "300", 3}, {"6004799503160661", "-54", 17}, {"1", "-1", 100}};
    std::vector<std::string> in_c;
    for (const auto &[m, e, digits] : cases)
        in_c.push_back(text([&](char *b, size_t n, size_t *l) { return cps_dyadic_format(end(m, e), digits, b, n, l); }));
    const char *used = nullptr;
    for (const char *name : {"de_DE.UTF-8", "de_DE.utf8", "fr_FR.UTF-8", "fr_FR.utf8", "es_ES.UTF-8", "it_IT.UTF-8",
                             "nl_NL.UTF-8", "ru_RU.UTF-8", "pt_BR.UTF-8"}) {
        if (std::setlocale(LC_ALL, name) && std::localeconv()->decimal_point[0] == ',') {
            used = name;
            break;
        }
    }
    if (!used) {
        std::setlocale(LC_ALL, "C");
        std::cout << "no comma-decimal locale installed: the locale check is skipped\n";
        return;
    }
    std::vector<std::string> in_comma;
    for (const auto &[m, e, digits] : cases)
        in_comma.push_back(text([&](char *b, size_t n, size_t *l) { return cps_dyadic_format(end(m, e), digits, b, n, l); }));
    std::setlocale(LC_ALL, "C");
    for (size_t i = 0; i < cases.size(); ++i)
        require(in_comma[i] == in_c[i], std::string("under ") + used + ": " + in_comma[i] + " for " + in_c[i]);
    require(in_c[0] == "1.5", "3/2 to 5 digits: " + in_c[0]);
    std::cout << "the same text under " << used << " as under the C locale\n";
}

void names() {
    // The name of a value outside the enumeration ("unknown") is tested from C (test_c_contract.c, enum_arguments):
    // in C++ the cast to the enum type is itself undefined for such a value.
    require(std::string(cps_termination_name(CPS_SUCCESS)) == "success" &&
                std::string(cps_termination_name(CPS_PRECISION_LIMIT)) == "precision limit",
            "termination names");
    require(std::string(cps_stop_reason_name(CPS_STOP_TIME_LIMIT)) == "time limit", "stop reason names");
}

void polynomial_reading() {
    cps_polynomial *p = nullptr;
    require(cps_polynomial_create(&p) == CPS_OK, "create");
    auto term_count = [&] {
        size_t count = 99;
        require(cps_polynomial_terms(p, nullptr, 0, &count) == CPS_OK, "term count");
        return count;
    };
    require(cps_polynomial_degree(p) == -1 && term_count() == 0, "the zero polynomial");
    // -1 + (1/3 - 2i) x^7 + 0.1 x^100 + (3/4 +/- 1e-9) x^101, and a coefficient set to zero at x^200
    require(cps_polynomial_set_coefficient(p, 0, "-1", "0") == CPS_OK && cps_polynomial_set_coefficient(p, 7, "1/3", "-2") == CPS_OK &&
                cps_polynomial_set_coefficient(p, 100, "0.1", "0") == CPS_OK &&
                cps_polynomial_set_ball_coefficient(p, 101, "3@-2", "0", "1e-9", "0") == CPS_OK &&
                cps_polynomial_set_coefficient(p, 200, "0", "0") == CPS_OK,
            "set");
    require(cps_polynomial_degree(p) == 101, "degree");
    require(term_count() == 4, "term count");
    std::vector<uint64_t> exponents(4);
    size_t count = 0;
    require(cps_polynomial_terms(p, exponents.data(), exponents.size(), &count) == CPS_OK && count == 4 &&
                exponents == std::vector<uint64_t>{0, 7, 100, 101},
            "term exponents");
    // a smaller array gets the first exponents and the full count
    std::vector<uint64_t> first(2, 99);
    require(cps_polynomial_terms(p, first.data(), first.size(), &count) == CPS_OK && count == 4 &&
                first == std::vector<uint64_t>{0, 7},
            "the first two exponents");
    require(cps_polynomial_terms(p, nullptr, 3, &count) == CPS_INVALID_ARGUMENT, "a null array with a capacity");
    auto part = [&](uint64_t index, cps_coefficient_part which) {
        return text([&](char *b, size_t n, size_t *l) { return cps_polynomial_get_coefficient(p, index, which, b, n, l); });
    };
    require(part(7, CPS_PART_REAL) == "1/3" && part(7, CPS_PART_IMAG) == "-2" && part(100, CPS_PART_REAL) == "1/10" &&
                part(101, CPS_PART_REAL) == "3/4" && part(5, CPS_PART_REAL) == "0" && part(1000, CPS_PART_IMAG) == "0",
            "coefficient text");
    // the text sets the same coefficient again: the radius included
    cps_polynomial *q = nullptr;
    require(cps_polynomial_create(&q) == CPS_OK, "create");
    for (uint64_t index : {0, 7, 100, 101})
        require(cps_polynomial_set_ball_coefficient(q, index, part(index, CPS_PART_REAL).c_str(), part(index, CPS_PART_IMAG).c_str(),
                                                   part(index, CPS_PART_REAL_RADIUS).c_str(),
                                                   part(index, CPS_PART_IMAG_RADIUS).c_str()) == CPS_OK,
                "set from text");
    for (uint64_t index : {0, 7, 100, 101})
        for (auto which : {CPS_PART_REAL, CPS_PART_IMAG, CPS_PART_REAL_RADIUS, CPS_PART_IMAG_RADIUS})
            require(part(index, which) ==
                        text([&](char *b, size_t n, size_t *l) { return cps_polynomial_get_coefficient(q, index, which, b, n, l); }),
                    "a coefficient did not round-trip");
    size_t unused = 0;
    require(cps_polynomial_get_coefficient(p, 0, static_cast<cps_coefficient_part>(7), nullptr, 0, &unused) == CPS_INVALID_ARGUMENT,
            "a bad part accepted");
    cps_polynomial_destroy(q);
    cps_polynomial_destroy(p);
}
} // namespace

int main() {
    try {
        dyadics();
        locale_free();
        names();
        polynomial_reading();
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        cps_cleanup_thread();
        return 1;
    }
    cps_cleanup_thread();
    std::cout << "exact text tests passed\n";
    return 0;
}
