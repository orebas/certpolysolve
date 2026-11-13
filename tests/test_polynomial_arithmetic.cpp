// Exact polynomial arithmetic on the input object (the 0.5 spec's C39): every operation of
// the C API against an independent exact reference (Boost rationals, exact_reference_api.hpp), over the Gaussian
// rationals; inputs written as integers, fractions, decimals and dyadics; the errors (a ball operand, a bad number,
// a degree beyond addressing, null arguments: *out stays NULL); and wilk30 built by from_roots and solved.
#include "exact_reference_api.hpp"
#include "internal.hpp"
#include <certpolysolve/certpolysolve.h>
#include <certpolysolve/certpolysolve.hpp>
#include <flint/fmpq.h>
#include <iostream>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

using namespace reference;

namespace {
using Owned = std::unique_ptr<cps_polynomial, decltype(&cps_polynomial_destroy)>;
using ComplexPolynomial = std::vector<Complex>;

Owned owned(cps_polynomial *p) { return Owned(p, cps_polynomial_destroy); }

std::string text(const Rational &q) {
    return boost::multiprecision::numerator(q).str() + "/" + boost::multiprecision::denominator(q).str();
}
Owned make(const ComplexPolynomial &c) {
    cps_polynomial *p = nullptr;
    require(cps_polynomial_create(&p) == CPS_OK, "create failed");
    for (size_t k = 0; k < c.size(); ++k)
        require(cps_polynomial_set_coefficient(p, k, text(c[k].real).c_str(), text(c[k].imag).c_str()) == CPS_OK,
                "set failed");
    return owned(p);
}
Rational rational(const certpolysolve::detail::Number &n) {
    fmpq_t q;
    fmpq_init(q);
    uint64_t budget = UINT64_MAX;
    require(n.rational_with_budget(q, budget, true), "a result coefficient does not expand");
    char *num = fmpz_get_str(nullptr, 10, fmpq_numref(q));
    char *den = fmpz_get_str(nullptr, 10, fmpq_denref(q));
    const Rational r{Integer{std::string(num)}, Integer{std::string(den)}};
    flint_free(num);
    flint_free(den);
    fmpq_clear(q);
    return r;
}
// A result's coefficients; the object keeps nonzero terms only (C40), and every coefficient is exact.
ComplexPolynomial read(const cps_polynomial *p) {
    ComplexPolynomial c;
    for (const auto &[e, k] : p->terms) {
        require(k.exact() && !k.zero(), "a result term is a ball or an exact zero");
        c.resize(static_cast<size_t>(e) + 1);
        c[static_cast<size_t>(e)] = {rational(k.real), rational(k.imag)};
    }
    return c;
}
ComplexPolynomial trimmed(ComplexPolynomial c) {
    while (!c.empty() && zero(c.back())) c.pop_back();
    return c;
}
ComplexPolynomial add(const ComplexPolynomial &a, const ComplexPolynomial &b, int sign = 1) {
    ComplexPolynomial r(std::max(a.size(), b.size()));
    for (size_t k = 0; k < r.size(); ++k) {
        const Complex x = k < a.size() ? a[k] : Complex(), y = k < b.size() ? b[k] : Complex();
        r[k] = {x.real + sign * y.real, x.imag + sign * y.imag};
    }
    return trimmed(r);
}
ComplexPolynomial multiply(const ComplexPolynomial &a, const ComplexPolynomial &b) {
    if (a.empty() || b.empty()) return {};
    ComplexPolynomial r(a.size() + b.size() - 1);
    for (size_t i = 0; i < a.size(); ++i)
        for (size_t j = 0; j < b.size(); ++j) r[i + j] = r[i + j] + a[i] * b[j];
    return trimmed(r);
}
bool equal(const ComplexPolynomial &a, const ComplexPolynomial &b) {
    if (a.size() != b.size()) return false;
    for (size_t k = 0; k < a.size(); ++k)
        if (a[k].real != b[k].real || a[k].imag != b[k].imag) return false;
    return true;
}
template <class Call> ComplexPolynomial result_of(Call call, const char *what) {
    cps_polynomial *out = nullptr;
    const cps_error e = call(&out);
    require(e == CPS_OK && out, (std::string(what) + " failed: " + cps_last_error()).c_str());
    const auto keep = owned(out);
    return read(out);
}

void sums_products_powers() {
    // a = 3/7 - 2x + (1/3 + 5i/2) x^3, b = -1 + i x + 4 x^2: both complex, of different degrees
    const ComplexPolynomial a = {{Rational(3, 7)}, {-2}, {0}, {Rational(1, 3), Rational(5, 2)}};
    const ComplexPolynomial b = {{-1}, {0, 1}, {4}};
    const auto pa = make(a), pb = make(b);
    require(equal(result_of([&](cps_polynomial **o) { return cps_polynomial_add(pa.get(), pb.get(), o); }, "add"), add(a, b)),
            "a + b");
    require(equal(result_of([&](cps_polynomial **o) { return cps_polynomial_subtract(pa.get(), pb.get(), o); }, "subtract"),
                  add(a, b, -1)),
            "a - b");
    require(equal(result_of([&](cps_polynomial **o) { return cps_polynomial_multiply(pa.get(), pb.get(), o); }, "multiply"),
                  multiply(a, b)),
            "a b");
    require(result_of([&](cps_polynomial **o) { return cps_polynomial_subtract(pa.get(), pa.get(), o); }, "a - a").empty(),
            "a - a is not the zero polynomial");
    // powers: (a b)^0 = 1, a^5 against repeated products, the zero polynomial to a positive power
    require(equal(result_of([&](cps_polynomial **o) { return cps_polynomial_power(pa.get(), 0, o); }, "power 0"), {{1}}),
            "a^0");
    ComplexPolynomial a5 = {{1}};
    for (int k = 0; k < 5; ++k) a5 = multiply(a5, a);
    require(equal(result_of([&](cps_polynomial **o) { return cps_polynomial_power(pa.get(), 5, o); }, "power 5"), a5), "a^5");
    const auto zero_polynomial = make({});
    require(result_of([&](cps_polynomial **o) { return cps_polynomial_power(zero_polynomial.get(), 3, o); }, "0^3").empty(),
            "0^3");
    // scaling, by a fraction, by i, and by a decimal and a dyadic written as such
    require(equal(result_of([&](cps_polynomial **o) { return cps_polynomial_scale(pb.get(), "-1/3", nullptr, o); }, "scale"),
                  multiply(b, {{Rational(-1, 3)}})),
            "b times -1/3");
    require(equal(result_of([&](cps_polynomial **o) { return cps_polynomial_scale(pb.get(), "0", "1", o); }, "scale by i"),
                  multiply(b, {{0, 1}})),
            "b times i");
    require(equal(result_of([&](cps_polynomial **o) { return cps_polynomial_scale(pb.get(), "0.1", "3@-2", o); }, "scale"),
                  multiply(b, {{Rational(1, 10), Rational(3, 4)}})),
            "b times 0.1 + 3/4 i");
    // a copy is equal and independent
    cps_polynomial *copy = nullptr;
    require(cps_polynomial_copy(pa.get(), &copy) == CPS_OK, "copy failed");
    const auto keep = owned(copy);
    require(equal(read(copy), a), "copy");
}

void from_roots() {
    // wilk30: (x - 1)(x - 2)...(x - 30), against the reference product
    std::vector<std::string> text_roots;
    ComplexPolynomial w = {{1}};
    for (int k = 1; k <= 30; ++k) {
        text_roots.push_back(std::to_string(k));
        w = multiply(w, {{-k}, {1}});
    }
    std::vector<const char *> roots;
    for (const auto &t : text_roots) roots.push_back(t.c_str());
    const auto wilk = result_of(
        [&](cps_polynomial **o) { return cps_polynomial_from_roots(roots.data(), nullptr, nullptr, roots.size(), o); }, "wilk30");
    require(equal(wilk, w), "from_roots(1..30)");
    // multiplicities and complex roots: (x - 1/2)^3 (x - i)^2 (x + i)^2, and a root with multiplicity 0 left out
    const char *real[] = {"1/2", "0", "0", "7"};
    const char *imag[] = {nullptr, "1", "-1", "0"};
    const uint64_t multiplicity[] = {3, 2, 2, 0};
    ComplexPolynomial expected = {{1}};
    for (int k = 0; k < 3; ++k) expected = multiply(expected, {{Rational(-1, 2)}, {1}});
    for (int k = 0; k < 2; ++k) expected = multiply(expected, {{0, -1}, {1}});
    for (int k = 0; k < 2; ++k) expected = multiply(expected, {{0, 1}, {1}});
    require(equal(result_of([&](cps_polynomial **o) { return cps_polynomial_from_roots(real, imag, multiplicity, 4, o); },
                            "from_roots with multiplicities"),
                  expected),
            "from_roots with multiplicities and complex roots");
    require(equal(result_of([&](cps_polynomial **o) { return cps_polynomial_from_roots(nullptr, nullptr, nullptr, 0, o); },
                            "from_roots of nothing"),
                  {{1}}),
            "from_roots of no roots is 1");
    // wilk30, solved: its 30 real roots, exactly the integers 1..30
    cps_polynomial *built = nullptr;
    require(cps_polynomial_from_roots(roots.data(), nullptr, nullptr, roots.size(), &built) == CPS_OK, "from_roots");
    const auto keep = owned(built);
    cps_options *options = nullptr;
    require(cps_options_create(&options) == CPS_OK && cps_options_set_real_only(options, 1) == CPS_OK, "options");
    cps_result *result = nullptr;
    require(cps_solve(built, options, &result) == CPS_OK, "solve");
    require(cps_result_termination(result) == CPS_SUCCESS && cps_result_size(result) == 30, "wilk30 was not solved");
    for (size_t k = 0; k < 30; ++k) {
        cps_root_view v{};
        require(cps_result_get_root(result, k, &v) == CPS_OK, "root view");
        const Rational lo = endpoint({v.real_lower.mantissa, v.real_lower.exponent});
        const Rational hi = endpoint({v.real_upper.mantissa, v.real_upper.exponent});
        int inside = 0;
        for (int r = 1; r <= 30; ++r) inside += lo <= r && r <= hi;
        require(inside == 1 && v.count == 1, "a wilk30 region does not hold exactly one of 1..30");
    }
    cps_result_destroy(result);
    cps_options_destroy(options);
}

void errors() {
    cps_polynomial *ball = nullptr, *out = reinterpret_cast<cps_polynomial *>(&ball);
    require(cps_polynomial_create(&ball) == CPS_OK && cps_polynomial_set_ball_coefficient(ball, 0, "1", "0", "1e-9", "0") == CPS_OK,
            "ball");
    const auto keep = owned(ball);
    require(cps_polynomial_add(ball, ball, &out) == CPS_INVALID_ARGUMENT && !out, "a ball operand was accepted");
    const auto x = make({{0}, {1}});
    out = reinterpret_cast<cps_polynomial *>(&ball);
    require(cps_polynomial_scale(x.get(), "abc", nullptr, &out) == CPS_INVALID_ARGUMENT && !out, "a bad number was accepted");
    out = reinterpret_cast<cps_polynomial *>(&ball);
    require(cps_polynomial_power(x.get(), UINT64_C(1) << 63, &out) == CPS_OUT_OF_RANGE && !out,
            "a degree beyond addressing was accepted");
    require(cps_polynomial_add(x.get(), nullptr, &out) == CPS_INVALID_ARGUMENT && !out, "a null operand was accepted");
    require(cps_polynomial_multiply(x.get(), x.get(), nullptr) == CPS_INVALID_ARGUMENT, "a null output was accepted");
    const char *roots[] = {"1", nullptr};
    require(cps_polynomial_from_roots(roots, nullptr, nullptr, 2, &out) == CPS_INVALID_ARGUMENT && !out,
            "a null root was accepted");
}
// C40: the object keeps its nonzero terms. Any exponent up to 2^63 - 1 can be set; an exact zero removes a term; a
// product of few terms is formed term by term and stays few terms at any degree; a solve that cannot hold the dense
// coefficient list refuses it.
using Term = std::tuple<uint64_t, Rational, Rational>;
std::vector<Term> terms_of(const cps_polynomial *p) {
    std::vector<Term> t;
    for (const auto &[e, c] : p->terms) t.emplace_back(e, rational(c.real), rational(c.imag));
    return t;
}
void sparse() {
    constexpr uint64_t largest = (UINT64_C(1) << 63) - 1; // a degree is an int64_t
    cps_polynomial *raw = nullptr;
    require(cps_polynomial_create(&raw) == CPS_OK, "create failed");
    const auto p = owned(raw);
    require(cps_polynomial_set_coefficient(p.get(), largest, "1", "0") == CPS_OK &&
                cps_polynomial_degree(p.get()) == static_cast<int64_t>(largest) && p.get()->terms.size() == 1,
            "the largest exponent");
    require(cps_polynomial_set_coefficient(p.get(), largest + 1, "1", "0") == CPS_OUT_OF_RANGE &&
                cps_polynomial_degree(p.get()) == static_cast<int64_t>(largest),
            "an exponent beyond 2^63 - 1 was accepted");
    require(cps_polynomial_set_coefficient(p.get(), largest, "0", "0") == CPS_OK && cps_polynomial_degree(p.get()) == -1 &&
                p.get()->terms.empty(),
            "an exact zero did not remove the term");
    require(cps_polynomial_set_ball_coefficient(p.get(), 5, "0", "0", "1/8", "0") == CPS_OK &&
                cps_polynomial_degree(p.get()) == 5,
            "a ball around zero is not a term");
    require(cps_polynomial_set_coefficient_d(p.get(), 5, 0.0, 0.0) == CPS_OK && cps_polynomial_degree(p.get()) == -1,
            "a double zero did not remove the term");
    // (x^1000000000 + 1)^2 by power and by multiplication: three terms
    const uint64_t big = 1000000000;
    cps_polynomial *b = nullptr;
    require(cps_polynomial_create(&b) == CPS_OK && cps_polynomial_set_coefficient(b, big, "1", "0") == CPS_OK &&
                cps_polynomial_set_coefficient(b, 0, "1", "0") == CPS_OK,
            "binomial");
    const auto binomial = owned(b);
    const std::vector<Term> square{{0, 1, 0}, {big, 2, 0}, {2 * big, 1, 0}};
    cps_polynomial *out = nullptr;
    require(cps_polynomial_power(binomial.get(), 2, &out) == CPS_OK, "power of a binomial");
    const auto by_power = owned(out);
    require(terms_of(by_power.get()) == square, "(x^1000000000 + 1)^2 by power");
    require(cps_polynomial_multiply(binomial.get(), binomial.get(), &out) == CPS_OK, "product of binomials");
    const auto by_product = owned(out);
    require(terms_of(by_product.get()) == square, "(x^1000000000 + 1)^2 by multiplication");
    require(cps_polynomial_derivative(by_power.get(), &out) == CPS_OK, "derivative");
    const auto derivative = owned(out);
    require(terms_of(derivative.get()) == std::vector<Term>{{big - 1, 2 * big, 0}, {2 * big - 1, 2 * big, 0}},
            "the derivative's terms");
    require(cps_polynomial_subtract(by_power.get(), by_product.get(), &out) == CPS_OK, "difference");
    const auto difference = owned(out);
    require(difference.get()->terms.empty(), "equal powers do not cancel");
    // the term-by-term path against the dense reference: a sparse factor times a dense one, and a sparse complex power
    const ComplexPolynomial dense_factor = {{Rational(3, 7)}, {-2, 1}, {0}, {Rational(1, 3), Rational(5, 2)}, {4},
                                            {0, -1}, {Rational(-7, 2)}, {1}};
    ComplexPolynomial sparse_factor(41);
    sparse_factor[0] = {-2};
    sparse_factor[40] = {0, 1}; // i x^40 - 2
    const auto pd = make(dense_factor), ps = make(sparse_factor);
    require(equal(result_of([&](cps_polynomial **o) { return cps_polynomial_multiply(ps.get(), pd.get(), o); }, "mixed"),
                  multiply(sparse_factor, dense_factor)),
            "a sparse times a dense polynomial");
    ComplexPolynomial seventh = {{1}};
    for (int k = 0; k < 7; ++k) seventh = multiply(seventh, sparse_factor);
    require(equal(result_of([&](cps_polynomial **o) { return cps_polynomial_power(ps.get(), 7, o); }, "sparse power"),
                  seventh),
            "(i x^40 - 2)^7");
    // a polynomial of degree above 2^62 to the power 1 is itself, and to the power 2 beyond the largest degree
    cps_polynomial *top = nullptr;
    require(cps_polynomial_create(&top) == CPS_OK && cps_polynomial_set_coefficient(top, (UINT64_C(1) << 62) + 1, "1", "0") == CPS_OK &&
                cps_polynomial_set_coefficient(top, 0, "3", "0") == CPS_OK,
            "a degree above 2^62");
    const auto keep_top = owned(top);
    require(cps_polynomial_power(top, 1, &out) == CPS_OK, "the power 1 of a degree above 2^62 was refused");
    const auto same = owned(out);
    require(terms_of(same.get()) == terms_of(top), "x^(2^62 + 1) + 3 to the power 1");
    out = nullptr;
    require(cps_polynomial_power(top, 2, &out) == CPS_OUT_OF_RANGE && !out, "a square beyond 2^63 - 1 was accepted");
    // a solve expands the terms; 2^62 + 1 coefficients are refused, not attempted
    cps_polynomial *huge = nullptr;
    require(cps_polynomial_create(&huge) == CPS_OK && cps_polynomial_set_coefficient(huge, UINT64_C(1) << 62, "1", "0") == CPS_OK &&
                cps_polynomial_set_coefficient(huge, 0, "-1", "0") == CPS_OK,
            "huge");
    const auto keep_huge = owned(huge);
    cps_options *options = nullptr;
    require(cps_options_create(&options) == CPS_OK, "options");
    cps_result *result = nullptr;
    require(cps_solve(huge, options, &result) == CPS_OUT_OF_RANGE && !result, "a solve of degree 2^62 was attempted");
    cps_options_destroy(options);
}
} // namespace

int main() {
    certpolysolve::ThreadScope lifetime;
    try {
        sums_products_powers();
        from_roots();
        errors();
        sparse();
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
    std::cout << "polynomial arithmetic tests passed\n";
    return 0;
}
