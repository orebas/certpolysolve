// Exact arithmetic on the input object's terms over the Gaussian rationals, behind the C API's polynomial functions
// (polynomial_add, polynomial_subtract, polynomial_multiply, polynomial_scale, polynomial_power, polynomial_derivative,
// polynomial_from_roots), and dense, the terms expanded into a coefficient list (declared in polynomial_input.hpp).
//
// Every coefficient, factor and root is read exactly (number: a ball is std::invalid_argument), an operand as Terms,
// its nonzero terms by exponent (from_input), and a result is written back without zero terms (to_input). A sum adds
// term by term and drops what cancels (accumulate). A product is formed term by term when it makes no more pairs than
// the dense result has coefficients (by_terms), and otherwise by FLINT's fmpq_poly_mul on the real and imaginary parts
// (dense_product); a power by square-and-multiply over product, or densely throughout when the operand's square is
// dense (dense_power); polynomial_from_roots raises each x - r to its multiplicity densely and multiplies the factors
// in a balanced tree. A degree above max_exponent, a dense length beyond WORD_MAX (require_addressable) and an
// exponent too large to expand exactly are std::out_of_range; dense throws ExpansionRefused for a degree a vector
// cannot address.
//
// State: none.
// Entry points: the polynomial_* functions, each called only by its C counterpart in c_api.cpp (cps_polynomial_add and
// the rest); dense, by session.cpp and real_only_route.cpp (dense_input, for a route that reads every coefficient).
// Not here: argument checks and the exceptions' cps_error codes (guarded, c_api.cpp); literals (Number, number.cpp).
// Design notes: docs/api-spec-0.5.md (section 10, items C39 and C40).
#include "numeric/arithmetic/polynomial_arithmetic.hpp"
#include <flint/fmpq_poly.h>
#include <algorithm>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace certpolysolve::detail {
namespace {
// An owned fmpq_poly_t (the raii.hpp pattern).
class RationalPolynomial {
    fmpq_poly_t value_;

  public:
    RationalPolynomial() { fmpq_poly_init(value_); }
    ~RationalPolynomial() { fmpq_poly_clear(value_); }
    RationalPolynomial(const RationalPolynomial &x) : RationalPolynomial() { fmpq_poly_set(value_, x.value_); }
    RationalPolynomial(RationalPolynomial &&x) noexcept : RationalPolynomial() { fmpq_poly_swap(value_, x.value_); }
    RationalPolynomial &operator=(RationalPolynomial x) noexcept {
        fmpq_poly_swap(value_, x.value_);
        return *this;
    }
    fmpq_poly_struct *value() & noexcept { return value_; }
    const fmpq_poly_struct *value() const & noexcept { return value_; }
    fmpq_poly_struct *value() && = delete;
    const fmpq_poly_struct *value() const && = delete;
};

// A dense polynomial over the Gaussian rationals, re + i im: FLINT's side of a product.
struct Gaussian {
    RationalPolynomial re, im;
    bool real() const { return fmpq_poly_is_zero(im.value()); }
    slong length() const { return std::max(fmpq_poly_length(re.value()), fmpq_poly_length(im.value())); }
};

// One exact coefficient re + i im, and a polynomial as its nonzero terms by exponent: the arithmetic's own form.
struct GaussianNumber {
    Rational re, im;
    bool zero() const { return fmpq_is_zero(re.value()) && fmpq_is_zero(im.value()); }
};
using Terms = std::map<uint64_t, GaussianNumber>;

// The largest degree a dense product can address: FLINT's lengths are slong.
void require_addressable(uint64_t degree) {
    if (degree >= static_cast<uint64_t>(WORD_MAX))
        throw std::out_of_range("the result's degree exceeds what a dense product can address");
}
uint64_t add_exponents(uint64_t a, uint64_t b) {
    if (a > max_exponent || b > max_exponent - a)
        throw std::out_of_range("the result's degree exceeds 2^63 - 1, the largest a polynomial can hold");
    return a + b;
}
// A number's exact value. The arithmetic expands exponents as written ("1@1000000" is a million-bit number), as an
// exact product must; an exponent beyond FLINT's word cannot be expanded at all.
void exact(fmpq_t out, const Number &n) {
    uint64_t budget = UINT64_MAX;
    if (!n.rational_with_budget(out, budget, true))
        throw std::out_of_range("a coefficient's exponent is too large to expand exactly");
}
GaussianNumber number(const Coefficient &c, uint64_t exponent) {
    if (!c.exact())
        throw std::invalid_argument("polynomial arithmetic needs exact coefficients; the coefficient of x^" +
                                    std::to_string(exponent) + " is a ball");
    GaussianNumber g;
    exact(g.re.value(), c.real);
    exact(g.im.value(), c.imag);
    return g;
}
Terms from_input(const cps_polynomial &p) {
    Terms t;
    for (const auto &[e, c] : p.terms) t.emplace_hint(t.end(), e, number(c, e));
    return t;
}
cps_polynomial to_input(Terms &&t) {
    cps_polynomial p;
    for (auto &[e, g] : t) {
        if (g.zero()) continue;
        Coefficient c;
        fmpq_swap(c.real.significand.value(), g.re.value());
        fmpq_swap(c.imag.significand.value(), g.im.value());
        p.terms.emplace_hint(p.terms.end(), e, std::move(c));
    }
    return p;
}
Terms one() {
    Terms t;
    GaussianNumber g;
    fmpq_one(g.re.value());
    t.emplace(0, std::move(g));
    return t;
}
uint64_t degree(const Terms &t) { return t.empty() ? 0 : t.rbegin()->first; }

// out[e] += g (or -= g); a term that cancels is removed.
void accumulate(Terms &out, uint64_t e, const GaussianNumber &g, bool subtract) {
    const auto it = out.try_emplace(e).first;
    auto &s = it->second;
    if (subtract) {
        fmpq_sub(s.re.value(), s.re.value(), g.re.value());
        fmpq_sub(s.im.value(), s.im.value(), g.im.value());
    } else {
        fmpq_add(s.re.value(), s.re.value(), g.re.value());
        fmpq_add(s.im.value(), s.im.value(), g.im.value());
    }
    if (s.zero()) out.erase(it);
}
// out = a b, (ar + i ai)(br + i bi) = (ar br - ai bi) + i (ar bi + ai br); out is neither a nor b.
void multiply(GaussianNumber &out, const GaussianNumber &a, const GaussianNumber &b, Rational &scratch) {
    fmpq_mul(out.re.value(), a.re.value(), b.re.value());
    fmpq_mul(scratch.value(), a.im.value(), b.im.value());
    fmpq_sub(out.re.value(), out.re.value(), scratch.value());
    fmpq_mul(out.im.value(), a.re.value(), b.im.value());
    fmpq_mul(scratch.value(), a.im.value(), b.re.value());
    fmpq_add(out.im.value(), out.im.value(), scratch.value());
}

Gaussian to_dense(const Terms &t) {
    Gaussian g;
    for (const auto &[e, c] : t) {
        if (!fmpq_is_zero(c.re.value())) fmpq_poly_set_coeff_fmpq(g.re.value(), static_cast<slong>(e), c.re.value());
        if (!fmpq_is_zero(c.im.value())) fmpq_poly_set_coeff_fmpq(g.im.value(), static_cast<slong>(e), c.im.value());
    }
    return g;
}
Terms from_dense(const Gaussian &g) {
    Terms t;
    Rational re, im;
    for (slong k = 0; k < g.length(); ++k) {
        fmpq_poly_get_coeff_fmpq(re.value(), g.re.value(), k);
        fmpq_poly_get_coeff_fmpq(im.value(), g.im.value(), k);
        if (fmpq_is_zero(re.value()) && fmpq_is_zero(im.value())) continue;
        GaussianNumber c;
        fmpq_swap(c.re.value(), re.value());
        fmpq_swap(c.im.value(), im.value());
        t.emplace_hint(t.end(), static_cast<uint64_t>(k), std::move(c));
    }
    return t;
}
Gaussian dense_product(const Gaussian &a, const Gaussian &b) {
    Gaussian r;
    fmpq_poly_mul(r.re.value(), a.re.value(), b.re.value());
    if (a.real() && b.real()) return r;
    RationalPolynomial t;
    fmpq_poly_mul(t.value(), a.im.value(), b.im.value());
    fmpq_poly_sub(r.re.value(), r.re.value(), t.value()); // re = ar br - ai bi
    fmpq_poly_mul(r.im.value(), a.re.value(), b.im.value());
    fmpq_poly_mul(t.value(), a.im.value(), b.re.value());
    fmpq_poly_add(r.im.value(), r.im.value(), t.value()); // im = ar bi + ai br
    return r;
}
Gaussian dense_one() {
    Gaussian r;
    fmpq_poly_one(r.re.value());
    return r;
}
Gaussian dense_power(const Gaussian &a, uint64_t exponent) {
    if (a.real()) {
        Gaussian r;
        fmpq_poly_pow(r.re.value(), a.re.value(), static_cast<ulong>(exponent));
        return r;
    }
    Gaussian result = dense_one(), base = a;
    for (uint64_t e = exponent; e; e >>= 1) {
        if (e & 1) result = dense_product(result, base);
        if (e > 1) base = dense_product(base, base);
    }
    return result;
}

// Term by term, a product pairs every term of a with every term of b; densely, FLINT forms all degree + 1
// coefficients. Term by term when it pairs no more terms than that (top is at most max_exponent, so top + 1 fits).
bool by_terms(size_t a_terms, size_t b_terms, uint64_t top) { return a_terms <= (top + 1) / b_terms; }
Terms product(const Terms &a, const Terms &b) {
    if (a.empty() || b.empty()) return {};
    const uint64_t top = add_exponents(degree(a), degree(b));
    if (!by_terms(a.size(), b.size(), top)) {
        require_addressable(top);
        return from_dense(dense_product(to_dense(a), to_dense(b)));
    }
    Terms r;
    GaussianNumber p;
    Rational scratch;
    for (const auto &[ea, ga] : a)
        for (const auto &[eb, gb] : b) {
            multiply(p, ga, gb, scratch);
            accumulate(r, ea + eb, p, false);
        }
    return r;
}
Terms power(const Terms &a, uint64_t exponent) {
    if (!exponent) return one();
    if (a.empty()) return {};
    const uint64_t d = degree(a);
    if (d && exponent > max_exponent / d)
        throw std::out_of_range("the result's degree exceeds 2^63 - 1, the largest a polynomial can hold");
    if (exponent == 1) return a;
    // From here exponent >= 2, so 2 d <= d exponent fits (the check above) and the square's degree can be formed.
    // A polynomial whose square is dense is raised densely (FLINT's power, or repeated dense squaring when complex).
    if (d && !by_terms(a.size(), a.size(), add_exponents(d, d))) {
        require_addressable(d * exponent);
        return from_dense(dense_power(to_dense(a), exponent));
    }
    Terms result = one(), base = a;
    for (uint64_t e = exponent; e; e >>= 1) {
        if (e & 1) result = product(result, base);
        if (e > 1) base = product(base, base);
    }
    return result;
}
Terms sum(Terms a, const Terms &b, bool subtract) {
    for (const auto &[e, g] : b) accumulate(a, e, g, subtract);
    return a;
}
} // namespace

DensePolynomial dense(const cps_polynomial &p) {
    DensePolynomial d;
    if (p.terms.empty()) return d;
    const uint64_t degree = p.terms.rbegin()->first;
    if (degree >= d.coefficients.max_size())
        throw ExpansionRefused("the polynomial's degree exceeds what a solve's coefficient list can address");
    d.coefficients.resize(static_cast<size_t>(degree) + 1);
    for (const auto &[exponent, c] : p.terms) d.coefficients[static_cast<size_t>(exponent)] = c;
    return d;
}

cps_polynomial polynomial_add(const cps_polynomial &a, const cps_polynomial &b) {
    return to_input(sum(from_input(a), from_input(b), false));
}
cps_polynomial polynomial_subtract(const cps_polynomial &a, const cps_polynomial &b) {
    return to_input(sum(from_input(a), from_input(b), true));
}
cps_polynomial polynomial_multiply(const cps_polynomial &a, const cps_polynomial &b) {
    return to_input(product(from_input(a), from_input(b)));
}
cps_polynomial polynomial_scale(const cps_polynomial &a, const Coefficient &factor) {
    const GaussianNumber f = number(factor, 0);
    Terms r;
    if (f.zero()) return to_input(std::move(r));
    Rational scratch;
    for (const auto &[e, g] : from_input(a)) {
        GaussianNumber p;
        multiply(p, g, f, scratch); // a nonzero term times a nonzero number is nonzero
        r.emplace_hint(r.end(), e, std::move(p));
    }
    return to_input(std::move(r));
}
cps_polynomial polynomial_power(const cps_polynomial &a, uint64_t exponent) {
    return to_input(power(from_input(a), exponent));
}
cps_polynomial polynomial_derivative(const cps_polynomial &a) {
    Terms r;
    for (const auto &[e, g] : from_input(a)) {
        if (!e) continue;
        GaussianNumber c;
        fmpq_mul_ui(c.re.value(), g.re.value(), static_cast<ulong>(e));
        fmpq_mul_ui(c.im.value(), g.im.value(), static_cast<ulong>(e));
        r.emplace_hint(r.end(), e - 1, std::move(c));
    }
    return to_input(std::move(r));
}
cps_polynomial polynomial_from_roots(std::span<const Coefficient> roots, std::span<const uint64_t> multiplicities) {
    uint64_t total = 0;
    for (const auto m : multiplicities) total = add_exponents(total, m);
    require_addressable(total);
    // The factors (x - r)^m, then a balanced product tree, so the work is that of a few products of half the size.
    std::vector<Gaussian> level;
    for (size_t k = 0; k < roots.size(); ++k) {
        if (!multiplicities[k]) continue;
        const GaussianNumber r = number(roots[k], 0);
        Gaussian linear;
        Rational negated;
        fmpq_neg(negated.value(), r.re.value());
        fmpq_poly_set_coeff_fmpq(linear.re.value(), 0, negated.value());
        fmpq_neg(negated.value(), r.im.value());
        fmpq_poly_set_coeff_fmpq(linear.im.value(), 0, negated.value());
        fmpq_poly_set_coeff_si(linear.re.value(), 1, 1);
        level.push_back(dense_power(linear, multiplicities[k]));
    }
    if (level.empty()) return to_input(one());
    while (level.size() > 1) {
        std::vector<Gaussian> next;
        for (size_t k = 0; k + 1 < level.size(); k += 2) next.push_back(dense_product(level[k], level[k + 1]));
        if (level.size() % 2) next.push_back(std::move(level.back()));
        level = std::move(next);
    }
    return to_input(from_dense(level.front()));
}
} // namespace certpolysolve::detail
