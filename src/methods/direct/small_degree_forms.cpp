// The closed forms of degree 3 for Auto's direct route: the contracts and the formulas are in small_degree_forms.hpp.
// The real forms and their cancellation-free pieces follow dev/docs/small-degree-closed-forms-20261026/cubic.md and its
// instruction-count probe (research/closed-forms/probe/closed_forms_probe.c); every formula here is an exact identity
// evaluated in balls, so a box holds its root whatever the rounding.
#include "methods/direct/small_degree_forms.hpp"
#include <algorithm>
#include <array>
#include <flint/fmpz_poly.h>
#include <flint/acb.h>
#include <flint/arb.h>

namespace certpolysolve::detail {
namespace {
// z = x y, z += k x y, z = k x on Gaussian integers.
void multiply(GaussianInteger &z, const GaussianInteger &x, const GaussianInteger &y) {
    Integer re, im;
    fmpz_mul(re.value(), x.re.value(), y.re.value());
    fmpz_submul(re.value(), x.im.value(), y.im.value());
    fmpz_mul(im.value(), x.re.value(), y.im.value());
    fmpz_addmul(im.value(), x.im.value(), y.re.value());
    z.re = std::move(re);
    z.im = std::move(im);
}
void add_multiple(GaussianInteger &z, const GaussianInteger &x, const GaussianInteger &y, slong k) {
    GaussianInteger t;
    multiply(t, x, y);
    fmpz_addmul_si(z.re.value(), t.re.value(), k);
    fmpz_addmul_si(z.im.value(), t.im.value(), k);
}
GaussianInteger scaled(const GaussianInteger &x, slong k) {
    GaussianInteger z;
    fmpz_mul_si(z.re.value(), x.re.value(), k);
    fmpz_mul_si(z.im.value(), x.im.value(), k);
    return z;
}
GaussianInteger product(const GaussianInteger &x, const GaussianInteger &y) {
    GaussianInteger z;
    multiply(z, x, y);
    return z;
}
void ball(acb_t out, const GaussianInteger &x, slong bits) {
    arb_set_round_fmpz(acb_realref(out), x.re.value(), bits);
    arb_set_round_fmpz(acb_imagref(out), x.im.value(), bits);
}
// A square root and a cube root whose balls do not straddle the cut: -z is rooted where Re z < 0 (as small_degree.cpp
// rotates the quadratic's discriminant); either branch serves the formulas below.
void square_root(acb_t out, const acb_t z, slong bits) {
    if (arb_is_negative(acb_realref(z))) {
        acb_neg(out, z);
        acb_sqrt(out, out, bits);
        acb_mul_onei(out, out);
    } else {
        acb_sqrt(out, z, bits);
    }
}
void cube_root(acb_t out, const acb_t z, slong bits) {
    if (arb_is_negative(acb_realref(z))) {
        acb_neg(out, z);
        acb_root_ui(out, out, 3, bits);
        acb_neg(out, out);
    } else {
        acb_root_ui(out, z, 3, bits);
    }
}
// The real Cardano pieces for integers D0, D1 and E = D1^2 - 4 D0^3 > 0: C = ((|D1| + sqrt E)/2)^(1/3) > 0,
// K = C + D0/C and L = C - D0/C, where the sum or difference would cancel written as |D1|/(C^2 - D0 + (D0/C)^2) and
// sqrt E/(C^2 + D0 + (D0/C)^2) (C^3 + D0^3/C^3 = |D1| and C^3 - D0^3/C^3 = sqrt E).
void real_cardano(arb_t K, arb_t L, const fmpz_t D0, const fmpz_t D1, const fmpz_t E, slong bits) {
    RealBall root_e, cube, d0, q, h;
    arb_set_round_fmpz(root_e.value(), E, bits);
    arb_sqrt(root_e.value(), root_e.value(), bits);
    arb_set_round_fmpz(cube.value(), D1, bits);
    arb_abs(cube.value(), cube.value());
    arb_add(cube.value(), cube.value(), root_e.value(), bits);
    arb_mul_2exp_si(cube.value(), cube.value(), -1);
    arb_root_ui(cube.value(), cube.value(), 3, bits);
    arb_set_round_fmpz(d0.value(), D0, bits);
    arb_div(q.value(), d0.value(), cube.value(), bits);
    const int d0_sign = fmpz_sgn(D0);
    if (d0_sign >= 0) {
        arb_add(K, cube.value(), q.value(), bits);
    } else {
        arb_sqr(h.value(), cube.value(), bits);
        arb_sub(h.value(), h.value(), d0.value(), bits);
        arb_addmul(h.value(), q.value(), q.value(), bits);
        arb_set_round_fmpz(K, D1, bits);
        arb_abs(K, K);
        arb_div(K, K, h.value(), bits);
    }
    if (d0_sign <= 0) {
        arb_sub(L, cube.value(), q.value(), bits);
    } else {
        arb_sqr(h.value(), cube.value(), bits);
        arb_add(h.value(), h.value(), d0.value(), bits);
        arb_addmul(h.value(), q.value(), q.value(), bits);
        arb_div(L, root_e.value(), h.value(), bits);
    }
}
Root simple_root(ComplexBall box, uint64_t multiplicity = 1) {
    Root root;
    root.box = std::move(box);
    root.count = root.multiplicity = multiplicity;
    return root;
}
// The roots of y^2 + B y + C for real balls B and C: two exactly real roots where the discriminant's ball is positive
// (the larger of -(B +- sqrt D)/2 in modulus, then C over it: no cancellation), a conjugate pair -B/2 +- i sqrt(-D)/2
// where it is negative; nullopt where its ball holds 0.
std::optional<std::array<ComplexBall, 2>> real_monic_quadratic(const arb_t B, const arb_t C, slong bits) {
    RealBall D, t;
    arb_sqr(D.value(), B, bits);
    arb_mul_2exp_si(t.value(), C, 2);
    arb_sub(D.value(), D.value(), t.value(), bits);
    std::array<ComplexBall, 2> out;
    if (arb_is_positive(D.value())) {
        RealBall root, plus, minus;
        arb_sqrt(root.value(), D.value(), bits);
        arb_add(plus.value(), B, root.value(), bits);
        arb_sub(minus.value(), B, root.value(), bits);
        Float plus_lower, minus_lower;
        arb_get_abs_lbound_arf(plus_lower.value(), plus.value(), bits);
        arb_get_abs_lbound_arf(minus_lower.value(), minus.value(), bits);
        if (arf_cmp(minus_lower.value(), plus_lower.value()) > 0) arb_swap(plus.value(), minus.value());
        if (arb_contains_zero(plus.value())) return std::nullopt;
        arb_neg(plus.value(), plus.value());
        arb_mul_2exp_si(plus.value(), plus.value(), -1);
        arb_set(acb_realref(out[0].value()), plus.value());
        arb_div(acb_realref(out[1].value()), C, plus.value(), bits);
        return out;
    }
    if (arb_is_negative(D.value())) {
        arb_neg(t.value(), D.value());
        arb_sqrt(t.value(), t.value(), bits);
        arb_mul_2exp_si(t.value(), t.value(), -1);
        arb_mul_2exp_si(acb_realref(out[0].value()), B, -1);
        arb_neg(acb_realref(out[0].value()), acb_realref(out[0].value()));
        arb_set(acb_imagref(out[0].value()), t.value());
        acb_conj(out[1].value(), out[0].value());
        return out;
    }
    return std::nullopt;
}
// The roots of y^2 + B y + C for complex balls, as small_degree.cpp's quadratic: the larger of -(B +- sqrt D)/2 in
// modulus, the square root rotated away from the cut, then C over it; nullopt where that root's ball holds 0.
std::optional<std::array<ComplexBall, 2>> complex_monic_quadratic(const acb_t B, const acb_t C, slong bits) {
    ComplexBall D, t, root, plus, minus;
    acb_sqr(D.value(), B, bits);
    acb_mul_2exp_si(t.value(), C, 2);
    acb_sub(D.value(), D.value(), t.value(), bits);
    square_root(root.value(), D.value(), bits);
    acb_add(plus.value(), B, root.value(), bits);
    acb_sub(minus.value(), B, root.value(), bits);
    Float plus_lower, minus_lower;
    acb_get_abs_lbound_arf(plus_lower.value(), plus.value(), bits);
    acb_get_abs_lbound_arf(minus_lower.value(), minus.value(), bits);
    if (arf_cmp(minus_lower.value(), plus_lower.value()) > 0) acb_swap(plus.value(), minus.value());
    if (acb_contains_zero(plus.value())) return std::nullopt;
    acb_neg(plus.value(), plus.value());
    acb_mul_2exp_si(plus.value(), plus.value(), -1);
    std::array<ComplexBall, 2> out;
    acb_set(out[0].value(), plus.value());
    acb_div(out[1].value(), C, plus.value(), bits);
    return out;
}
// The roots of the square-free real factor f2 y^2 + f1 y + f0 (integers), its discriminant's sign exact: two exactly
// real roots (the stable formula) or a conjugate pair.
std::array<ComplexBall, 2> exact_real_quadratic(const fmpz_t f0, const fmpz_t f1, const fmpz_t f2, slong bits) {
    Integer D, t;
    fmpz_mul(D.value(), f1, f1);
    fmpz_mul(t.value(), f2, f0);
    fmpz_submul_ui(D.value(), t.value(), 4);
    RealBall root, a, b, c;
    arb_set_round_fmpz(a.value(), f2, bits);
    arb_set_round_fmpz(b.value(), f1, bits);
    arb_set_round_fmpz(c.value(), f0, bits);
    std::array<ComplexBall, 2> out;
    if (fmpz_sgn(D.value()) > 0) {
        arb_set_round_fmpz(root.value(), D.value(), bits);
        arb_sqrt(root.value(), root.value(), bits);
        // Q = -(f1 + sign(f1) sqrt D)/2, the roots Q/f2 and f0/Q.
        if (fmpz_sgn(f1) < 0) arb_neg(root.value(), root.value());
        arb_add(root.value(), root.value(), b.value(), bits);
        arb_neg(root.value(), root.value());
        arb_mul_2exp_si(root.value(), root.value(), -1);
        arb_div(acb_realref(out[0].value()), root.value(), a.value(), bits);
        arb_div(acb_realref(out[1].value()), c.value(), root.value(), bits);
        return out;
    }
    // -f1/(2 f2) +- i sqrt(-D)/(2 f2)
    fmpz_neg(D.value(), D.value());
    arb_set_round_fmpz(root.value(), D.value(), bits);
    arb_sqrt(root.value(), root.value(), bits);
    arb_mul_2exp_si(a.value(), a.value(), 1);
    arb_div(acb_imagref(out[0].value()), root.value(), a.value(), bits);
    arb_div(acb_realref(out[0].value()), b.value(), a.value(), bits);
    arb_neg(acb_realref(out[0].value()), acb_realref(out[0].value()));
    acb_conj(out[1].value(), out[0].value());
    return out;
}
} // namespace

ExactCubic::ExactCubic(std::span<const GaussianInteger, 4> coefficients)
    : a_(coefficients[3]), b_(coefficients[2]), c_(coefficients[1]), d_(coefficients[0]) {
    real_ = a_.real() && b_.real() && c_.real() && d_.real();
    if (real_ && fmpz_sgn(a_.re.value()) < 0)
        for (auto *x : {&a_, &b_, &c_, &d_}) fmpz_neg(x->re.value(), x->re.value());
    // D0 = b^2 - 3ac, D1 = 2b^3 - 9abc + 27a^2 d, E = D1^2 - 4 D0^3.
    multiply(D0_, b_, b_);
    add_multiple(D0_, a_, c_, -3);
    multiply(D1_, product(b_, b_), b_);
    D1_ = scaled(D1_, 2);
    add_multiple(D1_, product(a_, b_), c_, -9);
    add_multiple(D1_, product(a_, a_), d_, 27);
    multiply(E_, D1_, D1_);
    add_multiple(E_, product(D0_, D0_), D0_, -4);
}

std::optional<std::vector<Root>> ExactCubic::roots(slong bits) const {
    std::vector<Root> out;
    ComplexBall a3, b;
    ball(a3.value(), a_, bits);
    acb_mul_ui(a3.value(), a3.value(), 3, bits);
    ball(b.value(), b_, bits);
    if (!square_free()) {
        // The exact multiple roots, each a quotient of Gaussian integers.
        const auto quotient = [&](const GaussianInteger &n, const GaussianInteger &m) {
            ComplexBall q, den;
            ball(q.value(), n, bits);
            ball(den.value(), m, bits);
            acb_div(q.value(), q.value(), den.value(), bits);
            return q;
        };
        if (D0_.zero()) {  // a triple root -b/(3a)
            out.push_back(simple_root(quotient(scaled(b_, -1), scaled(a_, 3)), 3));
            return out;
        }
        // The double root (9ad - bc)/(2 D0) and the simple root (4abc - 9a^2 d - b^3)/(a D0).
        GaussianInteger n2 = scaled(product(a_, d_), 9);
        add_multiple(n2, b_, c_, -1);
        GaussianInteger n1 = scaled(product(product(a_, b_), c_), 4);
        add_multiple(n1, product(a_, a_), d_, -9);
        add_multiple(n1, product(b_, b_), b_, -1);
        out.push_back(simple_root(quotient(n2, scaled(D0_, 2)), 2));
        out.push_back(simple_root(quotient(n1, product(a_, D0_))));
        return out;
    }
    if (real_ && fmpz_sgn(E_.re.value()) < 0) {
        // Three real roots: (2 sqrt(D0) cos((theta + 2 pi k)/3) - b)/(3a), theta = atan2(sqrt(-E), -D1). D0 > 0 here,
        // since 4 D0^3 = D1^2 - E > 0.
        RealBall y, x, theta, turn, radius, angle;
        arb_set_round_fmpz(y.value(), E_.re.value(), bits);
        arb_neg(y.value(), y.value());
        arb_sqrt(y.value(), y.value(), bits);
        arb_set_round_fmpz(x.value(), D1_.re.value(), bits);
        arb_neg(x.value(), x.value());
        arb_atan2(theta.value(), y.value(), x.value(), bits);
        arb_const_pi(turn.value(), bits);
        arb_mul_2exp_si(turn.value(), turn.value(), 1);
        arb_set_round_fmpz(radius.value(), D0_.re.value(), bits);
        arb_sqrt(radius.value(), radius.value(), bits);
        arb_mul_2exp_si(radius.value(), radius.value(), 1);
        for (slong k = 0; k < 3; ++k) {
            arb_set(angle.value(), theta.value());
            arb_addmul_si(angle.value(), turn.value(), k, bits);
            arb_div_ui(angle.value(), angle.value(), 3, bits);
            arb_cos(angle.value(), angle.value(), bits);
            arb_mul(angle.value(), angle.value(), radius.value(), bits);
            arb_sub(angle.value(), angle.value(), acb_realref(b.value()), bits);
            ComplexBall root;
            arb_div(acb_realref(root.value()), angle.value(), acb_realref(a3.value()), bits);
            out.push_back(simple_root(std::move(root)));
        }
        return out;
    }
    if (real_) {
        // One real root and a conjugate pair, from the real Cardano pieces K and L.
        RealBall K, L, h;
        real_cardano(K.value(), L.value(), D0_.re.value(), D1_.re.value(), E_.re.value(), bits);
        const slong d1_sign = fmpz_sgn(D1_.re.value()) >= 0 ? 1 : -1;
        ComplexBall real_root, pair;
        // (-b - sign(D1) K)/(3a)
        arb_mul_si(h.value(), K.value(), -d1_sign, bits);
        arb_sub(h.value(), h.value(), acb_realref(b.value()), bits);
        arb_div(acb_realref(real_root.value()), h.value(), acb_realref(a3.value()), bits);
        // (-b + sign(D1) K/2)/(3a) +- i sqrt(3) L/(6a)
        arb_mul_si(h.value(), K.value(), d1_sign, bits);
        arb_mul_2exp_si(h.value(), h.value(), -1);
        arb_sub(h.value(), h.value(), acb_realref(b.value()), bits);
        arb_div(acb_realref(pair.value()), h.value(), acb_realref(a3.value()), bits);
        arb_sqrt_ui(h.value(), 3, bits);
        arb_mul(h.value(), h.value(), L.value(), bits);
        arb_mul_2exp_si(h.value(), h.value(), -1);
        arb_div(acb_imagref(pair.value()), h.value(), acb_realref(a3.value()), bits);
        ComplexBall conjugate;
        acb_conj(conjugate.value(), pair.value());
        out.push_back(simple_root(std::move(real_root)));
        out.push_back(simple_root(std::move(pair)));
        out.push_back(simple_root(std::move(conjugate)));
        return out;
    }
    // Complex coefficients: C a cube root of the larger of (D1 +- sqrt E)/2, the roots -(b + xi^k C + D0/(xi^k C))/(3a).
    ComplexBall root_e, plus, minus, cube, d0, xi, term;
    ball(root_e.value(), E_, bits);
    square_root(root_e.value(), root_e.value(), bits);
    ball(plus.value(), D1_, bits);
    acb_sub(minus.value(), plus.value(), root_e.value(), bits);
    acb_add(plus.value(), plus.value(), root_e.value(), bits);
    Float plus_lower, minus_lower;
    acb_get_abs_lbound_arf(plus_lower.value(), plus.value(), bits);
    acb_get_abs_lbound_arf(minus_lower.value(), minus.value(), bits);
    if (arf_cmp(minus_lower.value(), plus_lower.value()) > 0) acb_swap(plus.value(), minus.value());
    if (acb_contains_zero(plus.value())) return std::nullopt;
    acb_mul_2exp_si(plus.value(), plus.value(), -1);
    cube_root(cube.value(), plus.value(), bits);
    if (acb_contains_zero(cube.value())) return std::nullopt;
    ball(d0.value(), D0_, bits);
    arb_set_si(acb_realref(xi.value()), -1);
    arb_sqrt_ui(acb_imagref(xi.value()), 3, bits);
    acb_mul_2exp_si(xi.value(), xi.value(), -1);
    for (int k = 0; k < 3; ++k) {
        if (k) acb_mul(cube.value(), cube.value(), xi.value(), bits);
        acb_div(term.value(), d0.value(), cube.value(), bits);
        acb_add(term.value(), term.value(), cube.value(), bits);
        acb_add(term.value(), term.value(), b.value(), bits);
        acb_neg(term.value(), term.value());
        ComplexBall root;
        acb_div(root.value(), term.value(), a3.value(), bits);
        out.push_back(simple_root(std::move(root)));
    }
    return out;
}
ExactQuartic::ExactQuartic(std::span<const GaussianInteger, 5> coefficients)
    : a_(coefficients[4]), b_(coefficients[3]), c_(coefficients[2]), d_(coefficients[1]), e_(coefficients[0]) {
    real_ = a_.real() && b_.real() && c_.real() && d_.real() && e_.real();
    if (real_ && fmpz_sgn(a_.re.value()) < 0)
        for (auto *x : {&a_, &b_, &c_, &d_, &e_}) fmpz_neg(x->re.value(), x->re.value());
    // D0 = c^2 - 3bd + 12ae, D1 = 2c^3 - 9bcd + 27b^2 e + 27ad^2 - 72ace, E = D1^2 - 4 D0^3.
    multiply(D0_, c_, c_);
    add_multiple(D0_, b_, d_, -3);
    add_multiple(D0_, a_, e_, 12);
    multiply(D1_, product(c_, c_), c_);
    D1_ = scaled(D1_, 2);
    add_multiple(D1_, product(b_, c_), d_, -9);
    add_multiple(D1_, product(b_, b_), e_, 27);
    add_multiple(D1_, product(a_, d_), d_, 27);
    add_multiple(D1_, product(a_, c_), e_, -72);
    multiply(E_, D1_, D1_);
    add_multiple(E_, product(D0_, D0_), D0_, -4);
    // P8 = 8ac - 3b^2, Q8 = b^3 - 4abc + 8a^2 d, R256 = -3b^4 + 16ab^2 c - 64a^2 bd + 256a^3 e.
    P8_ = scaled(product(a_, c_), 8);
    add_multiple(P8_, b_, b_, -3);
    multiply(Q8_, product(b_, b_), b_);
    add_multiple(Q8_, product(a_, b_), c_, -4);
    add_multiple(Q8_, product(a_, a_), d_, 8);
    const auto b2 = product(b_, b_);
    R256_ = scaled(product(b2, b2), -3);
    add_multiple(R256_, product(a_, b2), c_, 16);
    add_multiple(R256_, product(product(a_, a_), b_), d_, -64);
    add_multiple(R256_, product(product(a_, a_), a_), e_, 256);
    if (square_free()) return;
    if (!real_) {
        declined_ = true;
        return;
    }
    // The exact square-free factorization over Z: factors of degree 1 or 2 (a multiple root leaves at most two).
    IntegerPolynomial f;
    const GaussianInteger *parts[] = {&e_, &d_, &c_, &b_, &a_};
    for (slong k = 0; k < 5; ++k) fmpz_poly_set_coeff_fmpz(f.value(), k, parts[k]->re.value());
    struct Factorization {
        fmpz_poly_factor_t value;
        Factorization() { fmpz_poly_factor_init(value); }
        ~Factorization() { fmpz_poly_factor_clear(value); }
        Factorization(const Factorization &) = delete;
        Factorization &operator=(const Factorization &) = delete;
    } factorization;
    fmpz_poly_factor_squarefree(factorization.value, f.value());
    for (slong i = 0; i < factorization.value->num; ++i) {
        const fmpz_poly_struct *g = factorization.value->p + i;
        const slong degree = fmpz_poly_degree(g);
        if (degree < 1 || degree > 2) {
            declined_ = true;
            return;
        }
        Factor factor;
        factor.multiplicity = static_cast<uint64_t>(factorization.value->exp[i]);
        factor.coefficients.resize(static_cast<size_t>(degree) + 1);
        for (slong k = 0; k <= degree; ++k) fmpz_poly_get_coeff_fmpz(factor.coefficients[static_cast<size_t>(k)].re.value(), g, k);
        factors_.push_back(std::move(factor));
    }
}

std::optional<std::vector<Root>> ExactQuartic::roots(slong bits) const {
    std::vector<Root> out;
    if (!square_free()) {
        for (const auto &factor : factors_) {
            const auto &f = factor.coefficients;
            if (f.size() == 2) {  // -f0/f1
                ComplexBall root, den;
                ball(root.value(), f[0], bits);
                ball(den.value(), f[1], bits);
                acb_div(root.value(), root.value(), den.value(), bits);
                acb_neg(root.value(), root.value());
                out.push_back(simple_root(std::move(root), factor.multiplicity));
            } else {
                for (auto &root : exact_real_quadratic(f[0].re.value(), f[1].re.value(), f[2].re.value(), bits))
                    out.push_back(simple_root(std::move(root), factor.multiplicity));
            }
        }
        return out;
    }
    // y = x + b/(4a): p = P8/(8a^2), q = Q8/(8a^3), the shift -b/(4a).
    ComplexBall a, shift, p, q, t;
    ball(a.value(), a_, bits);
    ball(shift.value(), b_, bits);
    acb_div(shift.value(), shift.value(), a.value(), bits);
    acb_mul_2exp_si(shift.value(), shift.value(), -2);
    acb_neg(shift.value(), shift.value());
    ball(p.value(), P8_, bits);
    acb_sqr(t.value(), a.value(), bits);
    acb_mul_2exp_si(t.value(), t.value(), 3);
    acb_div(p.value(), p.value(), t.value(), bits);
    ball(q.value(), Q8_, bits);
    acb_mul(t.value(), t.value(), a.value(), bits);
    acb_div(q.value(), q.value(), t.value(), bits);
    std::vector<ComplexBall> ys;
    const auto take = [&](const std::optional<std::array<ComplexBall, 2>> &pair) {
        if (!pair) return false;
        for (const auto &y : *pair) ys.push_back(y);
        return true;
    };
    if (Q8_.zero()) {
        // The biquadratic: w^2 + p w + r with r = R256/(256 a^4), then y = +- sqrt(w).
        ComplexBall r;
        ball(r.value(), R256_, bits);
        acb_sqr(t.value(), a.value(), bits);
        acb_sqr(t.value(), t.value(), bits);
        acb_mul_2exp_si(t.value(), t.value(), 8);
        acb_div(r.value(), r.value(), t.value(), bits);
        std::optional<std::array<ComplexBall, 2>> w = real_ ? real_monic_quadratic(acb_realref(p.value()), acb_realref(r.value()), bits)
                                                            : complex_monic_quadratic(p.value(), r.value(), bits);
        if (!w) return std::nullopt;
        for (const auto &wk : *w) {
            ComplexBall y;
            if (real_ && arb_is_zero(acb_imagref(wk.value()))) {
                // A real w: real roots +- sqrt(w) where w > 0, the pair +- i sqrt(-w) where w < 0.
                if (arb_is_positive(acb_realref(wk.value()))) {
                    arb_sqrt(acb_realref(y.value()), acb_realref(wk.value()), bits);
                } else if (arb_is_negative(acb_realref(wk.value()))) {
                    arb_neg(acb_imagref(y.value()), acb_realref(wk.value()));
                    arb_sqrt(acb_imagref(y.value()), acb_imagref(y.value()), bits);
                } else {
                    return std::nullopt;
                }
            } else {
                square_root(y.value(), wk.value(), bits);
            }
            ComplexBall minus_y;
            acb_neg(minus_y.value(), y.value());
            ys.push_back(std::move(y));
            ys.push_back(std::move(minus_y));
        }
    } else if (real_) {
        // The largest real root of the resolvent: z = ((Q + D0/Q)/a - 2p)/3 with Q + D0/Q = 2 sqrt(D0) cos(theta/3),
        // theta = atan2(sqrt(-E), D1), where E < 0, and sign(D1) K where E > 0.
        RealBall sum, z, s, ratio, T, V;
        if (fmpz_sgn(E_.re.value()) < 0) {
            RealBall y, x, theta;
            arb_set_round_fmpz(y.value(), E_.re.value(), bits);
            arb_neg(y.value(), y.value());
            arb_sqrt(y.value(), y.value(), bits);
            arb_set_round_fmpz(x.value(), D1_.re.value(), bits);
            arb_atan2(theta.value(), y.value(), x.value(), bits);
            arb_div_ui(theta.value(), theta.value(), 3, bits);
            arb_cos(theta.value(), theta.value(), bits);
            arb_set_round_fmpz(sum.value(), D0_.re.value(), bits);
            arb_sqrt(sum.value(), sum.value(), bits);
            arb_mul_2exp_si(sum.value(), sum.value(), 1);
            arb_mul(sum.value(), sum.value(), theta.value(), bits);
        } else {
            RealBall L;
            real_cardano(sum.value(), L.value(), D0_.re.value(), D1_.re.value(), E_.re.value(), bits);
            if (fmpz_sgn(D1_.re.value()) < 0) arb_neg(sum.value(), sum.value());
        }
        arb_div(z.value(), sum.value(), acb_realref(a.value()), bits);
        arb_submul_ui(z.value(), acb_realref(p.value()), 2, bits);
        arb_div_ui(z.value(), z.value(), 3, bits);
        if (!arb_is_positive(z.value())) return std::nullopt;
        arb_sqrt(s.value(), z.value(), bits);
        // t, v = (p + z -+ q/s)/2; the factors y^2 + s y + t and y^2 - s y + v.
        arb_div(ratio.value(), acb_realref(q.value()), s.value(), bits);
        arb_add(T.value(), acb_realref(p.value()), z.value(), bits);
        arb_add(V.value(), T.value(), ratio.value(), bits);
        arb_sub(T.value(), T.value(), ratio.value(), bits);
        arb_mul_2exp_si(T.value(), T.value(), -1);
        arb_mul_2exp_si(V.value(), V.value(), -1);
        RealBall minus_s;
        arb_neg(minus_s.value(), s.value());
        if (!take(real_monic_quadratic(s.value(), T.value(), bits)) ||
            !take(real_monic_quadratic(minus_s.value(), V.value(), bits)))
            return std::nullopt;
    } else {
        // Complex coefficients: the resolvent root of largest modulus, z_k = ((xi^k Q + D0/(xi^k Q))/a - 2p)/3 for a cube
        // root Q of the larger of (D1 +- sqrt E)/2.
        ComplexBall root_e, plus, minus, cube, d0, xi, zk, z, s, ratio, T, V;
        ball(root_e.value(), E_, bits);
        square_root(root_e.value(), root_e.value(), bits);
        ball(plus.value(), D1_, bits);
        acb_sub(minus.value(), plus.value(), root_e.value(), bits);
        acb_add(plus.value(), plus.value(), root_e.value(), bits);
        Float plus_lower, minus_lower, best, lower;
        acb_get_abs_lbound_arf(plus_lower.value(), plus.value(), bits);
        acb_get_abs_lbound_arf(minus_lower.value(), minus.value(), bits);
        if (arf_cmp(minus_lower.value(), plus_lower.value()) > 0) acb_swap(plus.value(), minus.value());
        if (acb_contains_zero(plus.value())) return std::nullopt;
        acb_mul_2exp_si(plus.value(), plus.value(), -1);
        cube_root(cube.value(), plus.value(), bits);
        if (acb_contains_zero(cube.value())) return std::nullopt;
        ball(d0.value(), D0_, bits);
        arb_set_si(acb_realref(xi.value()), -1);
        arb_sqrt_ui(acb_imagref(xi.value()), 3, bits);
        acb_mul_2exp_si(xi.value(), xi.value(), -1);
        for (int k = 0; k < 3; ++k) {
            if (k) acb_mul(cube.value(), cube.value(), xi.value(), bits);
            acb_div(zk.value(), d0.value(), cube.value(), bits);
            acb_add(zk.value(), zk.value(), cube.value(), bits);
            acb_div(zk.value(), zk.value(), a.value(), bits);
            acb_submul_ui(zk.value(), p.value(), 2, bits);
            acb_div_ui(zk.value(), zk.value(), 3, bits);
            acb_get_abs_lbound_arf(lower.value(), zk.value(), bits);
            if (k == 0 || arf_cmp(lower.value(), best.value()) > 0) {
                acb_set(z.value(), zk.value());
                arf_set(best.value(), lower.value());
            }
        }
        square_root(s.value(), z.value(), bits);
        if (acb_contains_zero(s.value())) return std::nullopt;
        acb_div(ratio.value(), q.value(), s.value(), bits);
        acb_add(T.value(), p.value(), z.value(), bits);
        acb_add(V.value(), T.value(), ratio.value(), bits);
        acb_sub(T.value(), T.value(), ratio.value(), bits);
        acb_mul_2exp_si(T.value(), T.value(), -1);
        acb_mul_2exp_si(V.value(), V.value(), -1);
        ComplexBall minus_s;
        acb_neg(minus_s.value(), s.value());
        if (!take(complex_monic_quadratic(s.value(), T.value(), bits)) ||
            !take(complex_monic_quadratic(minus_s.value(), V.value(), bits)))
            return std::nullopt;
    }
    for (auto &y : ys) {
        ComplexBall x;
        if (real_) {
            // Exactly real roots stay exactly real: the real shift is added to the real part alone.
            arb_add(acb_realref(x.value()), acb_realref(y.value()), acb_realref(shift.value()), bits);
            arb_set(acb_imagref(x.value()), acb_imagref(y.value()));
        } else {
            acb_add(x.value(), y.value(), shift.value(), bits);
        }
        out.push_back(simple_root(std::move(x)));
    }
    return out;
}
} // namespace certpolysolve::detail
