#include "exact_reference_api.hpp"
#include "numeric/evaluate/midpoint_polynomial_evaluator.hpp"
#include "numeric/evaluate/blocked_polynomial_evaluator.hpp"
#include "numeric/candidates/candidates.hpp"
#include <iostream>

using namespace certpolysolve::detail;
namespace {
struct Point {
    acf_t value;
    Point() { acf_init(value); }
    ~Point() { acf_clear(value); }
};
reference::Rational exact(const arf_t value) {
    Rational rational;
    arf_get_fmpq(rational.value(), value);
    std::unique_ptr<char, decltype(&flint_free)> text(fmpq_get_str(nullptr, 10, rational.value()), flint_free);
    return reference::Rational(text.get());
}
void uncertain_balls() {
    for (slong length : {0, 1, 2, 8, 17, 63}) for (slong precision : {80, 192}) {
        ComplexPolynomial p;
        ComplexBall coefficient, point, value, derivative;
        for (slong j = 0; j < length; ++j) {
            acb_set_si_si(coefficient.value(), j%7-3, j%5-2);
            arb_add_error_2exp_si(acb_realref(coefficient.value()), -9);
            arb_add_error_2exp_si(acb_imagref(coefficient.value()), -9);
            acb_poly_set_coeff_acb(p.value(), j, coefficient.value());
        }
        BlockedPolynomialEvaluator evaluator(p, precision);
        acb_set_d_d(point.value(), .25, .5);
        arb_add_error_2exp_si(acb_realref(point.value()), -6);
        arb_add_error_2exp_si(acb_imagref(point.value()), -6);
        evaluator.evaluate(value.value(), derivative.value(), point.value());
        for (int re : {-1,1}) for (int im : {-1,1}) for (int side : {-1,0,1}) {
            reference::Complex z(reference::Rational(1)/4+reference::Rational(re)/64,
                                 reference::Rational(1)/2+reference::Rational(im)/64), y, d;
            for (slong j = length; j-- > 0;) {
                const reference::Rational delta = reference::Rational(side*(j%2 ? -1 : 1))/512;
                d = d*z + y;
                y = y*z + reference::Complex(reference::Rational(j%7-3)+delta,
                                             reference::Rational(j%5-2)-delta);
            }
            const auto contains = [](const acb_t ball, const reference::Complex &expected) {
                Rational real, imag;
                if (fmpq_set_str(real.value(), expected.real.str().c_str(), 10) ||
                    fmpq_set_str(imag.value(), expected.imag.str().c_str(), 10))
                    throw std::runtime_error("invalid independent test rational");
                reference::require(arb_contains_fmpq(acb_realref(ball),real.value()) &&
                                   arb_contains_fmpq(acb_imagref(ball),imag.value()),
                                   "blocked ball evaluation lost an exact uncertainty member");
            };
            contains(value.value(), y); contains(derivative.value(), d);
        }
    }
}
void check(slong length, slong precision, bool real, int scale) {
    ComplexPolynomial input;
    std::vector<reference::Complex> coefficients;
    ComplexBall coefficient;
    for (slong j = 0; j < length; ++j) {
        const slong re = (j * 17 + 3) % 23 - 11, im = real ? 0 : (j * 13 + 1) % 19 - 9;
        acb_set_si_si(coefficient.value(), re, im);
        acb_mul_2exp_si(coefficient.value(), coefficient.value(), scale);
        acb_poly_set_coeff_acb(input.value(), j, coefficient.value());
        coefficients.emplace_back(reference::Rational(re) * reference::power2(scale),
                                  reference::Rational(im) * reference::power2(scale));
    }
    MidpointPolynomialEvaluator evaluator(input.value(), precision);
    for (const auto &z : {reference::Complex(0), reference::Complex(1),
                          reference::Complex(reference::Rational(7)/8, reference::Rational(-3)/16)}) {
        Point point, value, derivative;
        const auto set = [](arf_t out, const reference::Rational &x) {
            Rational exact;
            const auto text = x.str();
            if (fmpq_set_str(exact.value(), text.c_str(), 10)) throw std::runtime_error("bad test rational");
            arf_set_fmpq(out, exact.value(), 1024, ARF_RND_NEAR);
        };
        set(acf_realref(point.value), z.real); set(acf_imagref(point.value), z.imag);
        evaluator.evaluate(value.value, derivative.value, point.value);
        reference::Complex y, d;
        reference::Rational envelope = 0, derivative_envelope = 0;
        const reference::Rational radius = reference::absolute(z.real) + reference::absolute(z.imag);
        for (auto i = coefficients.rbegin(); i != coefficients.rend(); ++i) {
            d = d * z + y; y = y * z + *i;
            derivative_envelope = derivative_envelope * radius + envelope;
            envelope = envelope * radius + reference::absolute(i->real) + reference::absolute(i->imag);
        }
        // Absolute scaled error, not relative error near cancellation. This
        // independent rational oracle catches value/derivative and block-edge
        // errors without treating candidate arithmetic as a proof enclosure.
        const auto verify = [&](const acf_t actual, const reference::Complex &expected,
                                const reference::Rational &majorant) {
            const reference::Rational error = reference::absolute(exact(acf_realref(actual)) - expected.real) +
                                              reference::absolute(exact(acf_imagref(actual)) - expected.imag);
            reference::require(error <= majorant * reference::power2(static_cast<int>(24 - precision)),
                               "midpoint block evaluation disagrees with exact rational arithmetic");
        };
        verify(value.value, y, envelope); verify(derivative.value, d, derivative_envelope);
    }
}
}
// The midpoint Aberth's own evaluation (candidates.cpp) has a sparse path,
// powers over the gaps between the nonzero coefficients, chosen when fewer
// than a quarter are nonzero. Until 30 September 2026 it never ran (its term
// list was built with a zero test that glibc's iszero template silenced), and
// it stopped at the lowest nonzero coefficient, returning p(z)/z^m: for a
// polynomial with a zero constant term the Newton correction of a polynomial
// without p's zero root. Checked against the rational oracle: sparse with
// m = 0, with m = 1 and m = 3, and a dense one through the same call.
void sparse_gaps(slong precision) {
    struct Shape { slong length; std::vector<slong> support; };
    for (const Shape &shape : {Shape{301, {300, 17, 0}}, Shape{301, {300, 1}}, Shape{257, {256, 40, 3}},
                               Shape{9, {8, 5, 3, 1, 0}}}) {
        ComplexPolynomial input;
        std::vector<reference::Complex> coefficients(static_cast<size_t>(shape.length));
        ComplexBall coefficient;
        for (const slong j : shape.support) {
            const slong re = (j * 17 + 3) % 23 - 11, im = (j * 13 + 1) % 19 - 9;
            acb_set_si_si(coefficient.value(), re, im);
            acb_poly_set_coeff_acb(input.value(), j, coefficient.value());
            coefficients[static_cast<size_t>(j)] = reference::Complex(reference::Rational(re), reference::Rational(im));
        }
        for (const auto &z : {reference::Complex(reference::Rational(7)/8, reference::Rational(-3)/16),
                              reference::Complex(reference::Rational(-33)/32, reference::Rational(1)/64)}) {
            ComplexBall point, value, derivative;
            const auto set = [](arb_t out, const reference::Rational &x) {
                Rational exact;
                if (fmpq_set_str(exact.value(), x.str().c_str(), 10)) throw std::runtime_error("bad test rational");
                arb_set_fmpq(out, exact.value(), 1024);
            };
            set(acb_realref(point.value()), z.real); set(acb_imagref(point.value()), z.imag);
            midpoint_polynomial_evaluation(value.value(), derivative.value(), input.value()->coeffs,
                                           acb_poly_length(input.value()), point.value(), precision);
            reference::Complex y, d;
            reference::Rational envelope = 0, derivative_envelope = 0;
            const reference::Rational radius = reference::absolute(z.real) + reference::absolute(z.imag);
            for (auto i = coefficients.rbegin(); i != coefficients.rend(); ++i) {
                d = d * z + y; y = y * z + *i;
                derivative_envelope = derivative_envelope * radius + envelope;
                envelope = envelope * radius + reference::absolute(i->real) + reference::absolute(i->imag);
            }
            const auto verify = [&](const acb_t actual, const reference::Complex &expected,
                                    const reference::Rational &majorant) {
                const reference::Rational error =
                    reference::absolute(exact(arb_midref(acb_realref(actual))) - expected.real) +
                    reference::absolute(exact(arb_midref(acb_imagref(actual))) - expected.imag);
                reference::require(error <= majorant * reference::power2(static_cast<int>(24 - precision)),
                                   "the sparse midpoint evaluation disagrees with exact rational arithmetic");
            };
            verify(value.value(), y, envelope); verify(derivative.value(), d, derivative_envelope);
        }
    }
}
int main() {
    certpolysolve::ThreadScope lifetime;
    try {
        uncertain_balls();
        for (slong precision : {64, 128, 340}) sparse_gaps(precision);
        for (slong length : {0, 1, 2, 3, 8, 16, 17, 63, 64, 65, 128, 257})
            for (slong precision : {64, 128, 340, 512})
                for (bool real : {false, true})
                    for (int scale : {-1000, 0, 1000}) check(length, precision, real, scale);
        std::cout << "midpoint values and derivatives agree with independent rational evaluations\n";
        return 0;
    } catch (const std::exception &e) { std::cerr << e.what() << '\n'; return 1; }
}
