#include <certpolysolve/certpolysolve.hpp>
#include "numeric/evaluate/disk_polynomial_evaluator.hpp"
#include "numeric/evaluate/disk_blocked_polynomial_evaluator.hpp"
#include "numeric/evaluate/polynomial_remainder.hpp"
#include <iostream>
namespace d = certpolysolve::detail;
template<class T>
concept BlockedStorage = requires { typename d::OwnedDiskBlockedPolynomialEvaluator<T>; };
static_assert(BlockedStorage<d::ComplexPolynomial> && BlockedStorage<d::PolynomialSnapshot>);
static_assert(!BlockedStorage<acb_poly_struct *> && !BlockedStorage<d::ConstComplexBallView>);
static_assert(!std::is_copy_constructible_v<d::SharedDiskBlockedPolynomialEvaluator>);
static_assert(std::is_nothrow_move_constructible_v<d::SharedDiskBlockedPolynomialEvaluator>);

static void
require(bool value, const char *why) {
    if (!value) throw std::runtime_error(why);
}
// Uses high-precision rectangular Horner evaluation of exact members of each
// coefficient/point family, independently of the blocked disk-error recurrence.
static void exact_member(d::ComplexBall &out, const acb_t family, int real_side, int imag_side) {
    acb_get_mid(out.value(), family);
    const auto end = [](arb_t value, const arb_t interval, int side) {
        if (side < 0) arb_get_lbound_arf(arb_midref(value), interval, 1024);
        if (side > 0) arb_get_ubound_arf(arb_midref(value), interval, 1024);
    };
    end(acb_realref(out.value()), acb_realref(family), real_side);
    end(acb_imagref(out.value()), acb_imagref(family), imag_side);
}

template<class Blocked>
static void check_blocked(Blocked &evaluate, const d::ComplexPolynomial &family, slong length) {
    d::ComplexBall point, value, derivative, only, exact_point, reference, reference_derivative;
    acb_set_d_d(point.value(), .5, .25);
    arb_add_error_2exp_si(acb_realref(point.value()), -12);
    arb_add_error_2exp_si(acb_imagref(point.value()), -13);
    require(evaluate.evaluate(value.value(), derivative.value(), point.value()),
            "blocked evaluation refused an uncertain point");
    require(evaluate.evaluate_value(only.value(), point.value()) && acb_equal(value.value(), only.value()),
            "blocked value-only path changed its enclosure");
    for (int side : {-1, 0, 1}) {
        d::ComplexPolynomial exact;
        for (slong j = 0; j < length; ++j) {
            d::ComplexBall coefficient;
            exact_member(coefficient, family.value()->coeffs + j, (j % 2 ? -side : side), side);
            acb_poly_set_coeff_acb(exact.value(), j, coefficient.value());
        }
        for (int re : {-1, 0, 1}) for (int im : {-1, 0, 1}) {
            exact_member(exact_point, point.value(), re, im);
            acb_poly_evaluate2(reference.value(), reference_derivative.value(), exact.value(),
                               exact_point.value(), 1024);
            require(acb_contains(value.value(), reference.value()) &&
                    acb_contains(derivative.value(), reference_derivative.value()),
                    "blocked enclosure misses an exact coefficient/point instance");
        }
    }
    d::ComplexBall serial_value(value), serial_derivative(derivative);
    {
        d::TeamScope workers(4);
        require(evaluate.evaluate(value.value(), derivative.value(), point.value()),
                "parallel blocked evaluation refused");
        require(acb_equal(value.value(), serial_value.value()) &&
                acb_equal(derivative.value(), serial_derivative.value()),
                "blocked serial/parallel enclosures differ");
        if (length >= 4096)
            require(d::split_jobs()[static_cast<size_t>(d::SplitSite::BlockedEvaluation)] != 0,
                    "blocked comparison missed the parallel kernel");
    }
    const d::ComplexBall saved_point(point), saved_value(value), saved_derivative(derivative);
    require(!evaluate.evaluate(value.value(), value.value(), point.value()) &&
            !evaluate.evaluate(point.value(), derivative.value(), point.value()) &&
            !evaluate.evaluate(value.value(), point.value(), point.value()) &&
            !evaluate.evaluate_value(point.value(), point.value()), "blocked output alias was accepted");
    require(acb_equal(point.value(), saved_point.value()) && acb_equal(value.value(), saved_value.value()) &&
            acb_equal(derivative.value(), saved_derivative.value()), "alias refusal changed an output");
}

static void blocked_ownership_and_enclosures() {
    for (slong length : {0, 1, 2, 17, 81, 4097}) {
        d::ComplexPolynomial source;
        for (slong j = 0; j < length; ++j) {
            d::ComplexBall coefficient;
            acb_set_si(coefficient.value(), j % 7 - 3);
            arb_set_si(acb_imagref(coefficient.value()), j % 5 - 2);
            arb_add_error_2exp_si(acb_realref(coefficient.value()), -100);
            arb_add_error_2exp_si(acb_imagref(coefficient.value()), -90);
            acb_poly_set_coeff_acb(source.value(), j, coefficient.value());
        }
        for (slong bits : {64, 128, 257}) {
            d::ComplexPolynomial owned;
            acb_poly_set(owned.value(), source.value());
            d::DiskBlockedPolynomialEvaluator original(std::move(owned), bits);
            d::DiskBlockedPolynomialEvaluator retained(std::move(original));
            // Moved-from evaluation refuses before borrowing its transferred storage.
            d::ComplexBall point, value, derivative;
            require(!original.evaluate(value.value(), derivative.value(), point.value()),
                    "moved blocked evaluator accepted a call");
            check_blocked(retained, source, length);
            // The snapshot and polynomial used to construct this evaluator
            // are gone before its first call. Only its immutable owner remains.
            auto shared_original = [&] {
                d::ComplexPolynomial local;
                acb_poly_set(local.value(), source.value());
                d::PolynomialSnapshot snapshot(std::move(local));
                return d::SharedDiskBlockedPolynomialEvaluator(snapshot, bits);
            }();
            d::SharedDiskBlockedPolynomialEvaluator shared(std::move(shared_original));
            require(!shared_original.evaluate(value.value(), derivative.value(), point.value()) &&
                    !shared_original.evaluate_value(value.value(), point.value()),
                    "moved shared evaluator borrowed transferred storage");
            check_blocked(shared, source, length);
            for (const auto &z : {std::pair{0., 0.}, std::pair{.5, .25}, std::pair{-1., 2.}}) {
                d::ComplexBall shared_value, shared_derivative;
                acb_set_d_d(point.value(), z.first, z.second);
                arb_add_error_2exp_si(acb_realref(point.value()), -12);
                arb_add_error_2exp_si(acb_imagref(point.value()), -13);
                require(retained.evaluate(value.value(), derivative.value(), point.value()) &&
                        shared.evaluate(shared_value.value(), shared_derivative.value(), point.value()),
                        "owning/shared comparison refused a finite point");
                require(acb_equal(value.value(), shared_value.value()) &&
                        acb_equal(derivative.value(), shared_derivative.value()),
                        "sharing immutable coefficients changed blocked evaluation bits");
            }
        }
    }
}

int
main() {
    certpolysolve::ThreadScope lifetime;
    try {
        blocked_ownership_and_enclosures();
        d::ComplexPolynomial p;
        d::ComplexBall coefficient, point, value, derivative, reference, reference_derivative;
        for (slong j = 0; j < 81; ++j) {
            acb_set_si(coefficient.value(), j % 7 - 3);
            arb_set_si(acb_imagref(coefficient.value()), j % 5 - 2);
            arb_add_error_2exp_si(acb_realref(coefficient.value()), -100);
            acb_poly_set_coeff_acb(p.value(), j, coefficient.value());
        }
        d::DiskPolynomialEvaluator evaluate(d::PolynomialSnapshot(p.value()), 128);
        for (const auto &z : { std::pair{ .5, .25 }, std::pair{ 1., 1. }, std::pair{ 0., 0. } }) {
            acb_set_d_d(point.value(), z.first, z.second);
            require(evaluate.evaluate(value.value(), derivative.value(), point.value()), "point evaluation refused");
            // A rectangular evaluation of uncertain coefficients may itself
            // wrap more widely than our disk enclosure. Compare exact members
            // of that coefficient family, not containment of the wider method.
            for (int side : {-1, 0, 1}) {
                d::ComplexPolynomial exact;
                for (slong j = 0; j < acb_poly_length(p.value()); ++j) {
                    acb_get_mid(coefficient.value(),p.value()->coeffs+j);
                    if (side < 0)
                        arb_get_lbound_arf(arb_midref(acb_realref(coefficient.value())),acb_realref(p.value()->coeffs+j),1024);
                    if (side > 0)
                        arb_get_ubound_arf(arb_midref(acb_realref(coefficient.value())),acb_realref(p.value()->coeffs+j),1024);
                    acb_poly_set_coeff_acb(exact.value(),j,coefficient.value());
                }
                acb_poly_evaluate2(reference.value(),reference_derivative.value(),exact.value(),point.value(),1024);
                require(acb_contains(value.value(),reference.value()) && acb_contains(derivative.value(),reference_derivative.value()),
                    "disk evaluation misses an independently evaluated exact coefficient instance");
            }
        }
        acb_set_d_d(point.value(), .5, .25);
        d::RealBall radius, remainder, separate_remainder, observed;
        arb_one(radius.value());
        arb_mul_2exp_si(radius.value(), radius.value(), -8);
        evaluate.linear_remainder(remainder.value(), point.value(), radius.value());
        // A remainder-only value must outlive the coefficients it borrowed at
        // construction. Its reference below is the independently shifted tail.
        const auto remainder_only = [&p] {
            d::ComplexPolynomial temporary;
            acb_poly_set(temporary.value(), p.value());
            return d::PolynomialLinearRemainder(temporary.coefficients());
        }();
        remainder_only.linear_remainder(separate_remainder.value(), point.value(), radius.value());
        d::ComplexPolynomial shifted;
        acb_poly_taylor_shift(shifted.value(), p.value(), point.value(), 1024);
        for (int re : { -1, 0, 1 })
            for (int im : { -1, 0, 1 }) {
                d::ComplexBall delta, term, residual;
                acb_set_d_d(delta.value(), re / 512., im / 512.);
                for (slong j = acb_poly_length(shifted.value()); j-- > 2;) {
                    acb_mul(residual.value(), residual.value(), delta.value(), 1024);
                    acb_add(residual.value(), residual.value(), shifted.value()->coeffs + j, 1024);
                }
                acb_mul(residual.value(), residual.value(), delta.value(), 1024);
                acb_mul(residual.value(), residual.value(), delta.value(), 1024);
                acb_abs(observed.value(), residual.value(), 1024);
                require(arb_le(observed.value(), remainder.value()),
                        "linear remainder misses exact shifted-polynomial tail");
                require(arb_le(observed.value(), separate_remainder.value()),
                        "owned remainder misses shifted tail after source destruction");
            }
        for (slong degree = -1; degree <= 1; ++degree) {
            d::ComplexPolynomial affine;
            if (degree >= 0)
                acb_poly_set_coeff_si(affine.value(), degree, 3);
            const d::PolynomialLinearRemainder affine_remainder(affine.coefficients());
            affine_remainder.linear_remainder(separate_remainder.value(), point.value(), radius.value());
            require(arb_is_zero(separate_remainder.value()), "affine polynomial has a nonlinear remainder");
        }
        arb_add_error_2exp_si(acb_realref(point.value()), -80);
        require(!evaluate.evaluate(value.value(), derivative.value(), point.value()),
                "disk point API accepted an uncertain input point");
        // Snapshot owns its coefficients across mutation of the source.
        acb_poly_zero(p.value());
        acb_zero(point.value());
        require(evaluate.evaluate(value.value(), derivative.value(), point.value()) && !acb_is_zero(value.value()),
                "evaluator retained a mutable coefficient borrow");
        std::cout << "disk evaluation and Taylor remainder enclosure checks passed\n";
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
