#pragma once
#include "internal.hpp"
#include <functional>
#include <map>

// Independent benchmark verifier: no native root solver, local model, Rouche
// certificate, or supplied-point polishing path is used here.
namespace giant_fft {
namespace d = certpolysolve::detail;
struct PointVerification {
    std::optional<d::Root> root;
    bool accurate = false;
    slong precision = 0;
    uint64_t evaluations = 0;
    std::string reason = "not attempted";
};

// F(z) = 2^(2b)|z-p|² - |z|² is convex for b > 0. Its maximum
// over a rectangle occurs at a corner, so four exact inequalities suffice.
inline bool
saved_point_accuracy(const acb_t point, const d::Root &root, uint64_t bits) {
    if (!acb_is_exact(point) || !acb_is_finite(point) || bits > static_cast<uint64_t>(WORD_MAX / 2)) return false;
    d::Rational px, py;
    arf_get_fmpq(px.value(), arb_midref(acb_realref(point)));
    arf_get_fmpq(py.value(), arb_midref(acb_imagref(point)));
    for (size_t i : { size_t(0), size_t(1) })
        for (size_t j : { size_t(2), size_t(3) }) {
            d::Rational x, y, dx, dy, left, right, term;
            uint64_t allowance = UINT64_MAX;
            const auto &ex = root.endpoints[i], &ey = root.endpoints[j];
            if (!d::Number(ex.mantissa_text() + "@" + ex.exponent_text()).rational_with_budget(x.value(), allowance) ||
                !d::Number(ey.mantissa_text() + "@" + ey.exponent_text()).rational_with_budget(y.value(), allowance))
                return false;
            fmpq_sub(dx.value(), x.value(), px.value());
            fmpq_sub(dy.value(), y.value(), py.value());
            fmpq_mul(left.value(), dx.value(), dx.value());
            fmpq_mul(term.value(), dy.value(), dy.value());
            fmpq_add(left.value(), left.value(), term.value());
            fmpq_mul_2exp(left.value(), left.value(), static_cast<ulong>(2 * bits));
            fmpq_mul(right.value(), x.value(), x.value());
            fmpq_mul(term.value(), y.value(), y.value());
            fmpq_add(right.value(), right.value(), term.value());
            if (fmpq_cmp(left.value(), right.value()) > 0) return false;
        }
    return true;
}

class OriginalCoefficientVerifier {
    certpolysolve::detail::DensePolynomial source_;
    struct Prepared {
        std::vector<d::ComplexBall> coefficients;
        std::vector<d::Magnitude> second;
    };
    std::map<slong, Prepared> prepared_;
    const Prepared &prepare(slong precision) {
        if (auto found = prepared_.find(precision); found != prepared_.end()) return found->second;
        Prepared p;
        p.coefficients.resize(source_.coefficients.size());
        p.second.resize(source_.coefficients.size() - 2);
        for (size_t i = 0; i < p.coefficients.size(); ++i) {
            source_.coefficients[i].enclose(p.coefficients[i].value(), precision);
            if (i >= 2) {
                acb_get_mag(p.second[i - 2].value(), p.coefficients[i].value());
                mag_mul_ui(p.second[i - 2].value(), p.second[i - 2].value(), i);
                mag_mul_ui(p.second[i - 2].value(), p.second[i - 2].value(), i - 1);
            }
        }
        return prepared_.emplace(precision, std::move(p)).first->second;
    }
    // Accumulate roundoff as a disk, rather than repeatedly rotating a
    // rectangle. Each center operation is enclosed by ACB, then recentered;
    // its radius is added to the propagated disk error. This avoids artificial
    // (|Re z|+|Im z|)^n wrapping at a fixed complex point.
    static bool evaluate(acb_t value,
                         acb_t derivative,
                         const std::vector<d::ComplexBall> &coefficients,
                         const acb_t z,
                         slong precision,
                         const std::function<bool()> &cancel) {
        if (!acb_is_exact(z)) return false;
        acb_zero(value);
        acb_zero(derivative);
        d::Magnitude norm, value_error, derivative_error, rounding;
        acb_get_mag(norm.value(), z);
        const auto recenter = [&](acb_t center, mag_t error) {
            mag_add(rounding.value(), arb_radref(acb_realref(center)), arb_radref(acb_imagref(center)));
            mag_add(error, error, rounding.value());
            acb_get_mid(center, center);
        };
        for (size_t i = coefficients.size(); i-- > 0;) {
            if (i % 1024 == 0 && cancel && cancel()) return false;
            acb_mul(derivative, derivative, z, precision);
            acb_add(derivative, derivative, value, precision);
            mag_mul(derivative_error.value(), derivative_error.value(), norm.value());
            mag_add(derivative_error.value(), derivative_error.value(), value_error.value());
            recenter(derivative, derivative_error.value());
            acb_mul(value, value, z, precision);
            acb_add(value, value, coefficients[i].value(), precision);
            mag_mul(value_error.value(), value_error.value(), norm.value());
            recenter(value, value_error.value());
        }
        acb_add_error_mag(value, value_error.value());
        acb_add_error_mag(derivative, derivative_error.value());
        return acb_is_finite(value) && acb_is_finite(derivative);
    }
    // On the containing disk, |p'(z)-p'(m)| <= R sup |p''|.
    // The coefficientwise nonnegative Horner bound has no complex wrapping.
    static bool jacobian_bound(acb_t jacobian,
                               const acb_t derivative,
                               const Prepared &p,
                               const acb_t box,
                               const acb_t midpoint,
                               slong precision,
                               const std::function<bool()> &cancel) {
        d::Magnitude radius, norm, bound;
        d::ComplexBall displacement;
        acb_sub(displacement.value(), box, midpoint, precision);
        acb_get_mag(radius.value(), displacement.value());
        acb_get_mag(norm.value(), box);
        for (size_t i = p.second.size(); i-- > 0;) {
            if (i % 1024 == 0 && cancel && cancel()) return false;
            mag_mul(bound.value(), bound.value(), norm.value());
            mag_add(bound.value(), bound.value(), p.second[i].value());
        }
        mag_mul(bound.value(), bound.value(), radius.value());
        acb_set(jacobian, derivative);
        acb_add_error_mag(jacobian, bound.value());
        return acb_is_finite(jacobian);
    }

  public:
    explicit OriginalCoefficientVerifier(const certpolysolve::detail::DensePolynomial &source)
      : source_(source) {
        while (!source_.coefficients.empty() && source_.coefficients.back().zero()) source_.coefficients.pop_back();
        if (source_.coefficients.size() < 2 || !std::all_of(source_.coefficients.begin(),
                                                            source_.coefficients.end(),
                                                            [](const auto &c) { return c.exact(); }))
            throw std::invalid_argument("verifier requires a nonconstant exact polynomial");
    }
    PointVerification verify(const acb_t saved_point,
                             uint64_t bits,
                             slong maximum_precision,
                             const std::function<bool()> &cancel = {}) {
        PointVerification result;
        if (!acb_is_exact(saved_point) || !acb_is_finite(saved_point) || bits < 2 ||
            bits > static_cast<uint64_t>(maximum_precision) || maximum_precision < 128) {
            result.reason = "invalid immutable point or verifier precision allowance";
            return result;
        }
        // The box is only a starting domain. Krawczyk contraction may move its
        // center; the saved point assessed by the final predicate never moves.
        for (slong precision = 128; precision <= maximum_precision;) {
            if (cancel && cancel()) {
                result.reason = "time limit";
                return result;
            }
            result.precision = precision;
            const auto &p = prepare(precision);
            d::ComplexBall box;
            acb_set(box.value(), saved_point);
            d::Float magnitude;
            acb_get_abs_ubound_arf(magnitude.value(), saved_point, precision);
            if (arf_is_zero(magnitude.value())) {
                // Strict relative accuracy at zero requires an exact zero root.
                if (source_.coefficients[0].zero() && !source_.coefficients[1].zero()) {
                    d::Root zero;
                    zero.multiplicity = 1;
                    zero.certificate = "exact simple zero from constant and linear coefficients";
                    zero.publish(precision);
                    result.root = std::move(zero);
                    result.accurate = true;
                    result.reason = "certified";
                } else
                    result.reason = "zero point lacks a simple exact zero proof";
                return result;
            }
            arf_mul_2exp_si(magnitude.value(), magnitude.value(), 1 - static_cast<slong>(bits));
            arb_add_error_arf(acb_realref(box.value()), magnitude.value());
            arb_add_error_arf(acb_imagref(box.value()), magnitude.value());
            for (unsigned contraction = 0; contraction < 8; ++contraction) {
                d::ComplexBall midpoint, value, derivative, jacobian, inverse, residual, displacement, image;
                acb_get_mid(midpoint.value(), box.value());
                result.evaluations += 2;
                if (!evaluate(value.value(), derivative.value(), p.coefficients, midpoint.value(), precision, cancel) ||
                    !jacobian_bound(
                      jacobian.value(), derivative.value(), p, box.value(), midpoint.value(), precision, cancel)) {
                    result.reason = cancel && cancel() ? "time limit" : "nonfinite original evaluation";
                    break;
                }
                if (acb_contains_zero(derivative.value())) {
                    result.reason = "midpoint derivative enclosure contains zero";
                    break;
                }
                acb_inv(inverse.value(), derivative.value(), precision);
                acb_get_mid(inverse.value(), inverse.value()); // Fixed real 2×2 preconditioner.
                if (!acb_is_finite(inverse.value()) || acb_is_zero(inverse.value())) {
                    result.reason = "invalid preconditioner";
                    break;
                }
                acb_mul(residual.value(), inverse.value(), jacobian.value(), precision);
                acb_neg(residual.value(), residual.value());
                acb_add_ui(residual.value(), residual.value(), 1, precision);
                acb_sub(displacement.value(), box.value(), midpoint.value(), precision);
                acb_mul(image.value(), residual.value(), displacement.value(), precision);
                acb_submul(image.value(), inverse.value(), value.value(), precision);
                acb_add(image.value(), image.value(), midpoint.value(), precision);
                d::RealBall norm, one;
                acb_abs(norm.value(), residual.value(), precision);
                arb_one(one.value());
                const bool unique = arb_lt(norm.value(), one.value());
                const bool interior = arb_contains_interior(acb_realref(box.value()), acb_realref(image.value())) &&
                                      arb_contains_interior(acb_imagref(box.value()), acb_imagref(image.value()));
                if (!unique || !interior) {
                    result.reason = !unique ? "Jacobian contraction not proved" : "Krawczyk inclusion not proved";
                    break;
                }
                d::Root root;
                root.box = image;
                root.multiplicity = 1;
                root.certificate = "original-coefficient 2D Krawczyk strict inclusion and contraction";
                root.publish(precision);
                result.accurate = saved_point_accuracy(saved_point, root, bits);
                result.root = std::move(root);
                if (result.accurate) {
                    result.reason = "certified";
                    return result;
                }
                result.reason = "unique root proved; saved-point accuracy remains inconclusive";
                box = std::move(image);
            }
            if (precision == maximum_precision) break;
            precision = precision > maximum_precision / 2 ? maximum_precision : precision * 2;
        }
        return result;
    }
};
} // namespace giant_fft
