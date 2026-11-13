// A reused evaluator cache must agree bit for bit with a fresh preparation and
// enclose an independent exact Horner value. Alternate point/origin, precision,
// sparse/dense coefficients and value/joint mode, including unsupported queries.
#include "numeric/evaluate/structured_source.hpp"
#include "numeric/evaluate/secular_disk_values.hpp"
#include <certpolysolve/certpolysolve.hpp>
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace d = certpolysolve::detail;
namespace {
void require(bool value, const std::string &why) {
    if (!value) throw std::runtime_error(why);
}

struct Outcome {
    bool available = false, limb = false, success = false;
    d::ComplexBall value, derivative;
};

template <class Owner, class Node>
concept CanPrepare = requires(Owner &&owner, Node &&node) {
    std::forward<Owner>(owner).prepare(64, std::forward<Node>(node), true);
};
static_assert(CanPrepare<d::DiskValues &, const d::ComplexBall &>);
static_assert(!CanPrepare<d::DiskValues, const d::ComplexBall &>);
static_assert(!CanPrepare<d::DiskValues &, d::ComplexBall>);
static_assert(!CanPrepare<d::DiskValues &, const d::ComplexBall>);
static_assert(!std::is_move_constructible_v<d::DiskValues> && !std::is_copy_constructible_v<d::DiskValues>);
static_assert(!std::is_copy_constructible_v<d::DiskValues::Prepared>);
static_assert(std::is_trivially_destructible_v<d::DiskValues::Prepared>);

// The arithmetic/cache checks below also passed with the original
// available()/value()/joint() protocol before replacing this adapter.
Outcome evaluate(d::DiskValues &disk, const d::ComplexBall &point, slong precision, bool value_only) {
    Outcome out;
    acb_set_si(out.derivative.value(), 17);
    if (auto prepared = disk.prepare(precision, point, value_only)) {
        out.available = true;
        out.limb = prepared->uses_limb_horner();
        out.success = prepared->evaluate(out.value.value(), out.derivative.value());
    }
    return out;
}

d::ComplexPolynomial polynomial(slong degree, bool dense) {
    d::ComplexPolynomial p;
    d::ComplexBall coefficient;
    acb_set_si(coefficient.value(), 2);
    for (slong k = 0; k <= degree; ++k)
        if (dense || k == 0 || k == 1 || k == degree)
            acb_poly_set_coeff_acb(p.value(), k, coefficient.value());
    return p;
}

std::vector<double> coefficient_sizes(const d::ComplexPolynomial &p) {
    std::vector<double> sizes;
    d::Magnitude size;
    for (slong k = 0; k < acb_poly_length(p.value()); ++k) {
        acb_get_mag(size.value(), p.value()->coeffs + k);
        sizes.push_back(mag_is_zero(size.value()) ? -HUGE_VAL : mag_get_d_log2_approx(size.value()));
    }
    return sizes;
}

d::PolynomialEncloser encloser(const d::ComplexPolynomial &p) {
    return [snapshot = d::PolynomialSnapshot(p.value())](acb_poly_t out, slong precision) {
        acb_poly_set_round(out, snapshot.value(), precision);
    };
}

// Exact Taylor translation by multiplying the current prefix by x + offset.
// Scalar dyadic operations keep this test independent of a rounded polynomial
// translation routine and give the source exact coefficients to enclose.
d::ComplexPolynomial translated(const d::ComplexPolynomial &p, const acb_t offset) {
    const slong length = acb_poly_length(p.value());
    std::vector<d::ComplexBall> coefficients(static_cast<size_t>(length));
    for (slong k = length; k-- > 0;) {
        for (slong j = length; --j > 0;) {
            acb_mul(coefficients[j].value(), coefficients[j].value(), offset, ARF_PREC_EXACT);
            acb_add(coefficients[j].value(), coefficients[j].value(), coefficients[j - 1].value(), ARF_PREC_EXACT);
        }
        acb_mul(coefficients[0].value(), coefficients[0].value(), offset, ARF_PREC_EXACT);
        acb_add(coefficients[0].value(), coefficients[0].value(), p.value()->coeffs + k, ARF_PREC_EXACT);
    }
    d::ComplexPolynomial result;
    for (slong j = 0; j < length; ++j) {
        require(acb_is_exact(coefficients[j].value()), "test translation was rounded");
        acb_poly_set_coeff_acb(result.value(), j, coefficients[j].value());
    }
    return result;
}

// Every coefficient is an integer and every test point a dyadic complex
// number. Exact scalar Horner arithmetic therefore gives an exact reference;
// the original leading coefficient is 2, so monic normalization is exact too.
void reference(acb_t value, acb_t derivative, const d::ComplexPolynomial &p,
               const acb_t point, const acb_t offset) {
    d::ComplexBall original_point;
    acb_add(original_point.value(), point, offset, ARF_PREC_EXACT);
    acb_zero(value);
    acb_zero(derivative);
    for (slong k = acb_poly_length(p.value()); k-- > 0;) {
        acb_mul(derivative, derivative, original_point.value(), ARF_PREC_EXACT);
        acb_add(derivative, derivative, value, ARF_PREC_EXACT);
        acb_mul(value, value, original_point.value(), ARF_PREC_EXACT);
        acb_add(value, value, p.value()->coeffs + k, ARF_PREC_EXACT);
    }
    acb_mul_2exp_si(value, value, -1);
    acb_mul_2exp_si(derivative, derivative, -1);
    require(acb_is_exact(value) && acb_is_exact(derivative), "reference arithmetic was rounded");
}

void selections(bool dense, bool centered) {
    constexpr slong degree = 128;
    const auto original = polynomial(degree, dense);
    d::ComplexBall offset;
    if (centered) acb_set_si(offset.value(), 8);
    const auto shifted = translated(original, offset.value());
    d::PolynomialSourceOptions options;
    options.input = d::InputSemantics::Exact;
    if (centered) {
        options.origin = encloser(original);
        acb_set(options.origin_offset.value(), offset.value());
        options.centred_log2 = coefficient_sizes(shifted);
        options.origin_log2 = coefficient_sizes(original);
    }
    d::PolynomialSource source(encloser(shifted), degree, true, options);
    require(source.rectangular_available(32) == dense, "test source has the wrong density classification");

    std::array<d::ComplexBall, 5> points;
    acb_set_d_d(points[1].value(), .5, .25);
    acb_sub(points[2].value(), points[1].value(), offset.value(), ARF_PREC_EXACT);
    acb_one(points[3].value());
    d::ComplexBall tiny;
    acb_one(tiny.value());
    // More than four machine limbs separates the nonzero mantissa bits: this
    // point is exact, but cannot use the short-multiplier Horner evaluator.
    acb_mul_2exp_si(tiny.value(), tiny.value(), -5 * FLINT_BITS);
    acb_add(points[3].value(), points[3].value(), tiny.value(), ARF_PREC_EXACT);
    acb_set(points[4].value(), points[1].value());
    arb_add_error_2exp_si(acb_realref(points[4].value()), -80);
    if (centered) {
        d::ComplexBall origin;
        require(!source.origin_point(origin.value(), points[1].value()) &&
                    source.origin_point(origin.value(), points[2].value()),
                "test does not exercise both origin choices");
    }

    // More than eight precisions evicts entries from the source's coefficient
    // cache. DiskValues' owned evaluators must still work when revisited.
    for (bool limbs : {false, true}) {
        d::DiskValues reused(source, limbs);
        for (slong precision : {32, 64, 96, 128, 160, 192, 224, 256, 288, 320, 32, 96}) {
            for (bool value_only : {true, false}) {
                for (size_t i = 0; i < points.size(); ++i) {
                    const auto point = points[i].value();
                    const std::string label = std::string(dense ? "dense" : "sparse") +
                        (centered ? " centered" : "") + ", precision " + std::to_string(precision) +
                        ", point " + std::to_string(i) + (value_only ? ", value" : ", joint") +
                        (limbs ? ", limbs" : ", blocks");
                    const bool short_point = i < 3;
                    const bool limb = limbs && dense && value_only && short_point;
                    const bool available = precision <= degree || limb;
                    auto actual = evaluate(reused, points[i], precision, value_only);
                    d::DiskValues fresh(source, limbs);
                    auto alone = evaluate(fresh, points[i], precision, value_only);
                    require(actual.available == available, label + ": preparation availability changed");
                    require(actual.limb == (available && limb), label + ": wrong evaluator chosen");
                    require(actual.available == alone.available && actual.limb == alone.limb &&
                                actual.success == alone.success,
                            label + ": reused cache changed preparation or evaluation status");
                    // The sparse joint evaluator recovers p' from z p'(z),
                    // so it refuses zero even though the schedule is available.
                    const bool sparse_joint_at_zero = !dense && !value_only && i == 0;
                    require(actual.success == (available && !sparse_joint_at_zero),
                            label + ": evaluation success/refusal changed");
                    if (!actual.success) continue;
                    require(acb_equal(actual.value.value(), alone.value.value()) &&
                                acb_equal(actual.derivative.value(), alone.derivative.value()),
                            label + ": reused selection differs from fresh evaluation");
                    if (value_only)
                        require(acb_equal_si(actual.derivative.value(), 17), label + ": value-only overwrote derivative");
                    // Blocked and sparse schedules also enclose uncertain
                    // points. Check exact members of the returned enclosure.
                    const int spread = i == 4 ? 1 : 0;
                    for (int side = -spread; side <= spread; ++side) {
                        d::ComplexBall sample, delta, value, derivative;
                        acb_get_mid(sample.value(), point);
                        if (side) {
                            acb_set_si(delta.value(), side);
                            acb_mul_2exp_si(delta.value(), delta.value(), -81);
                            acb_add(sample.value(), sample.value(), delta.value(), ARF_PREC_EXACT);
                        }
                        reference(value.value(), derivative.value(), original, sample.value(), offset.value());
                        require(acb_contains(actual.value.value(), value.value()),
                                label + ": monic value misses an exact member");
                        if (!value_only)
                            require(acb_contains(actual.derivative.value(), derivative.value()),
                                    label + ": monic derivative misses an exact member");
                    }
                }
            }
        }
    }
}
} // namespace

int main() {
    certpolysolve::ThreadScope lifetime;
    try {
        selections(true, false);
        selections(false, false);
        selections(true, true);
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
