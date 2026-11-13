// DiskSparsePolynomialEvaluator: value and slope of a polynomial with few
// terms. Checked against an independent evaluation (every term's power by
// FLINT's own powering at four times the precision, no chain and no table;
// the dense evaluation as well where the degree allows it), for containment
// and for width: the enclosures must be within a few units of rounding at the
// caller's precision, relative to the sum of the terms' sizes. Ball Horner
// over the terms, which this evaluator replaces in the certificate attempt,
// loses about half a bit a term and fails the width test from a few dozen
// terms on.
#include <certpolysolve/certpolysolve.hpp>
#include "numeric/evaluate/disk_polynomial_evaluator.hpp"
#include "numeric/evaluate/disk_sparse_polynomial_evaluator.hpp"
#include <iostream>
#include <set>
namespace d = certpolysolve::detail;
static void require(bool value, const char *why) {
    if (!value) throw std::runtime_error(why);
}
namespace {
struct Generator {
    uint64_t state;
    uint64_t next() {
        state = state * 6364136223846793005ULL + 1442695040888963407ULL;
        return state >> 33;
    }
};
// A polynomial of the given degree with `terms` nonzero coefficients, small
// Gaussian integers over a power of two; `with_constant` keeps exponent 0.
d::ComplexPolynomial sparse(Generator &random, slong degree, slong terms, bool with_constant, bool uncertain) {
    std::set<slong> support{degree};
    if (with_constant) support.insert(0);
    while (static_cast<slong>(support.size()) < terms)
        support.insert(static_cast<slong>(random.next() % static_cast<uint64_t>(degree)) + (with_constant ? 0 : 1));
    d::ComplexPolynomial p;
    d::ComplexBall c;
    for (const slong e : support) {
        const auto re = static_cast<slong>(random.next() % 19) - 9, im = static_cast<slong>(random.next() % 19) - 9;
        acb_set_si_si(c.value(), re ? re : 1, im);
        acb_mul_2exp_si(c.value(), c.value(), -static_cast<slong>(random.next() % 5));
        if (uncertain)
            arb_add_error_2exp_si(acb_realref(c.value()), -90);
        acb_poly_set_coeff_acb(p.value(), e, c.value());
    }
    return p;
}
// sum c_k z^k and sum k c_k z^(k-1) term by term at `precision`, and the sums
// of the terms' sizes.
void reference(acb_t value, acb_t slope, arb_t size, arb_t slope_size, const d::ComplexPolynomial &p,
               const acb_t z, slong precision) {
    d::ComplexBall power, term;
    d::RealBall magnitude;
    acb_zero(value);
    acb_zero(slope);
    arb_zero(size);
    arb_zero(slope_size);
    for (slong k = 0; k < acb_poly_length(p.value()); ++k) {
        const acb_srcptr c = p.value()->coeffs + k;
        if (acb_is_zero(c)) continue;
        acb_pow_ui(power.value(), z, static_cast<ulong>(k), precision);
        acb_mul(term.value(), c, power.value(), precision);
        acb_add(value, value, term.value(), precision);
        acb_abs(magnitude.value(), term.value(), precision);
        arb_add(size, size, magnitude.value(), precision);
        if (!k) continue;
        acb_pow_ui(power.value(), z, static_cast<ulong>(k - 1), precision);
        acb_mul(term.value(), c, power.value(), precision);
        acb_mul_si(term.value(), term.value(), k, precision);
        acb_add(slope, slope, term.value(), precision);
        acb_abs(magnitude.value(), term.value(), precision);
        arb_add(slope_size, slope_size, magnitude.value(), precision);
    }
}
// The enclosure's l1 radius is at most 2^(lost - precision) times `size`.
bool narrow(const acb_t enclosure, const arb_t size, slong precision, slong lost) {
    d::RealBall radius, allowed;
    arb_set_arf(allowed.value(), arb_midref(size));
    arb_mul_2exp_si(allowed.value(), allowed.value(), lost - precision);
    d::Magnitude total;
    mag_add(total.value(), arb_radref(acb_realref(enclosure)), arb_radref(acb_imagref(enclosure)));
    arf_set_mag(arb_midref(radius.value()), total.value());
    return arb_le(radius.value(), allowed.value());
}
} // namespace
int main() {
    certpolysolve::ThreadScope lifetime;
    try {
        Generator random{20261002};
        // The plan: its count is the least over the widths, and it is what an
        // evaluation performs.
        {
            const std::vector<slong> support{0, 3, 4, 19, 260, 1000, 70000};
            const auto plan = d::DiskSparsePolynomialEvaluator::plan(support);
            uint64_t least = UINT64_MAX;
            const unsigned bits = 17; // the largest gap, 69,000, has 17 bits
            for (unsigned width = 1; width <= bits; ++width) {
                const uint64_t windows = (bits + width - 1) / width;
                uint64_t count = windows * ((uint64_t(1) << width) - 2) + windows - 1;
                slong previous = 0;
                for (const slong e : support) {
                    for (auto gap = static_cast<uint64_t>(e - previous); gap; gap >>= width)
                        count += (gap & ((uint64_t(1) << width) - 1)) != 0;
                    previous = e;
                }
                least = std::min(least, count);
            }
            require(plan.multiplications == least, "the plan does not use the cheapest table width");
            require(d::DiskSparsePolynomialEvaluator::plan(std::vector<slong>{0}).multiplications == 0,
                    "a constant needs no multiplication");
            bool refused = false;
            try {
                (void)d::DiskSparsePolynomialEvaluator::plan(std::vector<slong>{4, 3});
            } catch (const std::invalid_argument &) {
                refused = true;
            }
            require(refused, "a descending support was accepted");
        }
        struct Case {
            slong degree, terms;
            bool with_constant, uncertain;
        };
        for (const auto &shape : {Case{40, 5, true, false}, Case{41, 12, false, false}, Case{1000, 30, true, false},
                                  Case{1200, 300, true, false}, Case{5000, 300, true, true},
                                  Case{5001, 2, false, false}, Case{70000, 259, true, false}}) {
            const auto p = sparse(random, shape.degree, shape.terms, shape.with_constant, shape.uncertain);
            for (const slong precision : {85, 170}) {
                d::DiskSparsePolynomialEvaluator evaluator(p, precision);
                require(evaluator.terms() == shape.terms, "the support was not read");
                for (const auto &z : {std::pair{0.6234375, 0.781494140625}, std::pair{0.25, -0.4375},
                                      std::pair{-1.0009765625, 0.03125}, std::pair{1.25, 0.75}}) {
                    d::ComplexBall point, value, slope, expected, expected_slope;
                    d::RealBall size, slope_size;
                    acb_set_d_d(point.value(), z.first, z.second);
                    const uint64_t before = evaluator.multiplications();
                    require(evaluator.evaluate(value.value(), slope.value(), point.value()), "an evaluation was refused");
                    require(evaluator.multiplications() - before == evaluator.plan().multiplications,
                            "an evaluation did not make the multiplications of its plan");
                    // Exact members of the coefficient family: the midpoints,
                    // and the ends of the uncertain real parts.
                    for (int side : {-1, 0, 1}) {
                        d::ComplexPolynomial member;
                        d::ComplexBall c;
                        for (slong k = 0; k < acb_poly_length(p.value()); ++k) {
                            const acb_srcptr source = p.value()->coeffs + k;
                            if (acb_is_zero(source)) continue;
                            acb_get_mid(c.value(), source);
                            if (side < 0)
                                arb_get_lbound_arf(arb_midref(acb_realref(c.value())), acb_realref(source), 1024);
                            if (side > 0)
                                arb_get_ubound_arf(arb_midref(acb_realref(c.value())), acb_realref(source), 1024);
                            acb_poly_set_coeff_acb(member.value(), k, c.value());
                        }
                        reference(expected.value(), expected_slope.value(), size.value(), slope_size.value(), member,
                                  point.value(), 4 * precision);
                        require(acb_contains(value.value(), expected.value()) &&
                                    acb_contains(slope.value(), expected_slope.value()),
                                "the sparse disk evaluation misses an independently evaluated instance");
                        if (shape.degree <= 1200) {
                            // The dense evaluation loses up to half a bit a
                            // term: give it the degree in bits to lose.
                            acb_poly_evaluate2(expected.value(), expected_slope.value(), member.value(), point.value(),
                                               4 * precision + shape.degree);
                            require(acb_contains(value.value(), expected.value()) &&
                                        acb_contains(slope.value(), expected_slope.value()),
                                    "the sparse disk evaluation misses FLINT's dense evaluation");
                        }
                    }
                    // Width, for exact coefficients: within 2^4 roundings of
                    // the terms' total size (the dot product's own rounding,
                    // the powers' with their guard bits, and the division).
                    if (!shape.uncertain)
                        require(narrow(value.value(), size.value(), precision, 4) &&
                                    narrow(slope.value(), slope_size.value(), precision, 4),
                                "the sparse disk evaluation is wider than a few roundings");
                    // An uncertain point: the enclosure holds the values at
                    // exact points of its ball.
                    d::ComplexBall wide = point, corner;
                    arb_add_error_2exp_si(acb_realref(wide.value()), -60);
                    arb_add_error_2exp_si(acb_imagref(wide.value()), -60);
                    require(evaluator.evaluate(value.value(), slope.value(), wide.value()),
                            "an uncertain point was refused");
                    for (int re : {-1, 1})
                        for (int im : {-1, 1}) {
                            acb_set(corner.value(), point.value());
                            d::RealBall shift;
                            arb_set_si(shift.value(), re);
                            arb_mul_2exp_si(shift.value(), shift.value(), -60);
                            arb_add(acb_realref(corner.value()), acb_realref(corner.value()), shift.value(), ARF_PREC_EXACT);
                            arb_set_si(shift.value(), im);
                            arb_mul_2exp_si(shift.value(), shift.value(), -60);
                            arb_add(acb_imagref(corner.value()), acb_imagref(corner.value()), shift.value(), ARF_PREC_EXACT);
                            d::ComplexPolynomial midpoints;
                            d::ComplexBall c;
                            for (slong k = 0; k < acb_poly_length(p.value()); ++k) {
                                if (acb_is_zero(p.value()->coeffs + k)) continue;
                                acb_get_mid(c.value(), p.value()->coeffs + k);
                                acb_poly_set_coeff_acb(midpoints.value(), k, c.value());
                            }
                            reference(expected.value(), expected_slope.value(), size.value(), slope_size.value(),
                                      midpoints, corner.value(), 4 * precision);
                            require(acb_contains(value.value(), expected.value()) &&
                                        acb_contains(slope.value(), expected_slope.value()),
                                    "the enclosure at an uncertain point misses a corner of the point's ball");
                        }
                }
                // Refusals: aliasing, a point that may be zero, a point that is not finite.
                d::ComplexBall point, value, slope;
                acb_set_d_d(point.value(), 0.5, 0.25);
                require(!evaluator.evaluate(value.value(), value.value(), point.value()) &&
                            !evaluator.evaluate(point.value(), slope.value(), point.value()),
                        "aliased outputs were accepted");
                acb_zero(point.value());
                require(!evaluator.evaluate(value.value(), slope.value(), point.value()), "the point zero was accepted");
                arb_add_error_2exp_si(acb_realref(point.value()), -80);
                require(!evaluator.evaluate(value.value(), slope.value(), point.value()),
                        "a point whose ball holds zero was accepted");
                acb_indeterminate(point.value());
                require(!evaluator.evaluate(value.value(), slope.value(), point.value()),
                        "an indeterminate point was accepted");
                // The remainder's majorant: the sum itself, to the magnitudes'
                // upward rounding, and never below it.
                acb_set_d_d(point.value(), 0.6234375, 0.781494140625);
                d::RealBall radius, bound, exact, term, rho;
                arb_one(radius.value());
                arb_mul_2exp_si(radius.value(), radius.value(), -57);
                evaluator.linear_remainder(bound.value(), point.value(), radius.value());
                acb_abs(rho.value(), point.value(), 4096);
                arb_add(rho.value(), rho.value(), radius.value(), 4096);
                for (slong k = 2; k < acb_poly_length(p.value()); ++k) {
                    if (acb_is_zero(p.value()->coeffs + k)) continue;
                    d::Magnitude upper;
                    acb_get_mag(upper.value(), p.value()->coeffs + k);
                    arf_set_mag(arb_midref(term.value()), upper.value());
                    mag_zero(arb_radref(term.value()));
                    arb_mul_si(term.value(), term.value(), k, 4096);
                    arb_mul_si(term.value(), term.value(), k - 1, 4096);
                    d::RealBall power;
                    arb_pow_ui(power.value(), rho.value(), static_cast<ulong>(k - 2), 4096);
                    arb_mul(term.value(), term.value(), power.value(), 4096);
                    arb_add(exact.value(), exact.value(), term.value(), 4096);
                }
                arb_mul(exact.value(), exact.value(), radius.value(), 4096);
                arb_mul(exact.value(), exact.value(), radius.value(), 4096);
                arb_mul_2exp_si(exact.value(), exact.value(), -1);
                require(arb_ge(bound.value(), exact.value()), "the sparse remainder bound is below the majorant");
                d::RealBall loose;
                arb_mul_2exp_si(loose.value(), exact.value(), -10);
                arb_add(loose.value(), loose.value(), exact.value(), 4096);
                require(arb_le(bound.value(), loose.value()), "the sparse remainder bound is far above the majorant");
                if (shape.degree <= 5001) {
                    d::DiskPolynomialEvaluator dense(d::PolynomialSnapshot(p.value()), precision);
                    d::RealBall dense_bound;
                    dense.linear_remainder(dense_bound.value(), point.value(), radius.value());
                    require(arb_ge(dense_bound.value(), exact.value()) && arb_le(dense_bound.value(), loose.value()),
                            "the dense remainder bound and the sparse one describe different sums");
                }
            }
        }
        std::cout << "sparse disk evaluation checks passed\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
