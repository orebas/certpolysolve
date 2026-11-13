// The certificate attempt's first tier (compensated_certificate.cpp), checked against 1,024-bit ball arithmetic:
// - the evaluation's bounds hold: |p(x) - (s + c)| <= value_error and |p'(x) - d| <= slope_error, at points near
//   roots, at roots that are binary64 numbers, and at arbitrary points, for real and complex, dense, sparse and wide
//   coefficients; the p'' majorant holds at points of its disk;
// - every rectangle the tier proves holds exactly one root: each root of a polynomial with planted dyadic roots is
//   tested for membership exactly;
// - the Newton retry: points too far for the first proof but in Newton's quadratic range are proved from the point
//   nearest the step's end, each rectangle holding exactly one planted root;
// - the smaller circles: at a 10-bit goal, where the goal's circles fail on the second-derivative term, circles of a
//   smaller radius from the same evaluation prove, each rectangle holding exactly one planted root and meeting the goal;
// - it declines where it must: a double root, roots closer than the goal's width, a goal beyond its precision;
// - the double-double repair (refine, proposal 1): from starts as poor as binary64 leaves the roots of a small cluster
//   (a real triple cluster 2^-18 wide; a real pair 2^-40 apart started as a conjugate pair, the symmetry trap), the
//   repaired points are each proved by the tier, each rectangle holding exactly one planted root and no two the same,
//   and the points not moved stay exactly where they were;
// - through the public API, Auto's solves use it (the attempt's count of its certificates) and every box holds
//   exactly one root of a 400-bit solve's boxes, which the ball certificate proves (the tier declines at that goal).
#include "certify/attempts/compensated_certificate.hpp"
#include "certify/chart_certificate.hpp"
#include "core/cover.hpp"
#include "core/raii.hpp"
#include <certpolysolve/certpolysolve.hpp>
#include <flint/acb_poly.h>
#include <flint/fmpq_poly.h>
#include "core/heuristics.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace d = certpolysolve::detail;

namespace {
constexpr slong reference_bits = 1024;

void require(bool ok, const std::string &why) {
    if (!ok) throw std::runtime_error(why);
}

// p = prod (x - root_j) with dyadic roots m 2^-shift (real and imaginary parts), as an acb_poly at the reference
// precision: exact where the coefficients fit it, and otherwise balls that hold the exact coefficients. The tier's
// proofs hold for every polynomial in the balls, so the exact one's planted roots can be checked against them.
struct Planted {
    d::ComplexPolynomial polynomial;
    std::vector<d::ComplexBall> roots;
};
Planted planted(const std::vector<std::pair<int64_t, int64_t>> &mantissas, slong shift) {
    Planted out;
    acb_poly_one(out.polynomial.value());
    d::ComplexPolynomial factor;
    for (const auto &[re, im] : mantissas) {
        d::ComplexBall root;
        acb_set_si_si(root.value(), re, im);
        acb_mul_2exp_si(root.value(), root.value(), -shift);
        acb_poly_zero(factor.value());
        acb_poly_set_coeff_si(factor.value(), 1, 1);
        d::ComplexBall negated;
        acb_neg(negated.value(), root.value());
        acb_poly_set_coeff_acb(factor.value(), 0, negated.value());
        acb_poly_mul(out.polynomial.value(), out.polynomial.value(), factor.value(), reference_bits);
        out.roots.push_back(std::move(root));
    }
    return out;
}

// The tier's bounds at x against the reference values (in the tier's units, the polynomial times 2^-scale).
void check_bounds(const d::CompensatedCertificate &tier, const acb_poly_t p, double xr, double xi,
                  const std::string &name) {
    const auto e = tier.evaluate(xr, xi);
    if (!e.ok) return; // declined: nothing claimed
    d::ComplexBall x, value, slope, computed, difference;
    acb_set_d_d(x.value(), xr, xi);
    acb_poly_evaluate2(value.value(), slope.value(), p, x.value(), reference_bits);
    acb_mul_2exp_si(value.value(), value.value(), -tier.scale());
    acb_mul_2exp_si(slope.value(), slope.value(), -tier.scale());
    // s + c, exactly.
    d::ComplexBall s, c;
    acb_set_d_d(s.value(), e.value_re, e.value_im);
    acb_set_d_d(c.value(), e.correction_re, e.correction_im);
    acb_add(computed.value(), s.value(), c.value(), reference_bits);
    acb_sub(difference.value(), value.value(), computed.value(), reference_bits);
    d::RealBall size, bound;
    acb_abs(size.value(), difference.value(), reference_bits);
    arb_set_d(bound.value(), e.value_error);
    require(arb_le(size.value(), bound.value()), name + ": the value's bound fails");
    acb_set_d_d(computed.value(), e.slope_re, e.slope_im);
    acb_sub(difference.value(), slope.value(), computed.value(), reference_bits);
    acb_abs(size.value(), difference.value(), reference_bits);
    arb_set_d(bound.value(), e.slope_error);
    require(arb_le(size.value(), bound.value()), name + ": the slope's bound fails");
}

// |p''| at points of the disk about x of radius `reach` against the majorant.
void check_second(const d::CompensatedCertificate &tier, const acb_poly_t p, double xr, double xi, double reach,
                  std::mt19937_64 &rng, const std::string &name) {
    const auto e = tier.evaluate(xr, xi);
    if (!e.ok) return;
    const double bound = tier.second_derivative_bound(e.modulus, reach);
    if (!std::isfinite(bound)) return;
    d::ComplexPolynomial second;
    acb_poly_derivative(second.value(), p, reference_bits);
    acb_poly_derivative(second.value(), second.value(), reference_bits);
    std::uniform_real_distribution<double> unit(-1, 1);
    for (int k = 0; k < 8; ++k) {
        d::ComplexBall z, value;
        acb_set_d_d(z.value(), xr + reach * unit(rng) * 0.7, xi + reach * unit(rng) * 0.7);
        acb_poly_evaluate(value.value(), second.value(), z.value(), reference_bits);
        acb_mul_2exp_si(value.value(), value.value(), -tier.scale());
        d::RealBall size, limit;
        acb_abs(size.value(), value.value(), reference_bits);
        arb_set_d(limit.value(), bound);
        require(arb_le(size.value(), limit.value()), name + ": the second derivative's majorant fails");
    }
}

cps_options goal(uint64_t bits) {
    cps_options o;
    o.bits = bits;
    o.goal = CPS_ROOT_APPROXIMATION;
    return o;
}

// Every root the tier proves: its published rectangle holds exactly one planted root.
uint64_t check_planted(const Planted &input, double perturbation, uint64_t bits, const std::string &name,
                       bool expect_all) {
    const d::CompensatedCertificate tier(input.polynomial.value());
    require(tier.available(), name + ": the tier is unavailable");
    const auto o = goal(bits);
    uint64_t proved = 0;
    for (const auto &root : input.roots) {
        d::ComplexBall point;
        const double xr = arf_get_d(arb_midref(acb_realref(root.value())), ARF_RND_NEAR) * (1 + perturbation);
        const double xi = arf_get_d(arb_midref(acb_imagref(root.value())), ARF_RND_NEAR) * (1 - perturbation);
        acb_set_d_d(point.value(), xr, xi);
        d::Root out;
        if (!tier.certify(out, point.value(), o, static_cast<slong>(bits + 32))) {
            require(!expect_all, name + ": a well-separated root was not proved");
            continue;
        }
        ++proved;
        int held = 0;
        for (const auto &other : input.roots) held += acb_contains(out.box.value(), other.value());
        require(held == 1, name + ": a proved rectangle holds " + std::to_string(held) + " planted roots");
        require(acb_contains(out.box.value(), root.value()), name + ": the rectangle misses its own root");
        require(d::meets_goal(out, o, static_cast<slong>(bits + 32)), name + ": the rectangle misses the goal");
    }
    return proved;
}
// The double-double repair on `input` from `starts` (one a planted root, in order): `moving` are repaired, the
// others held. Every moving point must then be proved, its rectangle holding exactly one planted root, no two the
// same; the held points must not move. Returns the points proved before the repair (the starts' own proofs).
uint64_t check_refine(const Planted &input, const std::vector<std::pair<double, double>> &starts,
                      const std::vector<int64_t> &moving, const std::string &name) {
    const d::CompensatedCertificate tier(input.polynomial.value());
    require(tier.available(), name + ": the tier is unavailable");
    const auto o = goal(53);
    const slong n = static_cast<slong>(starts.size());
    d::RootVector points(n);
    for (slong i = 0; i < n; ++i) acb_set_d_d(points[i], starts[static_cast<size_t>(i)].first, starts[static_cast<size_t>(i)].second);
    uint64_t before = 0;
    for (const auto i : moving) {
        d::Root out;
        before += tier.certify(out, points[i], o, 85);
    }
    d::RootVector held(n);
    _acb_vec_set(held.data(), points.data(), n);
    const unsigned steps = certpolysolve::detail::heuristics::double_double_refine_steps;
    uint64_t evaluated = 0;
    require(tier.refine(points.data(), n, moving, 53, evaluated, steps), name + ": refine moved nothing");
    // The charge: a step at each moving point at least, every point at every step at most.
    require(evaluated >= moving.size() && evaluated <= moving.size() * steps,
            name + ": refine counted " + std::to_string(evaluated) + " steps at its points");
    std::vector<int> owner(input.roots.size(), -1);
    for (const auto i : moving) {
        d::Root out;
        require(tier.certify(out, points[i], o, 85), name + ": a repaired point was not proved");
        int held_roots = 0;
        for (size_t r = 0; r < input.roots.size(); ++r)
            if (acb_contains(out.box.value(), input.roots[r].value())) {
                ++held_roots;
                require(owner[r] < 0, name + ": two repaired points proved the same root");
                owner[r] = static_cast<int>(i);
            }
        require(held_roots == 1, name + ": a repaired rectangle holds " + std::to_string(held_roots) + " planted roots");
    }
    for (slong i = 0; i < n; ++i)
        if (std::find(moving.begin(), moving.end(), i) == moving.end())
            require(acb_equal(points[i], held[i]), name + ": a held point moved");
    return before;
}
} // namespace

int main() {
    try {
        std::mt19937_64 rng(20261026);
        // 1. Bounds at random points, near roots and at roots, on planted polynomials of degree 1 to 64.
        uint64_t points = 0;
        for (const int n : {1, 2, 3, 5, 8, 13, 21, 32, 48, 64}) {
            std::uniform_int_distribution<int64_t> mantissa(-(int64_t(1) << 40), int64_t(1) << 40);
            std::vector<std::pair<int64_t, int64_t>> roots;
            for (int j = 0; j < n; ++j) roots.emplace_back(mantissa(rng), j % 3 ? mantissa(rng) : 0);
            const auto input = planted(roots, 40);
            const d::CompensatedCertificate tier(input.polynomial.value());
            require(tier.available(), "planted degree " + std::to_string(n) + ": unavailable");
            std::uniform_real_distribution<double> unit(-2, 2);
            for (const auto &root : input.roots) {
                const double rr = arf_get_d(arb_midref(acb_realref(root.value())), ARF_RND_NEAR);
                const double ri = arf_get_d(arb_midref(acb_imagref(root.value())), ARF_RND_NEAR);
                for (const double t : {0.0, 1e-15, 1e-9, 1e-3}) {
                    check_bounds(tier, input.polynomial.value(), rr * (1 + t), ri * (1 - t), "planted near root");
                    check_second(tier, input.polynomial.value(), rr * (1 + t), ri * (1 - t), 1e-10 * std::hypot(rr, ri),
                                 rng, "planted near root");
                    ++points;
                }
                check_bounds(tier, input.polynomial.value(), unit(rng), unit(rng), "planted, a random point");
                ++points;
            }
        }
        // Dense random coefficients, wide coefficients and three terms, real and complex, at random points.
        for (const int n : {2, 8, 32, 64, 128, 600}) {
            for (const int kind : {0, 1, 2, 3}) {
                d::ComplexPolynomial p;
                std::uniform_real_distribution<double> unit(-1, 1);
                std::uniform_int_distribution<int> wide(-64, 64);
                for (int k = 0; k <= n; ++k) {
                    if (kind == 3 && k != 0 && k != n / 2 && k != n) continue;
                    d::ComplexBall c;
                    const double re = unit(rng), im = kind == 1 ? unit(rng) : 0;
                    acb_set_d_d(c.value(), re, im);
                    if (kind == 2) acb_mul_2exp_si(c.value(), c.value(), wide(rng));
                    acb_poly_set_coeff_acb(p.value(), k, c.value());
                }
                const d::CompensatedCertificate tier(p.value());
                if (!tier.available()) continue;
                for (int k = 0; k < 64; ++k) {
                    const double radius = std::exp2(std::uniform_real_distribution<double>(-3, 3)(rng));
                    const double angle = std::uniform_real_distribution<double>(0, 6.283185307179586)(rng);
                    check_bounds(tier, p.value(), radius * std::cos(angle), radius * std::sin(angle), "random coefficients");
                    ++points;
                }
            }
        }
        // Ball coefficients: the bounds hold for the midpoints' polynomial and for an extreme member.
        {
            d::ComplexPolynomial p, extreme;
            for (int k = 0; k <= 16; ++k) {
                d::ComplexBall c;
                acb_set_si(c.value(), (k * 7919) % 23 - 11);
                acb_poly_set_coeff_acb(extreme.value(), k, c.value());
                d::Magnitude rad;
                mag_set_d(rad.value(), 0x1p-60);
                acb_add_error_mag(c.value(), rad.value());
                acb_poly_set_coeff_acb(p.value(), k, c.value());
            }
            const d::CompensatedCertificate tier(p.value());
            require(tier.available(), "ball coefficients: unavailable");
            for (int k = 0; k < 32; ++k) {
                std::uniform_real_distribution<double> unit(-1.5, 1.5);
                const double xr = unit(rng), xi = unit(rng);
                check_bounds(tier, extreme.value(), xr, xi, "ball coefficients");
                ++points;
            }
        }
        // 2. Proved rectangles hold exactly one planted root.
        uint64_t proved = 0;
        for (const int n : {1, 2, 4, 8, 16, 32, 64}) {
            std::uniform_int_distribution<int64_t> mantissa(-(int64_t(1) << 30), int64_t(1) << 30);
            std::vector<std::pair<int64_t, int64_t>> roots;
            for (int j = 0; j < n; ++j) roots.emplace_back(mantissa(rng), mantissa(rng));
            const auto input = planted(roots, 30);
            for (const double perturbation : {0.0, 0x1p-50, 0x1p-40})
                for (const uint64_t bits : {uint64_t(53), uint64_t(80)})
                    proved += check_planted(input, perturbation, bits, "planted degree " + std::to_string(n), false);
        }
        require(proved > 0, "the tier proved nothing");
        // 2b. The Newton retry (heuristics::compensated_newton_retries): from points too far for the first proof (its
        //     Taylor remainder over the step) but in Newton's quadratic range, certify proves at the point nearest the
        //     step's end. Each such rectangle holds exactly one planted root (check_planted), and some point that
        //     prove alone refuses is proved.
        uint64_t retried = 0;
        for (const int n : {8, 16, 32, 64}) {
            std::uniform_int_distribution<int64_t> mantissa(-(int64_t(1) << 30), int64_t(1) << 30);
            std::vector<std::pair<int64_t, int64_t>> roots;
            for (int j = 0; j < n; ++j) roots.emplace_back(mantissa(rng), mantissa(rng));
            const auto input = planted(roots, 30);
            const d::CompensatedCertificate tier(input.polynomial.value());
            for (const double perturbation : {0x1p-26, 0x1p-32, 0x1p-40}) {
                check_planted(input, perturbation, 53, "the retry at degree " + std::to_string(n), false);
                for (const auto &root : input.roots) {
                    const double xr = arf_get_d(arb_midref(acb_realref(root.value())), ARF_RND_NEAR) * (1 + perturbation);
                    const double xi = arf_get_d(arb_midref(acb_imagref(root.value())), ARF_RND_NEAR) * (1 - perturbation);
                    d::ComplexBall point;
                    acb_set_d_d(point.value(), xr, xi);
                    d::Root first, retry;
                    retried += !tier.prove(first, xr, xi, goal(53), 85) && tier.certify(retry, point.value(), goal(53), 85);
                }
            }
        }
        require(retried > 0, "the retry proved no point that the first proof refuses");
        // 2c. The smaller circles (heuristics::compensated_circle_curvature_share): at a 10-bit goal the goal's circles
        //     fail on the second-derivative term, whose majorant grows with the degree, and circles of a smaller radius
        //     from the same evaluation prove. Every rectangle holds exactly one planted root and meets the goal
        //     (check_planted, at the goal's precision and at the attempt's 85 bits), and some are narrower than half
        //     the goal's half-width, which only the smaller circles give.
        uint64_t narrower = 0;
        for (const int n : {8, 16, 32, 64}) {
            std::uniform_int_distribution<int64_t> mantissa(-(int64_t(1) << 30), int64_t(1) << 30);
            std::vector<std::pair<int64_t, int64_t>> roots;
            for (int j = 0; j < n; ++j) roots.emplace_back(mantissa(rng), mantissa(rng));
            const auto input = planted(roots, 30);
            for (const double perturbation : {0.0, 0x1p-40})
                check_planted(input, perturbation, 10, "the smaller circles at degree " + std::to_string(n), false);
            const d::CompensatedCertificate tier(input.polynomial.value());
            const auto o = goal(10);
            for (const auto &root : input.roots) {
                d::ComplexBall point;
                acb_set_d_d(point.value(), arf_get_d(arb_midref(acb_realref(root.value())), ARF_RND_NEAR),
                            arf_get_d(arb_midref(acb_imagref(root.value())), ARF_RND_NEAR));
                d::Root out;
                if (!tier.certify(out, point.value(), o, 85))
                    continue;
                int held = 0;
                for (const auto &other : input.roots) held += acb_contains(out.box.value(), other.value());
                require(held == 1 && acb_contains(out.box.value(), root.value()),
                        "the smaller circles: a rectangle at 85 bits does not hold exactly its own root");
                require(d::meets_goal(out, o, 85), "the smaller circles: a rectangle at 85 bits misses the goal");
                d::ComplexBall center;
                acb_get_mid(center.value(), out.box.value());
                d::RealBall goal_half;
                d::goal_half_width(goal_half.value(), center.value(), o, 85);
                d::Float twice, lower;
                arf_set_mag(twice.value(), arb_radref(acb_realref(out.box.value())));
                arf_mul_2exp_si(twice.value(), twice.value(), 1);
                arb_get_lbound_arf(lower.value(), goal_half.value(), 85);
                narrower += arf_cmp(twice.value(), lower.value()) < 0;
            }
        }
        require(narrower > 0, "the smaller circles: no rectangle narrower than half the goal's half-width");
        // 3. Where it must decline.
        {
            // A double root at 1: (x - 1)^2 (x + 2).
            const auto input = planted({{1, 0}, {1, 0}, {-2, 0}}, 0);
            const d::CompensatedCertificate tier(input.polynomial.value());
            d::ComplexBall point;
            acb_set_d_d(point.value(), 1 + 0x1p-30, 0);
            d::Root out;
            require(!tier.certify(out, point.value(), goal(53), 85), "a double root was proved simple");
            // Roots 2^-40 apart, which the 53-bit goal's rectangles separate: each one proved holds one.
            const auto near = planted({{(int64_t(1) << 40), 0}, {(int64_t(1) << 40) + 1, 0}, {-3, 0}}, 40);
            check_planted(near, 0, 53, "a pair 2^-40 apart", false);
            // Roots 2^-62 apart, closer than the 53-bit rectangle's half-width 2^-58: no rectangle may be proved
            // about either (it would hold both, and check_planted would fail).
            const auto close = planted({{(int64_t(1) << 62), 0}, {(int64_t(1) << 62) + 1, 0}, {-3, 0}}, 62);
            require(check_planted(close, 0, 53, "a pair 2^-62 apart", false) <= 1, "a pair 2^-62 apart was separated");
            // A goal beyond the two-word precision.
            const auto simple = planted({{3, 1}, {-5, 2}}, 2);
            const d::CompensatedCertificate easy(simple.polynomial.value());
            acb_set_d_d(point.value(), 0.75, 0.25);
            require(!easy.certify(out, point.value(), goal(120), 152), "a 120-bit goal was taken");
        }
        // 3b. may_certify, the necessary condition before the evaluation.
        uint64_t refused = 0, kept = 0;
        {
            // sum (k + 1) (-2x)^k to degree 1100: coefficients from 1 to 1101 2^1100, beyond what binary64 holds
            // at once. Its roots lie near |x| = 1/2, where every term has about the same size, so the coefficients the
            // words leave out defeat Rouché's test at every root: may_certify refuses, and certify with it. The same
            // sum without the scale (1 to 1101, held exactly) is never refused.
            constexpr int n = 1100;
            d::ComplexPolynomial scaled, plain;
            for (int k = 0; k <= n; ++k) {
                d::ComplexBall c;
                acb_set_si(c.value(), (k % 2 ? -1 : 1) * (k + 1));
                acb_poly_set_coeff_acb(plain.value(), k, c.value());
                acb_mul_2exp_si(c.value(), c.value(), k);
                acb_poly_set_coeff_acb(scaled.value(), k, c.value());
            }
            const d::CompensatedCertificate wide(scaled.value()), held(plain.value());
            require(wide.available() && held.available(), "the scaled sums: unavailable");
            // The words hold the exact sum's coefficients, so its points skip the test; the wide sum's carry errors.
            require(!held.carries_errors() && wide.carries_errors(), "the scaled sums: errors misread");
            for (int j = 0; j < 64; ++j) {
                const double angle = 6.283185307179586 * (j + 0.5) / 64;
                require(!wide.may_certify(0.5 * std::cos(angle), 0.5 * std::sin(angle), goal(53)),
                        "the 2^1100-wide sum: a point on its roots' circle was not refused");
                d::ComplexBall point;
                acb_set_d_d(point.value(), 0.5 * std::cos(angle), 0.5 * std::sin(angle));
                d::Root out;
                require(!wide.certify(out, point.value(), goal(53), 85), "the 2^1100-wide sum: a point was proved");
                require(held.may_certify(std::cos(angle), std::sin(angle), goal(53)), "the exact sum: a point was refused");
                ++refused;
            }
        }
        for (const int n : {4, 8, 16, 32, 64}) {
            // Planted roots with 30-bit parts: the coefficients need more than two words, so every one carries an
            // error. may_certify must keep every point that prove proves.
            std::uniform_int_distribution<int64_t> mantissa(-(int64_t(1) << 30), int64_t(1) << 30);
            std::vector<std::pair<int64_t, int64_t>> roots;
            for (int j = 0; j < n; ++j) roots.emplace_back(mantissa(rng), mantissa(rng));
            const auto input = planted(roots, 30);
            const d::CompensatedCertificate tier(input.polynomial.value());
            require(tier.available(), "planted for may_certify: unavailable");
            require(tier.carries_errors(), "planted for may_certify: coefficients beyond two words carry no error");
            for (const auto &root : input.roots)
                for (const double perturbation : {0.0, 0x1p-50, 0x1p-40, 0x1p-20})
                    for (const uint64_t bits : {uint64_t(53), uint64_t(80)}) {
                        const double xr = arf_get_d(arb_midref(acb_realref(root.value())), ARF_RND_NEAR) * (1 + perturbation);
                        const double xi = arf_get_d(arb_midref(acb_imagref(root.value())), ARF_RND_NEAR) * (1 - perturbation);
                        d::Root out;
                        if (!tier.prove(out, xr, xi, goal(bits), static_cast<slong>(bits + 32)))
                            continue;
                        require(tier.may_certify(xr, xi, goal(bits)), "may_certify refused a point prove proves");
                        ++kept;
                    }
        }
        require(kept > 0, "may_certify: nothing proved to keep");
        // 4. Through the public API: Auto's attempt uses the tier, and each box holds exactly one root of a 400-bit
        //    solve's boxes (proved by the ball certificate).
        uint64_t checked = 0, tier_boxes = 0;
        const auto contains =[](const certpolysolve::Box &outer, const certpolysolve::Box &inner) {
            return outer.real.lo <= inner.real.lo && inner.real.hi <= outer.real.hi && outer.imag.lo <= inner.imag.lo &&
                   inner.imag.hi <= outer.imag.hi;
        };
        const auto overlaps = [](const certpolysolve::Box &a, const certpolysolve::Box &b) {
            return !(a.real.hi < b.real.lo || b.real.hi < a.real.lo || a.imag.hi < b.imag.lo || b.imag.hi < a.imag.lo);
        };
        for (const int n : {3, 8, 24, 64}) {
            certpolysolve::Polynomial p;
            std::uniform_int_distribution<int> coefficient(-1000, 1000);
            for (int k = 0; k <= n; ++k)
                p.set(static_cast<uint64_t>(k), certpolysolve::Number(std::to_string(coefficient(rng) | 1)));
            const auto fast = certpolysolve::solve(p);
            const auto reference = certpolysolve::solve(p, {.bits = 400, .goal = certpolysolve::Goal::accuracy()});
            require(fast.goal_met() && reference.goal_met(), "degree " + std::to_string(n) + ": a solve failed");
            for (const auto &box : fast.roots()) {
                int inside = 0, meeting = 0;
                for (const auto &r : reference.roots()) {
                    if (contains(box.box, r.box)) ++inside;
                    else if (overlaps(box.box, r.box)) ++meeting;
                }
                require(inside == 1 && meeting == 0,
                        "degree " + std::to_string(n) + ": a box does not hold exactly one root");
                ++checked;
                tier_boxes += box.certificate.rfind("original-coefficient binary64 compensated", 0) == 0;
            }
        }
        require(tier_boxes > 0, "no public answer was proved by the tier");
        // 5. The double-double repair (proposal 1).
        {
            // A real triple cluster 2^-18 wide at 1 (1 - 2^-18, 1, 1 + 2^-18), with -2 and 3: its starts off by about
            // 2^-17 (binary64's u / s^2 for it), two of them a conjugate pair.
            const auto triple = planted({{(int64_t(1) << 18) - 1, 0}, {int64_t(1) << 18, 0}, {(int64_t(1) << 18) + 1, 0},
                                         {-(int64_t(2) << 18), 0}, {int64_t(3) << 18, 0}}, 18);
            const uint64_t triple_before =
                check_refine(triple, {{1 - 0x1p-18 + 0x1p-17, 0x1p-17}, {1 + 0x1p-17, -0x1p-17}, {1 + 0x1p-18 - 0x1p-17, 0},
                                      {-2, 0}, {3, 0}},
                             {0, 1, 2}, "a triple cluster 2^-18 wide");
            require(triple_before < 3, "the triple cluster's poor starts were all proved without the repair");
            // A real pair 2^-40 apart at 1, with -3, started as a conjugate pair 2^-30 off the axis.
            const auto pair = planted({{int64_t(1) << 40, 0}, {(int64_t(1) << 40) + 1, 0}, {-(int64_t(3) << 40), 0}}, 40);
            const uint64_t pair_before = check_refine(pair, {{1, 0x1p-30}, {1, -0x1p-30}, {-3, 0}}, {0, 1}, "a pair 2^-40 apart");
            require(pair_before == 0, "the pair's conjugate starts were proved without the repair");
        }
        std::cout << "compensated certificate: " << points << " evaluations bounded, " << proved
                  << " planted roots proved, each rectangle with exactly one root; " << retried
                  << " points proved by the Newton retry alone; declines held; may_certify refused "
                  << refused << " points of the 2^1100-wide sum and kept all " << kept << " points prove proves; " << checked
                  << " public boxes checked against 400-bit roots, " << tier_boxes << " of them the tier's\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << "test_compensated_certificate: " << e.what() << '\n';
        return 1;
    }
}
