#include "numeric/arithmetic/polynomial_arithmetic.hpp"
#include "internal.hpp"
#include "exact_winding_api.hpp"
#include "api_test_support.hpp"
#include "numeric/candidates/candidates.hpp"
#include "numeric/exact/centering.hpp"
#include "numeric/evaluate/evaluator.hpp"
#include "secular_evaluator.hpp"
#include "input/secular/secular.hpp"
#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
#include "solve/factor_state.hpp"
#endif
#include <cmath>
#include <complex>
#include <iostream>
#include <numbers>
#include <utility>

using namespace reference;
namespace d = certpolysolve::detail;
static Complex operator/(const Complex &a, const Complex &b) {
    const Rational norm = b.real * b.real + b.imag * b.imag;
    return {(a.real * b.real + a.imag * b.imag) / norm, (a.imag * b.real - a.real * b.imag) / norm};
}

static Complex value(const ComplexPolynomial &p, const Complex &z) {
    Complex out;
    for (auto i = p.rbegin(); i != p.rend(); ++i)
        out = out * z + *i;
    return out;
}
static void enclose(acb_t out, const Complex &z, slong prec = 128) {
    d::Number(z.real.str()).enclose(acb_realref(out), prec);
    d::Number(z.imag.str()).enclose(acb_imagref(out), prec);
}
static bool contains(const acb_t box, const Complex &z) {
    d::Rational re, im;
    require(fmpq_set_str(re.value(), z.real.str().c_str(), 10) == 0, "bad reference literal");
    require(fmpq_set_str(im.value(), z.imag.str().c_str(), 10) == 0, "bad reference literal");
    return arb_contains_fmpq(acb_realref(box), re.value()) && arb_contains_fmpq(acb_imagref(box), im.value());
}
static void polynomial(acb_poly_t out, const ComplexPolynomial &p, slong prec = 128) {
    d::ComplexBall c;
    for (size_t i = 0; i < p.size(); ++i) {
        enclose(c.value(), p[i], prec);
        acb_poly_set_coeff_acb(out, static_cast<slong>(i), c.value());
    }
}

// Deliberately test-only: public oracle solving also needs counts, bounds and
// multiplicity contracts, which this evaluation interface alone cannot supply.
class MandelbrotEvaluator final : public d::Evaluator {
    unsigned levels_;

  public:
    explicit MandelbrotEvaluator(unsigned levels) : levels_(levels) {}
    slong degree() const override { return (slong(1) << levels_) - 1; }
    bool evaluate(acb_t v, acb_t derivative, const acb_t z, slong prec) const override {
        acb_one(v);
        acb_zero(derivative);
        d::ComplexBall square, next;
        for (unsigned i = 0; i < levels_; ++i) {
            acb_mul(next.value(), v, derivative, prec);
            acb_mul(next.value(), next.value(), z, prec);
            acb_mul_2exp_si(next.value(), next.value(), 1);
            acb_sqr(square.value(), v, prec);
            acb_add(derivative, next.value(), square.value(), prec);
            acb_mul(v, z, square.value(), prec);
            acb_add_ui(v, v, 1, prec);
        }
        return acb_is_finite(v) && acb_is_finite(derivative);
    }
};
static void evaluators() {
    for (bool sparse : {false, true}) {
        ComplexPolynomial p(33);
        p[0] = {1, 2};
        p[1] = {-2, 1};
        p[17] = {3, -4};
        p[32] = {1};
        d::ComplexPolynomial enclosed;
        polynomial(enclosed.value(), p);
        d::CoefficientEvaluator evaluator(enclosed.value(), sparse);
        require(evaluator.degree() == 32, "evaluator lost degree");
        ComplexPolynomial derivative(32);
        for (size_t i = 1; i < p.size(); ++i)
            derivative[i - 1] = p[i] * Complex(Rational(i));
        for (const Complex center : {Complex(0), Complex(Rational(1) / 2, Rational(1) / 4), Complex(-1)}) {
            d::ComplexBall z, v, dv;
            enclose(z.value(), center);
            arb_add_error_2exp_si(acb_realref(z.value()), -8);
            arb_add_error_2exp_si(acb_imagref(z.value()), -9);
            require(evaluator.evaluate(v.value(), dv.value(), z.value(), 128), "finite evaluation failed");
            for (int i : {-1, 0, 1})
                for (int j : {-1, 0, 1}) {
                    Complex sample = center + Complex(Rational(i) / 256, Rational(j) / 512);
                    require(contains(v.value(), value(p, sample)), "value enclosure lost exact sample");
                    require(contains(dv.value(), value(derivative, sample)),
                            "derivative enclosure lost exact sample");
                }
        }
    }
    Polynomial p{1};
    for (unsigned level = 0; level <= 5; ++level) {
        MandelbrotEvaluator evaluator(level);
        require(evaluator.degree() == degree(p), "recurrence degree mismatch");
        for (const Rational z : {Rational(0), Rational(-1), Rational(1) / 3}) {
            d::ComplexBall x, v, dv;
            enclose(x.value(), Complex(z));
            require(evaluator.evaluate(v.value(), dv.value(), x.value(), 192), "recurrence evaluation failed");
            require(contains(v.value(), Complex(evaluate(p, z))),
                    "recurrence value differs from exact expansion");
            require(contains(dv.value(), Complex(evaluate(derivative(p), z))),
                    "recurrence derivative differs from exact expansion");
        }
        p = multiply(p, p);
        p.insert(p.begin(), 1);
    }
    d::ComplexPolynomial zero;
    d::CoefficientEvaluator evaluator(zero.value());
    d::ComplexBall x, v, dv;
    require(evaluator.evaluate(v.value(), dv.value(), x.value(), 64) && acb_is_zero(v.value()) &&
                acb_is_zero(dv.value()),
            "zero polynomial evaluation");
}
static void secular_identity() {
    d::ComplexPolynomial p;
    polynomial(p.value(), {{2}, {-3}, {1}});
    std::vector<d::ComplexBall> nodes(2), weights;
    acb_set_si(nodes[1].value(), 3);
    cps_options options;
    certpolysolve::detail::CandidateStatistics stats{};
    require(d::regenerate_secular(weights, p.value(), nodes, 128, options, stats),
            "secular regeneration failed");
    require(contains(weights[0].value(), Complex(Rational(2) / 3)) &&
                contains(weights[1].value(), Complex(Rational(-2) / 3)),
            "wrong secular sign or normalization");
    d::SecularEvaluator evaluator(weights, nodes);
    for (const Complex sample :
         {Complex(1), Complex(2), Complex(0, 1), Complex(Rational(1) / 2, Rational(1) / 3)}) {
        d::ComplexBall z, v, dv;
        enclose(z.value(), sample);
        require(evaluator.evaluate(v.value(), dv.value(), z.value(), 128),
                "secular evaluation outside poles failed");
        const Complex a(Rational(2) / 3), b(Rational(-2) / 3);
        const auto expected = a / sample + b / (sample - Complex(3)) - Complex(1);
        const auto slope =
            Complex(0) - a / (sample * sample) - b / ((sample - Complex(3)) * (sample - Complex(3)));
        require(contains(v.value(), expected) && contains(dv.value(), slope),
                "secular value/derivative differs from rational reference");
    }
    d::ComplexBall z, v, dv;
    require(!evaluator.evaluate(v.value(), dv.value(), z.value(), 128), "a pole was treated as a finite value");
    acb_zero(nodes[1].value());
    require(!d::regenerate_secular(weights, p.value(), nodes, 128, options, stats),
            "coincident interpolation nodes accepted");
    // Secular packets start exactly on their regeneration nodes. Successful
    // convergence here exercises the removable polynomial Newton limit there.
    acb_set_si(nodes[1].value(), 3);
    d::RootVector initial(2), candidates(2);
    acb_zero(initial.data());
    acb_set_si(initial.data() + 1, 3);
    auto progress = d::candidate_packet(candidates.data(), p.value(), initial.data(), 128, 16,
                                        d::CandidateEngine::Secular, options, stats);
    require(progress.finite && stats.secular_regenerations > 0, "pole-limit candidate step failed");
    require(_acb_poly_validate_roots(candidates.data(), p.value()->coeffs, 3, 128) == 2,
            "regenerated candidates failed original-polynomial certification");
}
static void known(const certpolysolve::Result &r, const std::vector<Complex> &roots, int bits) {
    if (!r.complete() || !r.goal_met())
        throw std::runtime_error("known roots incomplete: " + r.message());
    check_disjoint(r.roots());
    size_t total = 0;
    for (const auto &region : r.roots()) {
        size_t count = 0;
        for (const auto &z : roots)
            if (endpoint(region.box.real.lo) <= z.real && z.real <= endpoint(region.box.real.hi) &&
                endpoint(region.box.imag.lo) <= z.imag && z.imag <= endpoint(region.box.imag.hi))
                ++count;
        require(count == region.count && count == region.multiplicity.value_or(0),
                "exact known-root multiplicity mismatch");
        total += count;
        check_accuracy(region, bits);
    }
    require(total == roots.size(), "known roots missing");
}
static void native_tree_candidates() {
    Polynomial exact{1};
    for (int j = 1; j <= 128; ++j)
        exact = multiply(exact, {-j, 1});
    ComplexPolynomial complex;
    for (const auto &c : exact)
        complex.emplace_back(c);
    d::ComplexPolynomial p;
    polynomial(p.value(), complex, 1024);
    d::RootVector seeds(128), roots(128);
    for (slong j = 0; j < 128; ++j)
        enclose(seeds[j], Complex(Rational(j + 1) + power2(-10), power2(-1200)), 2048);
    cps_options options;
    certpolysolve::detail::CandidateStatistics statistics{};
    const auto progress = d::candidate_packet(roots.data(), p.value(), seeds.data(), 1024, 16,
                                              d::CandidateEngine::Aberth, options, statistics, {}, true);
    require(progress.finite, "tree proposals became nonfinite after coordinate underflow");
    require(statistics.root_interactions < statistics.candidate_sweeps * 128 * 127,
            "native tree did not reduce direct interactions");
    require(_acb_poly_validate_roots(roots.data(), p.value()->coeffs, 129, 1024) == 128,
            "native tree candidates did not certify against the original polynomial");
    for (int j = 1; j <= 128; ++j) {
        unsigned matches = 0;
        for (slong k = 0; k < 128; ++k)
            matches += contains(roots[k], Complex(j));
        require(matches == 1, "tree certificate disagrees with independent exact integer roots");
    }
}
static void stalled_exact_source() {
    // Two close pairs near +/-sqrt(3/a), plus five distant roots. The
    // coefficients are exact; Boost's winding calculation checks the exported
    // counts without using the solver's arithmetic or inclusion formula.
    const Rational a = power2(64);
    Polynomial exact(10);
    exact[0] = 9;
    exact[2] = -6 * a;
    exact[4] = a * a;
    exact[9] = a;
    certpolysolve::Options options;
    options.bits = 53;
    options.max_bits = 4096;
    const auto result = certpolysolve::solve(input(exact), options);
    require(result.complete() && result.goal_met(), "stalled exact source did not recover");
    ComplexPolynomial complex;
    for (const auto &coefficient : exact)
        complex.emplace_back(coefficient);
    check_disjoint(result.roots());
    uint64_t total = 0;
    for (const auto &root : result.roots()) {
        const auto count = rectangle_count(complex, root);
        require(count && *count == static_cast<int>(root.count), "stalled source winding disagreement");
        check_accuracy(root, 53);
        total += root.count;
    }
    require(total == 9, "stalled source lost roots");
}
static void candidates() {
#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
    // The DoubleSeed test engine's binary64 start (research build).
    {
        d::ComplexPolynomial p;
        acb_poly_set_coeff_si(p.value(), 0, -2);
        acb_poly_set_coeff_si(p.value(), 2, 1);
        d::RootVector seeds(2);
        cps_options o;
        certpolysolve::detail::CandidateStatistics statistics{};
        require(d::binary64_seeds(seeds.data(), p.value(), o, statistics), "binary64 initialization failed");
        require(_acb_poly_validate_roots(seeds.data(), p.value()->coeffs, 3, 128) == 2,
                "binary64 proposals failed original coefficient validation");
        d::ComplexBall tiny;
        acb_one(tiny.value());
        acb_mul_2exp_si(tiny.value(), tiny.value(), -4000);
        acb_poly_set_coeff_acb(p.value(), 1, tiny.value());
        require(!d::binary64_seeds(seeds.data(), p.value(), o, statistics),
                "binary64 range loss did not request fallback");
    }
#endif
    // Legacy, DoubleSeed and Hybrid are test engines of the research build.
    for (auto engine : {
#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
             d::CandidateEngine::Legacy, d::CandidateEngine::DoubleSeed,
#endif
             d::CandidateEngine::Aberth, d::CandidateEngine::Secular, d::CandidateEngine::Centered,
#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
             d::CandidateEngine::Hybrid,
#endif
         }) {
        for (const std::vector<Complex> roots :
             {std::vector<Complex>{{-2}, {Rational(1) / 3}, {2, 1}, {2, -1}},
              std::vector<Complex>{{0}, {0}, {1}, {1 + power2(-120)}, {0, 1}, {0, 1}},
              std::vector<Complex>{{power2(300)}, {power2(-300)}, {-power2(300)}}}) {
            // A one-unit solve stops honestly, and one with the 2,000 units the one-unit calls of a session used
            // to share meets the goal (one advance a solve state since 18 October 2026).
            auto p = input(from_roots(roots));
            cps_options o;
            o.bits = 53;
            o.max_work_bits = 8192;
            o.work_units = 1;
            d::SolveState unit(certpolysolve::detail::dense(*p.get()), {.backend = true, .candidates = engine});
            const auto stopped = unit.advance(&o);
            require(stopped.goal_met || stopped.termination == CPS_WORK_LIMIT, "a one-unit solve stopped for another reason");
            o.work_units = 2000;
            d::SolveState whole(certpolysolve::detail::dense(*p.get()), {.backend = true, .candidates = engine});
            known(d::make_result(CPS_OK, d::published_result(whole.advance(&o))), roots, 53);
        }
    }
}
// A region's exact endpoint, the dyadic m 2^e.
static void exact_endpoint(arf_t out, const certpolysolve::Dyadic &x) {
    d::Integer mantissa, exponent;
    require(fmpz_set_str(mantissa.value(), x.mantissa().c_str(), 10) == 0 &&
                fmpz_set_str(exponent.value(), x.exponent().c_str(), 10) == 0,
            "bad dyadic endpoint");
    arf_set_fmpz_2exp(out, mantissa.value(), exponent.value());
}
// Whether the closed region holds every point of the ball z: its exact endpoints against z's outward bounds.
static bool holds(const certpolysolve::Box &box, const acb_t z) {
    const auto within = [](const certpolysolve::Interval &interval, const arb_t x) {
        d::Float lower, upper, bound;
        exact_endpoint(lower.value(), interval.lo);
        exact_endpoint(upper.value(), interval.hi);
        arb_get_lbound_arf(bound.value(), x, 256);
        if (arf_cmp(bound.value(), lower.value()) < 0) return false;
        arb_get_ubound_arf(bound.value(), x, 256);
        return arf_cmp(bound.value(), upper.value()) <= 0;
    };
    return within(box.real, acb_realref(z)) && within(box.imag, acb_imagref(z));
}
static void analytic_and_overlap() {
#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
    // The research build's Hybrid test engine takes binomials in closed form at every degree.
    for (int n : {3, 4, 7, 12})
        for (const Complex constant : {Complex(-2), Complex(2), Complex(1, 2)}) {
            ComplexPolynomial p(static_cast<size_t>(n + 1));
            p[0] = constant;
            p.back() = Complex(1);
            auto input_poly = input(p);
            d::SolveState state(certpolysolve::detail::dense(*input_poly.get()),
                                {.backend = true, .candidates =d::CandidateEngine::Hybrid});
            cps_options o;
            o.bits = 80;
            auto result = d::make_result(CPS_OK, d::published_result(state.advance(&o)));
            require(result.goal_met() && result.roots().size() == static_cast<size_t>(n),
                    "binomial cover did not complete");
            check_disjoint(result.roots());
            for (const auto &root : result.roots()) {
                check_accuracy(root, 80);
                if (root.realness == certpolysolve::Realness::real) {
                    Polynomial real(static_cast<size_t>(n + 1));
                    real[0] = constant.real;
                    real.back() = 1;
                    require(constant.imag == 0 &&
                                roots_closed(real, endpoint(root.box.real.lo), endpoint(root.box.real.hi)) == 1,
                            "binomial real root count mismatch");
                } else if (endpoint(root.box.real.lo) == 0 && endpoint(root.box.real.hi) == 0 && n % 2 == 0 &&
                           constant.imag == 0) {
                    Polynomial imaginary(static_cast<size_t>(n + 1));
                    imaginary[0] = constant.real;
                    imaginary.back() = n % 4 == 0 ? 1 : -1;
                    require(roots_closed(imaginary, endpoint(root.box.imag.lo), endpoint(root.box.imag.hi)) ==
                                1,
                            "imaginary-axis Sturm check failed");
                } else if (!rouche_inside_rectangle(p, root)) {
                    const auto count =
                        rectangle_count(p, endpoint(root.box.real.lo), endpoint(root.box.real.hi),
                                        endpoint(root.box.imag.lo), endpoint(root.box.imag.hi));
                    require(count && *count == 1, "independent binomial winding check failed");
                }
            }
        }
    // Exercise interval coefficients and rotation away from the log branch cut.
    certpolysolve::Polynomial uncertain;
    uncertain.set(0, certpolysolve::ball(certpolysolve::Number("2"), certpolysolve::Number("1/100")));
    uncertain.set(7, certpolysolve::Number("1"));
    d::SolveState family(certpolysolve::detail::dense(*uncertain.get()),
                         {.backend = true, .candidates =d::CandidateEngine::Hybrid});
    cps_options o;
    o.max_work_bits = 128;
    o.absolute_tolerance = d::Number("1/100");
    auto result = d::make_result(CPS_OK, d::published_result(family.advance(&o)));
    require(result.goal_met(), "uncertain binomial cover failed");
#endif
    // The maintained route to the closed form (binomial_enclosures, isolate_complex): Auto sends inexact
    // coefficients to the D5 recipe (production_policy()), whose classical pass takes a binomial of degree 128 or
    // more as an nth root of the ratio times the nth roots of unity (an exact binomial is a power composition, solved
    // through a linear factor). The constants' signs give the unrotated root with two real roots and the rotated one
    // with one. As c ranges over the ball, the root |c|^(1/n) e^(i pi (a + 2k) / n) of x^n = -c (a = 1 for c > 0)
    // moves along a radial segment, and a region (convex) holds it when it holds the segment's two ends, enclosed in
    // ball arithmetic: each of the n disjoint regions must hold its own one. (Until W8's S2 the default suite reached
    // the closed form only through the Hybrid case above.)
    const slong precision = 256;
    for (const auto &[n, c] : {std::pair<slong, slong>{128, -2}, std::pair<slong, slong>{129, 2}}) {
        certpolysolve::Polynomial binomial;
        binomial.set(0, certpolysolve::ball(certpolysolve::Number(std::to_string(c)),
                                        certpolysolve::Number("1/1267650600228229401496703205376"))); // 2^-100
        binomial.set(static_cast<size_t>(n), certpolysolve::Number("1"));
        const auto answer = certpolysolve::solve(binomial);
        require(answer.goal_met() && answer.complete() && answer.roots().size() == static_cast<size_t>(n),
                ("the inexact binomial was not solved: " + answer.message()).c_str());
        check_disjoint(answer.roots());
        d::RealBall low, high; // the segments' ends' moduli, (|c| -+ 2^-100)^(1/n)
        for (auto *end : {&low, &high}) {
            arb_one(end->value());
            arb_mul_2exp_si(end->value(), end->value(), -100);
            if (end == &low) arb_neg(end->value(), end->value());
            arb_add_si(end->value(), end->value(), c < 0 ? -c : c, precision);
            arb_root_ui(end->value(), end->value(), static_cast<ulong>(n), precision);
        }
        const slong a = c > 0;
        std::vector<bool> taken(static_cast<size_t>(n));
        for (const auto &root : answer.roots()) {
            require(root.certificate == "binomial-nth-root-and-unity-cover",
                    ("the inexact binomial did not take the closed form: " + root.certificate).c_str());
            require(root.count == 1, "a binomial region claims more than one root");
            check_accuracy(root, 53);
            const auto k = ((std::lround(std::arg(root.centre()) / std::numbers::pi * static_cast<double>(n) / 2 -
                                         static_cast<double>(a) / 2) %
                             n) + n) % n;
            require(!taken[static_cast<size_t>(k)], "two regions for one root of the binomial");
            taken[static_cast<size_t>(k)] = true;
            d::Rational angle; // in units of pi
            fmpq_set_si(angle.value(), a + 2 * k, static_cast<ulong>(n));
            d::RealBall sine, cosine;
            arb_sin_cos_pi_fmpq(sine.value(), cosine.value(), angle.value(), precision);
            for (const auto *modulus : {&low, &high}) {
                d::ComplexBall end;
                arb_mul(acb_realref(end.value()), cosine.value(), modulus->value(), precision);
                arb_mul(acb_imagref(end.value()), sine.value(), modulus->value(), precision);
                require(holds(root.box, end.value()), "a binomial region misses its root");
            }
        }
    }
    std::vector<d::Root> boxes(160);
    for (size_t i = 0; i < boxes.size(); ++i) {
        acb_set_si(boxes[i].box.value(), static_cast<slong>(3 * i));
        arb_add_error_2exp_si(acb_realref(boxes[i].box.value()), -1);
    }
    require(d::disjoint(boxes), "sweep rejected separated boxes");
    acb_set(boxes.back().box.value(), boxes.front().box.value());
    arb_add_ui(acb_realref(boxes.back().box.value()), acb_realref(boxes.back().box.value()), 1, 64);
    require(!d::disjoint(boxes), "closed touching boxes were pruned by the sweep");
    arb_add_ui(acb_realref(boxes.back().box.value()), acb_realref(boxes.back().box.value()), 1, 64);
    require(!d::disjoint(boxes), "touching next neighbor was missed");
}
static void exact_centering() {
    const Polynomial p{2, -3, 5, 0, 7};
    d::IntegerPolynomial source, translated;
    for (size_t i = 0; i < p.size(); ++i)
        fmpz_poly_set_coeff_si(source.value(), static_cast<slong>(i), p[i].convert_to<slong>());
    for (const Rational c : {Rational(-1), Rational(2) / 7, Rational(-13) / 8}) {
        d::Rational center;
        fmpq_set_str(center.value(), c.str().c_str(), 10);
        require(d::center_integer(translated.value(), source.value(), center.value(), 100000),
                "exact centering refused a small input");
        Polynomial q;
        for (slong i = 0; i <= fmpz_poly_degree(translated.value()); ++i)
            q.emplace_back(d::integer_string(translated.value()->coeffs + i));
        // Compare values up to the removed scalar, using independent rationals.
        for (const Rational x : {Rational(-2), Rational(0), Rational(1) / 3, Rational(7)})
            require(evaluate(q, x) * p.back() == evaluate(p, x + c) * q.back(),
                    "exact shift changed the polynomial");
        require(!d::center_integer(translated.value(), source.value(), center.value(), 1),
                "centering ignored its preflight budget");
    }
    const ComplexPolynomial cp{{Rational(2) / 7, Rational(3) / 2}, {-3, 1}, {5, -2}, {0, 3}, {7, 2}};
    d::GaussianPolynomial source_gaussian(cp.size()), shifted;
    for (size_t i = 0; i < cp.size(); ++i) {
        fmpq_set_str(source_gaussian[i].real.value(), cp[i].real.str().c_str(), 10);
        fmpq_set_str(source_gaussian[i].imag.value(), cp[i].imag.str().c_str(), 10);
    }
    for (const Complex c : {Complex(-1), Complex(Rational(2) / 7, Rational(-13) / 8)}) {
        d::GaussianRational center;
        fmpq_set_str(center.real.value(), c.real.str().c_str(), 10);
        fmpq_set_str(center.imag.value(), c.imag.str().c_str(), 10);
        require(d::center_gaussian(shifted, source_gaussian, center, 100000),
                "Gaussian centering refused small polynomial");
        ComplexPolynomial q;
        for (const auto &z : shifted)
            q.emplace_back(Rational(d::integer_string(fmpq_numref(z.real.value()))),
                           Rational(d::integer_string(fmpq_numref(z.imag.value()))));
        for (const Complex x :
             {Complex(-2), Complex(0), Complex(Rational(1) / 3, Rational(2) / 7), Complex(7)})
            require(zero(value(q, x) * cp.back() - value(cp, x + c) * q.back()),
                    "Gaussian shift changed exact polynomial");
        require(!d::center_gaussian(shifted, source_gaussian, center, 1),
                "Gaussian centering ignored budget");
    }
}
static void sparse_translation_budget() {
    // The degree-383 sparse cluster polynomial has a large highest shifted
    // coefficient but graded storage growth. Its full exact shift fits 8M
    // bits; charging that largest coefficient 384 times wrongly rejects it.
    constexpr slong degree = 383;
    d::IntegerPolynomial source, pair, translated;
    d::Integer scale, coefficient;
    fmpz_one(scale.value());
    fmpz_mul_2exp(scale.value(), scale.value(), 64);
    fmpz_poly_set_coeff_si(pair.value(), 0, -1);
    fmpz_mul(coefficient.value(), scale.value(), scale.value());
    fmpz_poly_set_coeff_fmpz(pair.value(), 2, coefficient.value());
    fmpz_poly_pow(source.value(), pair.value(), 6);
    fmpz_poly_scalar_mul_si(source.value(), source.value(), -2);
    fmpz_poly_set_coeff_si(source.value(), degree, 1);
    d::Rational center;
    fmpq_one(center.value());
    fmpq_div_2exp(center.value(), center.value(), 64);
    require(d::center_integer(translated.value(), source.value(), center.value(), 8000000),
            "graded exact translation budget rejected a fitting polynomial");
    require(fmpz_is_one(translated.value()->coeffs), "local constant lost exact cancellation");
    uint64_t bits = 0;
    Integer at_one = 0;
    for (slong j = 0; j <= degree; ++j) {
        bits += static_cast<uint64_t>(fmpz_bits(translated.value()->coeffs + j));
        at_one += Integer(d::integer_string(translated.value()->coeffs + j));
    }
    require(bits <= 8000000, "translation exceeded its admitted coefficient budget");
    const Integer b = Integer(1) << 64;
    const Integer expected =
        boost::multiprecision::pow(Integer(b + 1), degree) -
        2 * boost::multiprecision::pow(b, degree + 6) * boost::multiprecision::pow(Integer(b + 2), 6);
    require(at_one == expected, "large exact translation differs from independent formula");
    require(!d::center_integer(translated.value(), source.value(), center.value(), 1024),
            "large exact translation ignored a small budget");
}
static void rejected_translation_cache() {
    d::IntegerPolynomial p;
    fmpz_poly_set_coeff_si(p.value(), 0, -1);
    fmpz_poly_set_coeff_si(p.value(), 3, 1);
    d::ExactTaylorCache cache(p, 100000);
    d::ComplexPolynomial out;
    d::ComplexBall center;
    acb_one(center.value());
    require(cache.enclose(out.value(), center.value(), 128), "small cached translation failed");
    acb_mul_2exp_si(center.value(), center.value(), 10000000);
    require(!cache.enclose(out.value(), center.value(), 128) && !cache.enclose(out.value(), center.value(), 256),
            "rejected translation reused a polynomial from another centre");
    acb_one(center.value());
    require(cache.enclose(out.value(), center.value(), 192) && acb_is_zero(out.value()->coeffs),
            "translation cache did not recover after a budget rejection");
}
#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
static void exact_translation_service_work() {
    d::IntegerPolynomial p, shifted, plain, scaled, scaled_shifted;
    for (slong j = 0; j < 4; ++j)
        fmpz_poly_set_coeff_si(p.value(), j, j + 1);
    d::Rational one;
    fmpq_one(one.value());
    double real_work = 0;
    require(d::center_integer(shifted.value(), p.value(), one.value(), 100000, nullptr, &real_work) &&
            d::center_integer(plain.value(), p.value(), one.value(), 100000) &&
            fmpz_poly_equal(shifted.value(), plain.value()) && real_work > 0,
            "research exact shift changed the polynomial or omitted its work");
    // Content cancellation leaves the same primitive output, but the actual
    // integer operands at the shift boundary are substantially larger.
    fmpz_poly_set(scaled.value(), p.value());
    for (slong j = 0; j < fmpz_poly_length(scaled.value()); ++j)
        fmpz_mul_2exp(scaled.value()->coeffs + j, scaled.value()->coeffs + j, 512);
    double scaled_work = 0;
    require(d::center_integer(scaled_shifted.value(), scaled.value(), one.value(), 100000, nullptr, &scaled_work) &&
            fmpz_poly_equal(shifted.value(), scaled_shifted.value()) && scaled_work > real_work,
            "exact service work ignored operand bits removed by later content division");
    d::Rational large;
    fmpq_one(large.value());
    fmpz_mul_2exp(fmpq_numref(large.value()), fmpq_numref(large.value()), 512);
    double growing_work = 0;
    require(d::center_integer(plain.value(), p.value(), large.value(), 100000, nullptr, &growing_work) &&
            growing_work > real_work,
            "exact service work ignored growth during a shift");
    // Two identical real components execute two distinct integer shifts.
    d::GaussianPolynomial gaussian(4), gaussian_shifted;
    for (size_t j = 0; j < gaussian.size(); ++j) {
        fmpq_set_fmpz(gaussian[j].real.value(), p.value()->coeffs + static_cast<slong>(j));
        fmpq_set_fmpz(gaussian[j].imag.value(), p.value()->coeffs + static_cast<slong>(j));
    }
    d::GaussianRational center;
    fmpq_one(center.real.value());
    double gaussian_work = 0;
    require(d::center_gaussian(gaussian_shifted, gaussian, center, 100000, nullptr, &gaussian_work) &&
            gaussian_work == 2 * real_work,
            "exact Gaussian service work did not count both executed real-component shifts");
    fmpq_zero(center.real.value());
    fmpq_one(center.imag.value());
    double imaginary_work = 0;
    require(d::center_gaussian(gaussian_shifted, gaussian, center, 100000, nullptr, &imaginary_work) &&
            imaginary_work > 0,
            "imaginary exact translation omitted its real-component shifts");

    d::ExactTaylorCache cache(p, 100000), unmetered(p, 100000);
    cache.research_service_meter(true);
    d::ComplexPolynomial out, reference;
    d::ComplexBall target;
    acb_one(target.value());
    require(cache.enclose(out.value(), target.value(), 128) &&
            unmetered.enclose(reference.value(), target.value(), 128) &&
            acb_poly_equal(out.value(), reference.value()), "metered cache changed its first enclosure");
    const double first_work = cache.research_work();
    require(first_work > 0 && unmetered.research_work() == 0, "exact service meter is not opt-in");
    require(cache.enclose(out.value(), target.value(), 256) &&
            unmetered.enclose(reference.value(), target.value(), 256) &&
            acb_poly_equal(out.value(), reference.value()), "metered precision reuse changed its enclosure");
    const double enclosure_work = cache.research_work() - first_work;
    require(enclosure_work > 0 && enclosure_work < first_work,
            "exact precision reuse omitted enclosure work or charged another construction");
    acb_mul_2exp_si(target.value(), target.value(), 10000000);
    const double before_refusal = cache.research_work();
    require(!cache.enclose(out.value(), target.value(), 128) &&
            !cache.enclose(out.value(), target.value(), 256) && cache.research_work() == before_refusal,
            "a preflight or cached refusal was charged as an executed exact shift");
    acb_one(target.value());
    require(cache.enclose(out.value(), target.value(), 128) &&
            cache.research_work() - before_refusal > enclosure_work,
            "returning to an evicted exact center did not charge its reconstruction");
}

static void exact_translation_reader_lifetime() {
    d::IntegerPolynomial p;
    fmpz_poly_set_coeff_si(p.value(), 0, -1);
    fmpz_poly_set_coeff_si(p.value(), 3, 1);
    d::PolynomialSourceOptions options;
    {
        d::FactorState factor;
        factor.kind = d::FactorState::Real;
        factor.integer = p;
        options.exact_taylor = d::owned_exact_taylor(factor, 100000, false, UINT64_MAX,
                                                    &options.exact_taylor_work);
    }
    require(options.exact_taylor_work && options.exact_taylor_work() == 0,
            "exact cache reader did not retain its own lifetime");
    const auto encloser = [p](acb_poly_t out, slong bits) { acb_poly_set_fmpz_poly(out, p.value(), bits); };
    d::PolynomialSource source(encloser, 3, true, options);
    options = {};
    source.research_service_meter(true);
    d::ComplexBall target, value, derivative;
    d::ComplexPolynomial expansion;
    acb_one(target.value());
    require(source.evaluate(value.value(), derivative.value(), target.value(), 128),
            "owned source lost exact cache after factor/options destruction");
    const double first_work = source.research_exact_taylor_work();
    require(first_work > 0 && source.evaluate(value.value(), derivative.value(), target.value(), 128) &&
            source.research_exact_taylor_work() == first_work,
            "prepared ball-cache reuse charged another exact operation");
    require(source.evaluate(value.value(), derivative.value(), target.value(), 256) &&
            source.research_exact_taylor_work() > first_work,
            "source precision reuse omitted the actual exact-cache enclosure");
    const double before_taylor = source.research_exact_taylor_work();
    require(source.taylor(expansion.value(), target.value(), 256) &&
            source.research_exact_taylor_work() > before_taylor && source.research_last_taylor_work() == 0,
            "explicit exact Taylor expansion was charged twice");

    d::FactorState limited;
    limited.kind = d::FactorState::Real;
    limited.integer = p;
    d::PolynomialSourceOptions refused;
    refused.exact_taylor = d::owned_exact_taylor(limited, 1, false, UINT64_MAX, &refused.exact_taylor_work);
    d::PolynomialSource fallback(encloser, 3, true, refused);
    fallback.research_service_meter(true);
    require(fallback.taylor(expansion.value(), target.value(), 128) &&
            fallback.research_exact_taylor_work() == 0 && fallback.research_last_taylor_work() > 0,
            "refused exact translation suppressed the executed Arb fallback charge");
}
#endif
static void public_secular() {
    using certpolysolve::Number;
    certpolysolve::SecularEquation p{{.residue = Number("2/3"), .pole = Number("0")},
                                 {.residue = Number("-2/3"), .pole = Number("3")}};
    known(certpolysolve::solve(p), {{1}, {2}}, 53);
    certpolysolve::SecularEquation empty{};
    auto r = certpolysolve::solve(empty);
    require(r.degree() == 0 && r.goal_met() && r.roots().size() == 0, "empty secular sum is not constant -1");
    certpolysolve::SecularEquation complex{
        {.residue = Number::complex(Number("1/3"), Number("2/7")), .pole = Number::complex(Number("2"), Number("-1"))}};
    known(certpolysolve::solve(complex), {{Rational(7) / 3, Rational(-5) / 7}}, 53);
    certpolysolve::SecularEquation repeated{{.residue = Number("1/2"), .pole = Number("0")},
                                        {.residue = Number("-1/2"), .pole = Number("2")}};
    known(certpolysolve::solve(repeated), {{1}, {1}}, 53);
    known(certpolysolve::real_roots(p), {{1}, {2}}, 53);
    // Native secular input: zero work converts nothing, one unit is the conversion, and the roots follow at 53
    // and 256 bits. (Until 18 October 2026 one session was advanced through these; one-shot solves since.)
    // A work budget of 0 selects the default, so zero work (no public call sets it since stage 3b-3) is the C
    // options' internal setting.
    const certpolysolve::SecularEquation source{{.residue = Number("1/3"), .pole = Number("0")}};
    auto zero_options = certpolysolve_test::raw_options();
    zero_options.get()->work_units = 0;
    zero_options.get()->work_units_default = false;
    const auto zero = certpolysolve_test::solve(source, zero_options.get());
    require(cps_result_termination(zero.get()) == CPS_WORK_LIMIT && zero.statistics().work_units == 0,
            "zero work performed conversion");
    certpolysolve::Options o;
    o.work_budget = 1;
    o.cell_budget = 65536;
    const auto converted = certpolysolve::solve(source, o);
    require(converted.statistics().work_units == 1 && converted.statistics().work_units == 1 &&
                !converted.goal_met(),
            "conversion was not one controller unit");
    o.work_budget = 100000;
    known(certpolysolve::solve(source, o), {{Rational(1) / 3}}, 53);
    o.bits = 256;
    known(certpolysolve::solve(source, o), {{Rational(1) / 3}}, 256);
    auto small = certpolysolve_test::raw_options();
    require(certpolysolve_test::set_exact_allowance(small.get(), 7) == CPS_OK, "exact allowance rejected");
    auto native = certpolysolve_test::solve(p, small.get());
    known(native, {{1}, {2}}, 53); // No coefficient expansion needs this budget.
    require(cps_options_set_real_only(small.get(), 1) == CPS_OK, "real-only option rejected");
    auto refused = certpolysolve_test::solve(p, small.get());
    require(cps_result_termination(refused.get()) == CPS_UNSUPPORTED_QUERY && !refused.goal_met(),
            "exact real fallback ignored its budget");
    known(certpolysolve::solve(p), {{1}, {2}}, 53); // Owned original remains usable.
}
int main() {
    certpolysolve::ThreadScope lifetime;
    try {
        evaluators();
        secular_identity();
        exact_centering();
        native_tree_candidates();
        stalled_exact_source();
        sparse_translation_budget();
        rejected_translation_cache();
#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
        exact_translation_service_work();
        exact_translation_reader_lifetime();
#endif
        candidates();
        analytic_and_overlap();
        public_secular();
        std::cout << "representation and candidate checks passed\n";
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
