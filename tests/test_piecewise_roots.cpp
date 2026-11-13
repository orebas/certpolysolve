#include "certify/attempts/newton_cover.hpp"
#include "core/engine_policy.hpp"
#include "core/engine_statistics.hpp"
#include "certify/newton_point.hpp"
#include "numeric/arithmetic/polynomial_arithmetic.hpp"
#include <algorithm>
#include <bit>
#include "exact_winding_api.hpp"
#include "numeric/evaluate/blocked_polynomial_evaluator.hpp"
#include "numeric/candidates/candidates.hpp"
#include "methods/piecewise/piecewise_roots.hpp"
#include "numeric/evaluate/structured_source.hpp"
#include <cmath>
#include <complex>
#include <numbers>
#include <cstdio>
#include <iostream>
#include <string>
#include <tuple>

namespace d = certpolysolve::detail;
using namespace reference;
static Rational endpoint(const d::Endpoint &e) {
    return reference::endpoint(certpolysolve::Dyadic{e.mantissa_text(), e.exponent_text()});
}
static void check_local_evaluation() {
    auto contains = [](const acb_t enclosure, const Complex &expected) {
        d::Rational real, imag;
        fmpq_set_str(real.value(), expected.real.str().c_str(), 10);
        fmpq_set_str(imag.value(), expected.imag.str().c_str(), 10);
        return arb_contains_fmpq(acb_realref(enclosure), real.value()) &&
               arb_contains_fmpq(acb_imagref(enclosure), imag.value());
    };
    // Exact rational value and derivative oracles exercise block boundaries,
    // a partial last block, zero/constant polynomials and uncertain inputs.
    for (const slong precision : {64, 192, 768}) for (const int degree : {-1, 0, 1, 2, 7, 37, 128}) {
        const Rational epsilon = power2(-40);
        ComplexPolynomial exact(static_cast<size_t>(degree + 1));
        d::ComplexPolynomial polynomial;
        d::ComplexBall coefficient, point, value, derivative;
        for (int j = 0; j <= degree; ++j) {
            exact[static_cast<size_t>(j)] = Complex(Rational(j % 7 - 3)/(j + 1), Rational(j % 5 - 2)/(j + 2));
            d::Number(exact[static_cast<size_t>(j)].real.str()).enclose(acb_realref(coefficient.value()), precision);
            d::Number(exact[static_cast<size_t>(j)].imag.str()).enclose(acb_imagref(coefficient.value()), precision);
            arb_add_error_2exp_si(acb_realref(coefficient.value()), -40);
            arb_add_error_2exp_si(acb_imagref(coefficient.value()), -40);
            acb_poly_set_coeff_acb(polynomial.value(), j, coefficient.value());
        }
        d::BlockedPolynomialEvaluator evaluator(polynomial, precision);
        for (const auto &center : {Complex(0), Complex(Rational(3)/8, Rational(1)/4), Complex(Rational(-5)/4)}) {
            d::Number(center.real.str()).enclose(acb_realref(point.value()), precision);
            d::Number(center.imag.str()).enclose(acb_imagref(point.value()), precision);
            arb_add_error_2exp_si(acb_realref(point.value()), -40);
            arb_add_error_2exp_si(acb_imagref(point.value()), -40);
            evaluator.evaluate(value.value(), derivative.value(), point.value());
            for (const int sign : {-1, 1}) {
                auto sample = exact;
                ComplexPolynomial differentiated(static_cast<size_t>(std::max(0, degree)));
                for (int j = 0; j <= degree; ++j) {
                    sample[static_cast<size_t>(j)] = sample[static_cast<size_t>(j)] +
                        Complex(sign * epsilon, (j % 2 ? -sign : sign) * epsilon);
                    if (j) differentiated[static_cast<size_t>(j - 1)] = sample[static_cast<size_t>(j)] * Complex(j);
                }
                const auto argument = center + Complex(sign * epsilon, -sign * epsilon);
                require(contains(value.value(), evaluate(sample, argument)) &&
                            contains(derivative.value(), evaluate(differentiated, argument)),
                        "blocked value/derivative excluded an exact uncertain-input sample");
            }
        }
    }
}
static void check() {
    check_local_evaluation();
    cps_options options;
    options.bits = 53;
    options.max_work_bits = 1024;
    d::RootModel model;
    acb_one(model.center.value());
    arb_set_d(model.radius.value(), 0.5);
    d::ComplexBall coefficient, point;
    acb_set_d(coefficient.value(), -0.5);
    acb_poly_set_coeff_acb(model.polynomial.value(), 1, coefficient.value());
    acb_set_d(coefficient.value(), 0.25);
    acb_poly_set_coeff_acb(model.polynomial.value(), 2, coefficient.value());
    d::Root root;
    require(d::certify_root_model(root, model, point.value(), options, 192),
            "linear Rouche certificate rejected exact root");
    ComplexPolynomial exact{{2}, {-3}, {1}};
    require(rectangle_count(exact, endpoint(root.endpoints[0]), endpoint(root.endpoints[1]),
                            endpoint(root.endpoints[2]), endpoint(root.endpoints[3])) == 1,
            "independent winding rejected published rectangle");
    arb_one(model.error.value());
    d::ChartFailure failure;
    require(!d::certify_root_model(root, model, point.value(), options, 192, &failure) &&
                failure == d::ChartFailure::ModelRemainder,
            "model remainder was ignored or misclassified");
    arb_zero(model.error.value());
    cps_options guard_options = options;
    guard_options.absolute_tolerance = d::Number("1/8");
    require(d::certify_root_model(root, model, point.value(), guard_options, 192) &&
                rectangle_count(exact, endpoint(root.endpoints[0]), endpoint(root.endpoints[1]),
                                endpoint(root.endpoints[2]), endpoint(root.endpoints[3])) == 1,
            "larger retirement rectangle failed independent winding");
    acb_set_si(point.value(), 2);
    require(!d::certify_root_model(root, model, point.value(), options, 192), "chart boundary was ignored");

    // Exact points cannot rescue a box wider than the pair separation. Check
    // both contracts independently on the exported dyadic rectangles.
    for (const int separation : {40, 100}) {
        const Rational epsilon = Rational(1) / Rational(Integer(1) << separation);
        const auto pair = from_roots({Complex(Rational(1) - epsilon), Complex(Rational(1) + epsilon)});
        const auto automatic_input = reference::input(pair);
        d::RootModel close;
        acb_one(close.center.value());
        arb_set_d(close.radius.value(), 0.5);
        d::Number(Rational(-epsilon * epsilon).str()).enclose(acb_realref(coefficient.value()), 512);
        arb_zero(acb_imagref(coefficient.value()));
        acb_poly_set_coeff_acb(close.polynomial.value(), 0, coefficient.value());
        acb_set_d(coefficient.value(), 0.25);
        acb_poly_set_coeff_acb(close.polynomial.value(), 2, coefficient.value());
        for (const auto goal : {CPS_ROOT_ACCURACY, CPS_ROOT_APPROXIMATION}) {
            auto query = options;
            query.goal = goal;
            d::PiecewiseRootStatistics automatic_stats;
            const auto automatic = d::solve_piecewise_roots(certpolysolve::detail::dense(*automatic_input.get()), query, automatic_stats);
            if (!automatic.goal_met || !d::disjoint(automatic.roots))
                throw std::runtime_error("automatic close-pair candidates k=" + std::to_string(separation) +
                    " goal=" + std::to_string(goal) + " count=" + std::to_string(automatic.roots.size()) +
                    " proposals=" + std::to_string(automatic_stats.proposals) +
                    " candidate=" + std::to_string(automatic_stats.candidate_limited) +
                    " rounding=" + std::to_string(automatic_stats.rounding_limited) +
                    " rejected=" + std::to_string(automatic_stats.seed_rejected));
            d::Number(Rational(2 * epsilon).str()).enclose(acb_realref(point.value()), 512);
            arb_zero(acb_imagref(point.value()));
            d::PiecewiseRootStatistics observed;
            require(d::certify_root_model(root, close, point.value(), query, 512, nullptr, &observed.certificate),
                    "separation-aware certificate failed exact close pair");
            require(rectangle_count(pair, endpoint(root.endpoints[0]), endpoint(root.endpoints[1]),
                                    endpoint(root.endpoints[2]), endpoint(root.endpoints[3])) == root.count,
                    "adaptive close-pair rectangle failed independent count");
            if (goal == CPS_ROOT_ACCURACY) {
                require(root.count == 1 && root.multiplicity == 1, "isolation accepted a cluster");
                if (separation == 100) require(observed.certificate.radius_shrinks > 0, "wide box was not shrunk");
            } else if (separation == 100)
                require(root.count == 2 && root.multiplicity == 0, "counted pair lost its contract");
        }
    }
    {
        // A translated complex triplet exercises m=3, not only quadratic
        // clusters. All original roots and exported counts are exact oracles.
        const Rational epsilon = power2(-100);
        const Complex center(-2, 3);
        const auto exact_triplet = from_roots({center - Complex(epsilon / 2), center,
                                               center + Complex(epsilon / 2)});
        d::RootModel triplet;
        acb_set_si_si(triplet.center.value(), -2, 3);
        arb_set_d(triplet.radius.value(), 0.5);
        d::Number(Rational(-epsilon * epsilon).str()).enclose(acb_realref(coefficient.value()), 512);
        arb_zero(acb_imagref(coefficient.value()));
        acb_poly_set_coeff_acb(triplet.polynomial.value(), 1, coefficient.value());
        acb_one(coefficient.value());
        acb_poly_set_coeff_acb(triplet.polynomial.value(), 3, coefficient.value());
        acb_zero(point.value());
        for (const auto goal : {CPS_ROOT_ACCURACY, CPS_ROOT_APPROXIMATION}) {
            auto query = options;
            query.goal = goal;
            require(d::certify_root_model(root, triplet, point.value(), query, 512),
                    "translated triplet certificate failed");
            require(root.count == (goal == CPS_ROOT_ACCURACY ? 1 : 3) &&
                        rectangle_count(exact_triplet, endpoint(root.endpoints[0]), endpoint(root.endpoints[1]),
                                        endpoint(root.endpoints[2]), endpoint(root.endpoints[3])) == root.count,
                    "translated triplet count failed independent winding");
        }
        arb_one(triplet.error.value());
        require(!d::certify_root_model(root, triplet, point.value(), options, 512),
                "triplet ignored an inadequate model remainder");
    }

    // Independently evaluate exact rational samples, including omitted bands
    // and truncated Taylor models. The FFT phase and envelope are both tested.
    d::ComplexPolynomial polynomial;
    ComplexPolynomial rational(257);
    for (slong j = 0; j <= 256; ++j) {
        rational[static_cast<size_t>(j)] = Complex(j % 5 - 2, j % 3 - 1);
        if (j == 0 || j == 256)
            rational[static_cast<size_t>(j)] = Complex(1);
        const auto &c = rational[static_cast<size_t>(j)];
        d::Number(c.real.str()).enclose(acb_realref(coefficient.value()), 192);
        d::Number(c.imag.str()).enclose(acb_imagref(coefficient.value()), 192);
        acb_poly_set_coeff_acb(polynomial.value(), j, coefficient.value());
    }
    d::PiecewiseRootStatistics stats;
    size_t seen = 0, samples = 0;
    bool omitted = false, truncated = false;
    require(d::piecewise_root_models(
                polynomial.value(), 8, 192, options, stats,
                [&](const d::RootModel &m) {
                    ++seen;
                    const bool needed = (acb_poly_degree(m.polynomial.value()) >= 96 && !truncated) ||
                                        (m.lower_index > 0 && !omitted);
                    if (!needed && (seen % 211 || samples >= 6))
                        return true;
                    ++samples;
                    omitted |= m.lower_index > 0;
                    truncated |= acb_poly_degree(m.polynomial.value()) >= 96;
                    d::ComplexBall z, w, value, power;
                    acb_set_d_d(w.value(), 0.25, 0.125);
                    acb_mul_arb(z.value(), w.value(), m.radius.value(), 192);
                    acb_add(z.value(), z.value(), m.center.value(), 192);
                    acb_get_mid(z.value(), z.value());
                    auto rational_component = [](const arb_t x) {
                        d::Integer a, b;
                        arf_get_fmpz_2exp(a.value(), b.value(), arb_midref(x));
                        return reference::endpoint(
                            certpolysolve::Dyadic{d::integer_string(a.value()), d::integer_string(b.value())});
                    };
                    const Complex q(rational_component(acb_realref(z.value())),
                                    rational_component(acb_imagref(z.value())));
                    acb_sub(w.value(), z.value(), m.center.value(), 192);
                    acb_div_arb(w.value(), w.value(), m.radius.value(), 192);
                    acb_poly_evaluate(value.value(), m.polynomial.value(), w.value(), 192);
                    arb_add_error(acb_realref(value.value()), m.error.value());
                    arb_add_error(acb_imagref(value.value()), m.error.value());
                    acb_pow_ui(power.value(), z.value(), static_cast<ulong>(m.lower_index), 192);
                    acb_mul(value.value(), value.value(), power.value(), 192);
                    acb_mul_2exp_si(value.value(), value.value(), m.scale);
                    const auto expected = evaluate(rational, q);
                    d::Rational re, im;
                    fmpq_set_str(re.value(), expected.real.str().c_str(), 10);
                    fmpq_set_str(im.value(), expected.imag.str().c_str(), 10);
                    require(arb_contains_fmpq(acb_realref(value.value()), re.value()) &&
                                arb_contains_fmpq(acb_imagref(value.value()), im.value()),
                            "piecewise tail or FFT excluded exact sample");
                    return true;
                }),
            "no piecewise models");
    require(samples > 0 && stats.fft_transforms > 0 && truncated && omitted,
            "bound test missed truncated or omitted terms");

    d::PiecewiseRootStatistics batched;
    size_t batch_models = 0;
    require(d::piecewise_root_models(polynomial.value(), 8, 192, options, batched,
                [&](const d::RootModel &) { ++batch_models; return true; }, {}, {size_t(1) << 20}) &&
                batch_models == seen && batched.batches > batched.rings,
            "bounded batching refused a feasible ring or lost models");

    const auto known =
        from_roots({Complex(1), Complex(2), Complex(0, 3), Complex(-1, 2), Complex(0), Complex(0)});
    const auto input = reference::input(known);
    {
        // Reusing a no-update observation must preserve every exported box,
        // not merely produce a second plausible approximation to this input.
        // A fresh proposal is polished once (proposals() leaves polishing to
        // the caller), so an uninterrupted solve has no repeated observation:
        // the switch changes neither the work nor a box here. It still spares
        // a resumed state the points already at their model's floor.
        certpolysolve::Polynomial circle;
        circle.set(0, certpolysolve::Number("-1"));
        circle.set(257, certpolysolve::Number("1"));
        d::PiecewiseRootConfig before, after;
        before.reuse_guards = after.reuse_guards = true;
        after.reuse_polish = true;
        d::PiecewiseRootStatistics old_work, new_work;
        const auto old_result = d::solve_piecewise_roots(certpolysolve::detail::dense(*circle.get()), options, old_work, before);
        const auto new_result = d::solve_piecewise_roots(certpolysolve::detail::dense(*circle.get()), options, new_work, after);
        require(old_result.goal_met && new_result.goal_met && new_work.polish_reused == 0 &&
                    new_work.polish_evaluations == old_work.polish_evaluations &&
                    old_result.roots.size() == new_result.roots.size(),
                "an uninterrupted solve polished a proposal twice, or polishing reuse lost completion");
        for (size_t j = 0; j < old_result.roots.size(); ++j)
            for (size_t k = 0; k < 4; ++k)
                require(old_result.roots[j].endpoints[k].mantissa_text() == new_result.roots[j].endpoints[k].mantissa_text() &&
                            old_result.roots[j].endpoints[k].exponent_text() == new_result.roots[j].endpoints[k].exponent_text(),
                        "polishing reuse changed a published rectangle");
    }
    {
        d::PolynomialSource source([&](acb_poly_t p, slong precision) {
            acb_poly_zero(p);
            for (const auto &[j, c] : input.get()->terms) {
                c.enclose(coefficient.value(), precision);
                acb_poly_set_coeff_acb(p, static_cast<slong>(j), coefficient.value());
            }
        }, static_cast<slong>(known.size()-1), false);
        for (const auto &z : {Complex(1), Complex(2), Complex(0, 3), Complex(-1, 2)}) {
            d::Number(z.real.str()).enclose(acb_realref(point.value()), 128);
            d::Number(z.imag.str()).enclose(acb_imagref(point.value()), 128);
            require(d::certify_newton_point(root, source, point.value(), options, 128) && root.count == 1 &&
                        rectangle_count(known, endpoint(root.endpoints[0]), endpoint(root.endpoints[1]),
                                        endpoint(root.endpoints[2]), endpoint(root.endpoints[3])) == 1,
                    "Newton derivative-variation certificate failed independent winding");
        }
        acb_set_si_si(point.value(), 12, 3);
        require(!d::certify_newton_point(root, source, point.value(), options, 128),
                "Newton certificate accepted a distant point");
        acb_zero(point.value());
        require(!d::certify_newton_point(root, source, point.value(), options, 128),
                "Newton certificate fabricated a relative certificate at zero");
    }
    for (bool optimized : {false, true}) {
        // The coefficient remainder must use the monic source's units, even
        // when the input is scaled by a tiny nonreal constant. Check all
        // returned counts using the independent exact winding implementation.
        const auto exact=from_roots({Complex(1),Complex(2),Complex(0,3),Complex(-1,2)});
        const Complex factor=optimized?Complex(power2(-100),power2(-102)):Complex(1);
        auto scaled=exact;for(auto &a:scaled)a=a*factor;
        const auto p=reference::input(scaled);
        d::PolynomialSource source([&](acb_poly_t out,slong bits) {
            d::ComplexBall a;
            acb_poly_zero(out);
            for(const auto &[j,c]:p.get()->terms) {
                c.enclose(a.value(),bits);
                acb_poly_set_coeff_acb(out,static_cast<slong>(j),a.value());
            }
        },4,false);
        cps_options query=options;query.goal=CPS_ROOT_APPROXIMATION;
        certpolysolve::detail::CandidateStatistics work{};std::vector<d::Root> cover;
        require(d::machine_newton_cover(cover,source,nullptr,query,192,true,work,true) && cover.size()==4,
                "coefficient Newton cover failed under nonreal scaling");
        for(const auto &r:cover)require(rectangle_count(scaled,endpoint(r.endpoints[0]),endpoint(r.endpoints[1]),
                endpoint(r.endpoints[2]),endpoint(r.endpoints[3]))==1,"coefficient Newton cover failed independent winding");
        d::PolynomialSource inconsistent([&](acb_poly_t out,slong bits) {
            const auto enclosed=source.coefficient_enclosure(bits);acb_poly_set(out,enclosed->polynomial.value());
        },3,false);
        require(!d::machine_newton_cover(cover,inconsistent,nullptr,query,192,true,work,true) && cover.size()==4,
                "inconsistent source degree reached candidate storage or replaced output");
        query.max_work_bits=128;
        require(!d::machine_newton_cover(cover,source,nullptr,query,192,true,work,true) && cover.size()==4,
                "coefficient certificate ignored precision ceiling");
        query.max_work_bits=4096;query.cancel=[](void *)->int{return 1;};
        require(!d::machine_newton_cover(cover,source,nullptr,query,192,true,work,true) && cover.size()==4,
                "cancelled coefficient cover replaced output");
    }
    {
        // A dense, well conditioned polynomial at complex points: 400 exact
        // rational points of the unit circle, ((1-t^2)/(1+t^2), 2t/(1+t^2))
        // and their mirror images, with t on a grid of 1/1024 nearest to the
        // roots of z^400 = -1. Ball Horner at degree 400 re-boxes its error
        // rectangle 400 times and loses up to half a bit a term, so its
        // 85-bit value has no correct bit at most of these points. The
        // attempt must evaluate by disk blocks and certify every root at its
        // first precision from one evaluation a point: a binary64 point is a
        // Newton step from its rectangle, and the certificate at the step's
        // end comes from the evaluation at the point (the linear model with
        // its remainder on a disk that holds the step), where a second
        // evaluation at the stepped point once gave it.
        std::vector<Complex> planted;
        constexpr int quarter = 100, grid = 1024;
        for (int k = 0; k < quarter; ++k) {
            const double angle = std::numbers::pi / 2 * (k + 0.5) / quarter;
            const Rational t = Rational(std::llround(std::tan(angle / 2) * grid)) / grid;
            const Rational x = (1 - t * t) / (1 + t * t), y = 2 * t / (1 + t * t);
            for (const auto &z : {Complex(x, y), Complex(-x, y), Complex(x, -y), Complex(-x, -y)})
                planted.push_back(z);
        }
        const auto p = reference::input(from_roots(planted));
        d::PolynomialSource source([&](acb_poly_t out, slong bits) {
            d::ComplexBall a;
            acb_poly_zero(out);
            for (const auto &[j, c] : p.get()->terms) {
                c.enclose(a.value(), bits);
                acb_poly_set_coeff_acb(out, static_cast<slong>(j), a.value());
            }
        }, static_cast<slong>(planted.size()), false);
        cps_options query = options;
        query.goal = CPS_ROOT_APPROXIMATION;
        d::EnginePolicy observed;
        observed.study_statistics = std::make_shared<d::StudyStatistics>();
        observed.compensated_certificates = false; // the ball certificate's own evaluations are under test
        certpolysolve::detail::CandidateStatistics work{};
        std::vector<d::Root> cover;
        require(d::machine_newton_cover(cover, source, nullptr, query, 85, true, work, true, &observed) &&
                    cover.size() == planted.size(),
                "the Newton-certificate attempt refused a dense polynomial with simple, separated complex roots");
        const auto &counts = *observed.study_statistics;
        require(counts.machine_cover_precision_retries == 0 && counts.machine_cover_peak_bits == 85,
                "the attempt's evaluation at complex points lost its precision (error-box wrapping)");
        // One evaluation a point, with a margin of an eighth of the points
        // for a step and a second one (a step and its certificate took two
        // at every point).
        const size_t margin = planted.size() / 8;
        require(counts.machine_cover_block_evaluations >= planted.size() &&
                    counts.machine_cover_block_evaluations <= planted.size() + margin,
                "the attempt did not evaluate by blocks, or evaluated its points a second time for their certificates");
        require(2 * counts.machine_cover_stepped_certificates > planted.size(),
                "the attempt did not certify its points at their Newton points from one evaluation");
        const auto holds_planted = [&](const std::vector<d::Root> &regions) {
            std::vector<uint8_t> found(planted.size(), 0);
            for (const auto &r : regions) {
                const Rational real_lower = endpoint(r.endpoints[0]), real_upper = endpoint(r.endpoints[1]);
                const Rational imag_lower = endpoint(r.endpoints[2]), imag_upper = endpoint(r.endpoints[3]);
                size_t inside = 0;
                for (size_t j = 0; j < planted.size(); ++j)
                    if (real_lower <= planted[j].real && planted[j].real <= real_upper &&
                        imag_lower <= planted[j].imag && planted[j].imag <= imag_upper) {
                        ++inside;
                        ++found[j];
                    }
                require(r.count == 1 && inside == 1,
                        "a rectangle of the attempt does not hold exactly its planted root");
            }
            require(std::all_of(found.begin(), found.end(), [](uint8_t f) { return f == 1; }),
                    "the attempt's rectangles do not cover every planted root once");
        };
        holds_planted(cover);
        require(counts.machine_cover_rung_evaluations == 0,
                "the Newton ladder ran below an 85-bit certificate precision (it starts at 128 bits)");
        // The same roots to 1,000 bits: the certificate's Newton steps from the
        // binary64 points climb the ladder (128, 256, 512, 1,024 bits) before
        // the certificate precision, where one evaluation remains: the point
        // the ladder leaves is 2^-749 of its modulus from the root, and the
        // certificate at its Newton point needs a disk of that radius around
        // a rectangle 2^255 times smaller (the disk carries a relative slack
        // for the certificate's own bounds on moduli, or it is refused).
        d::EnginePolicy high;
        high.study_statistics = std::make_shared<d::StudyStatistics>();
        cps_options thousand = query;
        thousand.bits = 1000;
        thousand.max_work_bits = 4096;
        std::vector<d::Root> fine;
        require(d::machine_newton_cover(fine, source, nullptr, thousand, 1032, true, work, true, &high) &&
                    fine.size() == planted.size(),
                "the attempt refused the dense polynomial at a 1,000-bit goal");
        const auto &ladder = *high.study_statistics;
        // Every point climbs the ladder, and at the certificate precision a
        // point needs one evaluation, with the margin above (from binary64
        // points, six steps there took six, and a step and its certificate
        // two).
        require(ladder.machine_cover_rung_evaluations >= planted.size(),
                "the attempt's Newton steps toward 1,000 bits did not use the ladder below the certificate precision");
        require(ladder.machine_cover_block_evaluations <= planted.size() + margin &&
                    2 * ladder.machine_cover_stepped_certificates > planted.size(),
                "the attempt evaluated its points more than once at the full certificate precision");
        for (const auto &r : fine)
            require(r.count == 1 && d::meets_goal(r, thousand, 1032), "a 1,000-bit rectangle misses its goal");
        holds_planted(fine);
    }
    {
        // An ill-conditioned dense polynomial: 100 rational points on each of
        // the circles of radius 1 and 5/4. At its worst roots |c| P''(|c|) /
        // |p'(c)|, with P'' the coefficient majorant of p'', is about 2^46.
        // The linear model's remainder passes only when it is bounded on the
        // certified rectangle's own disk (two half-widths: up to 2^55); on a
        // disk of 512 half-widths, as the attempt once used, it is 2^16 times
        // larger and the limit is 2^39. Binary64 cannot resolve these roots
        // (their evaluation loses up to 32 bits), so the attempt starts from
        // the planted points rounded to binary64; it needs a precision retry
        // and a Newton step.
        std::vector<Complex> planted;
        for (const auto &radius : {Rational(1), Rational(5) / 4})
            for (int k = 1; k <= 25; ++k) {
                const Rational t = Rational(k) / 26, scale = radius / (1 + t * t);
                const Rational x = (1 - t * t) * scale, y = 2 * t * scale;
                for (const auto &z : {Complex(x, y), Complex(-x, y), Complex(x, -y), Complex(-x, -y)})
                    planted.push_back(z);
            }
        const auto p = reference::input(from_roots(planted));
        d::PolynomialSource source([&](acb_poly_t out, slong bits) {
            d::ComplexBall a;
            acb_poly_zero(out);
            for (const auto &[j, c] : p.get()->terms) {
                c.enclose(a.value(), bits);
                acb_poly_set_coeff_acb(out, static_cast<slong>(j), a.value());
            }
        }, static_cast<slong>(planted.size()), false);
        d::RootVector seeds(static_cast<slong>(planted.size()));
        for (size_t j = 0; j < planted.size(); ++j) {
            d::Number(planted[j].real.str()).enclose(acb_realref(seeds[static_cast<slong>(j)]), 53);
            d::Number(planted[j].imag.str()).enclose(acb_imagref(seeds[static_cast<slong>(j)]), 53);
            acb_get_mid(seeds[static_cast<slong>(j)], seeds[static_cast<slong>(j)]);
        }
        cps_options query = options;
        query.goal = CPS_ROOT_APPROXIMATION;
        d::EnginePolicy observed;
        observed.study_statistics = std::make_shared<d::StudyStatistics>();
        observed.compensated_certificates = false; // the ball certificate's retries and steps are under test
        certpolysolve::detail::CandidateStatistics work{};
        std::vector<d::Root> cover;
        require(d::machine_newton_cover(cover, source, seeds.data(), query, 85, true, work, true, &observed) &&
                    cover.size() == planted.size(),
                "the Newton-certificate attempt refused simple roots whose remainder majorant is large");
        require(observed.study_statistics->machine_cover_precision_retries != 0 &&
                    observed.study_statistics->machine_cover_newton_steps != 0,
                "the ill-conditioned fixture needed no precision retry or no Newton step: it tests nothing");
        std::vector<uint8_t> found(planted.size(), 0);
        for (const auto &r : cover) {
            const Rational real_lower = endpoint(r.endpoints[0]), real_upper = endpoint(r.endpoints[1]);
            const Rational imag_lower = endpoint(r.endpoints[2]), imag_upper = endpoint(r.endpoints[3]);
            size_t inside = 0;
            for (size_t j = 0; j < planted.size(); ++j)
                if (real_lower <= planted[j].real && planted[j].real <= real_upper &&
                    imag_lower <= planted[j].imag && planted[j].imag <= imag_upper) {
                    ++inside;
                    ++found[j];
                }
            require(r.count == 1 && inside == 1, "a rectangle of the ill-conditioned attempt does not hold its root");
        }
        require(std::all_of(found.begin(), found.end(), [](uint8_t f) { return f == 1; }),
                "the ill-conditioned attempt's rectangles do not cover every planted root once");
    }
    {
        // The attempt goes on past a refusal. Shared by the fixtures below:
        // the source of a polynomial with planted roots, and the check that
        // every rectangle of a cover holds exactly its planted root.
        const auto source_of = [](const certpolysolve::Polynomial &p, slong fixed_bits = 0) {
            return d::PolynomialSource([&p, fixed_bits](acb_poly_t out, slong bits) {
                d::ComplexBall a;
                acb_poly_zero(out);
                for (const auto &[j, c] : p.get()->terms) {
                    c.enclose(a.value(), fixed_bits ? fixed_bits : bits);
                    acb_poly_set_coeff_acb(out, static_cast<slong>(j), a.value());
                }
            }, static_cast<slong>(cps_polynomial_degree(p.get())), false);
        };
        const auto holds_planted = [](const std::vector<d::Root> &cover, const std::vector<Complex> &planted) {
            std::vector<uint8_t> found(planted.size(), 0);
            for (const auto &r : cover) {
                const Rational real_lower = endpoint(r.endpoints[0]), real_upper = endpoint(r.endpoints[1]);
                const Rational imag_lower = endpoint(r.endpoints[2]), imag_upper = endpoint(r.endpoints[3]);
                size_t inside = 0;
                for (size_t j = 0; j < planted.size(); ++j)
                    if (real_lower <= planted[j].real && planted[j].real <= real_upper &&
                        imag_lower <= planted[j].imag && planted[j].imag <= imag_upper) {
                        ++inside;
                        ++found[j];
                    }
                if (r.count != 1 || inside != 1) return false;
            }
            return cover.size() == planted.size() &&
                   std::all_of(found.begin(), found.end(), [](uint8_t f) { return f == 1; });
        };
        const auto seeds_of = [](const std::vector<Complex> &planted) {
            auto seeds = std::make_unique<d::RootVector>(static_cast<slong>(planted.size()));
            for (size_t j = 0; j < planted.size(); ++j) {
                const auto z = (*seeds)[static_cast<slong>(j)];
                d::Number(planted[j].real.str()).enclose(acb_realref(z), 53);
                d::Number(planted[j].imag.str()).enclose(acb_imagref(z), 53);
                acb_get_mid(z, z);
            }
            return seeds;
        };
        std::vector<Complex> circle;
        constexpr int quarter = 100, grid = 1024;
        for (int k = 0; k < quarter; ++k) {
            const double angle = std::numbers::pi / 2 * (k + 0.5) / quarter;
            const Rational t = Rational(std::llround(std::tan(angle / 2) * grid)) / grid;
            const Rational x = (1 - t * t) / (1 + t * t), y = 2 * t / (1 + t * t);
            for (const auto &z : {Complex(x, y), Complex(-x, y), Complex(x, -y), Complex(-x, -y)})
                circle.push_back(z);
        }
        cps_options query = options;
        query.goal = CPS_ROOT_APPROXIMATION;
        {
            // Three pairs of roots 2^-23 apart beside 400 well separated
            // ones. Binary64 cannot tell a pair's roots apart (its points
            // stop about 2^-15 from the pair, where Newton's iteration
            // contracts by one half), so the first pass refuses six points
            // and certifies the rest; multiprecision packets on the refused
            // points after the pass (one round on all six, then on any left),
            // with the others held fixed, separate them, and they certify
            // after a precision raise that the rounds keep.
            auto planted = circle;
            const Rational gap = power2(-24);
            for (const auto &z : {Complex(Rational(3) / 2 + gap), Complex(Rational(3) / 2 - gap),
                                  Complex(0, 2 + gap), Complex(0, 2 - gap), Complex(0, -2 - gap),
                                  Complex(0, -2 + gap)})
                planted.push_back(z);
            const auto p = reference::input(from_roots(planted));
            auto source = source_of(p);
            d::EnginePolicy observed;
            observed.study_statistics = std::make_shared<d::StudyStatistics>();
            // The rounds and the raise are checked without the double-double repair (proposal 1), which proves the
            // pairs in the first pass (checked after).
            observed.double_double_repair = false;
            certpolysolve::detail::CandidateStatistics work{};
            std::vector<d::Root> cover;
            slong certified_at = 0;
            require(d::machine_newton_cover(cover, source, nullptr, query, 85, true, work, true, &observed, nullptr,
                                            nullptr, &certified_at) &&
                        holds_planted(cover, planted),
                    "the attempt did not finish a polynomial with three close pairs among separated roots");
            require(observed.study_statistics->machine_cover_rounds >= 1 &&
                        observed.study_statistics->machine_cover_rounds <= 6 &&
                        observed.study_statistics->machine_cover_precision_retries != 0 && certified_at > 85,
                    "the close pairs needed no round or no precision raise, or the raise was not reported");
            // With the repair: one double-double repair of the six refused points, then the compensated tier proves
            // them, with no multiprecision round and no precision raise.
            d::EnginePolicy repaired;
            repaired.study_statistics = std::make_shared<d::StudyStatistics>();
            std::vector<d::Root> repaired_cover;
            require(d::machine_newton_cover(repaired_cover, source, nullptr, query, 85, true, work, true, &repaired) &&
                        holds_planted(repaired_cover, planted),
                    "the repaired attempt did not finish the three close pairs");
            require(repaired.study_statistics->machine_cover_double_double_repairs == 1 &&
                        repaired.study_statistics->machine_cover_rounds == 0,
                    "the close pairs were not finished by one double-double repair and no round");
        }
        {
            // The cluster start (research/d8-clusters/README.md): four roots 2^-45 from 7/10, (z - 7/10)^4 = 2^-180,
            // among four separated ones. Binary64 leaves the cluster's points about 2^-13 from it; the tiers refuse
            // them and the double-double repair cannot reach them (4 x 45 bits beyond the leading ones). Their
            // Newton disks all meet; one packet from a circle about the centre, at the precision the cluster needs,
            // resolves it, and the attempt finishes. Without the start it does not.
            const Rational gap = power2(-45), centre = Rational(7) / 10;
            const std::vector<Complex> planted = {
                Complex(centre + gap), Complex(centre - gap), Complex(centre, gap), Complex(centre, -gap),
                Complex(Rational(-3) / 2), Complex(Rational(6) / 5), Complex(Rational(-1) / 3, Rational(4) / 5),
                Complex(Rational(-1) / 3, Rational(-4) / 5)};
            const auto p = reference::input(from_roots(planted));
            auto source = source_of(p);
            certpolysolve::detail::CandidateStatistics work{};
            d::EnginePolicy started;
            started.study_statistics = std::make_shared<d::StudyStatistics>();
            std::vector<d::Root> cover;
            require(d::machine_newton_cover(cover, source, nullptr, query, 85, true, work, true, &started) &&
                        holds_planted(cover, planted),
                    "the attempt did not finish a 4-fold cluster 2^-45 across with the cluster start");
            require(started.study_statistics->machine_cover_cluster_starts == 1,
                    "the 4-fold cluster was not certified by one cluster start");
            d::EnginePolicy plain;
            plain.study_statistics = std::make_shared<d::StudyStatistics>();
            plain.cluster_start = false;
            std::vector<d::Root> plain_cover;
            require(!d::machine_newton_cover(plain_cover, source, nullptr, query, 85, true, work, true, &plain),
                    "the attempt finished the 4-fold cluster without the cluster start: the fixture no longer tests it");
        }
        {
            // Inclusion disks: 100 rational points on each of four circles of
            // radius 1 to 5/2. On the outer circles |c| P''(|c|) / |p'(c)|
            // reaches 2^65, past what the linear model's remainder can carry
            // (2^55), so those roots are certified by the disk of radius
            // n |p/p'| and the closing count. From planted seeds (binary64
            // cannot resolve the inner roots).
            std::vector<Complex> planted;
            for (const auto &radius : {Rational(1), Rational(3) / 2, Rational(2), Rational(5) / 2})
                for (int k = 1; k <= 25; ++k) {
                    const Rational t = Rational(k) / 26, scale = radius / (1 + t * t);
                    const Rational x = (1 - t * t) * scale, y = 2 * t * scale;
                    for (const auto &z : {Complex(x, y), Complex(-x, y), Complex(x, -y), Complex(-x, -y)})
                        planted.push_back(z);
                }
            const auto p = reference::input(from_roots(planted));
            auto source = source_of(p);
            const auto seeds = seeds_of(planted);
            d::EnginePolicy observed;
            observed.study_statistics = std::make_shared<d::StudyStatistics>();
            certpolysolve::detail::CandidateStatistics work{};
            std::vector<d::Root> cover;
            require(d::machine_newton_cover(cover, source, seeds->data(), query, 85, true, work, true, &observed) &&
                        holds_planted(cover, planted),
                    "the attempt did not certify roots whose linear model is refused by its remainder");
            require(observed.study_statistics->machine_cover_disks != 0,
                    "the four-circle fixture needed no inclusion disk: it tests nothing");
            // The same polynomial through a source that returns its 85-bit
            // enclosure at every precision: rounding that precision does not
            // reduce ends the attempt at the first point that needs more.
            auto snapshot = source_of(p, 85);
            d::EnginePolicy limited;
            limited.study_statistics = std::make_shared<d::StudyStatistics>();
            std::vector<d::Root> kept = cover;
            require(!d::machine_newton_cover(kept, snapshot, seeds->data(), query, 85, true, work, true, &limited) &&
                        kept.size() == cover.size(),
                    "the attempt certified from a source whose rounding it cannot reduce, or replaced the output");
            require(limited.study_statistics->machine_cover_rounds == 0 &&
                        limited.study_statistics->machine_cover_block_evaluations < planted.size(),
                    "a refusal for irreducible rounding did not end the attempt at once");
        }
        {
            // Outside binary64's range: 48 rational points on each of the
            // circles of radius 2^-300 and 2^300. No scaling of the variable
            // brings the coefficients (spread over 2^28800) into binary64, so
            // the attempt's packet runs with a separate exponent, and the
            // certificates follow as for any other polynomial.
            std::vector<Complex> planted;
            constexpr int twelfth = 12;
            for (const auto &radius : {power2(-300), power2(300)})
                for (int k = 0; k < twelfth; ++k) {
                    const double angle = std::numbers::pi / 2 * (k + 0.5) / twelfth;
                    const Rational t = Rational(std::llround(std::tan(angle / 2) * grid)) / grid;
                    const Rational x = radius * (1 - t * t) / (1 + t * t), y = radius * 2 * t / (1 + t * t);
                    for (const auto &z : {Complex(x, y), Complex(-x, y), Complex(x, -y), Complex(-x, -y)})
                        planted.push_back(z);
                }
            const auto p = reference::input(from_roots(planted));
            auto source = source_of(p);
            d::EnginePolicy observed;
            observed.study_statistics = std::make_shared<d::StudyStatistics>();
            certpolysolve::detail::CandidateStatistics work{};
            std::vector<d::Root> cover;
            require(d::machine_newton_cover(cover, source, nullptr, query, 85, true, work, true, &observed) &&
                        holds_planted(cover, planted),
                    "the attempt did not solve a polynomial whose coefficients are outside binary64's range");
            require(observed.study_statistics->scaled_packets == 1,
                    "the out-of-range fixture did not need the scaled packet: it tests nothing");
        }
        {
            // A cluster the attempt cannot finish: sixteen roots 2^-70 apart
            // beside 400 separated ones. The first pass certifies the 400 and
            // refuses the sixteen; one round moves none of them to a
            // certificate, and the attempt ends there, without a precision
            // raise for points that are refused for their accuracy.
            auto planted = circle;
            for (int k = 0; k < 16; ++k)
                planted.push_back(Complex(Rational(3) / 2 + k * power2(-70)));
            const auto p = reference::input(from_roots(planted));
            auto source = source_of(p);
            d::EnginePolicy observed;
            observed.study_statistics = std::make_shared<d::StudyStatistics>();
            certpolysolve::detail::CandidateStatistics work{};
            std::vector<d::Root> cover;
            require(!d::machine_newton_cover(cover, source, nullptr, query, 85, true, work, true, &observed) &&
                        cover.empty(),
                    "the attempt claimed a cluster of sixteen roots 2^-70 apart, or left a partial output");
            require(observed.study_statistics->machine_cover_rounds <= 1 &&
                        observed.study_statistics->machine_cover_block_evaluations <= 4 * planted.size(),
                    "a cluster the attempt cannot finish cost more than its allowance");
        }
    }
    for (bool optimized : {false, true}) {
        d::PiecewiseRootConfig experiment;
        experiment.reuse_guards = true;
        experiment.recover_discovery = optimized;
        experiment.reuse_polish = optimized;
        experiment.conjugate_symmetry = optimized;
        experiment.blocked_polish = optimized;
        d::PiecewiseRootStatistics measured;
        const auto checked = d::solve_piecewise_roots(certpolysolve::detail::dense(*input.get()), options, measured, experiment);
        require(checked.goal_met && d::disjoint(checked.roots), "local-work ablation lost known complex roots");
        uint64_t total = 0;
        for (const auto &r : checked.roots) {
            total += r.count;
            require(rectangle_count(known, endpoint(r.endpoints[0]), endpoint(r.endpoints[1]),
                                    endpoint(r.endpoints[2]), endpoint(r.endpoints[3])) == r.count,
                    "local-work ablation failed independent winding");
        }
        require(total == known.size() - 1, "local-work ablation lost zero multiplicity");
        require(measured.reflected_roots == 0, "complex coefficients were treated as conjugate symmetric");
    }
    {
        const auto symmetric = from_roots({Complex(0), Complex(0), Complex(-2), Complex(1),
            Complex(0, 1), Complex(0, -1), Complex(1, 2), Complex(1, -2)});
        const auto p = reference::input(symmetric);
        d::PiecewiseRootConfig experiment;
        experiment.conjugate_symmetry = true;
        experiment.reuse_polish = experiment.recover_discovery = true;
        experiment.blocked_polish = true;
        d::PiecewiseRootStatistics observed;
        const auto r = d::solve_piecewise_roots(certpolysolve::detail::dense(*p.get()), options, observed, experiment);
        require(r.goal_met && d::disjoint(r.roots) && observed.reflected_roots > 0,
                "conjugate symmetry lost real, nonreal or zero roots");
        for (const auto &box : r.roots)
            require(rectangle_count(symmetric, endpoint(box.endpoints[0]), endpoint(box.endpoints[1]),
                                   endpoint(box.endpoints[2]), endpoint(box.endpoints[3])) == box.count,
                    "reflected dyadic rectangle failed independent winding");
    }
    {
        // A conjugate-filtered candidate pool used to retain only 79 roots.
        // Coordinate movement must not permanently suppress missing discovery.
        std::vector<Complex> planted;
        for (int k = 1; k <= 80; ++k) planted.emplace_back(k);
        const auto p = reference::input(from_roots(planted));
        d::PiecewiseRootConfig experiment;
        experiment.reuse_guards = true;
        experiment.conjugate_symmetry = experiment.reuse_polish = experiment.recover_discovery = true;
        experiment.blocked_polish = true;
        auto query = options;
        query.goal = CPS_ROOT_APPROXIMATION;
        d::PiecewiseRootState state(certpolysolve::detail::dense(*p.get()), experiment);
        query.work_units = 1;
        cps_result r;
        do {
            r = state.advance(query);
            require(r.goal_met || r.progress.work_units_this_call > 0,
                    "planted-root continuation exhausted precision before recovery");
        }
        while (!r.goal_met && state.statistics().recovery_rounds == 0);
        require(!r.goal_met, "planted-root case did not exercise discovery recovery");
        query.cancel = [](void *) -> int { return 1; };
        require(state.advance(query).termination == CPS_CANCELLED,
                "recovery cancellation was ignored");
        query.cancel = nullptr;
        query.work_units = 100000;
        r = state.advance(query);
        const auto &observed = state.statistics();
        require(r.goal_met && d::disjoint(r.roots) && r.roots.size() == planted.size() &&
                    observed.reseeded_models > 0,
                "stalled local discovery did not recover the missing root");
        for (int k = 1; k <= 80; ++k) {
            size_t matches = 0;
            for (const auto &box : r.roots)
                if (box.count == 1 && endpoint(box.endpoints[0]) <= k && k <= endpoint(box.endpoints[1]) &&
                    endpoint(box.endpoints[2]) <= 0 && 0 <= endpoint(box.endpoints[3])) ++matches;
            require(matches == 1, "recovered cover failed independent exact planted-root containment");
        }
    }
    {
        // p(z) = z^64 + (100*z+1)^19 has a tightly grouped set of
        // roots whose candidates need several rounds before certification.
        // Lack of final coverage must not repeatedly regenerate the same local
        // proposals. Coefficients and the cluster containment
        // check use the independent rational backend.
        ComplexPolynomial clustered(65);
        Integer coefficient = 1;
        for (int k = 0; k <= 19; ++k) {
            clustered[static_cast<size_t>(k)] = Complex(Rational(coefficient));
            coefficient *= 100 * (19 - k);
            coefficient /= k + 1;
        }
        clustered[64] = Complex(1);
        const auto p = reference::input(clustered);
        d::PiecewiseRootConfig experiment;
        experiment.reuse_guards = true;
        experiment.conjugate_symmetry = experiment.reuse_polish = experiment.recover_discovery = true;
        experiment.blocked_polish = true;
        auto query = options;
        query.goal = CPS_ROOT_APPROXIMATION;
        d::PiecewiseRootStatistics observed;
        const auto r = d::solve_piecewise_roots(certpolysolve::detail::dense(*p.get()), query, observed, experiment);
        require(r.goal_met && d::disjoint(r.roots) && observed.recovery_rounds > 0 &&
                    observed.repeated_recoveries_avoided > 0,
                "unchanged local discovery was repeated instead of retaining its candidates");
        uint64_t total = 0;
        for (const auto &box : r.roots) total += box.count;
        require(total == 64, "continued refinement lost root accounting");
        // On |z+1/100|=1/200, |z|^64 <= (3/200)^64 < (1/2)^19.
        // Rouche therefore gives exactly 19 roots in this disk. Check that
        // 19 published roots, with their complete rectangles, lie inside it.
        uint64_t grouped = 0;
        for (const auto &box : r.roots) {
            auto squared_extent = [](const Rational &a, const Rational &b) {
                return std::max(Rational(a*a), Rational(b*b));
            };
            const Rational squared_radius =
                squared_extent(endpoint(box.endpoints[0]) + Rational(1)/100,
                               endpoint(box.endpoints[1]) + Rational(1)/100) +
                squared_extent(endpoint(box.endpoints[2]), endpoint(box.endpoints[3]));
            if (squared_radius < Rational(1)/40000) grouped += box.count;
        }
        require(grouped == 19, "continued refinement failed the independent cluster count");
    }
    {
        d::PiecewiseRootState resumed(certpolysolve::detail::dense(*input.get()));
        auto query = options;
        query.work_units = 1;
        while (resumed.statistics().fft_transforms < 2) resumed.advance(query);
        const auto before = resumed.statistics().fft_transforms;
        query.cancel = [](void *) -> int { return 1; };
        const auto stopped = resumed.advance(query);
        require(stopped.termination == CPS_CANCELLED && stopped.progress.work_units_this_call == 0 &&
                    resumed.statistics().fft_transforms == before,
                "cancellation advanced or discarded a suspended FFT batch");
        query.cancel = nullptr;
        query.work_units = 100000;
        const auto complete = resumed.advance(query);
        require(complete.goal_met && stopped.roots.size() == 1,
                "resumption lost the retained batch or changed an old snapshot");
    }
    d::PiecewiseRootStatistics solving;
    auto answer = d::solve_piecewise_roots(certpolysolve::detail::dense(*input.get()), options, solving);
    require(answer.complete && answer.goal_met && d::disjoint(answer.roots),
            "standalone solve incomplete or overlapping");
    for (const auto &r : answer.roots) {
        const auto count = rectangle_count(known, endpoint(r.endpoints[0]), endpoint(r.endpoints[1]),
                                           endpoint(r.endpoints[2]), endpoint(r.endpoints[3]));
        require(count && static_cast<uint64_t>(*count) == r.count,
                "standalone count failed independent winding");
    }
    // Supplied points are only seeds. Duplicate seeds must not multiply counts;
    // an empty pool must never silently call a candidate generator.
    std::vector<d::ComplexBall> supplied;
    for (const auto &r : answer.roots)
        for (int copy = 0; copy < 3; ++copy)
            supplied.push_back(r.box);
    for (const auto goal : {CPS_ROOT_ACCURACY, CPS_ROOT_APPROXIMATION}) {
        cps_options higher = options;
        higher.goal = goal;
        higher.bits = 128;
        d::PiecewiseRootStatistics matching;
        const auto matched = d::certify_piecewise_points(certpolysolve::detail::dense(*input.get()), higher, matching, supplied);
        require(matched.complete && matched.goal_met && d::disjoint(matched.roots),
                "supplied points failed higher precision or deduplication");
        require(matching.guards > 0 && matching.retired_points > 0,
                "overlapping-chart retirement was not exercised");
        for (const auto &r : matched.roots)
            require(rectangle_count(known, endpoint(r.endpoints[0]), endpoint(r.endpoints[1]),
                                    endpoint(r.endpoints[2]), endpoint(r.endpoints[3])) == r.count,
                    "matched certificate failed independent winding");
    }
    cps_options bounded = options;
    bounded.work_units = 8;
    d::PiecewiseRootStatistics empty_stats;
    const auto empty = d::certify_piecewise_points(certpolysolve::detail::dense(*input.get()), bounded, empty_stats, {});
    require(!empty.goal_met && empty_stats.proposals == 0 && empty.roots.size() == 1,
            "empty supplied pool invented roots or invoked a hidden solver");
    d::PiecewiseRootStatistics memory_stats;
    const auto memory_limited =
        d::solve_piecewise_roots(certpolysolve::detail::dense(*input.get()), options, memory_stats, d::PiecewiseRootConfig{0});
    require(!memory_limited.goal_met && memory_limited.termination == CPS_WORK_LIMIT &&
                memory_stats.storage_refusals == 1 && memory_limited.roots.size() == 1,
            "model memory budget was ignored or lost retained certificates");
    const auto repeated_input = reference::input(from_roots({Complex(1), Complex(1), Complex(2)}));
    d::PiecewiseRootStatistics repeated_stats;
    const auto repeated = d::solve_piecewise_roots(certpolysolve::detail::dense(*repeated_input.get()), bounded, repeated_stats);
    require(!repeated.goal_met && std::none_of(repeated.roots.begin(), repeated.roots.end(),
                                             [](const d::Root &root) { return root.multiplicity > 1; }),
            "simple-root certificate invented nonzero repeated-root multiplicity");
    certpolysolve::Polynomial uncertain;
    uncertain.set(0, certpolysolve::ball(certpolysolve::Number("-1"), certpolysolve::Number("1/100"), certpolysolve::Number("0")));
    uncertain.set(1, certpolysolve::ball(certpolysolve::Number("1"), certpolysolve::Number("0"), certpolysolve::Number("0")));
    d::PiecewiseRootStatistics uncertain_stats;
    require(!d::solve_piecewise_roots(certpolysolve::detail::dense(*uncertain.get()), options, uncertain_stats).goal_met,
            "exact piecewise method accepted an uncertain polynomial");

    // Uniformly scaled roots exercise coordinate conversion beyond every
    // hardware floating-point exponent range, with exact input dyadics.
    certpolysolve::Polynomial large;
    large.set(0, certpolysolve::ball(certpolysolve::Number("2@8192"), certpolysolve::Number("0"), certpolysolve::Number("0")));
    large.set(1, certpolysolve::ball(certpolysolve::Number("-3@4096"), certpolysolve::Number("0"), certpolysolve::Number("0")));
    large.set(2, certpolysolve::ball(certpolysolve::Number("1"), certpolysolve::Number("0"), certpolysolve::Number("0")));
    d::PiecewiseRootStatistics large_stats;
    const auto large_answer = d::solve_piecewise_roots(certpolysolve::detail::dense(*large.get()), options, large_stats);
    require(large_answer.complete && large_answer.roots.size() == 2, "large exponent roots did not complete");
    for (const auto &r : large_answer.roots) {
        d::ComplexBall scaled = r.box;
        acb_mul_2exp_si(scaled.value(), scaled.value(), -4096);
        require(arb_contains_zero(acb_imagref(scaled.value())) &&
                    (arb_contains_si(acb_realref(scaled.value()), 1) ||
                     arb_contains_si(acb_realref(scaled.value()), 2)),
                "large exponent coordinate conversion excluded exact root");
    }
    options.cancel = [](void *) -> int { return 1; };
    auto cancelled = d::solve_piecewise_roots(certpolysolve::detail::dense(*input.get()), options, solving);
    require(!cancelled.complete && cancelled.termination == CPS_CANCELLED && cancelled.roots.size() == 1,
            "cancellation lost exact zero certificate");
}
// Chart exclusion: a dense polynomial of degree 300 with integer coefficients
// from a fixed generator (its roots crowd the unit circle). Many charts of its
// rings hold no root; with the exclusion (Pellet's test after root-squaring
// steps) they run no proposal iteration, and every published rectangle is the
// one the iteration finds without it.
static void chart_exclusion_keeps_the_answer() {
    constexpr int degree = 300;
    certpolysolve::Polynomial dense;
    uint64_t state = 20261005;
    for (int j = 0; j <= degree; ++j) {
        state = state * 6364136223846793005ULL + 1442695040888963407ULL;
        const auto value = static_cast<long>((state >> 33) % 2001) - 1000;
        dense.set(static_cast<size_t>(j), certpolysolve::Number(std::to_string(value == 0 ? 1 : value)));
    }
    cps_options options;
    options.goal = CPS_ROOT_APPROXIMATION;
    options.bits = 53;
    auto with = d::maintained_piecewise_config(), without = with;
    without.chart_exclusion = false;
    d::PiecewiseRootStatistics excluded, iterated;
    const auto kept = d::solve_piecewise_roots(certpolysolve::detail::dense(*dense.get()), options, excluded, with);
    const auto plain = d::solve_piecewise_roots(certpolysolve::detail::dense(*dense.get()), options, iterated, without);
    require(kept.goal_met && plain.goal_met && kept.roots.size() == plain.roots.size(),
            "the dense control was not solved");
    require(excluded.excluded_charts > 0 && excluded.machine_sweeps < iterated.machine_sweeps,
            "no chart without a root was excluded from the proposal iteration");
    for (size_t j = 0; j < kept.roots.size(); ++j)
        for (size_t k = 0; k < 4; ++k)
            require(kept.roots[j].endpoints[k].mantissa_text() == plain.roots[j].endpoints[k].mantissa_text() &&
                        kept.roots[j].endpoints[k].exponent_text() == plain.roots[j].endpoints[k].exponent_text(),
                    "chart exclusion changed a published rectangle");
}
// A root lies in the kept disks of several charts, and every chart after the first proposes it again. Those
// copies were polished in ball arithmetic and then retired; they are dropped on the iteration's own point, before
// polishing, and the answer is the same set of roots.
static void copies_are_dropped_before_polishing() {
    constexpr int degree = 300;
    certpolysolve::Polynomial dense;
    uint64_t state = 20261008;
    for (int j = 0; j <= degree; ++j) {
        state = state * 6364136223846793005ULL + 1442695040888963407ULL;
        const auto value = static_cast<long>((state >> 33) % 2001) - 1000;
        dense.set(static_cast<size_t>(j), certpolysolve::Number(std::to_string(value == 0 ? 1 : value)));
    }
    cps_options options;
    options.goal = CPS_ROOT_APPROXIMATION;
    options.bits = 53;
    d::PiecewiseRootStatistics stats;
    const auto result = d::solve_piecewise_roots(certpolysolve::detail::dense(*dense.get()), options, stats, d::maintained_piecewise_config());
    require(result.goal_met && result.roots.size() == static_cast<size_t>(degree), "the dense control was not solved");
    require(stats.unpolished_copies > 0, "copies of held roots were polished before they were retired");
    require(stats.unpolished_copies <= stats.retired_points + stats.owned_proposals_retired,
            "a dropped copy was not counted where a retired proposal is");
    // A candidate's coarse guard and its final certificate are proved at one point of one model: the final
    // certificate takes its short Taylor jet from the guard's longer one.
    require(stats.certificate.jets_reused > 0, "a final certificate recomputed the jet its guard had computed");
}
// A chart's iteration starts one point for each root its model has within twice the kept disk (the Newton
// polygon's count at radius two) and still returns every root of the kept disk: six roots inside the unit disk,
// four at radius 3/2 and twenty-six at radius three make ten points, where twice the dominant index plus sixteen
// made sixteen or more.
static void seeds_are_counted_within_twice_the_kept_disk() {
    constexpr slong precision = 128;
    const std::vector<std::complex<double>> inner{{0.1, 0.05}, {0.75, 0}, {-0.6, 0.3},
                                                  {0.2, -0.7}, {-0.3, -0.55}, {0.5, 0.5}};
    d::RootModel model;
    d::ComplexPolynomial factor;
    d::ComplexBall value;
    acb_poly_one(model.polynomial.value());
    for (const auto z : inner) {
        acb_poly_zero(factor.value());
        acb_set_d_d(value.value(), -z.real(), -z.imag());
        acb_poly_set_coeff_acb(factor.value(), 0, value.value());
        acb_poly_set_coeff_si(factor.value(), 1, 1);
        acb_poly_mul(model.polynomial.value(), model.polynomial.value(), factor.value(), precision);
    }
    // w^4 - (3/2)^4 e^(1.2 i) and w^26 - 3^26 e^i: roots at radius 3/2 and at radius three.
    for (const auto &[power, radius, angle] : {std::tuple{4, 1.5, 1.2}, std::tuple{26, 3.0, 1.0}}) {
        acb_poly_zero(factor.value());
        const auto constant = std::polar(std::pow(radius, power), angle);
        acb_set_d_d(value.value(), -constant.real(), -constant.imag());
        acb_poly_set_coeff_acb(factor.value(), 0, value.value());
        acb_poly_set_coeff_si(factor.value(), power, 1);
        acb_poly_mul(model.polynomial.value(), model.polynomial.value(), factor.value(), precision);
    }
    arb_one(model.radius.value());
    arb_one(model.error.value());
    arb_mul_2exp_si(model.error.value(), model.error.value(), -60);
    d::RootVector points(acb_poly_degree(model.polynomial.value()));
    d::PiecewiseRootStatistics stats;
    cps_options options;
    options.goal = CPS_ROOT_APPROXIMATION;
    options.bits = 53;
    const slong kept = d::propose_piecewise_candidates(points, model, options, precision, stats, d::PiecewiseRootConfig{});
    require(stats.machine_seeds == 10,
            "a chart's iteration did not start one point for each root within twice its kept disk");
    require(kept == static_cast<slong>(inner.size()), "the iteration did not keep one point for each root of the kept disk");
    for (const auto z : inner) {
        bool found = false;
        for (slong j = 0; j < kept && !found; ++j)
            found = std::abs(std::complex<double>(arf_get_d(arb_midref(acb_realref(points[j])), ARF_RND_NEAR),
                                                  arf_get_d(arb_midref(acb_imagref(points[j])), ARF_RND_NEAR)) - z) < 1e-9;
        require(found, "a root of the kept disk has no proposal");
    }
}
// A guard is proved from the model's Taylor jet at the candidate, by rungs: the final certificate's order, then
// twice the entries, up to guard_jet_order. A root far from the others (here 0.75, the nearest root 0.56 away, a
// rectangle of half-width 2^-8) is proved by the first rung, and its rectangle is the one the full order proves.
static void a_guard_takes_the_jet_it_needs() {
    constexpr slong precision = 128;
    const std::vector<std::complex<double>> inner{{0.1, 0.05}, {0.75, 0}, {-0.6, 0.3},
                                                  {0.2, -0.7}, {-0.3, -0.55}, {0.5, 0.5}};
    d::RootModel model;
    d::ComplexPolynomial factor;
    d::ComplexBall value, point;
    acb_poly_one(model.polynomial.value());
    for (const auto z : inner) {
        acb_poly_zero(factor.value());
        acb_set_d_d(value.value(), -z.real(), -z.imag());
        acb_poly_set_coeff_acb(factor.value(), 0, value.value());
        acb_poly_set_coeff_si(factor.value(), 1, 1);
        acb_poly_mul(model.polynomial.value(), model.polynomial.value(), factor.value(), precision);
    }
    for (const auto &[power, radius, angle] : {std::tuple{4, 1.5, 1.2}, std::tuple{26, 3.0, 1.0}}) {
        acb_poly_zero(factor.value());
        const auto constant = std::polar(std::pow(radius, power), angle);
        acb_set_d_d(value.value(), -constant.real(), -constant.imag());
        acb_poly_set_coeff_acb(factor.value(), 0, value.value());
        acb_poly_set_coeff_si(factor.value(), power, 1);
        acb_poly_mul(model.polynomial.value(), model.polynomial.value(), factor.value(), precision);
    }
    arb_one(model.radius.value());
    arb_one(model.error.value());
    arb_mul_2exp_si(model.error.value(), model.error.value(), -60);
    cps_options query;
    query.goal = CPS_ROOT_ACCURACY;
    query.bits = 53;
    query.absolute_tolerance = d::Number("1@-5");
    acb_set_d_d(point.value(), 0.75, 0);
    d::Root whole, laddered;
    require(d::certify_root_model(whole, model, point.value(), query, precision) && whole.count == 1,
            "the full jet did not prove the guard");
    d::CertificateJet jet;
    require(d::certify_root_model(laddered, model, point.value(), query, precision, nullptr, nullptr, nullptr, nullptr,
                                  0, true, nullptr, &jet) && laddered.count == 1,
            "the ladder did not prove the guard");
    require(jet.jet.size() == static_cast<size_t>(d::certificate_jet_order) + 1,
            "a guard its first rung proves computed further jet orders");
    for (size_t k = 0; k < 4; ++k)
        require(whole.endpoints[k].mantissa_text() == laddered.endpoints[k].mantissa_text() &&
                    whole.endpoints[k].exponent_text() == laddered.endpoints[k].exponent_text(),
                "the ladder proved another rectangle than the full jet");
    // The jet kept is the head of the full one, entry for entry: a final certificate at the same point takes it.
    d::Root final_rectangle;
    cps_options goal;
    goal.goal = CPS_ROOT_APPROXIMATION;
    goal.bits = 53;
    d::PiecewiseRootStatistics stats;
    require(d::certify_root_model(final_rectangle, model, point.value(), goal, precision, nullptr, &stats.certificate, nullptr,
                                  nullptr, 0, false, nullptr, &jet) && stats.certificate.jets_reused == 1,
            "the final certificate did not take the guard's jet");
    d::Root alone;
    require(d::certify_root_model(alone, model, point.value(), goal, precision), "the final certificate failed alone");
    for (size_t k = 0; k < 4; ++k)
        require(alone.endpoints[k].mantissa_text() == final_rectangle.endpoints[k].mantissa_text() &&
                    alone.endpoints[k].exponent_text() == final_rectangle.endpoints[k].exponent_text(),
                "a certificate from the shared jet differs from one from its own");
}
static void mixed_scale_cluster_regression() {
    // MPSolve's kam3_2 and kam2_2: A*x^9 + (A*x^2-3)^2,
    // or i*A*x^9 + (A*x^2-3)^2, with A=10^40.
    // Two extremely close pairs near zero coexist with five large roots.
    // Deferred polishing lost completion here; the retained optimizations
    // must preserve the bounded solve and its counted approximation contract.
    Integer magnitude = 1;
    for (int i = 0; i < 40; ++i) magnitude *= 10;
    const Rational a(magnitude);
    ComplexPolynomial polynomial(10);
    polynomial[0] = Complex(9);
    polynomial[2] = Complex(-6*a);
    polynomial[4] = Complex(a*a);
    cps_options options;
    options.goal = CPS_ROOT_APPROXIMATION;
    options.bits = 53;
    options.max_work_bits = 1024;
    options.work_units = 50000;
    for (const auto &leading : {Complex(a), Complex(0, a)}) {
        polynomial[9] = leading;
        const auto input = reference::input(polynomial);
        for (const bool optimized : {false, true}) {
            d::PiecewiseRootConfig config;
            config.reuse_guards = optimized;
            d::PiecewiseRootStatistics statistics;
            const auto answer = d::solve_piecewise_roots(certpolysolve::detail::dense(*input.get()), options, statistics, config);
            require(answer.goal_met && d::disjoint(answer.roots), "mixed-scale clustered solve lost completion");
            uint64_t total = 0;
            for (const auto &r : answer.roots) {
                total += r.count;
                require(rectangle_count(polynomial, endpoint(r.endpoints[0]), endpoint(r.endpoints[1]),
                                        endpoint(r.endpoints[2]), endpoint(r.endpoints[3])) == r.count,
                        "mixed-scale output failed independent exact winding");
            }
            require(total == 9, "mixed-scale output lost a root or multiplicity");
        }
    }
}
// Deflated discovery: the clustered Mahler polynomial prod_(k<5) (1 + a x^(2^k)), a = -(1 + 2^-20), of
// degree 31, times 2^100: the coefficient of x^m is (-(2^20 + 1))^s 2^(20 (5 - s)), s the binary digit
// sum of m. Each primitive 2^r-th root of unity gathers 5 - r roots at distances about 2^-(20 + k)
// from it. The first round's binary64 proposals leave two members of such clusters without a seed,
// the charts there cannot count their roots, and without the deflated seeding the solve climbed the
// precision ladder discovering nothing more (29 of 31 roots, 32,885 bits after 60 s). The regime exit
// is off: it would hand the stall to a global method.
static void deflated_discovery_finds_cluster_members() {
    constexpr int levels = 5, scale_bits = 20;
    const Integer scale = Integer(1) << scale_bits, factor = -(scale + 1);
    ComplexPolynomial polynomial(size_t{1} << levels);
    for (size_t m = 0; m < polynomial.size(); ++m) {
        const int digits = std::popcount(m);
        Integer coefficient = 1;
        for (int k = 0; k < digits; ++k) coefficient *= factor;
        for (int k = digits; k < levels; ++k) coefficient *= scale;
        polynomial[m] = Complex(Rational(coefficient));
    }
    const auto input = reference::input(polynomial);
    cps_options options;
    options.goal = CPS_ROOT_APPROXIMATION;
    options.bits = 53;
    options.max_work_bits = 1024;
    auto config = d::maintained_piecewise_config();
    config.regime_exit = false;
    d::PiecewiseRootStatistics statistics;
    const auto answer = d::solve_piecewise_roots(certpolysolve::detail::dense(*input.get()), options, statistics, config);
    require(answer.goal_met && d::disjoint(answer.roots),
            "a cluster member was never discovered (deflated seeding)");
    require(statistics.deflation_seeds > 0, "the clustered solve finished without deflated seeding");
    // The roots in closed form: factor k vanishes at (1 + 2^-20)^(-2^-k) exp(2 pi i j 2^-k). Every one
    // lies in exactly one published rectangle and every rectangle counts one root; the rectangles are
    // disjoint, so they isolate all 31. (Exact winding at these coefficient sizes takes minutes.)
    constexpr slong bits = 256;
    std::vector<size_t> hits(answer.roots.size(), 0);
    arb_t base, lo, hi;
    acb_t root, turn;
    arb_init(base); arb_init(lo); arb_init(hi);
    acb_init(root); acb_init(turn);
    fmpz_t mantissa, exponent;
    fmpz_init(mantissa); fmpz_init(exponent);
    const auto bound = [&](arb_t out, const d::Endpoint &e) {
        fmpz_set_str(mantissa, e.mantissa_text().c_str(), 10);
        fmpz_set_str(exponent, e.exponent_text().c_str(), 10);
        arb_set_fmpz(out, mantissa);
        arb_mul_2exp_fmpz(out, out, exponent);
    };
    for (int k = 0; k < levels; ++k)
        for (slong j = 0; j < (slong{1} << k); ++j) {
            arb_one(base);
            arb_mul_2exp_si(base, base, -scale_bits);
            arb_add_ui(base, base, 1, bits);
            arb_inv(base, base, bits);
            arb_root_ui(base, base, ulong{1} << k, bits); // ulong: 1UL is 32 bits on Windows
            acb_set_si(turn, j);
            acb_mul_2exp_si(turn, turn, 1 - k);
            acb_exp_pi_i(turn, turn, bits);
            acb_mul_arb(root, turn, base, bits);
            size_t containing = 0;
            for (size_t r = 0; r < answer.roots.size(); ++r) {
                const auto &e = answer.roots[r].endpoints;
                bound(lo, e[0]); bound(hi, e[1]);
                if (!arb_le(lo, acb_realref(root)) || !arb_le(acb_realref(root), hi)) continue;
                bound(lo, e[2]); bound(hi, e[3]);
                if (!arb_le(lo, acb_imagref(root)) || !arb_le(acb_imagref(root), hi)) continue;
                ++containing;
                ++hits[r];
            }
            require(containing == 1, "a Mahler root is not in exactly one published rectangle");
        }
    fmpz_clear(mantissa); fmpz_clear(exponent);
    arb_clear(base); arb_clear(lo); arb_clear(hi);
    acb_clear(root); acb_clear(turn);
    for (size_t r = 0; r < answer.roots.size(); ++r)
        require(answer.roots[r].count == 1 && hits[r] == 1, "a rectangle does not hold exactly one Mahler root");
}
// The certificate attempt's budget (machine_newton_cover's fallback_operations):
// the attempt as a whole spends at most what the caller's route costs if it fails. Sixteen rational points on the
// unit circle and a real pair 2^-20 apart at 2: binary64 leaves the pair unresolved, the certificates refuse it, and
// one repair round at the certificate's precision separates it (a pair 2^-24 apart needs more than one round's
// four sweeps, and the attempt ends there by its halving rule). Under the default budget (the global recipes'
// minimum) the round runs and every root is certified; under a budget smaller than the pass itself no round is
// paid for and the attempt refuses.
static void certificate_attempt_budget() {
    std::vector<Complex> roots;
    for (int k = 1; k <= 8; ++k) {
        const Rational t = Rational(k) / 9, x = (1 - t * t) / (1 + t * t), y = 2 * t / (1 + t * t);
        roots.push_back(Complex(x, y));
        roots.push_back(Complex(x, -y));
    }
    roots.push_back(Complex(Rational(2)));
    roots.push_back(Complex(Rational(2) + power2(-20)));
    const auto p = reference::input(from_roots(roots));
    d::PolynomialSource source([&](acb_poly_t out, slong bits) {
        d::ComplexBall a;
        acb_poly_zero(out);
        for (const auto &[j, c] : p.get()->terms) {
            c.enclose(a.value(), bits);
            acb_poly_set_coeff_acb(out, static_cast<slong>(j), a.value());
        }
    }, static_cast<slong>(roots.size()), true);
    cps_options query;
    query.goal = CPS_ROOT_APPROXIMATION;
    query.bits = 53;
    query.max_work_bits = 4096;
    // The rounds first, without the double-double repair (which proves the pair with no round); then the repair,
    // which the default budget pays for and a budget of one operation does not.
    for (const bool repair : {false, true})
        for (const uint64_t budget : {uint64_t{0}, uint64_t{1}}) {
            d::EnginePolicy observed;
            observed.double_double_repair = repair;
            observed.study_statistics = std::make_shared<d::StudyStatistics>();
            certpolysolve::detail::CandidateStatistics work{};
            std::vector<d::Root> cover;
            const bool certified = d::machine_newton_cover(cover, source, nullptr, query, 85, true, work, true,
                                                           &observed, nullptr, nullptr, nullptr, true, nullptr, {},
                                                           budget);
            const auto rounds = observed.study_statistics->machine_cover_rounds;
            const auto repairs = observed.study_statistics->machine_cover_double_double_repairs;
            if (!repair && budget == 0)
                require(certified && cover.size() == roots.size() && rounds >= 1 && repairs == 0,
                        "the attempt did not repair the close pair within its default budget");
            else if (!repair)
                require(!certified && rounds == 0, "a repair round ran beyond the attempt's budget");
            else if (budget == 0)
                require(certified && cover.size() == roots.size() && rounds == 0 && repairs == 1,
                        "the double-double repair did not prove the close pair within the default budget");
            else
                require(!certified && rounds == 0 && repairs == 0,
                        "the double-double repair ran beyond the attempt's budget");
        }
    // At a 3,322-bit goal, four roots 2^-8 from 2 beside the sixteen on the
    // unit circle need their certificates at twice the attempt's precision,
    // 6,708 bits. In a retry after the staged solve (every point 2^-12 from
    // its root, none settled) the fallback is the secular refinement at that
    // precision, and the budget is priced there: one round on the refused
    // points completes the cover. Priced at the attempt's 3,354 bits, the
    // certificate pass spent the budget, no round ran, and 4 of 20 points
    // went to the refinement (jacobi-characteristic-2 at 3,322 bits).
    roots.resize(16);
    for (const Complex direction : {Complex(1), Complex(0, 1), Complex(-1), Complex(0, -1)})
        roots.push_back(Complex(Rational(2)) + direction * Complex(power2(-8)));
    const auto q = reference::input(from_roots(roots));
    d::PolynomialSource clustered([&](acb_poly_t out, slong bits) {
        d::ComplexBall a;
        acb_poly_zero(out);
        for (const auto &[j, c] : q.get()->terms) {
            c.enclose(a.value(), bits);
            acb_poly_set_coeff_acb(out, static_cast<slong>(j), a.value());
        }
    }, static_cast<slong>(roots.size()), true);
    const slong n = static_cast<slong>(roots.size());
    d::RootVector staged(n);
    for (slong i = 0; i < n; ++i) {
        const Complex z = roots[static_cast<size_t>(i)] +
                          Complex(power2(-12) * Rational(1 + i % 3), power2(-12) * Rational(1 + i % 5));
        d::Number(z.real.str()).enclose(acb_realref(staged.data() + i), 3354);
        d::Number(z.imag.str()).enclose(acb_imagref(staged.data() + i), 3354);
        acb_get_mid(staged.data() + i, staged.data() + i);
    }
    cps_options high = query;
    high.bits = 3322;
    high.max_work_bits = 65536;
    {
        d::EnginePolicy observed;
        observed.study_statistics = std::make_shared<d::StudyStatistics>();
        certpolysolve::detail::CandidateStatistics work{};
        std::vector<d::Root> cover;
        d::PartialNewtonCover partial;
        const std::vector<uint8_t> none(static_cast<size_t>(n), 0);
        d::machine_newton_cover(cover, clustered, staged.data(), high, 3354, true, work, true, &observed, nullptr,
                                nullptr, nullptr, true, &partial, none, 0);
        require(partial.held() == static_cast<size_t>(n) && observed.study_statistics->machine_cover_rounds >= 1,
                "in a retry at a high goal the budget stopped the round that completes the cover");
    }
    // A first attempt at the same goal leaves its refusals to the staged
    // solve, which repairs them at 85-170 bits: its budget stays priced at
    // its first precision and pays for a round on one of the four. Priced at
    // the certificate precision it bought rounds at 6,708 bits that cost
    // more than the staged solve (resultant-769 scaled at 3,322 bits:
    // 14 -> 27 s).
    {
        d::EnginePolicy observed;
        observed.study_statistics = std::make_shared<d::StudyStatistics>();
        certpolysolve::detail::CandidateStatistics work{};
        std::vector<d::Root> cover;
        const bool certified = d::machine_newton_cover(cover, clustered, nullptr, high, 3354, true, work, true,
                                                       &observed, nullptr, nullptr, nullptr, true, nullptr, {}, 0);
        require(!certified && observed.study_statistics->machine_cover_rounds == 1,
                "a first attempt at a high goal spent beyond its first precision's budget");
    }
}
// The certificate attempt's binary64 packet goes on past its first quantum (heuristics::machine_cover_sweeps) while
// its quanta settle points and the attempt's budget pays (machine_newton_cover). From the equally spaced start,
// x^8192 - 3 x^4097 + 1 (two rings of roots a fraction of a spacing apart) needs about 110 sweeps: 1,698 points
// were still moving after 32 sweeps, 1,080 after 64, 556 after 96. Cut at one quantum the packet left them to the
// repair rounds, which refused the attempt after 12.7 s; with the quanta the packet converges and every point is
// certified (about 2 s).
static void packet_quanta() {
    const slong n = 8192, k = 4097;
    d::ComplexPolynomial p;
    acb_poly_set_coeff_si(p.value(), 0, 1);
    acb_poly_set_coeff_si(p.value(), k, -3);
    acb_poly_set_coeff_si(p.value(), n, 1);
    d::RootVector start(n);
    require(d::initial_candidates(start.data(), p.value(), 64), "no starting points");
    d::PolynomialSourceOptions sparse;
    sparse.sparse_coefficients = true;
    sparse.input = d::InputSemantics::Exact;
    d::PolynomialSource source([&](acb_poly_t out, slong) { acb_poly_set(out, p.value()); }, n, true, sparse);
    cps_options query;
    query.goal = CPS_ROOT_APPROXIMATION;
    query.bits = 53;
    query.max_work_bits = 4096;
    certpolysolve::detail::CandidateStatistics work{};
    std::vector<d::Root> cover;
    bool converged = false;
    const bool certified = d::machine_newton_cover(cover, source, start.data(), query, 85, true, work, true, nullptr,
                                                   nullptr, &converged);
    require(converged, "the packet ended unconverged (cut at its first quantum)");
    require(certified && cover.size() == static_cast<size_t>(n), "the attempt did not certify the trinomial");
}
// A sparse polynomial's certificate attempt starts from its Newton-polygon edges' binomial roots (initial_candidates,
// heuristics::binomial_seed_spacing_fraction). x^32768 - 3 x^10923 + 1 has two rings of roots half a root spacing
// apart; from the equally spaced start the attempt was refused after 16.6 s, from the binomial start it certifies
// (about 1.3 s). A real edge's binomial start is turned off the
// conjugation symmetry a Jacobi sweep of a real source would keep for ever.
static void sparse_binomial_seeds() {
    {
        d::ComplexPolynomial p;
        acb_poly_set_coeff_si(p.value(), 0, 1);
        acb_poly_set_coeff_si(p.value(), 5, -3);
        acb_poly_set_coeff_si(p.value(), 16, 1);
        d::RootVector seeds(16);
        require(d::initial_candidates(seeds.data(), p.value(), 64, true), "no binomial start for a trinomial");
        const auto at = [&](slong i) {
            return std::complex<double>(arf_get_d(arb_midref(acb_realref(seeds[i])), ARF_RND_NEAR),
                                        arf_get_d(arb_midref(acb_imagref(seeds[i])), ARF_RND_NEAR));
        };
        bool closed = true;
        for (slong i = 0; i < 16 && closed; ++i) {
            bool mirrored = false;
            for (slong j = 0; j < 16; ++j)
                mirrored = mirrored || std::abs(std::conj(at(i)) - at(j)) < 1e-9;
            closed = mirrored;
        }
        require(!closed, "a real edge's binomial start is closed under conjugation");
    }
    const slong n = 32768, k = 10923;
    d::PolynomialSource source([&](acb_poly_t out, slong) {
        acb_poly_zero(out);
        acb_poly_set_coeff_si(out, 0, 1);
        acb_poly_set_coeff_si(out, k, -3);
        acb_poly_set_coeff_si(out, n, 1);
    }, n, true);
    cps_options query;
    query.goal = CPS_ROOT_APPROXIMATION;
    query.bits = 53;
    query.max_work_bits = 4096;
    certpolysolve::detail::CandidateStatistics work{};
    std::vector<d::Root> cover;
    require(d::machine_newton_cover(cover, source, nullptr, query, 85, true, work, true) &&
                cover.size() == static_cast<size_t>(n),
            "the attempt did not certify the trinomial from its binomial start");
}
int main() {
    certpolysolve::ThreadScope lifetime;
    try {
        check();
        chart_exclusion_keeps_the_answer();
        copies_are_dropped_before_polishing();
        seeds_are_counted_within_twice_the_kept_disk();
        a_guard_takes_the_jet_it_needs();
        mixed_scale_cluster_regression();
        deflated_discovery_finds_cluster_members();
        certificate_attempt_budget();
        packet_quanta();
        sparse_binomial_seeds();
        std::cout << "piecewise original-input certificates and standalone solve passed\n";
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
