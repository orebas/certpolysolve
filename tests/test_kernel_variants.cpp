// A kernel's variants against its reference (platform/kernel_variant.hpp): a variant must return the reference's
// bits, so that one binary gives the same answers on every CPU it runs on. Every variant this machine
// runs is compared; on a machine without AVX2 and FMA the comparisons are of the reference with itself and say so.
//
// 1. The choice: the reference is always available; set_kernel_variant returns the choice before and refuses a
//    variant the machine does not run.
// 2. The compensated certificate's loop: compensated_kernel::avx2_fma against compensated_kernel::reference, both the
//    library's own objects, on all eleven numbers, for coefficients of every size (two-word, wide exponents,
//    subnormal, zero, real) and points inside, on and outside the unit circle, on the axes, tiny and so large that
//    the sums overflow. A result that is not a number must be not a number in both.
// 3. CompensatedCertificate's evaluations: evaluate_many on 1 to 9 points under each variant against evaluate a
//    point under the reference, every field; certify with what evaluate_ahead made (the admission, the evaluation,
//    the proof's geometry and its bound on the second derivative) against certify alone: the same verdict, the same
//    rectangle and the same kept evaluation, at two goals and precisions, also on a polynomial whose points the
//    admission refuses; values made for another point or precision are not used.
// 4. The iteration's dense loop: horner_kernel::avx2 against PreparedJets<double, 1>::evaluate (this file's build of
//    the header, with the library's floating-point flags), every field and the floating-point flags raised.
// 5. PreparedNewton<double>::evaluate_range against evaluate a point: points on both sides of the unit circle, points
//    left out, a polynomial with roots at zero (a shorter reversed polynomial) and a sparse one (no dense loop).
// 6. The machine engine itself, the library's own build of the sweep: packets under each variant end on the same
//    points, with the same progress and the same work counters.
// 8. The certificate's majorants: compensated_kernel::majorants_avx2 against compensated_kernel::majorant, for every
//    lower end the certificate uses, coefficients of every size and points whose sums overflow.
// 7. The direct interaction sum: interaction_kernel::add_groups, with a term at a time for what it leaves (as the
//    sweep calls it), against a term at a time throughout (this file's build of machine::Complex's operators): the
//    sum, the count of terms, the verdict at a point equal to z, and the floating-point flags, for every point of
//    sets of every size up to 13 and larger ones, with differences on the axes, of equal parts, tiny and huge, and a
//    coincident pair at every place.
#include "certify/attempts/compensated_certificate.hpp"
#include "certify/attempts/compensated_kernel.hpp"
#include "core/cover.hpp"
#include "core/engine_statistics.hpp"
#include "core/heuristics.hpp"
#include "core/raii.hpp"
#include "numeric/candidates/candidate_jets.hpp"
#include "numeric/candidates/candidates.hpp"
#include "numeric/candidates/aberth_step.hpp"
#include "numeric/candidates/horner_kernel.hpp"
#include "numeric/candidates/interaction_kernel.hpp"
#include "numeric/candidates/machine_candidates.hpp"
#include "platform/kernel_variant.hpp"
#include <flint/acb_poly.h>
#include <bit>
#include <cfenv>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace d = certpolysolve::detail;
namespace m = certpolysolve::detail::machine;

namespace {
void require(bool ok, const std::string &why) {
    if (!ok)
        throw std::runtime_error(why);
}
// The same binary64 number: the same bits, or both not a number (which operand's payload a NaN keeps is not part
// of what a kernel promises, and no caller reads it: a NaN ends the evaluation or the packet).
bool same(double a, double b) {
    return std::bit_cast<uint64_t>(a) == std::bit_cast<uint64_t>(b) || (std::isnan(a) && std::isnan(b));
}
bool same(m::Complex<double> a, m::Complex<double> b) { return same(a.re, b.re) && same(a.im, b.im); }

// The variants this machine runs, the reference first.
std::vector<d::KernelVariant> variants() {
    std::vector<d::KernelVariant> out{d::KernelVariant::Reference};
    if (d::kernel_variant_available(d::KernelVariant::Avx2Fma))
        out.push_back(d::KernelVariant::Avx2Fma);
    return out;
}
// A scope that runs the kernels in one variant and gives the choice back.
class Under {
    d::KernelVariant before_;

  public:
    explicit Under(d::KernelVariant variant) : before_(d::set_kernel_variant(variant)) {
        require(d::kernel_variant() == variant, "the variant asked for is not in force");
    }
    ~Under() { d::set_kernel_variant(before_); }
};

// Points of every kind a kernel meets: moduli from tiny to overflowing, both sides of the unit circle, the axes.
std::vector<m::Complex<double>> sample_points(std::mt19937_64 &rng, size_t random_points) {
    std::vector<m::Complex<double>> points = {{1, 0},      {0, 1},      {-1, 0},     {0, -1},   {0.5, 0},
                                              {0, -0.25},  {2, 0},      {0, 3},      {1e-200, 1e-201},
                                              {1e150, -1e150}, {0x1p-1040, 0}, {0.70710678118654757, 0.70710678118654757},
                                              {1, 0x1p-52}, {0.99999999999999989, 0}};
    std::uniform_real_distribution<double> unit(-1, 1);
    for (const double modulus : {1e-3, 0.5, 0.97, 1.0, 1.03, 2.0, 1e3})
        for (size_t k = 0; k < random_points; ++k) {
            const double angle = 3.141592653589793 * unit(rng), r = modulus * (1 + 0.01 * unit(rng));
            points.push_back({r * std::cos(angle), r * std::sin(angle)});
        }
    return points;
}

// 1.
void choice() {
    require(d::kernel_variant_available(d::KernelVariant::Reference), "the reference is not available");
    const d::KernelVariant start = d::kernel_variant();
    require(d::kernel_variant_available(start), "the variant in force is not available");
    const bool wide = d::kernel_variant_available(d::KernelVariant::Avx2Fma);
    require(start == (wide ? d::KernelVariant::Avx2Fma : d::KernelVariant::Reference),
            "the first choice is not the best available variant");
    require(d::set_kernel_variant(d::KernelVariant::Reference) == start, "set_kernel_variant lost the choice before");
    require(d::kernel_variant() == d::KernelVariant::Reference, "the reference was not taken");
    require(d::set_kernel_variant(d::KernelVariant::Avx2Fma) == d::KernelVariant::Reference,
            "set_kernel_variant lost the choice before");
    require(d::kernel_variant() == (wide ? d::KernelVariant::Avx2Fma : d::KernelVariant::Reference),
            "a variant the machine does not run was taken, or one it runs was refused");
    d::set_kernel_variant(start);
}

// 2.
using Term = d::compensated_kernel::Term;
using CompensatedSums = d::compensated_kernel::Sums;
bool same(const CompensatedSums &a, const CompensatedSums &b) {
    return same(a.sr, b.sr) && same(a.si, b.si) && same(a.cr, b.cr) && same(a.ci, b.ci) && same(a.dr, b.dr) &&
           same(a.di, b.di) && same(a.terms, b.terms) && same(a.terms_slope, b.terms_slope) &&
           same(a.values, b.values) && same(a.errors, b.errors) && same(a.errors_slope, b.errors_slope);
}
std::vector<Term> random_terms(std::mt19937_64 &rng, size_t degree, int kind) {
    std::uniform_real_distribution<double> unit(-1, 1);
    std::uniform_int_distribution<int> exponent(-200, 200);
    std::vector<Term> terms(degree + 1);
    for (Term &t : terms) {
        double scale = 1;
        if (kind == 1)
            scale = 0x1p-300;
        else if (kind == 2)
            scale = 0x1p300;
        else if (kind == 3)
            scale = std::ldexp(1.0, exponent(rng));
        else if (kind == 4)
            scale = 0x1p-1060; // subnormal coefficients
        t.re_hi = unit(rng) * scale;
        t.im_hi = kind == 5 ? 0 : unit(rng) * scale;
        t.re_lo = unit(rng) * scale * 0x1p-54;
        t.im_lo = kind == 5 ? 0 : unit(rng) * scale * 0x1p-54;
        t.error = std::abs(unit(rng)) * scale * 0x1p-108;
        if (kind == 6 && (rng() & 3) == 0)
            t = {};
    }
    if (terms.back().re_hi == 0)
        terms.back().re_hi = 0.75;
    return terms;
}
uint64_t compensated_loop() {
#ifdef CERTPOLYSOLVE_HAVE_AVX2_FMA_KERNELS
    if (!d::kernel_variant_available(d::KernelVariant::Avx2Fma))
        return 0;
    std::mt19937_64 rng(2026111701);
    uint64_t compared = 0;
    constexpr size_t group = d::compensated_kernel::group;
    for (const size_t degree : {1, 2, 3, 4, 5, 8, 13, 64, 257})
        for (int kind = 0; kind <= 6; ++kind) {
            const auto terms = random_terms(rng, degree, kind);
            const auto points = sample_points(rng, 2);
            for (size_t first = 0; first < points.size(); ++first) {
                // Four different points a call: a lane's result must not depend on its neighbours.
                double xr[group], xi[group], X[group];
                for (size_t l = 0; l < group; ++l) {
                    const auto z = points[(first + l * 5) % points.size()];
                    xr[l] = z.re, xi[l] = z.im;
                    X[l] = std::hypot(z.re, z.im) * (1 + 0x1p-50);
                }
                CompensatedSums wide[group];
                d::compensated_kernel::avx2_fma(terms.data(), degree, xr, xi, X, wide);
                for (size_t l = 0; l < group; ++l) {
                    const auto one = d::compensated_kernel::reference(terms.data(), degree, xr[l], xi[l], X[l]);
                    require(same(one, wide[l]),
                            "the compensated loop's four-point variant differs from the reference: degree " +
                                std::to_string(degree) + ", coefficients of kind " + std::to_string(kind) +
                                ", point " + std::to_string(first) + ", lane " + std::to_string(l));
                    ++compared;
                }
            }
        }
    return compared;
#else
    return 0;
#endif
}

// 3.
bool same(const d::CompensatedEvaluation &a, const d::CompensatedEvaluation &b) {
    return same(a.x_re, b.x_re) && same(a.x_im, b.x_im) && same(a.value_re, b.value_re) &&
           same(a.value_im, b.value_im) && same(a.correction_re, b.correction_re) &&
           same(a.correction_im, b.correction_im) && same(a.slope_re, b.slope_re) && same(a.slope_im, b.slope_im) &&
           same(a.value_error, b.value_error) && same(a.slope_error, b.slope_error) && same(a.modulus, b.modulus) &&
           same(a.amplification, b.amplification) && a.ok == b.ok;
}
cps_options goal(uint64_t bits) {
    cps_options o;
    o.bits = bits;
    o.goal = CPS_ROOT_APPROXIMATION;
    return o;
}
// prod (x - root_j) at 1,024 bits, the roots dyadic: m 2^-shift.
struct Planted {
    d::ComplexPolynomial polynomial;
    std::vector<d::ComplexBall> roots;
};
Planted planted(std::mt19937_64 &rng, int degree, slong shift, bool real_axis) {
    Planted out;
    acb_poly_one(out.polynomial.value());
    std::uniform_int_distribution<int64_t> mantissa(-(int64_t(1) << 40), int64_t(1) << 40);
    d::ComplexPolynomial factor;
    for (int j = 0; j < degree; ++j) {
        d::ComplexBall root, negated;
        acb_set_si_si(root.value(), mantissa(rng), real_axis || j % 3 == 0 ? 0 : mantissa(rng));
        acb_mul_2exp_si(root.value(), root.value(), -shift);
        acb_neg(negated.value(), root.value());
        acb_poly_zero(factor.value());
        acb_poly_set_coeff_si(factor.value(), 1, 1);
        acb_poly_set_coeff_acb(factor.value(), 0, negated.value());
        acb_poly_mul(out.polynomial.value(), out.polynomial.value(), factor.value(), 1024);
        out.roots.push_back(std::move(root));
    }
    return out;
}
// The points a later try proved in certify_ahead's comparisons of the tries made one at a time: the Newton retry
// must have been exercised.
uint64_t proved_on_a_later_try = 0;
// certify with what evaluate_ahead made, under every variant, against certify alone under the reference: the same
// verdict, the same kept evaluation and, where it proves, the same rectangle. Returns the points certify alone proves.
uint64_t certify_ahead(const d::CompensatedCertificate &tier, const std::vector<m::Complex<double>> &points,
                       const cps_options &o, slong precision, uint64_t &compared) {
    const auto n = static_cast<slong>(points.size());
    d::RootVector balls(n);
    std::vector<slong> indices;
    for (slong i = 0; i < n; ++i) {
        acb_set_d_d(balls.data() + i, points[static_cast<size_t>(i)].re, points[static_cast<size_t>(i)].im);
        indices.push_back(i);
    }
    std::vector<d::Root> alone(points.size());
    std::vector<d::CompensatedEvaluation> alone_seen(points.size());
    std::vector<uint8_t> alone_proved(points.size());
    uint64_t proved = 0;
    {
        const Under reference(d::KernelVariant::Reference);
        for (size_t i = 0; i < points.size(); ++i) {
            alone_proved[i] = tier.certify(alone[i], balls.data() + i, o, precision, &alone_seen[i]);
            proved += alone_proved[i];
        }
    }
    for (const auto variant : variants()) {
        const Under scope(variant);
        std::vector<d::CompensatedAhead> ahead(points.size());
        tier.evaluate_ahead(balls.data(), indices, o, precision, ahead);
        for (size_t i = 0; i < points.size(); ++i) {
            d::Root root;
            d::CompensatedEvaluation seen;
            const bool ok = tier.certify(root, balls.data() + i, o, precision, &seen, &ahead[i]);
            require(ok == bool(alone_proved[i]), "certify with the values made ahead gives another verdict");
            require(same(seen, alone_seen[i]), "certify with the values made ahead kept another evaluation");
            if (ok)
                require(acb_equal(root.box.value(), alone[i].box.value()) && root.count == alone[i].count &&
                            root.certificate == alone[i].certificate,
                        "certify with the values made ahead proves another rectangle");
            ++compared;
        }
        // Values made for another point, or at another precision, are not used: certify decides alone.
        if (points.size() >= 2) {
            d::Root root;
            d::CompensatedEvaluation seen;
            bool ok = tier.certify(root, balls.data(), o, precision, &seen, &ahead[1]);
            require(ok == bool(alone_proved[0]) && same(seen, alone_seen[0]),
                    "certify used values made for another point");
            d::Root again;
            ok = tier.certify(again, balls.data(), o, precision + 7, &seen, &ahead[0]);
            d::Root expected;
            d::CompensatedEvaluation expected_seen;
            const bool expected_ok = tier.certify(expected, balls.data(), o, precision + 7, &expected_seen);
            require(ok == expected_ok && same(seen, expected_seen) &&
                        (!ok || acb_equal(again.box.value(), expected.box.value())),
                    "certify used a geometry made at another precision");
        }
        // The tries one at a time, the later tries of all the points made together, as the attempt's first pass
        // makes them: certify's verdict and rectangle at every point.
        {
            std::vector<d::Root> roots(points.size());
            std::vector<uint8_t> ok(points.size());
            std::vector<std::pair<size_t, d::CompensatedRetry>> again, next;
            for (size_t i = 0; i < points.size(); ++i) {
                d::CompensatedRetry retry;
                d::CompensatedEvaluation seen;
                ok[i] = tier.certify_first(roots[i], balls.data() + i, o, precision, &seen, &ahead[i], retry);
                require(same(seen, alone_seen[i]), "certify_first kept another evaluation");
                require(!(ok[i] && retry.available), "certify_first proved a point and asked for another try");
                if (!ok[i] && retry.available)
                    again.emplace_back(i, retry);
            }
            for (unsigned k = 0; k < d::heuristics::compensated_newton_retries && !again.empty(); ++k) {
                std::vector<double> re, im;
                for (const auto &entry : again) {
                    re.push_back(entry.second.re);
                    im.push_back(entry.second.im);
                }
                std::vector<d::CompensatedAhead> later(again.size());
                tier.evaluate_ahead_at(re, im, o, precision, later);
                next.clear();
                for (size_t j = 0; j < again.size(); ++j) {
                    d::CompensatedRetry retry;
                    if (tier.certify_again(roots[again[j].first], again[j].second, o, precision, &later[j], retry)) {
                        ok[again[j].first] = 1;
                        ++proved_on_a_later_try;
                    } else if (retry.available) {
                        next.emplace_back(again[j].first, retry);
                    }
                }
                again.swap(next);
            }
            for (size_t i = 0; i < points.size(); ++i) {
                require(bool(ok[i]) == bool(alone_proved[i]),
                        "the tries made one at a time give another verdict than certify");
                if (ok[i])
                    require(acb_equal(roots[i].box.value(), alone[i].box.value()) &&
                                roots[i].certificate == alone[i].certificate,
                            "the tries made one at a time prove another rectangle than certify");
                ++compared;
            }
        }
    }
    return proved;
}
uint64_t certificate_evaluations() {
    std::mt19937_64 rng(2026111702);
    uint64_t compared = 0;
    for (const int degree : {1, 2, 5, 9, 33, 100})
        for (const slong shift : {40, 39}) { // roots up to 1 and up to 2 in modulus: both sides of the unit circle
            const auto input = planted(rng, degree, shift, degree == 9);
            const d::CompensatedCertificate tier(input.polynomial.value());
            require(tier.available(), "the tier is unavailable on a planted polynomial");
            // Points: near every root (close enough to prove, a little further for the retry, far for a refusal), and
            // the sample of all kinds, non-finite and zero included.
            std::vector<m::Complex<double>> points;
            for (const auto &root : input.roots) {
                const double rr = arf_get_d(arb_midref(acb_realref(root.value())), ARF_RND_NEAR);
                const double ri = arf_get_d(arb_midref(acb_imagref(root.value())), ARF_RND_NEAR);
                for (const double t : {0.0, 1e-15, 1e-9, 1e-7, 1e-6, 1e-5, 1e-4})
                    points.push_back({rr * (1 + t), ri * (1 - t)});
            }
            for (const auto &z : sample_points(rng, 1))
                points.push_back(z);
            points.push_back({0, 0});
            points.push_back({std::numeric_limits<double>::infinity(), 1});
            points.push_back({1, std::numeric_limits<double>::quiet_NaN()});
            // evaluate a point, under the reference.
            std::vector<d::CompensatedEvaluation> expected(points.size());
            {
                const Under reference(d::KernelVariant::Reference);
                for (size_t i = 0; i < points.size(); ++i)
                    expected[i] = tier.evaluate(points[i].re, points[i].im);
            }
            for (const auto variant : variants()) {
                const Under scope(variant);
                for (size_t i = 0; i < points.size(); ++i)
                    require(same(expected[i], tier.evaluate(points[i].re, points[i].im)),
                            "evaluate under a variant differs from the reference");
                for (size_t count = 1; count <= 9; ++count)
                    for (size_t first = 0; first + count <= points.size(); first += 7) {
                        std::vector<double> xr(count), xi(count);
                        for (size_t k = 0; k < count; ++k)
                            xr[k] = points[first + k].re, xi[k] = points[first + k].im;
                        std::vector<d::CompensatedEvaluation> many(count);
                        tier.evaluate_many(xr, xi, many);
                        for (size_t k = 0; k < count; ++k) {
                            require(same(expected[first + k], many[k]),
                                    "evaluate_many differs from evaluate: degree " + std::to_string(degree) +
                                        ", " + std::to_string(count) + " points, point " + std::to_string(k));
                            ++compared;
                        }
                    }
            }
            // certify with what evaluate_ahead found, against certify alone under the reference.
            const uint64_t proved = certify_ahead(tier, points, goal(53), 85, compared);
            require(proved >= input.roots.size(), "the tier proved fewer points than the polynomial has roots");
            // A coarse goal, where the smaller circles decide, and another precision.
            certify_ahead(tier, points, goal(10), 42, compared);
        }
    // A polynomial whose coefficients binary64 cannot hold at once (test_compensated_certificate.cpp, 3b): the sum of
    // (k + 1) (-2x)^k to degree 1,100. may_certify refuses the points on its roots' circle |x| = 1/2 and admits others;
    // the answers made ahead must be may_certify's own.
    {
        constexpr int n = 1100;
        d::ComplexPolynomial scaled;
        for (int k = 0; k <= n; ++k) {
            d::ComplexBall c;
            acb_set_si(c.value(), (k % 2 ? -1 : 1) * (k + 1));
            acb_mul_2exp_si(c.value(), c.value(), k);
            acb_poly_set_coeff_acb(scaled.value(), k, c.value());
        }
        const d::CompensatedCertificate wide(scaled.value());
        require(wide.available() && wide.carries_errors(), "the wide sum: unavailable, or its errors misread");
        std::vector<m::Complex<double>> points;
        for (int j = 0; j < 40; ++j) {
            const double angle = 6.283185307179586 * (j + 0.5) / 40;
            for (const double radius : {0.5, 0.25, 0.75, 1.5})
                points.push_back({radius * std::cos(angle), radius * std::sin(angle)});
        }
        uint64_t refused = 0;
        const auto o = goal(53);
        for (const auto variant : variants()) {
            const Under scope(variant);
            const auto count = static_cast<slong>(points.size());
            d::RootVector balls(count);
            std::vector<slong> indices;
            for (slong i = 0; i < count; ++i) {
                acb_set_d_d(balls.data() + i, points[static_cast<size_t>(i)].re, points[static_cast<size_t>(i)].im);
                indices.push_back(i);
            }
            std::vector<d::CompensatedAhead> ahead(points.size());
            wide.evaluate_ahead(balls.data(), indices, o, 85, ahead);
            for (size_t i = 0; i < points.size(); ++i) {
                const bool admitted = wide.may_certify(points[i].re, points[i].im, o);
                require((ahead[i].state == d::CompensatedAhead::State::Refused) == !admitted,
                        "evaluate_ahead's admission differs from may_certify's");
                refused += !admitted;
                ++compared;
            }
        }
        require(refused >= 40, "the wide sum: the points on its roots' circle were not refused");
        certify_ahead(wide, points, o, 85, compared);
    }
    return compared;
}

// 4.
bool same(const m::Jet<double> &a, const m::horner_kernel::Sums &b) {
    return same(a.value.re, b.value_re) && same(a.value.im, b.value_im) && same(a.derivative.re, b.derivative_re) &&
           same(a.derivative.im, b.derivative_im) && same(a.absolute_sum, b.absolute_sum);
}
constexpr int every_flag = FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW | FE_INEXACT;
std::vector<m::Complex<double>> random_coefficients(std::mt19937_64 &rng, size_t count, int kind) {
    std::uniform_real_distribution<double> unit(-1, 1);
    std::uniform_int_distribution<int> exponent(-300, 300);
    std::vector<m::Complex<double>> c(count);
    for (auto &a : c) {
        const double scale = kind == 1 ? std::ldexp(1.0, exponent(rng)) : kind == 2 ? 0x1p-1060 : 1;
        a = {unit(rng) * scale, kind == 3 ? 0 : unit(rng) * scale};
        if (kind == 4 && (rng() & 1))
            a = {};
    }
    if (c.back() == m::Complex<double>{})
        c.back() = {0.5, 0.25};
    return c;
}
uint64_t horner_loop() {
#ifdef CERTPOLYSOLVE_HAVE_AVX2_FMA_KERNELS
    if (!d::kernel_variant_available(d::KernelVariant::Avx2Fma))
        return 0;
    std::mt19937_64 rng(2026111703);
    uint64_t compared = 0;
    constexpr size_t group = m::horner_kernel::group;
    for (const size_t count : {1, 2, 3, 4, 5, 9, 17, 40, 257})
        for (int kind = 0; kind <= 3; ++kind) { // every kind is dense: PreparedJets' dense loop is the reference
            const auto coefficients = random_coefficients(rng, count, kind);
            const m::PreparedJets<double, 1> jets(coefficients);
            require(jets.dense(), "the reference would not take its dense loop");
            const auto points = sample_points(rng, 2);
            for (size_t first = 0; first < points.size(); ++first) {
                double zr[group], zi[group], radius[group];
                m::Complex<double> z[group];
                for (size_t l = 0; l < group; ++l) {
                    z[l] = points[(first + l * 3) % points.size()];
                    zr[l] = z[l].re, zi[l] = z[l].im;
                }
                // The norms first, outside both windows: the kernel is given them.
                for (size_t l = 0; l < group; ++l)
                    radius[l] = m::norm(z[l]);
                m::Jet<double> expected[group];
                std::feclearexcept(FE_ALL_EXCEPT);
                for (size_t l = 0; l < group; ++l)
                    expected[l] = jets.evaluate(z[l]);
                const int raised = std::fetestexcept(every_flag);
                // evaluate computes the norm again inside its window; its flags are taken apart.
                std::feclearexcept(FE_ALL_EXCEPT);
                for (size_t l = 0; l < group; ++l)
                    radius[l] = m::norm(z[l]);
                const int norm_flags = std::fetestexcept(every_flag);
                m::horner_kernel::Sums wide[group];
                std::feclearexcept(FE_ALL_EXCEPT);
                m::horner_kernel::avx2(coefficients.data(), jets.magnitudes().data(), count, zr, zi, radius, wide);
                const int kernel_flags = std::fetestexcept(every_flag);
                require((kernel_flags | norm_flags) == raised,
                        "the dense loop's four-point variant raises other floating-point flags: " +
                            std::to_string(kernel_flags | norm_flags) + " against " + std::to_string(raised));
                for (size_t l = 0; l < group; ++l) {
                    require(same(expected[l], wide[l]),
                            "the dense loop's four-point variant differs from the reference: " +
                                std::to_string(count) + " coefficients of kind " + std::to_string(kind) + ", point " +
                                std::to_string(first) + ", lane " + std::to_string(l));
                    ++compared;
                }
            }
        }
    return compared;
#else
    return 0;
#endif
}

// 5.
bool same(const m::NewtonEvaluation<double> &a, const m::NewtonEvaluation<double> &b) {
    return same(a.value, b.value) && same(a.denominator, b.denominator) &&
           same(a.coordinate_factor, b.coordinate_factor) && same(a.absolute_sum, b.absolute_sum) &&
           a.terms == b.terms;
}
uint64_t newton_ranges() {
    std::mt19937_64 rng(2026111704);
    uint64_t compared = 0;
    for (const size_t count : {2, 3, 6, 21, 130})
        for (int kind = 0; kind <= 6; ++kind) {
            // Kinds 0 to 3 as above; 4 with zeros among the coefficients; 5 with roots at zero (the reversed
            // polynomial is shorter); 6 sparse (three terms: neither side takes a dense loop).
            auto coefficients = random_coefficients(rng, count, std::min(kind, 4));
            if (kind == 5 && count > 3)
                coefficients[0] = coefficients[1] = {};
            if (kind == 6)
                for (size_t k = 1; k + 1 < coefficients.size(); ++k)
                    if (k != coefficients.size() / 2)
                        coefficients[k] = {};
            const m::PreparedNewton<double> newton(coefficients);
            const auto points = sample_points(rng, 3);
            std::vector<m::NewtonEvaluation<double>> expected(points.size());
            std::feclearexcept(FE_ALL_EXCEPT);
            for (size_t i = 0; i < points.size(); ++i)
                expected[i] = newton.evaluate(points[i]);
            const int raised = std::fetestexcept(every_flag);
            const m::NewtonEvaluation<double> untouched{{7, 7}, {7, 7}, {7, 7}, 7, 7};
            for (const auto variant : variants()) {
                const Under scope(variant);
                for (const size_t period : {1, 2, 3, 5}) { // every point; then every second, third, fifth left out
                    std::vector<m::NewtonEvaluation<double>> out(points.size(), untouched);
                    const auto active = [&](size_t j) { return period == 1 || j % period != 0; };
                    for (size_t first = 0; first < points.size(); first += 11) {
                        const size_t last = std::min(points.size(), first + 11);
                        if (period == 1 && first == 0)
                            std::feclearexcept(FE_ALL_EXCEPT);
                        newton.evaluate_range(points, first, last, active, out);
                    }
                    if (period == 1)
                        require(std::fetestexcept(every_flag) == raised,
                                "evaluate_range raises other floating-point flags than evaluate a point");
                    for (size_t i = 0; i < points.size(); ++i) {
                        require(same(out[i], active(i) ? expected[i] : untouched),
                                "evaluate_range differs from evaluate: " + std::to_string(count) +
                                    " coefficients of kind " + std::to_string(kind) + ", point " + std::to_string(i));
                        ++compared;
                    }
                }
            }
        }
    return compared;
}

// 6.
struct Packet {
    d::EngineAdvance progress;
    d::CandidateStatistics statistics{};
    std::vector<m::Complex<double>> points;
};
Packet packet(const acb_poly_t polynomial, slong degree, d::KernelVariant variant) {
    const Under scope(variant);
    d::RootVector seeds(degree), output(degree);
    require(d::initial_candidates(seeds.data(), polynomial, 64), "no starting points");
    // Underflow allowed and the non-dyadic scale, as in the certificate attempt's packet.
    auto engine = d::prepare_machine_engine(d::IterationMethod::Aberth, d::ArithmeticTier::Binary64, polynomial,
                                            seeds.view(), false, false, true, UINT64_MAX, true, true);
    require(bool(engine), "binary64 preparation failed");
    Packet out;
    cps_options options;
    out.progress = engine->advance(d::heuristics::machine_cover_sweeps, options, out.statistics);
    engine->export_candidates(output.data());
    for (slong i = 0; i < degree; ++i)
        out.points.push_back({arf_get_d(arb_midref(acb_realref(output.data() + i)), ARF_RND_NEAR),
                              arf_get_d(arb_midref(acb_imagref(output.data() + i)), ARF_RND_NEAR)});
    return out;
}
uint64_t engine_packets() {
    std::mt19937_64 rng(2026111705);
    std::uniform_real_distribution<double> unit(-1, 1);
    uint64_t compared = 0;
    for (const slong degree : {5, 8, 21, 64, 150, 301})
        for (int kind = 0; kind <= 3; ++kind) {
            // 0: real coefficients in [-1, 1]; 1: complex; 2: real, the roots pushed outside the unit circle
            // (coefficient k times 2^-k: p(x / 2)); 3: real with a wide spread of sizes (a harder packet, with
            // points that stop at the rounding of p).
            d::ComplexPolynomial polynomial;
            for (slong k = 0; k <= degree; ++k) {
                double re = unit(rng), im = kind == 1 ? unit(rng) : 0;
                if (kind == 2)
                    re = std::ldexp(re, -static_cast<int>(k));
                if (kind == 3)
                    re = std::ldexp(re, static_cast<int>((k * 37) % 61) - 30);
                d::ComplexBall c;
                acb_set_d_d(c.value(), re, im);
                acb_poly_set_coeff_acb(polynomial.value(), k, c.value());
            }
            d::ComplexBall lead;
            acb_set_d_d(lead.value(), 1, 0);
            acb_poly_set_coeff_acb(polynomial.value(), degree, lead.value());
            const Packet expected = packet(polynomial.value(), degree, d::KernelVariant::Reference);
            for (const auto variant : variants()) {
                const Packet got = packet(polynomial.value(), degree, variant);
                const std::string where = "degree " + std::to_string(degree) + ", kind " + std::to_string(kind);
                require(got.points.size() == expected.points.size(), "a packet lost points: " + where);
                for (size_t i = 0; i < got.points.size(); ++i)
                    require(same(got.points[i], expected.points[i]), "a packet ends on other points: " + where);
                require(got.progress.arithmetic_failure == expected.progress.arithmetic_failure &&
                            got.progress.arithmetic_limited == expected.progress.arithmetic_limited &&
                            got.progress.exhausted == expected.progress.exhausted &&
                            got.progress.active == expected.progress.active &&
                            got.progress.work == expected.progress.work &&
                            got.progress.candidates.converged == expected.progress.candidates.converged &&
                            got.progress.candidates.finite == expected.progress.candidates.finite &&
                            same(got.progress.candidates.correction_log2, expected.progress.candidates.correction_log2),
                        "a packet's progress differs: " + where);
                require(got.statistics.candidate_sweeps == expected.statistics.candidate_sweeps &&
                            got.statistics.evaluations == expected.statistics.evaluations &&
                            got.statistics.evaluation_terms == expected.statistics.evaluation_terms &&
                            got.statistics.root_interactions == expected.statistics.root_interactions,
                        "a packet's work counters differ: " + where);
                compared += got.points.size();
            }
        }
    return compared;
}

// 7.
struct DirectSum {
    m::Complex<double> sum;
    uint64_t terms = 0;
    bool finite = true; // false at a point equal to z
};
// The sweep's direct sum at point i (machine_candidates.cpp), a term at a time.
DirectSum direct_reference(const std::vector<m::Complex<double>> &points, size_t i) {
    DirectSum out;
    const m::Complex<double> z = points[i];
    for (size_t j = 0; j < points.size(); ++j)
        if (j != i) {
            const m::Complex<double> difference = z - points[j];
            if (difference == m::Complex<double>{}) {
                out.finite = false;
                break;
            }
            out.sum += m::interaction_term(difference);
            ++out.terms;
        }
    return out;
}
#ifdef CERTPOLYSOLVE_HAVE_AVX2_FMA_KERNELS
// The same with the kernel, as the sweep calls it: the neighbours below i, then those above.
DirectSum direct_grouped(const std::vector<m::Complex<double>> &points, size_t i) {
    DirectSum out;
    const m::Complex<double> z = points[i];
    const auto add_terms = [&](size_t first, size_t last) {
        size_t j = m::interaction_kernel::add_groups(points.data(), first, last, z.re, z.im, out.sum.re, out.sum.im);
        out.terms += j - first;
        for (; j < last; ++j) {
            const m::Complex<double> difference = z - points[j];
            if (difference == m::Complex<double>{}) {
                out.finite = false;
                return false;
            }
            out.sum += m::interaction_term(difference);
            ++out.terms;
        }
        return true;
    };
    if (add_terms(0, i))
        add_terms(i + 1, points.size());
    return out;
}
#endif
uint64_t interaction_sums() {
#ifdef CERTPOLYSOLVE_HAVE_AVX2_FMA_KERNELS
    if (!d::kernel_variant_available(d::KernelVariant::Avx2Fma))
        return 0;
    std::mt19937_64 rng(2026111706);
    std::uniform_real_distribution<double> unit(-1, 1);
    uint64_t compared = 0;
    const auto compare = [&](const std::vector<m::Complex<double>> &points, const std::string &where) {
        std::vector<DirectSum> expected(points.size());
        std::feclearexcept(FE_ALL_EXCEPT);
        for (size_t i = 0; i < points.size(); ++i)
            expected[i] = direct_reference(points, i);
        const int raised = std::fetestexcept(every_flag);
        std::feclearexcept(FE_ALL_EXCEPT);
        for (size_t i = 0; i < points.size(); ++i) {
            const DirectSum got = direct_grouped(points, i);
            require(got.finite == expected[i].finite && got.terms == expected[i].terms &&
                        (!got.finite || same(got.sum, expected[i].sum)),
                    "the interaction sum's four-term variant differs from the reference: " + where + ", point " +
                        std::to_string(i));
            ++compared;
        }
        require(std::fetestexcept(every_flag) == raised,
                "the interaction sum's four-term variant raises other floating-point flags: " + where);
    };
    for (size_t n = 2; n <= 70; n = n < 13 ? n + 1 : n + 19) {
        for (int kind = 0; kind <= 4; ++kind) {
            std::vector<m::Complex<double>> points(n);
            for (auto &p : points) {
                const double r = 1 + 0.3 * unit(rng), angle = 3.141592653589793 * unit(rng);
                p = {r * std::cos(angle), r * std::sin(angle)};
                if (kind == 1) // a lattice: differences on the axes and of equal parts, zeros of both signs
                    p = {std::round(4 * p.re) / 4, std::round(4 * p.im) / 4};
                if (kind == 2) // tiny differences: the quotients overflow
                    p = {1 + p.re * 0x1p-1000, p.im * 0x1p-1000};
                if (kind == 3) // huge ones: the quotients underflow
                    p = {p.re * 0x1p1000, p.im * 0x1p1000};
                if (kind == 4) // all real
                    p.im = 0;
            }
            if (kind == 1) { // the lattice repeats points: keep them distinct, a coincident pair is its own case below
                for (size_t a = 0; a < n; ++a)
                    points[a] = {points[a].re + static_cast<double>(a) * 0x1p-20, points[a].im};
            }
            compare(points, std::to_string(n) + " points of kind " + std::to_string(kind));
            // A coincident pair at every place: the sum fails at both points of the pair, after the same terms.
            if (kind == 0 && n <= 13)
                for (size_t a = 0; a < n; ++a)
                    for (size_t b = a + 1; b < n; ++b) {
                        auto with_pair = points;
                        with_pair[b] = with_pair[a];
                        compare(with_pair, std::to_string(n) + " points, " + std::to_string(a) + " equal to " +
                                               std::to_string(b));
                    }
        }
    }
    return compared;
#else
    return 0;
#endif
}

// 8.
uint64_t majorant_sums() {
#ifdef CERTPOLYSOLVE_HAVE_AVX2_FMA_KERNELS
    if (!d::kernel_variant_available(d::KernelVariant::Avx2Fma))
        return 0;
    std::mt19937_64 rng(2026111707);
    std::uniform_real_distribution<double> unit(0, 1);
    std::uniform_int_distribution<int> exponent(-300, 300);
    uint64_t compared = 0;
    constexpr size_t group = d::compensated_kernel::majorant_group;
    for (const uint64_t top : {0, 1, 2, 3, 7, 40, 300})
        for (int kind = 0; kind <= 3; ++kind) {
            std::vector<double> coefficients(top + 1);
            for (double &c : coefficients) {
                c = unit(rng);
                if (kind == 1)
                    c = std::ldexp(c, exponent(rng));
                if (kind == 2)
                    c *= 0x1p-1060; // subnormal
                if (kind == 3 && (rng() & 1))
                    c = 0;
            }
            for (uint64_t bottom = 0; bottom <= std::min<uint64_t>(top, 2); ++bottom)
                for (int round = 0; round < 3; ++round) {
                    double x[group], wide[group];
                    for (size_t l = 0; l < group; ++l) {
                        // Sixteen different points: below 1, about 1, above it (the sums overflow), tiny, zero.
                        const double choices[] = {unit(rng), 1 + 0x1p-30 * unit(rng), 1.5 + unit(rng), 1e3 * unit(rng),
                                                  0x1p-600 * unit(rng), 0, 1, 1e300};
                        x[l] = choices[(l + static_cast<size_t>(round) * 3) % 8];
                    }
                    d::compensated_kernel::majorants_avx2(coefficients.data(), top, bottom, x, wide);
                    for (size_t l = 0; l < group; ++l) {
                        require(same(wide[l], d::compensated_kernel::majorant(coefficients.data(), top, bottom, x[l])),
                                "the majorants' sixteen-point variant differs from the reference: top " +
                                    std::to_string(top) + ", bottom " + std::to_string(bottom) + ", kind " +
                                    std::to_string(kind) + ", lane " + std::to_string(l));
                        ++compared;
                    }
                }
        }
    return compared;
#else
    return 0;
#endif
}
} // namespace

int main() {
    try {
        choice();
        const bool wide = d::kernel_variant_available(d::KernelVariant::Avx2Fma);
        const uint64_t loop = compensated_loop();
        const uint64_t evaluations = certificate_evaluations();
        require(proved_on_a_later_try > 0, "no point was proved by a later try: the Newton retry was not exercised");
        const uint64_t horner = horner_loop();
        const uint64_t ranges = newton_ranges();
        const uint64_t packets = engine_packets();
        const uint64_t sums = interaction_sums();
        const uint64_t majorants = majorant_sums();
        require(!wide || (loop > 0 && horner > 0 && sums > 0 && majorants > 0),
                "the machine runs AVX2 and FMA, and their kernels were not compared");
        std::cout << "kernel variants: " << (wide ? "AVX2 and FMA compared with the reference" :
                                                    "this machine runs the reference only")
                  << "; compensated loop " << loop << " points, certificate evaluations " << evaluations << " ("
                  << proved_on_a_later_try << " proofs on a later try)"
                  << ", dense loop " << horner << " points, ranges " << ranges << ", packet points " << packets
                  << ", interaction sums " << sums << ", majorant sums " << majorants << "\n";
    } catch (const std::exception &e) {
        std::cerr << "test_kernel_variants: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
