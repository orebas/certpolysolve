// Conjugate pairs in the certificate attempt (src/certify/attempts/conjugate_pairs.hpp, newton_cover.cpp).
//
// On a polynomial with real coefficients the attempt certifies one point of each conjugate pair of points and gives
// the other the reflection of the region it proved; a reflection stands only in a cover the attempt completes. Which
// points are a pair is a guess from their places. This file checks the guess as a function, and then the attempt:
//
// 1. conjugate_partners on point sets: exact mirrors, coordinates on every kind of round number (the grid of an
//    earlier draft missed pairs with a coordinate on a cell's edge: 1, -2, 1/2), a lattice, a vertical and a
//    horizontal line of points, mirrors that are off by less and by more than the tolerance, a point with two
//    candidates (no pair on either side), points on the axis and points that are not finite; that the answer does
//    not depend on the order of the points; that the sorted search gives what comparing every two points gives; and
//    the scan limit (a lattice with more points a row than the limit is left unpaired, never mispaired).
// 2. machine_newton_cover on polynomials with planted roots, coordinates dyadic so that the coefficients are exact:
//    every region holds exactly one planted root, every pair gives one reflected region, each reflected region is
//    the exact mirror image of a region that was proved, and the kernel variants agree.
// 3. Nothing proved depends on the guess: with one deliberately wrong pair (EnginePolicy::conjugate_pair_fault) the
//    cover with the reflections is not disjoint, every reflection is replaced by the point's own certificate, and
//    the attempt certifies the cover it certifies without pairs, region for region.
// 4. A sparse source, where the ball certificate works alone: the reflected cover against the cover with the
//    reflection off, region for region.
// 5. Where the attempt leaves a point (a double real root; a goal beyond binary64 with one), no reflection is left
//    standing: what it hands back is the partial cover and the points it hands back without pairs, bit for bit.
#include "certify/attempts/conjugate_pairs.hpp"
#include "certify/attempts/newton_cover.hpp"
#include "core/cover.hpp"
#include "core/engine_policy.hpp"
#include "core/engine_statistics.hpp"
#include "core/heuristics.hpp"
#include "core/query.hpp" // cps_options
#include "core/raii.hpp"
#include "numeric/evaluate/structured_source.hpp"
#include "platform/kernel_variant.hpp"
#include <certpolysolve/certpolysolve.hpp>
#include <flint/acb_poly.h>
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>
#include <iostream>
#include <memory>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace d = certpolysolve::detail;

namespace {
void require(bool condition, const std::string &message) {
    if (!condition)
        throw std::runtime_error(message);
}

const int bits = d::heuristics::conjugate_pair_bits;
const std::size_t limit = d::heuristics::conjugate_scan_limit;
const std::string reflected_mark = "reflected in the real axis";

using Points = std::vector<d::PlanePoint>;
using Pairs = std::set<std::pair<std::size_t, std::size_t>>; // (below, above)

Pairs pairs_of(const Points &points, int tolerance_bits = bits, std::size_t scan_limit = limit) {
    const auto partner = d::conjugate_partners(points, tolerance_bits, scan_limit);
    require(partner.size() == points.size(), "the partner list has another length");
    Pairs out;
    for (std::size_t i = 0; i < points.size(); ++i) {
        if (partner[i] < 0)
            continue;
        const auto j = static_cast<std::size_t>(partner[i]);
        require(j < points.size() && partner[j] == static_cast<std::ptrdiff_t>(i), "a partner is not mutual");
        require((points[i].im > 0) != (points[j].im > 0) && points[i].im != 0 && points[j].im != 0,
                "partners on one side of the axis");
        if (points[i].im < 0)
            out.insert({i, j});
    }
    return out;
}

// The definition, by comparing every point below with every point above.
Pairs pairs_by_definition(const Points &points, int tolerance_bits) {
    const auto mirrors = [&](std::size_t a, std::size_t b) {
        const d::PlanePoint &p = points[a], &q = points[b];
        const double size = std::max(std::max(std::abs(p.re), std::abs(p.im)), std::max(std::abs(q.re), std::abs(q.im)));
        const double tolerance = std::ldexp(size, -tolerance_bits);
        return std::abs(p.re - q.re) <= tolerance && std::abs(std::abs(p.im) - std::abs(q.im)) <= tolerance;
    };
    std::vector<std::vector<std::size_t>> candidates(points.size());
    for (std::size_t a = 0; a < points.size(); ++a)
        for (std::size_t b = 0; b < points.size(); ++b)
            if (std::isfinite(points[a].re) && std::isfinite(points[a].im) && std::isfinite(points[b].re) &&
                std::isfinite(points[b].im) && points[a].im < 0 && points[b].im > 0 && mirrors(a, b)) {
                candidates[a].push_back(b);
                candidates[b].push_back(a);
            }
    Pairs out;
    for (std::size_t a = 0; a < points.size(); ++a)
        if (points[a].im < 0 && candidates[a].size() == 1 && candidates[candidates[a][0]].size() == 1)
            out.insert({a, candidates[a][0]});
    return out;
}

void pairing_tests() {
    // Exact mirrors on round numbers, every one on an edge of a power-of-two grid; each pair found.
    {
        Points points;
        for (const double x : {1.0, -2.0, 0.5, 0.0, 3.0, -0.75, 1024.0, 0x1p-3})
            for (const double y : {1.0, 0.5, 2.0, 0x1p-30, 7.0}) {
                points.push_back({x, y});
                points.push_back({x, -y});
            }
        require(pairs_of(points).size() == points.size() / 2, "round numbers: a pair was missed");
        require(pairs_of(points) == pairs_by_definition(points, bits), "round numbers: not the definition's pairs");
    }
    // Mirrors off by 2^-40 of their size either side of round numbers: found. Off by 2^-20: not pairs.
    {
        Points close, far;
        for (const double x : {1.0, -2.0, 0.5, 3.0}) {
            close.push_back({x * (1 + 0x1p-40), x == 0.5 ? 1.0 : 2.0 * (1 - 0x1p-41)});
            close.push_back({x * (1 - 0x1p-40), -(x == 0.5 ? 1.0 : 2.0) * (1 + 0x1p-41)});
            far.push_back({x * (1 + 0x1p-20), 2.0});
            far.push_back({x, -2.0});
        }
        require(pairs_of(close).size() == 4, "mirrors 2^-40 apart: a pair was missed");
        require(pairs_of(far).empty(), "mirrors 2^-20 apart were paired");
    }
    // The tolerance's edge: a mismatch of exactly 2^-bits of the size is a pair, the next number up is not.
    {
        const double size = 1.5, step = std::ldexp(size, -bits);
        require(pairs_of({{size, 1.0}, {size - step, -1.0}}).size() == 1, "a mismatch of the tolerance is no pair");
        require(pairs_of({{size, 1.0}, {std::nextafter(size - step, 0.0), -1.0}}).empty(),
                "a mismatch above the tolerance is a pair");
    }
    // Two candidates: neither they nor the point are paired, whichever side has the two.
    {
        const double nudge = 0x1p-30;
        require(pairs_of({{1.0, 1.0}, {1.0 + nudge, 1.0}, {1.0, -1.0}}).empty(), "two above, one below: paired");
        require(pairs_of({{1.0, 1.0}, {1.0, -1.0}, {1.0 + nudge, -1.0}}).empty(), "one above, two below: paired");
        // A third pair, clear of them, is found beside the cluster.
        require(pairs_of({{1.0, 1.0}, {1.0 + nudge, 1.0}, {1.0, -1.0}, {1.0 + nudge, -1.0}, {3.0, 1.0}, {3.0, -1.0}})
                        .size() == 1,
                "a pair beside a cluster was missed, or the cluster was paired");
    }
    // The axis, zero, and what is not finite: no partner.
    {
        const double inf = std::numeric_limits<double>::infinity(), nan = std::numeric_limits<double>::quiet_NaN();
        const Points points = {{1.0, 0.0}, {1.0, -0.0}, {0.0, 0.0}, {inf, 1.0}, {inf, -1.0}, {2.0, nan},
                               {2.0, -1.0}, {nan, 1.0}, {5.0, 1.0}, {5.0, -1.0}};
        const Pairs pairs = pairs_of(points);
        require(pairs.size() == 1 && pairs.count({9, 8}) == 1, "the axis or a point that is not finite was paired");
    }
    // A lattice (many points share each coordinate), a vertical and a horizontal line, and random points with
    // mirrors off by random amounts around the tolerance: the definition's pairs, in any order of the points, by the
    // comparison of every two (few points) and by the sorted search (many).
    std::mt19937_64 random(20261117);
    const auto check = [&](Points points, const std::string &name, std::size_t expected) {
        const Pairs want = pairs_by_definition(points, bits);
        if (expected != std::size_t(-1))
            require(want.size() == expected, name + ": the definition gives another number of pairs");
        require(pairs_of(points) == want, name + ": not the definition's pairs");
        // Another order: the same pairs, as pairs of places.
        std::vector<std::size_t> order(points.size());
        for (std::size_t i = 0; i < order.size(); ++i)
            order[i] = i;
        std::shuffle(order.begin(), order.end(), random);
        Points shuffled(points.size());
        for (std::size_t i = 0; i < order.size(); ++i)
            shuffled[i] = points[order[i]];
        Pairs back;
        for (const auto &[below, above] : pairs_of(shuffled))
            back.insert({order[below], order[above]});
        require(back == want, name + ": the pairs depend on the order of the points");
    };
    for (const int side : {3, 12, 40}) {
        Points lattice, vertical, horizontal;
        for (int a = 0; a < side; ++a) {
            for (int b = 1; b <= side; ++b) {
                lattice.push_back({a - side / 2.0, static_cast<double>(b)});
                lattice.push_back({a - side / 2.0, -static_cast<double>(b)});
            }
            vertical.push_back({0.25, a + 1.0});
            vertical.push_back({0.25, -(a + 1.0)});
            horizontal.push_back({a - side / 2.0, 3.0});
            horizontal.push_back({a - side / 2.0, -3.0});
        }
        check(lattice, "a lattice of side " + std::to_string(side), static_cast<std::size_t>(side) * side);
        check(vertical, "a vertical line of " + std::to_string(side), static_cast<std::size_t>(side));
        check(horizontal, "a horizontal line of " + std::to_string(side), static_cast<std::size_t>(side));
    }
    for (const int count : {5, 30, 700}) {
        Points points;
        std::uniform_real_distribution<double> place(-4, 4), exponent(-45, -18);
        for (int k = 0; k < count; ++k) {
            const double x = place(random), y = std::abs(place(random)) + 0x1p-10;
            points.push_back({x, y});
            // Its mirror, off by a random power of two between 2^-45 and 2^-18 of its size; some points get none
            // and some get two.
            const int copies = k % 11 == 0 ? 0 : k % 13 == 0 ? 2 : 1;
            for (int c = 0; c < copies; ++c)
                points.push_back({x * (1 + std::exp2(exponent(random))), -y * (1 - std::exp2(exponent(random)))});
        }
        check(points, std::to_string(count) + " random points and their mirrors", std::size_t(-1));
    }
    // The scan limit: a lattice with more points a row and a column than the limit is searched nowhere and pairs
    // nothing; with the limit above its side every pair is found. (With the library's limit that is a lattice of
    // over 256 points a side.)
    {
        Points lattice;
        for (int a = 0; a < 20; ++a)
            for (int b = 1; b <= 20; ++b) {
                lattice.push_back({static_cast<double>(a), static_cast<double>(b)});
                lattice.push_back({static_cast<double>(a), -static_cast<double>(b)});
            }
        require(pairs_of(lattice, bits, 8).empty(), "a lattice beyond the scan limit was paired");
        require(pairs_of(lattice, bits, 32).size() == 400, "a lattice within the scan limit: a pair was missed");
    }
}

// The monic polynomial with the given roots, exactly: the coordinates are dyadic and the precision is far above
// what their products need.
d::ComplexPolynomial from_roots(const std::vector<std::complex<double>> &roots) {
    const slong n = static_cast<slong>(roots.size());
    d::RootVector points(n);
    for (slong i = 0; i < n; ++i)
        acb_set_d_d(points[i], roots[static_cast<size_t>(i)].real(), roots[static_cast<size_t>(i)].imag());
    d::ComplexPolynomial p;
    acb_poly_product_roots(p.value(), points.data(), n, 16384);
    for (slong k = 0; k <= n; ++k) {
        const acb_srcptr c = acb_poly_get_coeff_ptr(p.value(), k);
        require(acb_is_exact(c), "a test polynomial's coefficient is not exact");
    }
    return p;
}

bool real_coefficients(const d::ComplexPolynomial &p) {
    for (slong k = 0; k <= acb_poly_degree(p.value()); ++k)
        if (!arb_is_zero(acb_imagref(acb_poly_get_coeff_ptr(p.value(), k))))
            return false;
    return true;
}

struct Outcome {
    bool certified = false;
    std::vector<d::Root> cover;
    d::PartialNewtonCover partial;
    d::StudyStatistics counts;
};

struct Settings {
    bool reflection = true, fault = false, sparse = false;
    uint64_t bits = 53;
};

Outcome attempt(const d::ComplexPolynomial &polynomial, const Settings &settings) {
    const slong degree = acb_poly_degree(polynomial.value());
    d::PolynomialSourceOptions source_options;
    source_options.sparse_coefficients = settings.sparse;
    source_options.input = d::InputSemantics::Exact;
    d::PolynomialSource source(
        [&](acb_poly_t out, slong precision) {
            d::ComplexBall a;
            acb_poly_zero(out);
            for (slong k = 0; k <= degree; ++k) {
                acb_set_round(a.value(), acb_poly_get_coeff_ptr(polynomial.value(), k), precision);
                acb_poly_set_coeff_acb(out, k, a.value());
            }
        },
        degree, real_coefficients(polynomial), source_options);
    cps_options options;
    options.goal = CPS_ROOT_APPROXIMATION;
    options.bits = settings.bits;
    d::EnginePolicy policy;
    policy.study_statistics = std::make_shared<d::StudyStatistics>();
    policy.conjugate_reflection = settings.reflection;
    policy.conjugate_pair_fault = settings.fault;
    Outcome out;
    d::CandidateStatistics work{};
    out.certified = d::machine_newton_cover(out.cover, source, nullptr, options, 85, true, work, true, &policy, nullptr,
                                            nullptr, nullptr, true, &out.partial);
    out.counts = *policy.study_statistics;
    return out;
}

bool holds(const d::Root &region, const std::complex<double> &root) {
    return arf_cmp_d(region.endpoints[0].value(), root.real()) <= 0 &&
           arf_cmp_d(region.endpoints[1].value(), root.real()) >= 0 &&
           arf_cmp_d(region.endpoints[2].value(), root.imag()) <= 0 &&
           arf_cmp_d(region.endpoints[3].value(), root.imag()) >= 0;
}

bool is_reflected(const d::Root &region) { return region.certificate.find(reflected_mark) != std::string::npos; }

// Whether `image` is `region` mirrored in the real axis, end for end.
bool mirror_image(const d::Root &image, const d::Root &region) {
    d::Float low, high;
    arf_neg(low.value(), region.endpoints[3].value());
    arf_neg(high.value(), region.endpoints[2].value());
    return arf_equal(image.endpoints[0].value(), region.endpoints[0].value()) &&
           arf_equal(image.endpoints[1].value(), region.endpoints[1].value()) &&
           arf_equal(image.endpoints[2].value(), low.value()) && arf_equal(image.endpoints[3].value(), high.value());
}

// Every region holds exactly one planted root and every planted root lies in exactly one region.
void require_cover(const std::vector<d::Root> &cover, const std::vector<std::complex<double>> &roots,
                   const std::string &name) {
    require(cover.size() == roots.size(), name + ": the cover has another number of regions than roots");
    std::vector<int> held(roots.size(), 0);
    for (const d::Root &region : cover) {
        int inside = 0;
        for (size_t j = 0; j < roots.size(); ++j)
            if (holds(region, roots[j])) {
                ++inside;
                ++held[j];
            }
        require(inside == 1 && region.count == 1, name + ": a region does not hold exactly one planted root");
    }
    for (const int h : held)
        require(h == 1, name + ": a planted root is not in exactly one region");
}

// Whether two covers are the same, region for region: the ball, the four ends, the count and the certificate.
bool same_region(const d::Root &a, const d::Root &b) {
    for (int k = 0; k < 4; ++k)
        if (!arf_equal(a.endpoints[k].value(), b.endpoints[k].value()))
            return false;
    return acb_equal(a.box.value(), b.box.value()) && a.count == b.count && a.certificate == b.certificate;
}
bool same_cover(const std::vector<d::Root> &a, const std::vector<d::Root> &b) {
    if (a.size() != b.size())
        return false;
    for (size_t i = 0; i < a.size(); ++i)
        if (!same_region(a[i], b[i]))
            return false;
    return true;
}

// Each reflected region is the mirror image of exactly one region that is not reflected.
uint64_t require_mirrors(const std::vector<d::Root> &cover, const std::string &name) {
    uint64_t reflected = 0;
    for (const d::Root &region : cover) {
        if (!is_reflected(region))
            continue;
        ++reflected;
        int sources = 0;
        for (const d::Root &other : cover)
            sources += !is_reflected(other) && mirror_image(region, other);
        require(sources == 1, name + ": a reflected region is not the mirror image of exactly one proved region");
    }
    return reflected;
}

std::vector<std::complex<double>> with_conjugates(const std::vector<std::complex<double>> &upper,
                                                  const std::vector<double> &real) {
    std::vector<std::complex<double>> roots;
    for (const auto &z : upper) {
        roots.push_back(z);
        roots.push_back(std::conj(z));
    }
    for (const double x : real)
        roots.emplace_back(x, 0.0);
    return roots;
}

void attempt_tests() {
    struct Case {
        std::string name;
        std::vector<std::complex<double>> upper;
        std::vector<double> real;
    };
    std::vector<Case> cases;
    // Round coordinates: i, 1 + i, Gaussian integers and halves, with real roots beside them.
    cases.push_back(Case{"round coordinates",
                         {{0.0, 1.0}, {1.0, 1.0}, {-2.0, 1.0}, {0.5, 0.5}, {2.0, 3.0}, {-1.0, 2.0}, {1.0, 0.25}, {-0.5, 4.0}},
                         {3.0, -3.0, 0.75}});
    // A vertical line of roots, a horizontal line, and a lattice.
    {
        Case vertical{"a vertical line", {}, {-1.5}};
        Case horizontal{"a horizontal line", {}, {}};
        Case lattice{"a lattice", {}, {0.125}};
        for (int k = 1; k <= 12; ++k) {
            vertical.upper.emplace_back(0.25, 0.375 * k);
            horizontal.upper.emplace_back(0.375 * k - 2.5, 1.0);
        }
        for (int a = -2; a <= 2; ++a)
            for (int b = 1; b <= 4; ++b)
                lattice.upper.emplace_back(0.5 * a, 0.5 * b);
        cases.push_back(vertical);
        cases.push_back(horizontal);
        cases.push_back(lattice);
    }
    uint64_t reflected_in_all = 0;
    const d::KernelVariant before = d::kernel_variant();
    for (const Case &c : cases) {
        const auto roots = with_conjugates(c.upper, c.real);
        const d::ComplexPolynomial p = from_roots(roots);
        require(real_coefficients(p), c.name + ": the test polynomial is not real");
        std::vector<Outcome> by_variant;
        for (const auto variant : {d::KernelVariant::Reference, d::KernelVariant::Avx2Fma}) {
            if (!d::kernel_variant_available(variant))
                continue;
            d::set_kernel_variant(variant);
            Outcome outcome = attempt(p, {});
            require(outcome.certified, c.name + ": the attempt did not certify its cover");
            require_cover(outcome.cover, roots, c.name);
            const uint64_t reflected = require_mirrors(outcome.cover, c.name);
            require(reflected == c.upper.size(), c.name + ": " + std::to_string(reflected) + " regions reflected for " +
                                                     std::to_string(c.upper.size()) + " pairs");
            require(outcome.counts.machine_cover_reflected_certificates == reflected,
                    c.name + ": the counter of reflected certificates differs from the cover");
            reflected_in_all += reflected;
            by_variant.push_back(std::move(outcome));
        }
        d::set_kernel_variant(before);
        if (by_variant.size() == 2) {
            require(by_variant[0].cover.size() == by_variant[1].cover.size(), c.name + ": the variants differ");
            for (size_t i = 0; i < by_variant[0].cover.size(); ++i)
                require(acb_equal(by_variant[0].cover[i].box.value(), by_variant[1].cover[i].box.value()) &&
                            by_variant[0].cover[i].certificate == by_variant[1].cover[i].certificate,
                        c.name + ": the kernel variants give different regions");
        }
        // With the reflection off: the same roots, no region reflected.
        const Outcome plain = attempt(p, {.reflection = false});
        require(plain.certified, c.name + ": the attempt without the reflection did not certify");
        require_cover(plain.cover, roots, c.name + ", reflection off");
        require(require_mirrors(plain.cover, c.name) == 0 && plain.counts.machine_cover_reflected_certificates == 0,
                c.name + ": a region was reflected with the reflection off");
        // One wrong pair: the cover with the reflections is not disjoint, none is taken, every point is asked
        // itself, and the attempt certifies the cover it certifies without pairs.
        const Outcome spoiled = attempt(p, {.fault = true});
        require(spoiled.certified, c.name + ": with one wrong pair the attempt did not certify");
        require(spoiled.counts.machine_cover_reflection_refusals > 0,
                c.name + ": the wrong pair was never tried, so the test shows nothing");
        require(spoiled.counts.machine_cover_reflected_certificates == 0 && require_mirrors(spoiled.cover, c.name) == 0,
                c.name + ": a reflection was taken from a cover with a wrong pair");
        require_cover(spoiled.cover, roots, c.name + ", one wrong pair");
        require(same_cover(spoiled.cover, plain.cover), c.name + ": with one wrong pair the cover is not the plain one");
    }
    require(reflected_in_all > 0, "no region was reflected in any attempt");
    // A complex polynomial: nothing is reflected.
    {
        const std::vector<std::complex<double>> roots = {{0.0, 1.0}, {1.0, 1.0},  {-2.0, 1.0}, {0.5, 0.5},
                                                         {2.0, 3.0}, {-1.0, -2.0}, {1.0, -0.25}};
        const d::ComplexPolynomial p = from_roots(roots);
        require(!real_coefficients(p), "the complex test polynomial is real");
        const Outcome outcome = attempt(p, {});
        require(outcome.certified, "complex coefficients: not certified");
        require_cover(outcome.cover, roots, "complex coefficients");
        require(outcome.counts.machine_cover_reflected_certificates == 0, "complex coefficients: a region was reflected");
    }
    // A sparse real polynomial, x^40 - 3 x^7 + 2 x^3 - 5: the ball certificate alone. Reflected, each region below
    // the axis is the mirror image of a proved one; and the cover agrees with the one made without the reflection,
    // region for region (each region of one meets exactly one region of the other).
    {
        d::ComplexPolynomial p;
        acb_poly_set_coeff_si(p.value(), 40, 1);
        acb_poly_set_coeff_si(p.value(), 7, -3);
        acb_poly_set_coeff_si(p.value(), 3, 2);
        acb_poly_set_coeff_si(p.value(), 0, -5);
        const Outcome reflected = attempt(p, {.sparse = true}), plain = attempt(p, {.reflection = false, .sparse = true});
        require(reflected.certified && plain.certified, "the sparse polynomial: not certified");
        require(reflected.counts.machine_cover_compensated_certificates == 0,
                "the sparse polynomial was asked of the compensated tier: the ball certificate's pass is not exercised");
        const uint64_t count = require_mirrors(reflected.cover, "the sparse polynomial");
        uint64_t below = 0;
        for (const d::Root &region : reflected.cover)
            below += arb_is_negative(acb_imagref(region.box.value()));
        require(count > 0 && count == below, "the sparse polynomial: " + std::to_string(count) + " regions reflected of " +
                                                 std::to_string(below) + " below the axis");
        require(plain.counts.machine_cover_reflected_certificates == 0, "the sparse polynomial: reflected when off");
        require(reflected.cover.size() == 40 && plain.cover.size() == 40, "the sparse polynomial: not 40 regions");
        for (const d::Root &a : reflected.cover) {
            int meets = 0;
            for (const d::Root &b : plain.cover)
                meets += acb_overlaps(a.box.value(), b.box.value());
            require(meets == 1, "the sparse polynomial: a reflected cover's region does not meet exactly one region "
                                "of the cover made without the reflection");
        }
    }
}
} // namespace

// Where the attempt leaves a point, it runs as it does without pairs: no reflection, and what it certifies, what it
// hands back as a partial cover and where it leaves every point are the same with the pairs on and off, bit for bit.
void open_point_tests() {
    struct Case {
        std::string name;
        std::vector<double> real;
        uint64_t bits;
    };
    const std::vector<std::complex<double>> upper = {{0.0, 1.0}, {1.0, 1.0}, {-2.0, 1.0}, {0.5, 0.5}, {2.0, 3.0}};
    const std::vector<Case> cases = {
        {"a double real root beside pairs", {0.75, 0.75, -3.0}, 53},
        {"two real roots 2^-30 apart beside pairs", {0.75, 0.75 + 0x1p-30, -3.0}, 53},
        {"a double real root, a goal beyond binary64", {0.75, 0.75, -3.0}, 100},
        {"no close roots, a goal beyond binary64", {0.75, -3.0}, 100},
    };
    uint64_t open = 0, reflected_beyond = 0;
    for (const Case &c : cases) {
        const d::ComplexPolynomial p = from_roots(with_conjugates(upper, c.real));
        require(real_coefficients(p), c.name + ": the test polynomial is not real");
        const Outcome on = attempt(p, {.bits = c.bits}), off = attempt(p, {.reflection = false, .bits = c.bits});
        require(on.certified == off.certified, c.name + ": certified with the pairs on and not off, or the reverse");
        if (on.certified) {
            // A cover: the same one, unless the pairs completed it (then the proved regions are the plain cover's).
            if (on.counts.machine_cover_reflected_certificates == 0)
                require(same_cover(on.cover, off.cover), c.name + ": the cover differs from the one without pairs");
            else {
                reflected_beyond += c.bits > 53;
                require(require_mirrors(on.cover, c.name) == upper.size(), c.name + ": not every pair was reflected");
                for (size_t i = 0; i < on.cover.size(); ++i)
                    require(is_reflected(on.cover[i]) || same_region(on.cover[i], off.cover[i]),
                            c.name + ": a proved region differs from the one without pairs");
            }
            continue;
        }
        ++open;
        require(on.counts.machine_cover_reflected_certificates == 0, c.name + ": a reflection in an attempt that left a point");
        require(static_cast<bool>(on.partial.points) == static_cast<bool>(off.partial.points) &&
                    on.partial.certified == off.partial.certified && on.partial.pass_complete == off.partial.pass_complete &&
                    on.partial.certificate_bits == off.partial.certificate_bits,
                c.name + ": the partial cover differs from the one without pairs");
        for (size_t i = 0; i < on.partial.certified.size(); ++i) {
            if (on.partial.certified[i])
                require(same_region(on.partial.roots[i], off.partial.roots[i]),
                        c.name + ": a region of the partial cover differs from the one without pairs");
            if (on.partial.points)
                require(acb_equal((*on.partial.points)[static_cast<slong>(i)], (*off.partial.points)[static_cast<slong>(i)]),
                        c.name + ": a point is left elsewhere than without pairs");
        }
        require(on.counts.machine_cover_compensated_certificates == off.counts.machine_cover_compensated_certificates &&
                    on.counts.machine_cover_double_double_repairs == off.counts.machine_cover_double_double_repairs &&
                    on.counts.machine_cover_rounds == off.counts.machine_cover_rounds,
                c.name + ": the attempt's counts differ from those without pairs");
    }
    require(open >= 2, "fewer than two of the fixtures left a point: the test shows too little");
    require(reflected_beyond >= 1, "no fixture with a goal beyond binary64 was completed by reflections");
}

int main() {
    certpolysolve::ThreadScope lifetime;
    try {
        pairing_tests();
        attempt_tests();
        open_point_tests();
        std::cout << "conjugate pairs tests passed\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "FAILED: " << error.what() << '\n';
        return 1;
    }
}
