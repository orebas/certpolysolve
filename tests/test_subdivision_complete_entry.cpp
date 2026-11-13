#include "reference_math.hpp"
#include "internal.hpp"
#include "methods/subdivision/subdivision.hpp"
#include <certpolysolve/certpolysolve.hpp>
#include <array>
#include <iostream>
#include <limits>

namespace d = certpolysolve::detail;
using reference::Integer;
using reference::Rational;
using reference::require;

namespace {
struct Point { Rational real, imag; };
constexpr slong fixture_bits = 2048;

std::vector<Point> ring() { return {{1, 0}, {0, 1}, {-1, 0}, {0, -1}}; }

d::ComplexPolynomial polynomial(const std::vector<Point> &points) {
    // Exact rational linear-factor expansion, independent of the solver.
    std::vector<Point> coefficients{{1, 0}};
    for (const auto &point : points) {
        std::vector<Point> next(coefficients.size() + 1, {0, 0});
        for (size_t i = 0; i < coefficients.size(); ++i) {
            const auto &a = coefficients[i];
            next[i].real -= a.real * point.real - a.imag * point.imag;
            next[i].imag -= a.real * point.imag + a.imag * point.real;
            next[i + 1].real += a.real;
            next[i + 1].imag += a.imag;
        }
        coefficients = std::move(next);
    }
    d::ComplexPolynomial result;
    for (size_t i = 0; i < coefficients.size(); ++i) {
        d::ComplexBall coefficient;
        d::Number(coefficients[i].real.str()).enclose(acb_realref(coefficient.value()), fixture_bits);
        d::Number(coefficients[i].imag.str()).enclose(acb_imagref(coefficient.value()), fixture_bits);
        require(acb_is_exact(coefficient.value()), "fixture coefficient was rounded");
        acb_poly_set_coeff_acb(result.value(), static_cast<slong>(i), coefficient.value());
    }
    return result;
}

std::vector<d::Root> point_cover(const std::vector<Point> &points, uint64_t multiplicity = 1) {
    std::vector<d::Root> cover;
    for (const auto &point : points) {
        d::Root root;
        d::Number(point.real.str()).enclose(acb_realref(root.box.value()), fixture_bits);
        d::Number(point.imag.str()).enclose(acb_imagref(root.box.value()), fixture_bits);
        require(acb_is_exact(root.box.value()), "fixture point was rounded");
        root.count = root.multiplicity = multiplicity;
        root.certificate = "exact-linear-factor";
        root.publish(fixture_bits);
        cover.push_back(std::move(root));
    }
    return cover;
}

d::Root global_bound(const d::ComplexPolynomial &p, uint64_t multiplicity = 1) {
    d::Root global;
    require(d::cauchy_cluster(global, p.value(), 64), "fixture needs a global bound");
    global.count *= multiplicity;
    global.multiplicity *= multiplicity;
    return global;
}

Rational endpoint(const d::Endpoint &x) {
    const auto text = x.text();
    return Rational(Integer(text.mantissa)) * reference::power2(std::stoi(text.exponent));
}

void check_cover(const std::vector<d::Root> &cover, const std::vector<Point> &points, uint64_t multiplicity = 1) {
    std::vector<unsigned> ownership(points.size(), 0);
    for (const auto &root : cover) {
        std::array<Rational, 4> box;
        for (size_t i = 0; i < box.size(); ++i) box[i] = endpoint(root.endpoints[i]);
        uint64_t contained = 0;
        for (size_t i = 0; i < points.size(); ++i)
            if (box[0] <= points[i].real && points[i].real <= box[1] &&
                box[2] <= points[i].imag && points[i].imag <= box[3]) {
                ++ownership[i];
                contained += multiplicity;
            }
        require(contained == root.count && contained != 0, "published count differs from exact planted roots");
    }
    for (const auto count : ownership) require(count == 1, "planted root ownership is not exactly once");
    require(d::disjoint(cover), "published regions overlap");
}

std::vector<d::Root> cover_of(const d::Subdivision &search) {
    std::vector<d::Root> cover;
    require(search.complete_cover(cover), "entry lost its complete certified cover");
    return cover;
}

void all_root_ring(uint64_t multiplicity) {
    const auto points = ring();
    const auto p = polynomial(points);
    const auto warm = point_cover(points, multiplicity);
    const auto global = global_bound(p, multiplicity);
    d::SearchStatistics statistics{};
    d::SubdivisionCompleteEntry entry;
    auto search = d::Subdivision::complete_start(global, warm, multiplicity, 32, statistics, &entry);
    require(search && search->complete_entry() && entry.used_warm_cover && !entry.exact_point,
            "complete warm ring did not enter by one normalized square");
    require(search->cells() == 1 && search->pending() == 1, "warm forest was imported instead of coalesced");
    require(search->min_bits() == 64 && search->max_bits() == 64, "entry inherited a high working precision");
    auto cover = cover_of(*search);
    require(cover.size() == 1 && cover[0].count == 4 * multiplicity, "all-root count was lost");
    check_cover(cover, points, multiplicity);
    std::vector<d::SubdivisionPrecisionState> state;
    search->precision_state(state);
    require(state.size() == 1 && state[0].newton_bits == 2 && state[0].trial_bits == 0,
            "complete entry retained a warm Newton speed or temporary precision");
    // Native entry must have a genuine root-free annulus, not merely a count
    // assigned to an imported support. Every exact root is strictly in R/2.
    const Rational radius = endpoint(cover[0].endpoints[1]); // this ring's center is zero
    for (const auto &point : points)
        require(point.real * point.real + point.imag * point.imag < radius * radius / 4,
                "initial all-root annulus lacks its containment margin");
    d::Subdivision generic(warm, multiplicity, 64, statistics, d::SubdivisionPrecision::TrialLocal);
    require(!generic.complete_entry(), "generic warm constructor claimed complete-route provenance");
    cps_options options;
    options.goal = CPS_ROOT_ISOLATION;
    options.max_work_bits = 1024;
    const d::Encloser enclose = [&](acb_poly_t out, slong bits) { acb_poly_set_round(out, p.value(), bits); };
    constexpr unsigned turns = 3000; // finite test limit, not a solver threshold
    unsigned used = 0;
    while (search->pending() && used++ < turns) {
        require(search->step(options, enclose, 4096) == CPS_STOP_NONE, "complete route unexpectedly stopped");
        std::vector<d::Root> next;
        if (search->complete_cover(next)) {
            check_cover(next, points, multiplicity);
            cover = std::move(next);
        }
    }
    require(!search->pending() && search->complete_entry(), "normalized ring failed to isolate persistently");
    require(cover.size() == points.size(), "isolated ring did not publish every singleton");
    options.goal = CPS_ROOT_ACCURACY;
    options.bits = 100;
    search->reopen(options);
    require(search->complete_entry(), "reopen lost admissible owned-state provenance");
    check_cover(cover_of(*search), points, multiplicity);
}

void translated_and_padded_geometry() {
    const Rational center = reference::power2(180), gap = reference::power2(-128);
    const std::vector<Point> points{{center - gap, 0}, {center + gap, 0}};
    const auto p = polynomial(points);
    const auto warm = point_cover(points);
    d::SearchStatistics statistics{};
    d::SubdivisionCompleteEntry entry;
    auto search = d::Subdivision::complete_start(global_bound(p), warm, 1, 512, statistics, &entry);
    require(search && entry.used_warm_cover && !entry.padded_width, "translated fine localization was discarded");
    const auto cover = cover_of(*search);
    check_cover(cover, points);
    require(endpoint(cover[0].endpoints[1]) - endpoint(cover[0].endpoints[0]) < 1,
            "normalization padded width to a fraction of the enormous global bound");
    require(search->max_bits() == 64, "geometry precision became the initial search precision");

    const Rational tiny = reference::power2(-256);
    const std::vector<Point> tiny_points{{1 - tiny, 0}, {1 + tiny, 0}};
    const auto tiny_p = polynomial(tiny_points);
    search = d::Subdivision::complete_start(global_bound(tiny_p), point_cover(tiny_points), 1, 64,
                                            statistics, &entry);
    require(search && entry.used_warm_cover && entry.padded_width, "tiny entry did not enforce source-derived floor");
    check_cover(cover_of(*search), tiny_points);
    std::vector<d::SubdivisionPrecisionState> state;
    search->precision_state(state);
    d::RealBall floor;
    arb_one(floor.value());
    arb_mul_2exp_si(floor.value(), floor.value(), -64);
    require(arb_ge(state[0].certificate_radius.value(), floor.value()), "normalization width is below its declared floor");
}

void unsuitable_covers_and_zero() {
    const auto points = ring();
    const auto p = polynomial(points);
    const auto global = global_bound(p);
    d::SearchStatistics statistics{};
    d::SubdivisionCompleteEntry entry;
    auto native = d::Subdivision::complete_start(global, {}, 1, 32, statistics, &entry);
    require(native && native->complete_entry() && !entry.used_warm_cover, "native bound alternative unavailable");
    const auto native_cover = cover_of(*native);
    auto oversized = global;
    acb_mul_2exp_si(oversized.box.value(), oversized.box.value(), 1);
    oversized.publish(64); // still an all-roots cover, but no useful bounded localization
    for (auto unsuitable : {std::vector<d::Root>{point_cover(points)[0]},
                            point_cover(std::vector<Point>(4, points[0])), std::vector<d::Root>{oversized}}) {
        auto search = d::Subdivision::complete_start(global, unsuitable, 1, 32, statistics, &entry);
        require(search && !entry.used_warm_cover, "unsuitable warm geometry accepted");
        require(acb_equal(cover_of(*search)[0].box.value(), native_cover[0].box.value()),
                "rejected warm cover changed the native alternative");
        check_cover(cover_of(*search), points);
    }
    require(!d::Subdivision::complete_start(global, {}, 0, 32, statistics), "zero multiplicity accepted");
    require(!d::Subdivision::complete_start(global, {}, 1, 0, statistics), "understated source scale accepted");
    auto shifted = global;
    arb_add_ui(acb_realref(shifted.box.value()), acb_realref(shifted.box.value()), 1, 64);
    require(!d::Subdivision::complete_start(shifted, {}, 1, 32, statistics), "non-global bound accepted");

    const std::vector<Point> zero{{0, 0}};
    const auto linear = polynomial(zero);
    auto exact = d::Subdivision::complete_start(global_bound(linear, 3), point_cover(zero, 3), 3, 32,
                                                statistics, &entry);
    require(exact && entry.exact_point && entry.used_warm_cover, "certified exact singleton was widened");
    check_cover(cover_of(*exact), zero, 3);
    cps_options options;
    options.goal = CPS_ROOT_ISOLATION;
    const d::Encloser enclose = [&](acb_poly_t out, slong bits) { acb_poly_set_round(out, linear.value(), bits); };
    require(exact->step(options, enclose, 0) == CPS_STOP_NONE && !exact->pending() && !exact->cells(),
            "exact singleton did not retire without a search cell");
    require(exact->complete_entry(), "exact singleton lost entry provenance");
}
} // namespace

int main() {
    certpolysolve::ThreadScope lifetime;
    try {
        all_root_ring(1);
        all_root_ring(3);
        translated_and_padded_geometry();
        unsuitable_covers_and_zero();
        std::cout << "complete subdivision entry: containment, normalization, ownership and continuation passed\n";
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
