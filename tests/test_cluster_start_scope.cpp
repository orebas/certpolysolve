// The cluster start's scope in the certificate attempt (src/certify/attempts/newton_cover.cpp;
// heuristics::cluster_start_points, cluster_start_pool and cluster_start_in_scope;
// dev/docs/cluster-start-limit-20261117.md).
//
// Since 17 November 2026 the count of 16 is of the points in the start's groups, and the points the compensated tier
// leaves altogether are bounded by a count that follows the degree, max(16, n / 64). Before, the 16 counted every
// point left, so that the points in no group, whose number grows with the degree, switched the start off.
//
// - The rule itself (cluster_start_in_scope), over points left, points in groups and degrees: where at most 16 points
//   are left it is the rule before the change at every degree, and up to degree 1,087 it is that rule altogether.
// - The attempt, on a dense polynomial with real coefficients times m pairs of roots 2^-45 apart, the pairs at the
//   m-th roots of unity times 1/2 (the factor (x^m - a^m)(x^m - b^m), a = 1/2, b = 1/2 + 2^-45, is of modulus
//   about 1 on the unit circle, so the dense polynomial's roots keep their conditioning). Binary64 cannot separate
//   a pair: the tier passes its two points on, the start makes a group of them, and nothing else in the attempt
//   finishes them.
//   - degree 4,112, eight pairs: the tier leaves their 16 points and two others, 18 in all. The start takes its
//     eight groups and the attempt certifies every root. Before the change it refused: the start did not run.
//   - degree 4,114, nine pairs: 18 points in groups, over the count. The start does not run; the attempt refuses.
//   - degree 1,040, eight pairs: the tier leaves the 16 points alone, and the start takes them, as before.
//   Each fixture's premise (how many points the tier left) is checked, so that a fixture that stops exercising its
//   case fails and says so. The sources are declared without real structure: no conjugate pairs are taken, and
//   the points the tier leaves are the points it does not certify.
// Not exercised by a fixture: more than 16 points left at a degree where the bound is 16. No polynomial of such a
// degree was found on which the tier leaves a point outside every group beside 16 points in groups (24 seeds leave
// none at degree 1,040; a further pair of roots changes the groups). The rule's own test covers that case.
#include "certify/attempts/newton_cover.hpp"
#include "core/engine_policy.hpp"
#include "core/engine_statistics.hpp"
#include "core/heuristics.hpp"
#include "core/query.hpp" // cps_options
#include "core/raii.hpp"
#include "numeric/evaluate/structured_source.hpp"
#include <certpolysolve/certpolysolve.hpp>
#include <flint/acb_poly.h>
#include <flint/fmpz.h>
#include <cstdint>
#include <iostream>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace d = certpolysolve::detail;
namespace h = certpolysolve::detail::heuristics;

namespace {
void require(bool condition, const std::string &message) {
    if (!condition)
        throw std::runtime_error(message);
}

// The dense polynomials' seed. With it the tier leaves two points outside the pairs at degree 4,112 and none at 1,040
// (of the seeds 1 to 24, sixteen leave one to three at degree 4,112 and none leaves any at 1,040).
constexpr uint64_t SEED = 1;

// Monic of degree n, the other coefficients m 2^-40 with |m| <= 2^40, from the generator's own 64-bit outputs (the
// same polynomial under every standard library).
d::ComplexPolynomial dense(uint64_t seed, slong n) {
    std::mt19937_64 rng(seed);
    d::ComplexPolynomial p;
    d::ComplexBall c;
    for (slong k = 0; k < n; ++k) {
        const int64_t mantissa = static_cast<int64_t>(rng() >> 23) - (int64_t(1) << 40);
        acb_set_si(c.value(), mantissa);
        acb_mul_2exp_si(c.value(), c.value(), -40);
        acb_poly_set_coeff_acb(p.value(), k, c.value());
    }
    acb_poly_set_coeff_si(p.value(), n, 1);
    return p;
}

// dense(seed, n) times (x^m - a^m)(x^m - b^m), a = 1/2 and b = 1/2 + 2^-45: m pairs of roots 2^-45 apart. Exact.
d::ComplexPolynomial with_pairs(uint64_t seed, slong n, slong pairs) {
    d::ComplexPolynomial p = dense(seed, n), factor;
    d::ComplexBall c;
    fmpz_t power;
    fmpz_init(power);
    for (const bool wider : {false, true}) {
        // a^m = 2^-m; b^m = (2^44 + 1)^m 2^(-45 m).
        fmpz_one(power);
        if (wider) {
            fmpz_set_ui(power, 1);
            fmpz_mul_2exp(power, power, 44);
            fmpz_add_ui(power, power, 1);
            fmpz_pow_ui(power, power, static_cast<ulong>(pairs));
        }
        acb_set_fmpz(c.value(), power);
        acb_mul_2exp_si(c.value(), c.value(), wider ? -45 * pairs : -pairs);
        acb_neg(c.value(), c.value());
        acb_poly_zero(factor.value());
        acb_poly_set_coeff_si(factor.value(), pairs, 1);
        acb_poly_set_coeff_acb(factor.value(), 0, c.value());
        acb_poly_mul(p.value(), p.value(), factor.value(), 8192);
    }
    fmpz_clear(power);
    for (slong k = 0; k <= acb_poly_degree(p.value()); ++k)
        require(acb_is_exact(acb_poly_get_coeff_ptr(p.value(), k)), "a test polynomial's coefficient is not exact");
    return p;
}

struct Outcome {
    bool certified = false;
    std::vector<d::Root> cover;
    d::PartialNewtonCover partial;
    d::StudyStatistics counts;
    uint64_t left = 0; // the points the compensated tier did not certify
};

// The attempt as a route calls it before its own method, on a source without real structure.
Outcome attempt(const d::ComplexPolynomial &polynomial) {
    const slong degree = acb_poly_degree(polynomial.value());
    d::PolynomialSourceOptions source_options;
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
        degree, false, source_options);
    cps_options options;
    options.goal = CPS_ROOT_APPROXIMATION;
    options.bits = 53;
    d::EnginePolicy policy;
    policy.study_statistics = std::make_shared<d::StudyStatistics>();
    Outcome out;
    d::CandidateStatistics work{};
    out.certified = d::machine_newton_cover(out.cover, source, nullptr, options, 85, true, work, true, &policy, nullptr,
                                            nullptr, nullptr, true, &out.partial);
    out.counts = *policy.study_statistics;
    out.left = static_cast<uint64_t>(degree) - out.counts.machine_cover_compensated_certificates;
    return out;
}

std::string describe(const Outcome &o, uint64_t degree) {
    return " (degree " + std::to_string(degree) + ": the tier left " + std::to_string(o.left) + " points, " +
           std::to_string(o.counts.machine_cover_cluster_starts) + " groups started, " +
           (o.certified ? "certified" : "refused") + ")";
}

void rule_tests() {
    require(h::cluster_start_points == 16, "the start's count of points in groups is not 16");
    for (const uint64_t degree : {uint64_t{1}, uint64_t{16}, uint64_t{64}, uint64_t{1024}, uint64_t{1087}})
        require(h::cluster_start_pool(degree) == h::cluster_start_points, "the bound is not the count up to degree 1,087");
    require(h::cluster_start_pool(1088) == 17 && h::cluster_start_pool(4112) == 64 && h::cluster_start_pool(32768) == 512 &&
                h::cluster_start_pool(uint64_t{1} << 40) == (uint64_t{1} << 34),
            "the bound is not one point in 64 above degree 1,087");
    for (const uint64_t degree : {uint64_t{2}, uint64_t{8}, uint64_t{256}, uint64_t{1087}, uint64_t{1088}, uint64_t{4112},
                                  uint64_t{32770}, uint64_t{1} << 20})
        for (size_t left = 0; left <= 600; ++left)
            for (size_t members = 0; members <= left + 4; ++members) {
                // The rule before the change: groups, among at most 16 points left.
                const bool before = members >= 2 && left <= 16;
                const bool now = h::cluster_start_in_scope(left, members, degree);
                if (left <= 16 || degree <= 1087)
                    require(now == before, "the rule differs from the rule before the change where it must not");
                else
                    require(now == (members >= 2 && members <= 16 && left <= degree / 64),
                            "the rule above 16 points left is not: at most 16 in groups among at most n / 64");
            }
}

void attempt_tests() {
    // Sixteen points in groups among more than sixteen left, above the bound's degree.
    {
        const d::ComplexPolynomial p = with_pairs(SEED, 4096, 8);
        const uint64_t n = static_cast<uint64_t>(acb_poly_degree(p.value()));
        const Outcome o = attempt(p);
        require(o.left > h::cluster_start_points && o.left <= h::cluster_start_pool(n),
                "the fixture with eight pairs no longer leaves more than 16 points within the bound" + describe(o, n));
        require(o.certified && o.cover.size() == n, "eight pairs among more than 16 points left are not certified" + describe(o, n));
        require(o.counts.machine_cover_cluster_starts == 8, "the start did not take the eight pairs" + describe(o, n));
    }
    // Eighteen points in groups: over the count.
    {
        const d::ComplexPolynomial p = with_pairs(SEED, 4096, 9);
        const uint64_t n = static_cast<uint64_t>(acb_poly_degree(p.value()));
        const Outcome o = attempt(p);
        require(o.left > h::cluster_start_points && o.left <= h::cluster_start_pool(n),
                "the fixture with nine pairs no longer leaves more than 16 points within the bound" + describe(o, n));
        require(!o.certified && o.counts.machine_cover_cluster_starts == 0,
                "eighteen points in groups are within the start's scope" + describe(o, n));
    }
    // Sixteen points left, all in groups, where the bound is 16: the start takes them, as before the change.
    {
        const d::ComplexPolynomial p = with_pairs(SEED, 1024, 8);
        const uint64_t n = static_cast<uint64_t>(acb_poly_degree(p.value()));
        const Outcome o = attempt(p);
        require(h::cluster_start_pool(n) == h::cluster_start_points && o.left == h::cluster_start_points,
                "the fixture of degree 1,040 no longer leaves the 16 points of its pairs" + describe(o, n));
        require(o.certified && o.cover.size() == n && o.counts.machine_cover_cluster_starts == 8,
                "eight pairs at degree 1,040 are not certified by their starts" + describe(o, n));
    }
}
} // namespace

int main() {
    certpolysolve::ThreadScope lifetime;
    try {
        rule_tests();
        attempt_tests();
        std::cout << "cluster start scope tests passed\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "FAILED: " << error.what() << '\n';
        return 1;
    }
}
