// The groups of the certificate attempt's first pass (src/certify/attempts/newton_cover.cpp).
//
// The first pass asks the compensated tier about every point in turn. Since 17 November 2026 it takes the points a
// group at a time (heuristics::compensated_ahead_points): what the tier computes first at a point is made for the
// group's points together, and so are the second tries of the points whose first try asks for one. A group of one
// point is the pass as it was before: each point's tries in turn. The kernels under those steps are compared with
// their references in tests/test_kernel_variants.cpp; this file checks the pass itself.
//
// machine_newton_cover runs with groups of 1, 2, 3, 5, 63 and 64 points under every kernel variant the machine has,
// and every run must give what the reference variant gives with one point a group: the same verdict, the same regions
// (ball, endpoints, count, certificate), the same partial cover where the attempt refuses, and the same counters.
// The polynomials:
// - x^n + (n x^(n-1) + ... + 2 x + 1) / (n + 1) for n = 511 and 1,023, the roots of MPSolve's "easy" polynomial,
//   where the binary64 iteration leaves many points a second try away from their certificate (63 and 231 of them;
//   required: the deferred tries must have been exercised), and for n = 127, where it leaves none;
// - random dense polynomials with real and with complex coefficients;
// - one with two roots 2^-45 apart, whose points the tier passes on to the tiers after it (required);
// - degrees 1, 2 and 65, around the group's ends.
// Checked by breaking the pass (17 November 2026): without the last group's second tries the test fails. Not
// exercised here: the order of the points passed on where one of them was passed on after a second try (the pass
// sorts them); no small polynomial was found that has one beside a point refused at once.
#include "certify/attempts/newton_cover.hpp"
#include "core/cover.hpp"
#include "core/engine_policy.hpp"
#include "core/engine_statistics.hpp"
#include "core/heuristics.hpp"
#include "core/query.hpp"
#include "core/raii.hpp"
#include "numeric/evaluate/structured_source.hpp"
#include "platform/kernel_variant.hpp"
#include <certpolysolve/certpolysolve.hpp>
#include <flint/acb_poly.h>
#include <cstdint>
#include <iostream>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace d = certpolysolve::detail;

namespace {
void require(bool condition, const std::string &message) {
    if (!condition)
        throw std::runtime_error(message);
}

// One run's outcome: everything the attempt returns or counts that is not a timer.
struct Outcome {
    bool certified = false;
    std::vector<d::Root> cover;
    d::PartialNewtonCover partial;
    slong certificate_bits = 0;
    bool machine_converged = false;
    d::CandidateStatistics work{};
    d::StudyStatistics counts;
};

bool same_root(const d::Root &a, const d::Root &b) {
    if (!acb_equal(a.box.value(), b.box.value()) || a.count != b.count || a.multiplicity != b.multiplicity ||
        a.realness != b.realness || a.certificate != b.certificate)
        return false;
    for (size_t k = 0; k < a.endpoints.size(); ++k)
        if (!arf_equal(a.endpoints[k].value(), b.endpoints[k].value()))
            return false;
    return true;
}

void require_same(const Outcome &a, const Outcome &b, const std::string &what) {
    require(a.certified == b.certified, what + ": another verdict");
    require(a.certificate_bits == b.certificate_bits && a.machine_converged == b.machine_converged,
            what + ": another precision or convergence report");
    require(a.cover.size() == b.cover.size(), what + ": another number of regions");
    for (size_t i = 0; i < a.cover.size(); ++i)
        require(same_root(a.cover[i], b.cover[i]), what + ": region " + std::to_string(i) + " differs");
    require(a.partial.certified == b.partial.certified && a.partial.pass_complete == b.partial.pass_complete &&
                a.partial.certificate_bits == b.partial.certificate_bits &&
                a.partial.roots.size() == b.partial.roots.size() && bool(a.partial.points) == bool(b.partial.points),
            what + ": another partial cover");
    for (size_t i = 0; i < a.partial.roots.size(); ++i) {
        if (a.partial.certified[i])
            require(same_root(a.partial.roots[i], b.partial.roots[i]),
                    what + ": the partial cover's region " + std::to_string(i) + " differs");
        if (a.partial.points)
            require(acb_equal(a.partial.points->data() + i, b.partial.points->data() + i),
                    what + ": the partial cover's point " + std::to_string(i) + " differs");
    }
    // The counters that are not timers.
    for (const auto field : {&d::CandidateStatistics::candidate_sweeps, &d::CandidateStatistics::evaluations,
                             &d::CandidateStatistics::evaluation_terms, &d::CandidateStatistics::root_interactions,
                             &d::CandidateStatistics::validation_calls, &d::CandidateStatistics::secular_regenerations,
                             &d::CandidateStatistics::double_seed_calls, &d::CandidateStatistics::candidate_stagnations,
                             &d::CandidateStatistics::real_refinement_steps, &d::CandidateStatistics::real_bisections})
        require(a.work.*field == b.work.*field, what + ": another work counter");
    for (const auto field :
         {&d::StudyStatistics::machine_cover_precision_retries, &d::StudyStatistics::machine_cover_peak_bits,
          &d::StudyStatistics::machine_cover_newton_steps, &d::StudyStatistics::machine_cover_block_evaluations,
          &d::StudyStatistics::machine_cover_stepped_certificates, &d::StudyStatistics::machine_cover_rung_evaluations,
          &d::StudyStatistics::machine_cover_rounds, &d::StudyStatistics::machine_cover_disks,
          &d::StudyStatistics::machine_cover_compensated_certificates,
          &d::StudyStatistics::machine_cover_compensated_later_tries,
          &d::StudyStatistics::machine_cover_double_double_repairs, &d::StudyStatistics::machine_cover_cluster_starts})
        require(a.counts.*field == b.counts.*field, what + ": another study counter");
}

// The attempt on an exact polynomial (the midpoints of `polynomial`'s coefficients) with `group` points a group.
Outcome attempt(const d::ComplexPolynomial &polynomial, size_t group, uint64_t bits) {
    const slong degree = acb_poly_degree(polynomial.value());
    d::PolynomialSource source(
        [&](acb_poly_t out, slong precision) {
            d::ComplexBall a;
            acb_poly_zero(out);
            for (slong k = 0; k <= degree; ++k) {
                acb_get_mid(a.value(), acb_poly_get_coeff_ptr(polynomial.value(), k));
                acb_set_round(a.value(), a.value(), precision);
                acb_poly_set_coeff_acb(out, k, a.value());
            }
        },
        degree, false);
    cps_options options;
    options.goal = CPS_ROOT_APPROXIMATION;
    options.bits = bits;
    d::EnginePolicy policy;
    policy.study_statistics = std::make_shared<d::StudyStatistics>();
    policy.compensated_ahead_points = group;
    Outcome out;
    out.certified = d::machine_newton_cover(out.cover, source, nullptr, options, 85, true, out.work, true, &policy,
                                            nullptr, &out.machine_converged, &out.certificate_bits, true, &out.partial);
    out.counts = *policy.study_statistics;
    return out;
}

struct Totals {
    uint64_t runs = 0, later_tries = 0, passed_on = 0;
};
// What a polynomial is in the test for: the run must have exercised it.
enum class Expect { Nothing, LaterTries, PassedOn };

// Every group size under every variant against one point a group under the reference.
void compare(const d::ComplexPolynomial &polynomial, uint64_t bits, const std::string &name, Totals &totals,
             Expect expect) {
    const d::KernelVariant before = d::kernel_variant();
    d::set_kernel_variant(d::KernelVariant::Reference);
    const Outcome reference = attempt(polynomial, 1, bits);
    const uint64_t degree = static_cast<uint64_t>(acb_poly_degree(polynomial.value()));
    const uint64_t passed_on = degree - reference.counts.machine_cover_compensated_certificates;
    if (expect == Expect::LaterTries)
        require(reference.counts.machine_cover_compensated_later_tries > 0,
                name + ": no point was certified on a second try, so the deferred tries were not exercised");
    if (expect == Expect::PassedOn)
        require(passed_on >= 2, name + ": the compensated tier passed on fewer than two points, so their order "
                                       "was not exercised");
    totals.later_tries += reference.counts.machine_cover_compensated_later_tries;
    totals.passed_on += passed_on;
    for (const auto variant : {d::KernelVariant::Reference, d::KernelVariant::Avx2Fma}) {
        if (!d::kernel_variant_available(variant))
            continue;
        d::set_kernel_variant(variant);
        for (const size_t group : {size_t(1), size_t(2), size_t(3), size_t(5), size_t(63), size_t(64), size_t(0)}) {
            const Outcome outcome = attempt(polynomial, group, bits);
            require_same(reference, outcome,
                         name + ", groups of " + std::to_string(group) + " points, variant " +
                             std::to_string(static_cast<int>(variant)));
            ++totals.runs;
        }
    }
    d::set_kernel_variant(before);
}

// x^n + (n x^(n-1) + ... + 2 x + 1) / (n + 1), n + 1 a power of two: the coefficients are exact.
d::ComplexPolynomial easy(slong n) {
    d::ComplexPolynomial p;
    d::ComplexBall c;
    for (slong k = 0; k <= n; ++k) {
        acb_set_si(c.value(), k + 1);
        acb_div_si(c.value(), c.value(), n + 1, 64);
        acb_poly_set_coeff_acb(p.value(), k, c.value());
    }
    return p;
}

// Monic, the other coefficients m 2^-40 with |m| <= 2^40 (and an imaginary part of that kind where asked).
d::ComplexPolynomial dense(std::mt19937_64 &rng, slong n, bool complex) {
    d::ComplexPolynomial p;
    d::ComplexBall c;
    std::uniform_int_distribution<int64_t> mantissa(-(int64_t(1) << 40), int64_t(1) << 40);
    for (slong k = 0; k < n; ++k) {
        acb_set_si_si(c.value(), mantissa(rng), complex ? mantissa(rng) : 0);
        acb_mul_2exp_si(c.value(), c.value(), -40);
        acb_poly_set_coeff_acb(p.value(), k, c.value());
    }
    acb_poly_set_coeff_si(p.value(), n, 1);
    return p;
}

// (x - 1/2)(x - 1/2 - 2^-45) times a dense polynomial: two roots that binary64 leaves about 2^-26 from their points,
// which the compensated tier passes on.
d::ComplexPolynomial close_pair(std::mt19937_64 &rng, slong n) {
    d::ComplexPolynomial p = dense(rng, n - 2, false), factor;
    d::ComplexBall root;
    for (const slong extra : {slong(0), slong(1)}) {
        acb_set_si(root.value(), -(slong(1) << 44) - extra);
        acb_mul_2exp_si(root.value(), root.value(), -45);
        acb_poly_zero(factor.value());
        acb_poly_set_coeff_si(factor.value(), 1, 1);
        acb_poly_set_coeff_acb(factor.value(), 0, root.value());
        acb_poly_mul(p.value(), p.value(), factor.value(), 4096);
    }
    return p;
}
} // namespace

int main() {
    certpolysolve::ThreadScope lifetime;
    try {
        std::mt19937_64 rng(20261117);
        Totals totals;
        compare(easy(511), 53, "the easy polynomial of degree 511", totals, Expect::LaterTries);
        compare(easy(1023), 53, "the easy polynomial of degree 1,023", totals, Expect::LaterTries);
        compare(easy(127), 53, "the easy polynomial of degree 127", totals, Expect::Nothing);
        compare(dense(rng, 130, false), 53, "a dense real polynomial of degree 130", totals, Expect::Nothing);
        compare(dense(rng, 97, true), 53, "a dense complex polynomial of degree 97", totals, Expect::Nothing);
        compare(close_pair(rng, 66), 53, "two roots 2^-45 apart in degree 66", totals, Expect::PassedOn);
        for (const slong n : {slong(1), slong(2), slong(65)})
            compare(dense(rng, n, n == 2), 53, "a dense polynomial of degree " + std::to_string(n), totals,
                    Expect::Nothing);
        std::cout << "attempt groups: " << totals.runs << " runs equal to one point a group under the reference ("
                  << totals.later_tries << " points certified on a later try, " << totals.passed_on
                  << " passed on by the compensated tier)\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << "test_attempt_groups: " << e.what() << "\n";
        return 1;
    }
}
