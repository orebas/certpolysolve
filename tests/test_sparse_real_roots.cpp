// B6's isolation (sparse_real_roots.hpp) against the reference Sturm count (exact_reference_api.hpp's rational Sturm
// sequence, cheap up to degree about 40 with small coefficients; past that its rationals grow without bound):
// random fewnomials, trinomials and binomials with counts known by construction up to degree 150, zero's
// multiplicity and the negative axis, the planted common root of consecutive derivatives that must stop
// honestly, and cancellation.
#include "exact_reference_api.hpp"
#include "methods/real/sparse_real_roots.hpp"
#include <certpolysolve/certpolysolve.hpp>
#include <flint/arf.h>
#include <flint/fmpq.h>
#include <algorithm>
#include <chrono>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

using namespace certpolysolve::detail;
namespace {
void require(bool value, const std::string &message) {
    if (!value) throw std::runtime_error(message);
}
// An interval end is a dyadic point or a tiny ball with the polynomial's sign proved over all of it (so no root
// in it): the closed interval of the outer bounds holds the same single root.
reference::Rational bound_of(const arb_t x, bool upper) {
    arf_t a;
    arf_init(a);
    if (upper) arb_get_ubound_arf(a, x, 4096);
    else arb_get_lbound_arf(a, x, 4096);
    fmpq_t q;
    fmpq_init(q);
    arf_get_fmpq(q, a);
    char *text = fmpq_get_str(nullptr, 10, q);
    reference::Rational out(text);
    flint_free(text);
    fmpq_clear(q);
    arf_clear(a);
    return out;
}
// The reference polynomial (dense, ascending) of sparse terms with integer coefficients.
reference::Polynomial dense_of(const std::vector<std::pair<slong, long long>> &terms) {
    reference::Polynomial p(static_cast<size_t>(terms.back().first) + 1, reference::Rational(0));
    for (const auto &[e, c] : terms) p[static_cast<size_t>(e)] = reference::Rational(c);
    return p;
}
std::vector<SparseRealTerm> sparse_of(const std::vector<std::pair<slong, long long>> &terms) {
    std::vector<SparseRealTerm> out;
    for (const auto &[e, c] : terms) out.push_back({e, Number(std::to_string(c))});
    return out;
}
// Every interval holds exactly one root of p by the reference Sturm count (closed interval, exact ends), the
// intervals are disjoint and ascending, and together with zero they hold every real root.
void check_against_sturm(const reference::Polynomial &dense, const SparseRealIsolation &iso, const std::string &label) {
    require(iso.complete, label + ": not complete: " + iso.stop);
    // The isolation is of the polynomial with its zero root stripped: an interval may start at the breakpoint 0,
    // which a closed-interval count over the unstripped polynomial would charge with the root at zero.
    require(iso.zero_multiplicity < dense.size(), label + ": zero multiplicity");
    const reference::Polynomial p(dense.begin() + static_cast<long>(iso.zero_multiplicity), dense.end());
    int found = 0;
    const auto check_side = [&](const std::vector<SparseRootInterval> &roots, bool negative) {
        reference::Rational previous_high = 0;
        for (const auto &root : roots) {
            reference::Rational lo = bound_of(root.low.value(), false), hi = bound_of(root.high.value(), true);
            require(lo <= hi && lo >= previous_high, label + ": intervals not ascending and disjoint");
            previous_high = hi;
            if (negative) std::swap(lo, hi), lo = -lo, hi = -hi;
            require(reference::roots_closed(p, lo, hi) == 1, label + ": an interval without exactly one root");
            ++found;
        }
    };
    check_side(iso.positive, false);
    check_side(iso.negative, true);
    // All real roots of the stripped polynomial: Sturm on a wide closed interval.
    const reference::Rational wide = reference::power2(200);
    const int all = reference::roots_closed(p, -wide, wide);
    require(all == found, label + ": " + std::to_string(all) + " real roots by Sturm, " + std::to_string(found) +
                              " isolated");
}

void trinomials_and_binomials() {
    cps_options o;
    for (const slong n : {4, 10, 31, 64, 120}) {
        // x^n - 2: one positive root; one negative for even n (Sturm up to degree 31; the count by construction).
        const std::vector<std::pair<slong, long long>> terms{{0, -2}, {n, 1}};
        SparseRealIsolator iso(sparse_of(terms), {});
        const auto result = iso.isolate(o);
        if (n <= 31) check_against_sturm(dense_of(terms), result, "x^" + std::to_string(n) + " - 2");
        require(result.complete && result.positive.size() == 1 &&
                    result.negative.size() == static_cast<size_t>(n % 2 == 0),
                "x^n - 2: root counts at n = " + std::to_string(n) + " (" + result.stop + ")");
    }
    // x^n - 3 x^e + 1 (the corpus's trinomials, e odd): two positive roots by Descartes' rule (two sign variations,
    // p(0) > 0 > p(1) and p(+inf) > 0) and no negative one (none for p(-x)).
    for (const slong n : {12, 60, 150}) {
        const slong e = 2 * (n / 6) + 1;
        const std::vector<std::pair<slong, long long>> terms{{0, 1}, {e, -3}, {n, 1}};
        SparseRealIsolator iso(sparse_of(terms), {});
        const auto result = iso.isolate(o);
        if (n <= 12) check_against_sturm(dense_of(terms), result, "trinomial-" + std::to_string(n));
        require(result.complete && result.positive.size() == 2 && result.negative.empty(),
                "the trinomial's two positive roots at n = " + std::to_string(n) + " (" + result.stop + ")");
    }
}

void random_fewnomials() {
    std::mt19937_64 rng(20261016);
    cps_options o;
    int roots = 0, inputs = 0;
    for (int t = 0; t < 60; ++t) {
        const slong degree = 2 + static_cast<slong>(rng() % 39);
        // At most degree + 1 distinct exponents exist.
        const size_t count = std::min<size_t>(2 + static_cast<size_t>(rng() % 6), static_cast<size_t>(degree) + 1);
        std::vector<std::pair<slong, long long>> terms{{0, 0}, {degree, 0}};
        while (terms.size() < count) {
            const slong e = static_cast<slong>(rng() % degree);
            bool seen = false;
            for (const auto &[f, c] : terms) seen = seen || f == e;
            if (!seen) terms.push_back({e, 0});
        }
        std::sort(terms.begin(), terms.end());
        for (auto &[e, c] : terms) {
            c = static_cast<long long>(rng() % 101) - 50;
            if (c == 0) c = 7;
        }
        std::string label = "random fewnomial " + std::to_string(t) + " (";
        for (const auto &[e, c] : terms) label += std::to_string(c) + " x^" + std::to_string(e) + " ";
        label += ")";
        SparseRealIsolator iso(sparse_of(terms), {});
        const auto result = iso.isolate(o);
        check_against_sturm(dense_of(terms), result, label);
        roots += static_cast<int>(result.positive.size() + result.negative.size());
        ++inputs;
    }
    std::cout << "random fewnomials: " << inputs << " inputs, " << roots << " real roots, all by Sturm\n";
}

void zero_and_negative() {
    // x^2 (x^4 - 2): zero twice, one positive and one negative root.
    const std::vector<std::pair<slong, long long>> terms{{2, -2}, {6, 1}};
    cps_options o;
    SparseRealIsolator iso(sparse_of(terms), {});
    const auto result = iso.isolate(o);
    require(result.complete && result.zero_multiplicity == 2 && result.positive.size() == 1 &&
                result.negative.size() == 1,
            "x^2 (x^4 - 2): zero twice and one root each side");
    check_against_sturm(dense_of(terms), result, "x^2 (x^4 - 2)");
}

void honest_stop() {
    // 45 (x^9/9 - 6 x^5/5 + 9 x + 1) = 5 x^9 - 54 x^5 + 405 x + 45: p_1 = p' = 45 (x^4 - 3)^2 has the double root
    // 3^(1/4), which is the root of p_2 too, so the sign of p_1 there is exactly 0: undecidable by approximation.
    const std::vector<std::pair<slong, long long>> terms{{0, 45}, {1, 405}, {5, -54}, {9, 5}};
    cps_options o;
    o.max_work_bits = 4096;
    SparseRealIsolator::Settings settings;
    settings.cap_bits = 4096;
    SparseRealIsolator iso(sparse_of(terms), settings);
    const auto result = iso.isolate(o);
    require(!result.complete, "a common root of consecutive derivatives was not noticed");
    require(result.stop.find("common root") != std::string::npos, "the stop does not name the cause: " + result.stop);
    require(result.deepest_bits == 4096, "the stop came before the cap");
    std::cout << "honest stop: " << result.stop << " after " << result.evaluations << " evaluations\n";
}

void cancellation() {
    cps_options o;
    int calls = 0;
    o.cancel = [](void *data) -> int { return ++*static_cast<int *>(data) > 3; };
    o.cancel_data = &calls;
    // Undecidable at any precision (a double root), so the loop runs until cancelled.
    const std::vector<std::pair<slong, long long>> stuck{{0, 45}, {1, 405}, {5, -54}, {9, 5}};
    SparseRealIsolator iso(sparse_of(stuck), {});
    const auto result = iso.isolate(o);
    require(!result.complete && result.stop == "cancelled", "cancellation ignored: " + result.stop);
}

// x^n - 2^(n k) has the simple positive root 2^k, and a negative root exactly when n is even.
// The exponents themselves fit in a few limbs; expanding the coefficients or endpoints is unnecessary.
// This checks both the isolation bound and refinement across binades beyond a signed machine word.
void compact_exponents() {
    const reference::Integer beyond = (reference::Integer(1) << (sizeof(slong) * 8 - 1)) + 1;
    for (const int direction : {1, -1}) {
        const reference::Integer k = direction * beyond;
        const std::string k_text = k.str();
        RealBall known;
        Number("1@" + k_text).enclose(known.value(), 64);
        require(arb_is_exact(known.value()), "known compact root was rounded");
        for (const slong degree : {1, 2, 127, 128}) {
            const reference::Integer exponent = k * degree;
            const std::vector<SparseRealTerm> terms{{0, Number("-1@" + exponent.str())}, {degree, Number("1")}};
            cps_options options;
            options.max_work_bits = 256;
            SparseRealIsolator::Settings settings;
            settings.cap_bits = 256;
            auto isolated = SparseRealIsolator(terms, settings).isolate(options);
            const std::string label = "compact binomial n=" + std::to_string(degree) + ", k=" + k_text;
            require(isolated.complete && isolated.positive.size() == 1 &&
                        isolated.negative.size() == static_cast<size_t>(degree % 2 == 0),
                    label + ": incorrect complete real-root count");
            SparseRealEvaluator evaluator(terms);
            for (auto *side : {&isolated.positive, &isolated.negative}) {
                for (auto &root : *side) {
                    // For even n, p(-x) = p(x), so both stored positive half-axis intervals use this evaluator.
                    require(arb_le(root.low.value(), known.value()) && arb_le(known.value(), root.high.value()),
                            label + ": Cauchy bound lost a root");
                    uint64_t evaluations = 0;
                    RealSignBudget budget{64, 256, 64, evaluations, true};
                    require(real_refine(evaluator, root.low, root.high, 32, budget, options),
                            label + ": compact-exponent refinement failed");
                    require(arb_le(root.low.value(), known.value()) && arb_le(known.value(), root.high.value()),
                            label + ": refinement lost the exact known root");
                    RealBall width, allowed;
                    arb_sub(width.value(), root.high.value(), root.low.value(), 128);
                    arb_mul_2exp_si(allowed.value(), known.value(), -32);
                    require(arb_le(width.value(), allowed.value()), label + ": refinement missed its width");
                }
            }
        }
    }

    // Public Auto must never call an empty search complete. A resource stop may remain incomplete;
    // any complete answer must include exactly the binomial's known real roots.
    for (const uint64_t degree : {127, 128}) {
        const reference::Integer exponent = beyond * degree;
        certpolysolve::Polynomial p;
        p.set(0, "-1@" + exponent.str()).set(degree, 1);
        certpolysolve::Options options;
        options.bits = 53;
        options.max_bits = 256;
        options.work_budget = 8192;
        options.time_limit = std::chrono::seconds(5);
        const auto result = certpolysolve::real_roots(p, options);
        require(result.get() && result.stop() != certpolysolve::Stop::error, "compact public input was rejected");
        if (!result.complete()) {
            require(!result.goal_met(), "incomplete compact answer claims success");
            continue;
        }
        require(result.roots().size() == (degree % 2 ? 1u : 2u), "complete compact answer omitted a real root");
        for (const int sign : {-1, 1}) {
            if (sign < 0 && degree % 2) continue;
            const certpolysolve::Dyadic known(std::to_string(sign), beyond.str());
            size_t held = 0;
            for (const auto &root : result.roots()) {
                held += root.box.real.contains(known) && root.box.imag.contains(certpolysolve::Dyadic(0));
                require(root.count == 1 && root.multiplicity.value_or(0) == 1 &&
                            root.realness == certpolysolve::Realness::real,
                        "compact root's count, multiplicity or realness is wrong");
            }
            require(held == 1, "complete compact answer does not enclose its exact known root");
        }
    }
}
} // namespace

int main() {
    try {
        compact_exponents();
        trinomials_and_binomials();
        random_fewnomials();
        zero_and_negative();
        honest_stop();
        cancellation();
    } catch (const std::exception &error) {
        std::cerr << "FAIL " << error.what() << '\n';
        return 1;
    }
    std::cout << "sparse real roots tests passed\n";
    return 0;
}
