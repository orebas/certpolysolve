#include "exact_reference_api.hpp"
#include "certify/counting.hpp"
#include "core/service_meter.hpp"
#include "internal.hpp"
#include <iostream>
#include <random>

namespace detail = certpolysolve::detail;
using reference::require;
using reference::Rational;
using reference::Complex;

namespace {
void set(arb_t out, const Rational& value) { detail::Number(value.str()).enclose(out, 512); }
void set(acb_t out, const Complex& value) { set(acb_realref(out), value.real); set(acb_imagref(out), value.imag); }
void set(acb_poly_t out, const reference::Polynomial& polynomial) {
    acb_poly_zero(out);
    detail::ComplexBall coefficient;
    for (size_t i = 0; i < polynomial.size(); ++i) {
        set(coefficient.value(), {polynomial[i]});
        acb_poly_set_coeff_acb(out, static_cast<slong>(i), coefficient.value());
    }
}
void set_from_roots(acb_poly_t out, const std::vector<Complex>& roots) {
    acb_poly_zero(out);
    std::vector<Complex> polynomial(1, {1});
    for (const auto& root : roots) {
        std::vector<Complex> next(polynomial.size() + 1);
        for (size_t j = 0; j < polynomial.size(); ++j) {
            next[j] = next[j] + polynomial[j] * Complex(-root.real, -root.imag);
            next[j + 1] = next[j + 1] + polynomial[j];
        }
        polynomial = std::move(next);
    }
    detail::ComplexBall coefficient;
    for (size_t j = 0; j < polynomial.size(); ++j) {
        set(coefficient.value(), polynomial[j]);
        acb_poly_set_coeff_acb(out, static_cast<slong>(j), coefficient.value());
    }
}
void predicate_examples() {
    // N = ceil(log2(1 + log2 n)) + 5, retaining the historical floor of 9.
    static_assert(detail::subdivision_graeffe_steps(0) == 9);
    static_assert(detail::subdivision_graeffe_steps(1) == 9);
    static_assert(detail::subdivision_graeffe_steps(32768) == 9);
    static_assert(detail::subdivision_graeffe_steps(32769) == 10);
    static_assert(std::numeric_limits<slong>::digits <= 31 ||
                  detail::subdivision_graeffe_steps(static_cast<slong>(2147483648LL)) == 10);
    static_assert(std::numeric_limits<slong>::digits <= 31 ||
                  detail::subdivision_graeffe_steps(static_cast<slong>(2147483649LL)) == 11);
    detail::ComplexPolynomial p;
    detail::ComplexBall center;
    detail::Float radius;
    set(p.value(), {-1, 0, 1});
    arf_one(radius.value());
    require(detail::pellet_count(p.value(), center.value(), radius.value(), 128, 8).count == -1,
            "roots on the unit circle must not produce a strict disk certificate");
    arf_set_ui(radius.value(), 2);
    require(detail::pellet_count(p.value(), center.value(), radius.value(), 128).count == 2, "count two failed");
    arf_set_ui_2exp_si(radius.value(), 1, -1);
    require(detail::pellet_count(p.value(), center.value(), radius.value(), 128).count == 0, "empty disk failed");
    require(detail::pellet_count(p.value(), center.value(), radius.value(), 128, 4, 1).count == -1, "wrong expected count accepted");
    // A positive child-disk count only ends an exclusion attempt. The disk
    // of radius 3/2 contains the root 5/4, but the square of half-width one
    // does not. An independent square certificate must still refuse it.
    set(p.value(), {-5, 4});
    arf_set_ui_2exp_si(radius.value(), 3, -1);
    require(detail::pellet_count(p.value(), center.value(), radius.value(), 128).count == 1,
            "positive exclusion-disk count was not recovered");
    arf_one(radius.value());
    detail::Root outside_square;
    require(!detail::certify_cluster_box(outside_square, p.value(), center.value(), radius.value(), 1, 1, 128),
            "positive containing-disk count was mistaken for square occupancy");
    set(p.value(), {-1, 0, 1});
    arf_zero(radius.value());
    require(detail::pellet_count(p.value(), center.value(), radius.value(), 128).count == -1, "zero radius accepted");
    arf_set_si(radius.value(), -1);
    require(detail::pellet_count(p.value(), center.value(), radius.value(), 128).count == -1, "negative radius accepted");
    arf_nan(radius.value());
    require(detail::pellet_count(p.value(), center.value(), radius.value(), 128).count == -1, "NaN radius accepted");
    arf_pos_inf(radius.value());
    require(detail::pellet_count(p.value(), center.value(), radius.value(), 128).count == -1, "infinite radius accepted");
    arf_one(radius.value());
    arb_add_error_2exp_si(acb_realref(center.value()), -10);
    require(detail::pellet_count(p.value(), center.value(), radius.value(), 128).count == -1, "uncertain disk center accepted");
    acb_zero(center.value());

    // A root on the boundary stays inconclusive. Honor the caller's full
    // allowance rather than silently clipping it at the former limit of 16.
    set(p.value(), {-1, 1});
    const auto boundary = detail::pellet_count(p.value(), center.value(), radius.value(), 128, 17);
    require(boundary.count == -1 && boundary.graeffe_steps == 17,
            "the caller's Graeffe allowance was silently shortened");

    // In z + ... + z^16 + [0, 64], every nonconstant modulus is 1.
    // The sum of the other exact terms decisively rejects the chosen 1,
    // while zero exclusion still needs a tighter constant enclosure.
    // An any-count attempt must retain that precision diagnosis.
    acb_poly_zero(p.value());
    for (slong j = 1; j <= 16; ++j)
        acb_poly_set_coeff_si(p.value(), j, 1);
    detail::ComplexBall wide_constant;
    acb_set_ui(wide_constant.value(), 32);
    arb_add_error_2exp_si(acb_realref(wide_constant.value()), 5);
    acb_poly_set_coeff_acb(p.value(), 0, wide_constant.value());
    for (unsigned steps : {0u, 9u}) {
        const auto undecided = detail::pellet_count(p.value(), center.value(), radius.value(), 85, steps);
        require(undecided.count < 0 && undecided.width_limited && undecided.exclusion_limited,
                "any-count comparison lost the direct exclusion precision diagnosis");
    }

    // The same ambiguity can hide a nonconstant dominant coefficient. The
    // exact constant is zero, so a separate zero-count diagnosis cannot see
    // it. One member has coefficient 64 at degree eight, larger than the
    // other fifteen moduli together: a soft any-count refusal must not
    // rule out that strong Pellet margin. No uniform count is asserted.
    // Unlike the wide constant above, this does not obstruct direct zero
    // exclusion and must not trigger its finite-cap split-abandonment rule.
    acb_poly_set_coeff_si(p.value(), 0, 0);
    acb_poly_set_coeff_acb(p.value(), 8, wide_constant.value());
    const auto hidden_count = detail::pellet_count(p.value(), center.value(), radius.value(), 85, 0);
    require(hidden_count.count < 0 && hidden_count.needs_precision && !hidden_count.width_limited &&
                !hidden_count.exclusion_limited,
            "all-index precision request was confused with the direct exclusion stop");
    acb_poly_set_coeff_si(p.value(), 8, 64);
    require(detail::pellet_count(p.value(), center.value(), radius.value(), 85, 0).count == 8,
            "the nonconstant uncertainty fixture has no strong counted member");
    // In z^2 + [0,2]z, uncertainty obstructs the selected positive count,
    // but the exact zero constant decisively rejects zero dominance. Keep
    // the rotated-shift opportunity without abandoning a split at the cap.
    set(p.value(), {0, 1, 1});
    arb_add_error_2exp_si(acb_realref(p.value()->coeffs + 1), 0);
    const auto selected = detail::pellet_count(p.value(), center.value(), radius.value(), 85, 0);
    require(selected.count < 0 && selected.needs_precision && selected.width_limited && !selected.exclusion_limited,
            "positive-count uncertainty was confused with direct zero-exclusion uncertainty");
    // Exact equality must be a finite soft refusal, rather than an
    // instruction to keep increasing precision at a root on the circle.
    set(p.value(), {-1, 1});
    for (slong precision : {64, 128, 256}) {
        const auto equality = detail::pellet_count(p.value(), center.value(), radius.value(), precision, 0);
        require(equality.count < 0 && !equality.needs_precision && !equality.width_limited,
                "exact Pellet equality did not have a finite soft refusal");
    }

    // These three positive roots give no dominant untransformed coefficient,
    // but root squaring makes the unit disk's one-root count decidable.
    set_from_roots(p.value(), {{Rational(3) / 4}, {Rational(3) / 2}, {Rational(7) / 4}});
    arf_one(radius.value());
    require(detail::pellet_count(p.value(), center.value(), radius.value(), 128, 0).count == -1, "Graeffe test does not exercise Graeffe");
    auto count = detail::pellet_count(p.value(), center.value(), radius.value(), 128, 8);
    require(count.count == 1 && count.graeffe_steps > 0, "Graeffe did not certify known count");

    // Uniformity over uncertainty, including changes in real/nonreal root type.
    acb_poly_zero(p.value());
    detail::ComplexBall uncertain;
    arb_add_error_2exp_si(acb_realref(uncertain.value()), -10);
    arb_add_error_2exp_si(acb_imagref(uncertain.value()), -10);
    acb_poly_set_coeff_acb(p.value(), 0, uncertain.value());
    acb_poly_set_coeff_si(p.value(), 2, 1);
    arf_set_ui_2exp_si(radius.value(), 1, -3);
    require(detail::pellet_count(p.value(), center.value(), radius.value(), 128).count == 2, "uncertain cluster not certified");
    acb_poly_set_coeff_acb(p.value(), 2, uncertain.value());
    require(detail::pellet_count(p.value(), center.value(), radius.value(), 128).count == -1, "variable degree accepted");
}

// Dovetailing bounds a disk test by a work allowance (CountContext::allowance):
// it stops before any Taylor shift or Graeffe step the allowance does not
// cover, keeps what it computed, and resumed with the same arguments gives the
// unbounded test's result with the same work in total (no operation repeated).
void bounded_tests_resume() {
    detail::ComplexPolynomial p;
    set_from_roots(p.value(), {{Rational(7) / 8}, {Rational(9) / 8}, {Rational(-15) / 16}, {Rational(0), Rational(17) / 16}});
    detail::ComplexBall center;
    arb_set_si(acb_realref(center.value()), 1);
    arb_mul_2exp_si(acb_realref(center.value()), acb_realref(center.value()), -10);
    detail::Float radius;
    arf_one(radius.value());
    certpolysolve::detail::SearchStatistics statistics{};
    detail::CountContext free(statistics);
    const auto full = detail::pellet_count(p.value(), center.value(), radius.value(), 128, 9, -1, &free);
    require(full.count == 2 && full.graeffe_steps >= 2, "the bounded-test fixture does not need two Graeffe steps");
    require(!free.suspension && free.completed.empty(), "an unbounded context kept test state");
    certpolysolve::detail::SearchStatistics bounded_statistics{};
    detail::CountContext bounded(bounded_statistics);
    // An allowance below the Taylor shift: the test stops at once, having done nothing.
    bounded.allowance = 1;
    auto result = detail::pellet_count(p.value(), center.value(), radius.value(), 128, 9, -1, &bounded);
    require(result.suspended && bounded.suspended && bounded.metered == 0 && bounded.pending > 1,
            "a test started a Taylor shift its allowance did not cover");
    // Each operation then allowed as it comes: one resumption per operation.
    unsigned rounds = 0;
    while (result.suspended && rounds <= full.graeffe_steps + 1) {
        bounded.allowance = bounded.metered + bounded.pending;
        result = detail::pellet_count(p.value(), center.value(), radius.value(), 128, 9, -1, &bounded);
        ++rounds;
    }
    require(!result.suspended && result.count == full.count && result.graeffe_steps == full.graeffe_steps &&
                result.width_limited == full.width_limited && result.needs_precision == full.needs_precision &&
                result.exclusion_limited == full.exclusion_limited,
            "a resumed test did not give the unbounded test's result");
    require(rounds == full.graeffe_steps + 1, "a resumed test repeated or skipped an operation");
    require(bounded.metered == free.metered && bounded_statistics.graeffe_steps == statistics.graeffe_steps &&
                bounded_statistics.disk_tests == 1 && bounded_statistics.taylor_shifts == 1,
            "a resumed test metered or counted other work than the unbounded one");
    // A test completed under an allowance is not computed again.
    result = detail::pellet_count(p.value(), center.value(), radius.value(), 128, 9, -1, &bounded);
    require(result.count == full.count && bounded.metered == free.metered && bounded_statistics.disk_tests == 1,
            "a completed bounded test was computed again");

    // A direct zero-exclusion obstruction must survive each suspended
    // Graeffe step, even though subsequent comparisons use different balls.
    set(p.value(), {1, 1});
    arb_add_error_2exp_si(acb_realref(p.value()->coeffs), 1);
    acb_zero(center.value());
    certpolysolve::detail::SearchStatistics uncertain_statistics{};
    detail::CountContext uncertain(uncertain_statistics);
    uncertain.allowance = 1;
    result = detail::pellet_count(p.value(), center.value(), radius.value(), 128, 9, -1, &uncertain);
    for (unsigned attempt = 0; result.suspended && attempt < 11; ++attempt) {
        uncertain.allowance = uncertain.metered + uncertain.pending;
        result = detail::pellet_count(p.value(), center.value(), radius.value(), 128, 9, -1, &uncertain);
    }
    require(!result.suspended && result.count < 0 && result.graeffe_steps == 9 && result.exclusion_limited &&
                uncertain_statistics.disk_tests == 1 && uncertain_statistics.taylor_shifts == 1,
            "suspension lost the direct zero-exclusion precision diagnosis");
}

// Admission and debit must agree after changing the Taylor price. All roots
// of x^128-1 lie in the unit disk, disjoint from this small off-axis disk.
void priced_taylor_resumes() {
    detail::ComplexPolynomial polynomial;
    acb_poly_set_coeff_si(polynomial.value(), 0, -1);
    acb_poly_set_coeff_si(polynomial.value(), 128, 1);
    detail::ComplexBall center;
    acb_set_d_d(center.value(), 2, 2);
    detail::Float radius;
    arf_set_ui_2exp_si(radius.value(), 1, -4);
    detail::EnginePolicy policy;
    policy.research_service_meter = true;
    const auto quote = detail::service_taylor_quote(polynomial.value(), center.value(), 128);
    require(quote.operation == detail::ServiceOperation::TaylorDivideConquer,
            "the admission fixture does not exercise the D&C price");
    detail::SearchStatistics free_statistics{}, bounded_statistics{};
    detail::CountContext free(free_statistics), bounded(bounded_statistics);
    free.research_service_policy = bounded.research_service_policy = &policy;
    const auto full = detail::pellet_count(polynomial.value(), center.value(), radius.value(), 128, 9, -1, &free);
    require(full.count == 0, "the independent unit-circle root bound disagrees with counting");
    bounded.allowance = std::nextafter(quote.price() / detail::operation_weight(128), 0.0);
    auto result = detail::pellet_count(polynomial.value(), center.value(), radius.value(), 128, 9, -1, &bounded);
    require(result.suspended && bounded.metered == 0 && bounded_statistics.taylor_shifts == 0,
            "a D&C shift ran before its quote was covered");
    for (unsigned round = 0; result.suspended && round < 12; ++round) {
        bounded.allowance = bounded.metered + bounded.pending;
        result = detail::pellet_count(polynomial.value(), center.value(), radius.value(), 128, 9, -1, &bounded);
    }
    require(!result.suspended && result.count == full.count && bounded.metered == free.metered &&
                bounded_statistics.taylor_shifts == 1 && free_statistics.taylor_shifts == 1,
            "resumption repeated or mischarged a D&C shift");
}

// The caller's cancellation is polled where the allowance is checked
// (CountContext::cancellation): a cancelled test suspends before its next
// Taylor shift or Graeffe step, and an uncancelled one is unchanged.
void cancelled_tests_stop() {
    detail::ComplexPolynomial p;
    set_from_roots(p.value(), {{Rational(7) / 8}, {Rational(9) / 8}, {Rational(-15) / 16}, {Rational(0), Rational(17) / 16}});
    detail::ComplexBall center;
    arb_set_si(acb_realref(center.value()), 1);
    arb_mul_2exp_si(acb_realref(center.value()), acb_realref(center.value()), -10);
    detail::Float radius;
    arf_one(radius.value());
    certpolysolve::detail::SearchStatistics statistics{};
    detail::CountContext free(statistics);
    const auto full = detail::pellet_count(p.value(), center.value(), radius.value(), 128, 9, -1, &free);
    require(full.count == 2 && full.graeffe_steps >= 2, "the cancellation fixture does not need two Graeffe steps");
    // The callback answers "cancelled" from its `after`-th poll on.
    struct Polls {
        unsigned seen = 0, after = 0;
        static int poll(void *data) {
            auto &self = *static_cast<Polls *>(data);
            return ++self.seen > self.after;
        }
    };
    for (unsigned after = 0; after <= 2; ++after) {
        Polls polls{0, after};
        cps_options options;
        options.cancel = Polls::poll;
        options.cancel_data = &polls;
        certpolysolve::detail::SearchStatistics cancelled_statistics{};
        detail::CountContext cancelled(cancelled_statistics);
        cancelled.cancellation = &options;
        const auto result = detail::pellet_count(p.value(), center.value(), radius.value(), 128, 9, -1, &cancelled);
        // Polls before the Taylor shift, then before each Graeffe step.
        require(result.suspended && cancelled.suspended && result.count == -1,
                "a cancelled test did not stop");
        require(cancelled_statistics.taylor_shifts == (after >= 1 ? 1u : 0u) &&
                    cancelled_statistics.graeffe_steps == (after >= 1 ? after - 1 : 0u),
                "a cancelled test ran an operation after the poll that saw the cancellation");
    }
    Polls never{0, std::numeric_limits<unsigned>::max()};
    cps_options options;
    options.cancel = Polls::poll;
    options.cancel_data = &never;
    certpolysolve::detail::SearchStatistics polled_statistics{};
    detail::CountContext polled(polled_statistics);
    polled.cancellation = &options;
    const auto result = detail::pellet_count(p.value(), center.value(), radius.value(), 128, 9, -1, &polled);
    require(!result.suspended && result.count == full.count && result.graeffe_steps == full.graeffe_steps &&
                polled.metered == free.metered && !polled.suspension && polled.completed.empty(),
            "an uncancelled poll changed the test");
}

// The stop cover (complete_cover_at_stop) stops its ladder of half-sides for the caller's cancellation,
// before a half-side and, inside its disk tests, before a Taylor shift or Graeffe step, and still returns
// a complete cover; a poll that never cancels changes nothing (review of 14 October 2026, S03). Sixteen
// roots near 1 whose coefficients reach 2^13 put Cauchy's bound about 2^14: the ladder has many rungs.
void stop_covers_stop_for_cancellation() {
    std::vector<Complex> roots;
    for (int j = 0; j < 16; ++j) roots.push_back({Rational(1) + Rational(j) / 64, Rational(j % 3) / 32});
    detail::ComplexPolynomial p;
    set_from_roots(p.value(), roots);
    const slong precision = 128;
    const auto counted = [](const std::vector<detail::Root> &cover) {
        uint64_t total = 0;
        for (const auto &root : cover) total += root.count;
        return total;
    };
    std::vector<detail::Root> reference;
    require(detail::complete_cover_at_stop(reference, p.value(), precision) && counted(reference) == 16,
            "the stop cover fixture has no complete cover");
    struct Polls {
        unsigned seen = 0, after = 0;
        static int poll(void *data) {
            auto &self = *static_cast<Polls *>(data);
            return ++self.seen > self.after;
        }
    };
    Polls never{0, std::numeric_limits<unsigned>::max()};
    cps_options options;
    options.cancel = Polls::poll;
    options.cancel_data = &never;
    std::vector<detail::Root> polled;
    require(detail::complete_cover_at_stop(polled, p.value(), precision, &options), "the polled cover failed");
    require(polled.size() == reference.size(), "an uncancelled poll changed the stop cover");
    for (size_t i = 0; i < polled.size(); ++i)
        require(acb_equal(polled[i].box.value(), reference[i].box.value()) && polled[i].count == reference[i].count,
                "an uncancelled poll changed the stop cover");
    const unsigned polls = never.seen;
    require(polls >= 8, "the stop cover fixture polls too rarely to show a stop");
    for (unsigned after : {0u, 1u, polls / 2, polls - 1}) {
        Polls cancel{0, after};
        options.cancel_data = &cancel;
        std::vector<detail::Root> cover;
        require(detail::complete_cover_at_stop(cover, p.value(), precision, &options) && counted(cover) == 16,
                "a cancelled stop cover is not complete");
        // The poll that sees the cancellation suspends the test in progress (or ends the ladder); at most the
        // outer test of the same rung polls once more before the ladder ends.
        require(cancel.seen <= after + 2, "a cancelled stop cover went on testing");
    }
}

void rectangle_sandwich() {
    detail::ComplexPolynomial p;
    detail::ComplexBall center;
    detail::Float radius;
    set_from_roots(p.value(), {{0}, {Rational(7) / 8, Rational(7) / 8}});
    arf_one(radius.value());
    require(detail::pellet_count(p.value(), center.value(), radius.value(), 256, 8).count == 1, "inner disk setup failed");
    detail::Root out;
    out.certificate = "unchanged";
    require(!detail::certify_cluster_box(out, p.value(), center.value(), radius.value(), 1, 1, 256, 8),
            "disk count was incorrectly claimed for a square containing an extra corner root");
    require(out.certificate == "unchanged", "failed certificate changed output");

    set_from_roots(p.value(), {{Rational(1) / 16}, {Rational(-1) / 16}, {8}});
    arf_set_ui_2exp_si(radius.value(), 1, -2);
    require(detail::certify_cluster_box(out, p.value(), center.value(), radius.value(), 2, 3, 128), "cluster box failed");
    require(out.count == 6 && out.multiplicity == 0 && out.realness == CPS_REAL_UNKNOWN,
            "cluster count invented multiplicity or realness");
    set_from_roots(p.value(), {{0}, {8}});
    require(detail::certify_cluster_box(out, p.value(), center.value(), radius.value(), 1, 3, 128), "single factor root failed");
    require(out.count == 3 && out.multiplicity == 3 && out.realness == CPS_REAL, "factor multiplicity lost");

    set_from_roots(p.value(), {{0}, {0}, {0}, {8}});
    require(detail::certify_cluster_box(out, p.value(), center.value(), radius.value(), 3, 2, 128),
            "exact multiple zero certificate failed");
    require(out.count == 6 && out.multiplicity == 6 && out.realness == CPS_REAL && acb_is_zero(out.box.value()),
            "exact zero valuation did not establish the point multiplicity");
    detail::ComplexBall uncertain;
    arb_add_error_2exp_si(acb_realref(uncertain.value()), -30);
    acb_poly_set_coeff_acb(p.value(), 0, uncertain.value());
    require(detail::certify_cluster_box(out, p.value(), center.value(), radius.value(), 3, 2, 128),
            "uncertain multiple-zero neighborhood failed");
    require(out.multiplicity == 0 && !acb_is_zero(out.box.value()),
            "uncertainty was incorrectly turned into an exact zero factor");

    // A compact center far from zero and an extremely small radius: export at
    // insufficient precision enlarges the rectangle outside the outer disk.
    const Rational huge = reference::power2(1000);
    set_from_roots(p.value(), {{huge}});
    set(center.value(), {huge});
    arf_set_ui_2exp_si(radius.value(), 1, -1000);
    bool needs_precision = false;
    require(!detail::certify_cluster_box(out, p.value(), center.value(), radius.value(), 1, 1, 128,
                                        detail::default_graeffe_steps, nullptr, &needs_precision) && needs_precision,
            "outward serialization lost its precision request");
    require(detail::certify_cluster_box(out, p.value(), center.value(), radius.value(), 1, 1, 4096,
                                       detail::default_graeffe_steps, nullptr, &needs_precision) && !needs_precision,
            "fresh precision did not resolve the outward serialization obstruction");
    require(acb_equal(out.box.value(), center.value()) && out.count == 1 && out.multiplicity == 1,
            "the exact linear root was not preserved after outward serialization");
}

void randomized_known_counts() {
    std::mt19937 generator(20260906);
    int certificates = 0, graeffe_certificates = 0;
    for (int sample = 0; sample < 240; ++sample) {
        std::vector<Complex> roots;
        const int degree = 1 + static_cast<int>(generator() % 14);
        for (int j = 0; j < degree; ++j)
            roots.emplace_back(Rational(static_cast<int>(generator() % 49) - 24) / 8,
                               Rational(static_cast<int>(generator() % 49) - 24) / 8);
        const Complex c{Rational(static_cast<int>(generator() % 9) - 4) / 4,
                        Rational(static_cast<int>(generator() % 9) - 4) / 4};
        const int exponent = static_cast<int>(generator() % 5) - 2;
        const Rational r = reference::power2(exponent);
        int exact_count = 0;
        bool boundary = false;
        for (const auto& z : roots) {
            const Rational distance_squared = (z.real - c.real) * (z.real - c.real)
                                            + (z.imag - c.imag) * (z.imag - c.imag);
            exact_count += distance_squared < r * r ? 1 : 0;
            boundary = boundary || distance_squared == r * r;
        }
        detail::ComplexPolynomial p;
        set_from_roots(p.value(), roots);
        detail::ComplexBall center;
        set(center.value(), c);
        detail::Float radius;
        arf_set_ui_2exp_si(radius.value(), 1, exponent);
        for (slong precision : {32, 64, 128, 256}) {
            const auto result = detail::pellet_count(p.value(), center.value(), radius.value(), precision, 6);
            if (result.count < 0) continue;
            require(!boundary && result.count == exact_count, "independent exact disk count disagrees");
            ++certificates;
            graeffe_certificates += result.graeffe_steps > 0 ? 1 : 0;
        }
    }
    require(certificates > 700 && graeffe_certificates > 100, "count tests did not exercise enough successful certificates");
    std::cout << certificates << " independently checked disk counts, " << graeffe_certificates << " using Graeffe\n";
}

// The rotated real shift behind CountContext::real_shift's retry: every
// coefficient modulus equals that of the exact shift p(c + y), so the Pellet
// comparisons and counts are those of the direct shift. Checked against exact
// Gaussian-rational shifts at centres off both axes; counts with the retry
// enabled are checked against exact counts as in randomized_known_counts.
void off_axis_real_shift() {
    std::mt19937 generator(20260929);
    int checked = 0;
    for (int sample = 0; sample < 120; ++sample) {
        std::vector<Complex> roots;
        const int degree = 1 + static_cast<int>(generator() % 24);
        for (int j = 0; j < degree; ++j)
            roots.emplace_back(Rational(static_cast<int>(generator() % 49) - 24) / 8,
                               Rational(static_cast<int>(generator() % 49) - 24) / 8);
        Complex c{Rational(1 + static_cast<int>(generator() % 8)) / 4, Rational(1 + static_cast<int>(generator() % 8)) / 4};
        if (generator() % 2) c.real = -c.real;
        if (generator() % 2) c.imag = -c.imag;
        std::vector<Complex> exact(1, {1});
        for (const auto& root : roots) {
            std::vector<Complex> next(exact.size() + 1);
            for (size_t j = 0; j < exact.size(); ++j) {
                next[j] = next[j] + exact[j] * Complex(-root.real, -root.imag);
                next[j + 1] = next[j + 1] + exact[j];
            }
            exact = std::move(next);
        }
        const auto shifted = reference::shifted(exact, c);
        detail::ComplexPolynomial p, rotated;
        set_from_roots(p.value(), roots);
        detail::ComplexBall center;
        set(center.value(), c);
        detail::rotated_real_shift(rotated.value(), p.value(), center.value(), 256);
        require(acb_poly_degree(rotated.value()) == degree, "rotated shift changed the degree");
        for (int k = 0; k <= degree; ++k) {
            detail::RealBall modulus, exact_modulus;
            acb_abs(modulus.value(), rotated.value()->coeffs + k, 256);
            arb_sqr(modulus.value(), modulus.value(), 256);
            set(exact_modulus.value(), shifted[k].real * shifted[k].real + shifted[k].imag * shifted[k].imag);
            require(arb_overlaps(modulus.value(), exact_modulus.value()), "rotated shift changed a coefficient modulus");
            ++checked;
        }
    }
    certpolysolve::detail::SearchStatistics statistics{};
    int certificates = 0;
    for (int sample = 0; sample < 240; ++sample) {
        std::vector<Complex> roots;
        const int degree = 1 + static_cast<int>(generator() % 14);
        for (int j = 0; j < degree; ++j)
            roots.emplace_back(Rational(static_cast<int>(generator() % 49) - 24) / 8,
                               Rational(static_cast<int>(generator() % 49) - 24) / 8);
        const Complex c{Rational(static_cast<int>(generator() % 9) - 4) / 4,
                        Rational(static_cast<int>(generator() % 9) - 4) / 4};
        const int exponent = static_cast<int>(generator() % 5) - 2;
        const Rational r = reference::power2(exponent);
        int exact_count = 0;
        bool boundary = false;
        for (const auto& z : roots) {
            const Rational distance_squared = (z.real - c.real) * (z.real - c.real)
                                            + (z.imag - c.imag) * (z.imag - c.imag);
            exact_count += distance_squared < r * r ? 1 : 0;
            boundary = boundary || distance_squared == r * r;
        }
        detail::ComplexPolynomial p;
        set_from_roots(p.value(), roots);
        detail::ComplexBall center;
        set(center.value(), c);
        detail::Float radius;
        arf_set_ui_2exp_si(radius.value(), 1, exponent);
        for (slong precision : {16, 32, 64, 128}) {
            detail::CountContext retrying(statistics);
            retrying.real_shift = true;
            const auto result = detail::pellet_count(p.value(), center.value(), radius.value(), precision, 6, -1, &retrying);
            if (result.count < 0) continue;
            require(!boundary && result.count == exact_count, "count with the rotated retry disagrees with the exact count");
            ++certificates;
        }
    }
    require(checked > 1000 && certificates > 600, "rotated real shift checks too few");
    std::cout << checked << " rotated-shift coefficient moduli and " << certificates
              << " counts with the rotated retry checked\n";
}

void untrusted_hints() {
    detail::ComplexPolynomial p;
    set_from_roots(p.value(), {{-2}, {2}});
    detail::RootVector hints(2);
    acb_set_si(hints.data(), -2);
    acb_set_si(hints.data() + 1, -2); // Both hints point at the SAME root.
    std::vector<detail::Root> output(1);
    output[0].certificate = "unchanged";
    cps_options options;
    if (detail::recover_clusters(output, p.value(), hints.data(), 2, 1, options, 128)) {
        require(output.size() == 2 && output[1].count == 2
                && arb_contains_si(acb_realref(output[1].box.value()), -2)
                && arb_contains_si(acb_realref(output[1].box.value()), 2),
                "duplicated hints invented a local cover missing the distant root");
    } else require(output.size() == 1, "failed recovery changed output");
    require(output[0].certificate == "unchanged", "recovery overwrote existing output");
    acb_indeterminate(hints.data());
    require(!detail::recover_clusters(output, p.value(), hints.data(), 2, 1, options, 128), "nonfinite hints accepted");
    options.cancel = [](void *) { return 1; };
    acb_set_si(hints.data(), -2);
    acb_set_si(hints.data() + 1, 2);
    require(!detail::recover_clusters(output, p.value(), hints.data(), 2, 1, options, 128), "cancelled recovery ran");
}

void coverage_quality() {
    std::vector<detail::Root> small(1), large(1);
    small[0].count = large[0].count = 2;
    small[0].multiplicity = large[0].multiplicity = 0;
    arb_add_error_2exp_si(acb_realref(small[0].box.value()), 0);
    arb_add_error_2exp_si(acb_realref(large[0].box.value()), 100);
    require(!detail::prefer_cover(large, small, 128), "larger equally coarse cover replaced a useful bound");
    require(detail::prefer_cover(small, large, 128), "narrower coarse cover rejected");
    std::vector<detail::Root> separated(2);
    separated[0].count = separated[0].multiplicity = 1;
    separated[1].count = separated[1].multiplicity = 1;
    require(detail::prefer_cover(separated, small, 128), "separated root partition not preferred");
}

void recovery_precision_handoff() {
    // Well-separated roots of (x+4)^48-1. Expanded-coordinate cancellation
    // makes 170-bit local counts inconclusive even with excellent proposals.
    // Re-enclosing at 340 bits permits all singleton certificates.
    constexpr slong degree = 48;
    detail::IntegerPolynomial exact;
    fmpz_poly_set_coeff_si(exact.value(), degree, 1);
    detail::Integer four;
    fmpz_set_ui(four.value(), 4);
    fmpz_poly_taylor_shift(exact.value(), exact.value(), four.value());
    fmpz_sub_ui(exact.value()->coeffs, exact.value()->coeffs, 1);
    detail::ComplexPolynomial polynomial;
    detail::RootVector hints(degree), known(degree);
    detail::ComplexBall unit;
    acb_unit_root(unit.value(), static_cast<ulong>(degree), 512);
    for (slong j=0; j<degree; ++j) {
        // Start at a nonreal root: an integer center can make the translation
        // exact even at the lower precision, avoiding the intended obstruction.
        acb_pow_ui(known[j], unit.value(), static_cast<ulong>(j+1), 512);
        acb_sub_ui(known[j], known[j], 4, 512);
        acb_get_mid(hints[j], known[j]);
    }
    cps_options options;
    options.goal = CPS_ROOT_APPROXIMATION;
    options.max_work_bits = 680;
    certpolysolve::detail::SearchStatistics statistics{};
    detail::CountContext context(statistics, true);
    std::vector<detail::Root> output(1);
    output[0].certificate = "retained";
    acb_poly_set_fmpz_poly(polynomial.value(), exact.value(), 170);
    bool needs_precision = false;
    require(!detail::recover_clusters(output, polynomial.value(), hints.data(), degree, 1, options,
                                     170, &context, &needs_precision), "uncertain recovery unexpectedly succeeded");
    require(needs_precision, "counting uncertainty was not returned to the source owner");
    require(output.size()==1 && output[0].certificate=="retained", "arithmetic refusal changed a retained cover");
    require(statistics.disk_tests==1 && statistics.graeffe_steps==0,
            "arithmetic refusal repeated counting instead of requesting precision");
    acb_poly_set_fmpz_poly(polynomial.value(), exact.value(), 340);
    output.clear();
    require(detail::recover_clusters(output, polynomial.value(), hints.data(), degree, 1, options,
                                    340, &context, &needs_precision), "re-enclosed source did not recover its roots");
    require(!needs_precision && output.size()==static_cast<size_t>(degree) && detail::disjoint(output),
            "successful recovery retained stale arithmetic refusal");
    for (slong j=0; j<degree; ++j)
        require(std::count_if(output.begin(), output.end(), [&](const auto &root) {
            return root.count==1 && acb_contains(root.box.value(), known[j]);
        })==1, "analytic translated roots-of-unity cover disagrees");

    // At a caller's precision ceiling, uncertainty may be intrinsic to the
    // coefficient family. Radius growth still proves its two-root cluster.
    acb_poly_zero(polynomial.value());
    detail::ComplexBall constant;
    acb_one(constant.value());
    arb_add_error_2exp_si(acb_realref(constant.value()), -160);
    acb_poly_set_coeff_acb(polynomial.value(), 0, constant.value());
    acb_poly_set_coeff_si(polynomial.value(), 1, -2);
    acb_poly_set_coeff_si(polynomial.value(), 2, 1);
    detail::RootVector pair(2);
    acb_one(pair[0]); acb_one(pair[1]);
    options.max_work_bits=170;
    detail::CountContext ceiling_context(statistics);
    output.clear();
    require(detail::recover_clusters(output, polynomial.value(), pair.data(), 2, 1, options,
                                    170, &ceiling_context, &needs_precision), "precision ceiling disabled geometric recovery");
    require(!needs_precision && output.size()==1 && output[0].count==2 && output[0].multiplicity==0,
            "uncertain pair lost its cluster semantics");
}

void compressed_exponent_certificates() {
    for (const char *encoded : {"100000000000000000000", "-100000000000000000000"}) {
        detail::Integer exponent, twice, linear_exponent, radius_exponent;
        require(fmpz_set_str(exponent.value(), encoded, 10) == 0, "test exponent parsing failed");
        fmpz_mul_ui(twice.value(), exponent.value(), 2);
        fmpz_add_ui(linear_exponent.value(), exponent.value(), 1);
        fmpz_sub_ui(radius_exponent.value(), exponent.value(), 80);
        detail::ComplexPolynomial p;
        detail::ComplexBall coefficient, center;
        detail::Number("1@" + detail::integer_string(twice.value())).enclose(acb_realref(coefficient.value()), 128);
        acb_poly_set_coeff_acb(p.value(), 0, coefficient.value());
        detail::Number("-1@" + detail::integer_string(linear_exponent.value())).enclose(acb_realref(coefficient.value()), 128);
        acb_poly_set_coeff_acb(p.value(), 1, coefficient.value());
        acb_poly_set_coeff_si(p.value(), 2, 1);
        detail::Number(std::string("1@") + encoded).enclose(acb_realref(center.value()), 128);
        detail::RealBall radius;
        detail::Number("1@" + detail::integer_string(radius_exponent.value())).enclose(radius.value(), 128);
        detail::Root root;
        require(detail::certify_cluster_box(root, p.value(), center.value(), arb_midref(radius.value()), 2, 1, 128),
                "local certificate failed for a compressed exponent beyond machine-integer range");
        require(root.count == 2 && root.multiplicity == 0, "compressed cluster count semantics failed");
        for (const auto& endpoint : root.endpoints)
            require(endpoint.mantissa_text().size() < 100 && endpoint.exponent_text().size() < 30,
                    "local certificate expanded a compressed exponent");
    }
}
}

int main() {
    int status = 0;
    try {
        predicate_examples();
        bounded_tests_resume();
        priced_taylor_resumes();
        cancelled_tests_stop();
        stop_covers_stop_for_cancellation();
        rectangle_sandwich();
        randomized_known_counts();
        off_axis_real_shift();
        untrusted_hints();
        coverage_quality();
        recovery_precision_handoff();
        compressed_exponent_certificates();
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        status = 1;
    }
    flint_cleanup();
    return status;
}
