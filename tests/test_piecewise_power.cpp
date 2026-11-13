#include "solve/prepared_piecewise.hpp"
#include "numeric/arithmetic/polynomial_arithmetic.hpp"
#include "exact_winding_api.hpp"
#include "methods/piecewise/piecewise_roots.hpp"
#include "numeric/exact/power_composition.hpp"
#include <iostream>

namespace d = certpolysolve::detail;
using namespace reference;

static void check_cover(const cps_result &result, const std::vector<Complex> &known,
                        const cps_options &options, bool complete) {
    require(d::disjoint(result.roots), "power images overlap");
    uint64_t covered = 0;
    for (const auto &root : result.roots) {
        const auto endpoint = [&](size_t i) {
            return reference::endpoint(certpolysolve::Dyadic{root.endpoints[i].mantissa_text(), root.endpoints[i].exponent_text()});
        };
        const auto left = endpoint(0), right = endpoint(1), bottom = endpoint(2), top = endpoint(3);
        uint64_t count = 0;
        for (const auto &z : known)
            count += left <= z.real && z.real <= right && bottom <= z.imag && z.imag <= top;
        require(count == root.count, "mapped region miscounted independently planted roots");
        if (root.multiplicity) require(count == root.multiplicity, "mapped multiplicity is wrong");
        if (complete) require(d::meets_goal(root, options, static_cast<slong>(result.work_bits)),
                              "mapped region lost original-coordinate accuracy");
        covered += count;
    }
    if (complete) require(result.complete && result.goal_met && covered == known.size(),
                           "mapped cover is incomplete");
}

static std::vector<Complex> orbit(const Rational &radius) {
    return {{radius}, {-radius}, {0, radius}, {0, -radius}};
}

static void solve_known(const std::vector<Complex> &known, const char *scale = "1") {
    const auto original = reference::input(from_roots(known));
    certpolysolve::detail::DensePolynomial input = certpolysolve::detail::dense(*original.get());
    const d::Number multiplier(scale);
    // Scaling is encoded without expanding its potentially million-bit exponent.
    for (auto &coefficient : input.coefficients) {
        coefficient.real.exponent = multiplier.exponent;
        coefficient.imag.exponent = multiplier.exponent;
    }
    auto config = d::maintained_piecewise_config();
    config.power_composition = true;
    for (const auto goal : {CPS_ROOT_APPROXIMATION, CPS_ROOT_ACCURACY}) {
        d::PreparedPiecewiseState state(input, config);
        cps_options options;
        options.goal = goal;
        options.max_work_bits = 1024;
        options.work_units = 1;
        cps_result result;
        for (unsigned step = 0; step < 1000; ++step) {
            result = state.advance(options);
            require(result.progress.work_units_this_call <= 1, "power preparation exceeded work allowance");
            check_cover(result, known, options, result.goal_met);
            if (result.goal_met) break;
        }
        check_cover(result, known, options, true);
        const auto prepared = state.prepared_structure();
        require(prepared && prepared->accounted_degree() == known.size() &&
                std::any_of(prepared->factors.begin(), prepared->factors.end(),
                    [](const auto &factor) { return factor.power > 1; }), "power identity was not retained");
        // Refine after a completed query; preserve both the exact zero count
        // and the original power map, rather than applying the map twice.
        options.bits = 170;
        options.work_units = 0;
        const auto old = result;
        result = state.advance(options);
        require(result.work_bits == old.work_bits && result.roots.size() == old.roots.size(),
                "zero-work query recomputed mapped roots");
        for (size_t i = 0; i < old.roots.size(); ++i)
            require(acb_equal(result.roots[i].box.value(), old.roots[i].box.value()),
                    "zero-work query changed the retained cover");
        options.work_units = 100;
        options.cancel = [](void *) { return 1; };
        result = state.advance(options);
        require((result.goal_met || result.termination == CPS_CANCELLED) &&
                    result.progress.work_units_this_call == 0 && result.roots.size() == old.roots.size(),
                "cancelled tightening discarded the retained power cover");
        for (size_t i = 0; i < old.roots.size(); ++i)
            require(acb_equal(result.roots[i].box.value(), old.roots[i].box.value()),
                    "cancelled tightening changed the retained power cover");
        options.cancel = nullptr;
        options.work_units = 1;
        for (unsigned step = 0; step < 1000; ++step) {
            result = state.advance(options);
            check_cover(result, known, options, result.goal_met);
            if (result.goal_met) break;
        }
        check_cover(result, known, options, true);
    }
}

static void cancelled_preparation() {
    auto known = orbit(1);
    auto second = orbit(2);
    known.insert(known.end(), second.begin(), second.end());
    auto input = reference::input(from_roots(known));
    auto config = d::maintained_piecewise_config();
    config.power_composition = true;
    for (unsigned checkpoint : {1u, 3u, 5u}) {
        struct Request { unsigned calls = 0, checkpoint; } request{0, checkpoint};
        cps_options options;
        options.goal = CPS_ROOT_APPROXIMATION;
        options.cancel_data = &request;
        options.cancel = [](void *p) -> int {
            auto &r = *static_cast<Request *>(p);
            return ++r.calls >= r.checkpoint;
        };
        d::PreparedPiecewiseState state(certpolysolve::detail::dense(*input.get()), config);
        const auto cancelled = state.advance(options);
        require(cancelled.termination == CPS_CANCELLED && cancelled.degree == 8,
                "power preparation cancellation lost input identity");
        check_cover(cancelled, known, options, false);
        options.cancel = nullptr;
        check_cover(state.advance(options), known, options, true);
    }
}

static void overlap_sweep() {
    // Compare against direct pairwise geometry, including the sweep threshold,
    // equal left endpoints, touching boxes and vertical stacks.
    for (size_t n : {size_t(17), size_t(128), size_t(257)}) {
        std::vector<d::Root> roots(n);
        for (size_t i = 0; i < n; ++i) {
            acb_set_si_si(roots[i].box.value(), static_cast<slong>(i % 7), static_cast<slong>(i / 7));
            if (i % 5 == 0) {
                arb_add_error_2exp_si(acb_realref(roots[i].box.value()), 0);
                arb_add_error_2exp_si(acb_imagref(roots[i].box.value()), 0);
            }
        }
        std::vector<bool> expected(n, false);
        for (size_t i = 0; i < n; ++i)
            for (size_t j = 0; j < i; ++j)
                if (acb_overlaps(roots[i].box.value(), roots[j].box.value())) expected[i] = expected[j] = true;
        require(d::overlapping_roots(roots) == expected && !d::disjoint(roots),
                "outward sweep missed a box intersection");
        // Which regions have another region on their reflection in the real
        // axis: the same sweep against the direct definition. The stacks
        // start at height zero, so boxes on the axis, boxes above it and
        // boxes that reach below it are all present.
        std::vector<bool> reflected(n, false);
        d::ComplexBall conjugate;
        for (size_t i = 0; i < n; ++i) {
            acb_conj(conjugate.value(), roots[i].box.value());
            for (size_t j = 0; j < n; ++j)
                if (i != j && acb_overlaps(conjugate.value(), roots[j].box.value())) reflected[i] = true;
        }
        require(d::reflection_met(roots) == reflected, "the reflection sweep disagrees with the pairwise definition");
    }
    d::ComplexBall conjugate_box;
    // A complete cover of a conjugation-invariant root set, as the driver's
    // final assembly sees it: real roots in boxes on the axis, conjugate
    // pairs in separate boxes, and pairs so close to the axis that one box
    // holds both. The assembly replaces a region whose reflection no other
    // region meets by its real segment. It used to look for the other region
    // with the boxes as they stood at that moment; it now decides from the
    // cover before any replacement. Both must make the same regions real.
    for (size_t n : {size_t(40), size_t(300)}) {
        std::vector<d::Root> cover;
        const auto add = [&](slong real, slong imaginary, slong radius_exponent, uint64_t count,
                             uint64_t multiplicity) {
            d::Root root;
            acb_set_si_si(root.box.value(), real, imaginary);
            arb_add_error_2exp_si(acb_realref(root.box.value()), radius_exponent);
            arb_add_error_2exp_si(acb_imagref(root.box.value()), radius_exponent);
            root.count = count;
            root.multiplicity = multiplicity;
            cover.push_back(std::move(root));
        };
        for (size_t k = 0; k < n; ++k) {
            const slong x = static_cast<slong>(4 * k);
            if (k % 4 == 0)
                add(x, 0, -2, 1, 1);                  // a real root
            else if (k % 4 == 1) {
                add(x, 5, -2, 1, 1);                  // a conjugate pair in two boxes,
                add(x, -5, -2, 1, 1);
                if (k % 8 == 1) add(x, 0, -3, 1, 1);  // sometimes with a real root between them
            } else if (k % 4 == 2)
                add(x, 0, -1, 2, 0);                  // a close pair in one box (two roots: not eligible)
            else if (k % 8 == 3)
                add(x, 0, -2, 2, 2);                  // a real double root (eligible: count = multiplicity)
            else {
                add(x, 3, -2, 2, 2);                  // a double root off the axis and its conjugate
                add(x, -3, -2, 2, 2);
            }
        }
        require(d::disjoint(cover), "the realness fixture's cover is not disjoint");
        // The sequential rule, as aggregate() had it.
        auto sequential = cover;
        for (size_t i = 0; i < sequential.size(); ++i) {
            auto &root = sequential[i];
            if (!(root.count == 1 || (root.multiplicity && root.count == root.multiplicity))) continue;
            acb_conj(conjugate_box.value(), root.box.value());
            bool other = false;
            for (size_t j = 0; j < sequential.size(); ++j)
                if (i != j && acb_overlaps(conjugate_box.value(), sequential[j].box.value())) {
                    other = true;
                    break;
                }
            if (!other) {
                arb_zero(acb_imagref(root.box.value()));
                root.realness = CPS_REAL;
            }
        }
        const auto met = d::reflection_met(cover);
        size_t real = 0;
        for (size_t i = 0; i < cover.size(); ++i) {
            const bool eligible = cover[i].count == 1 || (cover[i].multiplicity && cover[i].count == cover[i].multiplicity);
            const bool made_real = eligible && !met[i];
            require(made_real == (sequential[i].realness == CPS_REAL),
                    "deciding realness from the cover before any replacement differs from the sequential rule");
            real += made_real;
        }
        require(real != 0 && real != cover.size(), "the realness fixture has no real or no nonreal region");
        // A cover without overlaps is left as it is by the merge.
        auto merged = cover;
        d::merge_overlapping(merged, 64);
        require(merged.size() == cover.size(), "the merge changed a cover without overlaps");
        for (size_t i = 0; i < cover.size(); ++i)
            require(acb_equal(merged[i].box.value(), cover[i].box.value()) && merged[i].count == cover[i].count,
                    "the merge changed a region of a cover without overlaps");
        // With one overlap added, the two regions become one whose count is their sum.
        const uint64_t before = cover[0].count + cover[1].count;
        acb_union(merged[1].box.value(), merged[1].box.value(), merged[0].box.value(), 64);
        d::merge_overlapping(merged, 64);
        uint64_t total = 0, expected_total = 0;
        for (const auto &root : merged) total += root.count;
        for (const auto &root : cover) expected_total += root.count;
        require(d::disjoint(merged) && merged.size() < cover.size() && total == expected_total &&
                    merged[0].count >= before,
                "the merge did not join two overlapping regions");
    }
}

int main() {
    certpolysolve::ThreadScope lifetime;
    try {
        auto mixed = orbit(1);
        auto second = orbit(2);
        mixed.insert(mixed.end(), second.begin(), second.end());
        solve_known(mixed);
        solve_known(mixed, "1@-1000000");
        solve_known({Complex(1,1), Complex(1,-1), Complex(-1,1), Complex(-1,-1)});
        solve_known(orbit(power2(-200)));
        mixed.insert(mixed.end(), second.begin(), second.end());
        mixed.insert(mixed.end(), 5, Complex(0));
        solve_known(mixed);
        // The simple unit orbit is close to an unfinished repeated factor.
        // Its partial boxes must exclude q(z^k), not q(z), before publication.
        std::vector<Complex> close = orbit(1);
        for (const auto &r : {Rational(1) + power2(-160), Rational(2), Rational(3)}) {
            auto points = orbit(r);
            close.insert(close.end(), points.begin(), points.end());
            close.insert(close.end(), points.begin(), points.end());
        }
        solve_known(close);
        cancelled_preparation();
        overlap_sweep();
        std::cout << "piecewise power maps preserve counts, multiplicities, partial covers and continuation\n";
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
