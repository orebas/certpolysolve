// The pair tests of covers (disjoint, overlapping_roots, reflection_met) against every pair compared directly, on
// random covers of every size below, at and above the 128 regions from which the real-axis sweep used to start
// (since 25 October 2026 it sweeps at every size: dev/docs/easy-inputs-profile-20261025.md). Each answer depends
// only on which pairs meet, so the sweep must give exactly the brute-force answer.
#include "core/cover.hpp"
#include "certify/cover_sweep.hpp"
#include <algorithm>
#include <bit>
#include <cstdint>
#include <iostream>
#include <random>
#include <stdexcept>
#include <vector>

namespace d = certpolysolve::detail;

namespace {
void require(bool ok, const char *why) {
    if (!ok) throw std::runtime_error(why);
}

// A cover of n boxes: centres in the unit square, half-widths scaled so that a fraction of the pairs meet. With
// `mirrored`, the boxes come in conjugate pairs and some sit on the real axis, as a real polynomial's cover does.
std::vector<d::Root> cover(std::mt19937_64 &rng, size_t n, double spread, bool mirrored) {
    std::uniform_real_distribution<double> unit(-1.0, 1.0), width(0.0, 1.0);
    std::vector<d::Root> roots(n);
    for (size_t i = 0; i < n; ++i) {
        double x = unit(rng), y = unit(rng);
        if (mirrored && i % 2 == 1) {
            x = arf_get_d(arb_midref(acb_realref(roots[i - 1].box.value())), ARF_RND_NEAR);
            y = -arf_get_d(arb_midref(acb_imagref(roots[i - 1].box.value())), ARF_RND_NEAR);
        } else if (mirrored && i % 5 == 0) {
            y = 0.0;
        }
        acb_set_d_d(roots[i].box.value(), x, y);
        mag_t rx, ry;
        mag_init(rx);
        mag_init(ry);
        mag_set_d(rx, spread * width(rng));
        mag_set_d(ry, spread * width(rng));
        mag_add(arb_radref(acb_realref(roots[i].box.value())), arb_radref(acb_realref(roots[i].box.value())), rx);
        mag_add(arb_radref(acb_imagref(roots[i].box.value())), arb_radref(acb_imagref(roots[i].box.value())), ry);
        mag_clear(rx);
        mag_clear(ry);
    }
    return roots;
}

bool reflection_meets(const d::Root &a, const d::Root &b) {
    acb_t conjugate;
    acb_init(conjugate);
    acb_conj(conjugate, a.box.value());
    const bool meet = acb_overlaps(conjugate, b.box.value());
    acb_clear(conjugate);
    return meet;
}

void check(const std::vector<d::Root> &roots) {
    const size_t n = roots.size();
    std::vector<bool> overlap(n, false), met(n, false);
    bool any = false;
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < i; ++j) {
            if (acb_overlaps(roots[i].box.value(), roots[j].box.value())) {
                overlap[i] = overlap[j] = true;
                any = true;
            }
            if (reflection_meets(roots[i], roots[j]))
                met[i] = met[j] = true;
        }
    require(d::disjoint(roots) == !any, "disjoint differs from the pairs compared directly");
    require(d::disjoint(roots, 128) == !any, "rectangle sweep differs from the exact pair oracle");
    require(d::overlapping_roots(roots) == overlap, "overlapping_roots differs from the pairs compared directly");
    require(d::reflection_met(roots) == met, "reflection_met differs from the pairs compared directly");
    require(d::reflection_met(roots, 128) == met, "reflection sweep differs from the exact pair oracle");
}

d::Root rectangle(slong x, slong y, ulong rx, ulong ry) {
    d::Root root;
    acb_set_si_si(root.box.value(), x, y);
    mag_set_ui(arb_radref(acb_realref(root.box.value())), rx);
    mag_set_ui(arb_radref(acb_imagref(root.box.value())), ry);
    return root;
}

void boundaries_and_actual_boxes() {
    // Edge/corner contact is overlap; a unit gap is disjoint. Include points,
    // real/imaginary segments, repeated starts and zero-width real projections.
    std::vector<d::Root> boxes;
    boxes.push_back(rectangle(0, 0, 1, 1));
    boxes.push_back(rectangle(2, 0, 1, 1));
    boxes.push_back(rectangle(2, 2, 1, 1));
    boxes.push_back(rectangle(5, 0, 1, 1));
    boxes.push_back(rectangle(0, 4, 0, 1));
    boxes.push_back(rectangle(0, 5, 0, 0));
    boxes.push_back(rectangle(0, 6, 0, 0));
    boxes.push_back(rectangle(8, 0, 1, 0));
    for (size_t shift = 0; shift < boxes.size(); ++shift) {
        check(boxes);
        std::rotate(boxes.begin(), boxes.begin() + 1, boxes.end());
    }
    std::vector<d::Root> points;
    points.push_back(rectangle(0, 3, 0, 0));
    points.push_back(rectangle(0, -3, 0, 0));
    points.push_back(rectangle(4, 0, 1, 0));
    check(points); // every real event is simultaneous for the conjugate point pair
    for (auto &root : points) {
        root.multiplicity = 1;
        root.certificate = "test";
        root.publish(128);
    }
    d::prove_real_by_reflection(points, 128);
    require(points[0].realness == CPS_REAL_UNKNOWN && points[1].realness == CPS_REAL_UNKNOWN &&
                points[2].realness == CPS_REAL && points[2].certificate == "test+conjugation",
            "realness and certificate records must retain reflection semantics");
    std::vector<d::Root> cache(2);
    acb_set_d_d(cache[0].box.value(), 1.1, 0);
    mag_set_d(arb_radref(acb_realref(cache[0].box.value())), 1.0 / 3);
    cache[0].publish(16);
    d::Float radius, actual_upper;
    const auto real = acb_realref(cache[0].box.value());
    arf_set_mag(radius.value(), arb_radref(real));
    arf_add(actual_upper.value(), arb_midref(real), radius.value(), ARF_PREC_EXACT, ARF_RND_NEAR);
    require(arf_cmp(actual_upper.value(), cache[0].endpoints[1].value()) > 0,
            "fixture must have an actual box wider than the cached exported ends");
    arf_add(arb_midref(acb_realref(cache[1].box.value())), actual_upper.value(),
            cache[0].endpoints[1].value(), ARF_PREC_EXACT, ARF_RND_NEAR);
    arb_mul_2exp_si(acb_realref(cache[1].box.value()), acb_realref(cache[1].box.value()), -1);
    require(acb_overlaps(cache[0].box.value(), cache[1].box.value()), "point lies in the actual box");
    require(d::CoverSweep::prepare(cache, 16).has_value(), "actual-box check must use the sweep");
    check(cache);
    // Mutation without publication must likewise use the actual box.
    mag_set_ui(arb_radref(acb_realref(cache[0].box.value())), 4);
    check(cache);
    std::vector<d::Root> bad(1);
    acb_indeterminate(bad[0].box.value());
    require(!d::disjoint(bad, 128), "nonfinite disjointness must still be false");
}

void precision_and_compact_exponents() {
    // Native certificates are published at their generating P. The bound must
    // scale with that P, including precision well above the usual small cases.
    for (const slong precision : {8, 32, 85, 340, 4096, 32768}) {
        std::vector<d::Root> roots(2);
        acb_one(roots[0].box.value());
        mag_set_ui_2exp_si(arb_radref(acb_realref(roots[0].box.value())), 1, -precision - 8);
        acb_set_si(roots[1].box.value(), 2);
        for (auto &root : roots) root.publish(precision);
        auto sweep = d::CoverSweep::prepare(roots, precision);
        require(sweep && sweep->disjoint(), "published certificates must fit a P-dependent sweep bound");
    }
    // A billion-bit alignment is compact input, not permission to allocate it.
    std::vector<d::Root> compact(2);
    acb_one(compact[0].box.value());
    mag_set_ui_2exp_si(arb_radref(acb_realref(compact[0].box.value())), 1, -1'000'000'000);
    acb_set_si(compact[1].box.value(), 2);
    require(!d::CoverSweep::prepare(compact, 128), "unpublished huge alignment must retain the compact path");
    check(compact);
    compact[0].publish(128);
    require(d::CoverSweep::prepare(compact, 128).has_value(), "publication must restore the sweep bound");
    check(compact);

    // Huge exponent encodings are distinct from huge alignment spans. A scale
    // shared by midpoint and radius needs no exponent expansion.
    d::Integer exponent;
    fmpz_one(exponent.value());
    fmpz_mul_2exp(exponent.value(), exponent.value(), 100);
    for (const slong sign : {-1, 1}) {
        if (sign < 0) fmpz_neg(exponent.value(), exponent.value());
        std::vector<d::Root> roots;
        roots.push_back(rectangle(4, 0, 1, 1));
        roots.push_back(rectangle(8, 0, 1, 1));
        for (auto &root : roots) {
            arb_mul_2exp_fmpz(acb_realref(root.box.value()), acb_realref(root.box.value()), exponent.value());
            root.publish(85);
        }
        require(d::CoverSweep::prepare(roots, 85).has_value(), "large encoded exponents must remain compact");
        check(roots);
        if (sign < 0) fmpz_neg(exponent.value(), exponent.value());
    }
    // Stress the publication lemma with independent, widely separated center
    // and radius exponents, cancellation, and more input mantissa bits than P.
    std::mt19937_64 rng(20261106);
    const slong scales[] = {-1'000'000'000, -5000, -40, 0, 40, 5000, 1'000'000'000};
    const slong precisions[] = {8, 32, 85, 340};
    for (size_t trial = 0; trial < 256; ++trial) {
        const slong precision = precisions[trial % 4];
        std::vector<d::Root> roots(2);
        for (auto &root : roots)
            for (auto part : {acb_realref(root.box.value()), acb_imagref(root.box.value())}) {
                d::Integer mantissa;
                for (int word = 0; word < 4; ++word) {
                    fmpz_mul_2exp(mantissa.value(), mantissa.value(), 64);
                    fmpz_add_ui(mantissa.value(), mantissa.value(), static_cast<ulong>(rng()));
                }
                if (rng() % 2) fmpz_neg(mantissa.value(), mantissa.value());
                if (trial % 7 == 0) fmpz_zero(mantissa.value());
                arf_set_fmpz(arb_midref(part), mantissa.value());
                arf_mul_2exp_si(arb_midref(part), arb_midref(part), scales[rng() % 7]);
                if (trial % 11 != 0)
                    mag_set_ui_2exp_si(arb_radref(part), 1 + static_cast<ulong>(rng() % 127), scales[rng() % 7]);
            }
        for (auto &root : roots) root.publish(precision);
        require(d::CoverSweep::prepare(roots, precision).has_value(),
                "published native-shaped inputs must not use the quadratic fallback");
        check(roots);
    }
}

void growth_controls() {
    for (const size_t n : {128u, 512u, 2048u}) {
        std::vector<d::Root> stack;
        for (size_t i = 0; i < n; ++i)
            stack.push_back(rectangle(0, static_cast<slong>(4 * i), 1, 1));
        d::CoverSweepWork work;
        auto sweep = d::CoverSweep::prepare(stack, 128, &work);
        require(sweep.has_value(), "vertical stack must enter the rectangle sweep");
        require(sweep->disjoint(&work), "vertical stack is disjoint");
        const auto met = sweep->reflection_met(&work);
        require(std::none_of(met.begin(), met.end(), [](bool x) { return x; }), "self reflections do not mark");
        // Generous linearithmic budget: sorting plus the fixed number of
        // balanced-set queries/inserts/erases per event or newly marked ID.
        const uint64_t budget = 64 * n * (1 + std::bit_width(n));
        require(work.comparisons < budget, "vertical-stack comparison count exceeds linearithmic budget");
        std::cout << "stack n=" << n << " comparisons=" << work.comparisons << '\n';

        // n/2 upper vertical bars and n/2 lower horizontal bars. The original
        // cover is disjoint, but reflection creates (n/2)^2 cross edges.
        const size_t m = n / 2;
        std::vector<d::Root> crossing;
        for (size_t i = 0; i < m; ++i)
            crossing.push_back(rectangle(static_cast<slong>(4 * i), static_cast<slong>(2 * m + 1), 1, 2 * m));
        for (size_t j = 0; j < m; ++j)
            crossing.push_back(rectangle(static_cast<slong>(2 * (m - 1)), -static_cast<slong>(4 * j + 2),
                                         2 * (m - 1) + 1, 1));
        work = {};
        sweep = d::CoverSweep::prepare(crossing, 128, &work);
        require(sweep.has_value(), "cross-edge control must enter the sweep");
        require(sweep->disjoint(&work), "bar cover must be internally disjoint");
        const auto flags = sweep->reflection_met(&work);
        require(std::all_of(flags.begin(), flags.end(), [](bool x) { return x; }), "every bar has another mirror");
        require(work.marked == n, "each original ID must be marked exactly once");
        require(work.comparisons < budget, "reflection work must not enumerate its quadratic cross edges");
        std::cout << "cross n=" << n << " comparisons=" << work.comparisons << " marks=" << work.marked << '\n';
        if (n == 128) {
            check(stack);
            check(crossing);
            std::reverse(crossing.begin(), crossing.end());
            check(crossing);
        }
    }
}

void mixed_precision_ranks() {
    // Four conjugate point pairs differ only far down their real mantissas.
    // The other pairs have short integer coordinates. Equal-value ranks must
    // preserve conjugacy, while tiny exact gaps must preserve disjointness.
    constexpr size_t count = 128;
    for (const slong precision : {256, 4096, 32768}) {
        std::vector<d::Root> roots;
        for (size_t i = 0; i < count; ++i)
            roots.push_back(rectangle(static_cast<slong>(4 * (i / 2 + 1)), i % 2 ? -1 : 1, 0, 0));
        d::Float step, tiny;
        arf_one(step.value());
        arf_mul_2exp_si(step.value(), step.value(), -precision);
        arf_mul_2exp_si(tiny.value(), step.value(), -precision);
        for (size_t i = 0; i < 8; ++i) {
            auto real = arb_midref(acb_realref(roots[i].box.value()));
            arf_one(real);
            if (i / 2 < 2) arf_add(real, real, step.value(), ARF_PREC_EXACT, ARF_RND_NEAR);
            if (i / 2 == 1) arf_add(real, real, tiny.value(), ARF_PREC_EXACT, ARF_RND_NEAR);
            if (i / 2 == 3) arf_sub(real, real, step.value(), ARF_PREC_EXACT, ARF_RND_NEAR);
        }
        for (auto &root : roots) root.publish(2 * precision);
        std::mt19937_64 rng(static_cast<uint64_t>(precision));
        std::shuffle(roots.begin(), roots.end(), rng);
        uint64_t input_bits = 0;
        for (const auto &root : roots) {
            const auto real = arb_midref(acb_realref(root.box.value()));
            const auto imaginary = arb_midref(acb_imagref(root.box.value()));
            // Points contribute two equal X ends and four original/reflected Y ends.
            input_bits += 2 * static_cast<uint64_t>(std::max<slong>(1, arf_bits(real))) +
                          4 * static_cast<uint64_t>(std::max<slong>(1, arf_bits(imaginary)));
        }
        d::CoverSweepWork work;
        auto sweep = d::CoverSweep::prepare(roots, 2 * precision, &work);
        require(sweep.has_value(), "mixed-precision points must enter the sweep");
        require(work.endpoint_mantissa_bits == input_bits, "ranking must charge actual endpoint lengths");
        // Each merge level charges each emitted item once; adjacent
        // deduplication adds one more level. The larger list has 4*count ends.
        const uint64_t rank_budget = (1 + std::bit_width(4 * count - 1)) * input_bits;
        require(work.rank_mantissa_bits <= rank_budget, "rank formation exceeds its summed-size budget");
        const auto rank_comparisons = work.rank_comparisons;
        require(sweep->disjoint(&work), "exact gaps between long mantissas must stay disjoint");
        const auto met = sweep->reflection_met(&work);
        require(std::all_of(met.begin(), met.end(), [](bool x) { return x; }),
                "equal-value ranks must preserve every conjugate partner");
        require(work.rank_comparisons == rank_comparisons, "sweep queries must use ranks only");
        check(roots); // independent all-pairs Arb oracle, including the published boxes
        std::cout << "mixed precision=" << precision << " endpoint_bits=" << input_bits
                  << " rank_bit_charge=" << work.rank_mantissa_bits << '\n';
    }
}
} // namespace

int main() {
    try {
        boundaries_and_actual_boxes();
        precision_and_compact_exponents();
        growth_controls();
        mixed_precision_ranks();
        std::mt19937_64 rng(20261025);
        // Sizes around the old threshold (128) and small covers, where the sweep is new.
        const size_t sizes[] = {0, 1, 2, 3, 4, 5, 8, 16, 33, 64, 127, 128, 129, 200};
        // Spreads from covers that rarely meet to covers that mostly do.
        const double spreads[] = {0.0, 1e-3, 1e-2, 1e-1, 0.5};
        size_t covers = 0;
        for (const size_t n : sizes)
            for (const double spread : spreads)
                for (const bool mirrored : {false, true})
                    for (int repeat = 0; repeat < 8; ++repeat) {
                        check(cover(rng, n, spread, mirrored));
                        ++covers;
                    }
        // Since 26 October 2026 the sweep runs on binary64 bounds where they are finite: covers that test where they
        // stop deciding. Clusters whose boxes lie 2^-58 apart near 1 (apart exactly, not in binary64), boxes that
        // share an edge exactly (they meet), and covers with one box beyond binary64's range (the 64-bit sweep).
        for (int repeat = 0; repeat < 16; ++repeat) {
            std::vector<d::Root> tight(12);
            for (size_t k = 0; k < tight.size(); ++k) {
                arb_set_si(acb_realref(tight[k].box.value()), 1);
                arb_t step;
                arb_init(step);
                arb_set_si(step, static_cast<slong>(4 * k + (repeat % 3)));
                arb_mul_2exp_si(step, step, -58);
                arb_add(acb_realref(tight[k].box.value()), acb_realref(tight[k].box.value()), step, 128);
                arb_clear(step);
                mag_set_ui_2exp_si(arb_radref(acb_realref(tight[k].box.value())), 1, -58 - (repeat % 2));
                arb_set_si(acb_imagref(tight[k].box.value()), static_cast<slong>(k % 2) - (repeat % 3 == 0));
                mag_set_ui_2exp_si(arb_radref(acb_imagref(tight[k].box.value())), 1, -1);
            }
            check(tight);
            ++covers;
            std::vector<d::Root> edges(6);
            for (size_t k = 0; k < edges.size(); ++k) {
                acb_set_si_si(edges[k].box.value(), static_cast<slong>(2 * k), static_cast<slong>(repeat % 4));
                mag_set_ui(arb_radref(acb_realref(edges[k].box.value())), 1);
                mag_set_ui(arb_radref(acb_imagref(edges[k].box.value())), 1);
            }
            check(edges);
            ++covers;
            auto wide = cover(rng, 40, 1e-2, repeat % 2 == 1);
            arb_mul_2exp_si(acb_realref(wide[repeat % 40].box.value()), acb_realref(wide[repeat % 40].box.value()), 5000);
            check(wide);
            ++covers;
        }
        std::cout << covers << " random covers: disjoint, overlapping_roots and reflection_met agree with every pair "
                     "compared\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << "test_cover_sweep: " << e.what() << '\n';
        return 1;
    }
}
