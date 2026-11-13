// The certificate attempt's binary64 Aberth packet with the interaction tree
// (machine_candidates.cpp, from heuristics::machine_interaction_tree_degree
// on). The tree only changes how the sums of 1/(z_i - z_j) are formed: far
// fields from multipole series, the nearby points directly, and the tree
// follows every point the sweep moves. Its sums must equal the direct sums to
// within binary64 rounding after points move; the proposals must converge
// where the direct sums take them, and the packet must do less work.
#include "internal.hpp"
#include <numbers>
#include <certpolysolve/certpolysolve.hpp>
#include "numeric/candidates/candidates.hpp"
#include "core/heuristics.hpp"
#include "numeric/candidates/interaction_tree.hpp"
#include "numeric/candidates/machine_candidates.hpp"
#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <limits>
#include <set>
#include <utility>
#include <vector>
namespace d = certpolysolve::detail;
static void require(bool value, const char *why) {
    if (!value) throw std::runtime_error(why);
}
// The direct sum of 1/(z_i - z_j) over j != i, and the sum of the terms' sizes, the references of the tests below:
// each term in binary64, the sums with Neumaier's compensation, so that summing 4,095 terms adds about one rounding
// to the terms' own few. Not in long double, which is binary64 under MSVC and on macOS arm64: there a direct
// long double sum's own rounding was of the size the tests bound (the Windows session's finding, 21 October 2026).
struct DirectSum {
    std::complex<double> sum;
    double size = 0;
};
static DirectSum direct_sum(const std::vector<std::complex<double>> &points, size_t i) {
    double re = 0, re_lost = 0, im = 0, im_lost = 0, size = 0, size_lost = 0;
    const auto add = [](double &total, double &lost, double x) {
        const double t = total + x;
        lost += std::abs(total) >= std::abs(x) ? (total - t) + x : (x - t) + total;
        total = t;
    };
    for (size_t j = 0; j < points.size(); ++j)
        if (j != i) {
            const std::complex<double> difference = points[i] - points[j];
            const std::complex<double> term = 1.0 / difference;
            add(re, re_lost, term.real());
            add(im, im_lost, term.imag());
            add(size, size_lost, 1.0 / std::abs(difference));
        }
    return {{re + re_lost, im + im_lost}, size + size_lost};
}
// The largest error, over every target, of the tree's sum of 1/(z_i - z_j)
// against the direct sum, in unit roundoffs of the sum of the terms'
// sizes, after a quarter of the points moved by up to 2^-11. `follow`: the tree
// is told of the moves (move()); otherwise it keeps their old positions.
static double worst_error(const std::vector<std::complex<double>> &start, bool follow) {
    auto points = start;
    d::MachineInteractionTree tree(points);
    uint64_t state = 5;
    const auto unit = [&] {
        state = state * 6364136223846793005ULL + 1442695040888963407ULL;
        return static_cast<double>(state >> 11) * 0x1p-53 - 0.5;
    };
    for (size_t k = 0; k < points.size(); k += 4) {
        const std::complex<double> to = points[k] + std::complex<double>(0x1p-10 * unit(), 0x1p-10 * unit());
        if (follow)
            tree.move(k, to);
        points[k] = to;
    }
    double worst = 0;
    std::vector<size_t> near;
    for (size_t i = 0; i < points.size(); ++i) {
        std::complex<double> far;
        double tail = 0;
        require(tree.evaluate(i, far, near, tail), "the tree refused a target");
        std::complex<double> sum = far;
        for (const size_t j : near)
            sum += 1.0 / (points[i] - points[j]);
        const auto reference = direct_sum(points, i);
        worst = std::max(worst, std::abs(sum - reference.sum) / (reference.size * std::numeric_limits<double>::epsilon() / 2));
    }
    return worst;
}
// The program route's binary64 stage sums through a coarser tree
// (heuristics::ProgramFarField: far fields from a ratio of 0.5, series of 14
// terms), whose truncation lies far above binary64's rounding. Over every
// target: the largest error against the direct sum, and the largest
// amount by which it passes the tail the tree reports (its bound on what its
// series leave), both in unit roundoffs of the sum of the terms' sizes.
static std::pair<double, double> program_tree_errors(const std::vector<std::complex<double>> &points) {
    const d::InteractionTreeOf<d::heuristics::program_tree_order, d::heuristics::ProgramFarField> tree(points);
    double worst = 0, excess = -INFINITY;
    std::vector<size_t> near;
    for (size_t i = 0; i < points.size(); ++i) {
        std::complex<double> far;
        double tail = 0;
        require(tree.evaluate(i, far, near, tail), "the program tree refused a target");
        std::complex<double> sum = far;
        for (const size_t j : near)
            sum += 1.0 / (points[i] - points[j]);
        const auto reference = direct_sum(points, i);
        const double unit = reference.size * std::numeric_limits<double>::epsilon() / 2;
        const double error = std::abs(sum - reference.sum);
        worst = std::max(worst, error / unit);
        excess = std::max(excess, (error - tail) / unit);
    }
    return {worst, excess};
}
int main() {
    certpolysolve::ThreadScope lifetime;
    try {
        // 4,096 points: 64 clusters of 64, each of radius 1/64 about a centre
        // of a grid of spacing 1/8, so that clusters are far fields near the
        // acceptance ratio. Followed, the worst error is 5.1 roundings against
        // the compensated direct sum (6 against the long double one it
        // replaced; a plain direct sum's own error is of that size); not
        // followed, the moves stay in the far fields: 3.1e13 roundings.
        {
            std::vector<std::complex<double>> points(4096);
            uint64_t state = 7;
            const auto unit = [&] {
                state = state * 6364136223846793005ULL + 1442695040888963407ULL;
                return static_cast<double>(state >> 11) * 0x1p-53;
            };
            for (size_t k = 0; k < points.size(); ++k) {
                const size_t c = k / 64;
                const std::complex<double> centre(static_cast<double>(c % 8) / 8, static_cast<double>(c / 8) / 8);
                points[k] = centre + std::polar(unit() / 64, 2 * std::numbers::pi * unit());
            }
            const double followed = worst_error(points, true), unfollowed = worst_error(points, false);
            std::cout << "machine tree: worst error " << followed << " roundings with the moves followed, " << unfollowed
                      << " without\n";
            require(followed <= 16, "the tree's sums left binary64 rounding after points moved");
            require(unfollowed > 1e6, "the moves did not reach the far fields (check the fixture)");
            // The program stage's tree on the same points: its far fields are
            // the coarser rule's (an error of 1.3e10 roundings, about 1.4e-6
            // of the terms' sizes, against the binary64 tree's 6), and every
            // sum stays within the tail it reports (by 1e9 roundings or more).
            const auto [worst, excess] = program_tree_errors(points);
            require(worst > 1e6, "the program tree's sums are binary64-accurate: its coarser rule is not in force");
            require(excess <= 16, "a program tree's sum left the tail it reported");
            std::cout << "program tree: worst error " << worst << " roundings, " << excess
                      << " roundings past its tail at most\n";
        }
        // A sweep builds the tree only when that costs less than the direct
        // sums at the points still moving (machine_candidates.cpp). Every
        // point at a root of x^4096 - 1 but three, moved off theirs by
        // 2^-10: the first sweep builds it (no history yet) and the other
        // points stop at once; the second has three points moving, whose
        // direct sums (3 x 4,095 pairs) cost far less than a construction
        // over 4,096 points.
        {
            constexpr slong n = 4096;
            static_assert(n >= static_cast<slong>(d::heuristics::machine_interaction_tree_degree));
            d::ComplexPolynomial p;
            acb_poly_set_coeff_si(p.value(), 0, -1);
            acb_poly_set_coeff_si(p.value(), n, 1);
            d::RootVector start(n);
            std::vector<std::complex<double>> positions(n);
            for (slong k = 0; k < n; ++k) {
                const bool moved = k % 1366 == 0; // k = 0, 1,366 and 2,732
                positions[static_cast<size_t>(k)] =
                    std::polar(moved ? 1 + 0x1p-10 : 1.0, 2 * std::numbers::pi * static_cast<double>(k) / static_cast<double>(n));
                acb_set_d_d(start.data() + k, positions[static_cast<size_t>(k)].real(),
                            positions[static_cast<size_t>(k)].imag());
            }
            auto engine = d::prepare_machine_engine(d::IterationMethod::Aberth, d::ArithmeticTier::Binary64, p.value(),
                                                    start.view(), false, false, true, UINT64_MAX, false, true);
            require(engine != nullptr, "the binary64 engine was refused");
            cps_options options;
            certpolysolve::detail::CandidateStatistics first{}, second{};
            const auto one = engine->advance(1, options, first);
            const auto two = engine->advance(1, options, second);
            require(!one.arithmetic_failure && !two.arithmetic_failure && one.active == 3,
                    "the fixture no longer stops every point but three in its first sweep");
            const d::MachineInteractionTree construction(positions);
            require(second.root_interactions < construction.build_terms(),
                    "a sweep with three points moving rebuilt the interaction tree");
        }
        // A sparse polynomial of degree 8,192 with 64 small integer
        // coefficients: its roots lie near the unit circle and binary64
        // settles every one of them.
        constexpr slong degree = 8192, count = 64;
        static_assert(degree >= static_cast<slong>(d::heuristics::machine_interaction_tree_degree));
        uint64_t state = 11;
        const auto next = [&] {
            state = state * 6364136223846793005ULL + 1442695040888963407ULL;
            return state >> 33;
        };
        std::set<slong> support{0, degree};
        while (static_cast<slong>(support.size()) < count)
            support.insert(1 + static_cast<slong>(next() % (degree - 1)));
        d::ComplexPolynomial p;
        for (const slong e : support) {
            const auto c = static_cast<slong>(next() % 19) - 9;
            acb_poly_set_coeff_si(p.value(), e, c ? c : 1);
        }
        d::RootVector initial(degree);
        require(d::initial_candidates(initial.data(), p.value(), 64), "no initial points");
        cps_options options;
        d::RootVector tree_points(degree), direct_points(degree);
        uint64_t tree_work = 0, direct_work = 0;
        for (const bool tree : {true, false}) {
            certpolysolve::detail::CandidateStatistics statistics{};
            auto engine = d::prepare_machine_engine(d::IterationMethod::Aberth, d::ArithmeticTier::Binary64, p.value(),
                                                    initial.view(), false, false, true, UINT64_MAX, true, tree);
            require(engine != nullptr, "the binary64 engine was refused");
            const auto step = engine->advance(d::heuristics::machine_cover_sweeps, options, statistics);
            require(!step.arithmetic_failure && step.active == 0 && step.candidates.converged,
                    "the packet did not settle every point");
            engine->export_candidates((tree ? tree_points : direct_points).data());
            (tree ? tree_work : direct_work) = statistics.root_interactions;
        }
        // The same roots: every proposal of one packet has a proposal of the
        // other within 2^-40 of its size, and no two share one.
        std::vector<std::pair<double, double>> direct(degree);
        for (slong i = 0; i < degree; ++i)
            direct[static_cast<size_t>(i)] = {arf_get_d(arb_midref(acb_realref(direct_points[i])), ARF_RND_NEAR),
                                              arf_get_d(arb_midref(acb_imagref(direct_points[i])), ARF_RND_NEAR)};
        std::sort(direct.begin(), direct.end());
        std::vector<uint8_t> used(direct.size(), 0);
        for (slong i = 0; i < degree; ++i) {
            const double x = arf_get_d(arb_midref(acb_realref(tree_points[i])), ARF_RND_NEAR);
            const double y = arf_get_d(arb_midref(acb_imagref(tree_points[i])), ARF_RND_NEAR);
            const double tolerance = std::ldexp(std::hypot(x, y), -40);
            auto j = static_cast<size_t>(std::lower_bound(direct.begin(), direct.end(), std::pair{x - tolerance, -INFINITY}) -
                                         direct.begin());
            bool matched = false;
            for (; j < direct.size() && direct[j].first <= x + tolerance && !matched; ++j)
                if (!used[j] && std::hypot(direct[j].first - x, direct[j].second - y) <= tolerance)
                    used[j] = matched = true;
            require(matched, "the tree's packet converged to a different point");
        }
        // The direct packet sums degree - 1 pairs at every point and sweep;
        // the tree's work (its pairs, its series' terms and its
        // construction) is a fraction of that at this degree.
        require(tree_work * 2 < direct_work, "the tree did not halve the packet's interaction work");
        std::cout << "machine tree test passed (interaction work " << tree_work << " against " << direct_work << ")\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "FAILED: " << error.what() << '\n';
        return 1;
    }
}
