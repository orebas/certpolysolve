// The node products for short nodes (node_products.hpp: binary64 parts, and parts of up to two limbs) against
// ball arithmetic at a far higher precision.
// A product becomes the denominator of a secular residue that certificates rest on, so the property tested is
// enclosure with the radius the header derives: the exact product within (n - 1) 2^-121 of the held one, and the
// separation a lower bound of the least distance that is not slack by more than a part in a million. The node
// sets are the cases of the derivation: well separated, clustered to the last bit, parts with exponents far
// apart, zero parts, sizes at both ends of the admitted range.
#include "numeric/evaluate/node_products.hpp"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
namespace d = certpolysolve::detail;
static void require(bool value, const std::string &why) {
    if (!value) throw std::runtime_error(why);
}
namespace {
struct Random {
    uint64_t state;
    uint64_t next() {
        state = state * 6364136223846793005ULL + 1442695040888963407ULL;
        return state >> 11;
    }
    double unit() { return static_cast<double>(next()) * 0x1p-53; } // [0, 1)
    double signed_unit() { return 2 * unit() - 1; }
};
std::vector<d::ComplexBall> balls(const std::vector<std::pair<double, double>> &points) {
    std::vector<d::ComplexBall> nodes(points.size());
    for (size_t j = 0; j < points.size(); ++j) acb_set_d_d(nodes[j].value(), points[j].first, points[j].second);
    return nodes;
}
void check(const std::vector<d::ComplexBall> &nodes, const std::string &what, bool binary64) {
    constexpr slong reference = 1024;
    const d::ShortNodeProducts products(nodes);
    require(products.usable(), what + ": the nodes were refused");
    require(products.binary64() == binary64, what + ": the nodes were not taken in the form expected");
    d::ComplexBall held, exact, factor, gap;
    d::Magnitude separation, least, least_upper, distance, allowed, error;
    for (size_t i = 0; i < nodes.size(); ++i) {
        require(products.product(held.value(), separation.value(), i), what + ": a product was refused");
        acb_one(exact.value());
        mag_inf(least.value());
        mag_inf(least_upper.value());
        for (size_t j = 0; j < nodes.size(); ++j) {
            if (j == i) continue;
            acb_sub(factor.value(), nodes[i].value(), nodes[j].value(), reference);
            acb_mul(exact.value(), exact.value(), factor.value(), reference);
            acb_get_mag_lower(distance.value(), factor.value());
            mag_min(least.value(), least.value(), distance.value());
            acb_get_mag(distance.value(), factor.value());
            mag_min(least_upper.value(), least_upper.value(), distance.value());
        }
        require(acb_contains(held.value(), exact.value()),
                what + ": the product of node " + std::to_string(i) + " does not enclose the exact one");
        // Its radius is the derived one: (n - 1) 2^-121 of its modulus, both coordinates.
        acb_get_mag(allowed.value(), exact.value());
        mag_mul_ui(allowed.value(), allowed.value(), static_cast<ulong>(nodes.size() - 1));
        mag_mul_2exp_si(allowed.value(), allowed.value(), d::ShortNodeProducts::step_error_log2 + 2);
        mag_max(error.value(), arb_radref(acb_realref(held.value())), arb_radref(acb_imagref(held.value())));
        require(mag_cmp(error.value(), allowed.value()) <= 0, what + ": a product is wider than the derived bound");
        require(acb_rel_accuracy_bits(held.value()) >= d::ShortNodeProducts::accuracy_bits(nodes.size()) - 2,
                what + ": a product is less accurate than accuracy_bits says");
        // The separation is below the least distance (against its upper bound; the lower one is a magnitude,
        // itself a part in 2^29 short), by no more than a part in a million.
        require(mag_cmp(separation.value(), least_upper.value()) <= 0,
                what + ": the separation exceeds the least distance");
        mag_mul_2exp_si(distance.value(), least.value(), -20);
        mag_add(allowed.value(), separation.value(), distance.value());
        require(mag_cmp(least.value(), allowed.value()) <= 0, what + ": the separation is slack");
    }
}
void check(const std::vector<std::pair<double, double>> &points, const std::string &what) {
    check(balls(points), what, true);
}
// A number of about 117 bits: a binary64 number and another one 2^-64 of its size added exactly; `low` zero
// leaves the binary64 number.
void long_number(arf_t out, double high, double low) {
    d::Float tail;
    arf_set_d(out, high);
    arf_set_d(tail.value(), low);
    arf_mul_2exp_si(tail.value(), tail.value(), -64);
    arf_add(out, out, tail.value(), ARF_PREC_EXACT, ARF_RND_DOWN);
}
// Nodes of up to two limbs: spread, clustered to their last bits, mixed with binary64 nodes, parts far apart
// in size, and mantissas of exactly 128 bits.
void long_sets() {
    Random random{20261011};
    for (size_t n : {2u, 3u, 17u, 64u, 300u}) {
        std::vector<d::ComplexBall> nodes(n);
        const auto set = [&](size_t j, double re_high, double re_low, double im_high, double im_low) {
            long_number(arb_midref(acb_realref(nodes[j].value())), re_high, re_low);
            long_number(arb_midref(acb_imagref(nodes[j].value())), im_high, im_low);
        };
        for (size_t j = 0; j < n; ++j) set(j, random.signed_unit(), random.signed_unit(), random.signed_unit(), random.signed_unit());
        check(nodes, "long, spread, " + std::to_string(n) + " nodes", false);
        // Clusters: groups that share their binary64 parts and differ 2^-64 and less below them, or in the
        // last of their 117 bits.
        for (size_t j = 0; j < n; ++j) {
            const double base = 1 + static_cast<double>(j / 4);
            set(j, base, static_cast<double>(j % 4) * (j % 8 < 4 ? 1.0 : 0x1p-52), -0.375, j % 2 ? 0.5 : 0.5 + 0x1p-52);
        }
        check(nodes, "long, clustered, " + std::to_string(n) + " nodes", false);
        // Every other node a binary64 one.
        for (size_t j = 0; j < n; ++j)
            set(j, random.signed_unit(), j % 2 ? random.signed_unit() : 0, random.signed_unit(), j % 2 ? random.signed_unit() : 0);
        check(nodes, "long, mixed with binary64, " + std::to_string(n) + " nodes", n < 2);
        // Parts far apart in size (beyond the exact difference's reach), zero parts.
        for (size_t j = 0; j < n; ++j) {
            const int scale = static_cast<int>(random.next() % 400) - 200;
            set(j, j % 5 == 0 ? 0 : std::ldexp(random.signed_unit(), scale), j % 5 == 0 ? 0 : std::ldexp(random.signed_unit(), scale),
                j % 7 == 0 ? (j % 5 == 0 ? static_cast<double>(j + 1) : 0) : std::ldexp(random.signed_unit(), -scale),
                j % 7 == 0 ? 0 : std::ldexp(random.signed_unit(), -scale));
        }
        check(nodes, "long, mixed scales, " + std::to_string(n) + " nodes", false);
        // Mantissas of exactly 128 bits: 2^127 + ... + 1 scaled, every limb pattern at its ends.
        for (size_t j = 0; j < n; ++j)
            for (arb_ptr part : {acb_realref(nodes[j].value()), acb_imagref(nodes[j].value())}) {
                d::Float low;
                arf_set_ui(arb_midref(part), (random.next() | (uint64_t{1} << 63)) | 1);
                arf_mul_2exp_si(arb_midref(part), arb_midref(part), 64);
                arf_set_ui(low.value(), random.next() | 1);
                arf_add(arb_midref(part), arb_midref(part), low.value(), ARF_PREC_EXACT, ARF_RND_DOWN);
                arf_mul_2exp_si(arb_midref(part), arb_midref(part), -128 + static_cast<slong>(random.next() % 5));
                if (random.next() & 1) arf_neg(arb_midref(part), arb_midref(part));
            }
        check(nodes, "long, 128-bit mantissas, " + std::to_string(n) + " nodes", false);
    }
}
void random_sets() {
    Random random{20261009};
    for (size_t n : {2u, 3u, 17u, 64u, 300u}) {
        std::vector<std::pair<double, double>> points(n);
        // Spread over the unit square.
        for (auto &z : points) z = {random.signed_unit(), random.signed_unit()};
        check(points, "spread, " + std::to_string(n) + " nodes");
        // Clusters: groups within a few units in the last place of one another, in one or both parts.
        for (size_t j = 0; j < n; ++j) {
            const double base_real = 1 + static_cast<double>(j / 4), base_imaginary = -0.375;
            points[j] = {base_real + static_cast<double>(j % 4) * 0x1p-52 * base_real,
                         j % 2 ? base_imaginary : std::nextafter(base_imaginary, 0.0)};
        }
        check(points, "clustered, " + std::to_string(n) + " nodes");
        // Parts whose exponents are far apart, zero parts, real and imaginary nodes.
        for (size_t j = 0; j < n; ++j) {
            const int scale = static_cast<int>(random.next() % 240) - 120;
            const double real = j % 5 == 0 ? 0 : std::ldexp(random.signed_unit(), scale);
            const double imaginary = j % 7 == 0 ? 0 : std::ldexp(random.signed_unit(), -scale);
            points[j] = {real, j % 5 == 0 && j % 7 == 0 ? static_cast<double>(j + 1) : imaginary};
        }
        check(points, "mixed scales, " + std::to_string(n) + " nodes");
        // The ends of the admitted range.
        for (size_t j = 0; j < n; ++j)
            points[j] = {std::ldexp(0.5 + random.unit() / 2, j % 2 ? 499 : -499) * (j % 3 ? 1 : -1),
                         std::ldexp(0.5 + random.unit() / 2, j % 4 < 2 ? 499 : -499)};
        check(points, "range ends, " + std::to_string(n) + " nodes");
    }
}
// The products over a subset against another set of points (an unmoved node's update): enclosure of the
// exact product and the derived radius for the subset's factors; a zero factor and an index outside the other
// set are refused.
void subsets() {
    Random random{20261010};
    for (size_t n : {5u, 64u, 300u}) {
        std::vector<std::pair<double, double>> now(n), before(n);
        for (size_t j = 0; j < n; ++j) {
            now[j] = {random.signed_unit(), random.signed_unit()};
            // The other set: the same points moved a little, a few far away, one on a point of this set.
            before[j] = {now[j].first + std::ldexp(random.signed_unit(), -30), now[j].second + std::ldexp(random.signed_unit(), j % 3 ? -30 : 200)};
        }
        auto nodes = balls(now), others = balls(before);
        // The second and third sizes with long points: in this set, and in the other one.
        if (n == 64)
            for (size_t j = 0; j < n; ++j) long_number(arb_midref(acb_realref(nodes[j].value())), now[j].first, random.signed_unit());
        if (n == 300)
            for (size_t j = 0; j < n; ++j)
                long_number(arb_midref(acb_imagref(others[j].value())), before[j].second,
                            random.signed_unit() * std::fabs(before[j].second));
        const d::ShortNodeProducts products(nodes), other_products(others);
        require(products.usable() && other_products.usable(), "short points were refused");
        d::ComplexBall held, exact, factor;
        d::Magnitude allowed, error;
        constexpr slong reference = 1024;
        for (size_t i = 0; i < n; i += 1 + n / 7) {
            std::vector<size_t> subset;
            for (size_t j = 0; j < n; ++j)
                if (j != i && random.next() % 3) subset.push_back(j);
            for (const auto *set : {&other_products, &products}) {
                const auto &points = set == &products ? nodes : others;
                require(products.product_over(held.value(), i, subset, *set), "a subset product was refused");
                acb_one(exact.value());
                for (const size_t j : subset) {
                    acb_sub(factor.value(), nodes[i].value(), points[j].value(), reference);
                    acb_mul(exact.value(), exact.value(), factor.value(), reference);
                }
                require(acb_contains(held.value(), exact.value()), "a subset product does not enclose the exact one");
                acb_get_mag(allowed.value(), exact.value());
                mag_mul_ui(allowed.value(), allowed.value(), static_cast<ulong>(subset.size() + 1));
                mag_mul_2exp_si(allowed.value(), allowed.value(), d::ShortNodeProducts::step_error_log2 + 2);
                mag_max(error.value(), arb_radref(acb_realref(held.value())), arb_radref(acb_imagref(held.value())));
                require(mag_cmp(error.value(), allowed.value()) <= 0, "a subset product is wider than the derived bound");
            }
            // The empty subset is one, exactly.
            require(products.product_over(held.value(), i, {}, other_products) && acb_is_one(held.value()),
                    "the empty product is not one");
        }
        // A zero factor: node 1 against the other set's point 2 moved onto it; an index outside the other set.
        auto moved_balls = balls(before);
        acb_set(moved_balls[2].value(), nodes[1].value());
        const d::ShortNodeProducts moved_products(moved_balls);
        const std::vector<size_t> with_zero{0, 2, 3}, outside{0, n};
        require(!products.product_over(held.value(), 1, with_zero, moved_products), "a zero factor was not refused");
        require(!products.product_over(held.value(), 1, outside, other_products), "an index outside the set was taken");
    }
}
void refusals() {
    // Coincident nodes: a zero factor.
    {
        const auto nodes = balls({{1, 2}, {0.5, -1}, {1, 2}});
        const d::ShortNodeProducts products(nodes);
        d::ComplexBall product;
        d::Magnitude separation;
        require(products.usable(), "binary64 nodes were refused");
        require(!products.product(product.value(), separation.value(), 0) &&
                    !products.product(product.value(), separation.value(), 2),
                "the product of a node that coincides with another was not refused");
        require(products.product(product.value(), separation.value(), 1), "the product of the distinct node was refused");
    }
    // An inexact node, a node of more than two limbs, a size outside the range, a single node.
    {
        auto nodes = balls({{1, 2}, {0.5, -1}});
        arb_add_error_2exp_si(acb_realref(nodes[0].value()), -80);
        require(!d::ShortNodeProducts(nodes).usable(), "an inexact node was taken");
        nodes = balls({{1, 2}, {0.5, -1}});
        d::Float tiny;
        arf_one(tiny.value());
        arf_mul_2exp_si(tiny.value(), tiny.value(), -70);
        arf_add(arb_midref(acb_imagref(nodes[1].value())), arb_midref(acb_imagref(nodes[1].value())), tiny.value(),
                ARF_PREC_EXACT, ARF_RND_DOWN);
        require(d::ShortNodeProducts(nodes).usable() && !d::ShortNodeProducts(nodes).binary64(),
                "a node of 71 bits was not taken in the long form");
        arf_mul_2exp_si(tiny.value(), tiny.value(), -60);
        arf_add(arb_midref(acb_imagref(nodes[1].value())), arb_midref(acb_imagref(nodes[1].value())), tiny.value(),
                ARF_PREC_EXACT, ARF_RND_DOWN);
        require(!d::ShortNodeProducts(nodes).usable(), "a node of 131 bits was taken");
        nodes = balls({{1, 2}, {std::ldexp(1.0, d::ShortNodeProducts::exponent_range + 2), -1}});
        require(!d::ShortNodeProducts(nodes).usable(), "a node outside the admitted range was taken");
        nodes = balls({{1, 2}});
        require(!d::ShortNodeProducts(nodes).usable(), "a single node was taken");
    }
    require(d::ShortNodeProducts::accuracy_bits(1023) == 111 && d::ShortNodeProducts::accuracy_bits(2) == 121,
            "accuracy_bits is not 121 less the bits of the factor count");
}
} // namespace
int main() {
    try {
        refusals();
        random_sets();
        long_sets();
        subsets();
    } catch (const std::exception &error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
    std::cout << "node products: ok\n";
    return 0;
}
