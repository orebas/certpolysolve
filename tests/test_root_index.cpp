#include <certpolysolve/certpolysolve.hpp>
#include "methods/piecewise/root_index.hpp"
#include <iostream>
#include <random>

namespace d = certpolysolve::detail;
int main() {
    certpolysolve::ThreadScope lifetime;
    std::mt19937_64 random(0x716273);
    for (slong resolution : {0, 4, 24, 48, FLINT_BITS - 3}) {
        d::RootIndex index(resolution);
        std::vector<d::Root> roots;
        for (unsigned i = 0; i < 300; ++i) {
            d::Root root;
            acb_set_si_si(root.box.value(), static_cast<slong>(random() % 17) - 8,
                          static_cast<slong>(random() % 17) - 8);
            const auto width = -static_cast<slong>(random() % 70);
            arb_add_error_2exp_si(acb_realref(root.box.value()), width);
            arb_add_error_2exp_si(acb_imagref(root.box.value()), width + 1);
            const auto scale = static_cast<slong>(random() % 17) - 8;
            acb_mul_2exp_si(root.box.value(), root.box.value(), scale);
            index.insert(root, roots.size());
            roots.push_back(std::move(root));
        }
        // Include coincident points, axis-crossing wide boxes and widths far
        // below machine relative precision. Compare exact backend geometry,
        // independently of the grid's bins, with and without an ID filter.
        for (const auto &base : roots) {
            for (const slong shift : {0, -1, 1}) {
                d::Root query = base;
                acb_mul_2exp_si(query.box.value(), query.box.value(), shift);
                for (unsigned parity = 0; parity < 3; ++parity) {
                    const auto accept = [&](size_t i) { return parity == 2 || i % 2 == parity; };
                    bool expected = false;
                    for (size_t i = 0; i < roots.size(); ++i)
                        expected |= accept(i) && acb_overlaps(query.box.value(), roots[i].box.value());
                    if (index.overlaps(query, roots, accept) != expected) {
                        std::cerr << "spatial overlap index differs from exhaustive scan\n";
                        return 1;
                    }
                }
            }
        }
    }
    // Compact exponents near the backend index boundary must use exhaustive
    // geometry without overflowing the neighbor-bin loop or expanding them.
    d::RootIndex extreme;
    std::vector<d::Root> huge;
    for (slong exponent : {WORD_MIN + 1, WORD_MIN + 2, WORD_MAX - 2, WORD_MAX - 1}) {
        d::Root root;
        acb_one(root.box.value());
        acb_mul_2exp_si(root.box.value(), root.box.value(), exponent);
        extreme.insert(root, huge.size());
        huge.push_back(std::move(root));
    }
    for (size_t i = 0; i < huge.size(); ++i)
        if (!extreme.overlaps(huge[i], huge, [i](size_t j) { return i == j; })) {
            std::cerr << "extreme exponent fallback lost an overlap\n";
            return 1;
        }
    std::cout << "13500 filtered overlap queries and compact exponent boundaries match exhaustive geometry\n";
}
