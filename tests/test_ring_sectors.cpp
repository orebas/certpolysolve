// The admission's binary64 sector count against the 96-bit one the construction uses (ring_sectors_binary64,
// ring_sectors in piecewise_models.cpp): wherever the binary64 count is given, its bit_ceil equals the 96-bit
// count's. Checked on random rings, widths from 2^-40 to 64 octaves at starts from -200 to 200, and on rings placed
// where 4 pi gamma / R is a power of two, the only places where the count can change (there the binary64 count
// must be declined, or agree). Since 26 October 2026 the estimate takes the binary64 count (dev/docs/
// easy-inputs-profile-20261025.md, proposal A2).
#include "methods/piecewise/piecewise_roots.hpp"
#include <algorithm>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <numbers>
#include <random>
#include <stdexcept>

namespace d = certpolysolve::detail;

namespace {
size_t count(double sectors) { return std::bit_ceil(static_cast<size_t>(std::max(1.0, sectors))); }

struct Tally {
    uint64_t given = 0, declined = 0;
};

// One ring from 2^start to 2^next, as for_each_ring forms it.
void check(double start, double next, Tally &tally) {
    const auto fast = d::ring_sectors_binary64(start, next);
    if (!fast) {
        ++tally.declined;
        return;
    }
    ++tally.given;
    if (count(*fast) != count(d::ring_sectors(start, next)))
        throw std::runtime_error("the binary64 sector count differs from the 96-bit count");
}
} // namespace

int main() {
    try {
        std::mt19937_64 rng(20261026);
        std::uniform_real_distribution<double> log_width(-40.0, 6.0), place(-200.0, 200.0);
        Tally random, boundary;
        constexpr int random_rings = 200000;
        for (int i = 0; i < random_rings; ++i) {
            const double start = place(rng), next = start + std::exp2(log_width(rng));
            check(start, next, random);
        }
        // The widths at which 4 pi gamma / R = (8 pi / 3)(2^w + 1)/(2^w - 1) equals 2^k (k from 4, the first power
        // above 8 pi / 3, to 44), and just beside them.
        for (int k = 4; k <= 44; ++k) {
            const double target = std::exp2(k) * 3 / (8 * std::numbers::pi);
            const double width = std::log2(1 + 2 / (target - 1));
            for (const double offset : {0.0, 0x1p-44, -0x1p-44, 0x1p-41, -0x1p-41, 0x1p-38, -0x1p-38, 0x1p-30, -0x1p-30})
                for (const double start : {-100.25, -1.0, 0.0, 0.5, 37.0}) {
                    const double w = width * (1 + offset);
                    check(start, start + w, boundary);
                }
        }
        if (random.given < random_rings / 2) throw std::runtime_error("the binary64 count is declined too often");
        std::cout << "ring sectors: random rings " << random.given << " given, " << random.declined
                  << " declined; at powers of two " << boundary.given << " given, " << boundary.declined
                  << " declined; every given count equal to the 96-bit count\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << "test_ring_sectors: " << e.what() << '\n';
        return 1;
    }
}
