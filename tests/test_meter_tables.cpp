// The meter's tables of logarithms (src/core/dovetail.hpp): operation_weight and horner_weight compute the logarithms of
// their tables' constants once, at the first call, by the same portable::log2 on the same values, so every weight must
// be the very bits of the formula that computed them at every call. The reference functions below are those formulas
// verbatim (as of d15a108, renamed); the sweep compares the bits over precisions and point mantissas on, beside and
// between the tables' points, past their ends, and at a fixed pseudo-random sample, built with the library's
// floating-point flags.
#include "core/dovetail.hpp"
#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <vector>

namespace {
namespace d = certpolysolve::detail;
namespace portable = certpolysolve::detail::portable;

int failures = 0;
void require(bool ok, const char *what, int64_t bits, int64_t point_bits) {
    if (!ok && failures++ < 20)
        std::fprintf(stderr, "FAIL: %s at bits %lld, point bits %lld\n", what, static_cast<long long>(bits),
                     static_cast<long long>(point_bits));
}
bool same_bits(double a, double b) { return std::bit_cast<uint64_t>(a) == std::bit_cast<uint64_t>(b); }

// A fixed sequence (SplitMix64), the same on every run.
uint64_t state = 0x9e3779b97f4a7c15ULL;
uint64_t next() {
    uint64_t z = (state += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

double reference_operation_weight(int64_t bits) {
    constexpr int64_t binary64_bits = 53;
    if (bits <= binary64_bits)
        return 1;
    constexpr std::array<std::array<double, 2>, 9> measured{{{1, 148.0},
                                                             {2, 188.2},
                                                             {4, 230.6},
                                                             {8, 312.4},
                                                             {16, 684.3},
                                                             {32, 2059.3},
                                                             {64, 4976.1},
                                                             {256, 46535.6},
                                                             {1024, 227839.8}}};
    const double limbs = std::ceil(static_cast<double>(bits) / 64);
    size_t upper = 1;
    while (upper + 1 < measured.size() && measured[upper][0] < limbs)
        ++upper;
    const auto &a = measured[upper - 1], &b = measured[upper];
    const double slope = portable::log2(b[1] / a[1]) / portable::log2(b[0] / a[0]);
    return a[1] * portable::exp2(slope * portable::log2(limbs / a[0]));
}

double reference_horner_weight(int64_t bits, int64_t point_bits) {
    constexpr int64_t binary64_bits = 53;
    if (bits <= binary64_bits)
        return 1;
    constexpr size_t sizes = 9;
    constexpr std::array<double, sizes> precisions{64, 128, 256, 512, 1024, 2048, 4096, 16384, 65536};
    constexpr std::array<double, sizes> mantissas{53, 128, 256, 512, 1024, 2048, 4096, 16384, 65536};
    constexpr std::array<std::array<double, sizes>, sizes> measured{{
        {107, 107, 107, 107, 107, 107, 107, 107, 107},
        {130, 141, 141, 141, 141, 141, 141, 141, 141},
        {137, 143, 152, 152, 152, 152, 152, 152, 152},
        {142, 159, 191, 231, 231, 231, 231, 231, 231},
        {178, 223, 272, 332, 537, 537, 537, 537, 537},
        {276, 346, 479, 751, 1356, 1835, 1835, 1835, 1835},
        {403, 522, 822, 1372, 2478, 3629, 4149, 4149, 4149},
        {1165, 1597, 2741, 5505, 9529, 12601, 17280, 37441, 37441},
        {4701, 6083, 10816, 19495, 37784, 47513, 67527, 127684, 203445},
    }};
    const auto bracket = [](const std::array<double, sizes> &axis, double x, size_t &upper) {
        upper = 1;
        while (upper + 1 < sizes && axis[upper] < x)
            ++upper;
        return portable::log2(x / axis[upper - 1]) / portable::log2(axis[upper] / axis[upper - 1]);
    };
    const double p = static_cast<double>(std::max<int64_t>(bits, 64));
    const double m = static_cast<double>(std::clamp<int64_t>(point_bits, binary64_bits, std::max<int64_t>(bits, binary64_bits)));
    size_t row = 1, column = 1;
    const double along_precision = bracket(precisions, p, row);
    const double along_mantissa = std::clamp(bracket(mantissas, m, column), 0.0, 1.0);
    const auto at = [&](size_t i) {
        return portable::log2(measured[i][column - 1]) +
               along_mantissa * (portable::log2(measured[i][column]) - portable::log2(measured[i][column - 1]));
    };
    return portable::exp2(at(row - 1) + along_precision * (at(row) - at(row - 1)));
}

// Values on, beside and between the tables' points and past their ends (the points are 53 and powers of two): every
// power of two to 2^limit with its two neighbours and its 3/2 multiple, and 53 with its neighbours.
std::vector<int64_t> landmarks(int limit) {
    std::vector<int64_t> v{0, 1, 2, 52, 53, 54};
    for (int k = 2; k <= limit; ++k) {
        const int64_t x = int64_t{1} << k;
        v.insert(v.end(), {x - 1, x, x + 1, x + x / 2});
    }
    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end()), v.end());
    return v;
}
} // namespace

int main() {
    // operation_weight: every precision to 2^20 bits (every limb count to 16,384, the table's last point 1,024 and
    // beyond), the landmarks to 2^40, and a sample.
    long compared = 0;
    for (int64_t bits = -2; bits <= (int64_t{1} << 20); ++bits, ++compared)
        require(same_bits(d::operation_weight(bits), reference_operation_weight(bits)), "operation_weight", bits, 0);
    for (const int64_t bits : landmarks(40)) {
        require(same_bits(d::operation_weight(bits), reference_operation_weight(bits)), "operation_weight", bits, 0);
        ++compared;
    }
    for (int i = 0; i < 100000; ++i, ++compared) {
        const auto bits = static_cast<int64_t>(next() % (uint64_t{1} << 36));
        require(same_bits(d::operation_weight(bits), reference_operation_weight(bits)), "operation_weight", bits, 0);
    }
    std::printf("operation_weight: %ld precisions compared bit for bit\n", compared);

    // horner_weight: precisions to 600 bits one by one, then every 13th to 70,000, the landmarks to 2^24 and a
    // sample to 2^22; at each, the mantissa landmarks to 2^17, the precision's neighbours and multiples, and a sample.
    std::vector<int64_t> precisions;
    for (int64_t bits = 0; bits <= 600; ++bits)
        precisions.push_back(bits);
    for (int64_t bits = 601; bits <= 70000; bits += 13)
        precisions.push_back(bits);
    for (const int64_t bits : landmarks(24))
        precisions.push_back(bits);
    for (int i = 0; i < 2000; ++i)
        precisions.push_back(static_cast<int64_t>(next() % (uint64_t{1} << 22)));
    const auto mantissa_landmarks = landmarks(17);
    long pairs = 0;
    for (const int64_t bits : precisions) {
        std::vector<int64_t> mantissas = mantissa_landmarks;
        mantissas.insert(mantissas.end(), {bits - 1, bits, bits + 1, 2 * bits, bits + bits / 2, bits / 2});
        for (int i = 0; i < 8; ++i)
            mantissas.push_back(static_cast<int64_t>(next() % static_cast<uint64_t>(2 * bits + 2)));
        for (const int64_t point_bits : mantissas) {
            require(same_bits(d::horner_weight(bits, point_bits), reference_horner_weight(bits, point_bits)),
                    "horner_weight", bits, point_bits);
            ++pairs;
        }
    }
    std::printf("horner_weight: %ld (precision, mantissa) pairs compared bit for bit\n", pairs);
    if (failures)
        std::fprintf(stderr, "%d weights differ from the formula computed at every call\n", failures);
    return failures ? 1 : 0;
}
