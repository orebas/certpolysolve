#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <optional>
#include <tuple>
#include <utility>
#include <vector>
#include <flint/acb_poly.h>
#include "core/portable_math.hpp"

namespace certpolysolve::detail {
// Dovetailing a factor's fast iteration (secular) with certified subdivision.
// The abstract motivation: if the latter is complete, both engines receive
// fair bounded actions, and accounting bounds all work, a fixed share s gives
// a finish within about min(T_iteration / (1 - s), T_subdivision / s) plus
// one action. Those hypotheses need a separate implementation argument;
// this scheduler alone does not establish termination or a bit-cost bound.
// See dev/docs/m2a-fallback-design-20261105.md for the outstanding engine,
// pre-engagement, transfer and accounting obligations.

// Work is counted in binary64 complex multiply-adds, deterministically (no
// clock). An operation at `bits` of Arb precision costs operation_weight(bits)
// of them. The weight is a scheduling quantum: it sets how the two methods
// share work, never the validity of a certificate. These calibrated weights
// are not proved upper bounds for all bit operations or unmetered work.
// The table comes from a calibration on the reference machine.
// History: docs/code-history-20261021.md, dovetail.hpp 1.
//
// The meter's interpolations and counts below are in base 2 by portable::log2
// and portable::exp2 (src/portable_math.hpp), never the C library's log, exp,
// pow or log2: the counted work feeds every work budget, the regime exit, the
// certificate attempt's budget and the dovetail shares, so it must be the same
// number on every platform for those decisions to be; it need not be exact.
// History: docs/code-history-20261021.md, dovetail.hpp 2.
//
// The meter runs at every counted operation, and portable::log2 is a series,
// so the logarithms that depend only on a table's constants are computed once,
// at the first call, by the same portable::log2 on the same values: the same
// bits as computing them at every call (tests/test_meter_tables.cpp compares
// the two over a sweep).
// History: docs/code-history-20261021.md, dovetail.hpp 3.
inline double operation_weight(int64_t bits) {
    constexpr int64_t binary64_bits = 53;
    if (bits <= binary64_bits)
        return 1;
    // acb_addmul at 64 * limbs bits, relative to a binary64 complex
    // multiply-add, measured on the reference machine. Fixed overhead
    // dominates at one limb, multiplication beyond; interpolate in log-log
    // and extend the last segment's slope past the table.
    // History: docs/code-history-20261021.md, dovetail.hpp 4.
    static constexpr std::array<std::array<double, 2>, 9> measured{{{1, 148.0},
                                                                    {2, 188.2},
                                                                    {4, 230.6},
                                                                    {8, 312.4},
                                                                    {16, 684.3},
                                                                    {32, 2059.3},
                                                                    {64, 4976.1},
                                                                    {256, 46535.6},
                                                                    {1024, 227839.8}}};
    // slopes[upper]: the log-log slope of the segment that ends at measured[upper].
    static const std::array<double, measured.size()> slopes = [] {
        std::array<double, measured.size()> slope{};
        for (size_t upper = 1; upper < measured.size(); ++upper) {
            const auto &a = measured[upper - 1], &b = measured[upper];
            slope[upper] = portable::log2(b[1] / a[1]) / portable::log2(b[0] / a[0]);
        }
        return slope;
    }();
    const double limbs = std::ceil(static_cast<double>(bits) / 64);
    size_t upper = 1;
    while (upper + 1 < measured.size() && measured[upper][0] < limbs)
        ++upper;
    const auto &a = measured[upper - 1];
    return a[1] * portable::exp2(slopes[upper] * portable::log2(limbs / a[0]));
}

// Subdivision's steps are Taylor shifts, products and Graeffe steps of its
// factor, which FLINT computes by its block multiplication
// (arb_poly/mullow_block.c): after removing the average exponent slope, the
// coefficients are cut into blocks whose exponent height is at most
// 3 prec + 512 bits, and every pair of blocks is multiplied as integer
// polynomials by a Kronecker/FFT product. The cost therefore follows the
// polynomial's exponent profile, not only its length: a flat profile is one
// block and one quasi-linear product, far below degree^2, while each Graeffe
// step doubles the profile's curvature, which adds blocks. The meter counts
// the blocks as FLINT would and charges every block pair its FFT's
// butterflies, 3 (l_i + l_j) log2(2 (l_i + l_j)), or the schoolbook l_i l_j
// when that is less, weighted at the precision like every other operation.
// History: docs/code-history-20261021.md, dovetail.hpp 5.
inline double block_pair_operations(double a, double b) {
    const double fft = 3 * (a + b) * portable::log2(2 * (a + b));
    return std::min(a * b, fft);
}
// FLINT's block boundaries for one operand: the lengths of the blocks.
inline std::vector<double> flint_block_lengths(const acb_poly_t f, slong prec) {
    std::vector<double> lengths;
    std::vector<std::pair<slong, double>> exponents; // (index, top exponent)
    for (slong j = 0; j < acb_poly_length(f); ++j) {
        double top = -std::numeric_limits<double>::infinity();
        for (const arf_struct *m : {arb_midref(acb_realref(f->coeffs + j)), arb_midref(acb_imagref(f->coeffs + j))})
            if (!arf_is_zero(m) && arf_is_finite(m) && fmpz_fits_si(ARF_EXPREF(m)))
                top = std::max(top, static_cast<double>(fmpz_get_si(ARF_EXPREF(m))));
        if (std::isfinite(top))
            exponents.emplace_back(j, top);
    }
    if (exponents.size() < 2) {
        lengths.push_back(static_cast<double>(std::max<slong>(1, acb_poly_length(f))));
        return lengths;
    }
    // The average slope between the first and last nonzero coefficients
    // (_arb_poly_get_scale), then blocks of height below 3 prec + 512.
    const double scale = std::round((exponents.back().second - exponents.front().second) /
                                    static_cast<double>(exponents.back().first - exponents.front().first));
    const double height = 3.0 * static_cast<double>(prec) + 512;
    double top = 0, bottom = 0;
    slong start = exponents.front().first;
    bool open = false;
    for (const auto &[index, exponent] : exponents) {
        const double v = exponent - scale * static_cast<double>(index);
        if (!open) {
            top = bottom = v;
            start = index;
            open = true;
            continue;
        }
        const double t = std::max(top, v), b = std::min(bottom, v);
        if (t - b < height) {
            top = t;
            bottom = b;
        } else {
            lengths.push_back(static_cast<double>(index - start));
            top = bottom = v;
            start = index;
        }
    }
    lengths.push_back(static_cast<double>(acb_poly_length(f) - start));
    return lengths;
}
// A product of f by a polynomial of its own block structure (a square, or a
// Graeffe step's two half-length squares, whose pair sum is the same to the
// quasi-linear count).
inline double block_product_operations(const acb_poly_t f, slong prec) {
    const auto lengths = flint_block_lengths(f, prec);
    double total = 0;
    for (const double a : lengths)
        for (const double b : lengths)
            total += block_pair_operations(a, b);
    return total;
}
// A Taylor shift (acb_poly/taylor_shift.c): Horner's n^2 up to length 30,
// one convolution (a product and n scalings) when the precision exceeds
// twice the length, else the divide-and-conquer recursion, log2(n)/2
// products' worth for a quasi-linear product.
inline double block_taylor_shift_operations(const acb_poly_t f, slong prec) {
    const double n = static_cast<double>(std::max<slong>(1, acb_poly_length(f)));
    if (n <= 30)
        return n * n;
    const double product = block_product_operations(f, prec);
    if (static_cast<double>(prec) > 2 * n)
        return product + n;
    return std::max(1.0, portable::log2(n) / 2) * product;
}

// A Horner step y <- y x + c at `bits` of precision when the point x is exact
// with a `point_bits`-bit mantissa, relative to a binary64 complex
// multiply-add. Secular regenerates each node's residue at that node's own
// precision, and its nodes come from lower-precision iterations, so their
// mantissas are short, and a step at a short point costs a small part of
// operation_weight(bits): charging the full weight would overstate secular's
// work and give Subdivision more than its share of the run time. Measured on
// the reference machine; a scheduling quantum like operation_weight.
// History: docs/code-history-20261021.md, dovetail.hpp 6.
inline double horner_weight(int64_t bits, int64_t point_bits) {
    constexpr int64_t binary64_bits = 53;
    if (bits <= binary64_bits)
        return 1;
    constexpr size_t sizes = 9;
    static constexpr std::array<double, sizes> precisions{64, 128, 256, 512, 1024, 2048, 4096, 16384, 65536};
    static constexpr std::array<double, sizes> mantissas{53, 128, 256, 512, 1024, 2048, 4096, 16384, 65536};
    // measured[i][j]: precision precisions[i], mantissa mantissas[j]; a
    // mantissa longer than the precision costs what the full one does.
    static constexpr std::array<std::array<double, sizes>, sizes> measured{{
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
    // The table's logarithms, computed once (operation_weight says why): of
    // every measured value, and of every step along each axis
    // (steps[k] = log2(axis[k] / axis[k - 1]), k >= 1).
    struct Logarithms {
        std::array<std::array<double, sizes>, sizes> measured{};
        std::array<double, sizes> precision_steps{}, mantissa_steps{};
    };
    static const Logarithms logarithms = [] {
        Logarithms l;
        for (size_t i = 0; i < sizes; ++i)
            for (size_t j = 0; j < sizes; ++j)
                l.measured[i][j] = portable::log2(measured[i][j]);
        for (size_t k = 1; k < sizes; ++k) {
            l.precision_steps[k] = portable::log2(precisions[k] / precisions[k - 1]);
            l.mantissa_steps[k] = portable::log2(mantissas[k] / mantissas[k - 1]);
        }
        return l;
    }();
    // Log-log interpolation along one axis; past the table, the last
    // segment's slope continues.
    const auto bracket = [](const std::array<double, sizes> &axis, const std::array<double, sizes> &steps, double x,
                            size_t &upper) {
        upper = 1;
        while (upper + 1 < sizes && axis[upper] < x)
            ++upper;
        return portable::log2(x / axis[upper - 1]) / steps[upper];
    };
    const double p = static_cast<double>(std::max<int64_t>(bits, 64));
    const double m = static_cast<double>(std::clamp<int64_t>(point_bits, binary64_bits, std::max<int64_t>(bits, binary64_bits)));
    size_t row = 1, column = 1;
    const double along_precision = bracket(precisions, logarithms.precision_steps, p, row);
    const double along_mantissa = std::clamp(bracket(mantissas, logarithms.mantissa_steps, m, column), 0.0, 1.0);
    const auto &logs = logarithms.measured;
    const auto at = [&](size_t i) {
        return logs[i][column - 1] + along_mantissa * (logs[i][column] - logs[i][column - 1]);
    };
    return portable::exp2(at(row - 1) + along_precision * (at(row) - at(row - 1)));
}

// The share schedule of one factor.
class DovetailSchedule {
  public:
    // (unresolved roots, largest unresolved count, -regions): prefer_cover's
    // structure tuple, smaller is better.
    using Structure = std::tuple<uint64_t, uint64_t, int64_t>;
    // Subdivision's share of the metered work is 2^-share_log2 (at least 1:
    // a share above one half would make the iteration the junior partner).
    // prior_iteration: the iteration's work before engagement.
    // first_step: the expected cost of Subdivision's first step.
    explicit DovetailSchedule(unsigned share_log2 = 1, double prior_iteration = 0, double first_step = 1)
        : share_(std::ldexp(1.0, -static_cast<int>(std::max(1u, share_log2)))), base_share_(share_),
          prior_(prior_iteration), next_(std::max(1.0, first_step) / step_growth) {}
    double share() const { return share_; }
    // At the precision budget the iteration's model is final: regeneration
    // cannot raise its source precision, so whether the budget suffices is
    // for Subdivision to decide, and its share stops being the junior one.
    // The two engines then alternate as equal partners (a share of one half:
    // Subdivision runs while its work, the next step included, is within the
    // iteration's). The iteration keeps every chance of finishing at the
    // budget; Subdivision's verdict is no longer delayed by a share meant for
    // the regime in which the iteration could still improve its model.
    // A raised budget puts the factor below it again, and the share back.
    void at_budget(bool at) { share_ = at ? std::max(base_share_, 0.5) : base_share_; }
    bool at_budget_share() const { return share_ > base_share_; }
    // Subdivision runs a step when its share of all the work, that step
    // included, covers it. The step's cost is not known in advance. It is
    // predicted as step_growth times the last step's cost (first_step for
    // the first): Subdivision's precision doubles when a test is undecided,
    // and a step at doubled precision costs up to four times the last
    // (schoolbook multiplication). Predicting high only delays a step, and
    // the allowance it did not use stays banked; predicting low lets
    // Subdivision run ahead of its share, which is then paid off before the
    // next step. A short solve that finishes before the share covers the
    // first step pays nothing.
    static constexpr double step_growth = 4;
    bool subdivision_due() const {
        const double predicted = step_growth * next_;
        return subdivision_ + predicted <= share_ * (iteration_ + subdivision_ + predicted);
    }
    double iteration() const { return iteration_; }
    double subdivision() const { return subdivision_; }
    // The largest step the share covers now: the cost x with
    // subdivision + x <= share (iteration + subdivision + x). Subdivision's
    // disk tests run under this allowance (Subdivision::step): each Taylor
    // shift and Graeffe step is priced before it runs, so no step overruns
    // the share, whatever the last step cost.
    double available() const {
        return std::max(0.0, (share_ * (iteration_ + subdivision_) - subdivision_) / (1 - share_));
    }
    // A step stopped before an operation of this cost: the next step is due
    // exactly when the share covers it (subdivision_due predicts step_growth
    // times next_).
    void wait_for(double cost) { next_ = std::max(1.0, cost) / step_growth; }
    // Every step is charged at least one operation, so charged work bounds the
    // number of steps and a run of steps that compute nothing still advances
    // the schedule (the guarantee needs every step to count).
    void charge_subdivision(double work) {
        next_ = std::max(1.0, work);
        subdivision_ += next_;
    }
    // Explicit forecast for the integrated research candidate and its replay
    // driver. Pending primitive quotes still override this via wait_for().
    void charge_subdivision(double work, double forecast_growth) {
        charge_subdivision(work);
        if (forecast_growth != step_growth)
            wait_for(forecast_growth * std::max(1.0, work));
    }
    void charge_iteration(double work) {
        iteration_ += std::max(1.0, work);
        ++steps_;
    }
    // The mean cost of an iteration step since engagement (before the first
    // one, all the iteration's earlier work stands in for a step). A step of
    // either method costs the caller one work unit, so Subdivision's turn is
    // sized to this: work units then follow the share as metered work does.
    double iteration_step() const { return steps_ ? iteration_ / static_cast<double>(steps_) : std::max(1.0, prior_); }
    // A record of the iteration is a strictly better cover structure, a
    // halved largest unresolved radius, or a halved largest correction of
    // its last packet (the iteration's own convergence, which certificates
    // report only at the end: without it an iteration whose nodes are still
    // converging can count as overdue). All three are bounded for a bounded
    // precision budget.
    // History: docs/code-history-20261021.md, dovetail.hpp 7.
    void observe(const Structure &structure, double radius_log2,
                 double correction_log2 = std::numeric_limits<double>::infinity()) {
        const bool better = !structure_ || structure < *structure_;
        // Halved relative to the value at the last such record, so a slow
        // decrease still records once per halving, never on every step.
        const bool halved = std::isfinite(radius_log2) && radius_log2 <= radius_log2_ - 1;
        const bool contracted = std::isfinite(correction_log2) && correction_log2 <= correction_log2_ - 1;
        if (better) structure_ = structure;
        if (halved) radius_log2_ = radius_log2;
        if (contracted) correction_log2_ = correction_log2;
        if (better || halved || contracted)
            record_work_ = iteration_;
    }
    // True while the iteration has gone longer without a record than all of
    // its work before that record, engagement counting as the first record.
    // The schedule itself never uses this: it only matters once Subdivision
    // has spent a budget and the share has nothing left to give (solver.cpp).
    bool overdue() const { return iteration_ - record_work_ > prior_ + record_work_; }

  private:
    double share_, base_share_, prior_, next_;
    double iteration_ = 0, subdivision_ = 0, record_work_ = 0;
    uint64_t steps_ = 0;
    std::optional<Structure> structure_;
    double radius_log2_ = std::numeric_limits<double>::infinity();
    double correction_log2_ = std::numeric_limits<double>::infinity();
};
} // namespace certpolysolve::detail
