#pragma once
// Piecewise's counted work: PiecewiseWork (operations by kind, each weighted by the precision it ran at), the
// factors that price the kinds (piecewise_work_factors), the binary64 operation and the estimate's overhead in that
// unit, and construction_work / counted_work. Scheduling quanta only: they set when piecewise leaves for a global
// method and what its source finisher may spend, never what is computed or certified. Part of piecewise_roots.hpp.
#include "certify/chart_statistics.hpp"
#include <cstdint>

namespace certpolysolve::detail {
// Counted work (queue item 3.1): operations by kind. Every entry is a number
// of operations times operation_weight (dovetail.hpp) of the precision they
// ran at, so entries made at different precisions add; binary64 work has
// weight one. No clock. The factors below price the kinds (counted_work).
// History: docs/code-history-20261021.md, piecewise_roots.hpp 5.
struct PiecewiseWork {
    double enclosures = 0;  // coefficient enclosures at a round's start
    double grid = 0;        // unit roots of an angular grid or a transform plan
    double weights = 0;     // chart centres and scaled batch coefficients
    double recurrence = 0;  // Taylor recurrence terms
    double butterflies = 0; // radix-2 butterflies at the working precision
    double count = 0;       // local count terms
    double restriction = 0; // coefficients of restricted models
    double machine = 0;     // binary64 proposal iteration: Horner terms
    double interactions = 0; // binary64 proposal iteration: Aberth divisions
    double polish = 0;      // Horner terms of polishing
    double polish_bounds = 0; // local model error bounds while polishing
    ChartCertificateWork certificate;
    double finisher = 0;    // the source finisher's evaluation terms
    double finisher_preparation = 0; // the finisher's source enclosure, per coefficient
    double walk = 0;        // owners examined when a finishing cohort is formed
    template <class Op> PiecewiseWork &combine(const PiecewiseWork &other, Op op) {
        op(enclosures, other.enclosures); op(grid, other.grid); op(weights, other.weights);
        op(recurrence, other.recurrence); op(butterflies, other.butterflies);
        op(count, other.count); op(restriction, other.restriction); op(machine, other.machine);
        op(interactions, other.interactions); op(polish, other.polish); op(polish_bounds, other.polish_bounds);
        op(certificate.bounds, other.certificate.bounds); op(certificate.jets, other.certificate.jets);
        op(certificate.shrinks, other.certificate.shrinks); op(finisher, other.finisher);
        op(finisher_preparation, other.finisher_preparation); op(walk, other.walk);
        return *this;
    }
    PiecewiseWork &operator+=(const PiecewiseWork &other) {
        return combine(other, [](double &a, double b) { a += b; });
    }
    PiecewiseWork &operator-=(const PiecewiseWork &other) {
        return combine(other, [](double &a, double b) { a -= b; });
    }
    double total() const {
        return enclosures + grid + weights + recurrence + butterflies + count + restriction + machine +
               interactions + polish + polish_bounds + certificate.bounds + certificate.jets + certificate.shrinks +
               finisher + finisher_preparation + walk;
    }
    bool operator==(const PiecewiseWork &) const = default;
};
// One unit for every kind. A factor is the nanoseconds one counted unit of
// its kind took on the reference machine, fitted phase by phase against the
// stage timers on the corpus. The decisions below use only ratios of counted
// work, so the unit is immaterial and another machine changes them only as
// far as its balance between the kinds differs. Five kinds were too small
// beside the others in their phase to be fitted; each is one complex
// operation at the working precision and takes the recurrence's factor. Like
// the dovetail weights these are scheduling quanta: they set when piecewise
// leaves for a global method and how much its source finisher may spend,
// never what is computed or certified; any positive factors keep both rules
// sound. History: docs/code-history-20261021.md, piecewise_roots.hpp 6.
struct PiecewiseWorkFactors {
    double enclosures = 0.434;
    double grid = 0.434;
    double weights = 0.434;
    double recurrence = 0.434;
    double butterflies = 0.972;
    double count = 1.18;
    double restriction = 0.318;
    double machine = 5.29;
    double interactions = 5.65;
    double polish = 0.704;
    double polish_bounds = 0.434;
    double bounds = 0.434;
    double jets = 0.96;
    double shrinks = 45.4;
    double finisher = 0.145;
    double finisher_preparation = 12.6;
    double walk = 156;
};
inline constexpr PiecewiseWorkFactors piecewise_work_factors{};
// A binary64 complex multiply-add in the factors' unit: the dovetail meter's operation (operation_weight,
// dovetail.hpp, 0.907 ns on the same reference machine), and so the price of one operation of the global
// recipes' minimum budget, whose sweeps begin in binary64.
inline constexpr double binary64_operation_work = 0.907;
// The calibrated counted work one operation of the selector's piecewise estimate stands for, in binary64
// operations. The estimate (auto_selection.cpp) counts model construction and local solves, its admission's work
// in operations; the ledger above also counts what the estimate leaves out (certificates,
// polishing, the finisher). Measured as the largest ratio of a finished solve's counted work to its estimate
// over the inputs Auto sends to piecewise and piecewise finishes in its own regime (the Mahler product with
// radial clusters, where the estimate fails, lies far above). PiecewiseRootConfig::regime_exit_estimate uses
// it. This is an observed calibration, not a universal envelope. Like the
// factors, it sets when Auto leaves piecewise; finite-budget completion can
// change, while the original-input certificate requirements remain fixed.
// History: docs/code-history-20261021.md, piecewise_roots.hpp 7.
inline constexpr double estimate_operation_overhead = 2.25;
// Model construction alone: what construction_units counts in operations.
inline double construction_work(const PiecewiseWork &w) {
    const auto &f = piecewise_work_factors;
    return f.grid * w.grid + f.weights * w.weights + f.recurrence * w.recurrence +
           f.butterflies * w.butterflies;
}
inline double counted_work(const PiecewiseWork &w) {
    const auto &f = piecewise_work_factors;
    return construction_work(w) + f.enclosures * w.enclosures + f.count * w.count +
           f.restriction * w.restriction + f.machine * w.machine + f.interactions * w.interactions +
           f.polish * w.polish + f.polish_bounds * w.polish_bounds + f.bounds * w.certificate.bounds + f.jets * w.certificate.jets +
           f.shrinks * w.certificate.shrinks + f.finisher * w.finisher + f.finisher_preparation * w.finisher_preparation +
           f.walk * w.walk;
}
} // namespace certpolysolve::detail
