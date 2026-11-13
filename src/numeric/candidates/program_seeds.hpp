#pragma once
// Program input: starting points at the roots' radii (docs/slp-seed-placement-note-20261016.md, method A). A
// program has no coefficients, so the Newton polygon that places the coefficient route's start (initialize,
// candidates.cpp) is replaced by
// root counts on circles about the origin from the program's own evaluations. On the circle of radius rho the
// Cauchy sum (1/q) sum_g z_g p'(z_g) / p(z_g) over q equally spaced points z_g equals
// sum_j 1 / (1 - (xi_j / rho)^q) over the roots xi_j: about 1 for a root inside the circle and about 0 for one
// outside, with an error of about |xi_j / rho|^(+-q) from a root near the circle (Pan 2022, Theorem 8;
// reference-library/text/pan2022_rootsquaring.txt). A ladder of circles from the root bound inward says how many
// roots lie in each annulus, and that many seeds go on the annulus's middle circle at the placement angles, as
// the polygon's edges place the coefficient route's. The counts are proposals: a wrong one misplaces seeds by an
// annulus and costs sweeps. Nothing here certifies anything; the refinement's certificates are unchanged.
// History: docs/slp-seed-placement-note-20261016.md, section 9.
#include "core/query.hpp"
#include "core/engine_policy.hpp"
#include "core/engine_statistics.hpp"
#include "numeric/evaluate/program_machine.hpp"
#include <vector>

namespace certpolysolve::detail {
class MonicSource;

// One count on a circle: the Cauchy sum's real and imaginary parts (the imaginary part is zero but for roots
// near the circle, so its size measures the count's ambiguity), the sum's ball radius, and what stopped it.
struct CircleCount {
    double real = 0, imag = 0, radius = 0;
    bool finite = false;       // the sum was formed from finite balls
    bool near_root = false;    // p's ball contains zero at a sample while p' has correct bits: a root near the sample
    bool cancellation = false; // p and p' both without correct bits at a sample: this precision does not evaluate
                               // the program there
    bool binary64 = false;     // read from a machine evaluator (program_machine.hpp: binary64 or its range tier), not
                               // from balls
    uint64_t evaluations = 0;
};
// The g-th of `samples` points on the circle of radius 2^log2_radius about the origin: the angle
// 2 pi (g + rotation) / samples plus the offset of placement_angle, as an exact midpoint. Exposed for the tests.
void circle_sample(acb_t out, double log2_radius, unsigned samples, double rotation, unsigned g);
// The count of roots inside that circle by the Cauchy sum at its sample points. With `machine`, the samples
// go through the binary64 evaluator first (one batched call), then, when a sample's value leaves binary64's
// range and `block` is given, through the range tier (machine::BlockDouble); a sample at the floor (a root at
// a sample included), or a sum whose error estimate exceeds the tolerance, sends this circle to the balls at
// `precision`, with every remedy as before.
CircleCount circle_count(MonicSource &source, double log2_radius, unsigned samples, double rotation, slong precision,
                         certpolysolve::detail::CandidateStatistics &statistics, ProgramMachine<double> *machine = nullptr,
                         ProgramMachine<machine::BlockDouble> *block = nullptr);

struct SeedAnnulus {
    double outer_log2 = 0, inner_log2 = 0; // log2 of the radii
    uint64_t count = 0;                    // seeds placed on its middle circle
};
struct ProgramSeedReport {
    uint64_t circles = 0, evaluations = 0, sample_doublings = 0, rotations = 0, raises = 0;
    // Counts taken at the most samples while still not within the tolerance of an integer (roots within about
    // 1/samples of the circle, in log radius): accepted rounded and clamped, as proposals; among `circles`.
    uint64_t unsharp = 0;
    // Counts read in binary64 (among `circles`).
    uint64_t binary64_circles = 0;
    uint64_t placed = 0; // seeds placed from counts; the rest lie on the circle of the root bound
    slong precision = 0; // the precision the counts ended at (never lowered)
    bool budget_spent = false, cancelled = false;
    std::vector<SeedAnnulus> annuli;
};
// Places `degree` seeds for the roots of `source`, whose moduli are at most root_bound, from the counts on a
// ladder of circles: the radius halves from the bound until no root is inside, then the annuli holding the most
// roots are bisected in log radius down to a resolution, within a budget of evaluations (heuristics::program_seed_*
// and program_count_*). Where the ladder stops short (the budget, a circle whose count cannot be taken), the seeds
// not placed go on the circle of the bound, the start before this placement existed. Returns false only when
// cancelled, with `out` unspecified. `policy`, when given, receives one trace event.
bool program_seeds(acb_ptr out, MonicSource &source, slong degree, const mag_t root_bound, slong precision,
                   const cps_options &options, certpolysolve::detail::CandidateStatistics &statistics, const EnginePolicy *policy,
                   ProgramSeedReport *report = nullptr, ProgramMachine<double> *machine = nullptr);

} // namespace certpolysolve::detail
