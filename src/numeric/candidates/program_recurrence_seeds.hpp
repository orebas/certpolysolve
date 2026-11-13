#pragma once
// Structural seeds for a composition and a squaring recurrence (docs/slp-structural-seeding-note-20261018.md,
// section 3). A composition f(g(x)) through a cut node is seeded first, by its exact preimages: f's roots s, then
// g(x) = s for each, all solved roughly. The output is F(x, q) for an earlier node q (squaring_recurrence in
// program_structure.hpp). q's roots, from a rough solve of q's cone (the same rule again below, then the binary64
// stage), each give d seeds z = r + y / q'(r), one for each root y of F(r, .): near a root r of q, q(z) is about
// q'(r) (z - r), so the output vanishes where that equals a root of F(r, .). The leftover seeds lie on the root
// bound's circle, where Aberth's sums push them toward the roots no seed holds. Proposals only: nothing here is
// certified or published, and certification stays with the separated-roots pass on the output.
// History: docs/slp-structural-seeding-note-20261018.md, section 12.
#include "core/query.hpp"
#include "core/engine_statistics.hpp"
#include "input/program/program.hpp"
#include <cstdint>

namespace certpolysolve::detail {

struct RecurrenceSeedReport {
    uint64_t levels = 0;        // the rules' levels, the output's included
    uint64_t compositions = 0;  // levels seeded by a composition's exact preimages (f(g(x)) through a cut node)
    uint64_t base_degree = 0;   // q's degree at the output's level
    uint64_t sub_sweeps = 0;    // binary64 sweeps of the rough solves below the output
    uint64_t seeded = 0, leftover = 0; // the output's seeds from q's roots, and on the root bound's circle
};
// Writes `degree` seeds for the program's output (its degree proved `degree`, its roots within `root_bound`) when the
// output is a squaring recurrence and `degree` exceeds `threshold`. The rough solves below take the circle counts
// at `precision` where the rule does not apply. False, with `out` untouched, otherwise, or when a rough solve below
// refuses or is cancelled; the statistics count every evaluation and sweep taken.
bool recurrence_seeds(acb_ptr out, const Program &program, uint64_t degree, const mag_t root_bound, slong precision,
                      const cps_options &options, certpolysolve::detail::CandidateStatistics &statistics, RecurrenceSeedReport &report,
                      uint64_t threshold = 64);

} // namespace certpolysolve::detail
