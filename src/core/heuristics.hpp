#pragma once
#include "core/portable_math.hpp"
#include <limits>
#include <flint/flint.h>
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>

// Private build-time ablations. Supply these through CMAKE_CXX_FLAGS so every
// translation unit sees the same values; never override only one target.
// These are cost hypotheses, not mathematical eligibility or global size caps.
#ifndef CERTPOLYSOLVE_EXPERIMENT_LOCAL_COMPONENT_DEGREE
#define CERTPOLYSOLVE_EXPERIMENT_LOCAL_COMPONENT_DEGREE 64
#endif
#ifndef CERTPOLYSOLVE_EXPERIMENT_SECULAR_RECONSTRUCTION_DEGREE
#define CERTPOLYSOLVE_EXPERIMENT_SECULAR_RECONSTRUCTION_DEGREE 2048
#endif
#ifndef CERTPOLYSOLVE_EXPERIMENT_SECULAR_BLOCK_SEED_DEGREE
#define CERTPOLYSOLVE_EXPERIMENT_SECULAR_BLOCK_SEED_DEGREE 2048
#endif
#if defined(CERTPOLYSOLVE_BUILD_RESEARCH) && !defined(CERTPOLYSOLVE_EXPERIMENT_BINARY64_SEED_DEGREE)
#define CERTPOLYSOLVE_EXPERIMENT_BINARY64_SEED_DEGREE 2048
#endif
#ifndef CERTPOLYSOLVE_EXPERIMENT_CLUSTER_RECOVERY_DEGREE
#define CERTPOLYSOLVE_EXPERIMENT_CLUSTER_RECOVERY_DEGREE 2048
#endif
// The program route's binary64 stage (program_machine.hpp): 0 runs the preliminary phase in balls only, for the
// ablation measurement. History: docs/code-history-20261021.md, heuristics.hpp 1.
#ifndef CERTPOLYSOLVE_EXPERIMENT_PROGRAM_MACHINE
#define CERTPOLYSOLVE_EXPERIMENT_PROGRAM_MACHINE 1
#endif
// The program route's separated-roots pass (program_certificate.hpp): 0 sends every solve to the refinement, for
// the ablation measurement. History: docs/code-history-20261021.md, heuristics.hpp 2.
#ifndef CERTPOLYSOLVE_EXPERIMENT_PROGRAM_DISKS
#define CERTPOLYSOLVE_EXPERIMENT_PROGRAM_DISKS 1
#endif
// The program route's power lift (program_structure.hpp): 0 solves a power output c q^k itself, for the ablation
// measurement. History: docs/code-history-20261021.md, heuristics.hpp 3.
#ifndef CERTPOLYSOLVE_EXPERIMENT_PROGRAM_LIFT
#define CERTPOLYSOLVE_EXPERIMENT_PROGRAM_LIFT 1
#endif
// The program route's product lift (program_structure.hpp, output_factors): 0 solves a product of distinct factors
// whole, for the ablation measurement. History: docs/code-history-20261021.md, heuristics.hpp 4.
#ifndef CERTPOLYSOLVE_EXPERIMENT_PROGRAM_PRODUCT_LIFT
#define CERTPOLYSOLVE_EXPERIMENT_PROGRAM_PRODUCT_LIFT 1
#endif
// The program route's recurrence seeds (program_recurrence_seeds.hpp): 0 places every program's seeds by the
// circle counts, for the ablation measurement. History: docs/code-history-20261021.md, heuristics.hpp 5.
#ifndef CERTPOLYSOLVE_EXPERIMENT_PROGRAM_RECURRENCE_SEEDS
#define CERTPOLYSOLVE_EXPERIMENT_PROGRAM_RECURRENCE_SEEDS 1
#endif

namespace certpolysolve::detail::heuristics {
// Algorithm/cache choices, not input admission limits. Every disabled optional
// path has the ordinary candidate, evaluation or subdivision fallback.
//
// Classes (principled audit, fix 15; docs/principled-audit-20260929.md):
// - Cost crossover: a degree or size above which an optional helper is
//   skipped because it is assumed to cost more than it saves. Some are
//   measured crossovers (each says what was measured against what), the
//   others round-number cost hypotheses; none names or detects an input.
//   The ordinary path remains, so a crossover changes speed and trajectories,
//   never the validity of a certificate. Finite-budget goal attainment and
//   the stop reason can change with the trajectory.
// - Scheduling quantum: work done between checkpoints. Certificate conditions
//   are independent of this price; eventual progress also requires fair
//   service and the algorithm's own progress conditions.
// - Proposal resolution: how finely a heuristic proposal (a centre) is
//   rounded; exact arithmetic follows, so it changes cost, not validity.

// Cost crossover: Aberth interaction sums use the interaction tree from this
// degree; below it the direct O(n^2) sum is cheaper.
inline constexpr slong interaction_tree_degree = 128;
// W9's split of the multiprecision secular packet's Jacobi sweep (candidates.cpp, parallel.hpp): from this many
// points, in blocks of secular_packet_block points. Both fix the blocks by the degree alone; a point costs one
// O(n) model sum and its interactions, so a block is many times a join's cost.
inline constexpr size_t secular_packet_split_points = 256;
inline constexpr size_t secular_packet_block = 32;
// W9's split of the binary64 engine's evaluations (machine_candidates.cpp, parallel.hpp): from this many points,
// in blocks of machine_split_block points, fixed by the degree alone. A point costs a Horner pass over its
// factor's terms, so below a few hundred points a sweep's evaluations are too short to share.
inline constexpr size_t machine_split_points = 512;
inline constexpr size_t machine_split_block = 64;
// Cost crossover: the certificate attempt's binary64 Aberth packet sums its
// root interactions with the same tree from this degree (machine_candidates.cpp):
// the measured crossover of the tree's cost a point, construction included,
// against the packet's direct sum, on points spread like a polynomial's roots
// near the unit circle, rounded up to a power of two.
// History: docs/code-history-20261021.md, heuristics.hpp 6.
inline constexpr size_t machine_interaction_tree_degree = 2048;
// Cost crossover for small factors: component_seeds never escalates a local
// model at or below it, and small exact factors whose ordinary candidates
// stall move to the structured secular path.
inline constexpr slong small_complex_degree = 64;
// Cost crossover for component-scale restarts: a degree-k local Taylor model
// for components of at most this many members. The source's degree does not
// bound it: a close pair refused the model would be put back on a circle at
// every precision and then walk in at about a bit a sweep, and on the other
// large inputs measured, such a bound changed no work.
// History: docs/code-history-20261021.md, heuristics.hpp 7.
inline constexpr slong local_component_degree = CERTPOLYSOLVE_EXPERIMENT_LOCAL_COMPONENT_DEGREE;
static_assert(local_component_degree >= 0 && local_component_degree < WORD_MAX - 1);
// Cost crossover: above it, sparse factors use the automatic secular
// reconstruction.
inline constexpr slong secular_reconstruction_degree = CERTPOLYSOLVE_EXPERIMENT_SECULAR_RECONSTRUCTION_DEGREE;
static_assert(secular_reconstruction_degree >= 0 && secular_reconstruction_degree < WORD_MAX - 1);
// Cost crossovers (work guards) for optional initializers; above each, the
// ordinary secular initialization, ordinary candidates, or the retained count
// cover with subdivision remain.
inline constexpr slong secular_block_seed_degree = CERTPOLYSOLVE_EXPERIMENT_SECULAR_BLOCK_SEED_DEGREE;
static_assert(secular_block_seed_degree >= 0 && secular_block_seed_degree < WORD_MAX - 1);
#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
// The research build's DoubleSeed test engine's start (binary64_seeds).
inline constexpr slong binary64_seed_degree = CERTPOLYSOLVE_EXPERIMENT_BINARY64_SEED_DEGREE;
static_assert(binary64_seed_degree >= 0 && binary64_seed_degree < WORD_MAX - 1);
#endif
inline constexpr slong cluster_recovery_degree = CERTPOLYSOLVE_EXPERIMENT_CLUSTER_RECOVERY_DEGREE;
static_assert(cluster_recovery_degree >= 0 && cluster_recovery_degree < WORD_MAX - 1);
// Cost crossover: native secular input of at most native_secular_fallback_degree leaves the native iteration for the
// expanded coefficients after native_secular_fallback_rounds rounds at native_secular_fallback_bits bits or more
// (solver.cpp, native_advance; the audit's fix 2 would remove it once a stalled native iteration reaches the precision
// budget on its own).
inline constexpr slong native_secular_fallback_degree = 256;
inline constexpr unsigned native_secular_fallback_rounds = 16;
inline constexpr uint64_t native_secular_fallback_bits = 512;
// Cost crossover: exact centering (an exact Taylor shift, O(n^2) exact
// operations) is attempted from this degree.
inline constexpr slong centering_min_degree = 64;
// The all-real seed pass (real_line_seeds.hpp) runs on factors of at least this degree: the measured crossover, with
// the cell check below, of real-only solves with the pass against those without it, on the inputs of the real lists
// and the corpora that Newton's inequalities do not refuse. Below it the crawl from the circles is short, and
// seeding costs more than it saves. History: docs/code-history-20261021.md, heuristics.hpp 8.
inline constexpr slong real_line_seed_min_degree = 64;
// The pass ends without seeds when the nodes' certified sign changes number fewer than this share of the degree.
// On the same inputs the shares fall into two groups with a gap between, where the value lies. Above it are the
// larger inputs whose seeds paid, the smallest shares among them those of random Jacobi matrices, whose roots are
// spaced like independent points, so that some cells hold two. Below it are the inputs with mostly non-real roots
// or with roots too dense or too widely scaled for the nodes (clustered pairs, geometric spreads), on which seeding
// costs far more than it saves; the few small gains there are given up. 1/2 sits below the independent-spacing case
// with a margin. History: docs/code-history-20261021.md, heuristics.hpp 9.
inline constexpr double real_line_seed_min_cell_fraction = 0.5;
// Proposal resolution: an unresolved cluster's centre is added to the old
// centre at min(this, working precision / 2) bits.
inline constexpr slong global_centering_bits = 1024;
// Proposal resolution: exact centres are rounded to this many significant
// bits, a low-denominator dyadic that keeps the exact translation cheap
// (coefficient growth is about degree x centre bits).
inline constexpr slong initial_center_bits = 8;
// Proposal resolution: the same for local restart centres.
inline constexpr slong local_center_bits = 32;
// Scheduling quantum and cost unit: the classical candidate stage runs at
// least max(this, 2n) sweeps at a precision before escalating; the selector
// and piecewise's regime exit price a global solve as this many sweeps of
// degree x (degree + support) coefficient operations.
inline constexpr uint64_t initial_sweep_allowance = 256;
// Scheduling quanta: the sweeps one packet makes before its owner looks again
// (validation, model stopping, cancellation, the caller's budgets). Any value
// is correct. No option changes them: the caller's iteration budget
// (cps_options::max_iterations) only cuts a call's last packet short
// (cps_options::packet_sweeps).
// One classical candidate packet.
inline constexpr unsigned classical_packet_sweeps = 8;
// One packet on the secular model, and the polish packet of secular's
// preliminary phase.
inline constexpr unsigned secular_packet_sweeps = 8;
// Secular's ordinary preliminary packet: MPSolve's float preliminary phase
// permits 21 sweeps (its do/while test is inclusive).
inline constexpr unsigned preliminary_packet_sweeps = 21;
// A quantum of the binary64 packet of the Newton-certificate attempt
// (machine_newton_cover), which goes on in further quanta while they settle
// points and the attempt's budget pays for them.
inline constexpr unsigned machine_cover_sweeps = 32;
// The certificate attempt's start on a sparse polynomial (initial_candidates with binomial seeds): each
// Newton-polygon edge's seeds at the roots of its binomial c_a + c_b x^m, turned by this fraction of their spacing
// 2 pi / m. A fraction strictly between 0 and 1/2 keeps every seed nearer its own binomial root than either
// neighbour; 0 and 1/2 leave a real edge's seeds closed under conjugation, which a Jacobi sweep of a real source
// preserves together with the number of real-axis points (placement_angle). A tenth lies well inside that range.
// From the equally spaced start a sparse polynomial's seeds can lie far off the roots of one edge, near another
// ring, and the attempt be refused where it certifies from the binomial start.
// History: docs/code-history-20261021.md, heuristics.hpp 10.
inline constexpr double binomial_seed_spacing_fraction = 0.1;
// The program route's start (program_seeds.cpp; docs/slp-seed-placement-note-20261016.md): root counts on
// circles by the Cauchy sum at this many sample points, doubled while a count is not within
// program_count_tolerance of an integer with an imaginary part as small, up to program_count_samples_max.
// A root at relative distance t from the circle perturbs the count by about (1 - t)^q: 64 samples resolve
// the roots more than 7% from the circle to 0.01, 512 samples those more than 1% from it.
inline constexpr unsigned program_count_samples = 64, program_count_samples_max = 512;
inline constexpr double program_count_tolerance = 0.25;
// The ladder: the radius halves from the root bound until no root is inside (at most program_ladder_max_rungs
// rungs: a root at the origin never empties a circle), then the annuli holding more than
// max(program_annulus_capacity_min, degree / program_annulus_capacity_fraction) roots are bisected in log
// radius down to program_ladder_resolution_log2, the fullest first.
inline constexpr double program_ladder_step_log2 = 1, program_ladder_resolution_log2 = 0.25;
inline constexpr unsigned program_ladder_max_rungs = 256;
inline constexpr uint64_t program_annulus_capacity_min = 4, program_annulus_capacity_fraction = 32;
// The placement's evaluations: this many sweeps' worth (a sweep evaluates the program at degree points), and at
// least the floor so that a small degree gets its whole ladder. From the bound's circle the preliminary phase
// crawls about (n / 2) ln(R / r) sweeps. History: docs/code-history-20261021.md, heuristics.hpp 11.
inline constexpr uint64_t program_seed_budget_sweeps = 16, program_seed_budget_floor = 4096;
// The program route's binary64 stage (program_machine.hpp, program_advance): the preliminary phase's sweeps in
// doubles before any ball packet, batched node by node. program_machine_stage: the ablation (above).
inline constexpr bool program_machine_stage = CERTPOLYSOLVE_EXPERIMENT_PROGRAM_MACHINE != 0;
// Sweeps a work unit: a double sweep's measured cost is a small fraction of a ball sweep's, so a quantum of 32
// against the ball packet's 8 keeps a work unit's cost of the same order.
// History: docs/code-history-20261021.md, heuristics.hpp 12.
inline constexpr unsigned program_machine_packet_sweeps = 32;
// A value or derivative within this factor of its noise estimate is unreadable (the floor): a judgment, the
// noise being a lower bound of the rounding; MachineEngine's Roundoff rule uses 32 against a majorant.
inline constexpr double program_machine_floor_factor = 16;
// A relative step below 2^this is a converged point in binary64 (MachineEngine's SmallStep rule: 8 - 53).
inline constexpr int program_machine_small_step_log2 = -45;
// A point whose relative step has not fallen by a bit against its smallest so far for this many sweeps, or
// for the window (n + 2) / 2 when that is more, is at the floor (a straggler whose noise the floor estimate
// undercounts); at a large degree the window keeps the bulk from being floored while Aberth's iteration sorts
// out the angles. History: docs/code-history-20261021.md, heuristics.hpp 13.
inline constexpr uint32_t program_machine_patience = 32;
// The ball packet that polishes points every one of which converged in binary64: from 2^-45 a sweep or two
// reach the packet's convergence, and the refinement then regenerates once, not twice.
inline constexpr unsigned program_machine_polish_sweeps = 2;
// The binary64 stage builds its interaction tree over all n points only while this many or more are active: the
// tree's measured build cost a point is more than this many terms of the direct sum, so fewer stragglers sum
// directly over the n points for less. History: docs/code-history-20261021.md, heuristics.hpp 14.
inline constexpr size_t program_machine_tree_active = 256;
// The binary64 stage's sums use their own, coarser tree: a node is a far field once its radius is at most this ratio
// of its distance (0.2 for the coefficient routes' trees), its series cut at this many terms, a truncation bound of
// about 1.2e-4 of the node's sum of magnitudes (binary64's trees: 1e-16). Aberth's step damps an error in its sum by
// the square of the Newton correction (Bini 2024, section 2.1), and the sums only steer uncertified points, so a
// coarser sum costs at most sweeps. The pair is the fastest of the ratios and orders measured on deep Mandelbrot
// programs, the coefficient routes' pair among them. History: docs/code-history-20261021.md, heuristics.hpp 15.
inline constexpr double program_tree_acceptance = 0.5;
inline constexpr size_t program_tree_order = 14;
// That rule as the tree's far-field policy (interaction_tree.hpp), a constant in the traversal's code. A ratio read
// from the tree object at each node measured slower on the default trees' inputs, with the same answers and
// counters, so each tree's rule is in its type. History: docs/code-history-20261021.md, heuristics.hpp 16.
struct ProgramFarField {
    static constexpr double acceptance = program_tree_acceptance;
};
// The program route's second version, phase 2a (program_certificate.hpp): the separated-roots pass's rounds before
// it hands over to the refinement, and the rounds it allows without fewer overlapping disks (a cluster, or two
// points at one root, keeps them). A round evaluates only the points that moved or need more precision.
inline constexpr uint64_t program_disk_rounds = 32, program_disk_stalled_rounds = 4;
// Up to this many points to step, the pass takes their Aberth sums in balls over every point (close pairs, which
// binary64 midpoints cannot tell apart); above it in binary64, through the interaction tree at large degree.
inline constexpr size_t program_disk_exact_sums = 64;
// The Aberth steps a round gives each lost point among those few (a culprit whose disk is wider than 2^-20 of its
// point), with the converged points fixed in its sum (Newton's step on the polynomial with their roots divided out).
inline constexpr unsigned program_disk_inner_steps = 64;
inline constexpr bool program_disk_pass = CERTPOLYSOLVE_EXPERIMENT_PROGRAM_DISKS != 0;
// Phase 2b (program_certificate.hpp): a root count on a circle by the argument principle. The circle starts as
// this many arcs, and an arc whose image's enclosures both meet 0 is bisected, down to a 2^program_winding_arc_log2
// of the circle (finer means the circle passes too close to a root at this precision).
inline constexpr unsigned program_winding_arcs = 16;
inline constexpr int program_winding_arc_log2 = -24;
// The arcs a count may evaluate in all: base + scale (k + 2) for a circle expected to hold k roots. An arc's change
// of argument is under a half-turn, so k roots need at least 2k arcs, and the enclosures' widening many times more
// near a root of a deep program; the bound keeps one count that cannot resolve from costing a degree's worth of
// evaluations. History: docs/code-history-20261021.md, heuristics.hpp 17.
inline constexpr uint64_t program_winding_arc_base = 512, program_winding_arc_scale = 16;
// A cluster's first circle, about its points' centroid: this many times the reach of its disks.
inline constexpr unsigned program_cluster_circle = 2;
// For a program with a root of order v at zero, the route solving q = p / x^v: a cluster's count is taken on p
// itself when the circle's centre lies at least this many times v, plus one, radii from zero. Then z^v turns by less
// than a quarter radian around the whole circle, so p reads on q's arcs, and its enclosures are tighter (q's are p's
// divided, and q' carries p's own ball). Nearer zero, p's root there bends the circle's image, and q is counted.
// History: docs/code-history-20261021.md, heuristics.hpp 19.
inline constexpr unsigned program_undivided_count_margin = 8;
// A power output c q^k is solved as its base q and its cover lifted (program_structure.hpp). The ablation build
// (CERTPOLYSOLVE_EXPERIMENT_PROGRAM_LIFT=0) solves the output itself.
// History: docs/code-history-20261021.md, heuristics.hpp 20.
inline constexpr bool program_power_lift = CERTPOLYSOLVE_EXPERIMENT_PROGRAM_LIFT != 0;
// A product of distinct factors c q_1^k_1 ... q_m^k_m is solved factor by factor and the covers lifted and united
// (output_factors). The ablation build (CERTPOLYSOLVE_EXPERIMENT_PROGRAM_PRODUCT_LIFT=0) solves it whole.
// History: docs/code-history-20261021.md, heuristics.hpp 21.
inline constexpr bool program_product_lift = CERTPOLYSOLVE_EXPERIMENT_PROGRAM_PRODUCT_LIFT != 0;
// A squaring recurrence's output is seeded from its base's roots (program_recurrence_seeds.hpp) above this degree,
// at every level of the recursion; below it the circle counts seed.
// History: docs/code-history-20261021.md, heuristics.hpp 22.
inline constexpr bool program_recurrence_seeds = CERTPOLYSOLVE_EXPERIMENT_PROGRAM_RECURRENCE_SEEDS != 0;
inline constexpr uint64_t program_recurrence_seed_degree = 64;
// Points an evaluation at once: the chunk's arrays (slots x chunk x 48 bytes, about 9 slots for a three-term
// recurrence) stay in the second-level cache.
inline constexpr size_t program_machine_chunk = 1024;
// Newton steps of the certificate attempt at one point (machine_newton_cover).
// A step is taken only when the certificate refuses because the point is not
// accurate enough, and a later step only while the correction shrank at least
// fourfold. At a simple root Newton's iteration contracts quadratically; at a
// root of multiplicity m, or a cluster that acts like one, it contracts
// linearly with ratio (m-1)/m >= 1/2. Any bound below 1/2 tells the two
// apart; 1/4 = (1/2)^2 leaves one bit of margin for the corrections' own
// enclosures.
inline constexpr slong newton_cover_contraction_log2 = -2;
// A step is taken only from a correction whose enclosure has this many
// accurate bits: the fourfold rule tests two bits of the corrections it
// compares, and two more keep the enclosure's own width out of that test.
// A less accurate correction is rounding, which the certificate reports.
inline constexpr slong newton_cover_correction_bits = 2 - newton_cover_contraction_log2;
// The first rung of the certificate's Newton ladder (machine_newton_cover): a
// point from the binary64 packet holds at most binary64's bits, and a Newton
// step from b correct bits needs about 2b bits of evaluation to give 2b; so the
// steps below the certificate precision start at twice binary64's mantissa,
// rounded up to whole limbs (FLINT's unit of precision): 128 bits.
inline constexpr slong newton_ladder_first_rung =
    ((2 * std::numeric_limits<double>::digits + FLINT_BITS - 1) / FLINT_BITS) * FLINT_BITS;
// The least raise of the certificate's precision on a rounding refusal (machine_newton_cover's certify), as a
// ratio: half an octave. Raised to the deficit its evaluation shows alone, a point whose need is far above that
// deficit (Newton steps still to take, or a goal of thousands of bits) climbs a limb at a time, a raise and an
// evaluation each; at least half an octave, crossing an octave takes at most two raises, and a raise overshoots
// the need by at most sqrt(2) where a doubling overshoots by up to 2 (an evaluation's cost grows about as the
// precision^1.6 there: 1.74 against 3.03).
inline constexpr double certificate_raise_ratio = 1.4142135623730951;
// Quadratic convergence from one correct bit doubles the correct bits with
// each step, so bit_width(goal bits) steps reach any goal.
constexpr unsigned newton_cover_steps(uint64_t goal_bits) { return static_cast<unsigned>(std::bit_width(goal_bits)); }
// The sweeps of a round of the certificate attempt (a multiprecision Aberth
// packet on the points it refused, every other point held). A point repaired
// alone moves by Newton's iteration on the polynomial deflated by the held
// points, which doubles its correct bits at a simple root; points that move
// together triple them (Aberth). The count below is the tripling one: the
// smallest k with 3^k >= goal bits. It is a quantum, not a guarantee: at a
// 53-bit goal it is four sweeps, which also turn four correct bits into 64
// by doubling; the certificate's own Newton steps do the rest, and a point
// that needs more stays for the next round.
// Proposal 1 (dev/docs/cluster-double-double-design-20261027.md): the certificate attempt's double-double repair of
// the points the compensated tier refuses (CompensatedCertificate::refine), measured on the easy panel's clusters on
// 27 October 2026 (research/cluster-dd/dd_aberth_probe.py).
// - The nudge: each moving point turned by m times this angle (radians, m its rank among them). It must exceed the
//   double-double resolution of the points' differences, so that a conjugate pair stops being one, and stay below
//   their binary64 error, so that the start is not spoiled: in the probe 2^-50 left every 2^-40 pair stalled and
//   2^-20 slowed the rounded double roots; 2^-30 took the most cells soonest.
inline constexpr double double_double_nudge = 0x1p-30;
// - A point has converged when its correction is below 2^-(goal + this) of its modulus: certify's rectangle is the
//   goal's half-width (about 2^-(goal + 3) relative) and the binary64 point it starts from is within 2^-53 of this one.
inline constexpr int double_double_margin_bits = 8;
// - At most this many steps: the probe's slowest converging cells (pairs 2^-40 apart, rounded double roots) took up
//   to 30, almost all of them separating two points that started together.
inline constexpr unsigned double_double_refine_steps = 40;
// - At most this many refused points: each step costs every moving point a double-double evaluation (2 n terms) and
//   a sum over all n points, so the repair is for a few clusters, not for a cover that failed wholesale (the ball
//   certificate and its rounds take those, as before).
inline constexpr size_t double_double_refine_points = 16;
// The cluster start's scope (certificate_attempt.cpp, machine_newton_cover): the double-double repair's, a few clusters,
// not a cover that failed wholesale. The start takes its groups one by one, each a packet whose members interact with
// all n points at the cluster's precision, then their certificates; where many points are refused the pass, the
// rounds and the caller's route take them together (the third gate of the start: clustered-256, 208 refused points,
// its groups certified one by one in 0.56 s where the secular route takes 0.14 s; clustered-512 465, sos-near-axis-512
// 512, against 3 to 6 on the degree-8 and degree-32 clusters the start was built for).
inline constexpr size_t cluster_start_points = double_double_refine_points;
// - Its price in the attempt's budget (machine_newton_cover's operation_weight units, a binary64 complex multiply-add
//   being 1): a double-double complex multiply-add is four double-double products (about 9 operations each, TwoProd
//   with an FMA and the cross terms) and four accurate double-double sums (about 20 each), about 115 binary64
//   operations against binary64's 8, so about 15 times; rounded up. A step at one point is a round's sweep at one
//   point (costs.round_sweep: the value, the slope and n interactions) at this weight.
inline constexpr double double_double_operation_weight = 16;
// - The repair may take at most this share of what the attempt's budget has left, its steps capped to fit and only the
//   steps it runs charged: a repair that fails leaves the rounds that would have run without it at least the other
//   half, and the rounds on a refused set cost at most twice the first (certificate_attempt.cpp).
inline constexpr double double_double_budget_share = 0.5;
// - The reach test. Double-double finds a root of a cluster of k roots and size s to about u^2 / s^(k-1), u = 2^-53
//   (dev/docs/cluster-double-double-design-20261027.md), so a point of such a cluster never meets the stop above when
//   s^(k-1) < 2^-(double_double_bits - goal - double_double_margin_bits): every step it takes is waste (on the easy
//   panel's triangles of radius 2^-30 and crosses of radius 2^-60, a quarter of the solve, 27 October). Binary64's
//   points cannot show such a cluster (they scatter at about u^(1/k)), the double-double ones can after a few steps:
//   from step double_double_reach_first_step, every double_double_reach_every steps, a moving point whose (k-1)-th
//   nearest point lies within s of it (relative to its modulus) with s^(k-1) below that bound by
//   double_double_reach_margin_bits more, for some k up to double_double_reach_group, stops. The margin keeps the
//   test off clusters the repair does reach: a rounded 6-fold root spreads to about 2^-8.8, s^5 = 2^-44 against the
//   test's 2^-48 at the default goal. Replayed on 64 cluster inputs (research/cluster-dd/reach_probe.py): every one
//   of the 81 points that never converge is stopped; the 8 others it stops are of pairs and triangles too close for
//   the tier's own proof (its double-word value error), none of rounded multiple roots or close pairs.
inline constexpr int double_double_reach_margin_bits = 3;
inline constexpr unsigned double_double_reach_first_step = 8; // the probe stopped no point before step 8
inline constexpr unsigned double_double_reach_every = 4;
inline constexpr size_t double_double_reach_group = 8;
// The compensated tier's retries (CompensatedCertificate::certify): where its test fails at a binary64 point x and
// Newton's iteration from x is in its quadratic range (Kantorovich's h <= 1/2), the tier evaluates again at the binary64
// point nearest its Newton step. At high degree binary64's own iteration leaves simple roots about n u cond from their
// points (unisolve_easy3200: a median 2.3e-11, 27 October 2026), too far for the test's Taylor remainder over the step,
// though one step lands close enough. A step in the quadratic range doubles the correct bits, so from a point with at
// least half of binary64's bits one retry starts within binary64's rounding of the root, the closest a binary64 point
// can be: one.
inline constexpr unsigned compensated_newton_retries = 1;
// The compensated tier's smaller circles (CompensatedCertificate::prove). Where the goal's circles fail because the
// second-derivative term, which grows as the radius squared, passes |D| times the radius (a coarse goal: at 10 bits
// the half-width is 2^-10 of the root), the test runs again, from the same evaluation, on circles of a radius rho
// below the goal's. Its upper end is where that term, on the outer circle 2 rho, takes this share of the circle's
// |D| 2 rho: 0.5 K (2 rho)^2 <= share |D| 2 rho, rho <= share |D| / K. A quarter leaves three quarters to the value's
// and slope's errors. rho is the power of two halfway, in exponent, between that end and the least radius those errors
// allow, so that both keep a margin of the same factor; the upper end itself made boxes of 2^-18 of the root on
// rand-uniform-2048 at 10 bits, and led eliminant-katsura-9 into a geometric restart that the middle does not
// (dev/docs/slowdowns-20261029.md, 29 October 2026).
inline constexpr double compensated_circle_curvature_share = 0.25;
// The cluster start's seeds (certificate_attempt.cpp, machine_newton_cover): k points on a circle about the cluster's
// centre at angles 2 pi j / k, turned by this many radians, so that no seed lies on the real axis through a real centre:
// for a real polynomial, at a real point with the other points real or in conjugate pairs, the Aberth step is real
// and the seed would stay on the axis. Half a radian is the probe's choice (research/d8-clusters/start_probe.py, 1 to 4
// sweeps on its 27 clusters), not tuned.
inline constexpr double cluster_start_turn = 0.5;
constexpr unsigned newton_cover_round_sweeps(uint64_t goal_bits) {
    unsigned sweeps = 1;
    for (uint64_t reached = 3; reached < goal_bits; reached *= 3)
        ++sweeps;
    return sweeps;
}
// The certificate attempt's cost in the selector's unit (auto_selection.cpp):
// global sweeps, each degree x (degree + support) of the operations that
// piecewise's estimate counts. The weights below are measured times in units
// of such an operation, its own time measured over large inputs.
//
// Its binary64 packet's first quantum makes at most machine_cover_sweeps
// sweeps, each an evaluation and the root interactions at every point; that
// quantum is priced here (further quanta, while the packet settles points,
// come off the attempt's own budget in machine_newton_cover). A term of a dense
// source's binary64 evaluation and a pair of the direct interaction sum each
// measured about one of piecewise's operations. A sparse source's binary64
// evaluation (at most a quarter of the coefficients nonzero, as
// candidate_jets.hpp decides) takes machine_sparse_term_weight operations a
// term, the measured time of the gap's powers and the jet. From
// machine_interaction_tree_degree on, the interactions come from the tree at
// machine_tree_level_weight operations for each of bit_width(degree) levels:
// the measured time a point and sweep of the packet's tree on large sparse
// inputs, shared among its levels.
//
// Its certificates evaluate value and slope twice at every point. A
// coefficient of a dense source's evaluations (a block dot product) costs
// certificate_coefficient_weight operations, its measured time on a large
// dense input. A sparse source's evaluation (DiskSparsePolynomialEvaluator)
// takes its two dot products over the terms at the same weight, and first
// builds the terms' powers: the products its plan counts,
// certificate_multiplication_weight each, measured on the certificates of
// large sparse inputs.
//
// The weights are unit conversions measured on the reference machine, like
// operation_weight's table (dovetail.hpp): they decide which of two methods
// is tried first, never what either computes. The packet's term is an upper
// bound: a packet that converges stops its points early. The certificates'
// term is what they cost when every point is certified from two evaluations
// at the attempt's first precision, as on every input measured; Newton steps,
// precision raises and repair rounds add to it within the attempt's own
// allowances. History: docs/code-history-20261021.md, heuristics.hpp 23.
inline constexpr uint64_t machine_sparse_term_weight = 10;
inline constexpr uint64_t machine_tree_level_weight = 242;
inline constexpr uint64_t certificate_coefficient_weight = 17;
inline constexpr uint64_t certificate_multiplication_weight = 70;
inline constexpr unsigned certificate_evaluations = 2; // at each point, of the value and of the slope
// The repair rounds of the certificate attempt (machine_newton_cover): a sweep at one point of a round's
// multiprecision Aberth packet makes degree - 1 interactions and evaluates the coefficient enclosure (value and
// slope: two midpoint multiply-adds a coefficient of a dense source; a sparse source's terms as its certificate
// does). In the selector's operations at the attempt's first precision (85 bits at the default goal; an
// operation's time as certificate_coefficient_weight implies it): an interaction's time measured on large sparse
// inputs, and a dense coefficient's multiply-add's, from the time a degree unit and sweep measured on dense inputs
// at a higher precision, less the interactions there, halved, and scaled to 85 bits by operation_weight. Measured
// on the reference machine like the weights above; they bound how much repair the attempt's budget pays for, never
// what is certified. History: docs/code-history-20261021.md, heuristics.hpp 24.
inline constexpr uint64_t newton_cover_interaction_weight = 118;
inline constexpr uint64_t newton_cover_coefficient_weight = 34;
// The packet's and the attempt's costs at one point, exact counts of the selector's operations
// (docs/same-routes-every-platform-note-20261021.md: counts exact, ratios in binary64): the packet's
// machine_cover_sweeps sweeps of an evaluation and the interactions, and the certificates' evaluations. A global
// sweep costs degree x (degree + support) operations, so the cost in global sweeps is the count at a point over
// degree + support, and the whole cost is the count at a point times the degree. For a degree of at least 1.
constexpr uint64_t machine_packet_point_operations(uint64_t degree, uint64_t support) {
    const uint64_t evaluation = 4 * support <= degree ? machine_sparse_term_weight * support : degree + 1;
    const uint64_t interactions = degree >= machine_interaction_tree_degree
                                      ? machine_tree_level_weight * static_cast<uint64_t>(std::bit_width(degree))
                                      : degree - 1;
    return machine_cover_sweeps * (evaluation + interactions);
}
// `sparse_products`: the products one evaluation of a sparse source makes
// (DiskSparsePolynomialEvaluator::plan); zero for a dense source.
constexpr uint64_t certificate_attempt_point_operations(uint64_t degree, uint64_t support,
                                                        uint64_t sparse_products = 0) {
    return machine_packet_point_operations(degree, support) +
           certificate_evaluations *
               (2 * certificate_coefficient_weight * support + certificate_multiplication_weight * sparse_products);
}
// The geometric local restart's hints (secular_refinement.cpp, group_model):
// a node is certified when its inclusion disk, n times its residue, meets the
// goal, so it must lie within about 2^-(goal + log2 n) |z| of its root; the
// residue exceeds that distance by a factor 1 + O(distance / spacing). The
// nested solve that gives the hints is asked for this many bits more, so that
// factor and the midpoints' rounding leave the hints inside the certificate.
inline constexpr double local_restart_hint_margin_bits = 8;
// The packet on a restarted component's local model (component_seeds).
inline constexpr unsigned local_restart_sweeps = 32;
// Cost crossover: that model from the head of the source's expansion (PolynomialSource::taylor_head) in place
// of FLINT's whole Taylor shift. The head's k + 1 coefficients take k + 1 passes of Horner's rule at the working
// precision, and the dropped terms at the hints' midpoints at most k more, where the model's roots are determined.
// FLINT's shift (acb_poly/taylor_shift.c), measured in such passes on dense polynomials at a full-precision centre:
// - its divide-and-conquer recursion, above taylor_shift_horner_length where the precision is at most twice the
//   length: at least a measured share of the length at the shorter lengths, and a measured number of passes at
//   the longer. The head is taken where its 2 k + 1 passes are at most half of both: an eighth of the length
//   (restart_head_length_share, under half of that share) and restart_head_passes (half of that number);
// - one convolution, where the precision is above twice the length: the passes fall with the precision and
//   hardly depend on the length (convolution_shift_passes: the least over lengths 64 to 1,024 at each
//   precision of the grid). The head is taken where the k + 1 passes of an undetermined model, the usual one,
//   are at most half of that; the 2 k + 1 of a determined one are then within it. On a sparse input the
//   convolution can take fewer passes than a short head.
// History: docs/code-history-20261021.md, heuristics.hpp 25.
inline constexpr slong taylor_shift_horner_length = 30;
inline constexpr slong restart_head_passes = 27;
inline constexpr slong restart_head_length_share = 8;
inline double convolution_shift_passes(slong precision) {
    constexpr size_t sizes = 8;
    constexpr std::array<double, sizes> bits{256, 512, 1024, 2048, 4096, 8192, 16384, 32768};
    constexpr std::array<double, sizes> passes{28.4, 25.8, 20.4, 13.6, 11.3, 8.6, 5.3, 3.8};
    const double p = static_cast<double>(precision);
    if (p <= bits.front())
        return passes.front();
    // Past the grid: its last value (at 65,536 bits the grid's passes are no fewer).
    if (p >= bits.back())
        return passes.back();
    size_t upper = 1;
    while (bits[upper] < p)
        ++upper;
    // Log-log interpolation in base 2 by portable::log2 and portable::exp2: the restart decision below compares this
    // value, so it is the same on every platform (src/dovetail.hpp's meter, the same rule).
    const double along = portable::log2(p / bits[upper - 1]) / portable::log2(bits[upper] / bits[upper - 1]);
    return portable::exp2(portable::log2(passes[upper - 1]) +
                          along * (portable::log2(passes[upper]) - portable::log2(passes[upper - 1])));
}
// The head's cost is given in passes (PolynomialSource::taylor_head_passes): of a model whose roots are not
// determined, its k + 1 coefficients; of a determined one, those and the dropped terms at up to k midpoints.
// For a dense source they are k + 1 and 2 k + 1. A source with few terms takes the head term by term, a small
// fraction of a pass, while FLINT's shift does not use the zeros; on a sparse input with a close pair the shifts
// at its precision rungs can be a large share of the solve's work.
// History: docs/code-history-20261021.md, heuristics.hpp 26.
inline bool restart_model_from_head(slong length, slong precision, double undetermined_passes,
                                    double determined_passes) {
    if (length <= taylor_shift_horner_length)
        return false;
    if (precision <= 2 * length)
        return determined_passes <= static_cast<double>(restart_head_passes) &&
               static_cast<double>(restart_head_length_share) * determined_passes <= static_cast<double>(length);
    return 2.0 * undetermined_passes <= convolution_shift_passes(precision);
}
// The Newton steps one point takes in a piecewise polish pass.
inline constexpr unsigned piecewise_polish_steps = 16;
// The default work budget a root (cps_options::work_units_default). A unit is a
// controller step: about one local model, certificate or square of the
// searches. 64 a root is about four times the most that Auto's solves of
// large dense inputs used: a default caller's large solve finishes, while a
// solve that runs away still meets a limit that grows with its input, not a
// fixed one. History: docs/code-history-20261021.md, heuristics.hpp 27.
inline constexpr uint64_t default_work_units_per_root = 64;
// The rounds of an axis-mode attempt (docs/c1-axis-charts-note-20261015.md) before a real-only query goes to the
// complex route. Measured against Auto on the real-roots corpus, every input the axis mode answered faster did so
// within two rounds, and every slower one needed more (close pairs, pairs beside the axis, all-real inputs whose
// coefficients cancel on the axis). Two rounds keep every gain; a failed attempt costs the models of 16 and 32 bits
// up to its first certain failure. History: docs/code-history-20261021.md, heuristics.hpp 28.
inline constexpr uint64_t axis_attempt_rounds = 2;
// B6, the sparse recursion (docs/b6-sparse-real-roots-note-20261015.md). Cost hypotheses from a probe of sparse
// inputs of large degree: evaluations about 0.6 k^2 log2 n, each about k (2 log2(n / k) + 3) counted
// multiplications. One counted multiplication is about 800 of the selector's units (piecewise_operations): a
// multiplication's measured time over the unit's, the unit priced by the axis charts' time. The rule "800 x the
// prediction below the estimate of the route the selector would take" reproduces the measured crossovers against
// the axis charts; where the selector chooses the global recipes the recursion is cheaper by far at any small
// support. In binary64, with portable::log2: the same bits on every platform.
// History: docs/code-history-20261021.md, heuristics.hpp 29.
inline constexpr double sparse_evaluations_per_term_squared_log = 0.6;
inline constexpr double sparse_selector_units_per_multiplication = 800;
inline double sparse_recursion_evaluations(uint64_t degree, uint64_t support) {
    const double n = static_cast<double>(degree), k = static_cast<double>(support);
    return sparse_evaluations_per_term_squared_log * k * k * portable::log2(std::max(n, 2.0));
}
inline double sparse_recursion_multiplications(uint64_t degree, uint64_t support) {
    const double n = static_cast<double>(degree), k = static_cast<double>(support);
    const double log_gap = portable::log2(std::max(n / k, 2.0));
    return sparse_recursion_evaluations(degree, support) * k * (2 * log_gap + 3);
}
// An attempt's ceiling: this many times its predicted evaluations, after which the isolation stops as a work
// limit and the query goes on to the other routes. The count, not the precision-weighted work, tells the cases
// apart: a repeated root or a pair beyond the cap sends a breakpoint's refinement to the cap and the evaluations
// far beyond the prediction, while a pair the cap allows costs a few times the prediction even at a high
// precision, and the fewnomials about it; this ratio lies between them, as measured on sparse inputs.
// History: docs/code-history-20261021.md, heuristics.hpp 30.
inline constexpr double sparse_attempt_evaluation_ratio = 8;
} // namespace certpolysolve::detail::heuristics
