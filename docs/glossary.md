# Glossary and notation

The words the documentation uses, with the plainer phrase that means the same thing and the name the library's
headers and source use. A page for users says the plain phrase; the guide, the specification and the mathematics
say the term, defined here once.

## Vocabulary

| Term | In plain words | In the headers and source | Meaning |
|---|---|---|---|
| dyadic | an exact binary number | `Dyadic`, `cps_dyadic_view` | m · 2^e with integers m and e; every region's ends are dyadic, given as the two decimal integers |
| binary64 | a double | `double` | the IEEE 754 64-bit floating-point format, 53 significant bits; the plain numbers of `values()` and `cps_result_get_approximation` |
| region, box | a box | `Root`, `box` | a closed rectangle in the complex plane with exact ends, holding a proved number of roots |
| complete cover | the boxes hold every root | `complete` | disjoint regions whose proved counts add up to the degree, so every root is in one of them and none twice |
| count | how many roots are in the box | `count` | the number of roots in a region, counted with multiplicity, always proved |
| singleton | a box with one root | `multiplicity` set | a region proved to hold one distinct root, with its multiplicity |
| counted cluster | a box with several roots | `multiplicity` empty | a region with a proved count only; the roots in it may be distinct or repeated |
| goal | what you asked for | `Goal`, `cps_goal` | accuracy (each root alone to the requested bits), approximation (clusters allowed when the centre approximates each root), isolation (real roots each alone, any width), or a cluster radius |
| stop | the solve stopped early | `Stop`, `cps_result_termination` | an early end that returns the proved regions, an unmet goal and the reason |
| honest stop | the solver does not guess | | the rule that a stop returns only what is proved |
| certificate | the proof for a box | `certificate` | the inequality, checked on the original input in exact or ball arithmetic, that proves a region's count |
| ball arithmetic | arithmetic that carries its error | FLINT's `arb`, `acb` | each value is a midpoint and a radius guaranteed to contain the true value |
| working precision | the bits the solver computes with | `work_bits`, `max_bits` | the precision of the ball arithmetic, raised when a proof fails for lack of bits, up to the cap |
| accuracy | how many correct bits you asked for | `bits` | the relative error allowed in each root, as a power of two |
| certificate attempt | the fast first pass | `machine_newton_cover` | the double-precision iteration followed by a certificate at each point, which finishes most solves |
| packet | a bounded run of the iteration | `packet` | a fixed number of sweeps of the iteration before the solver looks at the result |
| secular representation | the polynomial written by its nodes | `SecularEquation`, `secular` | p written as Q·S with nodes bⱼ and residues wⱼ (after MPSolve); exact for any distinct nodes |
| regeneration | new residues at the current nodes | `regeneration` | recomputing the secular representation at the current candidates |
| restart | new positions for a group of candidates | `restart` | moving a cluster's candidates from a model built at the cluster's own scale |
| the secular count | counting roots in a disk by a winding number | `secular_disk_count` | a root count by the winding number of S around a circle, in double precision |
| the two-disk test | a count in a disk becomes a count in a square | | the same count in disks of radius r and 2r proves the count in the square between |
| ring, chart, local model | pieces of the plane with a short stand-in for p | `piecewise` | the piecewise method's geometry: rings around the origin, charts on each, and on each chart a short model of p with a proved error bound (after PWPoly) |
| rung | one step up a ladder | `rung` | one step of a precision or model-accuracy ladder |
| owner | a proved region in the ledger | `PiecewiseOwnership` | a region whose count is proved and kept while the models and the precision change |
| guard | a candidate's first, coarse proof | `guard` | the proof that reserves a root for a candidate before refinement |
| regime exit, hand-off | the piecewise route gives the input to the secular one | `Handoff` | the piecewise route leaving on counted work, storage, range or precision |
| subdivision fallback | a search that splits and contracts root regions | subdivision | a persistent search using proved root counts; its progress and complexity limitations are described in [the mathematics](mathematics.md#persistent-certified-subdivision-and-contraction) |
| dovetail, the share | subdivision shares work with the iteration | `dovetail_share` | estimated arithmetic work schedules subdivision; global Auto engages after precision growth leaves a certified cluster unresolved; the share is not a wall-time fraction |
| the selector | how the solver chooses a method | Auto, `Method::automatic` | the routing by a cost estimate in exact integers |
| controller work units | public counted work | `work_units`, `work_budget` | units charged for controller operations such as a preparation pass, backend call, local attempt or split; distinct from the scheduler's arithmetic estimates and internal allowance |
| the budgets | limits a caller may set | `work_budget`, `cell_budget`, `iteration_budget` | counted work; live cells of the subdivision search; sweeps of the iteration; each stops the solve with its own reason when spent |
| checkpoint | where a solve looks at the clock | | the points between stages where a solve checks its time limit and its cancellation token; a stop lands at the next one |
| Mahler's bound | how close two roots can be | | a lower bound on the distance between distinct roots, computed from the degree and the coefficients; the solver quotes it when a precision cap is too low |
| the three methods | classical, secular, piecewise | `Method` | the Aberth–Ehrlich iteration on the coefficients; the iteration on the secular representation; the local-model route |
| program | a polynomial given by how to compute it | `Program`, `cps_program` | input as a recorded list of operations, never expanded into coefficients |
| real-only query | only the real roots | `real_roots` | the query that returns every real root, each proved real |
| square-free | no repeated roots | `Squarefree` | proved or disproved by exact arithmetic on the input, the only way a multiplicity is ever learned |

## Notation (the mathematics page)

| Symbol | Meaning |
|---|---|
| p, n, aⱼ | the input polynomial, its degree and its coefficients, p = a₀ + a₁z + … + aₙzⁿ with aₙ ≠ 0 |
| ζⱼ | the roots of p, listed with multiplicity |
| D(c, r), D̄(c, r) | the open and the closed disk of centre c and radius r |
| Sq(c, r) | the closed square of centre c and half-width r |
| N(p, U) | the number of roots of p in U, with multiplicity |
| wind(γ, 0) | the winding number of a closed curve about 0 |
| b | the requested accuracy in bits: relative error at most 2⁻ᵇ |
| P | the working precision in bits |
| [m ± r] | a ball with midpoint m and radius r |
| \|A\| lower, \|A\| upper | a lower and an upper bound on \|A\| over its ball |
| bⱼ, wⱼ | the nodes and residues of the secular representation |
| Q, S | Q(z) = ∏ⱼ (z − bⱼ) and S(z) = 1 − Σⱼ wⱼ/(z − bⱼ), with p = Q·S |
| hⱼ, eⱼ | log₂\|aⱼ\|, and an integer with \|aⱼ\| < 2^eⱼ |
| F, E, ε | a local model, its error, and the bound \|E\| ≤ ε |
| R, γ, M | a chart's radius, a ring's middle radius, a ring's number of charts |
| W, A | the piecewise work estimate and the certificate attempt's cost |
| φ | the golden ratio (1 + √5)/2 |
