# The mathematics: what each certificate proves

This page states what CertPolySolve proves about its answers, by which certificate, and where the trust boundary lies.
It is written for a reader who wants the precise claims; [how it works](how-it-works.md) tells the story first, and
[the glossary](glossary.md) defines the words and the notation.

The trust boundary: FLINT supplies the arithmetic (exact integers and rationals, and ball arithmetic with rigorous
error bounds), the exact polynomial algebra and the low-level inclusion predicates. The root finding, the real
isolation and the certificates are CertPolySolve's own, and no FLINT root solver is called. A numerical iteration, a
deflation or an eigensolve only proposes candidates; every certificate is checked against the original input, and a
candidate that stops iterating is a proposal, never a proof. The code is not formally verified.

The sections follow a solve: the input domain; why coverage and accuracy are separate claims; how the solver
chooses its method; the certificates for one root; the first pass; local counts and recovery; the certified
subdivision and the conservative fallback; what is exported and with what accuracy; the independent validation;
secular input; the piecewise method; programs; real roots; and the limits of what is promised.

## Input domain

Coefficient input is a polynomial in the monomial basis, given by its nonzero
terms in ascending order; the routes that read every coefficient expand it to the
dense list. A coefficient is an exact real/complex value, or a Cartesian product
of real intervals specified by exact centers and nonnegative radii. All
arithmetic enclosures must contain those exact values/intervals. A certificate
for uncertain input applies uniformly to every polynomial in the coefficient
family, not just to its midpoint polynomial.
The other two inputs are native secular equations (below) and programs:
straight-line programs over exact constants, never expanded, whose degree,
leading coefficient, root bound and square-freeness verdict are proved from the
program itself (the program-input specification,
sections 4-7).

Numbers are kept as `rational_significand * radix^integer_exponent`, with radix 2
or 10. This separates exponent encoding from coefficient magnitude. A value such
as `1@100000000000000000000` need not allocate that many bits. Decimal powers are
evaluated with outward ball arithmetic and recreated at each working precision.

Exact algebraic preprocessing is used only for real exact inputs whose conversion
fits an explicit budget. Both coefficient expansion and common-denominator
storage are checked. The budget does NOT bound subsequent gcd/factorization
intermediate growth. If conversion is ineligible, the ball path remains available.

## Coverage is distinct from accuracy

For a nonzero polynomial of proved degree n, complete all-complex coverage means:

1. Every returned region has a certified root count, counting multiplicity.
2. Returned regions are pairwise disjoint, including after endpoint export.
3. Their counts sum to n.

No residual threshold establishes any of these facts by itself.

For `CPS_ROOT_ACCURACY`, the default goal for coefficient and secular input,
`goal_met` additionally requires separation of each distinct root and the
requested accuracy. An exact repeated root can meet this goal when its
multiplicity is proved. A small cluster of distinct or unresolved roots cannot
be relabeled as one repeated root.

`CPS_ROOT_APPROXIMATION` (opt-in for coefficient and secular input, the default
for programs) retains the same complete disjoint cover but
allows counted clusters to meet the requested accuracy. If a rectangle has
midpoint `c` and component radii `r_x,r_y`, then for **every** point `z` in it,
`|z-c| <= r_x+r_y`. The larger coordinate absolute lower bound `L` satisfies
`L <= |z|`. Thus `r_x+r_y <= 2^-bits L` proves the requested true-relative
error for every root in the region, independently of root separation or
multiplicity. The explicit absolute tolerance is an alternative sufficient
bound. A nontrivial box containing zero cannot pass the relative check merely
by being small in absolute units. Exact zero factors remain exact point roots.
The same predicate covers input uncertainty uniformly. A real-only query answers
the approximation goal with individual roots, as for accuracy.

An explicit cluster-radius query relaxes distinct-root separation: `goal_met`
then requires complete coverage and a requested absolute upper bound on each
exported region's Euclidean radius about its midpoint. It does not establish a
cluster's multiplicity or promise the 3-disk separation of a Ccluster natural
cluster. The result records which goal was requested. Default root queries retain
the original separation-and-accuracy meaning.

A real-only query concerns only the real roots, so its counts need not sum to the
full degree. It is answered for exact real input; ball and nonreal input are
refused. Every region of its answer is proved real, and when its goal is met each
holds one distinct real root with its multiplicity.

- **Coefficient input, Auto.** Three routes answer, in order, each stopping
  where it cannot certify and handing the query on, with its work counted. The
  sparse recursion goes first on an input with few terms, where its predicted
  cost is below the selector's estimate for the complex query: certified signs
  along the chain of fractional derivatives, in time counted in terms, not
  degree; it stops where consecutive derivatives may share a root, over its
  evaluation ceiling or on a cancellation (the B6
  note). Then, where the selector would
  choose the piecewise method, its axis charts: only the charts on the real axis
  are built, and the real roots are isolated from their real models in one
  dimension (the axis-chart note); an input with
  repeated or power factors ends it at once. Both attempts need a work budget of
  at least 64 units a root. Last, the complex route: a complex solve at the
  accuracy goal, whose regions proved real by conjugation symmetry are the
  answer. A region on the axis that is not proved real makes the route solve
  again at higher accuracy, until it is proved real or not real, or the cap, the
  work or sweep budget or a cancellation stops it. Where the complex route stops
  at its precision cap the query stops too; for exact real input the stop names
  a cap that suffices (Mahler's root-separation bound). What each route proves
  is in "Real roots of real input".
- **Secular input, Auto.** The equation is converted to its exact numerator
  polynomial, within the exact limits, and its real roots are isolated by the
  exact real engine below ("Exact real isolation and refinement").
- **The isolation goal** asks for each real root alone, with its multiplicity, in a
  closed interval of any width; the complex route then solves at the library's
  default 53 bits.
- **A closed region** [a, b] keeps the real roots x with a <= x <= b. A root at an
  end is decided by an exact test of the input there. An interval of a root of odd
  multiplicity that holds an end it is not equal to is placed on its side by the
  input's signs at the end and at the interval's lower end (the input changes sign
  once in the interval, at the root), and cut at a dyadic point between the end and
  the root; a root of even multiplicity is refined until its interval lies on one
  side.

## Routing: the selector

Under Auto, query dispatch first handles exact-coefficient complex isolation
through the complete prepared subdivision route, and tries eligible closed
forms for accuracy/approximation at degrees one through four. Real-only queries
use the separate real route described above. For the remaining coefficient
queries, the selector (`select_auto_route`, `src/solve/auto_selection.cpp`)
chooses the maintained secular or piecewise recipe, or the family-aware D5
recipe (`production_policy`). It reads the query and input alone, never a name,
hash or family label. Native secular input uses D5, and program Auto uses the
classical recipe without the coefficient selector. The
[method table](production-methods.md#the-selector) gives every supported
combination. Every route certifies its output against the original input.

- **Coefficient families.** Inexact coefficient input uses D5's family-aware
  certificates (`inexact-coefficients`). Isolation requires exact coefficients;
  unsupported combinations stop with the input/query outcome in the method table.
- **The precision cap.** Every other query takes the secular recipe (`global`)
  unless piecewise is admitted. Piecewise starts at `max(128, 2*16+b+64)` bits
  for a goal of `b` bits (`piecewise_start_precision`); under a lower cap
  `max_work_bits` it does no root work and its precision stop hands the query
  to secular, so such a cap takes the secular route at once
  (`precision-cap-below-piecewise-start`): a shortcut, not a choice.
- **Piecewise's admission** (`admit_piecewise`, `src/piecewise_admission.cpp`)
  runs with its development caps on the envelope's bend and the work ratio
  lifted (`lifted_policy`), at the goal raised to the default if below it,
  `b' = max(b, 53)` (`judged_query`: the estimate drops model rungs as the goal
  falls, while neither route's cost does). It refuses the cluster-radius goal,
  an unusable precision or a cancelled call, a constant nonzero part, a power
  composite (`p=x^v q(x^k)` with `k>1`, read from `support_stride`) and
  geometry its estimate cannot plan; none of these is a statement about the
  roots, and each leaves the query to secular.
- **The work estimate.** For each model rung `m = 16, 32, ...` up to
  `max(16, bit_ceil(2b'))` bits, `estimate_piecewise_models` walks piecewise's
  rings without building a model (the geometry: `src/piecewise_models.cpp`'s
  header). A ring with band `l..u`, `d = u-l`, `K` charts (a power of two) and
  Taylor degree `L = min(d, 4(m+16))` adds `K(L+1)^2 + (L+1)(d+1+K log2 K)`
  complex coefficient operations (local solves, model construction). Their
  total `W`, exact and saturating at `2^64-1`, is a cost estimate, not a bound.
- **The comparison.** With `n` the degree of `p/x^v` and `s` the number of
  its terms, a global Aberth sweep is priced at `n(n+s)` operations (`n`
  evaluations of `s` terms, `n^2` interactions) and a global solve at
  `initial_sweep_allowance = 256` sweeps, the least the classical candidate
  stage allows at one precision (`max(256, 2n)`, `src/classical_isolation.cpp`).
  The estimate passes when `W < 256 n(n+s)`, a cost rule (the selector
  campaign, "Offline analysis: the rule"). The
  test is exact: for positive integers `floor(floor(W/n)/(n+s))` is
  `floor(W/(n(n+s)))`, at least 256 exactly when `W >= 256 n(n+s)`, and no
  product that could pass `2^64` is formed. The walk stops once its total
  reaches the bound (`estimate_above_global_sweeps`); as the totals only grow,
  stopping early changes no verdict.
- **The square-free proof.** Last, for Auto's own route, the nonzero part must
  be proved square-free (`piecewise_squarefree_modular`) by a one-sided gcd
  test like that of "Individual root certificates", but modulo the split prime
  `998244353` (1 modulo 4, so `i` maps to a square root of `-1` in `F_p`;
  `src/modular_squarefree.cpp`). A failed test proves nothing and leaves the
  route secular; a proof chooses piecewise
  (`piecewise-estimate-below-global-sweep-budget`).
- **The certificate attempt first.** The global recipes begin with a
  certificate attempt (`machine_newton_cover`: a binary64 Aberth packet, then a
  certificate at each point; it succeeds when the `n` regions, each certified
  to hold a root, are pairwise disjoint). Where the selector chooses piecewise
  it may put that attempt first (`attempt_first`,
  `PiecewiseRootConfig::certificate_attempt`). Its price `A`, in the
  estimate's unit, is `n` times `certificate_attempt_point_operations`
  (`src/heuristics.hpp`): the packet's first quantum (`machine_cover_sweeps`
  sweeps, an upper bound) and the certificates' two evaluations of value and
  slope at each point. The attempt goes first when `W >= phi A`,
  `phi = (1+sqrt(5))/2` (`certificate-attempt-before-piecewise`).

  In an idealized two-attempt model, this order minimizes the worst-case ratio.
  Suppose the attempt costs `A` and
  succeeds or fails, which is not known beforehand, and piecewise costs `W` and
  answers. Compare each order with the cheaper order in hindsight, which costs
  `min(A,W)` if the attempt succeeds and `W` if it fails, and put `x = W/A`.
  The attempt first costs `A` or `A+W`: its worst ratio is
  `max(1, 1/x, 1+1/x) = 1+1/x`. Piecewise first costs `W`: its worst ratio is
  `max(1, x)`. Now `1+1/x <= max(1, x)` exactly when `x > 1` and
  `x^2-x-1 >= 0`, that is when `x >= phi`. So the rule takes, for every `x`, an
  order with the smaller worst case, and that worst case,
  `min(1+1/x, max(1, x))`, never exceeds `phi`, its value at `x = phi`.

  This argument concerns that model: `A` and `W` are estimates, and admitted
  operations can overrun an allowance. The implementation passes `W` as the
  attempt allowance (`certificate_attempt_fallback`, `machine_newton_cover`'s
  `fallback_operations`); it does not thereby prove a whole-call competitive
  ratio or a wall-time bound. A forced piecewise solve makes no such preliminary
  attempt. The comparison is in binary64;
  `tests/test_auto_selection.cpp` restates it exactly as `W^2 >= W A + A^2` in
  128-bit integers (`x >= phi` is `x^2 >= x+1` for `x > 0`).
- **All-real input.** On exact real coefficients, where the selector chose
  piecewise for the accuracy or approximation goal, the solver runs the
  all-real seed pass on the primitive integer polynomial (`every_root_seeded`,
  `src/real_line_routing.cpp`; degree 64 and up, `real_line_seed_min_degree`).
  If it gives a seed for every root, the secular route is taken instead, from
  those seeds (`every-root-seeded-before-piecewise`). The verdict is a cost
  signal: nothing is certified from it, and the secular route certifies every
  root (the note).
- **Real-only queries.** `RealOnlyRoute` (`src/real_only_route.cpp`) asks the
  selector for the complex query's route, without the square-free proof, and
  tries the sparse recursion where its predicted cost is below that route's
  estimate (`W` for piecewise, `256 n(n+s)` otherwise) and the axis charts
  where that route is piecewise, before the complex route, a fresh solve
  under the selector: cost rules too ("Coefficient input, Auto", above).
- **The same route on every platform.** Every value that decides a route, an
  order or a cap is an exact integer (saturating `resources::add` and
  `resources::multiply`, exact floor and ceiling quotients, 128-bit products)
  or a binary64 value of `+ - * /` (the golden ratio, measured prices, and
  `portable::log2` in the sparse recursion's prediction), built without
  contraction; nothing uses `long double` or the platform's transcendental
  functions. So a query takes the same route on every platform
  (the note); a binary64
  comparison departs from the exact statement only at a near-tie.

A wrong route therefore costs time only, and piecewise limits that time: its
regime exit (`regime_exit_due`, `src/piecewise_roots.cpp`) hands the original
input to secular, within what is left of the call's budget, once its counted
work passes what `W` stands for plus a global solve's minimum budget.

## Individual root certificates

The exact zero roots come first. The `v` lowest coefficients that are exact zeros
are set apart as the point 0, with count and multiplicity `v`, and the rest of
the preparation works on the nonzero part `p/x^v`: in the global driver
(`SolveState::Impl::preprocess`) as in piecewise's preparation. On exact real
input the one-sided modular test below (`squarefree_modular`) then runs first, as
on exact nonreal input: when it proves the nonzero part square-free, that part's
primitive form with a positive leading coefficient is its one factor, the factor
the integer decomposition returns for a square-free input, and no gcd over Z is
computed. Otherwise exact real polynomials are squarefree-decomposed. Each
squarefree factor is solved and its roots receive that factor's proven
multiplicity. Regions belonging to different factors, and the point 0, are
checked for overlap before claiming overall isolation.

Exact nonreal coefficients have a budgeted Q(i) path (`exact_complex_factors`).
Two one-sided squarefreeness tests come first: the image modulo the split prime
998244353, which is 1 modulo 4, so that `i` maps to a square root of -1 in F_p
(`squarefree_modular`); then the image in F_p[i], p = 2^31-1, which is 3 modulo
4, so that x^2+1 is irreducible. With denominators and the leading coefficient
nonzero after reduction, gcd(f,f')=1 modulo p proves squarefreeness in
characteristic zero. A failed modular test does not establish a repeated root;
it selects the exact characteristic-zero squarefree decomposition instead, whose
gcds over Z[i] are multimodular (Chinese remaindering and rational
reconstruction, a candidate accepted only when it divides both inputs).

Every exact polynomial division checks its remainder. The resulting monic
factors, raised to their multiplicities, must reconstruct the normalized input
exactly before they can be used. The modular shortcut retains the original
polynomial itself with multiplicity one. As on the real path, certified factor
regions must be disjoint and their weighted counts must sum to the original degree.
Tests reconstruct factors using separate Boost rational arithmetic, check bad
prime/denominator cases, and independently count repeated algebraic roots.

The Q(i) path uses separate configured exact-materialization and exact-operation
allowances, plus a per-polynomial intermediate-storage estimate. Unspecified
allowances have no fixed policy ceiling; compact exponent expansion separately
requires admission relative to owned input size (growth of at most one precision
cap a number) or explicit caller permission.
See [the resource policy](resource-policy.md) for the estimates and zero semantics.
Individual rational operations are preflighted for uncancelled operand growth.
These are conservative work guards, not hard process time or RAM limits. A guard
or cancellation discards the entire decomposition; no partially computed factor
information is promoted. Uncertain inputs never enter this exact path. Compressed
exponents that cannot fit the budget remain on the ball path.

For complex solving, native bounded Aberth packets compute numerical candidates.
Low-level FLINT inclusion predicates validate those proposals against the
coefficient balls. We accept a complete vector only when its isolated-root count
equals its degree. A smaller count does not certify completeness. Real-root
classification in a complex solve uses validated conjugation symmetry, never an
imaginary tolerance; the real-only route's sparse recursion and axis charts
prove realness by a sign change in a one-root interval ("Real roots of real
input").
The relevant pinned contract is in [FLINT 3.6 acb_poly documentation](https://raw.githubusercontent.com/flintlib/flint/v3.6.0/doc/source/acb_poly.rst).

Linear equations use ball division. Provably zero initial coefficients are set
apart before the exact preparation, as the point 0 (above); no approximate
deflation is performed. The all-root
count and disjointness checks also prevent double counting a zero root that may
appear in an uncertain reduced polynomial.

Candidate centers are only hints. Finite, distinct centers can be reused at higher
precision; invalid/coincident centers trigger fresh initialization. Iteration
allowances scale with working precision unless the user supplies an explicit
limit. Neither a candidate center nor its convergence is a certificate.

## The certificate attempt

`machine_newton_cover` (`src/certificate_attempt.cpp`) is the first pass that the
maintained recipes share: classical's isolation makes it once on each factor
that is neither linear nor an analytic binomial (`isolate_complex`), secular's
refinement once per goal (`newton_cover`), and the selector's piecewise route
once, before any chart, on an input it solves whole (`certificate_attempt` in
`src/prepared_piecewise.cpp`), where piecewise's estimate is at least the golden
ratio times the attempt's cost bound; a forced piecewise solve never makes it.
A binary64 Aberth-Ehrlich packet proposes `n` points, and a certificate at each
point proves a region about it. It serves a complex query's accuracy and
approximation goals on an expanded coefficient source (no recurrence, no
program) of degree `n>=1` whose leading coefficient excludes zero; `p` below is
the monic quotient it evaluates, and every bound holds uniformly over the
coefficient family. A refusal claims nothing: the caller's route continues.

### The packet and the certificate at a point

The points start as the caller's or on the Newton polygon's circles
(`initial_candidates`; near each edge's binomial roots on a sparse source).
After `prepare_machine_engine` scales them into binary64's range, a Gauss-Seidel
sweep moves each active point by Aberth's correction `N/(1-N*S)`, `N=p/p'`,
`S=sum_{j!=i} 1/(z_i-z_j)`, or by the same correction as `1/(p'/p-S)` where `p'`
is zero or `N` or `N*S` leaves binary64's range (`src/aberth_step.hpp`, with the
table of every kernel that takes the step). A point stops when its relative step
is below `2^-45`, or `|p(z)|` below the estimate `32(n+1)2^-53 sum_j |a_j||z|^j`
of Horner's rounding. Sweeps come in quanta of `machine_cover_sweeps` (32), cut
short by the call's iteration budget; another runs only while the last stopped a
point, others still move, and the budget (below) pays for it beside one
certificate pass. An arithmetic failure (an invalid, divide-by-zero or overflow
flag, a disallowed underflow, a non-finite or coincident point) restores the
quantum's starting points; in the first quantum, or with no engine, the packet
then runs with a separate exponent (`ScaledBinary64`), but before piecewise the
attempt is refused; a later failure ends the packet. Nothing here is certified.

A point first tries the compensated tier (`src/compensated_certificate.cpp`). It proves the same statement more
cheaply, wherever the attempt evaluates the dense coefficients (not a sparse source). The coefficients are held once
as two binary64 words each, scaled by a power of two, and the rest of each ball becomes an error bound. At the
binary64 point `x` nearest the point, `p(x)` comes from the compensated Horner scheme. Each step
`s_(k+1)x+a_k=s_k+t_k` is split exactly by error-free transformations (TwoProd by a fused multiply-add, TwoSum), so
`p(x)=s_0+E(x)` with `E(x)=sum_j t_j x^j`, and a second Horner recurrence evaluates `E`. `p'(x)` is the derivative
recurrence on the rounded `s_k`. Every error is bounded in binary64, for every polynomial in the coefficient balls.
In the normal range an operation is within a relative `2^-53`; below it a sum is exact and a product is off by at most
`2^-1075`. No floating-point flag is read.

One Newton step gives the centre `z1=x-q`, exactly as two words a part. Rouche's theorem then compares `p` with the
line `V+D(z-x)`, `V` and `D` the computed value and slope, on the circles `|z-z1|=r` and `|z-z1|=2r`. Here `r` is the
goal's half-width, and `|p''|` is bounded by the next paragraph's majorant `M` on the disk that holds both circles.
This proves exactly one root in the disk of radius `r` and no other in the disk of radius `2r`. The published square
of half-width `r` about `z1` must lie inside the square of half-width `5r/4` and meet the goal. The majorant's term
grows as the radius squared, so at a coarse goal it can fail these circles however close the point is. The test then
runs once more, from the same evaluation, on circles of a radius `rho<r`: the power of two halfway, in exponent,
between the least radius the errors allow and `|D|/(4M)` (`compensated_circle_curvature_share`). `M` still bounds
`|p''|` on their disk, and the square of half-width `rho` meets the goal as well. The tier declines a
goal of 106 bits or more, a half-width below binary64's normal range, a value or bound that overflows, and every
failed test.

Where it declines, the ball certificate runs. Each point is replaced by its exact midpoint `c`. Above 128 bits of certificate
precision, Newton steps first run on rungs doubling from 128
(`newton_ladder_first_rung`); each needs a correction with 4 accurate bits and,
after the first, at most a quarter of the last (a simple root's quadratic
convergence, not a cluster's linear one); a `b`-bit goal allows `bit_width(b)`
steps (`newton_cover_steps`), which propose. At the certificate precision
`p(c)` and `p'(c)` are enclosed. Let `c'` be `c` or the midpoint of its Newton
point `c-p(c)/p'(c)`, and `h` an upper bound of `|c'|` times `2^-(b+5)`, or an
eighth of the absolute tolerance if larger. The square of centre `c'` and
half-width `h` is published first, with exact outward endpoints ("Export and
accuracy"); `d` bounds from below the distance from `c'` to its sides. It must
lie inside a disk `|z-c|<=R`: `R=2h` for `c'=c`, and for the Newton point the
step's length (with a slack of `2^-32`) plus twice the half-width at `c`,
rounded up. There Taylor's integral remainder bounds `|p(z)-p(c)-p'(c)(z-c)|`
by `E=(R^2/2)M(|c|+R)`, with `M(s)=sum_{j>=2} j(j-1)|a_j|s^(j-2)` the
coefficient majorant of `p''` (`linear_remainder`, rounded up; over a sparse
source's support). The test, left side bounded below and right above, is

`|p'(c)| d > |p(c)+p'(c)(c'-c)| + E`.

On the boundary `|z-c'|>=d`, and `p(z)-p'(c)(z-c')` is `p(c)+p'(c)(c'-c)` plus
a remainder of modulus at most `E`, so it is smaller than `|p'(c)(z-c')|`. By
Rouche's theorem `p` has as many zeros inside as `p'(c)(z-c')`, exactly one, and
none on the boundary. At the Newton point the constant term vanishes but for
rounding, so the same evaluation certifies there when `E`, which grows with the
square of the step, is below `|p'(c)|d` (`certify_at_step`). The test is
conservative: `E` holds on the whole disk and ignores the cancellation in `p''`.

A failed test is classified by the first of `E`, the constant term's rounding
radius and the constant term to reach a quarter of the left side. A large
constant term (or `|p(c)|>=|p'(c)|h`, seen before any remainder) sends the point
to its Newton point's square, then to a step and a new evaluation. Rounding
doubles the certificate precision; at the cap, or where doubling no longer
shrinks the value's radius, no later round can help and the attempt ends. A
large `E` (or a square that never fits its disk) falls back on a disk needing no
model: since `p'/p=sum_j 1/(z-xi_j)` over the `n` roots, a root lies within
`n|p(c)/p'(c)|` of `c` ("Certificates" under "Program input"). The square about
`c` of half-side `n` times an upper bound of `|p(c)|` over a lower bound of
`|p'(c)|`, published, is taken when it meets the goal. The rectangle proves
exactly one root; the disk at least one. Without the coefficient bound, which
every maintained recipe sets, `certify_newton_point` serves instead.

### The cover, its repairs and its partial certificates

A complete cover is `n` regions pairwise disjoint as closed published
rectangles (`disjoint`). Each holds at least one of the `n` roots counted with
multiplicity, so each holds exactly one, which is simple, and none lies outside
them: the complete-vector rule of "Individual root certificates". The attempt
claims no realness. Classical and secular attach the factor's exact
multiplicity and recheck what they publish; piecewise adds any exact zero root
and, for real input, proves realness by reflection.

The pass tries the points in turn; a certified point moves to its region's
centre and is held there. Refused points are repaired in rounds: a
multiprecision Aberth packet (`candidate_packet`) of `newton_cover_round_sweeps`
sweeps (the least `k` with `3^k` at least the goal's bits) at the certificate
precision moves them, the others fixed, and each is certified again. A round on
two or more points must certify at least half of them or the attempt ends; a
lone point may stay. A raise of precision is kept only where it led to a
certificate. Refusals may spend one evaluation a point at the first precision
in all, and, once a point is certified, a point's raises only what is left.

A refusal after some certificate returns the certificates to secular
(`PartialNewtonCover`): the points, the certified regions less those that meet
another (two that meet may hold one root), the certificate precision, and
whether the pass reached every point (`pass_complete`). Each region holds at
least one root. Secular keeps them only from a complete pass
(`keep_certificates`), each meeting the goal and no other kept region; the node
stays in its region, inactive, and the validation (`secular_inclusion_cover`,
"Persistent interpolation certificates") gives a component whose members are
all kept their regions in place of its disk, unless a region then meets another:
disjoint regions each holding at least its count, the counts summing to the
degree, hold exactly their counts. After a goal change with a cover held, the
attempt runs without a packet on the points whose regions miss the new goal
(`settled`) and returns only such certificates. A refusal leaves classical's
packets the packet's points, and secular's preliminary phase too when the
packet converged and nothing was kept (`newton_point_handover`).

The spending is counted in the selector's coefficient operations, priced by
`newton_cover_costs` for one certificate evaluation of value and slope and one
sweep at one point of a round. The budget is `fallback_operations`, the cost of
the route a failure leaves the solve to: piecewise's estimate before piecewise,
otherwise `initial_sweep_allowance` (256) sweeps of `n(n+t)` operations for `t`
terms. The packet's quanta come off it, and the balance, as evaluations at the
first precision weighted by `operation_weight`, is charged with the certificate
evaluations at theirs; a round takes only the refused points it still pays for.
The prices decide what is tried, never what is proved.

## Local counts and candidate-guided recovery

When individual isolation fails at the outer precision limit, or during an
explicit cluster-radius query, a bounded local
recovery stage tries to certify groups of candidate roots. Candidates and their
Newton inclusion radii are only search hints. Duplicate or misleading hints do
not establish counts or completeness.

For an exact dyadic center c and positive dyadic radius r, form the ball-enclosed
Taylor polynomial q(w) = p(c + r w). If one coefficient satisfies

\[
\inf |q_k| > \sum_{j \ne k} \sup |q_j|,
\]

Rouche's theorem gives exactly k roots, counting multiplicity, inside the disk
|z-c| < r and none on its boundary. All inequalities are strict, outward ball
comparisons. Failure is **inconclusive**, never a zero-root claim. Up to four
Graeffe transforms are tried; root squaring preserves the number inside the unit
circle and the absence of roots on it. This is a sufficient Pellet/Graeffe test,
not a new theorem or a reproduction of the full Ccluster algorithm. See
[Imbach–Pan–Yap (2018), section 3](https://arxiv.org/pdf/1806.10584) and the
[FLINT 3.6 transform contract](https://raw.githubusercontent.com/flintlib/flint/v3.6.0/doc/source/acb_poly.rst).

A disk count alone does **not** establish the count in its bounding square:
additional roots may be in the corners. The implementation certifies the same
count in disks of radii r and 2r, publishes an outward square containing the inner
disk, and checks that the published containing ball is strictly inside the outer
disk. Thus the actual exported rectangle has exactly the asserted count.

Candidate groups merge when their search rectangles overlap; failed tests may
increase a search radius. Lower Taylor terms propose radius changes, with
neighbor-aware growth to permit nearby hints to merge before swallowing distant
roots. These are heuristics whose outputs must face the complete count predicate.
Only a disjoint certified partition whose weighted counts sum to the factor's
degree is accepted. A multi-root cluster has unknown multiplicity; exact factor
multiplicities are attached only when a factor region contains one root.

This candidate-guided recovery has explicit work guards: the private, overridable
`cluster_recovery_degree` crossover (default 2048), 16 geometric
rounds, and at most min(4096, max(64, 8*degree)) box attempts. It is not a
subdivision algorithm. The persistent subdivision fallback below can continue
from certified coverage if candidate recovery fails.
Among equally structured complete covers, a candidate with a worse largest
enclosure radius is not allowed to replace an earlier useful bound. This quality
ordering is a presentation choice, not an extra mathematical certificate.

## Persistent certified subdivision and contraction

`subdivision.cpp` implements a separate search that can run without numerical
root hints. Exact preprocessing may still provide factors and analytic linear
solutions. Each factor starts from a complete certified cover; dyadic squares
contain its rectangles. Expanded seed squares that overlap, including at corners,
are merged before separate root sets are assigned. Factor multiplicities remain
immutable, established only by exact preprocessing.

Splitting replaces a closed square by its four closed dyadic children. A child
is discarded only when a **strict** Pellet/Graeffe zero-count test succeeds on a
disk of radius 3/2 times its half-width, which strictly contains its corners.
Inconclusive tests retain the child. Eight-neighbor connectivity merges both
edge-touching and corner-touching cells. Exact dyadic centers must be representable
at the current working precision plus a guard; otherwise the search increases
precision or stops without discarding cells. Thus roots on subdivision boundaries
cannot fall through a gap.

The retained frontier covers every root of the factor. A component receives an
arbitrary positive count only after a containing disk count and the two-disk
rectangle certificate succeed, with the outer square disjoint from all other
component supports. That separation ensures its certificate cannot silently
include another component's root set. Whole-factor publication additionally
requires disjoint exported rectangles and the full weighted degree count.

For a component with certified count k, modified Newton proposes
`c - k*p(c)/p'(c)` and a smaller dyadic square. These approximate values supply
only a proposal. The proposal is snapped to an exact dyadic grid and must lie
inside the square backing the parent's certificate. Strict inner/outer disk
tests must prove the same count k. Equal counts in nested regions establish
that the same root set has been retained; only then may contraction replace the
parent's cells. Geometric component bounds and the square backing its certificate
are stored separately. Failure falls back to splitting or higher local precision.
An exact zero evaluation at an exact center, together with a strict unique-root
count, can collapse that root to an exact point, including zero.
For a strict count k around the exact center zero, exact vanishing of the first
k coefficients proves divisibility by x^k. All k counted roots must therefore
be the same zero, so the cluster can be published as an exact multiple point
(with any already established factor multiplicity). This inexpensive proof
also applies when squarefree preprocessing was unavailable. A coefficient ball
merely containing zero does not prove divisibility.

The controller schedules factors round-robin and complex components through a
queue. Each component can request a fresh enclosure from the immutable input or
exact factor at its own precision. A soft magnitude comparison distinguishes a
clear geometric failure from unresolved interval width; it is only a scheduling
hint. Only strict dominance establishes a mathematical count. For exact coprime
factors, overlapping exported rectangles are merged with summed counts; a tighter
local target separates the conflicting factor covers when the requested root
goal requires it. An incomplete replacement frontier never overwrites an older
complete cover.

The design follows the subdivision/Pellet/Graeffe/Newton ideas of
[Becker et al., root isolation](https://arxiv.org/abs/1509.06231),
[Becker et al., root clustering](https://arxiv.org/abs/2105.05183), and
[Imbach, Pan and Yap's Ccluster implementation](https://arxiv.org/abs/1806.10584).
The complete execution uses guarded Newton proposals, temporary local trial
precision and retained geometric progress. Under standard bit-arithmetic and
convergent enclosure assumptions, the reviewed correspondence gives isolation
cost Õ(n³+n²τ) for one prepared square-free Gaussian-integer factor, and an
additional Õ(n²b) for direct refinement of its retained isolating state. Here
n is the factor degree, τ its coefficient height and b the requested accuracy.
Exact preparation, cross-factor separation, lifting, assembly and output carry
separate charges. The model and correspondence
state the hypotheses and material deviations from the paper.

For direct exact-coefficient all-complex Auto queries on the global route,
a single bounded numerical attempt is followed by exclusive complete service.
With successful exact preparation and nonbinding caller and implementation
limits, this composes the complete engine's progress into eventual attainment.
The whole-call account is the bounded attempt plus the complete continuation,
including the surrounding charges. Standalone exact complex isolation uses
complete service immediately. The [method table](production-methods.md)
distinguishes these routes from explicit methods and staged piecewise, program
and real wrappers, which retain their own scheduling rules.

These are algorithm/model statements, not a bound on seconds or a formal proof
of the implementation. The coefficient-family contract and arbitrary warm
forests do not acquire the prepared-factor bound. A binding precision, work,
cell or storage limit can still prevent attainment. Counted output promises
its count and requested accuracy; it does not promise a natural-cluster
separation condition.

## Conservative global fallback

If the leading coefficient excludes zero, Cauchy's bound gives

\[
|z| < 1 + \max_{k<n}|a_k/a_n|
\]

for every root. To see this, put M = max |a_k/a_n|. For |z| > 1 + M,
the sum of the lower terms' magnitudes divided by |a_n z^n| is strictly less
than M/(|z|-1), hence less than 1; cancellation to zero is impossible.

We compute an outward upper bound R and return the containing square
`[-R,R] + i[-R,R]` with count n. Because every admissible polynomial has degree n
and every root is in this square, the count is uniform over coefficient uncertainty.

This gives safe coverage when separation and local recovery fail. It remains a
**coarse global cluster**, distinct from local Pellet/Graeffe certificates. Its
multiplicity is unknown when n > 1.

The global driver publishes this square once the degree is proved and assembles
the factors' covers once every factor has one. At a stop short of the goal that
is not a cancellation (a budget, the precision cap, a stall), a factor that has
had no step yet is completed on its own polynomial in `z` (a power composite's
reduced polynomial spread back to `z^k`): `complete_cover_at_stop`, the
smallest square about the origin, of a power-of-two half-side, whose count a
disk test certifies, or the factor's own Cauchy square. Its counts are
multiplied by the factor's multiplicity, the factors' regions are merged, and the
cover replaces the published one when it ranks ahead (`prefer_cover`). A
cancellation, a time limit included, completes nothing. Piecewise's factor loop
completes an unfinished factor the same way ("The axis mode and the prepared
pipeline").

## Export and accuracy

Public enclosures are closed rectangles whose four endpoints are exact dyadics.
Endpoint conversion rounds outwards. We then form containing balls and repeat
disjointness and accuracy checks, so serialization cannot silently invalidate
a certificate. Finite endpoint precision avoids expanding a very large exponent
gap into an enormous integer.

For a containing rectangle with midpoint c and component radii r_re, r_im,
let E = r_re + r_im. Then |z-c|_2 <= E for every point in the rectangle.
Let L be the larger certified lower bound for |Re z| and |Im z|. Then
|z|_2 >= L. Therefore

\[
E \le 2^{-b} L
\]

is a sufficient true-relative b-bit error test. This is conservative but avoids
confusing relative-to-midpoint accuracy with relative-to-root accuracy. Exact
zero roots are handled separately. If supplied, a positive absolute tolerance
can alternatively establish E <= tolerance. Precision refers to root accuracy,
not simply the mantissa width used by a backend.

`cps_result_get_approximation` explicitly produces an ordinary binary64 midpoint.
It is not an arbitrary-precision export and does not itself carry the enclosure's
error bound. The exact endpoints remain available. Unresolved clusters, unmet
goals, isolation answers, overflow, and underflow of a nonzero component are
rejected. The C++ `Root::value` is the region's centre rounded to the nearest
double for every region, refused in no case; it is a root's value only where the
goal is met and the multiplicity is established.

## Independent validation

The test-only reference implementation uses Boost exact integers/rationals and
does not call FLINT. It provides polynomial Euclidean arithmetic, squarefree
decomposition and Sturm counts. It verifies counts in the actual exported closed
intervals, including endpoint roots and repeated factors.

Complex tests additionally use an exact sufficient Rouche/Pellet predicate:
Taylor-expand p(c+w), and prove one term strictly dominates the sum of all others
on a disk contained in the exported rectangle. Exact 1-norm/max-norm bounds avoid
square roots. Combined with disjointness and the total degree, successful checks
independently establish coverage for those test cases.

These are independent mathematical checks, not formal proofs of the C++ program.
The test checker is intentionally unsuitable as a production high-degree solver.

The internal disk-count tests construct complex rational roots, form their
polynomial using independent Boost arithmetic, and compare each certified disk
count with exact squared-distance comparisons for those known roots. They cover
boundary roots, Graeffe-dependent counts, uncertain coefficients, invalid disks,
duplicate hints and the disk-versus-square corner trap.

A separate exact Sturm–Cauchy checker uses signed Euclidean remainder sequences
along rectangle edges to compute the winding number of the polynomial image.
It does not use Pellet/Graeffe or numerical roots. Boundary zeros give an
inconclusive result; degenerate rectangles are treated as exact points or line
segments. See [validation and replay](testing.md) for its domain and the completed
302,418-case property campaign. Sampled coefficient families are tests of the
implementation, not a replacement for uniform ball inequalities.

## Secular inputs and evaluation groundwork

For distinct exact poles `b_j` and nonzero residues `a_j`, let
`Q(x)=product_j(x-b_j)` and `P(x)=sum_j a_j Q(x)/(x-b_j)-Q(x)`.
The public secular input solves `P=0`. The leading coefficient is -1 and the
degree is the number of terms. At a pole,
`P(b_i)=a_i product_{j!=i}(b_i-b_j)` is nonzero. Thus clearing denominators
introduces no roots at excluded poles. Empty input gives `P=-1`.

The coefficient fallback uses exact Gaussian rationals and a streaming numerator
recurrence. A second reconstruction divides the finished denominator by each
linear factor, checks the exact remainder, and compares every numerator
coefficient. Term-count and arithmetic-growth limits are explicit; failed
conversion publishes no invented cover. Independent tests use the direct
product-of-all-except-one identity and known roots rather than that recurrence.

The internal `Evaluator` contract returns outward value and derivative enclosures
over the entire supplied complex rectangle and coefficient family. False means
unresolved evaluation (including a pole-containing secular argument), never
root exclusion. Dense/sparse coefficient and secular evaluators implement it.
The internal `MonicSource` contract additionally proves a fixed monic polynomial
degree. Native secular input evaluates `-P=Q(1-sum a/(x-b))` and its derivative
with a streaming product/numerator identity, including polynomial limits at
poles. It does not treat the rational function itself as a polynomial oracle.
The benchmark adapter still expands the audited MPSolve recurrence exactly at
the input boundary. The recurrence evaluator once enabled by recognizing
`p0=1`, `p(k+1)=1+x*p(k)^2` from the coefficients is gone; a recurrence is given as a
program. Program (straight-line
program) input is public and has its
own `MonicSource`, `ProgramSource`, which evaluates the program over dual numbers of
complex balls (see Program input, below); there is no public black-box
evaluation input.

### Persistent interpolation certificates

For exact, distinct dyadic interpolation nodes `c_i`, form outward enclosures of
`w_i=-p(c_i)/product_{j!=i}(c_i-c_j)` for the original monic polynomial. The
characteristic polynomial of `diag(c)+w*1^T` is exactly `p`. Gershgorin row
disks, centered at `c_i+w_i` with radius `(n-1)|w_i|`, lie inside the exported
disks centered at `c_i` with radius `n|w_i|`. Each connected component retains
its exact eigenvalue count: along `diag(c)+t*w*1^T`, for `0<=t<=1`, every
Gershgorin disk stays inside the same fixed enlarged cover, so no eigenvalue
can cross between its disjoint components. At `t=0` each node contributes one.
We merge the outward exported rectangles to a fixed
point, including contacts on closed boundaries, before publishing a disjoint
complete cover. A multi-node component does not prove root multiplicity.
Conjugation and the complete cover can prove realness of a unique root.

Zero regenerated residues retain the corresponding exact roots. They must not
be removed like invalid zero residues in the original rational input. A source
evaluation, interval bound or proposed node failure produces no new certificate.
For uncertain coefficients the interpolation identity and interval arithmetic
hold for every represented polynomial, so numerical midpoints cannot shrink
input uncertainty.

When at most one quarter of the nodes move, an unchanged node can use
`w'_i=w_i*product_changed (c_i-old_c_j)/(c_i-new_c_j)`. Insufficient accuracy
requests a full original-source evaluation. Updates commit transactionally.
Candidate precision and regeneration precision are independent; all are capped.
Candidates are certified before precision changes or restarts, and a changed
node set is regenerated before its new secular representation is used.

Component restarts may locate a centre through `p^(k-1)`, then initialize a
truncated local degree-k Taylor polynomial with coefficient-scale radii. Those
local roots are proposals only: every resulting count still uses the original
full source. Finished components remain frozen while unresolved roots refine.

### Separate axis-chart experiments

The local polynomial evaluator caches dyadic-centred Taylor polynomials and
discards terms only with outward bounds for both `sum |a_j| R^j` and
`sum j |a_j| R^(j-1)`. These tail bounds are added to value and derivative balls;
nonfinite evaluations fall back to the direct coefficient evaluator. It is the
research build's, with the Series test engine that selects it.

The candidate-only interaction tree uses direct arbitrary-precision neighbours,
self exclusion, and finite far-field Cauchy expansions. Range loss falls back to
direct sums. Its error estimates are not root certificates. Direct/tree A/B
uses identical Jacobi ordering; the native secular path can use Gauss--Seidel.
Native secular initialization can use
eigenvalues of strongly coupled 2-by-2 principal blocks of `diag(b)+a*1^T`.
Both initializers are followed by original-source certification.

## The piecewise method

The third maintained method for exact coefficients follows Imbach and Moroz's
piecewise approximation, without their complexity bounds; the derivations are in
the certification note. Charts cover where
roots can lie, `p` is a Taylor model with a rigorous error bound on each, and
roots are proposed and certified chart by chart, in rounds; geometry, proposals
and scheduling only steer. Exact zero roots, of order `v`, are split off as the
point 0 with count and multiplicity `v`; below, `p` is the nonzero part, of
degree `n`, coefficients `a_j`.

### Charts and models

The envelope (`PiecewiseCoefficientEnvelope`) is the Newton polygon, the concave
upper hull `h_j` of the integers `e_j` with `|a_j| < 2^e_j`. `for_each_ring`
cuts a window of root moduli into rings of at most an octave, each with the band
`l..u` of terms within `m = model bits + 16` bits of the largest, and makes none
where one term stays `m` bits above the rest. A ring `r_1 <= |z| <= r_2` has `S`
charts of radius `R = 3(r_2-r_1)/4` centred at `c = gamma zeta^k`,
`gamma = (r_1+r_2)/2`, `zeta = exp(2 pi i/S)`, `S` the least power of two at or
above `4 pi gamma/R`; none meets 0. Planned in binary64, this geometry sets cost
only: outside the axis mode no skipped stretch is certified root-free, nor needs
to be. Each Taylor row of a ring's models is one length-`S` DFT in balls. For a
scale `2^s` and `lo..hi` on one side of `l`, `bound(lo,hi,r)` bounds
`2^-s sum_(lo<=j<=hi) |a_j| r^(j-l)` by the number of terms times a power of two
at or above their largest `2^(h_j-s) r^(j-l)`, `log2 r` rounded outward. With
`2^s` above the band's unscaled bound at `gamma` and `F` the Taylor polynomial
at `c` of the band's terms over `z^l 2^s`, to the order kept,
`p(z)/(z^l 2^s) = F(w) + E(w)` on `z = c + R w`, `|w| <= 1`, with
`|E| <= bound(0,l-1,gamma-R) + bound(u+1,n,gamma+R) + T`, `T` bounding the
Taylor tail. Rows run to `L = min(u-l, 4m)`; if `L < u-l`, `T = M 2^-L` by
Cauchy's estimate, `M = bound(l,u,gamma+2R)` bounding the band on `|w| = 2`,
else `T = 0`. Under `adaptive_taylor` (maintained) the rows stop after the first
`t < L` with `q_t = rho (u-l-t-1)/(t+2) < 1` proved, `rho = R/gamma`, and
`T = sigma_t/(1-q_t) <= 2^-m bound(l,u,gamma-R)`, `sigma_t` the magnitude sum of
row `t+1`'s terms, each later row's terms being at most `q_t` times the last's.
Under the maintained flags `root_model_error` takes on a rectangle, where
smaller, the omitted terms at its own radii (`local_omitted_band_bound`) and
`T u^k`, the tail starting at order `k` and `|w| <= u <= 1` there.

### Certificates, counts and proposals

`certify_root_model` takes a chart point `w`, `|w| < 1`, and works about `z0`,
the exact midpoint of `c + R w`, at `w0 = (z0-c)/R`. The square about `z0`
starts at the goal's half-width (`goal_half_width`: `|z0| 2^-(b+5)` for a goal
of `b` bits, or an eighth of an absolute tolerance if larger), is halved on
failure at most `precision` times, and must meet the goal after export, lie in
the chart and, if `l > 0`, exclude 0. Let `r_min` be the distance from `z0` to
the exported edges and `r_max` that to the farthest point of the rebuilt ball,
over `R`, `f_k = F^(k)(w0)/k!`, `F_j` the coefficients of `F` and `J` the jet
order. A count `m` is proved by
`inf|f_m| r_min^m > sum_(k<=J,k!=m) sup|f_k| r_max^k + rem + eps(box)`, with
`rem = r_max^(J+1) sum_(j>J) C(j,J+1) sup|F_j| (|w0|+r_max)^(j-J-1)` Taylor's
integral remainder for the terms past the jet and `eps(box)` the model's error
on the rectangle. On the boundary `p/(z^l 2^s)` then differs from `f_m (w-w0)^m`
by less than its modulus, so by Rouché `p` has exactly `m` roots inside. `J` is
`min(deg F, 3)` (`certificate_jet_order`), or `min(deg F, 16)`
(`guard_jet_order`) under an absolute tolerance, as for a guard. The accuracy
goal tries `m = 1`, one simple root; the approximation goal `m = 1..J`, a
counted cluster for `m > 1`.

A chart's count (`count_root_model`) is Pellet's test at radius one: if the
largest coefficient has `inf|F_d| > sum_(j!=d) sup|F_j| + eps`, `eps` the
uniform error, `p` has exactly `d` roots in the open chart. No certificate uses
it: a zero count skips the chart, and one made up by owned regions inside
(`contained_in`) stops fresh proposals there. Two steps only propose:
`proposals`, a partial Ehrlich–Aberth iteration in binary64 on the model's
midpoints from `seed_count` points (closed forms in degrees one and two; skipped
where `chart_excluded` finds no root by Pellet after root-squaring; under
`interior_proposals`, on an inner disk still covering the ring), and `polish`,
at most `heuristics::piecewise_polish_steps` Newton steps on `F` in balls.

### Ownership, rounds and the exit

A candidate's first proof is a guard, the same test under `guard_options`
(accuracy goal, absolute tolerance `2R/(deg F + 1)`). `PiecewiseOwnership` keeps
the proofs disjoint with stable identities: a new one `identifies` an owner
(nested, same count), encloses all it meets with at least their counts (a
cluster), or is refused, changing nothing. Proposals in another candidate's
count-one region or a published root's guard are dropped. A final certificate
meeting a published one is withheld: published regions are disjoint, each holds
exactly its count, and the cover is complete when the counts and `v` sum to the
degree. Once owned singletons count `n` (`discovered`) every root is owned;
unowned candidates are dropped. For real input a nonreal certificate also proves
its conjugate; charts below the axis are skipped.

A round of `PiecewiseRootState` streams the models at the current model bits
(from `piecewise_initial_model_bits`, 16) and precision (from
`piecewise_start_precision`, `max(128, 2 model bits + b + 64)`). At a batch
boundary the source finisher (`begin_finishing`, `advance_finisher`) may spend
the suspended batch's construction work refining owned singletons on the input
`P = z^v p` itself: a Newton step kept in the owner's region, then the same test
on `P(x+Rw) = P(x) + R P'(x) w + E` at the new point `x`, `R = |x| 2^(4-b)`,
`|E| <= (R^2/2) sum_(j>=2) j(j-1)|A_j|(|x|+R)^(j-2)`, `A_j` the coefficients of
`P`; the proof must identify its owner. At a round's end deflated discovery
(`seed_missing_roots`) may seed missing roots by Aberth's correction against the
owned roots (Newton's step on `p` divided by them); then come a recovery (every
chart proposes again), a round more while candidates improve, or the next rung:
the precision doubles if more certificates failed on rounding than on the
model's error, else the model bits double and the precision follows, up to the
cap. The regime exit (`regime_exit_due`) is a cost rule with no clock. With `W`
the selector's estimate ("Routing: the selector"; 0 where it would not choose
piecewise) and `U = initial_sweep_allowance N (N + support)` a global solve's
minimum budget, `N` the degree of the whole input (the selector's `n` is that of
`p/x^v`) and `support` its number of nonzero coefficients (the selector's `s`;
`src/prepared_piecewise.cpp`), the state leaves once its counted work
exceeds `binary64_operation_work (estimate_operation_overhead W + U)`. Then, as
after a storage or range refusal or a precision stop, a fresh secular solve
takes the original input and the budget left; certificates do not transfer. A
spent budget or a cancellation stays with piecewise; in the complex mode every
stop but a cancellation returns a complete cover where it can (`cover_at_stop`).

### The axis mode and the prepared pipeline

Auto's real-only attempt on exact real input (`axis_only`, for
`heuristics::axis_attempt_rounds` rounds) builds only the charts centred on the
axis (the axis-chart note). There
`f = p/(z^l 2^s)` is real for real `w`, within `eps` of `Re F` on a piece of the
segment `|w| <= 13/16`, `eps` its `root_model_error`; the slope's error is that
of the piece widened by `1/8` (in the chart), times 8 (Cauchy). With `Re F` in
the Bernstein basis of its degree `K` at `precision + K` bits (`axis_classify`),
a piece holds no root if all coefficients exceed `eps` in one sign, one simple
root if every scaled difference exceeds the slope's error in one sign and the
end signs differ; where second differences fix the sign of `f''`, a bound at
Newton's minimum proves it root-free or counts its roots by the end signs;
others are bisected to depth `precision`. Coverage is proved: neighbours are cut
where both models decide one sign; a gap is an annulus where one term `a_i`
dominates at both radii, hence between (`annulus_root_free`: `|a_j| r^(j-i)` is
convex in `log r`); the constant and leading terms dominate inside and out. Each
root is refined by the model's signs, then the input's, until its segment lies
in its one-root region (`axis_finish_root`): count 1, proved real. No finisher,
deflated discovery or exit runs; a factored input ends it.

`PreparedPiecewiseState` first runs `prepare_squarefree_input`: the zero order,
the nonzero part as `q(z^k)` (`k` the exponents' gcd, `q(0) != 0`), and a
one-sided modular proof that it is square-free (`squarefree_modular`) or its
exact square-free factors, each `q_j(z^(k_j))` to a multiplicity `mu_j`. A
proved part of degree 3 or more with `k = 1`, or an unavailable preparation
(where a repeated root can only be a cluster), goes whole to one state, after
"The certificate attempt" (`machine_newton_cover`) where the selector asks.
Otherwise each factor of degree 3 or more gets a state, smaller ones the outward
closed forms. Only a complete cover of `q_j` is lifted (`lift_power_cover`),
each region, off 0, to its `k_j` branch images in squares of half-side four
times their larger radius; counts and proved multiplicities are multiplied by
`mu_j`. A lifted region holds at least its count and a complete factor's counts
sum to its degree, so disjoint regions hold exactly theirs. When every factor's
cover is complete (an unfinished factor is completed at a stop,
`complete_cover_at_stop`), regions that meet merge, counts added and the
multiplicity dropped (`merge_overlapping`, as the global driver merges: such a
counted cluster can meet the approximation goal), and the cover is published; an
assembly with a factor still unfinished (a cancellation, a refused lift)
publishes nothing, so the last complete cover stays, the input's Cauchy square
before any. When every factor meets its goal and the mapped cover does not, only
the factor whose region is the wider of the first meeting pair gets a tighter
goal, its accuracy doubled by at least 16 bits (`separation.hpp`, the global
driver's rule), or, where no covers meet, a power factor whose lifted regions
miss the goal 32 more bits; that factor alone is advanced again. Last come the
square-free verdict and, on a complete cover of real input,
`prove_real_by_reflection`.

## Program input

A program is the third input kind: a straight-line program over exact constants,
which the library never expands over the integers or the rationals. Its degree,
leading coefficient and root bound are proved from the program itself, its roots
are found from its own evaluations, and every certificate is checked against
those evaluations in ball arithmetic. Stages that find roots supply proposals
only. The program-input specification states the
contract; the seed-placement note, the
second-version note and the
structural-seeding note give each
stage's design and measurements.

### Programs and their arithmetics

A program is a list of nodes. Each node is the variable `x`, an exact constant
(real and imaginary parts exact numbers, as for coefficients), or `a+b`, `a-b`,
`a*b`, `-a` or `a^e` with `e>=2`, whose operands `a` and `b` are earlier nodes;
one node is the output, and its value is a polynomial `p` in `Q(i)[x]`. Each
node carries an upper bound of its degree: 1 for `x`, 0 for a constant, the
larger for a sum or a difference, the sum for a product and `e` times for a
power, saturating. The interpreter evaluates the nodes in order and releases
each value after its last reader, so a three-term recurrence holds a few values
at once, not one for every node.

Every operation is a ring operation, so running the program in another
commutative ring computes the image of `p` there (dual numbers, residues modulo
a prime ideal); other arithmetics truncate or bound it. The library runs it in
five:

- **Dual complex balls.** Pairs `(u,u')`, with `(u,u')*(w,w')=(uw,u'w+uw')` and
  `(u^e)'=e*u^(e-1)*u'`, in complex balls at the precision asked, starting from
  `(z,1)` for a ball `z`: outward enclosures of `p(z)` and `p'(z)` over every
  point of `z`. This is the only evaluation a certificate reads.
  `ProgramSource`, the program's `MonicSource`, divides both by the exact
  leading coefficient.
- **Top coefficients.** Each value as its degree and its highest `m`
  coefficients, exactly over `Q(i)`, with whether they are all of them. A sum
  whose leading coefficients cancel loses that many; a value with none left is
  undetermined, and the run is repeated with `m` doubled.
- **A modular image.** The coefficients reduced modulo a prime ideal of `Z[i]`
  (below).
- **A majorant.** Every constant replaced by an upper bound of its modulus and
  every subtraction by an addition, evaluated at a real `r>=0`. By induction
  over the nodes the value `M(r)` bounds `sum_i |a_i| r^i` over the node's
  coefficients `a_i`: the coefficients of a sum are bounded by the sums of the
  operands', and the majorant of a product by the product of the operands'
  majorants.
- **Machine arithmetic.** binary64, or binary64 with a separate exponent for
  values outside its range, at many points at once, with a noise estimate per
  point. It serves proposals only.

### What a solve proves before it starts

- **The degree and the leading coefficient** come from the top coefficients at
  widths 1, 2, 4, ... up to a limit, exactly. The program is refused, with a
  diagnostic, when it is longer than the exact-operation allowance, when a
  constant exceeds the exact storage allowance, when it is identically zero, or
  when its leading terms still cancel at the widest series. A program proved a
  nonzero constant has degree 0 and no roots.
- **A root bound.** With `n` the degree and `a_n` the leading coefficient, every
  root has `|z|<1+max_{i<n}|a_i|/|a_n|<=M(1)/|a_n|` (Cauchy). Below that, the
  least power of two `R` with `M(R)<2|a_n|R^n` is sought by halving steps and
  bisection. That inequality gives `sum_{i<n}|a_i|R^i<|a_n|R^n`, and since
  `sum_{i<n}|a_i|r^(i-n)` decreases in `r` it holds at every larger radius, so
  no root has modulus `R` or more. Where the leading terms cancel the test may
  hold nowhere, and Cauchy's bound stands. The square about the origin of
  half-side the bound, with count `n`, is the global cover: the answer of a
  solve stopped before anything finer.
- **Every root at zero.** Let `D` bound the bits of a common denominator `Q` of
  the coefficients, computed up the program from the constants (a sum's or a
  product's `Q` is the product of the operands', a power's the operand's to the
  exponent), so that a nonzero coefficient has `|a_i|>=1/Q>2^-D`, and let
  `|a_n|<2^L`. If `M(R)<2|a_n|R^n` at `R=2^-(D+L+1)`, then for `i<n`,
  `|a_i|<|a_n|R^(n-i)<=|a_n|R<2^-(D+1)`, so every such `a_i` is zero and
  `p=a_n x^n`. The answer is the point 0 with count and multiplicity `n`.
- **The order of a root at zero.** The program run over truncated power series
  at zero, with complex-ball coefficients, encloses the lowest coefficients of
  `p`. A coefficient whose ball lies inside the disk of radius `2^-D` about zero
  is zero, since a nonzero one is larger; one whose ball excludes zero is
  nonzero. The order `v` is proved when the first `v` coefficients are proved
  zero and the next nonzero. The precision doubles from 64 bits while a
  coefficient is undecided, far enough to cover the terms it cancels, and the
  coefficients are widened 1, 2, 4, ... at a time up to the widest series;
  undecided at the ceiling or past that width, the preparation claims nothing.
  `v=n` is the case above.
- **The square-free verdict.** For a prime `l=1 (mod 4)`, `-1` has a square root
  `iota` modulo `l`, and `(l,i-iota)` is a prime ideal of `Z[i]` with residue
  field `F_l`. Running the program there gives the image `f` of `p`, for a fixed
  sequence of primes just above `2^62`, skipping a prime that divides a
  constant's denominator or the leading coefficient, so that the degree is kept.
  If `gcd(f,f')=1` then `p` is square-free: a repeated factor of `p` over
  `Q(i)`, its denominators cleared, would reduce to a repeated factor of `f` of
  the same degree. A test that fails proves nothing; one more prime is tried,
  then the verdict is unknown. The image is priced against the exact-operation
  and temporary-storage allowances before it runs and skipped beyond them.
  Nothing but the diagnostics depends on the verdict.

### Structural identities

The program's graph is read syntactically, from the nodes the caller wrote, with
no rule for any named input. Exact identities then let the route certify a
smaller program and publish the output's roots from it:

- **Powers.** When the output is `c*q^k` with `k>=2`, read through powers,
  squares `q*q`, negations and products with constants, the zero sets of `p` and
  `q` coincide, each root of `p` having `k` times its multiplicity in `q`. The
  route certifies `q` and publishes its cover with every count, and every proved
  multiplicity, multiplied by `k`.
- **Products.** When the output is `c*q_1^k_1*...*q_m^k_m` with `m>=2` distinct
  factors of positive degree, the roots of `p` are the factors' roots, each
  `k_j` times as often. Each factor is solved alone. The lifted covers are
  published together only when every factor's cover is complete, the proved
  degrees add up to `n`, and regions of different factors are pairwise disjoint:
  the union is then a complete disjoint cover of `p`. Otherwise (a shared root,
  for instance) the output is solved whole. A solve stopped before every
  factor's cover is complete answers with the global cover, since a factor not
  yet solved may have roots inside another factor's regions.
- **A root of order `v` at zero.** The route solves `q=p/x^v` through the
  evaluations `q=p/z^v` and `q'=(p'-v*p/z)/z^v`, which are not finite on a ball
  holding zero, and publishes the point 0 with count and multiplicity `v` beside
  `q`'s cover; when a region of that cover holds zero, it counts the `v` roots
  too, and its multiplicity is unknown.

### Starting points

Every starting point is a proposal: a wrong one costs sweeps, never a wrong
answer.

- **Circle counts.** For `s` equally spaced samples `z_g=z_0*w^g` with
  `w=exp(2*pi*i/s)` on a circle about the origin, none of them a root,
  `(1/s) sum_g z_g p'(z_g)/p(z_g) = sum_j 1/(1-(xi_j/z_0)^s)` over the roots
  `xi_j`. A term is 1 up to about `|xi_j/z_0|^s` for a root inside the circle
  and 0 up to about `|z_0/xi_j|^s` for one outside, so the sum is about the
  number of roots inside, with doubt only from roots near the circle (Pan 2022,
  Theorem 8). A ladder of such counts, the radius halving from the root bound
  until no root is inside and the fullest annuli then bisected in log radius,
  says how many roots each annulus holds, and that many seeds go on its middle
  circle at the placement angles; seeds that no count places go on the bound's
  circle. The counts are rounded and clamped to be monotone in the radius, and
  never certified. They replace the Newton polygon that places the coefficient
  route's start.
- **A composition.** When every path from `x` to the output passes through one
  node, whose value is `g(x)`, the output is `f(g(x))` for the program `f` above
  that node, and `p=c*product_s(g(x)-s)` over the roots `s` of `f` counted with
  multiplicity: the roots of `p` are exactly the preimages under `g` of the
  roots of `f`. `f` and `g` are solved roughly, and each `g(x)=s` is started by
  a Newton step from `g`'s roots.
- **A squaring recurrence.** When the output is `F(x,q)` for an earlier node
  `q`, with `F` of degree `d>=2` in `q` and coefficients of low degree in `x`
  (Mandelbrot's `p_(k+1)=x*p_k^2+1` has `d=2`), near a simple root `r` of `q` we
  have `q(z)=q'(r)(z-r)+O((z-r)^2)`. The output therefore vanishes near the
  points where `q'(r)(z-r)` is a root `y` of `F(r,.)`, and each such root gives
  the seed `r+y/q'(r)`. `q` is solved roughly, by the same rules or by circle
  counts, and the seeds its roots leave unplaced go on the root bound's circle.

### The preliminary phase

From the seeds, Aberth's iteration runs on the program's own evaluations. It
runs in machine arithmetic first, as Jacobi sweeps
`z_i <- z_i - N_i/(1-N_i*S_i)`, with `N_i=p(z_i)/p'(z_i)` from one batched
evaluation of every active point and `S_i=sum_{j!=i} 1/(z_i-z_j)` over all
points. A point stops as converged when its relative step falls below `2^-45`;
at a floor when its value or derivative cannot be read above its noise estimate
(the rounding of the last sum on each path, a lower bound of the error by
design); and as held when its values leave the range of both machine
arithmetics. Over each window of about `n/2` sweeps the largest step must halve,
or every point still moving is put at its floor. Stopped points stay in the
others' sums. Packets of Aberth sweeps in balls then continue the points that
did not converge, the converged ones fixed in the sums. A derivative ball
without correct bits at a point means the working precision cannot evaluate the
program there, and the precision rises. Nothing in this phase certifies.

### Certificates

Every certificate reads the dual-ball evaluations of `p`, or of `q` for a root
at zero, at exact points or on balls. None depends on how the points were found.

- **Newton inclusion disks.** Since `p'(z)/p(z)=sum_j 1/(z-xi_j)` over the `n`
  roots counted with multiplicity, some root lies within `n|p(z)/p'(z)|` of `z`.
  At each exact point the radius is `n` times an upper bound of `|p(z)|` over a
  lower bound of `|p'(z)|`, which needs a derivative ball excluding zero
  (otherwise the precision rises, up to the cap). Each disk is published as a
  square rounded outward. When the `n` squares are pairwise disjoint
  (`disjoint`), each holds at least one root and together they hold `n`, so each
  holds exactly one and every root is simple: a complete cover with every count
  and multiplicity 1. A disk disjoint from the others proves nothing alone: two
  points converged to one root leave another root free to lie in some disjoint
  disk.
- **Root counts by the argument principle.** For a circle that no root meets,
  the number of roots inside is `1/(2*pi)` times the change of `arg p` around
  it. The circle is cut into arcs with dyadic ends. On an arc with midpoint
  `z_m`, half-length `h` and a ball `B` about `z_m` holding the arc, the image
  of the arc lies in the closed disk about `p(z_m)` of radius `|p'(B)|*h` (the
  mean-value inequality along segments in `B`), in the box `p(z_m)+p'(B)(B-z_m)`
  and in the enclosure `p(B)`. Each is compact and convex. If one misses 0, `p`
  has no zero on the arc and the argument changes along it by strictly less than
  a half-turn, so the change is the principal value of `arg(p(z_b)/p(z_a))` at
  the arc's ends; the implementation accepts it when its radius is below `2^-3`
  radians. An arc not resolved is bisected, down to a smallest arc. When every
  arc is resolved, the circle meets no root and the summed changes divided by
  `2*pi` hold a unique integer, the count. The count fails, rather than guesses,
  when a value at an arc's end has no argument at the working precision or the
  arcs run out.
- **Counted clusters.** Points whose disks keep overlapping are grouped into
  components. Each component gets a circle about its points' centroid that holds
  its members' squares, and a certified count `m>=1` on that circle; the circle
  is then shrunk about the same centre as far as the goal requires, keeping the
  count `m` (a smaller concentric circle with the same count holds the same
  roots). A component's region is a square holding its circle's closed disk, so
  it holds at least `m` roots. When the regions and the remaining Newton squares
  are pairwise disjoint and their counts add up to `n`, every region holds
  exactly its count and no root lies outside them. A cluster's region claims its
  count, and a multiplicity only when the count is 1. Counted clusters answer
  the approximation and cluster-radius goals; the accuracy goal needs each
  region to hold one distinct root with its multiplicity proved, which a counted
  cluster does not.
- **Counts for a root at zero.** With a root of order `v` at zero, a circle
  whose centre lies more than `8v+1` radii from zero is counted on `p` itself
  rather than on `q`. Its disk excludes 0, so `p` and `q=p/x^v` have the same
  zeros in it, and the count is `q`'s. Any disk excluding 0 would do for
  correctness; the margin keeps to circles where the enclosures of `p` are
  tighter than those of `q`, whose derivative carries `p`'s value divided by
  `z`.
- **The fallback.** What the separated-roots pass does not certify within its
  rounds goes to the first version's refinement (`SecularRefinement`) on the
  same source, from the preliminary phase's points. Its covers are the
  persistent interpolation certificates above, with residues regenerated from
  `n` program evaluations and `n^2` node products.

A region off the real axis is reported non-real, which holds for any
polynomial. A program whose constants are all real has real coefficients: they
are sums, differences and products of real numbers. Its non-real roots then come
in conjugate pairs, and the conjugation argument of the coefficient route
applies to the complete disjoint cover the finish publishes (after the point 0
and the power's lift): a region holding one distinct root (count 1, or a proved
multiplicity equal to its count) whose reflection in the real axis meets no
other region holds that root's conjugate as well, so the root is real
(`prove_real_by_reflection`). Such a region is reported real, as its real
segment: its imaginary ends become 0 and its real ends stay as published, so
the regions whose goal the solve judged only shrink. The proof runs in the
outermost solve only: a factor that the product lift solves alone keeps its
regions as they were, so the lift unites, merges and judges the factors' covers
as before, and the finish labels the cover the lift publishes. Any other region
meeting the axis is left unknown, and so is every region of a program with a
non-real constant, even where a root is real. A multiplicity above one is
claimed only where it is proved: at the point 0 (the order `v`, or `n` for a
monomial) and through a power's lift (`q`'s proved multiplicities times `k`). A
program's default goal is approximation. Under the accuracy goal every region
must hold one distinct root with its multiplicity proved, and success proves `p`
square-free when every count is 1. A real-only query, the isolation goal, a
region and the forced piecewise method are usage errors for a program.

## Candidate experiments and certification

The default complex path uses native candidate packets and certified subdivision,
with structural routes for exact power compositions,
useful exact centering, and sparse or small exact factors that need further work. Persistent interpolation covers
validate those routes against the original source. Small direct polynomial
candidates may evaluate coefficient midpoints; this never changes residue or
root certification. Internal experiments also provide bounded hardware-binary64 initialization, scaled
Gauss--Seidel Ehrlich--Aberth packets, full secular regeneration, exact variable
centering and a degree/observed-progress portfolio. Selection never depends on
fixture names or hashes. Binary64 initialization rejects range loss and limits
degree and sweeps; it only supplies starting points to higher-precision work.

Aberth uses `N/(1-N sum_{j!=i}1/(z_i-z_j))`, where `N=p/p'`. Every kernel that
takes this step (in binary64, scaled binary64, multiprecision or balls) is
listed in `src/aberth_step.hpp` with what it keeps: its `N`, its sums,
Gauss-Seidel or Jacobi order, its guards and stops; the operations they share
are there. Its scaling is an exact power-of-two variable change and coefficient
normalization. Midpoint arithmetic and Newton-polygon initialization are
heuristics. Finite, distinct candidates persist between packets and precision
changes. Adequately isolated roots can be frozen as candidate centers while
still participating in interactions and whole-vector validation. Tightened
accuracy clears that mask. No partial candidate vector is promoted to a complete
certified root result.

Full regeneration uses
`a_i=-p(b_i)/(c product_{j!=i}(b_i-b_j))`, for leading coefficient `c`, giving
`p/c=Q(1-sum_i a_i/(x-b_i))`. Balanced outward products reduce enclosure wrapping.
The removable Newton limit at `b_k` is
`a_k/(sum_{j!=k}(a_j+a_k)/(b_k-b_j)-1)`; ordinary off-pole evaluation includes the
logarithmic derivative of `Q`. Zero-containing node differences reject
regeneration. The weights and iteration remain proposals, validated against
the original coefficient balls by FLINT inclusion predicates. Residuals and
agreement with MPSolve do not certify counts.

Centering replaces a factor by a nonzero exact scalar multiple of `p(x+c)`.
Denominator clearing and integer/Gaussian-integer Taylor shift are bounded before
large powers are formed. Validation uses the exactly translated coefficients;
certified boxes are translated outward by `c` before ordinary coverage,
disjointness and accuracy checks. Independent exact evaluations check both
translation implementations.

The portfolio also has an analytic binomial branch for `a_n x^n+a_0`, with
provably zero intermediate coefficients and a nonzero enclosure of
`q=-a_0/a_n`. An outward nth root of `q`, multiplied by all nth roots of unity,
exhausts the `n` simple roots for every permitted coefficient choice. A fixed
rotation avoids a negative-real branch cut when possible; unresolved arithmetic
falls back. Only proved real branches have their imaginary interval set to zero.
The ordinary outward export and disjointness ledger must still pass before
coverage and accuracy are reported. Independent rational winding/axis-Sturm
checks cover the new branch, including complex and uncertain inputs.

For 128 or more candidate boxes, disjointness uses a sweep of outward real-axis
bounds. A pair is discarded only when those closed intervals are strictly
separated; all remaining pairs use the original complex overlap predicate.
Outward rounding can retain extra comparisons, never discard an overlapping
pair. Touching endpoints remain overlaps.

## Real roots of real input

A real-only query under Auto on exact real coefficients, at any goal but the
cluster radius, is answered by `RealOnlyRoute` as "Coverage is distinct from
accuracy" states; here is what each route proves. The exact real engine of
"Exact real isolation and refinement" (`RealRefinement`) is not on this
route: under public routing it answers only secular input's real-only queries.

- **The sparse recursion** (`SparseRealAttempt`) follows Sagraloff (2014)
  (the B6 note). With `x^i_0` the
  input's lowest term, 0 is a root of multiplicity `i_0`. On the positive
  half-line the chain `p_0=p/x^i_0`, `p_(j+1)=x^(1-e)*p_j'` (`e` the lowest
  positive exponent of `p_j`) loses a term a level down to a constant; a power
  of two `R` above Cauchy's bound for `p_0` bounds every level's roots
  (Gauss-Lucas). Bottom up, the breakpoints of level `j` are 0, `R` and the
  roots of `p_(j+1)` in brackets of proved opposite signs. The sign of `p_j`, a
  ball evaluation term by term (`SparseRealEvaluator`), is decided over each
  breakpoint's bracket, so no root of `p_j` lies there; between consecutive
  breakpoints `p_j'=x^(e-1)*p_(j+1)` has no zero, so `p_j` has one simple root
  if the two signs differ and none otherwise, and by induction none is missed.
  An undecided sign has its bracket refined and the precision raised, up to the
  cap, where it may mark a common root of `p_j` and `p_(j+1)` (as a repeated
  root of `p` is) and the recursion stops, as at its evaluation ceiling.
  Negative roots come from `p(-x)`. Each root is refined by the input
  (`real_refine`) and published as real, with count and multiplicity 1, once its
  box lies strictly inside its gap between breakpoints and meets the goal.
- **The axis charts** (`axis_only`; see "The piecewise method" and
  the axis-chart note) build only the charts centred
  on the axis. Each real model, with proved bounds on its errors for `p` and
  `p'`, splits its chart's segment into parts proved root-free and parts proved
  to hold one simple root; dominant terms prove the rest of the axis root-free,
  and each root is refined by the input until its box lies strictly inside its
  part. Repeated or power factors end the attempt.
- **The complex route** solves `RealOnlyRoute::complex_query` in a fresh
  `SessionState`: every root at the accuracy goal, with the seed pass below.
  Realness is proved by conjugation: in a complete disjoint cover of a real
  polynomial the conjugate of a region's one distinct root `z` lies in the
  region's mirror image, so in the region itself when no other region meets that
  image, and `z` is real (`prove_real_by_reflection`,
  `classify_cover_realness`). The route keeps these regions; the answer is
  complete when the solve met its goal, every kept root meets the caller's
  (`meets_goal`) and no other region meets the axis. Such a region, in an answer
  otherwise finished, makes the route solve afresh at twice the bits, until it
  is proved real or clear of the axis, the accuracy asked reaches the cap, or
  the work or sweep budget or a cancellation stops it.
- **A region** is applied to every answer by `RealRegionFilter`
  (the D1 note). A root's interval
  holding an end `e` is first tested exactly: for `e=r/s` in lowest terms the
  rational root theorem (`s` divides the integer form's leading coefficient, `r`
  its lowest nonzero one), then the exact value by Horner's rule along the
  support. If `p(e)!=0`, a root of odd multiplicity is placed by `p`'s signs and
  its interval cut (`place`, `cut`: "Coverage is distinct from accuracy"), the
  signs from balls up to the cap, then exact; one of even multiplicity, with no
  sign change, has the route solved afresh at twice the bits until its interval
  lies on one side. A counted cluster at a stop is dropped only when wholly
  outside. Under the isolation goal (`meets_goal`: one root, its multiplicity,
  any width) the sparse recursion publishes from a 1-bit goal.

### The all-real seed pass

`real_line_seeds` (its design)
proposes starting points on the real line for a polynomial whose roots are all
real. The internal flag `cps_options::real_line_seeds`, set by `complex_query`
and by the solver for Auto's queries on exact real coefficients, starts each
exact real factor of degree at least `heuristics::real_line_seed_min_degree`
from them, on the secular route and on piecewise's hand-off to it
(`seed_from_real_line`); the pass runs on the factor's integer polynomial
(centred where the start is), a root at zero divided out and seeded by 0.

- **Exact refusals.** If every root is real, Newton's inequalities
  `a_k^2*k*(n-k)>=a_(k-1)*a_(k+1)*(k+1)*(n-k+1)` hold for `0<k<n`; with
  `p(0)!=0` a zero coefficient lies between two of opposite signs, since the
  sign variations of `p(x)` and `p(-x)` must add up to `n` (Descartes' rule);
  and `(n-1)*p'^2>=n*p*p''` on the real line (Laguerre's inequality; off the
  roots, `n*sum_i 1/(x-x_i)^2>=(sum_i 1/(x-x_i))^2`, Cauchy-Schwarz).
  `newton_inequality_failure` checks the first two in integers, and each
  Laguerre iterate the third in balls (it is the step's discriminant). A certain
  failure proves a non-real root and ends the pass; passing proves nothing.
- **The interval.** Laguerre's method starts at `-B` and `B`, a power of two
  above `2*max_k |a_(n-k)/a_n|^(1/k)` (of Fujiwara's kind); for a polynomial
  with only real roots its iterates from outside move monotonically toward the
  extreme roots and never pass them. Its ends, once the steps fall below `2^-35`
  of the width, widened by `2^-30` of it, bound the interval.
- **Why an envelope.** A binary64 interpolant of `p` resolves values to about
  `2^-53` of its largest, but `log2|p|` can range over more than a thousand bits
  across the interval: it loses the sign changes where `|p|` is small
  (the earlier note). So the pass
  interpolates `f=p*2^-phi`: `2^-phi>0` keeps `p`'s sign changes and adds none,
  and `phi`, a smooth envelope of the local maxima of `log2|p|`, scales the
  peaks of `|f|` between roots to about one. Smooth, so that `f` keeps a short
  Chebyshev series: that of `2^-phi` has about `ln 2` times `phi`'s range in
  bits significant terms, which `f` needs as nodes beyond `n`, and a range of
  `log2|p|` over 4 bits a degree ends the pass.
- **The interpolant.** One change of basis and one DCT-I, in balls at precisions
  set by the cancellation and the range, give `p`'s values `v_j` at the nodes
  `cos(pi*j/N)`; fewer certain sign changes between nodes than half the degree
  end the pass. `phi` is a running median of `log2|v_j|`, smoothed by a Gaussian
  in the Chebyshev angle and raised until no local maximum lies above it; the
  `f_j`, rounded to binary64 once each `v_j` is known to `2^-54*2^phi_j`, give
  `f`'s coefficients by a binary64 DCT, doubling `N` to a small tail.
- **The seeds** are the interpolant's sign changes above its noise on a grid 16
  times finer, refined by the Illinois method, and its near-touching minima
  (more than `n` sign changes is noise and ends the pass), and the secant's zero
  in each cell between nodes with certain opposite signs and no seed. A seed
  alone in such a cell is kept, any other where `p`'s certain signs differ on
  its two sides (or, unproved, where a neighbour cut that test short); at most
  `n` are kept, the proved ones first.
- **Only starting points.** The pass proves nothing an answer uses: its certain
  sign changes bracket roots of odd multiplicity, but a seed only replaces a
  point of the circles' start, and the secular route certifies every region
  against the input. A wrong or missing seed costs sweeps, never an answer, and
  every failure of the pass (`RealLineSeedStop`) leaves no seeds.
- **Routing.** `every_root_seeded` (`real_line_routing.cpp`) runs the pass
  first, on the primitive part of the integer import, where the selector sends a
  query for all roots on exact real coefficients in dense form, at the accuracy
  or approximation goal, to piecewise (the real-only route's complex solve among
  them). If it seeds all `n` roots, `n` at least the degree floor, the secular
  route is taken from those seeds, which the factor equal to that polynomial
  reuses (`hint_real_line_seeds`). `n` seeds suggest, not prove, that every root
  is real; a misrouting costs time only.

## Exact real isolation and refinement

This exact real engine answers real-only queries on native secular input,
through the exact numerator. On coefficient input the routes of "Real roots of
real input" answer them; a public real-only query on inexact or non-real
coefficients is refused before any engine runs, so this engine answers no public
coefficient query (D5 routing, which only tests build, still reaches it).

Each exact squarefree real factor starts inside a strict global dyadic bound.
A persistent depth-first queue stores intervals and their exact primitive
integer affine transforms to [0,1]. Descartes sign variations after the positive
half-line transformation discard zero-root cells and isolate one-root cells.
Exact split-point zeros are separate rational components. A one-variation cell
is accepted only with nonzero endpoints, so boundary roots cannot be counted
twice. Cells with inconclusive variation counts remain in the queue.

Guarded modified-Newton proposals can contract a cell to one quarter of its
width. Both discarded sides must have zero Descartes variations and the new
endpoints must be nonroots. Failed proposals use exact dyadic subdivision.
These contractions are proposals validated by exact predicates; the implementation
does not claim the full complexity bounds of ANEWDSC. Queue steps, cancellations,
cell limits and precision changes retain already isolated roots.

Every retained irrational-root interval has one residual simple root and
non-root endpoints. A refinement step tries an interval Newton step first: its
image, with the derivative enclosed on the whole interval, is intersected with
its parent and accepted only when the represented result is contained in that
parent and strictly smaller. Then Abbott's quadratic interval refinement
(`refine_by_cells`): the secant through the ends places the root in one of
`N=2^g` equal cells (four at first), and the signs at that cell's two ends
confirm it, shrinking the interval `N`-fold; a confirmed cell squares `N` for
the next step and a refuted one halves `g`. Otherwise signs select a bisected
bracket. Every sign is certified (`sign_at`): balls at doubling precision, from
the precision that last decided a sign near the root, then the exact sign where
no ball decides, so a wrong guess costs time, never an answer; an exact zero at
a tested point becomes a rational component immediately. Bisection radii are
represented exactly, with a checked roundtrip, so a generic outward interval
constructor cannot enlarge the child beyond its parent. Precision increases
and rational-root re-enclosures consume later controller steps.

Factor multiplicities are retained. Publication checks outward exported
intervals for disjointness, true-relative accuracy and inter-factor overlap;
incomplete replacement work preserves a previous valid cover. Within a solve
the isolation keeps its intervals and local precision from one controller step
to the next. Individual exact arithmetic calls remain atomic.

## Source and representation contracts

`MonicSource` identifies one fixed monic polynomial and its known degree. Owned
joint evaluations enclose both p and p'; separate Newton and logarithmic-derivative
ratios fail when their denominator cannot be proved nonzero. Ratio failure is
not a root-count conclusion. Batch evaluation and optional Taylor/coefficient/
exact capabilities retain this richer interface. Public black-box inputs remain
deferred. Prepared handles own their polynomial and evaluator, surviving cache
eviction and source destruction; callbacks own their input coefficients.

For exact coefficients, cached translations construct a nonzero scalar multiple
of p(c+x) and normalize it to monic. Integer translation preflights every
coefficient of B^n p(x+A/B):
`sum(k>=j) a_k binom(k,j) A^(k-j) B^(n-k+j)`.
Checked bounds cover the full sum before large powers are materialized. Failed
cache preparation cannot expose a model from an earlier centre. Nearby queries
reuse a coarse dyadic centre and evaluate the complete shifted polynomial at
z-c; this is exact substitution, not an unbounded local approximation.

Automatic global centering considers dense exact real factors of degree at least
64. Its initial centroid is rounded to eight bits. The exact shifted monic
coefficient height must fall by about log2(degree) bits before the shift is used;
otherwise the original coordinates remain. This avoids damaging useful radial
scaling. When small exact real iteration stalls after a precision increase, the
controller can retain the certified cover and start a translated structured
source. Coefficient arithmetic quality is sampled separately from iteration
progress before increasing precision on large factors.

Exact `p(z)=q(z^k)` reduction uses the greatest common divisor of nonzero exponents
after squarefree preprocessing. A complete q-cover whose rectangles avoid zero
lifts through the k analytic root branches. Rotation removes a negative-real
branch cut; outward branch images carry the source count and multiplicity.
Balanced rectangle margins make subsequent contour checks possible. Any output
overlaps are merged before complete coverage is published. Individual roots
still require disjointness and accuracy in the original coordinate. Insufficient
lifting accuracy tightens the reduced goal and resumes retained state. The complex
solves of a real-only query on coefficient input use the reduction like any
complex query; the exact real engine (secular input) isolates on the original
polynomial.

## Resource and implementation limits

- The outer precision controller is bounded and records its last requested work
  precision. It does not claim to observe FLINT's internal peak precision.
- Exact gcd/squarefree work and individual exact arithmetic operations are atomic backend
  calls. They may allocate larger intermediates or take substantial time.
- The time limit and the cancellation token are checked between stages/attempts/factors and subdivision
  actions, not within those calls. A late cancellation can coexist with valid completed certificates.
- No hard time/RAM limit is promised. Use process isolation for strict job limits.
- Library-level allocation failure can terminate inside third-party C arithmetic;
  catching a C++ `bad_alloc` is not a universal recovery mechanism.
- Floating-point, compiler and third-party arithmetic correctness are trusted.
  The library is built without fast-math and does not alter global thread settings.
- Thread-local FLINT caches are released at explicit host-controlled boundaries
  (`ThreadScope` / `cps_cleanup_thread`). Objects made on the thread stay valid
  after it, which relies on FLINT's thread-aware memory manager: the library
  requires a FLINT with thread-local storage and thread support, and does not
  force global cache teardown.
- Caller work/cell ceilings can stop subdivision with an unmet goal. Bounded
  heuristic searches also have local stagnation stops; the complete policy uses
  its progress mechanism instead. A resumable frontier is not a portable proof object.
- There is no public evaluation-oracle interface, general native basis solver, GPU kernel,
  or formally verified checker yet.
