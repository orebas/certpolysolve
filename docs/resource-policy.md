# Resource policy and scaling limits

This is the current policy; the preceding audit
describes the old
state and is retained as historical evidence.

## Caller controls

Advanced controls: most callers set none of them. The controller work budget is a
soft limit (checked between units; see [the API guide](api.md#6-options)).

Eligible prepared exact-coefficient global Auto also has an internal allowance
for its numerical attempt and accompanying subdivision. Exhausting it starts
exclusive complete subdivision. The caller's remaining budget stays in force;
reaching this internal boundary does not itself report a caller resource stop.
The internal allowance, the scheduler's estimated arithmetic share and the public
controller work units are distinct quantities. None is a hard wall-time or
peak-memory bound. See
[the scheduling scope](production-methods.md) and
[the progress and cost account](mathematics.md#persistent-certified-subdivision-and-contraction).

The C API's `cps_resource_limits` separates three **per-stage estimates**:

| Setting | Meaning | Unspecified default |
|---|---|---|
| `exact_materialization_bits` | Admitted expanded coefficients, exact operand growth and applicable exact-preparation storage | `UINT64_MAX` |
| `exact_operations` | Exact preparation effort: imported-element admission for both modular square-free proofs (below the coefficient count Q(i) preprocessing is refused before any stage), degree-square estimates for bulk algebra/translation, explicit guarded arithmetic counts in Q(i), a conservative operation estimate for secular conversion | `UINT64_MAX` |
| `temporary_bytes` | Estimated engine/preparation/construction storage; also bounds applicable exact bit storage | `UINT64_MAX` |
| `flags` | `CPS_RESOURCE_ALLOW_EXPONENT_EXPANSION` explicitly allows expansion beyond the input-relative rule | `0` |

`UINT64_MAX` means no fixed policy ceiling. Zero is an explicit zero allowance.
These are not a global allocation counter, a bit-operation counter or an RSS/time
guarantee. Estimates restart at stage boundaries; a FLINT GCD/square-free call
is not preempted internally. Cancellation remains cooperative. External process
limits belong to the host/benchmark launcher.

Closed forms of degree two through four obey the same exact-effort policy.
Their admission uses the squared coefficient count, weighted by the word size
of their discriminant or invariant operands. Exact storage also admits the
quadratic discriminant and, before constructing them, each common-denominator
product and scaled coefficient for cubics and quartics. These conservative
operand bounds can refuse an expression whose cancellations would make the
final value smaller. A refused direct preparation
leaves the general numerical route available. Linear ball division needs no
exact algebra and remains available with zero exact-operation allowance.

```cpp
certpolysolve::Options options;
certpolysolve::ResourceLimits limits; // all allowances unspecified; flags = 0
limits.temporary_bytes = available_scratch_bytes; // optional caller policy
options.resources = limits;
```

C callers initialize with `cps_resource_limits_default`, then use
`cps_options_set_resource_limits` / `cps_options_get_resource_limits` (the `_v1`
names went in stage 3b-3). `cps_secular_create` takes these settings at
construction time (NULL: the defaults); the C++ `SecularEquation` uses the
defaults. Unknown flag bits are rejected. Existing ABI structures are not
resized.

An exact allowance alone is `exact_materialization_bits`, with the expansion
flag if wanted (the exact-budget setter of earlier versions went in 0.5). Zero
disables applicable exact preprocessing; a positive value with the flag admits
exponent expansion within that bit allowance.
Each solve reads its resource settings once, even with zero work; a solve
under other settings is another solve (the API has one-shot solves
only).

## Compact exponent rule

A short literal can denote an enormous expanded integer. Ordinary default input
`1@1000000000000` remains compact. Exact preprocessing compares a conservative
expanded-bit estimate against the existing owned coefficient representation,
including its stored significands/exponents and value-object storage. It writes
the coefficients out when the expanded estimate exceeds the owned representation
by at most one precision cap (`max_work_bits`, 65,536 bits by default) for each
nonzero number written out, unless the caller explicitly permits more. The
numerical solve holds each coefficient at no more than its precision cap, so
exact preparation may write each number out with as much again. A product of
roots written as m@e coefficients needs this room: the exponents a root adds
grow with the degree. Already supplied ordinary large integers do not encounter
a small fixed bit cap. A decimal exponent uses the conservative bound
4 × |exponent| bits; this is ceil(log₂10), not a benchmark-tuned parameter.

When the write-out is refused, the one-prime modular square-free test still runs
on the coefficients as written: it reduces each m@e number modulo the prime
directly, the exponent modulo p − 1. A pass keeps the polynomial's square-free
proof on the ball path; any other outcome says nothing, and multiplicities stay
unknown.

A single number that must be exact on its own (a region's end, a secular term)
is expanded only when the expansion is no larger than its owned representation.
Secular input is built before any solve, so no precision cap applies to it.

A refused polynomial transformation leaves numerical input available. Current
secular input construction requires exact normalized terms, so a refused
expansion returns `CPS_OUT_OF_RANGE` with a diagnostic. Proposed exact centers
receive the same materialization preflight before conversion from compact ball
midpoints into rational coordinates. This protects that conversion itself,
not just the subsequent Taylor shift.

## Disposition of the preceding audit

| Old restriction | Current disposition / remaining boundary |
|---|---|
| Public degree/index 1,000,000 | Removed; check backend index, container maximum and addressable storage before resize |
| Literal length 1,000,000 | Removed; valid grammar and string/backend representability remain |
| Secular normalization 8M bits | Removed; constructor uses explicit resource settings and compact-expansion policy |
| Requested/working precision 1,000,000 | Removed; backend signed precision representation remains, excluding FLINT's exact-precision sentinel |
| Exact default 1,048,576 / setter maximum 1B | Removed; unspecified default, uint64 caller allowance |
| Iteration maximum 1,000,000 | Removed; uint64 option, checked/clamped conversions at backend boundaries |
| Cell maximum 10,000,000 | Removed; positive uint64 caller allowance |
| Tolerance-comparison private 1M limits | Removed; current caller resource/working-precision limits and cancellation apply |
| Q(i) derived operation cap 1M | Removed; separate caller exact-operation setting |
| Shared engine 256 MiB / degree262144 | Removed; checked factory estimate versus caller temporary bytes |
| Structured companion degree349525 | Removed; checked storage estimate; optional method preconditions remain |
| Dense eigen n³ ≤ 2²⁹ | Retained as a named private cubic-work admission heuristic; temporary bytes are caller-controlled |
| JIT 16384 coefficients | Removed; checked addressable/caller preparation estimate |
| JIT retained cache64MiB | Retained as eviction policy; active modules keep ownership and ordinary evaluation remains possible |
| Piecewise construction256MiB | Removed; unspecified default, optional caller/private construction allowance determines batch capacity |
| Angular sectors2²⁸ | Removed; checked addressable sector/index representation remains |
| Piecewise coefficient log-magnitudes outside ±2⁵⁰ | Retained private geometry working range, with an explicit diagnostic; protects the current integer hull/geometry adapter, and can refuse this standalone route. It is an exponent restriction, not a degree or RAM limit. |
| FFT evaluator batches256MiB | Removed; checked caller storage estimate; ordinary evaluation fallback remains |
| Default 65536 working bits /100000 work units /65536 live cells | Retained as explicit caller-adjustable solve controls; none is an input degree ceiling |
| Controller iterations0 / work0 / exact0 | Automatic iterations / the automatic work budget (0 selects it since 0.5; no public setting means zero work) / disabled exact preprocessing |
| Benchmark256MiB,1GiB, time/precision/read-size limits | Explicit historical measurement recipes; not library defaults |

Overflow-sensitive precision growth, scaling products, exact bit accumulation,
allocation estimates and angular sizing use checked/saturating arithmetic.
Saturated storage estimates refuse admission rather than wrapping. The former
internal million-bit separation clamp masked a power-composition stopping defect:
the controller now stops when the caller's working precision is exhausted.

## Optional algorithm crossovers

These choices remain hypotheses about cost. A fallback's existence is not evidence
that its performance is competitive. This campaign does not indiscriminately
raise quadratic/cubic work thresholds.

There is no bound on the source's degree for the small cluster reconstruction: a former one left close pairs
above it to a slow walk and protected nothing measured (the comment in `src/heuristics.hpp`).

| Private name / default | Role and fallback |
|---|---|
| `local_component_degree`64 | Small cluster reconstruction; otherwise component-scale restart |
| `secular_reconstruction_degree`2048 | Selected automatic reconstruction choices; coefficient/direct or explicit fixed-method routes remain |
| `secular_block_seed_degree`2048 | Older block initializer; ordinary secular initialization |
| `binary64_seed_degree`2048 | Legacy double seed helper (the DoubleSeed test engine's, research build only); does **not** cap the modern machine engine |
| `cluster_recovery_degree`2048 | Bounded candidate-guided cluster recovery; retain count cover/subdivision |
| `annulus_estimate_degree`2048 | Optional three-center annulus probe; ordinary box counting |

Each listed constant has a private `CERTPOLYSOLVE_EXPERIMENT_<UPPERCASE_NAME>` macro
in `src/heuristics.hpp`. For a separate research build, pass the definition through
`CMAKE_CXX_FLAGS` so **all translation units agree**. Do not define it for only one
target. The dense work guard uses `CERTPOLYSOLVE_EXPERIMENT_DENSE_WORK_LIMIT`.
For example, `-DCERTPOLYSOLVE_EXPERIMENT_LOCAL_COMPONENT_DEGREE=32` changes that
one crossover; it grants no proof or performance guarantee. No override is enabled
in this campaign's qualified build.

Private factor traces emit `algorithm_degree_policy` eligibility records containing
the effective degree, cutoff and fallback description. They explicitly describe
eligibility, not a claim that every optional helper was attempted. Component
restart traces separately show actual active component sizes. Record compiler
flags and source/binary hashes with experimental observations.

Other retained cost/geometry choices: interaction-tree crossover128 and order20;
small-complex64; centering scheduling256/64 with8/32/1024-bit center heuristics;
Mandelbrot recognition65535 (expanded input still supported, recurrence deferred);
optional Taylor/annulus exponent probes±1M and machine proposal scaling±10⁹;
Pellet≤16 Graeffe steps per attempt; Cauchy sample/exclusion allowances; cluster
recovery≤4096 tests/16 rounds; annulus≤40 sampled exponents; local initializer
retry/iteration packets; Taylor cache8 centers and coefficient cache64 × (n+1).
These bound optional attempts or evict caches, then return to the controller.
Machine representable ranges, nonzero denominators, rigorous error inequalities,
integer dimensions and proof-count consistency are mathematical/safety checks
and remain mandatory.

The PSCASE readers'10001/2²⁰ rows, seed readers'10000/2²⁰ points, per-worker65536
working bits and86400-second argument maximum remain benchmark-format/protocol
limits. Historical recipes in the current piecewise/classical workers explicitly
spell out legacy exact settings; `resources-v1` selects current defaults. Frozen
old binaries are the authority for reproducing historical timings exactly.
The current machine's32GiB/1GiB aggregate guard and campaign deadlines are launcher
policy, not limits a larger-machine library user inherits.
