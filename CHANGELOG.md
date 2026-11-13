# Changelog

## 0.5.1 (12 November 2026)

- Renamed from PolySolve to CertPolySolve: namespace `certpolysolve`, headers `<certpolysolve/certpolysolve.hpp>`
  and `<certpolysolve/certpolysolve.h>`, CMake package `CertPolySolve` with target `CertPolySolve::certpolysolve`, C
  prefixes `cps_` and `CPS_` (were `ps_` and `PS_`). The entries below 0.5.0 keep the names they shipped under.
- Licence: GNU LGPL, version 3 or later (was MIT).
- `CPS_ROOT_ISOLATION` (`Goal::isolation()`) under Auto also accepts exact complex coefficient input. It returns
  each distinct root in a disjoint closed rectangle with its multiplicity, without requesting a width.
- Auto's global route for exact coefficient input starts subdivision alongside numerical refinement after a
  precision increase leaves a certified cluster unresolved. A shared internal work and precision allowance
  eventually gives subdivision exclusive service. Certified regions are retained, and caller limits still apply.

## 0.5.0 (29 October 2026)

This entry covers the changes since the transfer checkpoint of 8 September 2026, the tag `transfer-2026-09-08`
(docs/transfer-release.md). That version called itself 0.1.0-dev, with ABI 1. It took
one kind of input, a polynomial given by its coefficients. It had two goals and a cancellation callback. Its
real-root isolation and its local models used FLINT's root solvers.

Code written for the checkpoint needs changes. Before 1.0.0 no release promises compatibility of any kind.
[docs/api.md](docs/api.md) is the guide to the API as built, and its section 13 maps each removed call to its
replacement. The specification's section 10 lists the changes item by item (C1 to C40).

### The API

#### C

- **A specification.** [docs/api-spec-0.5.md](docs/api-spec-0.5.md) specifies every C and C++ entry point: inputs,
  options, goals, results, errors, lifetimes and thread guarantees. A note at its top lists where the built library
  departs from it on purpose. There the headers and docs/api.md describe the library.
- **Errors and outcomes.** A function that can fail returns a `ps_error`. `PS_OK` from a solve means that the call
  was valid and ran. The result says whether the goal was met. A precision stop, a spent budget, a time limit, a
  cancellation and a refused input are outcomes: each returns `PS_OK` and a result that says what happened.
- **Usage errors.** A combination of options that no input of the entry point's kind can answer now fails at the
  solve with `PS_INVALID_ARGUMENT`, and no result is made. Before 0.5 it returned `PS_OK` with
  `PS_UNSUPPORTED_QUERY`. A refusal that depends on the input stays an outcome. Section 7 of docs/api.md has one
  table per input kind.
- **Failed calls change nothing.** Setters and builders are transactional. `ps_last_error()` describes the calling
  thread's last failed call, and a successful call clears it.
- **NULL.** A NULL imaginary part or radius means 0 in every input function. Section 2 of docs/api.md lists every
  NULL rule.
- **Options.** Every setting but the cancellation token has a getter, and `ps_options_copy` copies a whole set. A
  setter checks only its own arguments. The solve checks how the settings combine.
  - `ps_options_set_precision(o, bits)` takes the requested accuracy alone. The precision cap has its own setter,
    `ps_options_set_max_work_bits`.
  - `ps_options_set_method` picks `PS_METHOD_AUTO` (the default), `PS_METHOD_CLASSICAL`, `PS_METHOD_SECULAR` or
    `PS_METHOD_PIECEWISE`.
  - `ps_options_set_goal` selects a goal (see "Queries and goals").
  - `ps_options_set_time_limit` bounds the wall time of a solve, in seconds on a monotonic clock. It covers everything
    the call does.
  - A cancellation token replaces the callback: `ps_cancel_token_create`, `ps_cancel_token_request`,
    `ps_cancel_token_requested`, `ps_cancel_token_destroy` and `ps_options_set_cancel_token`. Any thread may request
    a token at any time. The solve only reads a flag, so no caller code runs inside it.
  - `ps_options_set_work_budget` and `ps_options_set_live_cells` set the search's budgets, and 0 selects the default.
    The automatic work budget is 100,000 units or 64 units a root, whichever is more.
  - `ps_options_set_max_iterations` is now a budget of sweeps for the whole call, with its own stop reason. A sweep
    is one pass of the simultaneous iteration over all its root approximations.
  - Resource limits replace `ps_options_set_exact_budget`: `ps_resource_limits`, `ps_resource_limits_default`,
    `ps_options_set_resource_limits` and `ps_options_get_resource_limits`, with the flag
    `PS_RESOURCE_ALLOW_EXPONENT_EXPANSION`. By default no allowance has a fixed ceiling. The checkpoint's exact budget
    defaulted to 1,048,576 bits.
  - `ps_options_set_region` limits a real-only query to a closed region.
- **Results.**
  - `ps_result_stop_reason` says why a solve stopped: a work, cell or iteration budget, the precision cap, the input,
    a cancellation or the time limit. `ps_termination` gained `PS_WORK_LIMIT`.
  - Facts that travelled only in the diagnostics text are now fields: `ps_result_sufficient_work_bits`,
    `ps_result_squarefree` and `ps_result_handoff`. `ps_result_preparation_limited` says when the resource limits
    limited exact preparation, so that multiplicities may be unknown.
  - Exact root ends can be used without a big-number library: `ps_dyadic_compare`, `ps_dyadic_to_double` (rounded
    down, to nearest or up), `ps_dyadic_format`, `ps_dyadic_rational`, `ps_dyadic_midpoint` and `ps_dyadic_width`.
  - `ps_statistics` and `ps_engine_counters` replace `ps_statistics_v1`. `ps_statistics` is a promise, and its size
    field lets it grow. The engine counters are diagnostics and may change in any release. `ps_result_component_count`
    and `ps_result_get_component` describe the search's state at the end.
  - `ps_termination_name` and `ps_stop_reason_name` give English names.
  - The regions come in increasing order of their exact midpoint, real part first, on every route.
  - Nine invariants hold for every result (docs/api.md section 8). The termination is `PS_SUCCESS` exactly when the
    goal is met. A met goal implies a complete answer. The regions are pairwise disjoint, and each carries a proved
    count. The test suite checks them, and the fuzz harness checks them on every result.
  - A solve whose goal is met reports `PS_SUCCESS`, even when a cancellation arrives afterwards.
- **Versions.** `PS_VERSION_MAJOR`, `PS_VERSION_MINOR`, `PS_VERSION_PATCH` and `PS_VERSION_STRING` give the headers'
  version, and `ps_version()` the library's. `PS_ABI_VERSION` and `ps_abi_version()` are 2. The library does not
  check them. A program may compare them to detect a mismatched library.
- **Memory.** No exception crosses the C boundary. An allocation failure inside the library gives `PS_OUT_OF_MEMORY`
  and no result. One inside FLINT or GMP still ends the process.

#### C++

- **A new header.** `<polysolve/polysolve.hpp>` replaces the earlier wrapper. It needs C++20, lives in namespace
  `polysolve` and is header-only on top of the C API.
- **Values.** The types are values that copy. Numbers, coefficients, polynomials, exact ends, intervals, roots and
  stops print with `std::format` and `std::println`. Exact ends compare exactly, also with integers and doubles.
- **A list either way round, by name.** `Polynomial{-2, 0, 1}` reads a list with the constant term first, and says
  so; `Polynomial::lowest_first` is the same by name, and `Polynomial::highest_first({1, 0, -2})` reads the highest
  power first, as MATLAB's `roots` and `numpy.roots` do. Printing the polynomial shows which one was built.
  `Result::values()` gives the roots as plain `std::complex<double>` numbers, each region's value repeated as many
  times as the region holds roots. After a stop or an error it is empty.
- **Numbers.** A `Number` is one exact number: an integer, a double at its exact binary value, text such as `"1/3"`,
  `"0.1"` or `"3@-2"`, or a complex number. A character or a `bool` is refused at compile time.
- **Inputs.** A `Polynomial` is built from a list (`Polynomial{-2, 0, 1}`), from its terms (`Polynomial::sparse`),
  from its roots (`Polynomial::from_roots`) or by arithmetic on `Polynomial::x()`. It is read back with `degree()`,
  `p[e]`, `coefficients()` and `terms()`. `ball(value, radius)` makes a coefficient known within a radius.
  `SecularEquation` and `Program` are the other inputs.
- **Options.** `Options` is a struct of named fields: `bits`, `max_bits`, `time_limit` (any `std::chrono` duration),
  `goal`, `method`, `stop` (a `std::stop_token`), `work_budget`, `cell_budget`, `iteration_budget`,
  `absolute_tolerance` and `resources`. Designated initializers set them in that order, as in
  `{.bits = 128, .time_limit = 10s}`.
- **Solving.** `solve(input, options)` gives every complex root. `real_roots(input, options)` and
  `real_roots(polynomial, region, options)` give the real ones, and `between`, `at_least` and `at_most` make a
  region.
- **Results.** A `Result` has `ok()`, `stop()` (a `Stop` value that prints its name), `message()`, `complete()`,
  `degree()`, `roots()` and the facts of the C result. A `Root` has its exact `box`, its `count`, its `multiplicity`
  (empty when unknown), its `realness`, its `certificate`, its `centre()` and its `value`. The value is the number
  the C call `ps_result_get_approximation` gives, bit for bit, and it is empty where that call refuses: an unmet
  goal, an unresolved cluster or an isolation answer.
- **No exceptions.** Nothing throws except `std::bad_alloc`. An input remembers its first bad value (`ok()`,
  `status()`). A solve of it returns a `Result` whose `stop()` is `Stop::error`.
- **Threads.** A `ThreadScope` at the start of a thread frees the thread's FLINT caches when it ends.
- **GMP.** The optional `<polysolve/gmp.hpp>` lets GMP's `mpz_class` and `mpq_class` serve as numbers. Its `to_mpq`
  converts a root end to `mpq_class`.

#### Command line

- `polysolve_cli` gained `--method auto|classical|secular|piecewise`, `--approximation` and `--secular FILE`. The
  last reads a secular equation in MPSolve's file format (its modern subset).
- `polysolve_cli` reads its coefficients with the constant term first, or with the highest power first after
  `--highest-first`, and echoes the polynomial it read on stderr (`p(x) = x^2 - 2`), so a reversed list shows.
- `polysolve_cli --clusters R` asks for boxes of radius at most R, each with its count (the cluster-radius goal).

### Inputs

- **Coefficients.** The number grammar and the coefficients with error bounds (balls) are as before. New:
  - A polynomial keeps only its nonzero terms. Any exponent up to 2^63 − 1 can be set, so x^1000000000 + 1 is two
    terms. `ps_polynomial_terms` lists the exponents. Every route except the real-only sparse recursion still reads a
    full coefficient list, so a solve needs room for degree + 1 coefficients.
  - `ps_polynomial_degree`, `ps_polynomial_get_coefficient` and `ps_polynomial_copy` read and copy a polynomial. A
    coefficient comes back as exact text that sets the same value again.
  - Exact polynomial arithmetic over the Gaussian rationals (complex numbers with rational parts):
    `ps_polynomial_add`, `ps_polynomial_subtract`, `ps_polynomial_multiply`, `ps_polynomial_scale`,
    `ps_polynomial_power`, `ps_polynomial_derivative` and `ps_polynomial_from_roots`, which takes multiplicities and
    complex roots. In C++ they are `+`, `-`, `*`, `pow`, `derivative()` and `Polynomial::from_roots`. Arithmetic on
    few terms stays few terms. Expanding a product or a recurrence gives up the structure that a program keeps.
- **Native secular equations** (new). `ps_secular_create` and `ps_solve_secular`, in C++ `SecularEquation`, take
  S(x) = Σ aⱼ/(x − bⱼ) − 1 with exact residues aⱼ and distinct poles bⱼ. The eigenvalues of a diagonal matrix plus a
  rank-one change are the roots of such an equation. A query for the complex roots iterates and certifies on the
  terms. Only a real-only query, or a complex query that the terms do not finish, forms the exact numerator
  polynomial.
- **Programs** (new). A program, or straight-line program, is a list of steps that computes the polynomial: x, exact
  constants, sums, differences, products, negation and powers. It suits polynomials that are easy to compute and
  enormous to write out, such as iterations, recurrences and products.
  - A program is built with the C builders (`ps_program_create`, `ps_program_variable`, `ps_program_constant`,
    `ps_program_add`, `ps_program_subtract`, `ps_program_multiply`, `ps_program_negate`, `ps_program_power`,
    `ps_program_set_output`), read from its text form (`PSPROGRAM 1`, by `ps_program_parse`), or recorded from a C++
    function by `Program::from`.
  - `ps_solve_program` never expands the program. It proves the degree and the leading coefficient from the program
    and bounds the roots. Then it iterates and certifies on evaluations of the program.
  - A power in the program certifies a multiplicity, and it is solved through its base. A product is solved factor
    by factor. Squaring recurrences and compositions give the iteration its starting points from the program's
    structure.
  - When every constant is real, a region that holds one root, and whose mirror image in the real axis meets no other
    region, is reported real.
  - A program answers queries for all complex roots. Its default goal is the approximation goal.

### Queries and goals

- **Real roots only**, on exact real coefficients (`ps_options_set_real_only`, C++ `real_roots`). Every region of the
  answer is proved real, and no real root is missed. Under Auto the first of three routes that answers gives the
  result:
  - On inputs with few terms, a sparse recursion whose cost grows with the number of terms.
  - Where Auto would choose the piecewise method, that method's charts on the real axis alone.
  - The complex route, whose regions proved real are the answer.

  Where the complex route stops at the precision cap, the query stops too. A secular equation's real roots come from
  its exact numerator polynomial, by exact real isolation. Coefficients with error bounds and complex coefficients
  are refused.
- **Regions.** On coefficient input, `ps_options_set_region` (C++ `between(a, b)`, `at_least(a)`, `at_most(b)`)
  limits a real-only query to the real roots in a closed interval. A root exactly at an end is decided by an exact
  test, and it is in the answer.
- **Goals.**
  - `PS_ROOT_ACCURACY` (C++ `Goal::accuracy()`): each region holds one distinct root, with its multiplicity, to the
    requested relative accuracy. It is the default for coefficients and secular equations.
  - `PS_ROOT_APPROXIMATION` (`Goal::approximation()`), new: a region may hold several roots, counted, when its
    midpoint approximates each of them to the requested accuracy. It is the default for programs.
  - `PS_ROOT_ISOLATION` (`Goal::isolation()`), new, for real-only queries on coefficient input: each real root alone
    in a closed interval, with its multiplicity. No width is promised.
  - `PS_CLUSTER_RADIUS` (`Goal::clusters(r)`): each region lies within distance r of its midpoint.
  - `PS_GOAL_DEFAULT`, new: the options' initial goal, resolved per input kind at the solve. `ps_result_goal` reports
    the goal the solve applied.
- **The precision cap.** The default cap on the working precision is 65,536 bits (`PS_DEFAULT_MAX_WORK_BITS`). It was
  4,096. A solve on exact real input that stops at the cap names a cap that suffices, from Mahler's root-separation
  bound (`ps_result_sufficient_work_bits`).
- **Methods.** Auto chooses a route by a fixed rule that reads only the input and the query. A forced method names
  the principal numerical method. Forced piecewise hands a solve to the secular method within the same call when it
  stops for a reason that is not the caller's, and `ps_result_handoff` reports it. Section 7 of docs/api.md lists
  which method answers which query.

### Reliability and correctness

- **All root solving is PolySolve's own code.** FLINT supplies arithmetic, exact polynomial algebra and low-level
  inclusion tests. No route calls a FLINT root solver, and a test checks that none is linked.
- **Every stop keeps its proof.** A solve stopped by the precision cap, a budget, the time limit or a cancellation
  returns certified regions and the reason. Every region carries a proved count, whatever the termination.
- **Completeness.** The secular iteration has no convergence proof. A subdivision search, which has one, runs beside
  it on a share of the work. So every factor finishes or stops honestly within its budgets.
- **Deterministic budgets.** The piecewise method decides only by counted work. The work and iteration budgets are
  deterministic. An answer under a time limit depends on the machine and its load.
- **Fixed:** a solve blocked by a budget could report the precision cap far below the cap. It now reports the budget.
- **Fixed:** a power factor whose region overlapped another factor's could stop at a false precision cap.
- **Fixed:** from degree 40 on, the accuracy goal failed on exact inputs with multiple roots, because exact
  preparation refused to write their large coefficients out. It now writes them in a compact form within an
  allowance ([docs/resource-policy.md](docs/resource-policy.md), "Compact exponent rule").
- **Testing.**
  - A seeded property fuzzer solves small polynomials with planted roots through the public API, under tight limits.
    It checks every returned statement exactly: counts, disjointness, labels, the accuracy contract and stop reasons.
    Its families add real roots, regions and the isolation goal, native secular input, coefficients with error bounds
    and programs.
  - An independent checker, built on FLINT alone, audits real-root answers.
  - Contract tests cover the invariants, the tables of what each input kind and method supports, the NULL rules, the
    option getters and the C API called from a C program.
  - A `sanitized` preset builds the suite with AddressSanitizer and UndefinedBehaviorSanitizer. A `tsan` preset runs
    the thread tests under ThreadSanitizer.

### Performance

The changes with the largest effect, by the final report of the complex-roots phase
(docs/final-report-20261011/README.md):

- **Auto chooses its route by an estimate.** It asks the piecewise method for an estimate of its own work. Then it
  sends the input to the secular route (an iteration in the style of MPSolve) or to the piecewise route (local Taylor
  models in the style of PWPoly). A piecewise solve whose estimate proves wrong moves to the secular route.
- **A fast path first.** The secular route begins with a certificate attempt: one run of the Aberth iteration in
  double precision, then a certificate at every point. The piecewise route makes the attempt first when it is cheap
  against the route's estimate. A well-conditioned input within the range of doubles finishes there.
- **Clustered roots.** The secular route restarts a group of approximations that crawls onto a cluster. The restart
  starts from a local model at the cluster.
- **Large dense inputs.** The piecewise method builds its local models more cheaply.

On 11 October 2026 (main `3a9905f`), Auto finished all 769 inputs of the comparison against nine other solvers. On
the 106 inputs where either side needed a second or more, it was faster than or level with the best other solver on
every one. Its time was 0.12 of the best other solver's by geometric mean. The real-root, program and API work of
11 to 24 October came after this measurement.

Since that report:

- **Real roots.** The three real-only routes above. On exact real input, a factor of degree 64 or more whose roots
  are all real starts its iteration from seeds that the all-real seed pass places along the real axis. This holds
  for real-only queries and for queries for all roots.
- **Sparse input.** A real-only query on a polynomial with few terms works on the terms and needs no coefficient list.
- **Programs.** The degree, the root bound and the starting points come from the program. Roots are certified from
  evaluations of the program, so the coefficients are never formed.
- **Inputs with the root 0.** When the library solves a coefficient input for all its roots, the exact zero roots
  are set apart first. The root 0 is published as an exact point with its multiplicity, and the other roots are
  found from the polynomial divided by that power of x, by the route such a polynomial takes on its own. A modular
  test then tries to prove square-freeness, and FLINT factors the input only when the test does not prove it.

**Easy inputs.** From 25 to 29 October 2026 the work went to the everyday case: degree under 100, ordinary
coefficients, sometimes multiple or clustered roots. On the easy panel's 850 inputs of degree 1 to 128, warm solves
became about 3.5 times as fast from 25 to 27 October
(dev/docs/easy-panel-rerun-20261027.md), and about a tenth faster again
from 28 to 29 October. There the fastest certified competitor (MPSolve) takes 3.6 times our time by geometric mean.
It is faster in 8 of the 175 cells, all of them perturbed clusters of degree 8
(dev/docs/comparisons-20261029.md).

- **A certificate in binary64 first.** The certificate attempt first tries each root with Rouché's theorem from one
  compensated evaluation in binary64. Error-free transformations bound every rounding, and the test costs about a
  hundredth of the multiprecision certificate. Where the coefficients' own errors already defeat the test, it
  declines before it evaluates. Where the test fails close to the root, it evaluates once more at the Newton step: a
  degree-3,200 input went from 17.5 to 6.2 billion instructions. At a coarse goal the goal's own circle is too wide
  for the test, so it tries a smaller circle from the same evaluation: random dense inputs of degree 2,048 to 4,096 at
  10 bits now take about a third of their former time.
- **Closed forms for degrees 1 to 4.** Exact input of degree 1 or 2 is solved directly. Degrees 3 and 4 use
  Cardano's and Viète's forms and Descartes' factors in ball arithmetic, and exact invariants decide the multiple and
  the real roots. The goals and queries they decline go to the general route.
- **Small clusters.** Up to 16 points that the binary64 certificate refuses take Aberth steps together in
  double-double, and the certificate tries them again. A group of refused points whose Newton disks meet restarts
  from seeds about the cluster's centre, at the precision the cluster needs.
- **Precision raised to the need.** A certificate refused for rounding raises its precision to what its evaluation
  shows is missing. It used to double it. laguerre-4000 at 6,400 bits went from 190 to 120 s.
- **Less bookkeeping.** The disjointness test sweeps. Region ends stay binary until their text is written, once. The
  selector counts its rings in binary64 where that provably agrees with 96 bits. The certificate's binomials fit a
  machine word.
- **The all-real seed pass** keeps the result of a pass that stops, for a later hand-off to the secular route. It
  decides Newton's inequality from the coefficients' leading bits.

### Thread safety

- Any number of threads may call the library at once, each on objects it does not share.
- An input, an options object or a result may be read by any number of threads at once, while no thread changes or
  destroys it. Solving with an input or with options is reading them.
- An object may be made on one thread and used or destroyed on another, also after the making thread has called
  `ps_cleanup_thread`. At the checkpoint a thread had to destroy its objects before its cleanup.
- `ps_cancel_token_request` may be called on any thread at any time.
- No caller code runs inside a solve. The library keeps no mutable global state apart from each thread's error
  message. It reads no environment variable and writes nothing to `stdout` or `stderr`. It starts a thread only
  inside a solve whose options ask for more than one worker, and it joins every such thread before the solve
  returns.
- A solve saves the calling thread's floating-point environment. It runs with round-to-nearest and every exception
  masked, and it restores the caller's environment, raised flags included, before it returns.
- The diagnostics text and `ps_dyadic_format`'s output do not depend on the C or C++ locale: the decimal point is
  always `.`.
- A thread that calls the library needs at least 512 KiB of stack.
- Moving objects between threads relies on FLINT's thread-aware memory manager. So the build requires a FLINT built
  with thread-local storage and thread support (`FLINT_USES_TLS`, `FLINT_USES_PTHREAD`).

### Multithreading

A solve can split some of its loops across worker threads. `ps_options_set_workers` (C++ `Options::workers`) sets
the count. The default is one worker: then no thread starts and every loop runs as before. This is a first cut.
Its design: docs/multithreading-design-20261023.md.

- The same regions, counters and diagnostics at every worker count. A loop is split into blocks fixed by the
  data, never by the count.
- The library's own workers, threads of the C++ standard library, started and joined inside each solve. Nothing
  outlives the call. A block never waits and never locks. FLINT's thread pool is never used.
- Five loops split: the fixed-point transform, a Taylor row, the sweep of the multiprecision secular packet, the
  finisher's blocked evaluation and the binary64 engine's evaluations of a sweep. Every other loop runs on the
  calling thread.
- On eight workers the largest inputs measured (degree 16,384 to 24,576) solve 1.28 to 1.39 times as fast as on one,
  and smaller ones 1.0 to 1.26 times. The split loops scale almost linearly. The rest of a solve, about two thirds
  on the largest inputs, runs on the calling thread.

### Platforms and builds

- **Linux.** PolySolve is developed on Linux x86-64 and tested with Clang 21 and GCC 15.
- **Windows.** The MSVC portability work is on main: 128-bit integers and checked arithmetic for compilers without
  `__int128` and GCC's overflow built-ins, and 64-bit exponents where `long` has 32 bits. An MSVC build of
  21 October passed every test.
- **macOS** is not tested yet.
- **The same route on every platform.** Auto's routing decisions compare exact integer counts and double-precision
  ratios. They use no `long double`, and their base-2 logarithm and exponential are the library's own. So an input
  takes the same route on every platform. The numerical stages call the platform's math library, so the last bits of
  an answer can differ between platforms. The contract is the same on all of them.
- **The same bits on every CPU.** The library is built with `-ffp-contract=off`, so a product and a sum are fused
  only where the source writes `std::fma`. On x86-64 CPUs with FMA and AVX2, one binary gives the same bits. On
  other CPUs glibc's math functions can round differently, and a solve can take another iteration path. Its answer
  keeps the contract.
- **Presets.** [CMakePresets.json](CMakePresets.json) has named configurations (`release`, `dev`, `sanitized`, `tsan`,
  `research`) and test presets (`quick`, `correctness`, `sanitized`, `tsan`). Each builds in its own folder under
  `build/`.
- **Options.** `POLYSOLVE_PORTABLE_INT128` builds the portable 128-bit integers on GCC and Clang too. The research
  code, the stress campaigns and the benchmarks stay off by default.
- **Install.** The installed package now includes `<polysolve/gmp.hpp>` and the third-party notices of `licenses/`.

### Documentation

- [docs/api.md](docs/api.md) is rewritten for 0.5. It covers every input, option, entry point and result, with the C
  name beside each C++ name.
- [docs/quickstart-cpp.md](docs/quickstart-cpp.md) teaches the C++ API through five short programs, one for each
  common task, and [docs/examples-cpp.md](docs/examples-cpp.md) has eight more by topic. They are in `examples/`,
  and the test suite builds and runs each one.
  C and in C++, is in `examples/guide/`, and the suite builds and runs it too.
- [docs/mathematics.md](docs/mathematics.md) states what each certificate proves. It gained sections on Auto's
  routing, the certificate attempt, the piecewise method, program input and real roots.
- Every `src/*.cpp` file opens with a header: what it does, its stages in order, its state, its entry points and
  their callers.

### Removed

- FLINT's root solvers, which the checkpoint used for real-root isolation and local models.
- `ps_cancel_callback` and `ps_options_set_cancel`. The time limit and the cancellation token replace them.
- `ps_options_set_exact_budget`. The resource limits replace it.
- `ps_statistics_v1` and `ps_result_get_statistics_v1`. `ps_statistics` and `ps_engine_counters` replace them.
- The three-argument `ps_options_set_precision(o, bits, cap)`. It is now `ps_options_set_precision(o, bits)` with
  `ps_options_set_max_work_bits(o, cap)`.
- The earlier C++ wrapper: `check()` and its exceptions, `RootResult`, `RootEnclosure`, `roots()`, `clusters()`,
  `enclosures()`, `approximations()`, `set_ball` and the chained `Options` setters.
- These came after the checkpoint and went before 0.5:
  - The public sessions (`ps_session_*`, C++ `SolveSession`), which resumed a solve with a new budget, goal or cap.
    Every solve is now one call. A different budget, goal or cap is a new solve.
  - Interim names: the `_v1`, `_v2` and `_v3` structures and calls, `ps_options_set_search_budget`,
    `ps_secular_create_with_resources_v1` and `PS_COMPONENT_STAGNATION`. Section 13 of docs/api.md gives each
    replacement.
  - The recognition of the Mandelbrot recurrence. No route keys on a particular family of polynomials.
  - Support for BEFFT, an optional FFT library that never shipped. Its adapter, its test and its vcpkg port are kept
    in research/rejected/befft/, with a patch that restores them.
  - Research engines in the default library. The Sturm and portfolio engines, the exact residues, the machine
    engine's research methods (Jenkins–Traub and Laguerre among them) and the test candidate engines now build only
    with `POLYSOLVE_BUILD_RESEARCH`. No public method used them.
