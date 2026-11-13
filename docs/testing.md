# Testing

How to run CertPolySolve's test suite and what it checks.

## Build and run

```sh
cmake -S . -B build/release -DCMAKE_BUILD_TYPE=Release -DFLINT_ROOT=/path/to/flint
cmake --build build/release --parallel 4
ctest --test-dir build/release --output-on-failure
```

The tests need Boost's headers (Boost.Multiprecision is the exact reference), and a few need Python 3. Tests carry
labels, set in `CMakeLists.txt`: `ctest -L quick` runs the short cycle, `-L correctness` every mathematical test,
`-L api` the tests of the public API, `-L threads` the tests that solve on several threads at once. Use `ctest -N`
to list the tests available in your build. `-DCERTPOLYSOLVE_BUILD_STRESS=ON` also builds the stress campaign.

## What the suite checks

Every solve is one call, so a test that asks for a different budget, goal or cap makes a fresh solve.

- **The public API through its own headers:** `c_abi` (the C API from a C program, the token and the time limit
  included), `cpp_api` (the C++ header: numbers, polynomials built and read every way, printing, exact root ends,
  every solve and every stop without an exception, programs and secular equations, the two list forms and
  `values()`), `cpp_gmp`, `exact_text` (the exact end operations, the stop names, reading a polynomial back),
  `polynomial_arithmetic` (every arithmetic call against exact Boost rationals, and its errors), `quickstart_*` and
  `topic_*` (each program of [the quick start](quickstart-cpp.md) and of [the examples](examples-cpp.md)),
  `readme_first_solve`, `guide_first_solve_c` and `guide_first_solve_cpp` (the first solves of the README and of
  [the guide](api.md)), and `examples_quoted` (the pages quote their programs as they are).
- **The contract:** `api_additions` (the statistics, engine counters, components, stop reason, preparation flag and
  facts against the engines' own counters; the size protocol; every option getter, the copy and the setters; the
  resource limits; the testing poll counter of `certpolysolve_testing.h`), `result_invariants` (the outcome invariants
  and the regions' midpoint order on solves stopped at chosen checkpoints by that counter, for every method, a
  real-only query, a program and a secular equation), `cancellation` (the time limit, the token, the version macros,
  `cps_last_error` cleared by a successful call), `iteration_budget`, `floating_point_environment` (a solve runs in
  the default floating-point environment and gives the caller's back), `allocation_faults` (an allocation failure
  gives `CPS_OUT_OF_MEMORY` and no result) and `threads`.
- **Inputs and routes:** `program` (program input against exact expansions and coefficient twins),
  `borrowed_terms` (a solve on the borrowed terms answers as one on the coefficient list), `real_interval` (the
  isolation goal and regions, against Sturm counts and certified signs), `real_handoff`, `real_only_realness`,
  `sparse_real_roots`, `method_handoff`, `auto_selection` and `methods` (each explicit method, a one-unit solve's
  partial cover and full solves at 53 and 160 bits). `seed_equivariance` checks that initial candidates rescale
  exactly under x → 2ˢx.
- **The thread tests:** `parallel` and `parallel_solves` (a solve's own workers: the same answer, counters and
  statistics at 1, 3 and 8 workers), and the tests labelled `threads` (`threads`, `cancellation`, `mathematics`,
  `exact_complex`, `solve_state`).
- **The documentation:** `doc_links` (every local link resolves and each documentation folder's index lists its
  pages), `examples_quoted`, and `doc_residue` (the pages for users carry no dates, hashes or words of the project's
  own process).

## ThreadSanitizer

`cmake --preset tsan`, `cmake --build --preset tsan` and `ctest --preset tsan` build the thread tests with GCC's
ThreadSanitizer in `build/tsan` and run them with `tests/tsan.supp`. Set `CERTPOLYSOLVE_TSAN_FLINT_ROOT` to the prefix
of a FLINT built with `-fsanitize=thread` by the same compiler; against an uninstrumented FLINT, TSan sees nothing
inside FLINT. In a fresh copy of FLINT's source tree:

```sh
./configure --prefix=PREFIX --enable-thread-safe --enable-avx2 --disable-static CC=gcc CXX=g++ \
    CFLAGS="-O2 -g -fsanitize=thread" LDFLAGS=-fsanitize=thread && make -j6 && make install
```

GMP and MPFR stay the system's, so races inside them are not seen. Cap the run with a cgroup
(`systemd-run --user --scope -p MemoryMax=…`), not with `prlimit --as`: TSan needs an unlimited address space. The
one suppression is FLINT's benign allocator race.

## Sanitizers

`cmake --preset sanitized` builds with AddressSanitizer and UndefinedBehaviorSanitizer (Clang), and
`ctest --preset sanitized` runs the suite under them; the tests' time limits are scaled for the slower build.
