# Examples

Complete programs, each compiled and run by the test suite. The quick start's five are in `quickstart/`, the rest by
topic in `topics/`; [the quick start](../docs/quickstart-cpp.md) and [the examples page](../docs/examples-cpp.md)
show each with its output. All are C++23 (they print with `std::println`); the C version of the first solve is in
[the guide](../docs/api.md#1-headers-and-a-first-solve).

| Program | What it shows |
|---|---|
| [readme/first_solve.cpp](readme/first_solve.cpp) | The README's first solve: x² − 2, printed as numbers |
| [quickstart/01_first_solve.cpp](quickstart/01_first_solve.cpp) | Every root of x³ − 2x − 5 |
| [quickstart/02_real_roots.cpp](quickstart/02_real_roots.cpp) | The real roots: all of them, those in an interval, the positive ones, and only how many |
| [quickstart/03_many_digits.cpp](quickstart/03_many_digits.cpp) | Exact coefficients, huge coefficients, and many digits |
| [quickstart/04_repeated_and_close.cpp](quickstart/04_repeated_and_close.cpp) | Repeated roots, and roots very close together |
| [quickstart/05_time_limit.cpp](quickstart/05_time_limit.cpp) | A time limit, and stopping a solve from another thread |
| [topics/building.cpp](topics/building.cpp) | Building a polynomial: with arithmetic, from its roots, by a recurrence, term by term |
| [topics/stability.cpp](topics/stability.cpp) | A proof that a closed loop is stable, or that it is not |
| [topics/precision_cap.cpp](topics/precision_cap.cpp) | Two roots 2⁻¹⁰⁰ apart, and the precision cap |
| [topics/noisy_coefficients.cpp](topics/noisy_coefficients.cpp) | Coefficients known only within bounds |
| [topics/eliminant.cpp](topics/eliminant.cpp) | An eliminant from a Gröbner basis or a rational univariate representation |
| [topics/programs.cpp](topics/programs.cpp) | A polynomial given by the steps that compute it, never expanded |
| [topics/secular.cpp](topics/secular.cpp) | Eigenvalues of a diagonal-plus-rank-one matrix, as a secular equation |
| [topics/threads.cpp](topics/threads.cpp) | Many solves at once, on four threads |
| [guide/first_solve.cpp](guide/first_solve.cpp), [guide/first_solve.c](guide/first_solve.c) | The guide's first solve to 128 bits, in C++ and in C |
| [solve.cpp](solve.cpp) | The command-line tool, `certpolysolve_cli` |
