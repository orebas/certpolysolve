# The CertPolySolve API

This guide describes CertPolySolve's interface: every input, option, entry point and result, with the C name beside
each C++ name. Read [the quick start](quickstart-cpp.md) first; it shows the common tasks as complete programs.
[The mathematics](mathematics.md) says what each certificate proves, and [the methods](production-methods.md)
how `CPS_METHOD_AUTO` chooses a route.

## Contents

1. [Headers and a first solve](#1-headers-and-a-first-solve)
2. [Conventions](#2-conventions)
3. [Coefficient polynomials](#3-coefficient-polynomials)
4. [Secular equations](#4-secular-equations)
5. [Programs](#5-programs)
6. [Options](#6-options)
7. [Solving](#7-solving)
8. [Results](#8-results)
9. [Statistics](#9-statistics)
10. [Threads, lifetimes and FLINT](#10-threads-lifetimes-and-flint)
11. [Versions](#11-versions)
12. [The command-line tool](#12-the-command-line-tool)
13. [Changes in 0.5](#13-changes-in-05)

## 1. Headers and a first solve

- `<certpolysolve/certpolysolve.h>` is the C API. It compiles as C99 or later and as C++ (its declarations are
  `extern "C"`). Its names start with `cps_` and `CPS_`.
- `<certpolysolve/certpolysolve.hpp>` is the C++ API, in namespace `certpolysolve`. It needs C++20 (`<format>`,
  `<stop_token>`, concepts) and is header-only on top of the C API. Its advanced parts are the C structs
  themselves: `cps_statistics`, `cps_engine_counters`, `cps_component_progress` and `cps_resource_limits`.
- `<certpolysolve/gmp.hpp>` (optional) makes GMP's `mpz_class` and `mpq_class` numbers and converts a root end to
  `mpq_class` (`to_mpq`). Link `gmpxx` and `gmp` to use it.
- With CMake: `find_package(CertPolySolve REQUIRED)`, then link `CertPolySolve::certpolysolve`.

Both headers are public. Other languages bind through the C API.

The roots of x² − 2 to 128 bits, within 10 seconds, in C (`examples/guide/first_solve.c`, which the suite builds
and runs):

```c
#include <certpolysolve/certpolysolve.h>
#include <stdio.h>

int main(void) {
    cps_polynomial *p = NULL;
    cps_options *o = NULL;
    cps_result *r = NULL;
    int status = 1;
    if (cps_polynomial_create(&p) != CPS_OK ||
        cps_polynomial_set_coefficient(p, 0, "-2", NULL) != CPS_OK || /* x^2 - 2 */
        cps_polynomial_set_coefficient(p, 2, "1", NULL) != CPS_OK)
        goto done;
    if (cps_options_create(&o) != CPS_OK || cps_options_set_precision(o, 128) != CPS_OK ||
        cps_options_set_time_limit(o, 10.0) != CPS_OK)
        goto done;
    if (cps_solve(p, o, &r) != CPS_OK) { /* a usage error: no result */
        fprintf(stderr, "%s\n", cps_last_error());
        goto done;
    }
    if (!cps_result_goal_met(r)) { /* a stop: the regions are proved, the goal is not met */
        fprintf(stderr, "%s: %s\n", cps_termination_name(cps_result_termination(r)), cps_result_diagnostics(r));
        goto done;
    }
    for (size_t i = 0; i < cps_result_size(r); ++i) {
        cps_root_view v;
        char lower[64], upper[64];
        if (cps_result_get_root(r, i, &v) != CPS_OK ||
            cps_dyadic_format(v.real_lower, 40, lower, sizeof lower, NULL) != CPS_OK ||
            cps_dyadic_format(v.real_upper, 40, upper, sizeof upper, NULL) != CPS_OK)
            goto done;
        printf("[%s, %s] holds %llu root(s)\n", lower, upper, (unsigned long long)v.count);
    }
    status = 0;
done:
    cps_result_destroy(r);
    cps_options_destroy(o);
    cps_polynomial_destroy(p);
    return status;
}
```

Its calls, as the header declares them (the header includes `<stddef.h>` and `<stdint.h>` itself):

```c
/* Declared in <certpolysolve/certpolysolve.h>. */
cps_error cps_polynomial_create(cps_polynomial **out);
cps_error cps_polynomial_set_coefficient(cps_polynomial *, uint64_t index, const char *real, const char *imag);
cps_error cps_options_create(cps_options **out);
cps_error cps_options_set_precision(cps_options *, uint64_t bits);
cps_error cps_options_set_time_limit(cps_options *, double seconds);
cps_error cps_solve(const cps_polynomial *, const cps_options *, cps_result **out);
const char *cps_last_error(void);
int cps_result_goal_met(const cps_result *);
cps_termination cps_result_termination(const cps_result *);
const char *cps_termination_name(cps_termination);
const char *cps_result_diagnostics(const cps_result *);
size_t cps_result_size(const cps_result *);
cps_error cps_result_get_root(const cps_result *, size_t index, cps_root_view *out);
cps_error cps_dyadic_format(cps_dyadic_view x, unsigned digits, char *buffer, size_t size, size_t *length);
void cps_result_destroy(cps_result *result);
void cps_options_destroy(cps_options *options);
void cps_polynomial_destroy(cps_polynomial *polynomial);
```

- A region's ends are *dyadic* numbers: m · 2^e with integers m and e, so they are exact. A `cps_dyadic_view` holds
  m and e as two decimal strings, and `cps_dyadic_format` writes the number in decimal with the digits you ask for.
- Text comes back into your buffer by one rule (section 3). A buffer too short is not an error: the text is cut to
  fit, and the length the call reports is longer than the buffer.
- [Approximations](#approximations) in section 8 gives the roots as plain doubles instead.

Without CMake, with the library installed under a prefix such as `$HOME/.local` (`/usr/local` needs no `-I`, `-L`
or `-rpath`):

```sh
cc -std=c99 -I "$HOME/.local/include" first_solve.c -o first_solve \
   -L "$HOME/.local/lib" -Wl,-rpath,"$HOME/.local/lib" -lcertpolysolve
```

The same in C++ (C++23 for `std::println`; `examples/guide/first_solve.cpp`):

```cpp
#include <certpolysolve/certpolysolve.hpp>
#include <chrono>
#include <print>

using namespace std::chrono_literals;

int main() {
    const certpolysolve::Polynomial p{-2, 0, 1}; // x^2 - 2, ascending
    const certpolysolve::Result result = certpolysolve::solve(p, {.bits = 128, .time_limit = 10s});
    if (!result.ok()) { // an error, or a stop with proved regions short of the goal
        std::println(stderr, "{}: {}", result.stop(), result.message());
        return 1;
    }
    for (const certpolysolve::Root &root : result.roots())
        std::println("{} holds {} root(s), centre {:.40}", root.box.real, root.count, root);
}
```

`CPS_OK` from `cps_solve` means only that the call was valid and ran. Whether the goal was met is in the result:
`cps_result_goal_met` and `cps_result_termination`. A solve that stops early still returns proved regions. In C++,
`ok()` says both: no error, and the goal met.

## 2. Conventions

### Errors and outcomes

**C.** A function that can fail returns `cps_error`:

| Code | Meaning |
|---|---|
| `CPS_OK` (0) | The call succeeded. |
| `CPS_INVALID_ARGUMENT` (1) | The call broke a precondition: a NULL required pointer, a malformed literal, a value outside its range, an operand that is not an earlier node, a program without an output at the solve, or a usage error of the solve (section 7). |
| `CPS_OUT_OF_MEMORY` (2) | An allocation of the library's own failed. Nothing was published. (An allocation failure inside FLINT or GMP aborts the process: section 10.) |
| `CPS_INTERNAL_ERROR` (3) | Anything else failed inside the library. No partial result is published. |
| `CPS_OUT_OF_RANGE` (4) | The request is well formed but beyond what the library represents or admits: an index at or past the end of a sequence (a region, a component); an exponent above 2^63 − 1, or arithmetic whose result's degree would be; at the solve, a degree beyond what memory or a coefficient list can address; a precision above 2^63 − 2; a literal longer than addressable storage; a secular equation over its resource limits at creation; a binary64 conversion that overflows or turns a nonzero value into zero; a decimal outside MPFR's exponent range (`cps_dyadic_format`). |
| `CPS_UNRESOLVED` (5) | A binary64 approximation was asked of a region the result does not resolve (section 8). |

- A function of another return type cannot fail. Given NULL, it returns a fallback (below).
- A call that fails changes nothing. A handle output (`T **out`) is set to NULL and a node output to
  `CPS_NO_NODE`; other outputs keep their values; the object operated on is unchanged, since setters and builders
  are transactional.
- `cps_last_error()` describes the calling thread's last failed call. A call that returns `cps_error` sets it: the
  reason when it fails, an empty string when it succeeds. The pointer is valid on the calling thread until that
  thread's next `cps_error`-returning call.
  The text is at most 511 characters, for people; its wording is not part of the contract.

**C++.** Nothing throws except `std::bad_alloc`, when memory runs out (a C call's `CPS_OUT_OF_MEMORY` included).
Printing never throws: a binary exponent too large to print as a fraction, or text the library has not checked
yet, prints as it was given. `to_mpq` (`gmp.hpp`) takes an end as a result gives it; an exponent beyond
`long long` is a number of more than 2^63 bits and gives `std::bad_alloc`, as any allocation that fails does, and
text that is not a dyadic gives `std::invalid_argument`, as gmpxx does for a malformed number.

- An input remembers its first bad value: `ok()` is false and `status()` says why. Arithmetic on a bad
  polynomial passes the mistake on.
- A solve of a bad input, a solve with a bad option, and a solve the C API rejects return a `Result` whose
  `stop()` is `Stop::error`, with the reason in `message()`, no regions and `goal_met()` false.
- An exception thrown by your own function inside `Program::from` propagates unchanged.

**Outcomes are not errors.** A precision stop, a spent budget, a time limit, a cancellation and a refused input
all return `CPS_OK` and a result that says what happened (C++: a `Result` whose `stop()` names it).

### NULL arguments in C

| Where NULL appears | Meaning |
|---|---|
| A handle argument of a `cps_error` function | `CPS_INVALID_ARGUMENT` |
| A handle argument of a function of another return type | That function's fallback (next table) |
| The argument of a destroy function or of `cps_cancel_token_request` | Nothing happens |
| An output pointer | `CPS_INVALID_ARGUMENT` |
| The options of a solve | The defaults |
| The resource limits of `cps_secular_create` | The defaults |
| A real part, a pole's real part, a program text, a tolerance, a cluster radius | `CPS_INVALID_ARGUMENT` |
| An imaginary part or a radius, in every input function | The value 0 |
| A region end | That side is unbounded (the empty string means the same) |
| The token of `cps_options_set_cancel_token` | The options' token is cleared |

| Function given a NULL handle | Returns |
|---|---|
| `cps_result_degree`, `cps_polynomial_degree` | −1 |
| `cps_result_termination` | `CPS_UNSUPPORTED_QUERY` |
| `cps_result_stop_reason` | `CPS_STOP_INPUT` |
| `cps_result_complete`, `cps_result_goal_met`, `cps_result_preparation_limited` | 0 |
| `cps_result_goal` | `CPS_ROOT_ACCURACY` |
| `cps_result_work_bits`, `cps_result_size`, `cps_result_component_count`, `cps_result_sufficient_work_bits` | 0 |
| `cps_result_diagnostics` | the static string "null result" |
| `cps_result_squarefree`, `cps_result_handoff` | `CPS_SQUAREFREE_UNKNOWN`, `CPS_HANDOFF_NONE` |
| `cps_program_size`, `cps_program_output` | 0, `CPS_NO_NODE` |
| an options getter | the setting's default |
| `cps_cancel_token_requested` | 0 |
| `cps_resource_limits_default` | does nothing |

`cps_termination_name` and `cps_stop_reason_name` return "unknown" for a value outside their enums.

### Numbers

Every number given as text uses one grammar: coefficients and radii, secular residues and poles, program
constants, the absolute tolerance, the cluster radius and the region ends.

```text
literal   = space* ( rational | dyadic | decimal ) space*
space     = " " | TAB | CR | LF                         (around the literal only)
rational  = integer "/" integer                         p/q, q not 0
dyadic    = integer "@" integer                         m * 2^e
decimal   = sign? mantissa ( ("e" | "E") integer )?     m * 10^e
mantissa  = digit+ ( "." digit* )? | "." digit+
integer   = sign? digit+
sign      = "+" | "-"
```

- Values are exact: `0.1` is one tenth and `3@-2` is 3/4. Integers and exponents may have any length, and an
  exponent stays compact: `1@100000000000000000000` takes a few bytes. Whether a compact value is ever expanded is
  the resource limits' decision (section 6).
- Rejected with `CPS_INVALID_ARGUMENT`: empty text, space inside a literal, hexadecimal, `inf`, `nan`, expressions,
  a second operator (`1@2@3`) and a zero denominator. Radii, the tolerance and the cluster radius must be ≥ 0.
- A literal longer than addressable storage gives `CPS_OUT_OF_RANGE`.
- A double passed as a double (`cps_polynomial_set_coefficient_d`, and every C++ double) enters at its exact binary
  value. NaN and ±infinity are rejected.

In C++ a `Number` is one exact number: an integer of any integral type, a double (its exact binary value), text in
the grammar above, a `std::complex<double>`, `Number::complex(re, im)` from two real numbers, or a value of a type
with a `number_traits` specialization (`<certpolysolve/gmp.hpp>` gives `mpz_class` and `mpq_class` one). A character
or a `bool` is refused at compile time (`'3'` would be 51). `real()` and `imag()` return the parts as text. Text
is checked when the library uses it, not when the `Number` is made.

### Units

| Quantity | Unit |
|---|---|
| Requested precision (`bits`) | Bits of relative accuracy: every root z in a region with midpoint c satisfies abs(z − c) ≤ 2^−bits · abs(z) |
| Precision cap (`max_bits`), working precision, `work_bits` | Bits of the ball arithmetic's precision |
| Absolute tolerance, cluster radius | A distance in the complex plane, in the roots' own units |
| Time limit | Seconds (C: a `double`; C++: a `std::chrono` duration), on a monotonic clock |
| Times in the statistics | Nanoseconds of wall time, on a monotonic clock |
| Work budget | Controller work units (section 6) |
| Iteration budget | Sweeps of the simultaneous iterations |
| Live cells | Search squares kept at once |
| Resource limits | Bits, word operations and bytes, each a per-stage estimate |
| Region ends | mantissa · 2^exponent, two decimal integer strings |

## 3. Coefficient polynomials

C++ `Polynomial`, C `cps_polynomial`: p(x) = Σ aᵢ xⁱ, coefficients in ascending order.

- A coefficient is an exact complex number, or a ball: the box [re ± r_re] + i [im ± r_im], componentwise, not a
  disk. A coefficient never set is exactly zero.
- The polynomial keeps its nonzero terms only. Any exponent up to 2^63 − 1 can be set, x¹⁰⁰⁰⁰⁰⁰⁰⁰⁰ + 1 is two
  terms, an exact zero removes a term, and a ball whose value is zero is a term.
- The degree is the highest exponent whose coefficient is not exactly zero; a ball counts as nonzero. If the top
  ball contains zero, the solve may end with `CPS_INDETERMINATE_DEGREE`.
- A certificate on ball input holds for every polynomial whose coefficients lie in the balls.
- Every route except the real-only route's sparse recursion reads the coefficients as a list of degree + 1
  entries, so a degree beyond what memory holds fails at the solve (`CPS_OUT_OF_MEMORY`, or `CPS_OUT_OF_RANGE`
  beyond what the list can address), not when the coefficient is set.

### Building a polynomial

| C++ | C | What it makes |
|---|---|---|
| `Polynomial{a0, a1, ...}` | `cps_polynomial_create(&p)`, then a setter per coefficient | the constant term first: `{-2, 0, 1}` is x² − 2 |
| `Polynomial::lowest_first({a0, a1, ...})` | | the same, named |
| `Polynomial::highest_first({an, ..., a0})` | | the highest power first, as MATLAB's `roots` and `numpy.roots` read a list: `{1, 0, -2}` is x² − 2. Print the polynomial to see which you built |
| `Polynomial(range)` (explicit), `lowest_first(range)`, `highest_first(range)` | | from a range of numbers, strings or coefficients |
| `Polynomial(c)`, `Polynomial{c}` | | the constant c |
| `p.set(e, c)` (returns `p`, so calls chain) | `cps_polynomial_set_coefficient(p, e, re, im)`, `cps_polynomial_set_coefficient_d(p, e, re, im)`, `cps_polynomial_set_ball_coefficient(p, e, re, im, re_radius, im_radius)` | sets the coefficient of xᵉ |
| `Polynomial::sparse({{e, c}, ...})`, `Polynomial::sparse(range)` | the setters | from its terms, in any order |
| `Polynomial::x()`, `Polynomial::monomial(e, c = 1)` | | x, and c·xᵉ |
| `ball(value, radius)`, `ball(value, re_radius, im_radius)` | `cps_polynomial_set_ball_coefficient` | a ball coefficient |
| copy construction and assignment | `cps_polynomial_copy(p, &out)` | an independent copy, balls included |
| destructor | `cps_polynomial_destroy(p)` | |

- The C setters read their strings during the call only. The real part must not be NULL; a NULL imaginary part or
  radius is 0. A negative radius or a malformed literal gives `CPS_INVALID_ARGUMENT`, an exponent above 2^63 − 1
  `CPS_OUT_OF_RANGE`; either way the polynomial is unchanged.
- `ball(value, radius)` bounds both parts of a complex value by the radius, and only the real part of a real
  value. `ball(value, re_radius, im_radius)` bounds them separately.
- C++ coefficients are `Coefficient` values: `value` (a `Number`), `real_radius` and `imag_radius` (`Number`s,
  0 for an exact coefficient), and `exact()`. Anything a `Number` is made from converts to one.

```cpp
using certpolysolve::Polynomial;
const Polynomial x = Polynomial::x();
const Polynomial p = pow(x - 1, 3) * (x + 2) - "1/3";      // exact: x^4 - x^3 - 3x^2 + 5x - 7/3
const auto s = Polynomial::sparse({{1000000, 1}, {7, -3}, {0, "1/2"}});
const Polynomial noisy{certpolysolve::ball(2, "1e-9"), certpolysolve::ball(-3, "1e-9"), 1};
const std::vector<std::string> lines{"-2", "0", "0", "1"};
const Polynomial cube(lines);                                // x^3 - 2
```

### Reading a polynomial

| C++ | C |
|---|---|
| `degree()` (−1 for the zero polynomial), `is_zero()` | `cps_polynomial_degree(p)` (−1 for the zero polynomial or NULL) |
| `p[e]`, `coefficient(e)`: a `Coefficient`, zero beyond the degree | `cps_polynomial_get_coefficient(p, e, part, buffer, size, &length)`, `part` one of `CPS_PART_REAL`, `CPS_PART_IMAG`, `CPS_PART_REAL_RADIUS`, `CPS_PART_IMAG_RADIUS` |
| `coefficients()`: degree + 1 coefficients, ascending (for few terms and a high degree, read `terms()`) | |
| `terms()`: the nonzero terms, ascending, as `Polynomial::Term{exponent, coefficient}` | `cps_polynomial_terms(p, exponents, capacity, &count)` |
| `a == b`: the same exact coefficients and radii | |
| `std::format("{}", p)`, `std::println` | |

- A coefficient part comes back as exact text in the literal grammar ("p/q", "p", or "m@e" for a binary
  exponent), so the text sets the same value again.
- Every C call that writes text follows one protocol: at most `size` bytes, the terminating NUL included, go into
  `buffer`; `*length` (if not NULL) receives the full length without the NUL, so `size` 0 asks for the length.
  A short buffer is not an error: the call returns `CPS_OK` with the text cut to `size` − 1 bytes and a NUL, and
  `*length` ≥ `size` says it was cut.
- `cps_polynomial_terms` writes at most `capacity` exponents, ascending; `*count` (if not NULL) receives the number
  of terms, so `capacity` 0 asks for it.
- A polynomial prints exactly, highest power first: `x^4 - x^3 - 3x^2 + 5x - 7/3`, `(3/4 +/- 1/1000)x + 1`, and a
  bad polynomial as `<bad polynomial: reason>`. A `Number` prints as `-1/3` or `1 - 2i`, a `Coefficient` as
  `3/4 +/- 1/1000`.

### Exact arithmetic

| C++ | C |
|---|---|
| `a + b`, `a - b`, `a * b` (numbers mix in: `2 * x - "1/3"`), `-a`, `+=`, `-=`, `*=` | `cps_polynomial_add(a, b, &out)`, `cps_polynomial_subtract`, `cps_polynomial_multiply`; `cps_polynomial_scale(a, re, im, &out)` for a times a number |
| `pow(a, n)` (any polynomial to the power 0 is 1) | `cps_polynomial_power(a, n, &out)` |
| `a.derivative()` | `cps_polynomial_derivative(a, &out)` |
| `Polynomial::from_roots(range)`: a range of numbers, or of pairs (root, multiplicity); also `from_roots({r1, r2, ...})` | `cps_polynomial_from_roots(real, imag, multiplicity, count, &out)` |

- The arithmetic is exact, over the Gaussian rationals. Every coefficient of an operand must be exact: a ball
  coefficient gives `CPS_INVALID_ARGUMENT` (in C++, a bad result polynomial).
- A coefficient's exponent is expanded as written: `"1@1000000"` is a million-bit integer.
- The arithmetic works on the terms. A product, or a step of a power, pairs terms when that pairs no more of them
  than the dense result has coefficients, and is dense otherwise: (x¹⁰⁰⁰⁰⁰⁰ + 1)² stays three terms.
- A result whose degree would exceed 2^63 − 1 gives `CPS_OUT_OF_RANGE`.
- In C each call makes a new polynomial in `*out` (NULL on error); the operands stay unchanged and may be the same
  object.
- `cps_polynomial_from_roots` forms the product of (x − r_k)^m_k with r_k = `real[k]` + i `imag[k]`. NULL `imag`
  makes every root real and a NULL `imag[k]` that root; NULL `multiplicity` means 1 each, and a multiplicity of 0
  leaves the root out. `count` 0 gives the polynomial 1.
- Expanding a product or a recurrence gives up the structure the program route uses: there a product's factors
  are solved one by one, and a recurrence is evaluated without its large, cancelling coefficients. A program
  (section 5) keeps that structure.

## 4. Secular equations

C++ `SecularEquation`, C `cps_secular`: S(x) = Σⱼ aⱼ / (x − bⱼ) − 1, with exact complex residues aⱼ and poles bⱼ.
The poles are excluded from the domain, and the roots are the zeros of S. The degree is the number of terms. The
empty sum is the nonzero constant −1, which has no roots.

| C++ | C |
|---|---|
| `SecularEquation{{.residue = r, .pole = b}, ...}`; `SecularEquation(range)` (explicit) of `SecularEquation::Term{residue, pole}` | `cps_secular_create(terms, count, limits, &out)` with `cps_secular_term_view {residue_real, residue_imag, pole_real, pole_imag}` |
| `ok()`, `status()` | `cps_last_error()` after a failed create |
| copies share the equation, which never changes | `cps_secular_destroy(s)` |

```cpp
const certpolysolve::SecularEquation s{{.residue = "2/3", .pole = 0}, {.residue = "-2/3", .pole = 3}}; // roots 1 and 2
```

```c
cps_secular_term_view terms[2] = {{"2/3", NULL, "0", NULL}, {"-2/3", NULL, "3", NULL}};
cps_secular *s = NULL;
if (cps_secular_create(terms, 2, NULL, &s) != CPS_OK) fprintf(stderr, "%s\n", cps_last_error());
```

- Creation copies the terms and normalizes them to exact rationals; the strings are read during the call only.
  `count` 0 with NULL `terms` is the empty sum. NULL imaginary parts are 0.
- `limits` are the resource limits of section 6 (NULL: the defaults). The C++ `SecularEquation` uses the defaults.
- `CPS_INVALID_ARGUMENT`: NULL `out`; NULL `terms` with `count` > 0; a NULL `residue_real` or `pole_real`; a
  malformed literal; a zero residue; two equal poles (compared exactly, so `1/2` equals `0.5`); an unknown flag.
- `CPS_OUT_OF_RANGE`: `count` beyond addressable storage; the estimated term storage over `temporary_bytes`;
  `count` above `exact_operations` / 4; an exponent expansion beyond the compact-exponent rule without
  `CPS_RESOURCE_ALLOW_EXPONENT_EXPANSION`; a normalization over min(`exact_materialization_bits`,
  8 · (`temporary_bytes` − the term storage estimate)).
- A complex query iterates and certifies on the terms, without expanding the equation into a polynomial, and can
  succeed with no exact allowance at all. A real-only query, and a complex query that this does not finish, use
  the exact numerator polynomial Π(x − bⱼ) · S(x), under the solve's resource limits; a refused conversion is the
  outcome `CPS_UNSUPPORTED_QUERY`, with `cps_result_preparation_limited` set.

## 5. Programs

C++ `Program`, C `cps_program`: a polynomial given by a straight-line program that computes it, for polynomials
that are easy to compute and enormous to write out (iterations, recurrences, products).

- A program is a list of nodes, numbered 0, 1, 2, ... in the order they are appended. A node is x, an exact complex
  constant, a + b, a − b, a · b, −a, or aᵉ with e ≥ 2; its operands are earlier nodes. Several x nodes all
  denote x.
- One node is the output, and the polynomial is its value.
- A solve never expands the program. It proves the degree and the leading coefficient from the program's top
  coefficients, bounds the roots by a majorant, and iterates on evaluations of the program.

### Building a program in C

| Call | Effect |
|---|---|
| `cps_program_create(&p)`, `cps_program_destroy(p)` | an empty program, with no output |
| `cps_program_variable(p, &node)` | appends x |
| `cps_program_constant(p, re, im, &node)` | appends the constant re + i im (literals; NULL `im` is 0) |
| `cps_program_add(p, a, b, &node)`, `cps_program_subtract`, `cps_program_multiply` | appends a + b, a − b, a · b |
| `cps_program_negate(p, a, &node)` | appends −a |
| `cps_program_power(p, a, e, &node)` | appends aᵉ; e must be at least 2 |
| `cps_program_set_output(p, node)` | makes `node` the output, replacing an earlier one |
| `cps_program_parse(text, &p)` | a program from the text form below |
| `cps_program_size(p)`, `cps_program_output(p)` | the number of nodes; the output node, or `CPS_NO_NODE` (`UINT64_MAX`) when none is set |

- Each builder writes the new node's index to `*node`, or `CPS_NO_NODE` when it fails, and a failed call leaves the
  program unchanged.
- `CPS_INVALID_ARGUMENT`: a NULL program or output pointer; a NULL constant real part or a malformed literal; an
  operand that is not an earlier node; an exponent below 2; an output that is not a node.

```c
/* (x^2 + x + 1)^5 - 1/3 */
uint64_t x, square, sum, one, base, power, third, out;
cps_program *q = NULL;
if (cps_program_create(&q) == CPS_OK && cps_program_variable(q, &x) == CPS_OK &&
    cps_program_multiply(q, x, x, &square) == CPS_OK && cps_program_add(q, square, x, &sum) == CPS_OK &&
    cps_program_constant(q, "1", NULL, &one) == CPS_OK && cps_program_add(q, sum, one, &base) == CPS_OK &&
    cps_program_power(q, base, 5, &power) == CPS_OK && cps_program_constant(q, "1/3", NULL, &third) == CPS_OK &&
    cps_program_subtract(q, power, third, &out) == CPS_OK && cps_program_set_output(q, out) == CPS_OK) {
    /* cps_solve_program(q, options, &result) */
}
cps_program_destroy(q);
```

### The text form

```text
PSPROGRAM 1
x            # node 0
mul 0 0      # node 1: x^2
add 1 0      # node 2: x^2 + x
c 1          # node 3: the constant 1 (an optional second field is the imaginary part)
add 2 3      # node 4
pow 4 5      # node 5
c 1/3        # node 6
sub 5 6      # node 7
out 7
```

- Lines end with LF; a trailing CR is ignored. `#` starts a comment that runs to the end of the line. Fields are
  separated by spaces or tabs. Blank lines are ignored.
- The first line that is not blank must be `PSPROGRAM 1`. Each later line is one operation: `x`;
  `c REAL [IMAG]`; `add A B`, `sub A B`, `mul A B`; `neg A`; `pow A E` (E ≥ 2); `out N`. A, B, N and E are unsigned
  decimal integers.
- Every operation but `out` appends one node, numbered in order. The text needs at least one `out` line; the last
  one wins.
- A NULL text, a syntax error, a bad operand or a missing `out` line gives `CPS_INVALID_ARGUMENT`;
  `cps_last_error()` names the line where there is one.

### The C++ recorder

```cpp
auto power = certpolysolve::Program::from([](auto x) { return pow(x * x + x + 1, 5) - "1/3"; });
auto mandelbrot = certpolysolve::Program::from([](auto c) {
    auto p = c;
    for (int k = 1; k < 8; ++k) p = p * p + c; // degree 128; its coefficients are never formed
    return p;
});
const certpolysolve::Program parsed("PSPROGRAM 1\nx\nmul 0 0\nout 1\n"); // the text form
```

- `Program::from(f)` calls f once with a `ProgramValue` for x and records the arithmetic f does. The value f
  returns becomes the output: a `ProgramValue` of this recording, or a plain number (a constant).
- A `ProgramValue` has `+`, `-`, `*`, unary `-`, `+=`, `-=`, `*=` and `pow(value, n)` (n = 0 gives the
  constant 1 and n = 1 the value itself; the text form's `pow` takes n ≥ 2). Numbers mix in
  (integers, doubles at their exact value, text, complex numbers), and each mixed operation appends a constant
  node.
- It has no comparison and no conversion to `bool`, so f cannot branch on the value of x; such code does not
  compile. Loops and branches over integers, helper functions and captured values are fine.
- Use a `ProgramValue` only inside the call of f that produced it: an operand recorded by another program is
  refused (its node numbers mean nothing here).
- The first failure (a malformed constant, an operand or a returned value that this recording did not produce) is
  remembered: `ok()` is false and `status()` says why, and a solve reports it as `Stop::error`.
- `Program(text)` parses the text form. A default-constructed `Program` is empty and has no output, and a solve of
  it is an error.
- Copies of a `Program` share it.

### Solving a program

- Methods: `CPS_METHOD_AUTO` and `CPS_METHOD_CLASSICAL` iterate in Aberth mode, `CPS_METHOD_SECULAR` in secular mode;
  `CPS_METHOD_PIECEWISE` is a usage error.
- All complex roots only: a real-only query, the isolation goal and a region are usage errors (section 7).
- Realness. A region off the real axis is `CPS_NONREAL`. When every constant of the program is real, a region holding
  one distinct root whose reflection in the real axis meets no other region is `CPS_REAL` (by conjugation, as for
  coefficient input; the region is a real interval). Any other region meeting the axis is `CPS_REAL_UNKNOWN`, and so
  is every such region of a program with a non-real constant.
- Goals. The default goal is `CPS_ROOT_APPROXIMATION`, counted clusters allowed. `CPS_CLUSTER_RADIUS` is answered
  under all three methods. `CPS_ROOT_ACCURACY` is a best effort: a cluster that does not separate within the
  precision cap ends the call with `CPS_PRECISION_LIMIT` and the counted cluster, and `cps_result_squarefree` (and
  the diagnostics) say whether the roots are known to be distinct: "distinct roots, unresolved at the cap" when
  the program's image modulo a word-size prime proved the polynomial square-free, "possibly a multiple root"
  otherwise.
- A multiplicity above one is certified only through a power in the program (the output c·qᵏ, or a factor qᵏ of a
  product): q's regions, each with the count and the multiplicity q's region proves, times k. Or when every root
  is zero (the program computes c·xⁿ, proved from its majorant and its constants' denominators): the point 0,
  with count and multiplicity n.
- Outcomes that depend on the program: `CPS_ZERO_POLYNOMIAL` for a program that computes zero;
  `CPS_INDETERMINATE_DEGREE` when its leading terms cancel beyond what the library resolves; `CPS_UNSUPPORTED_QUERY`
  with `CPS_STOP_INPUT` when its length or the solve's storage exceeds the resource limits. A nonzero constant is
  answered as degree 0 with no roots.
- The budgets, the time limit, the token and the iteration budget apply as for coefficients. The solve copies the
  program, so it may be destroyed after the call.

## 6. Options

C++ `Options` is a struct of named fields, each with its default; set them by name (`o.bits = 64;`) or with
designated initializers in declaration order (`{.bits = 64, .time_limit = 60s}`). C `cps_options` is made by
`cps_options_create(&o)` with every default, copied by `cps_options_copy(o, &out)` (an independent copy of every
setting, naming the same token) and freed by `cps_options_destroy(o)`.

### Settings and defaults

| Setting | C++ field | C setter / getter | Default | Valid values |
|---|---|---|---|---|
| Requested precision | `bits` | `cps_options_set_precision` / `cps_options_precision` | `CPS_DEFAULT_BITS` = 53 | 1 to 2^63 − 2 |
| Precision cap | `max_bits` | `cps_options_set_max_work_bits` / `cps_options_max_work_bits` | `CPS_DEFAULT_MAX_WORK_BITS` = 65,536 | `CPS_MIN_MAX_WORK_BITS` (32) to 2^63 − 2; may be below the requested precision |
| Time limit | `time_limit` (`std::chrono::duration<double>`) | `cps_options_set_time_limit` (seconds) / `cps_options_time_limit` | none | 0 or +∞ for none, else a positive number of seconds |
| Goal | `goal` (`std::optional<Goal>`; empty: the input kind's default) | `cps_options_set_goal` / `cps_options_goal` | `CPS_GOAL_DEFAULT` | a `cps_goal` |
| Cluster radius | `Goal::clusters(r)` | `cps_options_set_cluster_radius` / `cps_options_cluster_radius` | none ("0") | a literal ≥ 0 |
| Method | `method` | `cps_options_set_method` / `cps_options_method` | `Method::automatic`, `CPS_METHOD_AUTO` | `automatic`, `classical`, `secular`, `piecewise` |
| Cancellation (advanced) | `stop` (`std::stop_token`) | `cps_options_set_cancel_token` | none | a token, or none |
| Work budget (advanced) | `work_budget` (empty or 0: automatic) | `cps_options_set_work_budget` / `cps_options_work_budget` | automatic: max(`CPS_DEFAULT_WORK_UNITS` = 100,000, `CPS_DEFAULT_WORK_UNITS_PER_ROOT` = 64 × degree) | any `uint64_t`; 0 selects the automatic budget |
| Live cells (advanced) | `cell_budget` (empty or 0: the default) | `cps_options_set_live_cells` / `cps_options_live_cells` | `CPS_DEFAULT_LIVE_CELLS` = 65,536 | any `uint64_t`; 0 selects the default |
| Iteration budget (advanced) | `iteration_budget` (empty: none) | `cps_options_set_max_iterations` / `cps_options_max_iterations` | none (0) | any `uint64_t`; 0 is none |
| Absolute tolerance (advanced) | `absolute_tolerance` (`std::optional<Number>`) | `cps_options_set_absolute_tolerance` / `cps_options_absolute_tolerance` | none ("0") | a literal ≥ 0; 0 disables it |
| Resource limits (advanced) | `resources` (`std::optional<ResourceLimits>`) | `cps_options_set_resource_limits` / `cps_options_get_resource_limits` | every allowance `UINT64_MAX`, no flag | any values; the only flag is `CPS_RESOURCE_ALLOW_EXPONENT_EXPANSION` |
| Real-only | `real_roots(...)` | `cps_options_set_real_only` / `cps_options_real_only` | off | on or off |
| Region | `real_roots(p, region, ...)` | `cps_options_set_region` / `cps_options_get_region` | none | four literals, each side may be unbounded |
| Workers (advanced) | `workers` | `cps_options_set_workers` / `cps_options_workers` | 1: no thread started | 1 or more; at most the machine's hardware threads are used |

- "Degree" in the work budget is the coefficient input's degree, a secular equation's number of terms, or a
  program's degree bound (which saturates).
- The default cap, 65,536 bits, covers every input of the 769-input benchmark comparison (the most any needed was
  59,904 bits, for a polynomial with four clusters of 192 roots). A solve that needs more stops at the cap with a
  complete cover, and for exact real input names a cap that suffices (`cps_result_sufficient_work_bits`).
- A setter checks only its own arguments: NULL, the literal, the sign and the range. A failing setter leaves the
  options unchanged. `cps_options_set_method` and `cps_options_set_goal` take an `int`, so a value outside the enum is
  rejected before any enum value exists. Checks that involve two settings, or a setting and the input kind, are
  the solve's (section 7).
- A getter returns the setting; given NULL, the default. `cps_options_work_budget` returns 0 while the budget is
  automatic, `cps_options_live_cells` the cells in force (65,536 by default), and `cps_options_time_limit` 0 when
  there is none (also after +∞ was set). A string getter returns the literal as last set, trimmed of surrounding
  space ("0" when none is set; NULL for an unbounded region side), valid until the options change or are
  destroyed. `cps_options_get_region` and `cps_options_get_resource_limits` give `CPS_INVALID_ARGUMENT` for NULL.
- C++ applies the fields at the solve: a bad value (`.bits = 0`, a negative time limit, `.max_bits = 16`) makes
  the solve's `Result` an error with the setter's reason.
- Workers: a count above 1 lets each solve call start up to that many threads less one and join them before it
  returns; the answer, its counters and its diagnostics are the same at every count. One worker is the default,
  because a program that runs many solves at once already fills the machine. Large inputs gain most (about 1.4 at
  eight workers at degree 24,576 on forced piecewise); inputs below degree 8,192 gain little. See section 10.

```cpp
certpolysolve::Options options;
options.bits = 128;                                  // relative accuracy of every root
options.max_bits = 4096;                             // the most working precision the solver may use
options.time_limit = std::chrono::seconds(30);
options.goal = certpolysolve::Goal::approximation();
auto result = certpolysolve::solve(p, options);
auto clusters = certpolysolve::solve(p, {.goal = certpolysolve::Goal::clusters("1e-8")});
```

### Goals

Each goal requires `complete` and the condition below. "The accuracy test" for a region with midpoint c: every
root z in it satisfies abs(z − c) ≤ 2^−bits · abs(z) or, with a positive absolute tolerance t, abs(z − c) ≤ t;
it is checked with outward bounds.

| C++ | C | Each region |
|---|---|---|
| `Goal::accuracy()` | `CPS_ROOT_ACCURACY` | holds one distinct root with its multiplicity established, and passes the accuracy test |
| `Goal::approximation()` | `CPS_ROOT_APPROXIMATION` | passes the accuracy test; counted clusters are allowed (count above 1, multiplicity unknown) |
| `Goal::clusters(r)` | `CPS_CLUSTER_RADIUS` with `cps_options_set_cluster_radius(o, "r")` | lies within the absolute distance r of its midpoint (checked by the sum of the half-widths); no separation or multiplicity is required, and Ccluster's 3-disk separation is not promised |
| `Goal::isolation()` | `CPS_ROOT_ISOLATION` | holds one distinct root alone, with its multiplicity established, in a disjoint closed rectangle (a real interval for a real-only query); exact coefficient input and Auto only; no width is promised and the requested precision is not used |
| empty `goal` | `CPS_GOAL_DEFAULT` (4) | the default for the input kind: `CPS_ROOT_ACCURACY` for coefficients and secular equations, `CPS_ROOT_APPROXIMATION` for programs |

- `cps_result_goal` reports the goal the solve applied, never `CPS_GOAL_DEFAULT`.
- Programs default to approximation because a program can never certify a multiplicity above one except through
  a power (section 5), so accuracy is a best effort there.
- A real-only query answers the approximation goal with individual roots, as for accuracy.
- Complex isolation uses exact preparation followed by subdivision. Caller resource limits can stop it with
  a counted cover whose roots are not yet isolated; the result distinguishes complete coverage from a met goal.
  If exact factor preparation does not fit its resource allowance, the query is refused with
  `CPS_UNSUPPORTED_QUERY` and `preparation_limited()` set.
- Under the cluster-radius goal the requested precision and the cap still steer the precision schedule; success is
  decided by the radius alone.
- The goal and the radius are coupled in C. A positive radius sets the goal to `CPS_CLUSTER_RADIUS`; a zero radius
  clears it and sets `CPS_GOAL_DEFAULT`. `cps_options_set_goal(o, CPS_CLUSTER_RADIUS)` needs a positive radius set
  first (`CPS_INVALID_ARGUMENT` otherwise), and every other goal clears the radius. In C++, `Goal::clusters(0)` is
  the default goal.
- Settings that the query does not use are ignored, not rejected: the tolerance under the cluster-radius goal, the
  requested precision under the isolation goal.

### Budgets

**Work budget.** A soft, deterministic limit on controller work units.
- A unit is a degree-proof or preparation pass, a backend call, a bounded local certification or refinement
  attempt, a disk-exclusion attempt, a square split or a frontier transition.
- It is checked between units, and the unit in progress is finished: at high precision one unit can take long.
  Nothing is promised about how far a solve gets within it.
- It belongs to the call: every engine of the call shares it, a method's hand-off and the solves a real-only
  query runs included.
- A solve that spends it returns its certified cover with `CPS_WORK_LIMIT` and `CPS_STOP_WORK_BUDGET`.
- The automatic budget, 100,000 units or 64 a root if that is more, gives a degree-8,192 input 524,288.

**Live cells.** The most search squares subdivision keeps at once, parents and untested children included. Not a
byte limit. A search it blocks ends with `CPS_WORK_LIMIT` and `CPS_STOP_CELL_BUDGET`.

**Iteration budget.** At most this many sweeps of the simultaneous iterations: the Aberth and secular packets of
the classical and secular methods, at every precision, the binary64 ones included, counted by
`cps_statistics.sweeps`.
- Not under it: subdivision and the exact real-root isolation, which make no sweeps; the real-only route's
  attempts (the sparse recursion and the axis charts); and the piecewise method's own work (its solves of small
  local models appear in `sweeps` but are not sweeps of the input). Their limits are the work and cell budgets,
  the time limit and the token. A solve the piecewise method hands to the secular method is under the budget from
  the hand-off on.
- A spent budget returns the last certified cover with `CPS_WORK_LIMIT` and `CPS_STOP_ITERATION_BUDGET`. The value
  changes no engine's own packet lengths.
- It is not a time bound: a sweep at tens of thousands of bits costs as much as a million binary64 sweeps. Use the
  time limit for a deadline.

### The time limit and the cancellation token

**Time limit.** A wall-clock limit for one solve, on a monotonic clock from the moment the solve call is entered.
- It covers everything the call does: every engine, a method's hand-off, and the solves of the real-only route.
- The solve checks it at its cooperative checkpoints and stops at the first one after the limit, with
  `CPS_CANCELLED`, `CPS_STOP_TIME_LIMIT` (C++ `Stop::time_limit`) and its certified partial answer. The overrun is
  bounded by the longest step that cannot be interrupted: one FLINT call, or one work unit.
- A solve whose goal is met before the limit is seen reports `CPS_SUCCESS`.
- An answer under a time limit depends on the machine and its load. The work and iteration budgets are
  deterministic.

**Cancellation (advanced).** In C++, `Options::stop` takes a `std::stop_token`; when its `std::stop_source` is asked
to stop, from any thread, the solve stops at its next checkpoint with `Stop::cancelled`. The C++ solve makes a
library token for the call and requests it from a `std::stop_callback`, which runs on the thread that asks to
stop, not inside the solve.

In C the token is a flag the library owns:

| C | |
|---|---|
| `cps_cancel_token_create(&t)` | a new token, not requested |
| `cps_cancel_token_request(t)` | requests it: from any thread, at any time, also while solves read it; lock-free; sticky |
| `cps_cancel_token_requested(t)` | 1 once requested |
| `cps_cancel_token_destroy(t)` | releases the caller's handle |
| `cps_options_set_cancel_token(o, t)` | names the token that solves with these options read; NULL clears it |

- Every solve whose options name the token stops at its next checkpoint with `CPS_CANCELLED`, `CPS_STOP_CANCELLED`
  and its certified partial answer. A token requested before the solve starts stops it at once.
- A solve whose goal is met reports `CPS_SUCCESS`, whatever is requested afterwards.
- Options that name a token share it, and so do their copies, so the caller may destroy its handle at any time.
  One token may serve any number of options and solves. To cancel again, make a new token.
- No caller code runs inside a solve: cancellation is the solver reading a flag.

### Resource limits

Advanced: most callers leave them unspecified. Each is a per-stage estimate, not a guarantee of peak memory or
wall time.

```c
typedef struct cps_resource_limits {
    uint64_t exact_materialization_bits; /* expanded coefficients, exact operand growth, exact preparation's storage */
    uint64_t exact_operations;           /* exact preparation's effort, in word operations */
    uint64_t temporary_bytes;            /* estimated engine, preparation and construction storage */
    uint64_t flags;                      /* CPS_RESOURCE_ALLOW_EXPONENT_EXPANSION, or 0 */
} cps_resource_limits;
```

- `cps_resource_limits_default(&limits)` sets every allowance to `UINT64_MAX` (unspecified: no fixed ceiling) and
  no flag. In C++, `ResourceLimits` derives from the struct and its constructor does the same.
- 0 is an explicit zero allowance. `exact_materialization_bits` = 0 disables exact preparation, so the solve works
  on balls.
- The exact storage policy: the exact storage limit is min(`exact_materialization_bits`, 8 · `temporary_bytes`);
  exact preparation writes compact coefficients out when the written-out numbers exceed their owned representation
  by at most one precision cap (`max_work_bits`) each, and a region's end is expanded only when the expanded value is
  no larger than its owned representation. `CPS_RESOURCE_ALLOW_EXPONENT_EXPANSION` permits any expansion within the
  limit.
- A refused exact path falls back to ball arithmetic, and `cps_result_preparation_limited` reports it, since
  multiplicities may then be unknown.
- `cps_options_set_resource_limits` copies the struct; NULL or an unknown flag gives `CPS_INVALID_ARGUMENT`.
- [The resource policy](resource-policy.md) details the estimates.

```cpp
certpolysolve::ResourceLimits no_exact;
no_exact.exact_materialization_bits = 0;
auto result = certpolysolve::solve(p, {.goal = certpolysolve::Goal::clusters("1e-8"), .resources = no_exact});
```

## 7. Solving

### Entry points

| C | C++ | Input |
|---|---|---|
| `cps_solve(p, options, &result)` | `solve(const Polynomial &, const Options & = {})`; `real_roots(const Polynomial &, const Options & = {})`; `real_roots(const Polynomial &, const Region &, const Options & = {})` | a coefficient polynomial |
| `cps_solve_secular(s, options, &result)` | `solve(const SecularEquation &, const Options & = {})`; `real_roots(const SecularEquation &, const Options & = {})` | a secular equation |
| `cps_solve_program(q, options, &result)` | `solve(const Program &, const Options & = {})` | a program, which must have an output |

- C: NULL options mean the defaults. On `CPS_OK`, `*result` is a new result that the caller owns, independent of
  the input and the options. On an error `*result` is NULL and nothing is published: `CPS_INVALID_ARGUMENT` for a
  NULL input or output, a program without an output and the usage errors below; `CPS_OUT_OF_MEMORY`;
  `CPS_OUT_OF_RANGE` for a degree beyond what a coefficient list can address; `CPS_INTERNAL_ERROR`.
- C++: `solve` and `real_roots` return a `Result` and throw nothing but `std::bad_alloc`. `real_roots` is the solve
  with real-only switched on; there is none for programs.
- The solve reads its input and options during the call. Neither may be changed or destroyed until it returns;
  other threads may read them meanwhile. The result does not refer to them.

### What each input kind and method supports

Three outcomes appear in the tables:
- **answered**: the solve works on the query;
- **usage error**: `CPS_INVALID_ARGUMENT` from the solve, no result, the reason in `cps_last_error()` (C++:
  `Stop::error`);
- **refusal**: `CPS_OK` with termination `CPS_UNSUPPORTED_QUERY` and stop reason `CPS_STOP_INPUT` (C++:
  `Stop::unsupported`).

A usage error is a combination of options that no input of the entry point's kind can answer; the solve finds it
before any work. A refusal depends on the input. "Ball input" is a coefficient polynomial with at least one nonzero
radius; "the default goal" is `CPS_GOAL_DEFAULT`.

**Coefficients** (`cps_solve`; C++ `solve` and `real_roots` of a `Polynomial`):

| Query | `CPS_METHOD_AUTO` | `CPS_METHOD_CLASSICAL`, `_SECULAR`, `_PIECEWISE` |
|---|---|---|
| All complex roots; accuracy, approximation or the default goal | answered (exact or ball input) | answered on exact input; refusal on ball input |
| All complex roots; cluster radius | answered (exact or ball input) | usage error |
| All complex roots; isolation | answered on exact input; refusal on ball input | usage error |
| Real-only; accuracy, approximation, isolation or the default goal; with or without a region | answered on exact real input; refusal on ball input or nonreal input | usage error |
| Real-only; cluster radius | usage error | usage error |
| A region without real-only | usage error | usage error |

**Secular equations** (`cps_solve_secular`; C++ `solve` and `real_roots` of a `SecularEquation`):

| Query | `CPS_METHOD_AUTO` | `CPS_METHOD_SECULAR` | `CPS_METHOD_CLASSICAL`, `_PIECEWISE` |
|---|---|---|---|
| All complex roots; accuracy, approximation or the default goal | answered | answered | usage error |
| All complex roots; cluster radius | answered | usage error | usage error |
| Real-only; accuracy, approximation or the default goal | answered when the numerator is real and within the exact limits; refusal otherwise | usage error | usage error |
| Real-only with the cluster radius; the isolation goal; a region | usage error | usage error | usage error |

**Programs** (`cps_solve_program`; C++ `solve` of a `Program`):

| Query | `CPS_METHOD_AUTO`, `_CLASSICAL`, `_SECULAR` | `CPS_METHOD_PIECEWISE` |
|---|---|---|
| All complex roots; approximation (the default), accuracy (a best effort) or cluster radius | answered | usage error |
| Real-only (any goal); the isolation goal; a region | usage error | usage error |

**Outcomes that depend on the input**, at every entry point, all with stop reason `CPS_STOP_INPUT`:

| Situation | Termination |
|---|---|
| The input is zero: every coefficient exactly zero, or a program that computes zero | `CPS_ZERO_POLYNOMIAL` |
| The degree is not proved: the top coefficient's ball still contains zero at the precision cap, or a program's leading terms cancel beyond what the library resolves | `CPS_INDETERMINATE_DEGREE` |
| A refusal cell above | `CPS_UNSUPPORTED_QUERY` |
| Over the resource limits: a program's length or solve storage; a secular equation's exact numerator; the exact test at a region's end, when it needs the input in exact integer form beyond the exact storage policy | `CPS_UNSUPPORTED_QUERY` |

A nonzero constant input is answered as degree 0, with no regions, complete, the goal met and `CPS_SUCCESS`: a
constant polynomial, the empty secular sum and a constant program alike.

### Methods

- `CPS_METHOD_AUTO` (C++ `Method::automatic`) chooses a route by a fixed rule that reads only the input and the
  query. The choice is not part of the contract; [the methods](production-methods.md) describe it.
- A forced method (`classical`, `secular`, `piecewise`) names the principal numerical method. Shared exact
  preparation, direct formulas and the method's internal stages remain part of it.
- Forced piecewise hands to the secular method, within the same call and budgets, every stop that is not the
  caller's: its regime exit, a storage refusal, a refusal because the coefficients' exponents exceed its working
  range, and its precision budget. The regime exit is Auto's where Auto would choose piecewise (piecewise leaves
  once its work overruns the selector's estimate by a global solve's minimum budget), and elsewhere piecewise leaves
  once its work passes that budget: so on an input Auto sends to piecewise, forcing piecewise gives Auto's piecewise
  route, and on most others the answer comes from secular after a short attempt. `cps_result_handoff` reports the
  hand-off and its reason. The caller's own stops stay with piecewise: a spent work budget, the time limit, the
  token.

### Real-only queries and regions

A real-only query asks for the real roots of an exact real input. When its goal is met, each comes alone, with its
multiplicity, in a region proved real.

- On coefficient input under Auto, the answer comes from the first of three routes that answers: the sparse
  recursion, on an input with few terms where its predicted cost is below that of the route Auto would take; the
  piecewise method's axis charts, where Auto would choose the piecewise method; the complex route, whose regions
  proved real are the answer. The first two are attempts, made only when the work budget is at least 64 units a
  root and the region (if any) meets the real axis; an attempt that stops hands the query on. The work of every
  route counts against the call's budgets. Which route answered is not part of the contract (the diagnostics name
  it).
- On a secular equation, the solve converts the equation to its exact numerator polynomial and isolates the real
  roots exactly.
- Under the isolation goal the complex route solves at the default 53 bits, so isolation costs what the default
  accuracy query costs.

A region is a closed rectangle [real_lower, real_upper] × [imag_lower, imag_upper].

| C++ | C |
|---|---|
| `real_roots(p, between(a, b))`, `at_least(a)`, `at_most(b)`: ends as numbers; the imaginary sides unbounded | `cps_options_set_real_only(o, 1)` and `cps_options_set_region(o, real_lower, real_upper, imag_lower, imag_upper)`: literals, NULL or "" for an unbounded side |

```cpp
auto window = certpolysolve::real_roots(p, certpolysolve::between("-1/2", "1/2")); // closed: an end that is a root counts
auto positive = certpolysolve::real_roots(p, certpolysolve::at_least(0));
auto intervals = certpolysolve::real_roots(p, {.goal = certpolysolve::Goal::isolation()}); // no width asked
```

In C, the real roots of 16x⁵ − 20x³ + 5x − 1/2 between 0 and 1 (`examples/guide/real_roots_region.c`):

```c
#include <certpolysolve/certpolysolve.h>
#include <stdio.h>

int main(void) {
    cps_polynomial *p = NULL;
    cps_options *o = NULL;
    cps_result *r = NULL;
    int status = 1;
    if (cps_polynomial_create(&p) != CPS_OK ||
        cps_polynomial_set_coefficient(p, 0, "-1/2", NULL) != CPS_OK ||
        cps_polynomial_set_coefficient(p, 1, "5", NULL) != CPS_OK ||
        cps_polynomial_set_coefficient(p, 3, "-20", NULL) != CPS_OK ||
        cps_polynomial_set_coefficient(p, 5, "16", NULL) != CPS_OK)
        goto done;
    if (cps_options_create(&o) != CPS_OK || cps_options_set_real_only(o, 1) != CPS_OK ||
        cps_options_set_region(o, "0", "1", NULL, NULL) != CPS_OK) /* 0 <= x <= 1; NULL: no bound */
        goto done;
    if (cps_solve(p, o, &r) != CPS_OK) {
        fprintf(stderr, "%s\n", cps_last_error());
        goto done;
    }
    if (!cps_result_goal_met(r)) {
        fprintf(stderr, "%s: %s\n", cps_termination_name(cps_result_termination(r)), cps_result_diagnostics(r));
        goto done;
    }
    for (size_t i = 0; i < cps_result_size(r); ++i) {
        cps_root_view v;
        char lower[32], upper[32];
        if (cps_result_get_root(r, i, &v) != CPS_OK ||
            cps_dyadic_format(v.real_lower, 20, lower, sizeof lower, NULL) != CPS_OK ||
            cps_dyadic_format(v.real_upper, 20, upper, sizeof upper, NULL) != CPS_OK)
            goto done;
        printf("[%s, %s]\n", lower, upper);
    }
    status = 0;
done:
    cps_result_destroy(r);
    cps_options_destroy(o);
    cps_polynomial_destroy(p);
    return status;
}
```

```
[0.10452846326765347104, 0.10452846326765347176]
[0.49999999999999999827, 0.50000000000000000173]
[0.97814760073380563453, 0.97814760073380564132]
```

```c
/* Declared in <certpolysolve/certpolysolve.h>. */
cps_error cps_options_set_real_only(cps_options *, int enabled);
cps_error cps_options_set_region(cps_options *, const char *real_lower, const char *real_upper, const char *imag_lower, const char *imag_upper);
```

The root 1/2 comes in an interval about 2⁻⁵³ wide around it, as every root does at the default precision; only a
root at an end of the region comes as an exact point.

- The query answers exactly the real roots x with real_lower ≤ x ≤ real_upper, when imag_lower ≤ 0 ≤
  imag_upper. A rectangle that misses the real axis gets an empty, complete answer, without work. `complete` and
  `goal_met` refer to the roots in the region.
- The ends are closed. A root exactly at an end is decided by an exact test and is in the answer; when the end is
  dyadic, its region is that exact point.
- A root of odd multiplicity whose interval holds an end it is not equal to is placed on its side by the input's
  signs, and its interval is cut at a dyadic point between the end and the root. A root of even multiplicity is
  refined until its interval lies on one side; on the complex route that refines every root.
- A region saves no time on the complex route.
- The setter checks each end's literal only; four unbounded sides clear the region. The solve checks the rest
  before any work, each failure a usage error: on each axis the lower end must not exceed the upper (decided
  exactly within the precision cap), and each real end must expand to an exact rational within the exact storage
  policy.

## 8. Results

C++ `Result`, C `cps_result`. A result never changes and owns everything it exposes. C++ results copy, and copies
share the C result.

### The outcome

| C++ | C | |
|---|---|---|
| `ok()` | `cps_result_goal_met(r)` and `CPS_SUCCESS` | the goal is met (C++: and no error) |
| `stop()` | `cps_result_termination(r)`, `cps_result_stop_reason(r)` | why the solve stopped |
| `complete()` | `cps_result_complete(r)` | the regions hold every root asked for |
| `goal_met()` | `cps_result_goal_met(r)` | the goal is met |
| `degree()` | `cps_result_degree(r)` | the input's degree |
| `work_bits()` | `cps_result_work_bits(r)` | the last outer working precision, not the backend's peak |
| `message()` | `cps_result_diagnostics(r)` | the solver's trace, or the error |
| `roots()` | `cps_result_size(r)`, `cps_result_get_root(r, i, &view)` | the regions |
| `values()` | `cps_result_get_approximation(r, i, &re, &im)` for each region | the roots as numbers: each region's `value`, repeated as many times as the region holds roots (a double root gives two copies), in the regions' order. Doubles, not enclosures. A region without a value is left out, so the list has `degree()` numbers exactly when every region has one; after a stop or an error it is empty |
| (`cps_result_goal(result.get())`) | `cps_result_goal(r)` | the goal applied |
| `get()` | | the C result, for the C calls; NULL for an error result |
| destructor | `cps_result_destroy(r)` | |

`degree` is −1 for the zero input, and an upper bound with `CPS_INDETERMINATE_DEGREE`. For a secular equation it is
the number of terms; for a program, the proved degree, or its degree bound when the program was refused. A
real-only query does not reduce it.

| C termination | C stop reason | C++ `Stop` | `complete` | Regions |
|---|---|---|---|---|
| `CPS_SUCCESS` | `CPS_STOP_NONE` | `none` (`ok()`) | 1 | the answer |
| `CPS_PRECISION_LIMIT` | `CPS_STOP_PRECISION` | `precision_cap` | 0 or 1 | certified, possibly partial |
| `CPS_PRECISION_LIMIT` | `CPS_STOP_INPUT`: the coefficients' own uncertainty blocks every pending region | `input_uncertainty` | 0 or 1 | certified, possibly partial |
| `CPS_WORK_LIMIT` | `CPS_STOP_WORK_BUDGET`, `CPS_STOP_CELL_BUDGET`, `CPS_STOP_ITERATION_BUDGET` | `work_budget`, `cell_budget`, `iteration_budget` | 0 or 1 | certified, possibly partial |
| `CPS_CANCELLED` | `CPS_STOP_CANCELLED` (the token) or `CPS_STOP_TIME_LIMIT` | `cancelled`, `time_limit` | 0 or 1 | certified, possibly partial |
| `CPS_UNSUPPORTED_QUERY` | `CPS_STOP_INPUT` | `unsupported` | 0 | none, or the regions certified before the refusal |
| `CPS_ZERO_POLYNOMIAL` | `CPS_STOP_INPUT` | `zero_polynomial` | 0 | none (degree −1) |
| `CPS_INDETERMINATE_DEGREE` | `CPS_STOP_INPUT` | `indeterminate_degree` | 0 | none |
| (a `cps_error` other than `CPS_OK`: no result) | | `error` | 0 | none |

`cps_termination_name` and `cps_stop_reason_name` give static English names ("success", "precision limit", ...;
"none", "work budget", "time limit", ...); in C++, `name(stop)` and `std::format("{}", stop)` ("precision cap",
"time limit", ...).

### Invariants

These hold for every result.

1. `goal_met` implies `complete`.
2. The termination is `CPS_SUCCESS` if and only if the goal is met. A cancellation or a time limit observed after
   the goal is met does not change it.
3. The stop reason is `CPS_STOP_NONE` if and only if the termination is `CPS_SUCCESS`.
4. Termination, stop reason, coverage and regions go together as the table above says.
5. Every region carries a proved count, whatever the termination.
6. The regions are pairwise disjoint as closed rectangles.
7. A complete all-complex answer has counts that sum to the degree. A complete real-only answer holds every real
   root of the requested domain, each in exactly one region.
8. A region's multiplicity is unknown (0) or equal to its count.
9. A region proved real is a real interval: both its imaginary ends are 0.

### Regions

```c
typedef struct cps_dyadic_view { const char *mantissa; const char *exponent; } cps_dyadic_view;
typedef struct cps_root_view {
    cps_dyadic_view real_lower, real_upper, imag_lower, imag_upper;
    uint64_t count;
    uint64_t multiplicity; /* 0 = not established */
    cps_realness realness;
    const char *certificate;
} cps_root_view;
```

| C++ `Root` | C `cps_root_view` | |
|---|---|---|
| `box.real.lo`, `box.real.hi`, `box.imag.lo`, `box.imag.hi` (`Dyadic`) | `real_lower`, `real_upper`, `imag_lower`, `imag_upper` | the closed rectangle; each end exactly mantissa · 2^exponent |
| `count` | `count` | the roots in the region, counted with multiplicity; proved |
| `multiplicity` (`std::optional`, empty when unknown) | `multiplicity` (0 when unknown) | when established, equal to `count`: the region holds one distinct root of that multiplicity |
| `realness`: `Realness::real`, `not_real`, `unknown` | `realness`: `CPS_REAL`, `CPS_NONREAL`, `CPS_REAL_UNKNOWN` | every root in the region is proved real (the region is a real interval), or the region misses the real axis, or neither |
| `value` (`std::optional<std::complex<double>>`) | `cps_result_get_approximation(r, i, &re, &im)` | the root as a double, the C call's number, and only where the C call gives one (see Approximations below); empty otherwise |
| `centre()` | `cps_dyadic_midpoint` of each part, then `cps_dyadic_to_double` | the box's centre rounded to doubles, for any region: what a `Root` prints |
| `certificate` | `certificate` | a short label naming how the count was proved |

- `cps_result_get_root` gives `CPS_INVALID_ARGUMENT` for NULL and `CPS_OUT_OF_RANGE` for an index at or past
  `cps_result_size`. The view's strings are borrowed until the result is destroyed; a C++ `Root` owns its data and
  outlives the result.
- The regions come in increasing order of their exact midpoint's real part, then its imaginary part, on every
  route.
- A count is never a claimed multiplicity: a region of count 3 with an unknown multiplicity may hold three distinct
  roots.
- Every region of a real-only answer is proved real.
- A program's region is `CPS_REAL` only by conjugation, when every constant of the program is real (section 5,
  "Solving a program").

### Exact ends

| C++ (`Dyadic`, `Interval`) | C | |
|---|---|---|
| `a <=> b`, `a == b`; with integers and doubles too (`root.box.real.hi < 0`) | `cps_dyadic_compare(a, b, &order)`: −1, 0 or 1 | exact comparison |
| `d.to_double(Rounding::down)`, `nearest` (the default), `up` | `cps_dyadic_to_double(d, CPS_ROUND_DOWN / _NEAREST / _UP, &out)` | the binary64 number rounded toward −∞, to nearest or toward +∞, correctly also where the exponent range overflows or underflows |
| `d.to_string()` | `cps_dyadic_rational(d, buffer, size, &length)` | the exact value, "p/q" or "p" in lowest terms |
| `d.format(n)`; `std::format("{:.N}", d)` (`{}`: 17 digits) | `cps_dyadic_format(d, n, buffer, size, &length)` | decimal with at most n significant digits |
| `interval.mid()`, `interval.width()` | `cps_dyadic_midpoint(a, b, ...)`, `cps_dyadic_width(a, b, ...)` | (a + b)/2 and b − a, exact; C writes "m@e" |
| `interval.contains(x)`, `interval.is_point()`; `std::format("{}", interval)` gives "[lo, hi]" | | |
| `certpolysolve::to_mpq(d)` (`<certpolysolve/gmp.hpp>`) | | the exact value as `mpq_class` |

- `cps_dyadic_format` writes as C's `"%.*g"` writes a double with n (≥ 1) as its precision, and as C++'s
  `std::format("{:.N}")` does: at most n significant digits, correctly rounded to nearest with ties to even,
  trailing zeros dropped, the exponent form (`7e-30`) below 10⁻⁴ and from 10ⁿ. The decimal point is `.` whatever the
  host's locale. A value beyond MPFR's exponent range (about 2^(±2^30)) gives `CPS_OUT_OF_RANGE`.
- A malformed or NULL end gives `CPS_INVALID_ARGUMENT`. The text calls follow the buffer protocol of section 3.
- A `Root` prints as its centre: `{}` to double precision, `{:.N}` with N digits from the exact box.

### Approximations

For the roots as plain numbers, convert each region's centre to a double. The roots of x² − 2 as doubles, as the
README's first program prints them in C++ (`examples/guide/plain_numbers.c`):

```c
#include <certpolysolve/certpolysolve.h>
#include <stdio.h>

int main(void) {
    cps_polynomial *p = NULL;
    cps_result *r = NULL;
    int status = 1;
    if (cps_polynomial_create(&p) != CPS_OK ||
        cps_polynomial_set_coefficient(p, 0, "-2", NULL) != CPS_OK || /* constant term first: x^2 - 2 */
        cps_polynomial_set_coefficient(p, 2, "1", NULL) != CPS_OK)
        goto done;
    if (cps_solve(p, NULL, &r) != CPS_OK) { /* NULL options: the defaults */
        fprintf(stderr, "%s\n", cps_last_error());
        goto done;
    }
    if (!cps_result_goal_met(r)) {
        fprintf(stderr, "%s: %s\n", cps_termination_name(cps_result_termination(r)), cps_result_diagnostics(r));
        goto done;
    }
    for (size_t i = 0; i < cps_result_size(r); ++i) {
        cps_root_view v;
        double re, im;
        if (cps_result_get_root(r, i, &v) != CPS_OK || cps_result_get_approximation(r, i, &re, &im) != CPS_OK)
            goto done;
        for (uint64_t k = 0; k < v.count; ++k) { /* a double root prints twice */
            if (im == 0)
                printf("%.17g\n", re);
            else
                printf("%.17g %+.17gi\n", re, im);
        }
    }
    status = 0;
done:
    cps_result_destroy(r);
    cps_polynomial_destroy(p);
    return status;
}
```

```c
/* Declared in <certpolysolve/certpolysolve.h>. */
cps_error cps_result_get_approximation(const cps_result *, size_t index, double *real, double *imag);
cps_error cps_dyadic_midpoint(cps_dyadic_view a, cps_dyadic_view b, char *buffer, size_t size, size_t *length);
cps_error cps_dyadic_to_double(cps_dyadic_view x, cps_rounding rounding, double *out);
```

`cps_result_get_approximation(r, i, &re, &im)` converts region i's midpoint to binary64, rounded to nearest. It is
an explicit, lossy conversion; the requested precision applies to the exact region, not to it. It fails, leaving
the outputs unchanged:

| Case | Code |
|---|---|
| NULL `r`, `re` or `im` | `CPS_INVALID_ARGUMENT` |
| `i` ≥ `cps_result_size(r)` | `CPS_OUT_OF_RANGE` |
| the result is incomplete or its goal is unmet; the region is a cluster (multiplicity ≠ count); the goal is `CPS_ROOT_ISOLATION`, which promises no width | `CPS_UNRESOLVED` |
| a part overflows binary64, or a nonzero part would round to zero | `CPS_OUT_OF_RANGE` |

A nonzero subnormal result is allowed, and may lose relative accuracy.

C and C++ give the same numbers. A C++ `Root`'s `value` is the C call's number for its region, bit for bit, and it
is empty where the C call refuses: after a stop, for a cluster whose multiplicity is not proved, for an isolation
answer, and outside binary64's range. `values()` lists the values, so it is empty after a stop.

A region's centre as a double, for any region and goal, is `Root::centre()` in C++. In C, write each part's exact
centre with `cps_dyadic_midpoint`, which gives "m@e"; split that text at the `@` into a `cps_dyadic_view` and convert
it with `cps_dyadic_to_double`. A centre is not a root's value unless the region has one: it may be a cluster's
centre, and it may overflow to infinity or underflow to zero.

### Facts

| C++ | C | |
|---|---|---|
| `sufficient_work_bits()` | `cps_result_sufficient_work_bits(r)` | 0 unless the solve stopped at its precision cap (`CPS_PRECISION_LIMIT` with `CPS_STOP_PRECISION`) on exact real input, whose distinct factors are then integer polynomials; then 2 · (separation bits + size bits + requested bits) from Mahler's root-separation bound: a cap at which the roots can be isolated at the requested precision. A bound, not a promise that a solve under it meets its goal within its other budgets. |
| `squarefree()`: `Squarefree::proved`, `disproved`, `unknown` | `cps_result_squarefree(r)`: `CPS_SQUAREFREE_PROVED`, `_DISPROVED`, `_UNKNOWN` | proved: every root is simple; disproved: a root of multiplicity above one is proved; unknown: neither (also where exact preparation was limited, and on ball input). A program's verdict comes from its image modulo a word-size prime, so it is proved or unknown. |
| `handoff()`: `Handoff::none`, `regime_exit`, `storage`, `range`, `precision` | `cps_result_handoff(r)`: `CPS_HANDOFF_NONE`, `_REGIME_EXIT`, `_STORAGE`, `_RANGE`, `_PRECISION` | the piecewise method, forced or chosen by Auto, handed the solve to the secular method within the call, and why: it spent more than a global solve's predicted cost, its storage was refused, the coefficients' exponents exceed its working range, or it reached its precision budget. The answer then comes from the secular method. |
| `preparation_limited()` | `cps_result_preparation_limited(r)` | exact preparation was refused or limited by the resource limits, so multiplicities may be unknown |

The diagnostics may repeat these facts; the getters are the reference.

### Diagnostics and certificates

- `cps_result_diagnostics` (C++ `message()`) is text for people: the preparation and adaptation trace. Its format
  is not part of the contract and may change in any release. It does not depend on the C or C++ locale. It is
  borrowed until the result is destroyed.
- `certificate` is a short label for people naming how a region's count was proved, for example
  `pellet-local-count`, `pellet-graeffe-local-count` or `cauchy-global-count`. The vocabulary is not fixed, and it
  is not a replayable proof object.

## 9. Statistics

Statistics are observations, not certificates; a field may be zero for a stage that did not run.

**What the call spent** is a promise: the fields have exactly these meanings.

| `cps_statistics` field | |
|---|---|
| `size` | in: `sizeof(cps_statistics)`; out: the bytes filled |
| `work_units` | controller work units the call spent, every engine and internal solve included |
| `sweeps` | sweeps of the simultaneous iterations: what the iteration budget counts |
| `precision_increases` | raises of the outer working precision |
| `final_work_bits` | the last outer working precision (`cps_result_work_bits`) |
| `highest_work_bits` | the highest outer working precision reached |
| `regions` | exported regions (`cps_result_size`) |
| `roots` | the sum of the regions' counts |
| `resolved_regions` | regions whose multiplicity is established |
| `exact_factors` | factors produced by exact preparation |
| `total_ns` | wall time of the call, on the solve's monotonic clock |
| `preparation_ns` | wall time of exact preparation, every engine included |

**The engines' counters** (`cps_engine_counters`) are diagnostics: fields may be added, removed or renamed in any
release. They hold stage times in nanoseconds (`input_enclosure_ns` to `regeneration_ns`; the stages overlap),
counts (precision attempts, isolation passes, disk tests, Graeffe steps, Taylor shifts and cache hits, squares
split and excluded, Newton attempts and successes, backend calls, evaluations and their terms, root interactions,
validations, secular regenerations, binary64 seed calls, candidate stagnations, real refinement steps and
bisections) and the search's state at the end (pending components, live cells, the pending precision range, input
uncertainty).

**The search's components** are diagnostics too: `cps_result_component_count(r)` and
`cps_result_get_component(r, i, &component)`, a `cps_component_progress` {`factor_index`, `work_bits`,
`live_cells`, `certified_count`, `state`}, the state one of `CPS_COMPONENT_ACTIVE`, `_FINISHED`,
`_PRECISION_LIMIT`, `_CELL_LIMIT`. An index names an entry of this result only, not a region and not a
persistent handle; an index at or past the count gives `CPS_OUT_OF_RANGE`.

| C++ | C |
|---|---|
| `result.statistics()` (a `cps_statistics`) | `cps_result_get_statistics(r, &statistics)` |
| `result.engine_counters()` (a `cps_engine_counters`) | `cps_result_get_engine_counters(r, &counters)` |
| `result.components()` (a `std::vector<cps_component_progress>`) | `cps_result_component_count`, `cps_result_get_component` |

- In C, set `size` to the struct's `sizeof` before the call. A smaller size gives `CPS_INVALID_ARGUMENT` and leaves
  the output unchanged; otherwise the library fills its layout, leaves the bytes beyond it untouched and writes
  the bytes it filled to `size`. A NULL result or output gives `CPS_INVALID_ARGUMENT`.
- The C++ calls set the size themselves; for an error result they return zeros and an empty vector.
- Times are on the solve's monotonic clock. They exclude the caller's own work before and after the call, and
  stage times need not sum to the call's time.

```c
cps_statistics statistics;
statistics.size = sizeof statistics;
if (cps_result_get_statistics(r, &statistics) == CPS_OK)
    printf("%llu work units in %llu ns\n", (unsigned long long)statistics.work_units,
           (unsigned long long)statistics.total_ns);
```

## 10. Threads, lifetimes and FLINT

### Guarantees

- **Concurrent solves.** Any number of threads may call the library at once, each on objects it does not share.
- **Read-only sharing.** An input, an options object or a result may be read by any number of threads at once,
  while no thread changes or destroys it. Solving with an input or options is reading them.
- **Moving between threads.** An object may be made on one thread and used or destroyed on another, also after
  the making thread has called `cps_cleanup_thread`. This relies on FLINT's thread-aware memory manager, so the
  build requires a FLINT built with thread-local storage and thread support (`FLINT_USES_TLS`,
  `FLINT_USES_PTHREAD`); it fails otherwise.
- **Tokens.** `cps_cancel_token_request` may be called on any thread at any time, also while solves read the token.
- **Changes are exclusive.** A setter, a builder or a destroy on an object must not run while another thread uses
  that object: that is a data race.
- **No hidden activity.** No caller code runs inside a solve. The library starts no threads, except that a solve
  whose options ask for more than one worker starts up to workers - 1 threads and joins them before it returns: no
  thread outlives the call. It keeps no mutable global state apart from each thread's error message, reads no
  environment variable, and writes nothing to `stdout` or `stderr`.
- **The same answer at every worker count.** A solve's answer, counters and diagnostics do not depend on the
  worker count. A stop by a time limit or a token depends on the machine, as it does with one worker.

### Per-thread state

- `cps_last_error` is per thread (section 2).
- FLINT keeps thread-local caches (constants, prime tables, FFT contexts, free integer lists), which grow with the
  precision a thread has used.
- A solve saves the calling thread's floating-point environment at entry, runs with round-to-nearest, every
  exception masked, and flush-to-zero and denormals-are-zero off, and restores the caller's environment, its
  raised flags included, before it returns.
- A call may change `errno`.
- A solve's workers run in the default floating-point environment, start with every signal blocked, never use
  FLINT's thread pool, and free their FLINT caches before they end. The calling thread's caches and environment
  are as without workers.

### FLINT, GMP and the host

- FLINT 3.6 or later is required, built with thread-local storage and thread support. `cps_backend_version()`
  returns FLINT's version string.
- An error inside FLINT or GMP, an allocation failure inside them included, calls `abort()`: it ends the process
  and every solve on every thread. The resource limits bound FLINT's allocations only as estimates.
- The host must not install GMP or FLINT memory functions that are not thread-safe, and must not change them
  while a solve runs.
- The library never calls `flint_set_num_threads`. A host that calls it while a solve uses FLINT's thread pool
  makes FLINT abort. A host thread that raised its FLINT worker count lets FLINT routines inside its solves use the
  shared pool.
- `cps_cleanup_thread()` frees the calling thread's FLINT caches, and MPFR's (which MPFR may share between threads).
  Call it before a thread that used the library exits, at a point where no other FLINT user on that thread holds
  pointers into those caches; a short-lived thread that skips it leaks them. It is never needed for correctness,
  and objects made on the thread stay valid after it. In C++, construct a `certpolysolve::ThreadScope` first at the
  thread's entry: its destructor calls it. Do not nest scopes. The main thread needs none.
- A solve runs on the caller's stack, and FLINT and GMP take temporaries from it. A thread that calls the library
  needs at least 512 KiB of stack: over the 6,894 solves of the gate corpus, every method and query, the deepest
  took 169 KiB, the thread's own 16 KiB included.
- Above 4,592 bits, FLINT's thread-local π and log 2 are rounded from a cache whose precision depends on earlier
  work on the same thread, so an exported end can differ in near-tie cases. Certificates are unaffected. A test of
  bit-identical output should run each solve on a cleaned thread.

### Lifetimes

| Object or pointer | Valid until |
|---|---|
| A handle (`cps_polynomial`, `cps_secular`, `cps_program`, `cps_options`, `cps_result`) | its destroy call; it may be used on any thread |
| A token | the last of: the caller's destroy, and every options object and solve that holds it |
| Strings and views from a result | the result's destroy |
| A string from an options getter | the next change to those options, or their destroy |
| `cps_last_error()` | the calling thread's next `cps_error`-returning call, on that thread only |
| `cps_version()`, `cps_backend_version()`, the name functions' strings | forever (static) |
| A `ProgramValue` | the end of the call of f that made it |
| A C++ `Root`, `Polynomial`, `Number`, statistics struct | its own lifetime: it owns its data |

The library never frees caller memory, and the caller never frees anything the library returns except through the
destroy functions.

## 11. Versions

| C | C++ | |
|---|---|---|
| `CPS_VERSION_MAJOR` 0, `CPS_VERSION_MINOR` 5, `CPS_VERSION_PATCH` 1, `CPS_VERSION_STRING` "0.5.1" | the same macros | the version of the headers |
| `cps_version()` | `certpolysolve::version()` | the library's version string: `CPS_VERSION_STRING` of the headers it was built with |
| `CPS_ABI_VERSION` (2), `cps_abi_version()` | `certpolysolve::abi_version()` | the binary interface's number, raised whenever it changes; informational |
| `cps_backend_version()` | `certpolysolve::backend_version()` | FLINT's version string |

Before 1.0.0 no compatibility of any kind is promised: any release may change the API, the ABI and the behaviour.
A program may compare `cps_abi_version()` with `CPS_ABI_VERSION`, and `cps_version()` with `CPS_VERSION_STRING`, to
detect a mismatched library; the library does not check.

## 12. The command-line tool

`certpolysolve_cli` (built from `examples/solve.cpp`) solves one polynomial or one secular equation:

```sh
certpolysolve_cli [--real] [--approximation | --clusters R] [--bits N] [--method auto|classical|secular|piecewise] [--highest-first] a0 a1 ... an
certpolysolve_cli [--real] [--approximation | --clusters R] [--bits N] [--method ...] --secular FILE
```

- Coefficients are exact, in the literal grammar, with the constant term first (`-2 0 1` is x² − 2); after
  `--highest-first` (or `--descending`) the highest power comes first (`1 0 -2`). The tool echoes the polynomial it
  read on stderr (`p(x) = x^2 - 2`), so a reversed list shows. `--bits N` (1 to 1,000,000) sets the requested
  precision and a cap of max(4,096, N + 64), at most 1,000,000. `--real` asks for the real roots,
  `--approximation` for the approximation goal, `--clusters R` for boxes of radius at most R, each with its count.
- It prints the degree, `complete`, `goal_met`, the termination and the last working precision, then each region
  with its exact ends `mantissa * 2^(exponent)`, count, multiplicity and certificate, then the diagnostics. It exits
  with 0 when the goal is met, 2 for a stop, 1 for an error.
- A secular file uses MPSolve's modern secular subset and needs an explicit `Secular; Degree=N;` declaration. It
  accepts `Dense`, `Real` or `Complex`, and `Integer`, `Rational` or `Floatingpoint` literals, read exactly.
  Comments start with `!`; trailing data and conflicting headers are rejected; an input precision must be omitted
  or zero. The data are the terms in order, each a residue and then a pole (with `Complex`, each as a real and an
  imaginary part):

```text
Secular;
Rational;
Real;
Degree=2;
2/3 0
-2/3 3
```

## 13. Changes in 0.5

The pre-0.5 interface (version 0.1.0-dev, ABI 1) changed as follows. The design is recorded in
[the 0.5 specification](api-spec-0.5.md), its decisions and
the C++ proposal.

| Removed or renamed | Replacement |
|---|---|
| `cps_cancel_callback`, `cps_options_set_cancel` | `cps_options_set_time_limit`; `cps_cancel_token_*` with `cps_options_set_cancel_token` (C++: `Options::stop`, a `std::stop_token`). No caller code runs inside a solve. |
| `cps_options_set_precision(o, bits, cap)` | `cps_options_set_precision(o, bits)` and `cps_options_set_max_work_bits(o, bits)` |
| `cps_options_set_search_budget(o, units, cells)` (0 units: no work) | `cps_options_set_work_budget` and `cps_options_set_live_cells` (0: the default) |
| `cps_options_set_exact_budget(o, bits)` | `exact_materialization_bits` in the resource limits, with `CPS_RESOURCE_ALLOW_EXPONENT_EXPANSION` if wanted |
| `cps_statistics_v1`, `_v2`, `_v3`, `cps_progress_v1` and their getters | `cps_statistics` (`cps_result_get_statistics`), `cps_engine_counters` (`cps_result_get_engine_counters`), `cps_result_stop_reason`, `cps_result_preparation_limited` |
| `cps_component_progress_v1`, `cps_result_get_component_progress_v1` | `cps_component_progress`, `cps_result_get_component` |
| `CPS_COMPONENT_STAGNATION` (3), which no solve produced | none |
| `cps_resource_limits_v1`, `cps_resource_limits_default_v1`, `cps_options_set_resource_limits_v1`, `cps_options_get_resource_limits_v1`, `cps_options_set_region_v1` | the same names without `_v1`, and `cps_options_get_region` |
| `cps_secular_create(terms, count, out)`, `cps_secular_create_with_resources_v1` | `cps_secular_create(terms, count, limits or NULL, out)` |
| `cps_program_output` returning `int64_t`, −1 when unset | `uint64_t`, `CPS_NO_NODE` when unset; a failed builder writes `CPS_NO_NODE` |
| sessions: `cps_session_create`, `_create_secular`, `_create_program`, `_advance`, `_destroy`; C++ `SolveSession` | one-shot solves; a different budget, goal or cap is a new solve |
| the pre-0.5 C++ header (`legacy.hpp` after stage 3a, not installed, deleted in stage 3b-3): `RootResult`, its exceptions, `enclosures()`, `approximations()`, `progress()`, `statistics_v2/v3()`; `roots(input)`, `clusters(input, radius)`; `RootEnclosure`; `SecularTerm`; `ProgramConstant`; `CancelToken`; chained `Options` setters | `Result` (`ok()`, `stop()`, `message()`; nothing thrown but `std::bad_alloc`), `roots()` of `Root` values with exact `Dyadic` ends and `value`, `statistics()`, `engine_counters()`; `solve(input, options)` and `real_roots`; `SecularEquation::Term`; plain numbers; `std::stop_token`; `Options` fields |
| `cps_version()` "0.1.0-dev", `CPS_ABI_VERSION` 1 | "0.5.0", 2, and the version macros |

Behaviour that changed:
- A query no input of the entry point's kind can answer is `CPS_INVALID_ARGUMENT` from the solve; it was `CPS_OK`
  with `CPS_UNSUPPORTED_QUERY`. The region setter checks literals only; the ends' order and exact form are checked
  at the solve.
- The options' initial goal is `CPS_GOAL_DEFAULT`, resolved per input kind, and a zero cluster radius restores it.
- A NULL imaginary part or radius means 0 in every input function.
- The regions come in exact midpoint order on every route.
- A solve whose goal is met reports `CPS_SUCCESS` even when a cancellation is seen afterwards, and
  `CPS_ZERO_POLYNOMIAL` carries `CPS_STOP_INPUT`.
- A call that succeeds clears `cps_last_error()`.
- A solve runs in the default floating-point environment and restores the caller's.
- Objects may move between threads and outlive their thread's cleanup.
- The coefficient object keeps its nonzero terms only, and exact polynomial arithmetic is in the library.
- The default precision cap is 65,536 bits.
