# CertPolySolve 0.5 API specification

The contract of the 0.5 API: the inputs, the options, the goals, the results, the errors, the lifetimes and the
thread guarantees. The C++ interface uses the C API and supplies owning values, options and results. [The guide](api.md) explains both by use;
[the glossary](glossary.md) defines the words. Where this text and the headers disagree, the headers and the guide
describe the library. The normative text uses the present tense: "returns", "must", "may".

---

## 1. Overview

### 1.1 Objects

| Object | C type | Made by | Freed by | Changed by |
|---|---|---|---|---|
| Coefficient polynomial | `cps_polynomial` | `cps_polynomial_create` | `cps_polynomial_destroy` | its setters |
| Native secular equation | `cps_secular` | `cps_secular_create` | `cps_secular_destroy` | nothing (immutable) |
| Program (straight-line) | `cps_program` | `cps_program_create`, `cps_program_parse` | `cps_program_destroy` | its builders (C only) |
| Options | `cps_options` | `cps_options_create`, `cps_options_copy` | `cps_options_destroy` | its setters |
| Cancellation token (advanced) | `cps_cancel_token` | `cps_cancel_token_create` | `cps_cancel_token_destroy` | `cps_cancel_token_request` |
| Result | `cps_result` | a solve | `cps_result_destroy` | nothing (immutable) |

All C handle types are opaque. The C++ input and result classes manage their handles automatically;
copies of immutable secular equations, programs and results share ownership.

### 1.2 Lifetimes

- The caller owns every object it makes and frees it with the matching destroy function. A C++ object does this in its destructor.
- The input and options must remain alive and unchanged throughout a solve. Once the solve returns, the caller may change or destroy them; the result does not depend on them.
- A result owns everything it exposes. Every string and view it hands out stays valid until the result is destroyed.
- A token is shared. Options that name a token, and solves that use it, each keep it alive. `cps_cancel_token_destroy` releases only the caller's own reference.
- The library never frees caller memory. Apart from the destroy functions, the caller never frees anything the library returns.

### 1.3 Call sequence

1. Make an input: a polynomial, a secular equation or a program.
2. Optionally, make options and set what differs from the defaults (§4).
3. Solve. The call returns a result, or an error code if the call itself was wrong (§2.1).
4. Read the outcome (`termination`, `complete`, `goal_met`) and the regions (§6).
5. Destroy the result, the options and the input, in any order.
6. Before a thread that used the library exits, call `cps_cleanup_thread()`, or let a `ThreadScope` do it (§8).

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

The same solve in C++ (C++23 for `std::println`):

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

`CPS_OK` means the C call ran and returned a result. Inspect `goal_met` and `termination` to learn whether
it reached the mathematical goal (§6.2). C++ `Result::ok()` means the goal was met.

### 1.4 What a simple caller needs

An input, `solve`, and the result's `goal_met`, `size` and regions. Options are needed only to change the
precision or to set a time limit. Everything else is advanced, and a caller may ignore it: the
cancellation token, the work, cell and iteration budgets, the resource limits, the statistics, the engine
counters and the diagnostics.

---

## 2. Conventions

### 2.1 Errors in C

- Every C function that can fail returns `cps_error`. `CPS_OK` (0) means success.
- A function that returns any other type cannot fail. Passing it NULL is a usage error, which it tolerates by returning the fallback value listed in §2.3.
- Nothing crosses the C boundary: no C++ exception, no `longjmp`. No C function calls user code.
- A call that fails changes nothing. Its outputs are:
  - a handle output (`T **out`) is set to NULL;
  - a node output is set to `CPS_NO_NODE`;
  - every other output is left as it was;
  - the object operated on is unchanged, because setters and builders are transactional.
- A mathematical non-success is not an error. A precision stop, a spent budget, a cancellation or a refused input all return `CPS_OK` with a result that says what happened (§6.2).
- Usage errors and refusals are told apart **per entry point**:
  - An option combination that no input of the entry point's kind can satisfy is a usage error. The solve returns `CPS_INVALID_ARGUMENT` and no result.
  - A refusal that depends on the input is an outcome: `CPS_OK` and a result with termination `CPS_UNSUPPORTED_QUERY` (or one of the other input terminations).
  - §5.2 lists every case.

| Code | Meaning |
|---|---|
| `CPS_OK = 0` | The call succeeded. |
| `CPS_INVALID_ARGUMENT = 1` | The call breaks a stated precondition. Examples: a NULL required pointer; a malformed literal; a value outside its range; an operand that is not an earlier node; a program without an output at the solve; an option combination that §5.2 marks as a usage error. |
| `CPS_OUT_OF_MEMORY = 2` | A C++ allocation inside the library failed. Nothing was published. A failure inside FLINT or GMP aborts the process instead (§8.3). |
| `CPS_INTERNAL_ERROR = 3` | Anything else failed inside the library. No partial result is published. This is a library bug unless the message says otherwise. |
| `CPS_OUT_OF_RANGE = 4` | The request is well formed but beyond what the library can represent or admit. Exactly these cases: an index at or past the end of a sequence (a root or a component); a coefficient index beyond the backend's limit or addressable storage; a precision above the backend's maximum; a secular input over its resource limits at creation (§3.2); a binary64 conversion that would overflow or flush a nonzero value to zero (§6.6). |
| `CPS_UNRESOLVED = 5` | A binary64 approximation was asked of a region the result does not resolve (§6.6). |


`const char *cps_last_error(void)`:
- Every function that returns `cps_error` sets the calling thread's message: "" when it returns `CPS_OK`, a description when it fails. Functions of other return types leave the message unchanged.
- The pointer is valid on the calling thread only, until that thread's next `cps_error`-returning call.
- The text is informational. Its format is not promised.
- In C++ a result's `message()` and an input's `status()` carry the text (§2.2); `cps_last_error` serves callers that mix C calls into C++ code.

### 2.2 Errors in C++

The C++ header is header-only over the C API and throws nothing but `std::bad_alloc`, when memory runs out. No
call throws for a usage error, a library error or a mathematical outcome:
- A result carries its status: `ok()` is true when the solve ran and met its goal, `stop()` says why a solve
  stopped (`Stop::error` for a usage or library error, with nothing solved), and `message()` holds the error's
  text or the solver's diagnostics.
- An input (a polynomial, a secular equation or a program) remembers its first bad value: `ok()` is false and
  `status()` says why. Options are plain fields, validated when solving. A bad input or invalid options produce
  a result whose `stop()` is `Stop::error` with the reason in `message()`.
- `CPS_INTERNAL_ERROR` is reported the same way. It is not thrown.
- An exception thrown by the caller's own code inside `Program::from(f)` propagates unchanged, and the recording is
  discarded (§3.3).

[The guide](api.md) describes the C++ types and their C counterparts. The tables below give the C signatures.

### 2.3 NULL rules

| Where NULL appears | Meaning |
|---|---|
| A handle argument of a `cps_error` function | `CPS_INVALID_ARGUMENT` |
| A handle argument of a function of another return type | That function's fallback value (table below) |
| The argument of a destroy function, or of `cps_cancel_token_request` | No effect |
| An output pointer | `CPS_INVALID_ARGUMENT` |
| The options of a solve | The defaults |
| The resource limits of `cps_secular_create` | The defaults |
| A real part, a pole's real part, a program text, a tolerance or a cluster radius | `CPS_INVALID_ARGUMENT` |
| An imaginary part or a radius, in every input function | The value 0 |
| A region end | That side is unbounded. The empty string means the same. |
| The token of `cps_options_set_cancel_token` | The options' token is cleared |

Fallback values for a NULL handle:

| Getter | Fallback |
|---|---|
| `cps_result_degree` | −1 |
| `cps_result_termination` | `CPS_UNSUPPORTED_QUERY` |
| `cps_result_stop_reason` | `CPS_STOP_INPUT` |
| `cps_result_complete`, `cps_result_goal_met`, `cps_result_preparation_limited` | 0 |
| `cps_result_goal` | `CPS_ROOT_ACCURACY` |
| `cps_result_work_bits`, `cps_result_size`, `cps_result_component_count`, `cps_result_sufficient_work_bits` | 0 |
| `cps_result_diagnostics` | the static string "null result" |
| `cps_result_squarefree` | `CPS_SQUAREFREE_UNKNOWN` |
| `cps_result_handoff` | `CPS_HANDOFF_NONE` |
| `cps_program_size` | 0 |
| `cps_program_output` | `CPS_NO_NODE` |
| options getters | the default value of that setting |
| `cps_cancel_token_requested` | 0 |
| `cps_resource_limits_default` | does nothing |

The rule for C++ objects:
- A moved-from C++ object may only be destroyed or assigned to.
- Using it in any other way is a usage error. Its getters return the fallback values above, and its fallible calls return `CPS_INVALID_ARGUMENT`.

### 2.4 Number literals

Every number the caller passes as text uses one grammar: coefficient parts and radii, secular terms,
program constants, the absolute tolerance, the cluster radius and the region ends. The same grammar applies to all of them.

```
literal   = space* ( rational | dyadic | decimal ) space*
space     = " " | TAB | CR | LF                       (only around the literal, never inside it)
rational  = integer "/" integer                       value p/q; q must not be 0
dyadic    = integer "@" integer                       value m * 2^e
decimal   = sign? mantissa ( ("e" | "E") integer )?   value m * 10^e
mantissa  = digit+ ( "." digit* )? | "." digit+
integer   = sign? digit+
sign      = "+" | "-"
```

- Values are exact. Exponents are integers of any size and stay compact: `1@100000000000000000000` is a few bytes. Whether a value is expanded is governed by the resource limits (§4.4).
- Rejected: empty text, whitespace inside a literal, hexadecimal, `inf`, `nan`, expressions, and a second operator (for example `1@2@3`).
- Radii, the tolerance and the cluster radius must also be ≥ 0.
- A double passed as a double (`cps_polynomial_set_coefficient_d`, the C++ double forms) enters at its exact binary value. NaN and ±infinity give `CPS_INVALID_ARGUMENT`.
- A literal longer than the backend can address gives `CPS_OUT_OF_RANGE`.

### 2.5 Units

| Quantity | Unit |
|---|---|
| Requested precision | Bits of relative accuracy: the region's half-diagonal is at most 2^−bits times sqrt(min abs(Re z)^2 + min abs(Im z)^2) over the region, a lower bound for the modulus of every root in it. So every point of the region is within 2^−bits · abs(z) of the root z it holds. An absolute tolerance, when set, is an alternative bound on the half-diagonal. |
| Working precision, the precision cap, `work_bits` | Bits of the ball arithmetic's binary precision |
| Absolute tolerance, cluster radius | An absolute distance in the complex plane, in the roots' own units |
| Time limit | Seconds (a `double`), on a monotonic clock |
| Durations in statistics | Nanoseconds of wall-clock time on a monotonic clock |
| Work budget | Controller work units (§4.4) |
| Iteration budget | Sweeps of the simultaneous iterations (§4.4) |
| Live cells | Retained search squares (§4.4) |
| `exact_materialization_bits` / `exact_operations` / `temporary_bytes` | Bits / word operations / bytes, each a per-stage estimate (§4.4) |
| Region endpoints | mantissa · 2^exponent, as two decimal integer strings |

### 2.6 Headers and language levels

- `include/certpolysolve/certpolysolve.h` is the C API. It compiles as C99 or later, and as C++ (it uses `extern "C"`).
- `include/certpolysolve/certpolysolve.hpp` is the C++ API, in namespace `certpolysolve`. It requires C++20 and is header-only on top of `certpolysolve.h`.
- Both headers are public. The C++ API uses the C enums, structs and macros directly, and they are part of the C++ contract.
- `CPS_API` marks exported functions. It means `dllexport` while the library is built (`CERTPOLYSOLVE_BUILDING_LIBRARY`), `dllimport` otherwise on Windows, and default visibility elsewhere. A static-library switch is not part of 0.5.

---

## 3. Inputs

### 3.1 Coefficient polynomial

What it means:
- p(x) = Σ aᵢ xⁱ, coefficients in ascending order.
- Each aᵢ is an exact complex number or a box: componentwise intervals [re ± r_re] + i·[im ± r_im].
- An unset coefficient is exactly zero.
- The degree is the highest index whose coefficient is not exactly zero. A box is never treated as zero. If the top box contains zero, the solve may end with `CPS_INDETERMINATE_DEGREE`.
- A certificate on box input holds for every polynomial whose coefficients lie in the boxes.
- There is no fixed limit on degree, literal length or precision. Representability and the resource limits are the only bounds.

| C signature | Contract |
|---|---|
| `cps_error cps_polynomial_create(cps_polynomial **out)` | Makes the zero polynomial. Errors: NULL `out` (`CPS_INVALID_ARGUMENT`); out of memory. |
| `void cps_polynomial_destroy(cps_polynomial *p)` | NULL does nothing. |
| `cps_error cps_polynomial_set_coefficient(cps_polynomial *p, uint64_t index, const char *real, const char *imag)` | Sets a_index = real + i·imag exactly, overwriting any earlier value. `real` must not be NULL. A NULL `imag` means 0. `CPS_INVALID_ARGUMENT`: NULL `p` or `real`, or a malformed literal. `CPS_OUT_OF_RANGE`: `index` > 2^63 − 1. Only nonzero terms are stored: a large index does not allocate all preceding coefficients. `CPS_OUT_OF_MEMORY` means a C++ allocation failed. Strings are read during the call only. |
| `cps_error cps_polynomial_set_coefficient_d(cps_polynomial *p, uint64_t index, double real, double imag)` | Sets the exact binary values. NaN or ±infinity gives `CPS_INVALID_ARGUMENT`. The index rules are as above. |
| `cps_error cps_polynomial_set_ball_coefficient(cps_polynomial *p, uint64_t index, const char *real, const char *imag, const char *real_radius, const char *imag_radius)` | Sets the box [real ± real_radius] + i·[imag ± imag_radius]. The radii are componentwise, not a disk. A negative radius gives `CPS_INVALID_ARGUMENT`. NULL `imag` or a NULL radius means 0. Otherwise as `set_coefficient`. |

### 3.2 Native secular equation

What it means:
- S(x) = Σⱼ aⱼ/(x − bⱼ) − 1, with residues aⱼ and poles bⱼ exact complex numbers.
- The poles are excluded from the domain. The roots are the zeros of S.
- The degree is the number of terms. The empty sum is the nonzero constant −1, which has no roots.

| C signature | Contract |
|---|---|
| `typedef struct cps_secular_term_view { const char *residue_real, *residue_imag, *pole_real, *pole_imag; } cps_secular_term_view;` | One term. The strings are literals (§2.4), read during the create call only. `residue_real` and `pole_real` must not be NULL. NULL imaginary parts mean 0. |
| `cps_error cps_secular_create(const cps_secular_term_view *terms, size_t count, const cps_resource_limits *limits, cps_secular **out)` | Copies the terms and normalizes them to exact rationals. `count` 0 with NULL `terms` is allowed. NULL `limits` means the defaults. **`CPS_INVALID_ARGUMENT`:** NULL `out`; NULL `terms` with `count` > 0; a NULL required literal; a malformed literal; a zero residue; two equal poles (compared exactly, so `1/2` equals `0.5`); unknown flag bits. **`CPS_OUT_OF_RANGE`:** `count` beyond addressable storage; the estimated storage over `temporary_bytes`; `count` > `exact_operations` / 4; exponent expansion beyond the input-relative rule without `CPS_RESOURCE_ALLOW_EXPONENT_EXPANSION`; normalization over min(`exact_materialization_bits`, 8 · (`temporary_bytes` − the term storage estimate)). |
| `void cps_secular_destroy(cps_secular *s)` | NULL does nothing. |

- Some queries need the exact numerator polynomial: real-only queries, and unresolved cases that fall back to it. For those, the solve applies the options' resource limits again, and a refusal there is an outcome (§5.2).

### 3.3 Program (straight-line program)

What it means:
- A program is a list of nodes, numbered 0, 1, 2, … in the order they are appended.
- Each node is one of: x; an exact complex constant; a + b; a − b; a · b; −a; or aᵉ with e ≥ 2. Operands a and b must be earlier nodes.
- One node is the output, and the polynomial is the output's value.
- Several x nodes all denote the same x.
- The program is never expanded. Its degree and leading coefficient are proved from the program itself, and its roots are bounded by a majorant.
- A failed builder call leaves the program unchanged.

`#define CPS_NO_NODE UINT64_MAX`

| C signature | Contract |
|---|---|
| `cps_error cps_program_create(cps_program **out)` | Makes an empty program with no output. NULL `out` gives `CPS_INVALID_ARGUMENT`. |
| `void cps_program_destroy(cps_program *p)` | NULL does nothing. A solve copies its program, so destroying the program does not affect results. |
| `cps_error cps_program_variable(cps_program *p, uint64_t *out)` | Appends an x node and writes its index to `*out`. |
| `cps_error cps_program_constant(cps_program *p, const char *real, const char *imag, uint64_t *out)` | Appends the exact constant real + i·imag. NULL `imag` means 0. NULL `real` or a malformed literal gives `CPS_INVALID_ARGUMENT`. |
| `cps_error cps_program_add(cps_program *p, uint64_t a, uint64_t b, uint64_t *out)`; the same for `_subtract` and `_multiply` | Appends a node. An operand that is not an earlier node gives `CPS_INVALID_ARGUMENT`. |
| `cps_error cps_program_negate(cps_program *p, uint64_t a, uint64_t *out)` | As above. |
| `cps_error cps_program_power(cps_program *p, uint64_t a, uint64_t exponent, uint64_t *out)` | Appends aᵉ. An exponent below 2 gives `CPS_INVALID_ARGUMENT`. In C++ only, `pow(a, 0)` is the constant 1 and `pow(a, 1)` is a itself; neither appends a power node. |
| `cps_error cps_program_set_output(cps_program *p, uint64_t node)` | Makes `node` the output, replacing any earlier one. A node at or past the size gives `CPS_INVALID_ARGUMENT`. |
| `cps_error cps_program_parse(const char *text, cps_program **out)` | Builds a program from the text form below. NULL text, a syntax error, a bad operand or a missing `out` line gives `CPS_INVALID_ARGUMENT`; the message names the line where there is one. |
| `uint64_t cps_program_size(const cps_program *p)` | The number of nodes. |
| `uint64_t cps_program_output(const cps_program *p)` | The output node, or `CPS_NO_NODE` when none is set. |

In every builder, `*out` receives the new node's index on success and `CPS_NO_NODE` on failure. A NULL `p` or `out` gives `CPS_INVALID_ARGUMENT`.

**Text form** (normative):
- Lines end with LF. A trailing CR is ignored. `#` starts a comment that runs to the end of the line.
- Fields are separated by spaces or tabs. Blank lines are ignored.
- The first non-blank line must be exactly `PSPROGRAM 1`.
- Each later line is one operation:
  - `x`;
  - `c REAL [IMAG]` (literals as in §2.4);
  - `add A B`, `sub A B`, `mul A B`;
  - `neg A`;
  - `pow A E` (E ≥ 2);
  - `out N`.
- A, B, N and E are unsigned decimal integers.
- Every operation except `out` appends one node, numbered in order. If there are several `out` lines, the last one wins. The text must contain at least one.

**The C++ recorder**:
- `Program::from(f)` calls f exactly once, with a `ProgramValue` for x. Its returned expression or plain number
  becomes the output.
- A `ProgramValue` supports `+`, `-`, `*`, unary `-`, `+=`, `-=`, `*=` and `pow(value, exponent)`. It mixes with
  the types accepted by `Number`: integers, exact binary values of doubles, text, complex numbers and types
  with a `number_traits` specialization. Each mixed operation appends a constant node.
- A `ProgramValue` has no comparisons or conversion to `bool`. Loops and branches over integers are fine.
- The first recording error is remembered: `ok()` is false and `status()` says why. A solve reports that error.
- f runs before any solve and may call library functions. Its exceptions propagate out of `from`.
- A `ProgramValue` is valid only during the call of f that produced it. Using one after `from` returns is undefined.
- `Program(text)` parses the text form. A default-constructed `Program` has no output, so solving it is an error.
  Copies share their program.

---

## 4. Options

### 4.1 Settings, defaults and ranges

| Setting | Default | Valid values | Used by |
|---|---|---|---|
| Method | `CPS_METHOD_AUTO` | `CPS_METHOD_AUTO`, `_CLASSICAL`, `_SECULAR`, `_PIECEWISE` | §5.2 |
| Requested precision | `CPS_DEFAULT_BITS` = 53 | 1 to 2^63 − 2 | The accuracy and approximation goals. The isolation goal ignores it. Under the cluster-radius goal it only steers the precision schedule. |
| Precision cap (`max_work_bits`) | `CPS_DEFAULT_MAX_WORK_BITS` = 65,536 | `CPS_MIN_MAX_WORK_BITS` (32) to 2^63 − 2. It may be below the requested precision. | Every solve |
| Goal | `CPS_GOAL_DEFAULT` | A `cps_goal` value | §6.4 |
| Cluster radius | none ("0") | A literal ≥ 0 | The `CPS_CLUSTER_RADIUS` goal |
| Absolute tolerance | none ("0") | A literal ≥ 0 | The accuracy and approximation goals |
| Real-only | off | on or off | §5.2 |
| Region | none (all four sides unbounded) | Four literals, each of which may be unbounded | Real-only queries under Auto on coefficient input |
| Work budget | 0: automatic, max(`CPS_DEFAULT_WORK_UNITS` = 100,000, `CPS_DEFAULT_WORK_UNITS_PER_ROOT` = 64 × degree) | Any `uint64_t`; 0 selects the automatic budget | Every solve |
| Live cells | 0, which means `CPS_DEFAULT_LIVE_CELLS` = 65,536 | Any `uint64_t`; 0 selects the default | Subdivision |
| Iteration budget | 0: no limit | Any `uint64_t` | The classical and secular iterations |
| Time limit | 0: no limit | 0 or +∞ for none; otherwise a finite positive number of seconds | Every solve |
| Cancellation token | none | A token, or none | Every solve |
| Resource limits | every allowance `UINT64_MAX` (unspecified), flags 0 | Any values. Flags may contain only `CPS_RESOURCE_ALLOW_EXPONENT_EXPANSION`. | Exact preparation and storage estimates |
| Workers | 1: no thread started | 1 or more; 0 is invalid. At most the machine's hardware threads are used. | Every solve: the threads it may use (§8.1, T6 and T8) |

"Degree" in the work budget means:
- for coefficients, the highest index whose coefficient is not exactly zero;
- for secular input, the number of terms;
- for a program, its degree bound (which saturates).

All the defaults are named macros in `certpolysolve.h`.

### 4.2 Functions

Setters validate only their own arguments (§4.3). A failing setter leaves the options unchanged. In C++ the options are a struct of fields, checked by the solve. Getters
that return a value cannot fail; given NULL they return the default. The two `get_*` functions that
return `cps_error` (`cps_options_get_region`, `cps_options_get_resource_limits`) give `CPS_INVALID_ARGUMENT`
for NULL.

| C signature | Contract |
|---|---|
| `cps_error cps_options_create(cps_options **out)` | Makes options with every default of §4.1. NULL `out` gives `CPS_INVALID_ARGUMENT`. |
| `cps_error cps_options_copy(const cps_options *source, cps_options **out)` | An independent copy of every setting. The copy names the same token as the source. |
| `void cps_options_destroy(cps_options *o)` | NULL does nothing. |
| `cps_error cps_options_set_method(cps_options *o, int method)` / `cps_method cps_options_method(const cps_options *o)` | A value that is not a `cps_method` gives `CPS_INVALID_ARGUMENT`. The integer parameter lets C reject bad values before an enum exists. |
| `cps_error cps_options_set_precision(cps_options *o, uint64_t bits)` / `uint64_t cps_options_precision(const cps_options *o)` | 0 gives `CPS_INVALID_ARGUMENT`. Above 2^63 − 2 gives `CPS_OUT_OF_RANGE`. |
| `cps_error cps_options_set_max_work_bits(cps_options *o, uint64_t bits)` / `uint64_t cps_options_max_work_bits(const cps_options *o)` | Below `CPS_MIN_MAX_WORK_BITS` gives `CPS_INVALID_ARGUMENT`. Above 2^63 − 2 gives `CPS_OUT_OF_RANGE`. |
| `cps_error cps_options_set_goal(cps_options *o, int goal)` / `cps_goal cps_options_goal(const cps_options *o)` | A value that is not a `cps_goal` gives `CPS_INVALID_ARGUMENT`. `CPS_CLUSTER_RADIUS` while no positive radius is set gives `CPS_INVALID_ARGUMENT`. Any other goal clears the radius. The getter may return `CPS_GOAL_DEFAULT`. |
| `cps_error cps_options_set_cluster_radius(cps_options *o, const char *radius)` / `const char *cps_options_cluster_radius(const cps_options *o)` | A positive radius sets the goal to `CPS_CLUSTER_RADIUS` with that radius. Zero clears the radius and sets the goal to `CPS_GOAL_DEFAULT`. NULL, a malformed literal or a negative value gives `CPS_INVALID_ARGUMENT`. This setter is coupled to the goal (§4.3). |
| `cps_error cps_options_set_absolute_tolerance(cps_options *o, const char *tolerance)` / `const char *cps_options_absolute_tolerance(const cps_options *o)` | A value ≥ 0; "0" disables it. A positive t also accepts a region when every root in it lies within t of the midpoint, as an alternative to the relative test. There is no implicit tolerance near zero. |
| `cps_error cps_options_set_real_only(cps_options *o, int enabled)` / `int cps_options_real_only(const cps_options *o)` | Nonzero means on. Asks for the real roots only (§5.2, §6.3). |
| `cps_error cps_options_set_region(cps_options *o, const char *real_lower, const char *real_upper, const char *imag_lower, const char *imag_upper)` / `cps_error cps_options_get_region(const cps_options *o, const char **real_lower, const char **real_upper, const char **imag_lower, const char **imag_upper)` | A closed rectangle. NULL or "" leaves a side unbounded, and four unbounded sides clear the region. The setter checks each end's literal only. The checks that relate two ends, or an end to the limits, are made at the solve (§4.3,). The getter writes NULL for an unbounded side. Meaning: §5.3. |
| `cps_error cps_options_set_work_budget(cps_options *o, uint64_t units)` / `uint64_t cps_options_work_budget(const cps_options *o)` | 0 selects the automatic budget of §4.1. Any other value is used as given. Meaning: §4.4. |
| `cps_error cps_options_set_live_cells(cps_options *o, uint64_t cells)` / `uint64_t cps_options_live_cells(const cps_options *o)` | 0 selects `CPS_DEFAULT_LIVE_CELLS`. |
| `cps_error cps_options_set_max_iterations(cps_options *o, uint64_t sweeps)` / `uint64_t cps_options_max_iterations(const cps_options *o)` | 0 means no limit. Meaning: §4.4. |
| `cps_error cps_options_set_time_limit(cps_options *o, double seconds)` / `double cps_options_time_limit(const cps_options *o)` | 0 or +∞ means no limit. NaN or a negative value gives `CPS_INVALID_ARGUMENT`. Meaning: §4.4. |
| `cps_error cps_options_set_cancel_token(cps_options *o, cps_cancel_token *token)` | Names the token that solves with these options poll. NULL clears it. The options keep the token alive (§1.2). |
| `cps_error cps_options_set_resource_limits(cps_options *o, const cps_resource_limits *limits)` / `cps_error cps_options_get_resource_limits(const cps_options *o, cps_resource_limits *out)` | NULL gives `CPS_INVALID_ARGUMENT`. Unknown flag bits give `CPS_INVALID_ARGUMENT`. The struct is copied. |
| `cps_error cps_options_set_workers(cps_options *o, uint32_t workers)` / `uint32_t cps_options_workers(const cps_options *o)` | 0 gives `CPS_INVALID_ARGUMENT`. A solve uses at most the machine's hardware threads; where it cannot start them it runs on its own thread, with the same answer. NULL options read 1. |

String getters return the literal as the caller last set it, trimmed of surrounding space. When nothing
was set they return "0" (the tolerance and the radius) or NULL (an unbounded region side). The pointer is
valid until the options are next changed or destroyed.

`cps_options_set_exact_budget` is removed. To get the same effect, set
`exact_materialization_bits`, and the expansion flag if wanted, in the resource limits.

**Resource limits:**

| C | C++ |
|---|---|
| `typedef struct cps_resource_limits { uint64_t exact_materialization_bits, exact_operations, temporary_bytes, flags; } cps_resource_limits;` | `struct ResourceLimits : cps_resource_limits { ResourceLimits() noexcept; };` (the constructor fills in the defaults) |
| `void cps_resource_limits_default(cps_resource_limits *out)` (NULL does nothing) | the constructor |
| `#define CPS_RESOURCE_ALLOW_EXPONENT_EXPANSION UINT64_C(1)` | the same macro |

### 4.3 Validation rules

- A setter checks only its own arguments: NULL, the literal grammar, signs, and ranges.
- There is one coupling. `set_cluster_radius` sets the goal, and `set_goal(CPS_CLUSTER_RADIUS)` needs a positive radius already set.
- Every other check that involves two settings, or a setting and the input kind, is made by the solve, before any work:
  - a combination that no input of that entry point's kind can satisfy is `CPS_INVALID_ARGUMENT` from the solve call (§5.2);
  - a refusal that depends on the input is an outcome (§5.2).
- At the solve, the region's ends must satisfy:
  - lower ≤ upper on each axis, decided exactly at no more than the precision cap;
  - each real end converts to an exact rational within the exact storage policy (§4.4).

  Otherwise the solve returns `CPS_INVALID_ARGUMENT`. All of these depend on the options alone.
- Settings that their query does not use are ignored, not rejected. Examples: a tolerance under the cluster-radius goal, or the requested precision under the isolation goal.

### 4.4 Meaning of the budgets and limits

**Work budget.** A soft, deterministic limit on controller work units.
- A unit is one of: a degree-proof or preparation pass; a backend call; a bounded local certification or refinement attempt; a disk-exclusion attempt; a square split; a frontier transition.
- It is checked between units, and the unit in progress always finishes. At high precision one unit can take a long time.
- The budget belongs to the call. Every engine of the call shares it, including a hand-off and the solves the library runs internally for a real-only query.
- When it is spent, the solve returns its certified cover with `CPS_WORK_LIMIT` and `CPS_STOP_WORK_BUDGET`.

**Live cells.** The most search squares that subdivision keeps at once, parents and untested children included. This is not a byte limit. A search blocked by it ends with `CPS_WORK_LIMIT` and `CPS_STOP_CELL_BUDGET`.

**Iteration budget.** At most this many sweeps of the classical and secular simultaneous iterations, at every precision, the binary64 sweeps included.
- They are counted in `cps_statistics.sweeps`.
- Not under this budget: subdivision; exact real-root isolation; the piecewise method's own work (its local-model solves appear in `sweeps` but do not count against the budget); and the axis charts. Their limits are the work and cell budgets, the time limit and the token. A solve that piecewise hands to secular is under the budget from the hand-off on.
- When it is spent, the solve returns its last certified cover with `CPS_WORK_LIMIT` and `CPS_STOP_ITERATION_BUDGET`. It is not a time bound.

**Time limit.**
- The solver measures it on a monotonic clock, from the moment the solve function is entered.
- It covers everything the call does, including every engine and every solve the library runs internally.
- It is checked at the same cooperative checkpoints as the token. The solve stops at the first checkpoint after the limit, with `CPS_CANCELLED`, `CPS_STOP_TIME_LIMIT` and its certified partial answer.
- The overrun is bounded only by the longest step that cannot be interrupted: one FLINT call, or one work unit.
- A time-limited answer depends on the machine and its load. The work and iteration budgets are deterministic.

**Cancellation token** (advanced).

| C | C++ |
|---|---|
| `typedef struct cps_cancel_token cps_cancel_token;` | a `std::stop_token`, named in `Options::stop` |
| `cps_error cps_cancel_token_create(cps_cancel_token **out)` | a `std::stop_source` of the caller's (`get_token()`) |
| `void cps_cancel_token_request(cps_cancel_token *t)` | `std::stop_source::request_stop()` |
| `int cps_cancel_token_requested(const cps_cancel_token *t)` | `std::stop_token::stop_requested()` |
| `void cps_cancel_token_destroy(cps_cancel_token *t)` | the `std::stop_source`'s destructor |

- A request is sticky: once requested, a token stays requested. To cancel again, make a new token.
- Every solve whose options name the token stops at its next checkpoint, with `CPS_CANCELLED`, `CPS_STOP_CANCELLED` and its certified partial answer. A solve whose token was requested before it started does no work.
- `cps_cancel_token_request` may be called from any thread at any time, including while solves poll the token. It is lock-free; the API does not promise async-signal safety.
- One token may serve any number of options and solves.
- The library runs no user code during a solve. Cancellation is the solver reading a flag.

**Resource limits.** Each is a per-stage estimate. None is a guarantee of peak RSS or wall time.
- `UINT64_MAX` means unspecified: no fixed ceiling. 0 is an explicit zero allowance.
- `exact_materialization_bits` bounds the admitted expanded coefficients, exact operand growth and the storage of applicable exact preparation.
- `exact_operations` bounds exact preparation effort, in word operations.
- `temporary_bytes` bounds estimated engine, preparation and construction storage. It also bounds exact bit storage, at 8 bits a byte.
- The *exact storage policy* is:
  - the exact storage limit is min(`exact_materialization_bits`, 8 · `temporary_bytes`);
  - exact preparation may write compact coefficients out when the written-out numbers exceed their owned representation by at most one precision cap (`max_work_bits`) each; a single number that must be exact on its own (a region's end, a secular term) may be expanded when the expanded value is no larger than its owned representation;
  - with `CPS_RESOURCE_ALLOW_EXPONENT_EXPANSION`, any expansion within the limit is allowed.
- `exact_materialization_bits` = 0 disables exact preparation, so the solve works on balls.
- A refused exact path falls back to ball solving, and the result reports `cps_result_preparation_limited`.

---

## 5. Solving

### 5.1 Entry points

| C signature | Contract |
|---|---|
| `cps_error cps_solve(const cps_polynomial *p, const cps_options *o, cps_result **out)` | One solve of a coefficient polynomial. |
| `cps_error cps_solve_secular(const cps_secular *s, const cps_options *o, cps_result **out)` | One solve of a native secular equation. |
| `cps_error cps_solve_program(const cps_program *p, const cps_options *o, cps_result **out)` | One solve of a program. The program must have an output. |

Common to all three C entry points:
- NULL options mean the defaults.
- Errors:
  - `CPS_INVALID_ARGUMENT`: NULL input or `out`; a program without an output; any usage error of §4.3 or §5.2.
  - `CPS_OUT_OF_MEMORY`;
  - `CPS_INTERNAL_ERROR`.

  On error, `*out` is NULL and nothing is published.
- On `CPS_OK`, `*out` is a new result that the caller owns. It is independent of the input and the options.
- The input and options must remain alive and unchanged while the call runs (§8).
- A C++ `solve` without explicit options uses a default-constructed `Options`.

### 5.2 What each input kind and method supports

Three outcomes appear in the tables:
- **answered**: the solve works on the query;
- **usage error**: `CPS_INVALID_ARGUMENT` from the solve and no result;
- **refusal**: `CPS_OK`, with termination `CPS_UNSUPPORTED_QUERY` and stop reason `CPS_STOP_INPUT`.

"Ball input" means a coefficient polynomial with at least one nonzero radius. "The default goal" is
`CPS_GOAL_DEFAULT`, resolved per input kind as §6.4 states.

**`cps_solve` (coefficients):**

| Query | `CPS_METHOD_AUTO` | `CPS_METHOD_CLASSICAL`, `_SECULAR`, `_PIECEWISE` |
|---|---|---|
| All complex roots; accuracy, approximation or default goal | answered (exact or ball input) | answered on exact input; refusal on ball input |
| All complex roots; cluster-radius goal | answered (exact or ball input) | usage error |
| All complex roots; isolation goal | answered on exact input; refusal on ball input | usage error |
| Real-only; accuracy, approximation, isolation or default goal; with or without a region | answered on exact real input; refusal on ball input or nonreal exact input | usage error |
| Real-only; cluster-radius goal | usage error | usage error |
| A region without real-only | usage error | usage error |

**`cps_solve_secular`:**

| Query | `CPS_METHOD_AUTO` | `CPS_METHOD_SECULAR` | `CPS_METHOD_CLASSICAL`, `_PIECEWISE` |
|---|---|---|---|
| All complex roots; accuracy, approximation or default goal | answered | answered | usage error |
| All complex roots; cluster-radius goal | answered | usage error | usage error |
| Real-only; accuracy, approximation or default goal | answered when the terms give a real numerator within the exact limits; refusal otherwise | usage error | usage error |
| Real-only with the cluster-radius goal; the isolation goal; a region | usage error | usage error | usage error |

**`cps_solve_program`:**

| Query | `CPS_METHOD_AUTO`, `_CLASSICAL`, `_SECULAR` | `CPS_METHOD_PIECEWISE` |
|---|---|---|
| All complex roots; approximation (the default for programs), accuracy (best effort, §6.4) or cluster-radius goal | answered | usage error |
| Real-only (any goal); the isolation goal; a region | usage error | usage error |

**Refusals and input terminations that depend on the input**, at every entry point:

| Situation | Termination | Stop reason |
|---|---|---|
| The input is zero: every coefficient exactly zero, or a program that computes zero | `CPS_ZERO_POLYNOMIAL` | `CPS_STOP_INPUT` |
| The degree is not proved: the top coefficient's box still contains zero at the precision cap, or a program's leading terms cancel beyond what the library resolves | `CPS_INDETERMINATE_DEGREE` | `CPS_STOP_INPUT` |
| The query is not available for this input (the "refusal" cells above) | `CPS_UNSUPPORTED_QUERY` | `CPS_STOP_INPUT` |
| Over the resource limits: a program's length or solve storage; a secular numerator's exact conversion; a real-only query whose exact real preparation exceeds the exact limits; an exact test at a region end that needs the input in exact integer form beyond the exact storage policy | `CPS_UNSUPPORTED_QUERY` | `CPS_STOP_INPUT` |

A nonzero constant input is answered as degree 0, with no regions, complete, the goal met and
`CPS_SUCCESS`. This holds for a constant polynomial, the empty secular sum and a constant program.

**Methods.**
- `CPS_METHOD_AUTO` chooses a route by a fixed, inexpensive rule. The choice is not part of the contract.
- A forced method names the principal numerical method. Shared exact preparation, direct formulas and the method's internal stages remain part of it.
- Forced piecewise hands to the secular method, within the same call and budgets, every stop that is not the caller's: its regime exit (Auto's where Auto would choose piecewise; elsewhere once its work passes a global solve's minimum budget), a storage refusal, a refusal because the coefficients' exponents exceed its range, and its precision budget. `cps_result_handoff` (§6.5) reports the hand-off and its reason. The caller's own stops stay with piecewise: a spent work budget, the time limit, the token.

### 5.3 Real-only queries and regions

- A real-only query asks for the real roots of an exact real input. Each comes alone in a region proved real (§6.3).
- Under Auto the answer comes from the first of three routes that answers:
  1. the sparse recursion, on inputs with few terms;
  2. the piecewise method's axis charts;
  3. the complex route.

  Each route stops honestly and hands the query on. The work of every route counts against the call's budgets. Which route answered is not part of the contract.
- A region is the closed rectangle [real_lower, real_upper] × [imag_lower, imag_upper].
  - The query answers exactly the real roots x with real_lower ≤ x ≤ real_upper, when imag_lower ≤ 0 ≤ imag_upper.
  - A rectangle that misses the real axis gets an empty, complete answer, without work.
  - `complete` and `goal_met` refer to the roots in the region.
- The ends are closed.
  - A root exactly at an end is decided by an exact test and is in the answer. When the end is dyadic, the region is that exact point.
  - A root of odd multiplicity whose interval holds an end, but is not equal to it, is placed on its side by the input's signs, and its interval is cut at a dyadic point between the end and the root.
  - A root of even multiplicity is refined until its interval lies on one side. On the complex route that refines every root.
- A region saves no time on the complex route.

---

## 6. Results

### 6.1 Functions

| C signature | Contract |
|---|---|
| `void cps_result_destroy(cps_result *r)` | NULL does nothing. |
| `int64_t cps_result_degree(const cps_result *r)` | The input's degree: −1 for the zero input. With `CPS_INDETERMINATE_DEGREE`, an upper bound. For secular input, the number of terms. For a program, the proved degree, or its degree bound when it was refused. Not reduced by real-only filtering. |
| `cps_termination cps_result_termination(const cps_result *r)` | How the solve ended (§6.2). |
| `cps_stop_reason cps_result_stop_reason(const cps_result *r)` | Why it stopped (§6.2). |
| `int cps_result_complete(const cps_result *r)` | Whether the regions cover every root asked for (§6.3). |
| `int cps_result_goal_met(const cps_result *r)` | Whether the goal was met (§6.4). |
| `cps_goal cps_result_goal(const cps_result *r)` | The goal the solve applied. Never `CPS_GOAL_DEFAULT`. |
| `uint64_t cps_result_work_bits(const cps_result *r)` | The last outer working precision. This is not the backend's peak precision. |
| `size_t cps_result_size(const cps_result *r)` | The number of regions. |
| `const char *cps_result_diagnostics(const cps_result *r)` | Text for people (§6.7). It is borrowed until the result is destroyed. |
| `cps_error cps_result_get_root(const cps_result *r, size_t index, cps_root_view *out)` | Region `index`. NULL gives `CPS_INVALID_ARGUMENT`; an index at or past `size` gives `CPS_OUT_OF_RANGE`. The view's strings are borrowed until the result is destroyed. The C++ copies own their strings and outlive the result. |
| `cps_error cps_result_get_approximation(const cps_result *r, size_t index, double *real, double *imag)` | §6.6. |
| `uint64_t cps_result_sufficient_work_bits(const cps_result *r)`; `cps_squarefree cps_result_squarefree(const cps_result *r)`; `cps_handoff cps_result_handoff(const cps_result *r)`; `int cps_result_preparation_limited(const cps_result *r)` | §6.5. |
| `cps_error cps_result_get_statistics(const cps_result *r, cps_statistics *out)`; `cps_error cps_result_get_engine_counters(const cps_result *r, cps_engine_counters *out)`; `size_t cps_result_component_count(const cps_result *r)`; `cps_error cps_result_get_component(const cps_result *r, size_t index, cps_component_progress *out)` | §7. |

### 6.2 Invariants

These hold for every result. Each is checked by a test.

- **I1.** `goal_met` implies `complete`.
- **I2.** `termination == CPS_SUCCESS` if and only if `goal_met`. A cancellation or time limit observed after the goal is met does not change the termination.
- **I3.** `stop_reason == CPS_STOP_NONE` if and only if `termination == CPS_SUCCESS`.
- **I4.** Termination, stop reason, coverage and regions go together as follows:

| `termination` | `complete` | `goal_met` | `stop_reason` | Regions |
|---|---|---|---|---|
| `CPS_SUCCESS` | 1 | 1 | `CPS_STOP_NONE` | the answer |
| `CPS_PRECISION_LIMIT` | 0 or 1 | 0 | `CPS_STOP_PRECISION`; or `CPS_STOP_INPUT` when the coefficients' own uncertainty, not the precision, blocks every pending region | certified; possibly partial |
| `CPS_WORK_LIMIT` | 0 or 1 | 0 | `CPS_STOP_WORK_BUDGET`, `CPS_STOP_CELL_BUDGET` or `CPS_STOP_ITERATION_BUDGET` | certified; possibly partial |
| `CPS_CANCELLED` | 0 or 1 | 0 | `CPS_STOP_CANCELLED` (token) or `CPS_STOP_TIME_LIMIT` | certified; possibly partial |
| `CPS_UNSUPPORTED_QUERY` | 0 | 0 | `CPS_STOP_INPUT` | none, or the regions certified before the refusal |
| `CPS_ZERO_POLYNOMIAL` | 0 | 0 | `CPS_STOP_INPUT` | none (degree −1) |
| `CPS_INDETERMINATE_DEGREE` | 0 | 0 | `CPS_STOP_INPUT` | none |

- **I5.** Every exported region carries a proved count, whatever the termination.
- **I6.** The regions are pairwise disjoint as closed rectangles.
- **I7.** A complete all-complex answer has counts that sum to the degree. A complete real-only answer holds every real root in the requested domain, each in exactly one region.
- **I8.** `multiplicity` is 0 or equal to `count`.
- **I9.** `realness == CPS_REAL` implies `imag_lower = imag_upper = 0`. A test checks this on every route.

### 6.3 Regions

```c
typedef struct cps_dyadic_view { const char *mantissa; const char *exponent; } cps_dyadic_view;
typedef struct cps_root_view {
    cps_dyadic_view real_lower, real_upper, imag_lower, imag_upper;
    uint64_t count;
    uint64_t multiplicity;
    cps_realness realness;
    const char *certificate;
} cps_root_view;
typedef enum cps_realness { CPS_REAL_UNKNOWN = 0, CPS_REAL = 1, CPS_NONREAL = 2 } cps_realness;
```

In C++ a region is a `Root`: its `box` (`Box { Interval real, imag; }`, each `Interval { Dyadic lo, hi; }` with exact
ends), its optional `value` (the approximation as a `std::complex<double>`), `count`, `multiplicity` (empty when not known), `realness`
and `certificate`.

- A region is the closed rectangle [real_lower, real_upper] × [imag_lower, imag_upper]. Each end is exactly mantissa · 2^exponent, with both written as decimal integer strings.
- `count` is the number of roots in the region, counted with multiplicity, and it is proved.
- `multiplicity` is 0 when not established. When it is established, it equals `count`, and the region holds one distinct root of that multiplicity. A count is never a claimed multiplicity.
- `realness`:
  - `CPS_REAL`: every root in the region is real, and the region is a real interval;
  - `CPS_NONREAL`: the region misses the real axis, so none of its roots is real;
  - `CPS_REAL_UNKNOWN`: neither was proved.

  A program's region is `CPS_REAL` only by conjugation, when every constant of the program is real: a region of the
  complete cover holding one distinct root whose reflection in the real axis meets no other region (the coefficient
  route's proof). A program with a non-real constant has no `CPS_REAL` region.
- Every region of a real-only answer is `CPS_REAL` and holds one distinct real root with its multiplicity.
- `certificate`: see §6.7.
- Root order: regions come in increasing order of their exact midpoint's real part, then its imaginary part, on every route.

### 6.4 Goals

```c
typedef enum cps_goal {
    CPS_ROOT_ACCURACY = 0, CPS_CLUSTER_RADIUS = 1, CPS_ROOT_APPROXIMATION = 2, CPS_ROOT_ISOLATION = 3,
    CPS_GOAL_DEFAULT = 4                                   /* */
} cps_goal;
```

Each goal below requires `complete`, plus the condition stated. "The accuracy test" for a region with
midpoint c is: every root z in the region satisfies abs(z − c) ≤ 2^−bits · abs(z), or, with a positive
absolute tolerance t, abs(z − c) ≤ t. It is checked with outward bounds.

- **`CPS_ROOT_ACCURACY`**: every region holds one distinct root with its multiplicity established, and passes the accuracy test.
- **`CPS_ROOT_APPROXIMATION`**: every region passes the accuracy test. Counted clusters are allowed (count > 1, multiplicity 0). Real-only queries answer it as `CPS_ROOT_ACCURACY`.
- **`CPS_CLUSTER_RADIUS`**: every region lies within the absolute radius r of its midpoint, checked by the sum of the half-widths. Neither separation nor multiplicity is required. Ccluster's 3-disk separation property is not promised.
- **`CPS_ROOT_ISOLATION`**: every distinct root lies alone in a disjoint closed rectangle with its multiplicity established. A real-only query returns closed real intervals. This goal requires Auto and exact coefficient input (exact real input for a real-only query). No width is promised, and the requested precision is not used.
- **`CPS_GOAL_DEFAULT`**: the initial value of the setting. The solve resolves it per input kind: `CPS_ROOT_ACCURACY` for coefficients and secular input, `CPS_ROOT_APPROXIMATION` for a program. `cps_result_goal` reports the resolved goal.
  - The reason programs differ is recorded in the program-input spec (`docs/slp-input-spec-20261015.md` §2): a program can never certify a multiplicity above one, so accuracy is only a best effort on a program.
- **Programs under `CPS_ROOT_ACCURACY`**: a cluster that does not separate within the precision cap ends the call with `CPS_PRECISION_LIMIT` and the counted cluster. `cps_result_squarefree` then says whether the roots are known to be distinct (§6.5).

### 6.5 Facts

```c
typedef enum cps_squarefree { CPS_SQUAREFREE_UNKNOWN = 0, CPS_SQUAREFREE_PROVED = 1, CPS_SQUAREFREE_DISPROVED = 2 } cps_squarefree;
typedef enum cps_handoff {
    CPS_HANDOFF_NONE = 0, CPS_HANDOFF_REGIME_EXIT = 1, CPS_HANDOFF_STORAGE = 2,
    CPS_HANDOFF_RANGE = 3, CPS_HANDOFF_PRECISION = 4
} cps_handoff;
```

- **`cps_result_sufficient_work_bits`** is 0 unless both of these hold:
  - the solve stopped at its precision cap (`CPS_PRECISION_LIMIT` with `CPS_STOP_PRECISION`);
  - the input is exact and real, so that its distinct factors are integer polynomials (Gaussian-rational and ball inputs give 0).

  When both hold, it is 2 · (separation bits + size bits + requested bits), from Mahler's root-separation bound. That is a cap at which the roots can be isolated at the requested precision. It is a bound from mathematics, not a promise that a solve under that cap meets its goal within its other budgets.
- **`cps_result_squarefree`**:
  - `PROVED` means the input is proved square-free: every root is simple.
  - `DISPROVED` means a root of multiplicity above one is proved.
  - `UNKNOWN` means neither.
  - For a program the verdict comes from its image modulo a word-size prime, so it is `PROVED` or `UNKNOWN`.
  - For coefficient input it comes from exact preparation when that ran, and is `UNKNOWN` when preparation was limited or the input is a ball input.
- **`cps_result_handoff`** reports that the piecewise method, forced or chosen by Auto, handed the solve to the secular method within the call, and why:
  - `REGIME_EXIT`: it spent more than a global solve's predicted cost;
  - `STORAGE`: its storage was refused;
  - `RANGE`: the coefficients' exponents exceed its working range;
  - `PRECISION`: it hit its precision budget.

  The answer then comes from the secular method.
- **`cps_result_preparation_limited`**: exact preparation was refused or limited by the resource limits, so multiplicities may be unknown.
- Certificates stay informational text (§6.7). The diagnostics may repeat these facts, but these getters are the reference.

### 6.6 Approximations

`cps_result_get_approximation(r, i, &re, &im)` is an explicit, lossy conversion of region i's midpoint to
binary64, rounded to nearest. It fails, leaving the outputs unchanged, in these cases:

| Case | Code |
|---|---|
| NULL `r`, `real` or `imag` | `CPS_INVALID_ARGUMENT` |
| `i` ≥ `size` | `CPS_OUT_OF_RANGE` |
| the result is incomplete or its goal is unmet; the region is a cluster (multiplicity ≠ count); the goal is `CPS_ROOT_ISOLATION`, which promises no width | `CPS_UNRESOLVED` |
| a component overflows binary64, or a nonzero component would round to zero | `CPS_OUT_OF_RANGE` |

Nonzero subnormal results are allowed, and they may lose relative accuracy. The requested precision
applies to the exact region, not to this conversion.

In C++ each `Root` carries its `value`, a `std::optional<std::complex<double>>` holding exactly the number
`cps_result_get_approximation` gives for its region and empty where that call refuses, and `Result::values()` lists
the values, each repeated as many times as its region's count (empty after a stop). `Root::centre()` is the box's
centre as doubles, for any region. C and C++ never give different values.

### 6.7 Diagnostics and certificates

- `cps_result_diagnostics` is text for people: the preparation and adaptation trace.
  - Its format is not part of the contract and may change in any release.
  - It does not depend on the C or C++ locale.
  - Every fact a caller needs is also available from a getter (§6.5).
- `certificate` is a short label for people, naming how the region's count was proved, for example `pellet-local-count`. The vocabulary is not fixed. It is not a replayable proof object.

---

## 7. Statistics and progress

Statistics are observations, not certificates. A field may be zero for a stage that did not run.

**The stable block** is a promise. The fields have exactly the meanings given here.

```c
typedef struct cps_statistics {
    size_t   size;               /* in: sizeof(cps_statistics); out: the bytes filled */
    uint64_t work_units;         /* controller work units the call spent, every engine and internal solve included */
    uint64_t sweeps;             /* sweeps of the simultaneous iterations; what the iteration budget counts */
    uint64_t precision_increases;/* raises of the outer working precision */
    uint64_t final_work_bits;    /* the last outer working precision (= cps_result_work_bits) */
    uint64_t highest_work_bits;  /* the highest outer working precision reached */
    uint64_t regions;            /* exported regions (= cps_result_size) */
    uint64_t roots;              /* the sum of the regions' counts */
    uint64_t resolved_regions;   /* regions whose multiplicity is established */
    uint64_t exact_factors;      /* factors produced by exact preparation */
    uint64_t total_ns;           /* wall time of the call, on the solve's monotonic clock */
    uint64_t preparation_ns;     /* wall time of exact preparation, every engine included */
} cps_statistics;
CPS_API cps_error cps_result_get_statistics(const cps_result *, cps_statistics *out);
```

How the size field works:
- The caller sets `out->size = sizeof(cps_statistics)`.
- A size smaller than the 0.5 layout gives `CPS_INVALID_ARGUMENT`, and `*out` is left unchanged.
- Otherwise the library fills the fields that fit, leaves the bytes beyond untouched, and sets `size` to the number of bytes it filled.
- NULL gives `CPS_INVALID_ARGUMENT`.
- The C++ `statistics()` sets the size itself.
- The stage times promised are `total_ns` and `preparation_ns`; the finer stages are diagnostics.

**The engine counters** are diagnostics. Fields may be added, removed or renamed in any release. No
meaning is promised beyond the comment on each field.

```c
typedef struct cps_engine_counters {
    size_t   size;               /* as in cps_statistics */
    /* stage times, ns, overlapping */
    uint64_t input_enclosure_ns, isolation_ns, coverage_validation_ns, backend_ns, taylor_shift_ns,
             counting_ns, refinement_ns, search_ns, candidate_ns, evaluation_ns, interaction_ns,
             validation_ns, regeneration_ns;
    /* counts */
    uint64_t precision_attempts, isolation_passes, disk_tests, graeffe_steps, taylor_shifts,
             taylor_cache_hits, boxes_split, boxes_excluded, newton_attempts, newton_successes,
             backend_calls, evaluations, evaluation_terms, root_interactions, validation_calls,
             secular_regenerations, double_seed_calls, candidate_stagnations, real_refinement_steps,
             real_bisections;
    /* the search's state at the end */
    uint64_t pending_components, live_cells, min_pending_bits, max_pending_bits, input_uncertainty;
} cps_engine_counters;
CPS_API cps_error cps_result_get_engine_counters(const cps_result *, cps_engine_counters *out);
```

A struct with a size field, so that fields can be added behind the ones a caller knows.

**Component progress** is also diagnostics. Its indices identify entries in this result only: they are
not region indices and not persistent handles.

```c
typedef enum cps_component_state {
    CPS_COMPONENT_ACTIVE = 0, CPS_COMPONENT_FINISHED = 1, CPS_COMPONENT_PRECISION_LIMIT = 2, CPS_COMPONENT_CELL_LIMIT = 4
} cps_component_state;                                     /* 3 (stagnation) is never produced */
typedef struct cps_component_progress {
    uint64_t factor_index, work_bits, live_cells, certified_count;
    cps_component_state state;
} cps_component_progress;
CPS_API size_t   cps_result_component_count(const cps_result *);
CPS_API cps_error cps_result_get_component(const cps_result *, size_t index, cps_component_progress *out);
```

`cps_result_get_component` gives `CPS_INVALID_ARGUMENT` for NULL and `CPS_OUT_OF_RANGE` for an index at or
past the count.

**Progress.** There is no separate progress snapshot, because a solve is one call. Today's progress fields
are placed as follows:

| Today's field | New place |
|---|---|
| the stop reason | `cps_result_stop_reason` |
| the call's work | `cps_statistics.work_units` |
| the highest precision | `cps_statistics.highest_work_bits` |
| the preparation flag | `cps_result_preparation_limited` |
| everything else | the engine counters |

```c
typedef enum cps_stop_reason {
    CPS_STOP_NONE = 0, CPS_STOP_WORK_BUDGET = 1, CPS_STOP_CELL_BUDGET = 2, CPS_STOP_PRECISION = 3,
    CPS_STOP_CANCELLED = 4, CPS_STOP_INPUT = 5, CPS_STOP_ITERATION_BUDGET = 6,
    CPS_STOP_TIME_LIMIT = 7                                /* */
} cps_stop_reason;
```

---

## 8. Threads, lifetimes and FLINT

### 8.1 Guarantees

- **T1. Concurrent solves.** Any number of threads may call the library at the same time, each on objects it does not share. (True in today's code: audit A1-A6.)
- **T2. Read-only sharing.** An input, an options object or a result may be read by any number of threads at once, provided no thread changes or destroys it meanwhile. Reading includes solving with the input or the options.
- **T3. Moving between threads.** An object may be created on one thread and used or destroyed on another, including after the creating thread has called `cps_cleanup_thread`.
  - This relies on FLINT's thread-aware memory manager (audit B2).
  - The build must therefore require a FLINT built with thread-local storage (enforced today) and with thread support, `FLINT_USES_PTHREAD` (to enforce: audit B9).
- **T4. Tokens.** `cps_cancel_token_request` may be called on any thread at any time, including concurrently with solves that poll the token and with other requests.
- **T5. Changes are exclusive.** A setter, a builder or a destroy on an object must not run while any other thread uses that object. Doing so is a data race and undefined behaviour.
- **T6. No hidden activity.**
  - No user code runs inside a solve. `Program::from` runs f before any solve, and f may call any library function.
  - The library starts no threads, except that a solve whose options ask for more than one worker starts up to
    workers − 1 threads and joins them before it returns. No thread outlives the call. The
    workers never use FLINT's thread pool and start with every signal blocked.
  - It keeps no mutable global state apart from the per-thread error message.
  - It reads no environment variables and writes nothing to `stdout` or `stderr`.
- **T7. Re-entrancy.** Any function may be called on any thread at any time, subject to T5. The API does not promise async-signal safety for any library function.
- **T8. The same answer at every worker count.** A solve's answer, counters and diagnostics do not
  depend on the worker count. A stop by a time limit or a token depends on the machine already, and still does.

### 8.2 Per-thread state

- `cps_last_error` is per thread (§2.1).
- FLINT keeps thread-local caches: arb constants, prime tables, FFT contexts, and free `mpz` lists. A thread's caches grow with the precision it has used.
- **Floating-point environment.** A solve saves the calling thread's floating-point environment at entry. It runs with round-to-nearest, all exceptions masked, and flush-to-zero and denormals-are-zero off. It restores the saved environment before it returns.
- A call may change `errno`.

### 8.3 FLINT, GMP and the host

- **FLINT** 3.6 or later is required, built with thread-local storage (enforced at build time) and, under T3, thread support. `cps_backend_version()` returns FLINT's version string.
- **Fatal failures.** An error inside FLINT or GMP, including an allocation failure inside them, calls `abort()`. That ends the process and every solve on every thread. The library does not intercept it, and installing a throwing FLINT error handler would not help, because it would unwind through C frames. The resource limits bound FLINT's allocations only as estimates.
- **Memory functions.** The host must not install GMP or FLINT memory functions that are not thread-safe. It must not change them while any solve runs.
- **FLINT's thread pool.** The library never calls `flint_set_num_threads`. A host that calls it while a solve uses FLINT's pool makes FLINT abort. A host thread that raised its FLINT worker count lets FLINT routines inside its solves use the shared pool. The library documents this as a host rule (read in FLINT 3.6.0's source): holding the count at one would not keep a solve's bits from depending on the pool, because FLINT's algorithm choices read the pool's idle threads, not the thread's count (under Reproducibility below). So the host rule is the one that matters.
- **`void cps_cleanup_thread(void)`**, or a C++ `ThreadScope` whose destructor calls it.
  - What it does:
    - frees the calling thread's FLINT caches;
    - frees MPFR's caches, which MPFR may share between threads.
  - When to call it:
    - before a thread that used the library exits;
    - at a point where no other FLINT user on that thread holds pointers into those caches;
    - never while that thread is inside a library call (which is impossible without callbacks).
  - It is never needed for correctness. Skipping it leaks the caches of a short-lived thread.
  - Under T3, objects created on the thread stay valid after it.
  - `ThreadScope` is constructed first at a worker's entry, and scopes are not nested.
- **Stack.** A solve runs on the caller's stack, and FLINT and GMP take temporaries from the stack. A host thread should give a solve at least 512 KiB of stack. Measured over the test set, every method and query, a solve's peak is at most 153 KiB above the 16 KiB an empty thread takes (median 23 KiB, 99th percentile 102 KiB); the deepest are high-precision inputs, where GMP and FLINT put temporaries on the stack.
  - **Proposed statement:** a thread that calls the library needs at least 512 KiB of stack (the default for secondary threads on macOS; glibc's default is 8 MiB, Windows' 1 MiB), three times the largest measured peak.
- **Reproducibility.**
  - Above 4,592 bits, FLINT's thread-local π and log 2 are rounded from a cache whose precision depends on earlier work on the same thread. An exported endpoint can then differ in near-tie cases. Certificates are unaffected (audit B3).
  - Above about 10,000 bits (20,000 in some cases), FLINT 3.6's arb exp, sin and cos, and atan choose their algorithm by `flint_get_num_available_threads()` (`arb/exp_arf.c:110`, `arb/sin_cos_arf_generic.c:295`, `arb/atan_arf_bb.c:317`, `arb/exp_arf_bb.c:231` and the `*_sum_bs_powtab.c` files). Once FLINT's process-wide pool exists, that function returns the pool's idle threads plus one, whatever the calling thread's own setting (`thread_support/get_num_available_threads.c`). So in a process whose host has created the pool (`flint_set_num_threads` above one), the endpoints a high-precision solve exports can change with the machine's load: each is still a valid enclosure, and certificates hold. The library never creates the pool, so its own answers do not depend on load.
  - A test of bit-identical output should run each solve on a cleaned thread.
  - Answers under a time limit depend on the machine's speed.

### 8.4 Lifetimes, all in one place

| Object or pointer | Valid until |
|---|---|
| A handle | Its destroy call. Under T3 it may be used on any thread. |
| A token | The last of: the caller's destroy, and every options object and solve that holds it |
| Strings and views from a result | The result's destroy |
| A string from an options getter | The next change to those options, or their destroy |
| `cps_last_error()` | The calling thread's next `cps_error`-returning call; only on that thread |
| `cps_version()`, `cps_backend_version()` | Forever (static strings) |
| A `ProgramValue` | The end of the call of f that made it |
| A C++ snapshot (`Root`, `std::string`, statistics structs) | Its own lifetime; it owns its data |

---

## 9. Versioning

| C | C++ | Meaning |
|---|---|---|
| `#define CPS_VERSION_MAJOR 0`, `CPS_VERSION_MINOR 5`, `CPS_VERSION_PATCH 1`, `CPS_VERSION_STRING "0.5.1"` | the same macros | The version of the headers. |
| `const char *cps_version(void)` | `const char *version()` | The library's version, the same string as `CPS_VERSION_STRING` in the headers it was built with. It is static. |
| `#define CPS_ABI_VERSION 2u`; `uint32_t cps_abi_version(void)` | `uint32_t abi_version()` | Informational. Raised whenever the binary interface changes. |
| `const char *cps_backend_version(void)` | `const char *backend_version()` | FLINT's version string. |

- Before 1.0.0 there is no compatibility promise of any kind. Any release may change the API, the ABI and the behaviour.
- A program may compare `cps_abi_version()` with `CPS_ABI_VERSION`, and `cps_version()` with `CPS_VERSION_STRING`, to detect a mismatched library. The library itself does not check.

---

---

## Appendix. The C declarations at a glance

```c
/* versions */           CPS_VERSION_MAJOR/MINOR/PATCH/STRING, CPS_ABI_VERSION; cps_abi_version, cps_version, cps_backend_version
/* errors, threads */    cps_error; cps_last_error; cps_cleanup_thread
/* resource limits */    cps_resource_limits, CPS_RESOURCE_ALLOW_EXPONENT_EXPANSION; cps_resource_limits_default
/* polynomial */         cps_polynomial_create, _destroy, _set_coefficient, _set_coefficient_d, _set_ball_coefficient
/* secular */            cps_secular_term_view; cps_secular_create, cps_secular_destroy
/* program */            CPS_NO_NODE; cps_program_create, _destroy, _variable, _constant, _add, _subtract, _multiply,
                         _negate, _power, _set_output, _parse, _size, _output
/* token */              cps_cancel_token_create, _request, _requested, _destroy
/* options */            cps_options_create, _copy, _destroy;
                         set/get: _method, _precision, _max_work_bits, _goal, _cluster_radius, _absolute_tolerance,
                         _real_only, _region (get_region), _work_budget, _live_cells, _max_iterations, _time_limit,
                         _workers, _resource_limits (get_resource_limits); set only: _cancel_token
                         CPS_DEFAULT_BITS, CPS_DEFAULT_MAX_WORK_BITS, CPS_MIN_MAX_WORK_BITS,
                         CPS_DEFAULT_WORK_UNITS, CPS_DEFAULT_WORK_UNITS_PER_ROOT, CPS_DEFAULT_LIVE_CELLS
/* solving */            cps_method, cps_goal; cps_solve, cps_solve_secular, cps_solve_program
/* results */            cps_termination, cps_stop_reason, cps_realness, cps_dyadic_view, cps_root_view;
                         cps_result_destroy, _degree, _termination, _stop_reason, _complete, _goal_met, _goal,
                         _work_bits, _size, _diagnostics, _get_root, _get_approximation
/* facts */              cps_squarefree, cps_handoff; cps_result_sufficient_work_bits, _squarefree, _handoff,
                         _preparation_limited
/* statistics */         cps_statistics, cps_engine_counters, cps_component_state, cps_component_progress;
                         cps_result_get_statistics, _get_engine_counters, _component_count, _get_component
/* removed */            cps_cancel_callback, cps_options_set_cancel, cps_options_set_exact_budget,
                         cps_options_set_search_budget, cps_secular_create_with_resources_v1, every _v1/_v2/_v3 name,
                         cps_result_get_progress_v1, CPS_COMPONENT_STAGNATION
```
