# CertPolySolve for C++: a quick start

Five short programs, one for each common task. Each is a complete file in
[`examples/quickstart/`](../examples/quickstart/) that the test suite builds and runs, and the output under it is
what it printed. [More examples](examples-cpp.md) cover the rest by topic, and [the guide](api.md) has every function.

1. [The first solve](#1-the-first-solve)
2. [Real roots](#2-real-roots)
3. [Many digits, exact coefficients](#3-many-digits-exact-coefficients)
4. [Repeated roots and close roots](#4-repeated-roots-and-close-roots)
5. [Time limits](#5-time-limits)
6. [When a solve stops early](#6-when-a-solve-stops-early)

Building: the header needs C++20 (GCC 15 and Clang 21 are tested) and these programs C++23 for `std::println`.
[The README](../README.md#install) says how to build and install the library; then, in CMake,
`find_package(CertPolySolve REQUIRED)` and `target_link_libraries(app PRIVATE CertPolySolve::certpolysolve)`. Everything is in one
header, `<certpolysolve/certpolysolve.hpp>`.

## 1. The first solve

The three roots of x³ − 2x − 5:

```cpp
#include <certpolysolve/certpolysolve.hpp>
#include <print>

int main() {
    certpolysolve::Polynomial p{-5, -2, 0, 1}; // constant term first: -5 - 2x + 0x^2 + x^3
    certpolysolve::Result result = certpolysolve::solve(p);

    for (const certpolysolve::Root &root : result.roots()) std::println("{}", root);
}
```

```text
-1.0472757407711633 - 1.1359398890889283i
-1.0472757407711633 + 1.1359398890889283i
2.0945514815423265
```

Three things to know:

- **The constant term comes first.** `Polynomial{-5, -2, 0, 1}` is −5 − 2x + x³. If your list has the highest
  power first, as MATLAB's and NumPy's `roots` do, write `Polynomial::highest_first({1, 0, -2, -5})`. Printing the
  polynomial (`std::println("{}", p)`) shows which one you built.
- **A `double` is used at its exact value.** `0.1` in C++ is the double nearest one tenth, and that is the number
  CertPolySolve uses. When you mean one tenth, write it as text: `"1/10"` or `"0.1"` (section 3).
- **Every root is correct to 15 digits** by default (53 bits of relative accuracy, a double's worth), and you can
  ask for more. `{}` prints a root's centre as a double. `result.values()` is all the roots as
  `std::complex<double>` once the goal is met, a repeated root repeated (a `std::complex<double>` has no `std::format` support before
  C++26: print `z.real()` and `z.imag()`, or print the `Root`). `root.box` is the box proved to hold the root, with
  exact ends.

## 2. Real roots

All the real roots, the real roots in an interval, the positive ones, and only how many there are:

```cpp
#include <certpolysolve/certpolysolve.hpp>
#include <print>
#include <string_view>

using certpolysolve::Rounding;

void print(std::string_view title, const certpolysolve::Result &result) {
    std::println("{}: {}", title, result.roots().size());
    for (const certpolysolve::Root &root : result.roots()) // the interval's ends as doubles, rounded outward
        std::println("  [{}, {}]", root.box.real.lo.to_double(Rounding::down), root.box.real.hi.to_double(Rounding::up));
}

int main() {
    certpolysolve::Polynomial quintic{-1, -1, 0, 0, 0, 1}; // x^5 - x - 1: one real root of five
    print("real roots of x^5 - x - 1", certpolysolve::real_roots(quintic));

    certpolysolve::Polynomial t5{-0.5, 5, 0, -20, 0, 16}; // 16x^5 - 20x^3 + 5x - 1/2 (the double -0.5 is exactly -1/2)
    print("positive roots of T5(x) - 1/2", certpolysolve::real_roots(t5, certpolysolve::at_least(0)));
    print("roots in [-1/2, 1/2]", certpolysolve::real_roots(t5, certpolysolve::between("-1/2", "1/2"))); // closed: ends count

    // each root alone in an interval, not refined: the cheapest query
    const certpolysolve::Result count = certpolysolve::real_roots(t5, {.goal = certpolysolve::Goal::isolation()});
    std::println("T5(x) - 1/2 has {} distinct real roots", count.roots().size());
}
```

```text
real roots of x^5 - x - 1: 1
  [1.1673039782614185, 1.1673039782614187]
positive roots of T5(x) - 1/2: 3
  [0.10452846326765347, 0.10452846326765348]
  [0.49999999999999994, 0.5000000000000001]
  [0.9781476007338056, 0.9781476007338057]
roots in [-1/2, 1/2]: 2
  [0.10452846326765347, 0.10452846326765348]
  [0.5, 0.5]
T5(x) - 1/2 has 5 distinct real roots
```

- Every root in the answer is proved real, and no real root is missed.
- `real_roots(p)` gives the real roots, `real_roots(p, region)` those in a closed region: `between(a, b)`,
  `at_least(a)` or `at_most(b)`, with ends given as any number, `0` or `"-1/2"`. Options follow the region:
  `real_roots(p, between(0, 1), {.bits = 200})`. `solve(p)` gives every complex root.
- `root.box.real.lo` and `.hi` are the interval's exact ends, numbers m · 2^e (`Dyadic`); `to_double(Rounding::down)`
  and `to_double(Rounding::up)` round them outward to doubles.
- Options are a struct with named fields, `{.bits = 64}` or `{.goal = certpolysolve::Goal::isolation()}`, in their
  declared order; the fields you leave out keep their defaults.
- Real-only queries need real, exact coefficients: doubles and text are exact, coefficients given within bounds
  are not.

## 3. Many digits, exact coefficients

Any coefficient can be given as text, exactly, and any number of digits can be asked for:

```cpp
#include <certpolysolve/certpolysolve.hpp>
#include <print>
#include <string>
#include <string_view>
#include <vector>

void print_root(std::string_view title, const certpolysolve::Polynomial &p, uint64_t bits, unsigned digits) {
    // the real roots, each to `bits` bits of relative accuracy, printed with `digits` significant digits
    for (const certpolysolve::Root &root : certpolysolve::real_roots(p, {.bits = bits}).roots())
        std::println("{}{}", title, root.box.real.mid().format(digits));
}

int main() {
    // A double is used at its exact binary value: the double 0.1 is 0.1000000000000000055511151231257827...
    certpolysolve::Polynomial from_double{-0.1, 1}; // x - (the double nearest 0.1)
    certpolysolve::Polynomial from_text{"-1/10", 1}; // x - 1/10, exactly
    print_root("root of x - 0.1 (double): ", from_double, 120, 30);
    print_root("root of x - 1/10 (text):  ", from_text, 120, 30);

    // Text: integers of any size, p/q, decimals, scientific notation, m@e for m * 2^e.
    certpolysolve::Polynomial mixed{"-123456789012345678901234567890", "2/3", "1.5e-3", "3@-2"};
    print_root("a real root of a mixed cubic: ", mixed, 53, 17);

    // A whole list of text coefficients at once, ascending, as when they are read from a file.
    const std::vector<std::string> lines{"-2", "0", "0", "1"};
    certpolysolve::Polynomial cube(lines); // x^3 - 2
    print_root("2^(1/3) to 100 digits: ", cube, 340, 100); // 340 bits > 100 * log2(10)
}
```

```text
root of x - 0.1 (double): 0.100000000000000005551115123126
root of x - 1/10 (text):  0.1
a real root of a mixed cubic: 5480471259.847853
2^(1/3) to 100 digits: 1.259921049894873164767210607278228350570251464701507980081975112155299676513959483729396562436255094
```

- Text forms: integers of any length, `p/q`, decimals such as `0.1` or `-2.5e-7` (exact: `"0.1"` is one tenth), and
  `m@e` for m · 2^e. A complex coefficient is `certpolysolve::Number::complex("1/3", -2)`.
- `.bits` asks for that many bits of relative accuracy in every root; a decimal digit is about 3.3 bits.
- `format(n)` and `{:.n}` write at most n significant digits, correctly rounded; an exact root prints exactly
  (`0.5`), and trailing zeros are left out.

## 4. Repeated roots and close roots

```cpp
#include <certpolysolve/certpolysolve.hpp>
#include <print>
#include <string_view>

void print(std::string_view title, const certpolysolve::Result &result) {
    std::println("{} ({})", title, result.ok() ? "goal met" : "stopped");
    for (const certpolysolve::Root &root : result.roots()) {
        std::print("  real {}, imaginary {}: count {}, ", root.box.real.mid().format(12), root.box.imag.mid().format(3),
                   root.count);
        if (root.multiplicity) std::println("multiplicity {}", *root.multiplicity);
        else std::println("multiplicity not known");
    }
}

int main() {
    // (x - 1)^3 (x + 2) with exact coefficients: the library proves that 1 is a triple root.
    certpolysolve::Polynomial exact{-2, 5, -3, -1, 1};
    print("(x - 1)^3 (x + 2)", certpolysolve::solve(exact));

    // (x - 0.1)^3 with its coefficients rounded to doubles: no longer a cube, but three distinct roots close to 0.1.
    const double a = 0.1;
    certpolysolve::Polynomial rounded{-a * a * a, 3 * a * a, -3 * a, 1};
    print("(x - 0.1)^3 in doubles, every root alone", certpolysolve::solve(rounded));

    // Often one box and its count is the useful answer: clusters of radius at most 1e-4.
    print("the same, as clusters of radius 1e-4", certpolysolve::solve(rounded, {.goal = certpolysolve::Goal::clusters("1e-4")}));
}
```

```text
(x - 1)^3 (x + 2) (goal met)
  real -2, imaginary 0: count 1, multiplicity 1
  real 1, imaginary 0: count 3, multiplicity 3
(x - 0.1)^3 in doubles, every root alone (goal met)
  real 0.0999997687897, imaginary 4e-07: count 1, multiplicity 1
  real 0.0999997687897, imaginary -4e-07: count 1, multiplicity 1
  real 0.100000462421, imaginary 0: count 1, multiplicity 1
the same, as clusters of radius 1e-4 (goal met)
  real 0.100002955943, imaginary -9.31e-07: count 3, multiplicity not known
```

- With exact coefficients a repeated root is found and its multiplicity proved.
- Rounding the coefficients turns a repeated root into a tight group of distinct roots. The default separates
  them all, which can take a lot of precision for very tight groups.
- `{.goal = certpolysolve::Goal::clusters("r")}` instead accepts boxes of radius at most r, each with its count. A box
  with count 3 says "three roots in here", not that they are one triple root.

## 5. Time limits

```cpp
#include <certpolysolve/certpolysolve.hpp>
#include <chrono>
#include <print>
#include <stop_token>
#include <thread>

using namespace std::chrono_literals;

certpolysolve::Polynomial slow_input() { // degree 1,000; at 1,000 bits its whole solve takes seconds
    certpolysolve::Polynomial p;
    for (uint64_t k = 0; k <= 1000; ++k) p.set(k, static_cast<int>((k * 7919 + 13) % 1000) - 500);
    return p;
}

void report(const certpolysolve::Polynomial &p, const certpolysolve::Options &options) {
    const auto start = std::chrono::steady_clock::now();
    const certpolysolve::Result result = certpolysolve::solve(p, options);
    const std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - start;

    uint64_t roots_in_boxes = 0;
    for (const certpolysolve::Root &root : result.roots()) roots_in_boxes += root.count;
    std::println("stopped ({}) after {:.3f} s: {} boxes holding {} of {} roots, every root covered: {}", result.stop(),
                 elapsed.count(), result.roots().size(), roots_in_boxes, result.degree(), result.complete() ? "yes" : "no");
}

int main() {
    report(slow_input(), {.bits = 1000, .time_limit = 100ms}); // at most 0.1 s of wall-clock time

    // Advanced: a std::stop_token; the solve stops at its next checkpoint after another thread asks.
    std::stop_source source;
    certpolysolve::Options stoppable;
    stoppable.bits = 1000;
    stoppable.stop = source.get_token();
    std::jthread stopper([&] {
        std::this_thread::sleep_for(50ms);
        source.request_stop();
    });
    report(slow_input(), stoppable);
}
```

```text
stopped (time limit) after 0.101 s: 1 boxes holding 1000 of 1000 roots, every root covered: yes
stopped (cancelled) after 0.051 s: 1 boxes holding 1000 of 1000 roots, every root covered: yes
```

- `.time_limit` takes any `std::chrono` duration and bounds one solve's wall-clock time. A solve that reaches it
  stops at its next check of the clock, between its stages, and returns what it has proved so far. The boxes are always correct, but early in a
  solve they may say little: here, one box holding all 1,000 roots.
- Advanced: `.stop` takes a `std::stop_token`, and the solve stops at its next check when its `std::stop_source` is
  asked to, from any thread. Use it for a Cancel button or to stop the rest of a batch.

## 6. When a solve stops early

A solve that does not reach its goal still returns proved boxes. `ok()` is false, `stop()` says why, and `message()`
is a readable account of what the solver did.

| `stop()` | What happened | What to do |
|---|---|---|
| `Stop::none` | The goal is met (`ok()`). | |
| `Stop::precision_cap` | Some roots need more working precision than the cap allows; the default cap is 65,536 bits. | Raise `.max_bits` (the message often names a cap that is enough), or ask for clusters (section 4). |
| `Stop::input_uncertainty` | Coefficients given within bounds are too uncertain to tell the roots further apart. | Narrow the bounds, or ask for clusters of a radius. |
| `Stop::time_limit`, `Stop::cancelled` | The time limit or your `std::stop_token` stopped it (section 5). | Allow more time. |
| `Stop::work_budget`, `Stop::cell_budget`, `Stop::iteration_budget` | A budget you set ran out: `.work_budget` (counted work), `.cell_budget` (live cells of the subdivision search) or `.iteration_budget` (sweeps of the iteration). | Raise it. |
| `Stop::unsupported` | The question cannot be asked of this input, for instance real roots of a polynomial with complex coefficients. | Change the question; the message says what was refused. |
| `Stop::zero_polynomial` | Every coefficient is zero: every number is a root. | |
| `Stop::indeterminate_degree` | The leading coefficient's bounds include zero. | Narrow its bounds, or leave out the term. |
| `Stop::error` | The input or an option was wrong, and nothing was solved. | Read `message()`. |

Mistakes do not throw. A polynomial remembers its first bad value (`p.ok()` is false and `p.status()` says why),
and a solve of a bad input, or with a bad option, returns a result whose `stop()` is `Stop::error`:

```cpp
certpolysolve::Polynomial p{1, "1/0"};
if (!p.ok()) std::println(stderr, "{}", p.status()); // zero denominator

const certpolysolve::Result result = certpolysolve::solve(p);
if (result.stop() == certpolysolve::Stop::error) std::println(stderr, "{}", result.message()); // the same reason
```

The only exception the library throws is `std::bad_alloc`, when memory runs out.

## Where next

- [More examples](examples-cpp.md): building polynomials with arithmetic and from roots, proving a control loop
  stable, the precision cap, coefficients with error bounds, eliminants, polynomials given by a program, secular
  equations, many solves on many threads, and the C++ traps to know.
- [The guide](api.md): every type, option and function, with the C name beside each C++ name.
- [How CertPolySolve works](how-it-works.md): the ideas in two pages.
