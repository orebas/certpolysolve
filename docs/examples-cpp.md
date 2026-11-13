# CertPolySolve for C++: examples by topic

Complete programs beyond [the quick start](quickstart-cpp.md), each a file in
[`examples/topics/`](../examples/topics/) that the test suite builds and runs; the output under each is what it
printed.

1. [Building polynomials](#1-building-polynomials)
2. [Reading the answer: a proof of stability](#2-reading-the-answer-a-proof-of-stability)
3. [When the default precision is not enough](#3-when-the-default-precision-is-not-enough)
4. [Coefficients with errors](#4-coefficients-with-errors)
5. [Eliminants from Gröbner bases and RUR](#5-eliminants-from-gröbner-bases-and-rur)
6. [Programs: polynomials given by how to compute them](#6-programs-polynomials-given-by-how-to-compute-them)
7. [Secular equations](#7-secular-equations)
8. [Threads](#8-threads)
9. [C++ traps](#9-c-traps)

## 1. Building polynomials

A list of coefficients is one way to give a polynomial. You can also build it the way you would write it: with
exact arithmetic, from its roots, by a recurrence, or term by term. And you can read it back:

```cpp
#include <certpolysolve/certpolysolve.hpp>
#include <print>
#include <vector>

int main() {
    using certpolysolve::Polynomial;
    const Polynomial x = Polynomial::x();

    // Arithmetic is exact: +, -, * and pow, with numbers mixed in. A fraction is text: 1/3 in C++ is 0.
    const Polynomial p = pow(x - 1, 3) * (x + 2) - "1/3";
    std::println("p = {}", p);
    std::println("p' = {}", p.derivative());

    // Wilkinson's polynomial (x - 1)(x - 2)...(x - 20), as a product and from its roots.
    Polynomial w = 1;
    std::vector<int> roots;
    for (int k = 1; k <= 20; ++k) {
        w *= x - k;
        roots.push_back(k);
    }
    std::println("the same both ways: {}", w == Polynomial::from_roots(roots));
    std::println("degree {}; the coefficient of x^19 is {}, the constant term {}", w.degree(), w[19], w[0]);

    // Chebyshev's T_6 by its recurrence T_(k+1) = 2x T_k - T_(k-1).
    Polynomial previous = 1, current = x;
    for (int k = 1; k < 6; ++k) {
        Polynomial next = 2 * x * current - previous;
        previous = current;
        current = next;
    }
    std::println("T_6 = {}", current);
    std::print("its coefficients, constant term first:");
    for (const certpolysolve::Coefficient &c : current.coefficients()) std::print(" {}", c);
    std::println("");

    // Term by term: coefficients you do not set are zero, and only the nonzero terms are stored, at any degree.
    Polynomial sparse;
    sparse.set(1000, 1).set(7, -3).set(0, "1/2");
    std::print("{}: degree {}, terms", sparse, sparse.degree());
    for (const Polynomial::Term &term : sparse.terms()) std::print(" ({}, {})", term.exponent, term.coefficient);
    std::println("");
}
```

```text
p = x^4 - x^3 - 3x^2 + 5x - 7/3
p' = 4x^3 - 3x^2 - 6x + 5
the same both ways: true
degree 20; the coefficient of x^19 is -210, the constant term 2432902008176640000
T_6 = 32x^6 - 48x^4 + 18x^2 - 1
its coefficients, constant term first: -1 0 18 0 -48 0 32
x^1000 - 3x^7 + 1/2: degree 1000, terms (0, 1/2) (7, -3) (1000, 1)
```

- `Polynomial::x()` is x. `+`, `-`, `*`, `pow(p, n)` and `derivative()` are exact (the coefficients are rational,
  or complex with rational parts), and numbers mix in: integers, doubles, and text such as `"1/3"`.
- `Polynomial::from_roots(r)` is the product of (x − r) over a range of numbers; a range of pairs (root,
  multiplicity) works too, and roots may be complex (`certpolysolve::Number::complex(0, 1)` is i).
- Reading: `degree()` (−1 for the zero polynomial), `p[i]` or `p.coefficient(i)` (zero beyond the degree),
  `coefficients()` (constant term first), `terms()` (the nonzero ones, lowest power first) and `==`. A
  polynomial, a coefficient and a number all print with `std::format` and `std::println`, exactly.
- `p.set(i, c)` sets the coefficient of xⁱ and returns `p`, so calls chain. A polynomial stores only its nonzero
  terms: x¹⁰⁰⁰⁰⁰⁰⁰⁰⁰ + 1 is two terms, and arithmetic on few terms stays few terms.
  `Polynomial::sparse({{1000000000, 1}, {0, 1}})` builds one from its terms. A solve still reads every coefficient
  up to the degree.
- A mistake does not throw: a bad number such as `"1/0"` is remembered (`p.ok()` is false and `p.status()` says
  why), arithmetic passes it on, and a solve reports it.

## 2. Reading the answer: a proof of stability

The answer is a list of boxes, `result.roots()`. Each box is proved to contain exactly `count` roots, and together
the boxes contain every root. That makes a proof possible. A control loop is stable when every root of its
characteristic polynomial has a negative real part, and the boxes can show that it is, or that it is not:

```cpp
#include <certpolysolve/certpolysolve.hpp>
#include <print>
#include <string_view>

std::string_view verdict(double gain) {
    certpolysolve::Polynomial p{1 + gain, 3, 3, 1}; // s^3 + 3s^2 + 3s + 1 + K
    certpolysolve::Result result = certpolysolve::solve(p);
    if (!result.ok()) return "no answer";

    // Each box holds roots, and together the boxes hold every root, so this is a proof. Its ends are exact numbers,
    // and comparing them with 0 is exact too.
    bool all_left = true;
    for (const certpolysolve::Root &root : result.roots()) {
        if (root.box.real.lo > 0) return "unstable (a root is proved right of the axis)";
        if (root.box.real.hi >= 0) all_left = false; // this box touches the axis
    }
    return all_left ? "stable (every root is proved left of the axis)" : "undecided (a root box touches the axis)";
}

int main() {
    for (double gain : {7.9, 7.999999, 8.0, 8.1}) std::println("K = {}: {}", gain, verdict(gain));
}
```

```text
K = 7.9: stable (every root is proved left of the axis)
K = 7.999999: stable (every root is proved left of the axis)
K = 8: undecided (a root box touches the axis)
K = 8.1: unstable (a root is proved right of the axis)
```

At K = 8 two roots lie exactly on the axis, at ±i√3, so no box can be on either side of it.

What a box holds (`certpolysolve::Root`):

| Field | Meaning |
|---|---|
| `box` | The box: `box.real` and `box.imag`, intervals with exact ends `lo` and `hi` (`certpolysolve::Dyadic`, numbers m · 2^e). Ends compare exactly with each other, with integers and with doubles; `to_double(Rounding::down)` and `to_double(Rounding::up)` round them outward; `to_string()` writes one as an exact fraction and `format(n)` with n digits. An interval's `mid()` and `width()` are exact too. |
| `value` | The root as a `std::complex<double>`, where the answer gives one (a `std::optional`): set when the result is complete with its goal met and the box holds one root of proved multiplicity; empty after a stop, for an unresolved cluster and for an isolation answer. The C call `cps_result_get_approximation` gives the same number. |
| `centre()` | The box's centre as a `std::complex<double>`, for any box: what printing shows. Where `value` is set, the two are equal. |
| `count` | How many roots are in the box, counted with multiplicity. Always proved. |
| `multiplicity` | A `std::optional`: for a box proved to hold one root, its multiplicity; empty when the box may hold several different roots or the multiplicity is not known. |
| `realness` | `Realness::real` when the roots are proved real, `Realness::not_real` when proved not real, `Realness::unknown` otherwise. |

A root prints as its centre: `{}` to double precision, `{:.30}` with 30 digits from the exact box. About the whole
result: `ok()`, `stop()`, `complete()` (the boxes hold every root; false only when a solve stopped before it had
proved a cover, which a time limit early in a solve can do), `roots()`, `values()` (the roots as
`std::complex<double>`, repeated by count; empty after a stop), `degree()` and `message()`.

## 3. When the default precision is not enough

To tell roots apart, the solver works at higher precision as needed, up to a cap of 65,536 bits. A root that needs
more stops the solve:

```cpp
#include <certpolysolve/certpolysolve.hpp>
#include <print>
#include <ranges>
#include <string_view>

void solve_with_cap(const certpolysolve::Polynomial &p, uint64_t cap) {
    // 53 bits of accuracy, using at most `cap` bits of working precision
    const certpolysolve::Result result = certpolysolve::solve(p, {.bits = 53, .max_bits = cap});

    if (result.ok()) std::print("cap {}: goal met; boxes with counts", cap);
    else std::print("cap {}: stopped ({}); boxes with counts", cap, result.stop());
    for (const certpolysolve::Root &root : result.roots()) std::print(" {}", root.count);
    std::println("");
    for (auto part : std::views::split(result.message(), ';')) { // a stop at the cap may say what cap would do
        const std::string_view clause(part.begin(), part.end());
        if (clause.contains("suffices")) std::println(" {}", clause);
    }
}

int main() {
    // (x - 1)(x - 1 - 2^-100) = x^2 - (2 + 2^-100) x + (1 + 2^-100), written exactly: "m@e" is m * 2^e.
    certpolysolve::Polynomial p{"1267650600228229401496703205377@-100",  // (2^100 + 1) * 2^-100
                            "-2535301200456458802993406410753@-100", // -(2^101 + 1) * 2^-100
                            1};
    solve_with_cap(p, 96);                           // too little to tell the two roots apart
    solve_with_cap(p, certpolysolve::Options{}.max_bits); // the default cap, 65,536 bits
}
```

```text
cap 96: stopped (precision cap); boxes with counts 2
  by Mahler's bound the roots are more than 2^-105 apart and less than 2^3 in size, so a cap of 322 bits suffices to isolate them (cps_options_set_max_work_bits)
cap 65536: goal met; boxes with counts 1 1
```

- `.bits` is the accuracy you want, `.max_bits` the most working precision the solver may use.
- A solve that stops at the cap still returns proved boxes: here one box with both roots in it. To go further,
  solve again with a larger cap (the message often names one that is enough, from Mahler's bound: a lower bound on
  the distance between distinct roots, computed from the coefficients), or accept the box as a cluster.

## 4. Coefficients with errors

When a coefficient is known only within bounds, from a measurement or a floating-point computation with an error
bound, give it as a value and a radius:

```cpp
#include <certpolysolve/certpolysolve.hpp>
#include <print>

int main() {
    using certpolysolve::ball;
    certpolysolve::Polynomial p{ball(2, "1e-9"), ball(-3, "1e-9"), ball(1, "1e-9")}; // value, then radius
    std::println("{}", p);

    // The boxes hold the roots of EVERY polynomial within those bounds, so ask for boxes of a size, not for digits.
    certpolysolve::Result result = certpolysolve::solve(p, {.goal = certpolysolve::Goal::clusters("1e-6")});
    std::println("{}", result.ok() ? "goal met" : "stopped");
    for (const certpolysolve::Root &root : result.roots())
        std::println("  real part in {:.12}, count {}", root.box.real, root.count);
}
```

```text
(1 +/- 1/1000000000)x^2 + (-3 +/- 1/1000000000)x + (2 +/- 1/1000000000)
goal met
  real part in [0.999999988686, 1.00000001131], count 1
  real part in [1.99999996889, 2.00000003111], count 1
```

- The answer is about the whole family of polynomials within the bounds: each box holds `count` roots of every
  one of them.
- `ball(value, radius)` bounds each part of a complex value by the radius; `ball(value, real_radius, imag_radius)`
  bounds them separately. The text `"1e-9"` is exactly 1/1,000,000,000, and prints so.
- If the leading coefficient's bounds include zero, the degree itself is unknown and the solve stops with
  `Stop::indeterminate_degree`.

## 5. Eliminants from Gröbner bases and RUR

An eliminant, the univariate polynomial a Gröbner basis or a rational univariate representation (RUR) of a
polynomial system hands over, arrives as integer coefficients of any size. Read them as text, then ask for the real roots:

```cpp
#include <certpolysolve/certpolysolve.hpp>
#include <chrono>
#include <print>
#include <ranges>
#include <string>
#include <string_view>
#include <vector>

using namespace std::chrono_literals;

std::vector<std::string> read_coefficients(std::string_view text) {
    std::vector<std::string> coefficients;
    for (auto line : std::views::split(text, '\n'))
        if (!line.empty()) coefficients.emplace_back(line.begin(), line.end());
    return coefficients;
}

int main() {
    // Stands for the text of a file. This one is (10^30 x - 7) times the minimal polynomial of
    // sqrt(2) + sqrt(3) + sqrt(5): nine real roots, one of them 7 * 10^-30.
    const std::string_view file = "-4032\n"
                                  "576000000000000000000000000000000\n"
                                  "6720\n"
                                  "-960000000000000000000000000000000\n"
                                  "-2464\n"
                                  "352000000000000000000000000000000\n"
                                  "280\n"
                                  "-40000000000000000000000000000000\n"
                                  "-7\n"
                                  "1000000000000000000000000000000\n";
    certpolysolve::Polynomial f(read_coefficients(file));

    // each real root to 64 relative bits; give up after 60 s
    const certpolysolve::Result result = certpolysolve::real_roots(f, {.bits = 64, .time_limit = 60s});
    if (!result.ok()) {
        std::println(stderr, "stopped: {}\n{}", result.stop(), result.message());
        return 1;
    }

    for (const certpolysolve::Root &root : result.roots()) {
        std::println("{:.20} (multiplicity {})", root, root.multiplicity.value_or(0)); // 0: not known
        std::println("  in [{},\n      {}]", root.box.real.lo.to_string(), root.box.real.hi.to_string()); // exact
    }

    // Only the intervals, each holding one real root, without refining them: the cheapest real query.
    std::println("{} real roots isolated", certpolysolve::real_roots(f, {.goal = certpolysolve::Goal::isolation()}).roots().size());

    // Every complex root, for the complex solutions of the system.
    std::println("{} complex roots", certpolysolve::solve(f).roots().size());
}
```

```text
-5.3823323474417620387 (multiplicity 1)
  in [-426432301928897787736773608925/79228162514264337593543950336,
      -426432301928897787735328800083/79228162514264337593543950336]
-2.5539052226955719411 (multiplicity 1)
  in [-404682436059506456122932461907/158456325028528675187087900672,
      -404682436059506456121561344477/158456325028528675187087900672]
-1.9182307323040074517 (multiplicity 1)
  in [-151977896198838192531337441303/79228162514264337593543950336,
      -151977896198838192530822520161/79228162514264337593543950336]
-0.91019639244218264592 (multiplicity 1)
  in [-72113187700306367143969936675/79228162514264337593543950336,
      -72113187700306367143725607691/79228162514264337593543950336]
7e-30 (multiplicity 1)
  in [43939712147706765346776089193/6277101735386680763835789423207666416102355444464034512896,
      43939712147706765346924962729/6277101735386680763835789423207666416102355444464034512896]
0.91019639244218264592 (multiplicity 1)
  in [36056593850153183571862803845/39614081257132168796771975168,
      36056593850153183571984968337/39614081257132168796771975168]
1.9182307323040074517 (multiplicity 1)
  in [151977896198838192530822520161/79228162514264337593543950336,
      151977896198838192531337441303/79228162514264337593543950336]
2.5539052226955719411 (multiplicity 1)
  in [404682436059506456121561344477/158456325028528675187087900672,
      404682436059506456122932461907/158456325028528675187087900672]
5.3823323474417620387 (multiplicity 1)
  in [426432301928897787735328800619/79228162514264337593543950336,
      426432301928897787736773609461/79228162514264337593543950336]
9 real roots isolated
9 complex roots
```

- Accuracy is relative: the tiny root 7 · 10⁻³⁰ gets 64 correct bits too.
- The input need not be square-free: a repeated factor shows up as a root with `multiplicity` above 1.
- `to_string()` writes an interval end exactly, a fraction with a power of two as denominator, ready for exact
  arithmetic downstream. With GMP, include [`<certpolysolve/gmp.hpp>`](../include/certpolysolve/gmp.hpp) and link `gmpxx`:
  `certpolysolve::to_mpq(root.box.real.lo)` is the end as an `mpq_class`, and `mpz_class` and `mpq_class` values can be
  coefficients.
- More bits later means a new solve with a higher `.bits`; nothing is kept between solves.

## 6. Programs: polynomials given by how to compute them

Some polynomials are easy to compute but enormous to write out: iterations, recurrences, products. You can give
CertPolySolve the computation instead of the coefficients. Write it as an ordinary C++ function of `x`:

```cpp
auto mandelbrot = certpolysolve::Program::from([](auto c) {
    auto p = c;
    for (int k = 1; k < 8; ++k) p = p * p + c;
    return p;
});
certpolysolve::Result result = certpolysolve::solve(mandelbrot);
```

`Program::from` runs your function once and records the arithmetic it does. That is all a program is. The
recorded steps are what the solver evaluates, so the coefficients are never formed: this polynomial has 129 of
them, large integers, and with 21 iterations it would have over a million.

More examples, all in one file:

```cpp
#include <certpolysolve/certpolysolve.hpp>
#include <algorithm>
#include <cmath>
#include <numbers>
#include <print>
#include <string_view>
#include <vector>

void summary(std::string_view title, const certpolysolve::Result &result) {
    std::println("{}: degree {}, {} boxes, {}", title, result.degree(), result.roots().size(),
                 result.ok() ? "goal met" : "stopped");
}

int main() {
    // The Mandelbrot polynomial p_8(c): p_1 = c, p_(k+1) = p_k^2 + c. Degree 128.
    auto mandelbrot = certpolysolve::Program::from([](auto c) {
        auto p = c;
        for (int k = 1; k < 8; ++k) p = p * p + c;
        return p;
    });
    summary("Mandelbrot p_8", certpolysolve::solve(mandelbrot));

    // Chebyshev's T_50 by its recurrence T_(k+1) = 2x T_k - T_(k-1).
    auto chebyshev = certpolysolve::Program::from([](auto x) {
        auto previous = x, current = 2 * x * x - 1; // T_1 and T_2
        for (int k = 2; k < 50; ++k) {
            auto next = 2 * x * current - previous;
            previous = current;
            current = next;
        }
        return current;
    });
    summary("Chebyshev T_50", certpolysolve::solve(chebyshev));

    // A tridiagonal matrix's characteristic polynomial det(T - x I), by the three-term recurrence. Here the 1D
    // Laplacian of size 40 (2 on the diagonal, -1 beside it), whose eigenvalues are 2 - 2 cos(k pi / 41).
    const int n = 40;
    const std::vector<double> diagonal(n, 2.0), beside(n - 1, -1.0);
    auto laplacian = certpolysolve::Program::from([&](auto x) {
        auto previous = diagonal[0] - x;                                     // the leading 1 x 1 minor
        auto current = (diagonal[1] - x) * previous - beside[0] * beside[0]; // the leading 2 x 2 minor
        for (int k = 2; k < n; ++k) {
            auto next = (diagonal[k] - x) * current - beside[k - 1] * beside[k - 1] * previous;
            previous = current;
            current = next;
        }
        return current;
    });
    const certpolysolve::Result eigenvalues = certpolysolve::solve(laplacian);
    summary("Laplacian of size 40", eigenvalues);
    std::vector<double> computed;
    for (const certpolysolve::Root &root : eigenvalues.roots())
        if (root.value) computed.push_back(root.value->real());
    std::ranges::sort(computed);
    double error = 0;
    for (int k = 1; k <= n; ++k)
        error = std::max(error, std::abs(computed[k - 1] - (2 - 2 * std::cos(k * std::numbers::pi / (n + 1)))));
    std::println("  largest difference from 2 - 2 cos(k pi / 41): {:.3g}", error);

    // Wilkinson's polynomial (x - 1)(x - 2)...(x - 30) as a product: its coefficients, up to 30!, never appear.
    auto wilkinson = certpolysolve::Program::from([](auto x) {
        auto p = x - 1;
        for (int k = 2; k <= 30; ++k) p *= x - k;
        return p;
    });
    // every root in its own box (a program's default allows close roots to share one)
    summary("Wilkinson W_30", certpolysolve::solve(wilkinson, {.goal = certpolysolve::Goal::accuracy()}));

    // Powers and exact fractions: (x^2 + x + 1)^5 - 1/3.
    auto power = certpolysolve::Program::from([](auto x) { return pow(x * x + x + 1, 5) - "1/3"; });
    summary("(x^2 + x + 1)^5 - 1/3", certpolysolve::solve(power));
}
```

```text
Mandelbrot p_8: degree 128, 128 boxes, goal met
Chebyshev T_50: degree 50, 50 boxes, goal met
Laplacian of size 40: degree 40, 40 boxes, goal met
  largest difference from 2 - 2 cos(k pi / 41): 4.44e-16
Wilkinson W_30: degree 30, 30 boxes, goal met
(x^2 + x + 1)^5 - 1/3: degree 10, 10 boxes, goal met
```

What you can write inside the function:

- `+`, `-`, `*`, `+=`, `-=`, `*=`, and `pow(value, n)` for a whole power n.
- Numbers mixed in freely: integers, doubles at their exact value, and text for an exact fraction or decimal, as
  in `- "1/3"` above (complex: `certpolysolve::Number::complex("1", "-2")`).
- Loops, `if`s on loop counters, helper functions, values captured from outside: anything that does not depend on
  the value of x. There is no `<` or `==` on x, so code that would branch on it does not compile. There is no
  division by x either; the result must be a polynomial.

What is different from coefficient input:

- The roots are all the complex roots: `real_roots` does not take a program. Take `root.value->real()` when you know
  the roots are real, as in the Laplacian example (`value` is set once the goal is met).
- By default a box may hold a few very close roots when its centre approximates all of them to the requested
  accuracy. Ask for `{.goal = certpolysolve::Goal::accuracy()}` to have every root alone, as in the Wilkinson example.

## 7. Secular equations

The eigenvalues of a diagonal matrix plus a rank-one change, D + ρ z zᵀ, are the roots of a secular equation,
Σ aⱼ / (x − bⱼ) = 1 with aⱼ = ρ zⱼ² and bⱼ = dⱼ. CertPolySolve takes it in that form:

```cpp
#include <certpolysolve/certpolysolve.hpp>
#include <print>
#include <string>
#include <vector>

int main() {
    const std::vector<std::string> d = {"1", "2", "3", "4"};
    const std::vector<std::string> a = {"1", "1/4", "1/9", "1/16"}; // rho z_j^2

    std::vector<certpolysolve::SecularEquation::Term> terms;
    for (size_t j = 0; j < d.size(); ++j) terms.push_back({.residue = a[j], .pole = d[j]});
    const certpolysolve::SecularEquation s(terms);

    const certpolysolve::Result result = certpolysolve::real_roots(s);
    std::println("{}", result.ok() ? "goal met" : "stopped");
    for (const certpolysolve::Root &root : result.roots()) std::println("  {:.17}", root);
}
```

```text
goal met
  1.5855831027681075
  2.4513655120803733
  3.2524927312691694
  4.1341697649934609
```

Each term is a residue aⱼ and a pole bⱼ, numbers as anywhere else (complex ones with `Number::complex`). The poles
must be distinct and the residues nonzero. A short equation can be written in place:
`certpolysolve::SecularEquation s{{.residue = 1, .pole = 1}, {.residue = "1/4", .pole = 2}};`. `solve(s)` and
`real_roots(s)` work as for polynomials.

## 8. Threads

Solves can run at the same time on any number of threads:

```cpp
#include <certpolysolve/certpolysolve.hpp>
#include <atomic>
#include <print>
#include <thread>
#include <vector>

int verdict(double gain) { // 1 proved stable, -1 proved unstable, 0 undecided
    const certpolysolve::Result result = certpolysolve::solve(certpolysolve::Polynomial{1 + gain, 3, 3, 1});
    if (!result.ok()) return 0;
    bool all_left = true;
    for (const certpolysolve::Root &root : result.roots()) {
        if (root.box.real.lo > 0) return -1;
        if (root.box.real.hi >= 0) all_left = false;
    }
    return all_left ? 1 : 0;
}

int main() {
    const int gains = 400, threads = 4;
    const double first = 6, step = 0.01; // K = 6.00, 6.01, ..., 9.99
    std::vector<int> verdicts(gains);
    std::atomic<int> next{0};

    {
        std::vector<std::jthread> pool;
        for (int t = 0; t < threads; ++t)
            pool.emplace_back([&] {
                certpolysolve::ThreadScope scope; // first thing in the thread: frees its FLINT caches when it ends
                for (int i; (i = next++) < gains;) verdicts[i] = verdict(first + step * i);
            });
    } // the jthreads join here

    for (int i = 0; i < gains; ++i)
        if (verdicts[i] != 1) {
            std::println("stable proved up to K = {:.2f}; at K = {:.2f} {}", first + step * (i - 1), first + step * i,
                         verdicts[i] < 0 ? "unstable" : "undecided");
            break;
        }
}
```

```text
stable proved up to K = 7.99; at K = 8.00 undecided
```

- Polynomials, options and results may be shared between threads for reading, and moved to and destroyed on any
  thread. Do not change an object while another thread uses it.
- Put a `certpolysolve::ThreadScope` first in every thread you start. FLINT, which the library uses for arithmetic,
  keeps caches per thread; the scope frees them when the thread ends. The main thread does not need one.
- One solve can also use several threads of its own: `{.workers = 4}` splits some of its loops across four, with
  the same answer at every count. One worker is the default, and then a solve starts no thread.

## 9. C++ traps

C++ has a few traps for numbers. The library removes the ones it can, and these are the others:

- **`1/3` is 0.** Integer division happens before the library sees the number. Write the fraction as text, `"1/3"`.
  (`1.0/3` is the double nearest one third, which is not one third.)
- **`0.1` is not one tenth.** A double literal is the nearest binary number, and the library uses that exact value.
  Write `"0.1"` or `"1/10"` for one tenth.
- **Large integers overflow before the library sees them.** `1000000 * 1000000` overflows an `int`, and 30! fits no
  integer type. Write big numbers as text, or build them with polynomial arithmetic, which is exact.
- **`certpolysolve::Polynomial p(certpolysolve::Number(x));` declares a function** (C++'s "most vexing parse"). Write
  `auto p = certpolysolve::Polynomial(x);` or use braces. Braces and parentheses mean the same thing for every type
  here: `Polynomial{3}` and `Polynomial(3)` are both the constant 3 (unlike `std::vector<int>{3}` and
  `std::vector<int>(3)`).
- **Designated initializers go in declaration order.** `{.bits = 64, .time_limit = 60s}` compiles and
  `{.time_limit = 60s, .bits = 64}` does not. The order is `bits`, `max_bits`, `time_limit`, `goal`, `method` (forces
  one of the three methods; the default lets the solver choose), `stop`, then `work_budget`, `cell_budget`,
  `iteration_budget`, `absolute_tolerance`, `resources` and `workers`. Fields can also be set one at a time: `certpolysolve::Options o;
  o.time_limit = 60s;`.

Two traps the library closes: a character such as `'3'` (which C++ reads as 51) or a `bool` is refused as a number
at compile time, and nothing in the interface needs a `new`, a `delete` or a `free`.
