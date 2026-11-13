// The 0.5 C++ interface (include/certpolysolve/certpolysolve.hpp; docs/cpp-api-proposal-20261019.md and the requests that
// followed: polynomial arithmetic, a polynomial type that gets, sets, prints and gives its coefficients, values
// without C++ traps): numbers of every kind, a polynomial built and read every way, braces and parentheses meaning
// the same, printing, exact root ends compared with integers and doubles, every way to solve, every stop reported
// without an exception, programs and secular equations.
#include <certpolysolve/certpolysolve.hpp>
#include <chrono>
#include <cmath>
#include <format>
#include <iostream>
#include <numbers>
#include <stdexcept>
#include <stop_token>
#include <string>
#include <thread>
#include <utility>
#include <vector>

using namespace std::chrono_literals;
using certpolysolve::Number;
using certpolysolve::Polynomial;

namespace {
void require(bool condition, const std::string &message) {
    if (!condition) throw std::runtime_error(message);
}
void numbers() {
    require(Number(3).real() == "3" && Number(-7LL).real() == "-7", "integers");
    require(Number(0.5).real() == "1@-1", "a double is its exact binary value: " + Number(0.5).real());
    require(Number("1/3").real() == "1/3" && Number(std::string("0.1")).real() == "0.1", "text");
    const Number z = std::complex<double>(1.5, -2);
    require(z.real() == "3@-1" && z.imag() == "-1@1" && !z.is_real(), "complex<double>");
    const Number w = Number::complex("1/2", 3);
    require(w.real() == "1/2" && w.imag() == "3", "Number::complex");
}
void polynomials() {
    // braces and parentheses mean the same thing: the constant 3
    require(Polynomial{3} == Polynomial(3) && Polynomial(3).degree() == 0, "{3} and (3)");
    // ascending lists mixing integers, text, doubles; a range; a constant; x and arithmetic
    const Polynomial a{1, "-1/3", 0.5};
    require(a.degree() == 2 && a[1].value.real() == "-1/3" && a[2].value.real() == "1/2", "a list");
    const std::vector<std::string> lines{"2", "0", "-1"};
    require(Polynomial(lines) == Polynomial{2, 0, -1}, "a range");
    const auto x = Polynomial::x();
    require(x * x - 2 == Polynomial{-2, 0, 1}, "x^2 - 2 from arithmetic");
    require(std::format("{}", Polynomial{"1/3", -2, 0, 1}) == "x^3 - 2x + 1/3", "printing: " + std::format("{}", Polynomial{"1/3", -2, 0, 1}));
    require(std::format("{}", Polynomial{"-7/3", "-1/2", 0, 1}) == "x^3 - (1/2)x - 7/3",
            "fractions: " + std::format("{}", Polynomial{"-7/3", "-1/2", 0, 1}));
    require(std::format("{}", Polynomial{0, Number::complex(1, -2), 0, 3}) == "3x^3 + (1 - 2i)x", "complex printing");
    require(std::format("{}", Polynomial{}) == "0" && std::format("{}", -x + 1) == "-x + 1", "zero and signs");
    // Printing never throws: an exponent beyond long long, or text not checked yet, prints as given.
    require(std::format("{}", Number("1@100000000000000000000")) == "1@100000000000000000000" &&
                std::format("{}", Number("1@x")) == "1@x",
            "printing a large or unchecked exponent");
    require(std::format("{}", Polynomial{"3@-2", "1@100000000000000000000"}) == "(1@100000000000000000000)x + 3/4",
            "printing a polynomial with a large exponent: " + std::format("{}", Polynomial{"3@-2", "1@100000000000000000000"}));
    // a number and a coefficient print on their own, exactly
    require(std::format("{}", Number(0.5)) == "1/2" && std::format("{}", Number("-1/3")) == "-1/3" &&
                std::format("{}", Number::complex(1, -2)) == "1 - 2i" && std::format("{}", a[1]) == "-1/3",
            "numbers and coefficients printed");
    // wilk30 three ways: a product loop, from_roots, and the expanded coefficients read back
    Polynomial w = 1;
    for (int k = 1; k <= 30; ++k) w *= x - k;
    std::vector<int> roots;
    for (int k = 1; k <= 30; ++k) roots.push_back(k);
    require(w == Polynomial::from_roots(roots), "a product loop and from_roots agree");
    require(w.degree() == 30 && w[0].value.real() == "265252859812191058636308480000000" && w[30].value.real() == "1",
            "wilk30's coefficients");
    require(Polynomial(w.coefficients()) == w, "coefficients() rebuilds the polynomial");
    // a recurrence: Chebyshev T_5 = 16x^5 - 20x^3 + 5x
    Polynomial t0 = 1, t1 = x;
    for (int k = 1; k < 5; ++k) {
        auto t2 = 2 * x * t1 - t0;
        t0 = std::move(t1);
        t1 = std::move(t2);
    }
    require(t1 == Polynomial{0, 5, 0, -20, 0, 16}, "Chebyshev T_5 by its recurrence");
    require(t1.derivative() == Polynomial{5, 0, -60, 0, 80}, "derivative");
    require(pow(x + 1, 3) == Polynomial{1, 3, 3, 1}, "pow");
    // sparse building and terms
    const auto s = Polynomial::monomial(1000) - 3 * Polynomial::monomial(7) + 1;
    const auto terms = s.terms();
    require(terms.size() == 3 && terms[0].exponent == 0 && terms[1].exponent == 7 && terms[2].exponent == 1000 &&
                terms[1].coefficient.value.real() == "-3",
            "terms");
    // sparse input (C40): a billion-degree polynomial of three terms is built, raised and read back as three terms
    const auto huge = Polynomial::monomial(1000000000) - 3 * Polynomial::monomial(7) + 1;
    require(huge.ok() && huge.degree() == 1000000000 && huge.terms().size() == 3, "a sparse polynomial: " + huge.status());
    const auto squared = pow(Polynomial::monomial(1000000000) + 1, 2);
    require(squared.terms().size() == 3 && squared[1000000000].value.real() == "2", "a sparse square");
    require(huge == 1 - Polynomial::monomial(7, 3) + Polynomial::monomial(1000000000) && !(huge == squared),
            "sparse equality");
    require(Polynomial{"3@-2", "0.5"} == Polynomial{"3/4", "1/2"}, "one value written two ways");
    require(Polynomial::sparse({{1000000000, 1}, {0, 1}, {7, -3}}) == huge, "Polynomial::sparse, in any order");
    const std::vector<Polynomial::Term> listed{{7, -3}, {1000000000, 1}, {0, 1}};
    require(Polynomial::sparse(listed) == huge, "Polynomial::sparse of a range");
    // from_roots with multiplicities and complex roots
    const std::vector<std::pair<Number, int>> with{{1, 2}, {Number::complex(0, 1), 1}, {Number::complex(0, -1), 1}};
    require(Polynomial::from_roots(with) == pow(x - 1, 2) * (x * x + 1), "from_roots with multiplicities");
    // a ball coefficient: kept, printed, refused by exact arithmetic without an exception
    Polynomial b{1, certpolysolve::ball("1/2", "1/1000")};
    require(!b[1].exact() && b[1].real_radius.real() == "1/1000", "a ball");
    require(std::format("{}", b) == "(1/2 +/- 1/1000)x + 1" && std::format("{}", b[1]) == "1/2 +/- 1/1000",
            "ball printing: " + std::format("{}", b));
    const auto bad = b * x;
    require(!bad.ok() && !bad.status().empty(), "arithmetic on a ball reports, not throws");
    // a bad number: remembered, and the solve reports it
    const Polynomial wrong{1, "abc"};
    require(!wrong.ok(), "a bad number is remembered");
    const auto r = certpolysolve::solve(wrong);
    require(!r.ok() && r.stop() == certpolysolve::Stop::error && !r.message().empty(), "a bad input's solve reports it");
}
void solving() {
    const auto x = Polynomial::x();
    // real roots of x^2 - 2: exact ends compared with integers and doubles; printing with digits
    const auto two = certpolysolve::real_roots(x * x - 2, {.bits = 64});
    require(two.ok() && two.roots().size() == 2, "real roots of x^2 - 2: " + two.message());
    // the positive root (the order of regions is not promised yet: the spec's C26), compared exactly with integers
    // and doubles: its 64-bit interval lies strictly between 1.41421356237309 and 1.4142135623731 (sqrt(2)'s
    // nearest double, about 1e-16 away, is outside it, as it should be)
    const auto &root = two.roots()[0].box.real.lo > 0 ? two.roots()[0] : two.roots()[1];
    require(root.box.real.lo > 1 && root.box.real.hi < 2 && root.box.real.lo > 1.41421356237309 &&
                root.box.real.hi < 1.4142135623731 && !root.box.real.contains(1.4142135623730951) &&
                root.realness == certpolysolve::Realness::real,
            "exact comparisons");
    require(std::format("{:.20}", root) == "1.4142135623730950488", "20 digits: " + std::format("{:.20}", root));
    require(root.box.real.width() < certpolysolve::Dyadic(1e-19), "64 bits of accuracy");
    // every complex root, options by plain assignment, a region
    certpolysolve::Options o;
    o.bits = 80;
    const auto all = certpolysolve::solve(pow(x, 2) + 1, o);
    require(all.ok() && all.roots().size() == 2 && all.roots()[0].realness == certpolysolve::Realness::not_real, "x^2 + 1");
    const auto positive = certpolysolve::real_roots(Polynomial::from_roots(std::vector<int>{-2, 1, 3}), certpolysolve::at_least(0));
    require(positive.ok() && positive.roots().size() == 2, "the real roots at least 0");
    // stops, each reported without an exception: the zero polynomial, a time limit, a cancellation from another thread
    require(certpolysolve::solve(Polynomial{}).stop() == certpolysolve::Stop::zero_polynomial, "the zero polynomial");
    Polynomial slow;
    for (uint64_t k = 0; k <= 1000; ++k) slow.set(k, static_cast<int>((k * 7919 + 13) % 1000) - 500);
    const auto limited = certpolysolve::solve(slow, {.bits = 1000, .time_limit = 50ms});
    require(!limited.ok() && limited.stop() == certpolysolve::Stop::time_limit && limited.complete(), "a time limit");
    // After a stop no region has a value, as the C call refuses one (C and C++ never give
    // different values); the boxes and their centres stay.
    require(limited.values().empty(), "no values after a stop");
    for (const certpolysolve::Root &root : limited.roots())
        require(!root.value && std::isfinite(root.centre().real()), "a stopped region has a centre and no value");
    std::stop_source source;
    std::jthread canceller([&] {
        std::this_thread::sleep_for(20ms);
        source.request_stop();
    });
    certpolysolve::Options cancellable;
    cancellable.bits = 1000;
    cancellable.stop = source.get_token();
    const auto cancelled = certpolysolve::solve(slow, cancellable);
    require(cancelled.stop() == certpolysolve::Stop::cancelled, std::format("a std::stop_token: {}", cancelled.stop()));
    // a goal of clusters, and a bad option reported
    const auto clusters = certpolysolve::solve(pow(x - 1, 3), {.goal = certpolysolve::Goal::clusters("1e-6")});
    require(clusters.ok() && clusters.roots().size() == 1 && clusters.roots()[0].count == 3, "a cluster of 3");
    require(certpolysolve::solve(x, {.bits = 0}).stop() == certpolysolve::Stop::error, "a bad option");
    // the advanced controls: an absolute tolerance, and the exact-arithmetic allowances
    certpolysolve::Options tolerant;
    tolerant.absolute_tolerance = "1e-10";
    require(certpolysolve::solve(x * x - 2, tolerant).ok(), "an absolute tolerance");
    tolerant.absolute_tolerance = "-1";
    require(certpolysolve::solve(x * x - 2, tolerant).stop() == certpolysolve::Stop::error, "a negative tolerance accepted");
    certpolysolve::ResourceLimits limits;
    cps_resource_limits defaults;
    cps_resource_limits_default(&defaults);
    require(limits.exact_operations == defaults.exact_operations && limits.flags == defaults.flags,
            "ResourceLimits starts at the defaults");
    certpolysolve::Options allowed;
    allowed.resources = limits;
    require(certpolysolve::solve(x * x - 2, allowed).ok(), "the default allowances");
    limits.flags = UINT64_C(1) << 7; // no such flag
    allowed.resources = limits;
    require(certpolysolve::solve(x * x - 2, allowed).stop() == certpolysolve::Stop::error, "an unknown resource flag accepted");
    // workers: the same roots at one and at four, and 0 reported as a bad option
    const auto one = certpolysolve::solve(slow, {.bits = 100});
    const auto four = certpolysolve::solve(slow, {.bits = 100, .workers = 4});
    require(one.ok() && four.ok() && one.roots().size() == four.roots().size(), "four workers");
    for (size_t i = 0; i < one.roots().size(); ++i)
        require(std::format("{:.40}", one.roots()[i]) == std::format("{:.40}", four.roots()[i]), "the same root at four workers");
    require(certpolysolve::solve(x, {.workers = 0}).stop() == certpolysolve::Stop::error, "zero workers accepted");
}
void other_inputs() {
    // a program: Mandelbrot's p_3 times T_4, constants as plain numbers
    const auto program = certpolysolve::Program::from([](auto x) {
        auto m = x * 0 + 1;
        for (int j = 0; j < 3; ++j) m = x * m * m + 1;
        auto previous = x * 0 + 1, chebyshev = x;
        for (int k = 1; k < 4; ++k) {
            auto next = 2 * x * chebyshev - previous;
            previous = chebyshev;
            chebyshev = next;
        }
        return m * chebyshev - "1/3";
    });
    require(program.ok(), "a program: " + program.status());
    const auto p = certpolysolve::solve(program, {.bits = 64});
    require(p.ok() && p.degree() == 11, "a program solved: " + p.message());
    // a program from its text form (cps_program_parse; C35): (x - 2)^2 + 1, roots 2 +/- i; a malformed text is
    // recorded, not thrown, and its solve is the error
    const certpolysolve::Program parsed("PSPROGRAM 1\nx\nc 2\nsub 0 1\nmul 2 2\nc 1\nadd 3 4\nout 5\n");
    require(parsed.ok(), "a program text was refused: " + parsed.status());
    const auto from_text = certpolysolve::solve(parsed);
    require(from_text.ok() && from_text.degree() == 2 && from_text.roots().size() == 2 &&
                from_text.roots()[0].box.real.contains(2) && from_text.roots()[1].box.real.contains(2),
            "the parsed program's roots: " + from_text.message());
    // pow(value, 0) is the constant 1 and pow(value, 1) the value (the text form's pow needs 2 or more): x (x - 2)
    const auto powers = certpolysolve::Program::from([](auto x) { return pow(x, 1) * pow(x - 5, 0) * (x - 2); });
    require(powers.ok(), "pow with 0 or 1: " + powers.status());
    const auto powered = certpolysolve::solve(powers);
    require(powered.ok() && powered.degree() == 2 && powered.roots().size() == 2 &&
                powered.roots()[0].box.real.contains(0) && powered.roots()[1].box.real.contains(2),
            "pow with 0 or 1: " + powered.message());
    // An operand recorded by another program is refused, not read as this program's node with that number.
    std::optional<certpolysolve::ProgramValue> leaked;
    const auto first = certpolysolve::Program::from([&](auto x) { leaked = x * x; return x - 1; });
    const auto second = certpolysolve::Program::from([&](auto x) { return x * *leaked + 1; });
    require(first.ok() && !second.ok() && second.status().find("another program") != std::string::npos,
            "a foreign operand: " + second.status());
    const certpolysolve::Program malformed("PSPROGRAM 1\nx\nfrobnicate 0\nout 1\n");
    require(!malformed.ok() && !malformed.status().empty(), "a malformed program text was accepted");
    require(certpolysolve::solve(malformed).stop() == certpolysolve::Stop::error, "a malformed program was solved");
    // a secular equation with named fields: roots 1, 2 and 4 (test_d3's)
    const certpolysolve::SecularEquation secular{{.residue = "36/5", .pole = -2}, {.residue = "-4/3", .pole = 0}, {.residue = "2/15", .pole = 3}};
    const auto s = certpolysolve::real_roots(secular);
    require(s.ok() && s.roots().size() == 3 && s.roots()[2].box.real.contains(4), "a secular equation: " + s.message());
}
// C42: a list either way round, named, and the roots as numbers.
void list_forms_and_values() {
    const auto text = [](const Polynomial &p) { return std::format("{}", p); };
    // x^2 + x - 2 is not a palindrome, so a reversed list is another polynomial.
    const Polynomial brace{-2, 1, 1};
    require(text(Polynomial::lowest_first({-2, 1, 1})) == text(brace), "lowest_first is the brace form");
    require(text(Polynomial::highest_first({1, 1, -2})) == text(brace), "highest_first reverses: " + text(Polynomial::highest_first({1, 1, -2})));
    require(text(Polynomial::highest_first({1, 1, -2})) != text(Polynomial::lowest_first({1, 1, -2})), "the two forms differ on an asymmetric list");
    require(text(Polynomial::highest_first({1, 0, -2})) == "x^2 - 2", "x^2 - 2 highest first: " + text(Polynomial::highest_first({1, 0, -2})));
    const auto padded = Polynomial::highest_first({0, 1, 0, -2});
    require(padded.degree() == 2 && text(padded) == "x^2 - 2", "a leading zero drops: " + text(padded));
    require(Polynomial::highest_first({}).is_zero() && Polynomial::lowest_first({}).is_zero(), "an empty list is the zero polynomial");
    const std::vector<Number> numbers{"1", "1", "-2"};
    require(text(Polynomial::highest_first(numbers)) == text(brace), "highest_first from a range");
    const std::vector<int> ints{-2, 1, 1};
    require(text(Polynomial::lowest_first(ints)) == text(brace), "lowest_first from a range");
    // values(): (x - 1)^2 (x - 2) has the double root 1, so 1 appears twice, then 2; doubles from the boxes.
    const auto result = certpolysolve::solve(Polynomial{-2, 5, -4, 1});
    require(result.ok(), "the cubic solved: " + result.message());
    const auto values = result.values();
    require(values.size() == 3, "three values for a cubic, the double root twice: " + std::to_string(values.size()));
    const auto near = [](std::complex<double> z, double x) { return std::abs(z - std::complex<double>(x, 0)) < 1e-9; };
    require(near(values[0], 1) && near(values[1], 1) && near(values[2], 2), "the values are 1, 1, 2");
    // Each value is the C call's, bit for bit, and the centre where there is one.
    for (size_t i = 0; i < result.roots().size(); ++i) {
        double re = 0, im = 0;
        require(cps_result_get_approximation(result.get(), i, &re, &im) == CPS_OK, "the C call gives a value");
        const auto &root = result.roots()[i];
        require(root.value && root.value->real() == re && root.value->imag() == im, "the C++ value is the C value");
        require(root.centre() == *root.value, "the centre and the value of a resolved region");
    }
    require(certpolysolve::solve(Polynomial{}).values().empty(), "no values after an error (the zero polynomial)");
}
} // namespace

int main() {
    certpolysolve::ThreadScope lifetime;
    try {
        numbers();
        polynomials();
        solving();
        other_inputs();
        list_forms_and_values();
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
    std::cout << "C++ interface tests passed\n";
    return 0;
}
