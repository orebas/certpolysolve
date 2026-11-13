// Exact coefficients, huge coefficients, and many digits.
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
