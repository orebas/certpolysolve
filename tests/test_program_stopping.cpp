#include "exact_reference_api.hpp"
#include <iostream>
#include <stdexcept>

namespace {
void require(bool condition, const char *message) {
    if (!condition) throw std::runtime_error(message);
}

// Fuzz-9 seed 17646018995326954843. Preserve the recorded operations: expanding
// or collecting the factors first would exercise a different solver route.
const char *repeated_program =
    "PSPROGRAM 1\nx\nc 29605/65536 -373/652\nsub 0 1\n"
    "c 29605/65536 -373/652\nsub 0 3\nmul 2 4\npow 5 2\n"
    "c 29605/65536 -373/652\nsub 0 7\npow 8 5\n"
    "c -75/56 421107/524288\nsub 0 10\npow 11 2\nmul 6 9\nmul 13 12\nout 14\n";

void check_cover(const certpolysolve::Result &result) {
    using reference::Rational;
    const std::vector<std::pair<reference::Complex, uint64_t>> known{
        {{Rational(29605) / 65536, -Rational(373) / 652}, 9},
        {{-Rational(75) / 56, Rational(421107) / 524288}, 2}};
    require(result.get() && result.complete(), "stopped program lost its complete cover");
    uint64_t total = 0;
    for (const auto &root : result.roots()) {
        uint64_t expected = 0;
        for (const auto &[point, multiplicity] : known)
            if (reference::endpoint(root.box.real.lo) <= point.real &&
                point.real <= reference::endpoint(root.box.real.hi) &&
                reference::endpoint(root.box.imag.lo) <= point.imag &&
                point.imag <= reference::endpoint(root.box.imag.hi))
                expected += multiplicity;
        require(root.count == expected, "program cover count disagrees with exact planted roots");
        total += root.count;
    }
    require(total == 11, "program cover lost a planted root");
}
} // namespace

int main() {
    certpolysolve::ThreadScope lifetime;
    try {
        const certpolysolve::Program program(repeated_program);
        require(program.ok(), "program fixture refused");
        certpolysolve::Options options;
        options.method = certpolysolve::Method::secular;
        options.goal = certpolysolve::Goal::approximation();
        options.bits = 100;
        options.max_bits = 72;
        options.work_budget = 512;
        const auto limited = certpolysolve::solve(program, options);
        check_cover(limited);
        require(!limited.goal_met() && limited.stop() == certpolysolve::Stop::precision_cap,
                "exhausted stationary program model repeats until the work budget");
        require(limited.statistics().work_units < *options.work_budget,
                "precision stop consumed the entire work allowance");

        options.max_bits = 144;
        const auto raised = certpolysolve::solve(program, options);
        check_cover(raised);
        require(raised.goal_met(), "raising the cap no longer resolves the repeated program");
        for (const auto &root : raised.roots()) reference::check_accuracy(root, 100);

        // Requested accuracy above the working cap is not a blanket refusal:
        // an exact root can still meet it, without an approximation error.
        const certpolysolve::Program exact("PSPROGRAM 1\nx\nout 0\n");
        options.bits = 256;
        options.max_bits = 64;
        const auto exact_result = certpolysolve::solve(exact, options);
        const auto exact_roots = exact_result.roots();
        require(exact_result.complete() && exact_result.goal_met() && exact_roots.size() == 1,
                "precision-cap handling rejected an exact root");
        const auto &zero = exact_roots[0];
        require(zero.count == 1 && reference::endpoint(zero.box.real.lo) == 0 &&
                    reference::endpoint(zero.box.real.hi) == 0 && reference::endpoint(zero.box.imag.lo) == 0 &&
                    reference::endpoint(zero.box.imag.hi) == 0,
                "exact-root control did not return the exact zero");
        std::cout << "stationary model: " << limited.statistics().work_units
                  << " work units; raised cap: " << raised.statistics().work_units << "; exact covers checked\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
