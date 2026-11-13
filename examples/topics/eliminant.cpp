// An eliminant, as a Groebner basis or a rational univariate representation hands it over: integer coefficients
// of any size, one a line, constant term first. Its real roots in intervals with exact rational ends.
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
