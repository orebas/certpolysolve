// When the default is not enough: two roots 2^-100 apart, and the precision cap.
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
