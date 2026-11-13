// word_binomial (chart_certificate.hpp) against fmpz_bin_uiui for every n up to word_binomial_degree and every m up to
// n + 1, and the remainder loop's recurrence C(k - 1, j) = C(k, j) (k - j) / k in 128 bits against fmpz's; also that
// word_binomial_degree is the largest degree whose binomials all fit in 63 bits.
#include "certify/chart_certificate.hpp"
#include "core/raii.hpp"
#include <flint/fmpz.h>
#include <iostream>
#include <stdexcept>
#include <string>

namespace d = certpolysolve::detail;

namespace {
void require(bool ok, const std::string &why) {
    if (!ok) throw std::runtime_error(why);
}
bool equal(const fmpz_t x, d::uint128 y) {
    d::Integer hi, lo;
    fmpz_set_ui(hi.value(), static_cast<ulong>(y >> 64));
    fmpz_mul_2exp(hi.value(), hi.value(), 64);
    fmpz_set_ui(lo.value(), static_cast<ulong>(y));
    fmpz_add(hi.value(), hi.value(), lo.value());
    return fmpz_equal(x, hi.value());
}
} // namespace

int main() {
    try {
        uint64_t checks = 0;
        d::Integer reference;
        const ulong top = static_cast<ulong>(d::word_binomial_degree);
        for (ulong n = 0; n <= top; ++n)
            for (ulong m = 0; m <= n + 1; ++m) {
                fmpz_bin_uiui(reference.value(), n, m);
                require(equal(reference.value(), d::word_binomial(n, m)), "C(" + std::to_string(n) + ", " +
                                                                              std::to_string(m) + ")");
                require(fmpz_bits(reference.value()) <= 63, "a binomial above 63 bits at n " + std::to_string(n));
                ++checks;
                // The remainder loop from k = n down to m: C(k - 1, m) from C(k, m).
                d::uint128 word = d::word_binomial(n, m);
                d::Integer exact;
                fmpz_bin_uiui(exact.value(), n, m);
                for (ulong k = n; k > m && m >= 1; --k) {
                    word = word * (k - m) / k;
                    fmpz_mul_ui(exact.value(), exact.value(), k - m);
                    fmpz_divexact_ui(exact.value(), exact.value(), k);
                    require(equal(exact.value(), word), "the recurrence at n " + std::to_string(n) + ", k " +
                                                            std::to_string(k) + ", m " + std::to_string(m));
                    ++checks;
                }
            }
        fmpz_bin_uiui(reference.value(), top + 1, (top + 1) / 2);
        require(fmpz_bits(reference.value()) > 63, "word_binomial_degree is not the largest that fits");
        std::cout << "word binomial: " << checks << " checks ok\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "test_word_binomial: " << error.what() << '\n';
        return 1;
    }
}
