// Two exact square-free factors whose covers overlap, one of them a power composite: p(z) = (z^2 - 2)(z - c)^2 with
// c = round(sqrt(2) 2^128) / 2^128, within 2^-130 of sqrt 2. The factors z^2 - 2 (solved as w - 2 in w = z^2 and
// lifted) and z - c (multiplicity 2) meet beside sqrt 2, and the solve tightens the factor with the larger region,
// the power composite. Until 21 October 2026 its finished test ignored the tighter goal: the factor was never
// stepped, the goal was raised again at no work up to its representable limit, and Auto stopped after 4 work units
// at 85 working bits "at the precision cap of 65536 bits" (for all roots and for the real-only query alike; found by
// replaying the fuzzer's precision stops). Checked: both
// queries meet the accuracy goal, find every root with its multiplicity, and work far below the cap.
#include "exact_winding_api.hpp"
#include <certpolysolve/certpolysolve.hpp>
#include <iostream>

using namespace reference;

namespace {
certpolysolve::Options query() {
    certpolysolve::Options options;
    options.method = certpolysolve::Method::automatic;
    options.goal = certpolysolve::Goal::accuracy();
    options.bits = 53;
    return options;
}
// The working precision the answers may use: far below the default cap of 65,536 bits that the false stop reported,
// and above what separating the two roots 2^-130 apart at a 53-bit goal needs.
constexpr uint64_t working_bits_bound = 1024;

void check(const certpolysolve::Result &result, uint64_t roots, uint64_t real_regions, const char *what) {
    require(result.goal_met() && result.complete() && cps_result_termination(result.get()) == CPS_SUCCESS &&
                cps_result_stop_reason(result.get()) == CPS_STOP_NONE,
            what);
    uint64_t count = 0, real = 0;
    for (const auto &region : result.roots()) {
        count += region.count;
        real += region.realness == certpolysolve::Realness::real;
    }
    require(count == roots, "the answer does not count every root with its multiplicity");
    require(real == real_regions, "the answer's real regions are not the three real roots");
    require(cps_result_work_bits(result.get()) <= working_bits_bound, "the answer worked far above what it needs");
}
} // namespace

int main() {
    certpolysolve::ThreadScope lifetime;
    try {
        const Rational c(Integer("481231938336009023090067544955250113854"), Integer(1) << 128);
        // (z^2 - 2)(z - c)^2 = z^4 - 2c z^3 + (c^2 - 2) z^2 + 4c z - 2c^2.
        const std::vector<Rational> coefficients{-2 * c * c, 4 * c, c * c - 2, -2 * c, Rational(1)};
        certpolysolve::Polynomial input;
        for (size_t j = 0; j < coefficients.size(); ++j)
            input.set(j, certpolysolve::Number(coefficients[j].str()));
        require(input.ok(), "the input recorded a bad value");
        // All roots: -sqrt 2, sqrt 2 and c twice, each region real.
        check(certpolysolve::solve(input, query()), 4, 3, "Auto for all roots did not meet the accuracy goal");
        // The real-only query: the same three real regions.
        check(certpolysolve::real_roots(input, query()), 4, 3, "Auto's real-only query did not meet the accuracy goal");
        std::cout << "power overlap tests passed\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "FAILED: " << error.what() << '\n';
        return 1;
    }
}
