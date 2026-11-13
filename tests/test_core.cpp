#include "api_test_support.hpp"
#include <certpolysolve/certpolysolve.hpp>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

#define CHECK(test) do { if (!(test)) throw std::runtime_error("FAIL line " + std::to_string(__LINE__) + ": " #test); } while (0)
using namespace certpolysolve;

int main() {
    certpolysolve::ThreadScope thread_lifetime;
    try {
        auto simple = solve(Polynomial{-2, 0, 1});
        CHECK(simple.complete() && simple.goal_met() && simple.roots().size() == 2);
        CHECK(simple.roots()[0].value && std::abs(simple.roots()[0].value->real() + std::sqrt(2.0)) < 1e-15);
        CHECK(simple.roots()[0].value->imag() == 0);
        auto repeated = solve(Polynomial{1, -2, 1});
        CHECK(repeated.goal_met() && repeated.roots().size() == 1);
        CHECK(repeated.roots()[0].multiplicity.value_or(0) == 2);
        // Each root counted as often as its multiplicity (what approximations(true) expanded to).
        uint64_t expanded = 0;
        for (const auto &root : repeated.roots()) expanded += root.multiplicity.value_or(0);
        CHECK(repeated.complete() && repeated.goal_met() && expanded == 2);
        auto zero = solve(Polynomial{});
        CHECK(cps_result_termination(zero.get()) == CPS_ZERO_POLYNOMIAL && !zero.complete());
        auto constant = solve(Polynomial{2, 0, 0});
        CHECK(constant.degree() == 0 && constant.goal_met() && constant.roots().size() == 0);
        auto monomial = solve(Polynomial{0, 0, 0, 7});
        CHECK(monomial.goal_met() && monomial.roots()[0].multiplicity.value_or(0) == 3);
        auto complex_roots = solve(Polynomial{1, 0, 2, 0, 1});
        CHECK(complex_roots.goal_met() && complex_roots.roots().size() == 2);
        CHECK(complex_roots.roots()[0].multiplicity.value_or(0) == 2);
        CHECK(complex_roots.roots()[0].realness == Realness::not_real);
        auto real_none = real_roots(Polynomial{1, 0, 2, 0, 1});
        CHECK(real_none.complete() && real_none.goal_met() && real_none.roots().size() == 0);
        auto real_repeated = real_roots(Polynomial{1, -2, 1});
        CHECK(real_repeated.goal_met() && real_repeated.roots()[0].multiplicity.value_or(0) == 2);
        Options high_precision;
        high_precision.bits = 1024;
        auto high = solve(Polynomial{-2, 0, 1}, high_precision);
        CHECK(high.goal_met());
        CHECK(high.roots()[0].box.real.lo.mantissa().size() > 300);
        // The default precision cap is 65,536 bits (CPS_DEFAULT_MAX_WORK_BITS;
        // 4,096 before 1 October 2026): a request for more working bits than
        // the old cap, under default options, is honoured rather than cut off.
        CHECK(CPS_DEFAULT_MAX_WORK_BITS == 65536);
        Options deep;
        deep.bits = 8192; // the cap is left at its default
        auto deep_roots = solve(Polynomial{-2, 0, 1}, deep);
        CHECK(deep_roots.goal_met() && deep_roots.work_bits() >= 8192 &&
              deep_roots.roots()[0].box.real.lo.mantissa().size() > 2400);
        Options low_precision;
        low_precision.bits = 256;
        low_precision.max_bits = 32;
        auto low = solve(Polynomial{-2, 0, 1}, low_precision);
        CHECK(low.complete() && !low.goal_met() && cps_result_termination(low.get()) == CPS_PRECISION_LIMIT);
        // An unmet goal has no approximation: the C getter refuses it.
        double approximation_real = 0, approximation_imag = 0;
        CHECK(cps_result_get_approximation(low.get(), 0, &approximation_real, &approximation_imag) == CPS_UNRESOLVED);
        Polynomial uncertain;
        uncertain.set(0, Number("-1")).set(1, ball(Number("0"), Number("1")));
        auto unknown_degree = solve(uncertain, low_precision);
        CHECK(cps_result_termination(unknown_degree.get()) == CPS_INDETERMINATE_DEGREE && !unknown_degree.complete());
        Polynomial wide;
        wide.set(0, ball(Number("-1"), Number("1/10"))).set(2, Number("1"));
        auto wide_roots = solve(wide, low_precision);
        CHECK(wide_roots.complete() && !wide_roots.goal_met());
        Polynomial huge;
        huge.set(0, Number("-1@10000000")).set(1, Number("1"));
        auto huge_root = solve(huge);
        CHECK(huge_root.goal_met() && huge_root.roots().size() == 1);
        CHECK(huge_root.roots()[0].box.real.lo.exponent() == "10000000");
        CHECK(cps_result_get_approximation(huge_root.get(), 0, &approximation_real, &approximation_imag) == CPS_OUT_OF_RANGE);
        Polynomial tiny;
        tiny.set(0, Number("-1@-10000000")).set(1, Number("1"));
        auto tiny_root = solve(tiny);
        CHECK(tiny_root.goal_met());
        CHECK(tiny_root.roots()[0].box.real.lo.exponent() == "-10000000");
        Polynomial subnormal;
        subnormal.set(0, Number("-3@-1075")).set(1, Number("1"));
        const auto subnormal_root = solve(subnormal);
        CHECK(subnormal_root.goal_met());
        CHECK(subnormal_root.roots()[0].value && subnormal_root.roots()[0].value->real() == 2 * std::numeric_limits<double>::denorm_min());
        // The exact 1.5-subnormal-unit root remains unchanged by the lossy
        // ties-to-even binary64 value above.
        CHECK(subnormal_root.roots()[0].box.real.lo.mantissa() == "3");
        CHECK(subnormal_root.roots()[0].box.real.lo.exponent() == "-1075");
        subnormal.set(0, Number("-1@-1076"));
        const auto underflow = solve(subnormal);
        CHECK(underflow.goal_met());
        CHECK(cps_result_get_approximation(underflow.get(), 0, &approximation_real, &approximation_imag) == CPS_OUT_OF_RANGE);
        auto ball_only = certpolysolve_test::raw_options();
        CHECK(certpolysolve_test::set_exact_allowance(ball_only.get(), 0) == CPS_OK);
        CHECK(cps_options_set_precision(ball_only.get(), 128) == CPS_OK && cps_options_set_max_work_bits(ball_only.get(), 256) == CPS_OK);
        auto cluster = certpolysolve_test::solve(Polynomial{1, -2, 1}, ball_only.get());
        CHECK(cluster.complete() && !cluster.goal_met());
        CHECK(cluster.roots()[0].count == 2 && cluster.roots()[0].multiplicity.value_or(0) == 0);
        std::cout << "C++ core tests passed (FLINT " << cps_backend_version() << ")\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n'; return 1;
    }
}
