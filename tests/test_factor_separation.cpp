// Covers of different exact square-free factors that overlap are tightened by doubling the affected factor's
// goal, and a factor whose working precision is already at the cap is tightened all the same: its own step reports
// the cap if it cannot meet the tighter goal there (solver.cpp, the overlap branch of the global driver; 18 October
// 2026). The input is the axis fuzz family's seed 5952222659831003457, degree 16: the triple real root -1/2 lies
// 2^-53 from the double pair -1/2 +- i 2^-53, and the double real root -251927/262144 lies 2^-64 from the simple
// pair -251927/262144 +- i 2^-64, each beside roots of another square-free factor. Under a 16-bit goal, Auto, forced
// secular and forced piecewise (which hands its stops to secular) stopped for the precision at every cap below 304
// bits, because the tighter goal was set from the factor's working precision and the solve stopped once that
// precision reached the cap; the old default route answered from 160 bits. Checked: each of the three answers at
// caps of 160 and 256 bits, every region one distinct root of its multiplicity and its count proved exactly on the
// input (rectangle_count), the regions disjoint and the counts summing to the degree; forced classical as well.
#include "exact_winding_api.hpp"
#include <certpolysolve/certpolysolve.hpp>
#include <iostream>
#include <string>

using namespace reference;

namespace {
ComplexPolynomial seed_input() {
    const Rational half = Rational(-1) / Rational(2), near = Rational(-251927) / Rational(262144);
    const Rational tiny53 = power2(-53), tiny64 = power2(-64);
    std::vector<Complex> roots{Complex(half, tiny53),  Complex(half, -tiny53), Complex(half, tiny53),
                               Complex(half, -tiny53), Complex(near, tiny64),  Complex(near, -tiny64),
                               Complex(half, 0),       Complex(half, 0),       Complex(half, 0),
                               Complex(near, 0),       Complex(near, 0),       Complex(Rational(5233529) / Rational(16777216), 0),
                               Complex(Rational(-5) / Rational(8), Rational(43) / Rational(512)),
                               Complex(Rational(-5) / Rational(8), Rational(-43) / Rational(512)),
                               Complex(Rational(-10033) / Rational(16384), Rational(5) / Rational(8)),
                               Complex(Rational(-10033) / Rational(16384), Rational(-5) / Rational(8))};
    return from_roots(roots);
}

void check(const ComplexPolynomial &p, certpolysolve::Method method, uint64_t cap) {
    const std::string what = "method " + std::to_string(static_cast<int>(method)) + ", cap " + std::to_string(cap);
    certpolysolve::Options options;
    options.method = method;
    options.goal = certpolysolve::Goal::accuracy();
    options.bits = 16;
    options.max_bits = cap;
    const auto result = certpolysolve::solve(input(p), options);
    require(cps_result_termination(result.get()) == CPS_SUCCESS && result.complete() && result.goal_met(),
            (what + ": " + result.message()).c_str());
    const auto regions = result.roots();
    check_disjoint(regions);
    uint64_t total = 0;
    for (const auto &region : regions) {
        require(region.multiplicity.value_or(0) == region.count, (what + ": a region is not one distinct root").c_str());
        const auto count = rectangle_count(p, region);
        require(count && static_cast<uint64_t>(*count) == region.count,
                (what + ": a region's count is not the input's").c_str());
        total += region.count;
    }
    require(total == 16, (what + ": the counts do not sum to the degree").c_str());
}
} // namespace

int main() {
    certpolysolve::ThreadScope lifetime;
    try {
        const auto p = seed_input();
        for (const certpolysolve::Method method : {certpolysolve::Method::automatic, certpolysolve::Method::secular,
                                                certpolysolve::Method::piecewise, certpolysolve::Method::classical})
            for (const uint64_t cap : {uint64_t{160}, uint64_t{256}}) check(p, method, cap);
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
    std::cout << "factor separation tests passed\n";
    return 0;
}
