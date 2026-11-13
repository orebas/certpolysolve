// At a stop short of the goal that is not a cancellation, the global driver (SolveState, solver.cpp) completes every
// factor that has had no step yet with its stop cover (complete_cover_at_stop on the factor's own polynomial in x), as
// piecewise's factor loop completes an unfinished factor (the S8 map's rule 3; 23 October 2026). Before, such a
// stop published the degree proof's Cauchy square of the whole input. Checked under forced classical and forced
// secular (both run the global driver) at work budgets of 2, 3 and 4 units (the degree proof, the preparation, then one
// factor step each): the answer is complete, its regions are disjoint, each region's count is the input's own count in
// it (rectangle_count, exact) and the counts sum to the degree; and unless the goal was met, every region lies in the
// square of half-side 1024 about the origin, far inside each input's Cauchy bound (5,001 to 24,000), so the stop did
// not keep the Cauchy square. The inputs: two square-free factors whose roots lie at very different distances from the
// origin, (x - 100)(x + 120)(x^2 + 1)^2; the square-free (x - 100)(x + 120)(x^2 + 1), whose one factor has had no step
// at a budget of 2; and a power composite beside a double root, (x^3 - 2)(x - 50)^2, whose factor x^3 - 2 is stored as
// w - 2 (its stop cover is computed on x^3 - 2 itself: a square about the origin has no lift from w - 2).
#include "exact_winding_api.hpp"
#include <certpolysolve/certpolysolve.hpp>
#include <iostream>
#include <string>

using namespace reference;

namespace {
// The square about the origin that every published region lies in, by its half-side.
Rational reach(const std::vector<certpolysolve::Root> &regions) {
    Rational largest = 0;
    for (const auto &region : regions)
        for (const auto &end : {region.box.real.lo, region.box.real.hi, region.box.imag.lo, region.box.imag.hi})
            largest = std::max(largest, absolute(endpoint(end)));
    return largest;
}

void check(const ComplexPolynomial &p, uint64_t degree, certpolysolve::Method method, uint64_t budget,
           const std::string &name) {
    const std::string what = name + ", method " + std::to_string(static_cast<int>(method)) + ", budget " +
                             std::to_string(budget);
    certpolysolve::Options options;
    options.method = method;
    options.work_budget = budget;
    const auto result = certpolysolve::solve(input(p), options);
    require(result.complete(), (what + ": not complete: " + result.message()).c_str());
    const auto regions = result.roots();
    check_disjoint(regions);
    uint64_t total = 0;
    for (const auto &region : regions) {
        const auto count = rectangle_count(p, region);
        require(count && static_cast<uint64_t>(*count) == region.count,
                (what + ": a region's count is not the input's").c_str());
        total += region.count;
    }
    require(total == degree, (what + ": the counts do not sum to the degree").c_str());
    if (!result.goal_met())
        require(reach(regions) <= Rational(1024),
                (what + ": the stop kept a region reaching the input's Cauchy bound: " + result.message()).c_str());
}
} // namespace

int main() {
    certpolysolve::ThreadScope lifetime;
    try {
        const Complex i(Rational(0), Rational(1)), minus_i(Rational(0), Rational(-1));
        const Complex far1(Rational(100), Rational(0)), far2(Rational(-120), Rational(0));
        const auto two_factors = from_roots({far1, far2, i, minus_i, i, minus_i});
        const auto square_free = from_roots({far1, far2, i, minus_i});
        // (x^3 - 2)(x - 50)^2 from its coefficients: x^5 - 100 x^4 + 2500 x^3 - 2 x^2 + 200 x - 5000.
        ComplexPolynomial power_composite;
        for (const int c : {-5000, 200, -2, 2500, -100, 1})
            power_composite.push_back(Complex(Rational(c), Rational(0)));
        for (const certpolysolve::Method method : {certpolysolve::Method::classical, certpolysolve::Method::secular})
            for (const uint64_t budget : {uint64_t{2}, uint64_t{3}, uint64_t{4}}) {
                check(two_factors, 6, method, budget, "(x - 100)(x + 120)(x^2 + 1)^2");
                check(square_free, 4, method, budget, "(x - 100)(x + 120)(x^2 + 1)");
                check(power_composite, 5, method, budget, "(x^3 - 2)(x - 50)^2");
            }
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
    std::cout << "stop cover tests passed\n";
    return 0;
}
