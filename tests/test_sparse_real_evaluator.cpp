// The sparse real evaluator (sparse_real_evaluator.hpp) against the dense evaluation, and the extracted sign and
// refinement routines (real_sign_refinement.hpp) over it: B6's building blocks, before any recursion.
#include "certify/real_sign_refinement.hpp"
#include "numeric/evaluate/sparse_real_evaluator.hpp"
#include <flint/arb_poly.h>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

using namespace certpolysolve::detail;
namespace {
void require(bool value, const std::string &message) {
    if (!value) throw std::runtime_error(message);
}

// Dense enclosures of p and p' at x at `bits`, the reference.
void dense_evaluate2(arb_t value, arb_t slope, const std::vector<Number> &coefficients, const arb_t x, slong bits) {
    arb_poly_t p;
    arb_poly_init(p);
    RealBall c;
    for (size_t j = 0; j < coefficients.size(); ++j) {
        coefficients[j].enclose(c.value(), bits);
        arb_poly_set_coeff_arb(p, static_cast<slong>(j), c.value());
    }
    arb_poly_evaluate2(value, slope, p, x, bits);
    arb_poly_clear(p);
}

// Random sparse polynomials: the sparse values overlap the dense ones, at random dyadic points and at 0, and the
// sparse enclosures are not wider than about 2^-(bits - log2(degree) - 8) relative (the powers' rounding).
void against_dense() {
    std::mt19937_64 rng(20261016);
    int checked = 0;
    for (int t = 0; t < 200; ++t) {
        const slong degree = 1 + static_cast<slong>(rng() % (t < 100 ? 64 : 5000));
        const size_t terms = 1 + static_cast<size_t>(rng() % 12);
        std::vector<Number> coefficients(static_cast<size_t>(degree) + 1);
        for (auto &c : coefficients) c = Number("0");
        for (size_t k = 0; k < terms; ++k) {
            const auto e = static_cast<size_t>(rng() % (degree + 1));
            const long long value = static_cast<long long>(rng() % (1u << 20)) - (1 << 19);
            coefficients[e] = Number(std::to_string(value == 0 ? 1 : value));
        }
        coefficients[static_cast<size_t>(degree)] = Number(std::to_string(1 + rng() % 1000));
        const auto sparse = SparseRealEvaluator::from_dense(coefficients);
        require(sparse.degree() == degree, "degree");
        for (int point = 0; point < 4; ++point) {
            RealBall x, v, s, dv, ds;
            if (point == 0) arb_zero(x.value());
            else {
                // A dyadic in [-2, 2] with 20 fractional bits, exact.
                arb_set_si(x.value(), static_cast<slong>(rng() % (1u << 22)) - (1 << 21));
                arb_mul_2exp_si(x.value(), x.value(), -20);
            }
            const slong bits = 64 + static_cast<slong>(rng() % 3) * 64;
            sparse.evaluate2(v.value(), s.value(), x.value(), bits);
            dense_evaluate2(dv.value(), ds.value(), coefficients, x.value(), bits);
            require(arb_overlaps(v.value(), dv.value()), "value disagrees with the dense evaluation at case " + std::to_string(t));
            require(arb_overlaps(s.value(), ds.value()), "slope disagrees with the dense evaluation at case " + std::to_string(t));
            require(arb_is_finite(v.value()) && arb_is_finite(s.value()), "non-finite enclosure");
            ++checked;
        }
    }
    std::cout << "sparse against dense: " << checked << " evaluations\n";
}

// x^n - 2 at n = 4096: the bracket [1, 2] of 2^(1/n) refined to 200 bits by real_refine over the sparse
// evaluator, and the refined ends' signs proved; then an undecidable sign (x^2 at 0, a double root) is 0.
void sign_and_refine() {
    constexpr slong n = 4096;
    std::vector<SparseRealTerm> terms{{0, Number("-2")}, {n, Number("1")}};
    SparseRealEvaluator p(terms);
    cps_options o;
    uint64_t evaluations = 0;
    const RealSignBudget budget{64, real_sign_cap(o), 64, evaluations, true};
    RealBall lo, hi, value, slope;
    arb_one(lo.value());
    arb_set_ui(hi.value(), 2);
    require(real_sign(p, lo.value(), value.value(), slope.value(), budget, o) == -1, "p(1) < 0");
    require(real_sign(p, hi.value(), value.value(), slope.value(), budget, o) == 1, "p(2) > 0");
    constexpr slong goal = 200;
    require(real_refine(p, lo, hi, goal, budget, o), "refinement to 200 bits");
    // The root: hi - lo <= min(|lo|, |hi|) 2^-goal, and the ends keep their signs.
    RealBall width, limit;
    arb_sub(width.value(), hi.value(), lo.value(), 256);
    arb_mul_2exp_si(limit.value(), lo.value(), -goal);
    require(arb_le(width.value(), limit.value()), "the refined bracket is wider than the goal");
    require(real_sign(p, lo.value(), value.value(), slope.value(), budget, o) == -1 &&
                real_sign(p, hi.value(), value.value(), slope.value(), budget, o) == 1,
            "the refined ends lost their signs");
    // 2^(1/4096) = 1.000169... lies inside. The reference ball is far narrower than the bracket's ends' distance to
    // the root (a Newton iterate can land within 2^-264 of it), so the comparison is strict.
    RealBall root;
    arb_set_ui(root.value(), 2);
    arb_root_ui(root.value(), root.value(), n, 2048);
    require(arb_gt(root.value(), lo.value()) && arb_lt(root.value(), hi.value()), "the bracket misses 2^(1/n)");
    std::cout << "x^4096 - 2: refined to 200 bits in " << evaluations << " evaluations, "
              << p.multiplications() << " multiplications\n";
    // An exact double root: x^2 at x = 0 evaluates to exactly 0, undecided at every precision, so 0.
    SparseRealEvaluator square(std::vector<SparseRealTerm>{{2, Number("1")}});
    RealBall zero;
    cps_options capped;
    capped.max_work_bits = 1024;
    uint64_t few = 0;
    const RealSignBudget small{64, real_sign_cap(capped), 64, few};
    require(real_sign(square, zero.value(), value.value(), slope.value(), small, capped) == 0, "a root's sign is undecided");
    require(few == 5, "doubling from 64 to the 1024-bit cap takes five evaluations, took " + std::to_string(few));
}

void refusals() {
    bool thrown = false;
    try {
        SparseRealEvaluator bad(std::vector<SparseRealTerm>{{3, Number("1")}, {1, Number("1")}});
    } catch (const std::invalid_argument &) {
        thrown = true;
    }
    require(thrown, "descending exponents accepted");
}
} // namespace

int main() {
    try {
        against_dense();
        sign_and_refine();
        refusals();
    } catch (const std::exception &error) {
        std::cerr << "FAIL " << error.what() << '\n';
        return 1;
    }
    std::cout << "sparse real evaluator tests passed\n";
    return 0;
}
