// Auto's direct route for exact cubics (small_degree.cpp, small_degree_forms.cpp): every root lies in exactly one box,
// a multiple root's box holds all its copies with the exact multiplicity, every box meets the goal, the realness of
// real input is exact (real roots in exactly real boxes, the pairs proved not real), and the square-free verdict is
// exact. Planted inputs: distinct real roots (Viete's form), one real root and a pair (the real Cardano form), exact
// double and triple roots, a negative leading coefficient, near-double roots 2^-40 apart and near-triple roots
// 2^-24 apart, nearly real pairs, roots from 2^-30 to 2^30 and from 2^-300 to 2^300, complex coefficients with
// simple, double and triple roots. Unplanted inputs (x^3 - 2, the casus irreducibilis x^3 - 3x + 1, x^3 - i, dense
// complex coefficients) against FLINT's roots of the exact input at 4,096 bits. At 53, 200 and 1,000 bits. Through the
// public API: the solve names the route, and a real-only query returns the real roots with the direct certificate.
#include "methods/direct/small_degree.hpp"
#include "core/cover.hpp"
#include "core/raii.hpp"
#include <flint/acb_poly.h>
#include <flint/fmpq.h>
#include <certpolysolve/certpolysolve.hpp>
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace d = certpolysolve::detail;

namespace {
void require(bool ok, const std::string &why) {
    if (!ok) throw std::runtime_error(why);
}
constexpr slong reference_bits = 4096;

struct Exact {  // a Gaussian rational
    d::Rational re, im;
};
Exact exact(const char *re, const char *im = "0") {
    Exact z;
    const d::Number a(re), b(im);
    uint64_t budget = UINT64_MAX / 4;
    require(a.rational_with_budget(z.re.value(), budget, true) && b.rational_with_budget(z.im.value(), budget, true),
            std::string("not a number: ") + re);
    return z;
}
std::string text(const fmpq_t q) {
    char *s = fmpq_get_str(nullptr, 10, q);
    std::string out(s);
    flint_free(s);
    return out;
}
struct Planted {
    Exact value;
    uint64_t multiplicity;
};
// lead times the product of (x - r) over the planted roots, each as often as its multiplicity, lowest degree first.
d::DensePolynomial expand(const Exact &lead, const std::vector<Planted> &roots) {
    std::vector<Exact> c(1);
    fmpq_set(c[0].re.value(), lead.re.value());
    fmpq_set(c[0].im.value(), lead.im.value());
    d::Rational t, u;
    for (const auto &root : roots)
        for (uint64_t m = 0; m < root.multiplicity; ++m) {
            std::vector<Exact> next(c.size() + 1);
            for (size_t k = 0; k < c.size(); ++k) {
                // next[k+1] += c[k]; next[k] -= r c[k]
                fmpq_add(next[k + 1].re.value(), next[k + 1].re.value(), c[k].re.value());
                fmpq_add(next[k + 1].im.value(), next[k + 1].im.value(), c[k].im.value());
                fmpq_mul(t.value(), root.value.re.value(), c[k].re.value());
                fmpq_mul(u.value(), root.value.im.value(), c[k].im.value());
                fmpq_sub(t.value(), t.value(), u.value());
                fmpq_sub(next[k].re.value(), next[k].re.value(), t.value());
                fmpq_mul(t.value(), root.value.re.value(), c[k].im.value());
                fmpq_mul(u.value(), root.value.im.value(), c[k].re.value());
                fmpq_add(t.value(), t.value(), u.value());
                fmpq_sub(next[k].im.value(), next[k].im.value(), t.value());
            }
            c = std::move(next);
        }
    d::DensePolynomial p;
    for (const auto &x : c) {
        d::Coefficient k;
        k.real = d::Number(text(x.re.value()));
        k.imag = d::Number(text(x.im.value()));
        p.coefficients.push_back(std::move(k));
    }
    return p;
}
d::DensePolynomial polynomial(std::initializer_list<std::pair<const char *, const char *>> coefficients) {
    d::DensePolynomial p;
    for (const auto &[re, im] : coefficients) {
        d::Coefficient c;
        c.real = d::Number(re);
        c.imag = d::Number(im);
        p.coefficients.push_back(std::move(c));
    }
    return p;
}
cps_options options(uint64_t bits) {
    cps_options o;
    o.bits = bits;
    o.goal = CPS_ROOT_ACCURACY;
    return o;
}
bool real_input(const d::DensePolynomial &p) {
    for (const auto &c : p.coefficients)
        if (!c.imag.zero()) return false;
    return true;
}
// Whether a box holds a reference root, given as a 4,096-bit ball, part by part: the ball's part inside the box's, or,
// where the box's part is one exact point (an exactly real or exactly imaginary root of the closed forms, which no ball
// of positive radius fits in), that point inside the ball's part.
bool holds(const acb_t box, const acb_t reference) {
    for (const bool imaginary : {false, true}) {
        const arb_struct *b = imaginary ? acb_imagref(box) : acb_realref(box);
        const arb_struct *r = imaginary ? acb_imagref(reference) : acb_realref(reference);
        if (!arb_contains(b, r) && !(arb_is_exact(b) && arb_contains(r, b))) return false;
    }
    return true;
}
std::string ball_text(const acb_t z) {
    char *re = arb_get_str(acb_realref(z), 30, 0), *im = arb_get_str(acb_imagref(z), 30, 0);
    std::string out = std::string(re) + " + i " + im;
    flint_free(re);
    flint_free(im);
    return out;
}
struct Expectation {
    uint64_t regions;
    std::vector<uint64_t> multiplicities;  // sorted descending
    cps_squarefree squarefree;
};
uint64_t checked = 0;
// The direct route's answer: the counts and multiplicities, the goal, realness for real input, and each reference
// root (an exact value, or a reference ball) in exactly one box whose multiplicity is the root's.
void check(const std::string &name, const d::DensePolynomial &p, const std::vector<std::pair<d::ComplexBall, uint64_t>> &reference,
           const Expectation &want) {
    for (const uint64_t bits : {uint64_t(53), uint64_t(200), uint64_t(1000)}) {
        const auto o = options(bits);
        const auto result = d::direct_small_degree(p, o);
        const std::string at = name + " at " + std::to_string(bits) + " bits";
        require(result.has_value(), at + ": the direct route declined");
        require(result->goal_met && result->complete && result->termination == CPS_SUCCESS, at + ": goal");
        require(result->roots.size() == want.regions, at + ": region count " + std::to_string(result->roots.size()));
        require(result->squarefree == want.squarefree, at + ": square-free verdict");
        std::vector<uint64_t> multiplicities;
        uint64_t total = 0;
        for (const auto &root : result->roots) {
            require(root.multiplicity == root.count, at + ": a multiplicity is unknown");
            require(d::meets_goal(root, o, static_cast<slong>(result->work_bits)), at + ": a box misses the goal");
            require(root.certificate.find("closed form") != std::string::npos, at + ": the certificate");
            if (real_input(p))
                require(root.realness == CPS_REAL ? arb_is_zero(acb_imagref(root.box.value()))
                                                 : root.realness == CPS_NONREAL,
                        at + ": realness not decided for real input");
            total += root.count;
            multiplicities.push_back(root.multiplicity);
        }
        std::sort(multiplicities.rbegin(), multiplicities.rend());
        require(multiplicities == want.multiplicities, at + ": multiplicities");
        require(total == static_cast<uint64_t>(result->degree), at + ": the counts do not add up");
        for (const auto &[r, multiplicity] : reference) {
            int holding = 0;
            for (const auto &root : result->roots)
                if (holds(root.box.value(), r.value())) {
                    ++holding;
                    require(root.multiplicity == multiplicity, at + ": a root's box has another multiplicity");
                    if (real_input(p) && arb_is_zero(acb_imagref(r.value())))
                        require(root.realness == CPS_REAL, at + ": a real root's box is not proved real");
                    if (real_input(p) && !arb_contains_zero(acb_imagref(r.value())))
                        require(root.realness == CPS_NONREAL, at + ": a nonreal root's box is not proved nonreal");
                }
            if (holding != 1) {
                std::string boxes;
                for (const auto &root : result->roots) boxes += "\n  box " + ball_text(root.box.value());
                require(false, at + ": a root is not in exactly one box (" + std::to_string(holding) + "): root " +
                                   ball_text(r.value()) + boxes);
            }
        }
        ++checked;
    }
}
std::vector<std::pair<d::ComplexBall, uint64_t>> exact_reference(const std::vector<Planted> &roots) {
    std::vector<std::pair<d::ComplexBall, uint64_t>> out;
    for (const auto &r : roots) {
        d::ComplexBall b;
        arb_set_fmpq(acb_realref(b.value()), r.value.re.value(), reference_bits);
        arb_set_fmpq(acb_imagref(b.value()), r.value.im.value(), reference_bits);
        out.emplace_back(std::move(b), r.multiplicity);
    }
    return out;
}
std::vector<std::pair<d::ComplexBall, uint64_t>> numerical_reference(const d::DensePolynomial &p) {
    d::ComplexPolynomial f;
    for (size_t k = 0; k < p.coefficients.size(); ++k) {
        d::ComplexBall c;
        p.coefficients[k].enclose(c.value(), reference_bits);
        acb_poly_set_coeff_acb(f.value(), static_cast<slong>(k), c.value());
    }
    const slong n = acb_poly_degree(f.value());
    acb_ptr roots = _acb_vec_init(n);
    require(acb_poly_find_roots(roots, f.value(), nullptr, 0, reference_bits) == n, "the reference roots");
    std::vector<std::pair<d::ComplexBall, uint64_t>> out;
    for (slong k = 0; k < n; ++k) {
        d::ComplexBall b;
        acb_set(b.value(), roots + k);
        // Real input: a reference root whose ball meets the real axis is real (the input's nonreal roots come in
        // conjugate pairs well apart from it here), so its imaginary part is set to exactly 0 for the realness check.
        if (real_input(p) && arb_contains_zero(acb_imagref(b.value()))) arb_zero(acb_imagref(b.value()));
        out.emplace_back(std::move(b), 1);
    }
    _acb_vec_clear(roots, n);
    return out;
}
void planted(const std::string &name, const Exact &lead, const std::vector<Planted> &roots, const Expectation &want) {
    check(name, expand(lead, roots), exact_reference(roots), want);
}
} // namespace

int main() {
    try {
        const cps_squarefree proved = CPS_SQUAREFREE_PROVED, disproved = CPS_SQUAREFREE_DISPROVED;
        const Expectation three{3, {1, 1, 1}, proved};
        // Real coefficients, three real roots (Viete's form), one real root and a pair (the real Cardano form).
        planted("(x-1)(x-2)(x-3)", exact("1"), {{exact("1"), 1}, {exact("2"), 1}, {exact("3"), 1}}, three);
        planted("2(x-1/3)(x+5/7)(x-11)", exact("2"), {{exact("1/3"), 1}, {exact("-5/7"), 1}, {exact("11"), 1}}, three);
        planted("-3(x-1)(x-2)(x+4)", exact("-3"), {{exact("1"), 1}, {exact("2"), 1}, {exact("-4"), 1}}, three);
        planted("(x+1)(x^2+1)", exact("1"), {{exact("-1"), 1}, {exact("0", "1"), 1}, {exact("0", "-1"), 1}}, three);
        planted("-5(x-2)(x-(1+3i))(x-(1-3i))", exact("-5"),
                {{exact("2"), 1}, {exact("1", "3"), 1}, {exact("1", "-3"), 1}}, three);
        // Exact multiple roots.
        planted("(x-1)^2(x+2)", exact("1"), {{exact("1"), 2}, {exact("-2"), 1}}, {2, {2, 1}, disproved});
        planted("-7(x+3/2)^2(x-5)", exact("-7"), {{exact("-3/2"), 2}, {exact("5"), 1}}, {2, {2, 1}, disproved});
        planted("(x-1/3)^3", exact("1"), {{exact("1/3"), 3}}, {1, {3}, disproved});
        // Near-multiple roots, wide and extreme scales.
        planted("near-double 2^-40", exact("1"),
                {{exact("1"), 1}, {exact("1099511627777@-40"), 1}, {exact("-2"), 1}}, three);
        planted("near-triple 2^-24", exact("1"),
                {{exact("1"), 1}, {exact("16777217@-24"), 1}, {exact("16777215@-24"), 1}}, three);
        planted("nearly real pair 2^-40", exact("1"),
                {{exact("2"), 1}, {exact("1", "1@-40"), 1}, {exact("1", "-1@-40"), 1}}, three);
        planted("roots 2^-30 to 2^30", exact("1"), {{exact("1@30"), 1}, {exact("1"), 1}, {exact("1@-30"), 1}}, three);
        planted("roots 2^-300 to 2^300", exact("3"),
                {{exact("1@300"), 1}, {exact("-1@-300"), 1}, {exact("3"), 1}}, three);
        // Complex coefficients.
        planted("(x-i)(x-2)(x+1+i)", exact("1"), {{exact("0", "1"), 1}, {exact("2"), 1}, {exact("-1", "-1"), 1}}, three);
        planted("(2+i)(x-(1+i))^2(x-2)", exact("2", "1"), {{exact("1", "1"), 2}, {exact("2"), 1}}, {2, {2, 1}, disproved});
        planted("(x-i)^3", exact("1"), {{exact("0", "1"), 3}}, {1, {3}, disproved});
        planted("complex near-double 2^-30", exact("1"),
                {{exact("0", "1"), 1}, {exact("1@-30", "1"), 1}, {exact("-3"), 1}}, three);
        // Unplanted: D0 = 0 (a pure cube after the shift), the casus irreducibilis, complex D0 = 0, dense complex.
        for (const auto &[name, p] : std::vector<std::pair<std::string, d::DensePolynomial>>{
                 {"x^3 - 2", polynomial({{"-2", "0"}, {"0", "0"}, {"0", "0"}, {"1", "0"}})},
                 {"x^3 - 3x + 1", polynomial({{"1", "0"}, {"-3", "0"}, {"0", "0"}, {"1", "0"}})},
                 {"x^3 - i", polynomial({{"0", "-1"}, {"0", "0"}, {"0", "0"}, {"1", "0"}})},
                 {"dense complex", polynomial({{"0.3", "0.7"}, {"-1.1", "0.25"}, {"0.6", "-0.9"}, {"-0.45", "0.15"}})},
                 {"dense real", polynomial({{"-0.13", "0"}, {"-0.86", "0"}, {"-0.82", "0"}, {"-0.15", "0"}})}})
            check(name, p, numerical_reference(p), three);
        // Through the public API: the route is named; a real-only query returns the three real roots, proved real,
        // with the direct certificate.
        certpolysolve::Polynomial cubic;
        cubic.set(0, certpolysolve::Number("1")).set(1, certpolysolve::Number("-3")).set(3, certpolysolve::Number("1"));
        const auto solved = certpolysolve::solve(cubic);
        require(solved.goal_met() && solved.message().find("auto-route=direct") != std::string::npos,
                "the public solve did not take the direct route");
        const auto real = certpolysolve::real_roots(cubic);
        require(real.goal_met() && real.roots().size() == 3, "the real-only query: three real roots");
        for (const auto &root : real.roots())
            require(root.realness == certpolysolve::Realness::real && root.certificate.find("closed form") != std::string::npos,
                    "the real-only query: a root not proved real by the direct route");
        certpolysolve::Polynomial one_real;
        one_real.set(0, certpolysolve::Number("-2")).set(3, certpolysolve::Number("1"));
        const auto single = certpolysolve::real_roots(one_real);
        require(single.goal_met() && single.roots().size() == 1 &&
                    single.roots()[0].realness == certpolysolve::Realness::real,
                "the real-only query of x^3 - 2: one real root");
        std::cout << "direct route for cubics: " << checked << " solves (22 inputs at three goals) checked against "
                  << "their exact or reference roots; the public solve and two real-only queries\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << "test_small_degree_cubic: " << e.what() << '\n';
        return 1;
    }
}
