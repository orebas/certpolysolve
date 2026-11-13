// Auto's direct route for exact degree 1 and 2 (small_degree.cpp): every box holds its root (computed apart by the
// textbook formula at 2,048 bits), the counts and multiplicities are right (an exact double root is one region of multiplicity 2), the goal is
// met at 53 and 200 bits, and the inputs it must decline return nullopt (inexact coefficients, a nonzero part of
// degree 0, degree 5 (degrees 3 and 4: test_small_degree_cubic and _quartic), a real-only query, a region, a value the exact storage policy will not expand, the cluster-radius
// and isolation goals). Through the public API the solve names the route.
#include "methods/direct/small_degree.hpp"
#include "core/cover.hpp"
#include "core/raii.hpp"
#include <certpolysolve/certpolysolve.hpp>
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

cps_options options(uint64_t bits, cps_goal goal = CPS_ROOT_ACCURACY) {
    cps_options o;
    o.bits = bits;
    o.goal = goal;
    return o;
}

// The roots, as 2,048-bit balls, of c0 + c1 x + c2 x^2 (c2 may be zero): the reference the boxes must hold. The
// textbook formula loses bits to cancellation where the roots' moduli differ greatly (400 for roots 2^200 apart), so
// it runs far above the 200-bit goal tested.
std::vector<d::ComplexBall> reference(const d::DensePolynomial &p) {
    const slong bits = 2048;
    std::vector<d::ComplexBall> c(3), out;
    for (size_t k = 0; k < p.coefficients.size() && k < 3; ++k) p.coefficients[k].enclose(c[k].value(), bits);
    if (p.coefficients.size() == 2 || acb_is_zero(c[2].value())) {
        d::ComplexBall r;
        acb_div(r.value(), c[0].value(), c[1].value(), bits);
        acb_neg(r.value(), r.value());
        out.push_back(std::move(r));
        return out;
    }
    d::ComplexBall disc, t, root;
    acb_mul(disc.value(), c[1].value(), c[1].value(), bits);
    acb_mul(t.value(), c[2].value(), c[0].value(), bits);
    acb_mul_2exp_si(t.value(), t.value(), 2);
    acb_sub(disc.value(), disc.value(), t.value(), bits);
    if (acb_contains_zero(disc.value())) {
        // A double root (the discriminant's ball holds 0 here only for the exact double roots below): -c1 / (2 c2),
        // without the square root of a ball about 0, whose imaginary radius no real box could hold.
        acb_div(root.value(), c[1].value(), c[2].value(), bits);
        acb_mul_2exp_si(root.value(), root.value(), -1);
        acb_neg(root.value(), root.value());
        out.push_back(root);
        out.push_back(root);
        return out;
    }
    acb_sqrt(disc.value(), disc.value(), bits);
    for (int sign : {1, -1}) {
        acb_mul_si(t.value(), disc.value(), sign, bits);
        acb_sub(root.value(), t.value(), c[1].value(), bits);
        acb_div(root.value(), root.value(), c[2].value(), bits);
        acb_mul_2exp_si(root.value(), root.value(), -1);
        out.push_back(root);
    }
    return out;
}

// Each reference root lies in exactly one box (the double root's two copies in the same one), the counts add up to
// the degree, and every box meets the goal.
void check(const char *name, const d::DensePolynomial &p, uint64_t bits, uint64_t regions, uint64_t multiplicity_two,
           cps_squarefree squarefree) {
    const auto o = options(bits);
    const auto result = d::direct_small_degree(p, o);
    require(result.has_value(), std::string(name) + ": the direct route declined");
    require(result->goal_met && result->complete && result->termination == CPS_SUCCESS, std::string(name) + ": goal");
    require(result->roots.size() == regions, std::string(name) + ": region count");
    require(result->squarefree == squarefree, std::string(name) + ": square-free verdict");
    uint64_t total = 0, doubles = 0;
    for (const auto &root : result->roots) {
        total += root.count;
        doubles += root.multiplicity == 2;
        require(root.multiplicity == root.count, std::string(name) + ": a multiplicity is unknown");
        require(d::meets_goal(root, o, static_cast<slong>(result->work_bits)), std::string(name) + ": a box misses the goal");
    }
    require(total + 0 == static_cast<uint64_t>(result->degree), std::string(name) + ": the counts do not add up");
    require(doubles == multiplicity_two, std::string(name) + ": double roots");
    // The nonzero part's reference roots (the exact zero roots are the exact point 0).
    size_t zeros = 0;
    while (zeros + 1 < p.coefficients.size() && p.coefficients[zeros].zero()) ++zeros;
    d::DensePolynomial nonzero;
    nonzero.coefficients.assign(p.coefficients.begin() + static_cast<ptrdiff_t>(zeros), p.coefficients.end());
    for (const auto &r : reference(nonzero)) {
        int holding = 0;
        for (const auto &root : result->roots) holding += acb_contains(root.box.value(), r.value());
        require(holding == 1, std::string(name) + ": a root is not in exactly one box");
    }
}
} // namespace

int main() {
    try {
        for (uint64_t bits : {uint64_t(53), uint64_t(200)}) {
            check("x^2 - 2", polynomial({{"-2", "0"}, {"0", "0"}, {"1", "0"}}), bits, 2, 0, CPS_SQUAREFREE_PROVED);
            check("x^2 + 1", polynomial({{"1", "0"}, {"0", "0"}, {"1", "0"}}), bits, 2, 0, CPS_SQUAREFREE_PROVED);
            check("(x - 1/3)^2", polynomial({{"1/9", "0"}, {"-2/3", "0"}, {"1", "0"}}), bits, 1, 1,
                  CPS_SQUAREFREE_DISPROVED);
            check("3x + 1", polynomial({{"1", "0"}, {"3", "0"}}), bits, 1, 0, CPS_SQUAREFREE_PROVED);
            check("x^2 + 3x", polynomial({{"0", "0"}, {"3", "0"}, {"1", "0"}}), bits, 2, 0, CPS_SQUAREFREE_PROVED);
            check("x^2 + i", polynomial({{"0", "1"}, {"0", "0"}, {"1", "0"}}), bits, 2, 0, CPS_SQUAREFREE_PROVED);
            check("(x - (1+i))^2", polynomial({{"0", "2"}, {"-2", "-2"}, {"1", "0"}}), bits, 1, 1,
                  CPS_SQUAREFREE_DISPROVED);
            check("roots 2^200 apart", polynomial({{"1", "0"}, {"1@200", "0"}, {"1", "0"}}), bits, 2, 0,
                  CPS_SQUAREFREE_PROVED);
            check("1@600 x^2 - 1", polynomial({{"-1", "0"}, {"0", "0"}, {"1@600", "0"}}), bits, 2, 0,
                  CPS_SQUAREFREE_PROVED);
            // Within the compact-exponent allowance (number.hpp: one precision cap of growth for each nonzero number
            // written out): 2^-100000 writes out to about 100,000 bits against two caps of 65,536.
            check("x^2 - 2^-100000", polynomial({{"-1@-100000", "0"}, {"0", "0"}, {"1", "0"}}), bits, 2, 0,
                  CPS_SQUAREFREE_PROVED);
        }
        // Declined: the general route answers these.
        auto inexact = polynomial({{"-2", "0"}, {"0", "0"}, {"1", "0"}});
        inexact.coefficients[0].real_radius = d::Number("1@-60");
        require(!d::direct_small_degree(inexact, options(53)), "an inexact coefficient was taken");
        require(!d::direct_small_degree(polynomial({{"0", "0"}, {"0", "0"}, {"1", "0"}}), options(53)),
                "a nonzero part of degree 0 was taken");
        require(!d::direct_small_degree(polynomial({{"1", "0"}, {"0", "0"}, {"0", "0"}, {"0", "0"}, {"0", "0"}, {"1", "0"}}),
                                        options(53)),
                "degree 5 was taken");
        auto real_only = options(53);
        real_only.real_only = true;
        require(!d::direct_small_degree(polynomial({{"-2", "0"}, {"0", "0"}, {"1", "0"}}), real_only),
                "a real-only query was taken");
        // Beyond the allowance: an exponent of three caps, where the two nonzero numbers allow two.
        const std::string beyond = "-1@-" + std::to_string(3 * uint64_t(CPS_DEFAULT_MAX_WORK_BITS));
        require(!d::direct_small_degree(polynomial({{beyond.c_str(), "0"}, {"0", "0"}, {"1", "0"}}), options(53)),
                "a value the compact-exponent rule refuses was expanded");
        // The goals that do not ask for one region a distinct root: the general route answers them.
        auto cluster = options(53, CPS_CLUSTER_RADIUS);
        cluster.cluster_radius = d::Number("100");
        require(!d::direct_small_degree(polynomial({{"-2", "0"}, {"0", "0"}, {"1", "0"}}), cluster),
                "a cluster-radius goal was taken");
        require(!d::direct_small_degree(polynomial({{"-2", "0"}, {"0", "0"}, {"1", "0"}}), options(53, CPS_ROOT_ISOLATION)),
                "an isolation goal was taken");
        // Through the public API: the route is named, and an inexact input takes the general route.
        certpolysolve::Polynomial exact;
        exact.set(0, certpolysolve::Number("-2")).set(2, certpolysolve::Number("1"));
        const auto solved = certpolysolve::solve(exact);
        require(solved.goal_met() && solved.message().find("auto-route=direct") != std::string::npos,
                "the public solve did not take the direct route");
        certpolysolve::Polynomial ball;
        ball.set(0, certpolysolve::ball(certpolysolve::Number("-2"), certpolysolve::Number("1@-60"))).set(2, certpolysolve::Number("1"));
        const auto general = certpolysolve::solve(ball);
        require(general.message().find("auto-route=direct") == std::string::npos, "an inexact input took the direct route");
        std::cout << "direct route for degree 1 and 2: 18 inputs proved (two goals), 7 declined, the public route named\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << "test_small_degree: " << e.what() << '\n';
        return 1;
    }
}
