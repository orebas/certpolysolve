#include "exact_winding_api.hpp"
#include <certpolysolve/certpolysolve.hpp>
#include <algorithm>
#include <iostream>
#include <random>
#include <tuple>

// Equivalent inputs must give equivalent answers.
//  A. An exact variable rescaling x -> 2^s x rescales every exported region
//     exactly: every floating-point decision is then exactly equivariant
//     (docs/secular-scaling-findings-20260928.md).
//  B. Planted roots with long dyadic ("irrational", "transcendental") and
//     complex values, rotated, shifted, conjugated or scaled: every exported
//     region holds exactly its count of the transformed planted roots.
//  C. Uncertain (ball) coefficients: the regions cover the roots of the
//     midpoint polynomial, a member of the input family.
// Checks use exact rational arithmetic only; no FLINT or solver output is trusted.
namespace {
using namespace reference;

// 64-bit dyadic stand-ins for pi and sqrt(2), used as exact input values:
// round(pi 2^64) and floor(sqrt(2) 2^64).
const Rational pi64 = Rational(Integer("57952155664616982739")) / power2(64);
const Rational sqrt2_64 = Rational(Integer("26087635650665564424")) / power2(64);
const Complex unit{Rational(3, 5), Rational(4, 5)}; // |unit| = 1, at an irrational angle
constexpr int bits = 53;

certpolysolve::Result solve(const certpolysolve::Polynomial &p, certpolysolve::Method method) {
    certpolysolve::Options options;
    options.method = method;
    options.goal = certpolysolve::Goal::approximation();
    options.bits = bits;
    return certpolysolve::solve(p, options);
}
const char *name(certpolysolve::Method method) {
    return method == certpolysolve::Method::classical ? "classical"
           : method == certpolysolve::Method::secular ? "secular" : "auto";
}

using Box = std::tuple<Rational, Rational, Rational, Rational, uint64_t>;
std::vector<Box> boxes(const certpolysolve::Result &result, const Rational &scale) {
    std::vector<Box> out;
    for (const auto &e : result.roots())
        out.emplace_back(endpoint(e.box.real.lo) * scale, endpoint(e.box.real.hi) * scale,
                         endpoint(e.box.imag.lo) * scale, endpoint(e.box.imag.hi) * scale, e.count);
    std::sort(out.begin(), out.end());
    return out;
}

// A: q(x) = p(2^s x) has roots r 2^-s; scaling q's regions by 2^s must give p's.
void exact_rescaling(const ComplexPolynomial &p, int s, certpolysolve::Method method, const char *label) {
    ComplexPolynomial q = p;
    Rational factor = 1;
    for (auto &c : q) {
        c = {c.real * factor, c.imag * factor};
        factor *= power2(s);
    }
    const auto original = solve(input(p), method), scaled = solve(input(q), method);
    require(original.complete() && original.goal_met() && scaled.complete() && scaled.goal_met(),
            "rescaling test solve did not meet its goal");
    if (boxes(original, 1) != boxes(scaled, power2(s))) {
        std::cerr << label << " " << name(method) << " s=" << s << '\n';
        require(false, "an exact 2^s rescaling changed the exported regions");
    }
}

bool contains(const certpolysolve::Root &region, const Complex &z) {
    return endpoint(region.box.real.lo) <= z.real && z.real <= endpoint(region.box.real.hi) &&
           endpoint(region.box.imag.lo) <= z.imag && z.imag <= endpoint(region.box.imag.hi);
}
// B, C: each planted root lies in exactly one region, and each region's count
// equals the planted roots inside it.
void verify_planted(const certpolysolve::Result &result, const std::vector<Complex> &roots, const char *label) {
    if (!(result.complete() && result.goal_met())) {
        std::cerr << label << ": " << result.message() << '\n';
        require(false, "planted-root solve did not meet its goal");
    }
    const auto &regions = result.roots();
    check_disjoint(regions);
    std::vector<uint64_t> found(regions.size());
    for (const auto &z : roots) {
        size_t inside = 0;
        for (size_t i = 0; i < regions.size(); ++i)
            if (contains(regions[i], z)) {
                ++found[i];
                ++inside;
            }
        if (inside != 1) std::cerr << label << '\n';
        require(inside == 1, "a planted root is not in exactly one exported region");
    }
    for (size_t i = 0; i < regions.size(); ++i) {
        require(found[i] == regions[i].count, "a region's count differs from its planted roots");
        // max(|re|,|im|) >= |z|/sqrt(2): allow that factor in the relative check.
        if (regions[i].count == 1) check_accuracy(regions[i], bits - 1);
    }
}

std::vector<Complex> planted() {
    return {{sqrt2_64}, {-pi64 / 4, Rational(1, 3)}, {Rational(3, 7), -sqrt2_64}, {pi64, pi64 / 8},
            {Rational(-5, 2)}, {Rational(1, 5), Rational(2)}, {Rational(-1, 3), Rational(-3, 4)},
            // A tight pair, 2^-40 apart: still resolvable at 53-bit relative accuracy.
            {Rational(7, 8), Rational(1, 8)}, {Rational(7, 8) + power2(-40), Rational(1, 8)},
            {Rational(0), Rational(-2)}};
}
Complex conj(const Complex &z) { return {z.real, -z.imag}; }
ComplexPolynomial scaled(ComplexPolynomial p, const Complex &c) {
    for (auto &x : p) x = x * c;
    return p;
}

// D: 96 rational points of the circle of radius 2^-16 about `centre`
// (t = tan(angle/2) = (2k+1)/32; k and -1-k give conjugate points).
std::vector<Complex> cloud(const Complex &centre) {
    std::vector<Complex> roots;
    for (int k = -48; k < 48; ++k) {
        const Rational t(2 * k + 1, 32), d = Rational(1) + t * t;
        roots.push_back(centre + Complex{(Rational(1) - t * t) / d * power2(-16), Rational(2) * t / d * power2(-16)});
    }
    return roots;
}
} // namespace

int main() {
    certpolysolve::ThreadScope lifetime;
    const certpolysolve::Method methods[] = {certpolysolve::Method::classical, certpolysolve::Method::secular,
                                         certpolysolve::Method::automatic};
    // A. Real and Gaussian random polynomials, and a near-multiple cluster.
    std::mt19937_64 random(0xe9a1e7);
    auto coefficient = [&] { return Rational(static_cast<long long>(random() % 2001) - 1000); };
    ComplexPolynomial real_random, gaussian_random;
    for (int k = 0; k <= 30; ++k) real_random.push_back({coefficient()});
    for (int k = 0; k <= 24; ++k) gaussian_random.push_back({coefficient(), coefficient()});
    std::vector<Complex> cluster_roots;
    // Four roots within 2^-20 of 3/4 + i/2, plus six spread roots.
    const Complex centre{Rational(3, 4), Rational(1, 2)};
    const Complex offsets[] = {{power2(-20)}, {0, power2(-20)}, {-power2(-20)}, {0, -power2(-20)}};
    for (const auto &o : offsets) cluster_roots.push_back(centre + o);
    for (int k = 1; k <= 6; ++k) cluster_roots.push_back({Rational(-k, 3), Rational(k % 3, 2)});
    const ComplexPolynomial cluster = from_roots(cluster_roots);
    for (const auto method : methods)
        for (const int s : {1, 3, -2}) {
            exact_rescaling(real_random, s, method, "real random degree 30");
            exact_rescaling(gaussian_random, s, method, "Gaussian random degree 24");
            exact_rescaling(cluster, s, method, "near-multiple cluster");
        }
    // B. Transformed planted roots.
    const auto roots = planted();
    struct Variant {
        const char *label;
        std::vector<Complex> roots;
        Complex scalar;
    };
    std::vector<Variant> variants;
    variants.push_back({"identity", roots, {1}});
    std::vector<Complex> rotated, shifted_roots, conjugated;
    for (const auto &r : roots) {
        rotated.push_back(r * conj(unit)); // roots of p(unit x)
        shifted_roots.push_back(r - Complex{sqrt2_64});
        conjugated.push_back(conj(r));
    }
    variants.push_back({"rotation by (3+4i)/5", rotated, {1}});
    variants.push_back({"shift by sqrt2", shifted_roots, {1}});
    variants.push_back({"conjugation", conjugated, {1}});
    variants.push_back({"scalar pi", roots, {pi64}});
    variants.push_back({"scalar (3+4i)/5", roots, unit});
    for (const auto &variant : variants)
        for (const auto method : methods) {
            const auto p = scaled(from_roots(variant.roots), variant.scalar);
            verify_planted(solve(input(p), method), variant.roots, variant.label);
        }
    // C. Ball coefficients (Auto: explicit methods require exact input).
    const std::vector<Complex> separated(roots.begin(), roots.begin() + 7);
    const auto midpoint = from_roots(separated);
    certpolysolve::Polynomial uncertain;
    for (size_t j = 0; j < midpoint.size(); ++j)
        uncertain.set(j, certpolysolve::ball(certpolysolve::Number::complex(certpolysolve::Number(midpoint[j].real.str()),
                                                                    certpolysolve::Number(midpoint[j].imag.str())),
                                         certpolysolve::Number("1@-200"),
                                         certpolysolve::Number(midpoint[j].imag == 0 ? "0" : "1@-200")));
    verify_planted(solve(uncertain, certpolysolve::Method::automatic), separated, "ball coefficients");
    // D. A tight cloud far from the origin, above centering_min_degree. The
    // explicit methods centre at the root centroid rounded to 8 bits. For the
    // 16-bit centres below that leaves an offset of about 2^-10, some 64 cloud
    // radii, which costs about 6 bits per root unless the centre is refined
    // (docs/equivalence-screen-findings-20260928.md, class B). The work-bit
    // cap makes the regression deterministic: unrefined, these need more.
    constexpr uint64_t cloud_work_bits = 512;
    const Complex centres[] = {{Rational(-3, 8)}, {Rational(14745, 65536)},
                               {Rational(-14745, 65536), Rational(19661, 65536)}};
    for (const auto &centre : centres) {
        const auto roots = cloud(centre);
        const auto p = from_roots(roots);
        for (const auto method : {certpolysolve::Method::classical, certpolysolve::Method::secular}) {
            certpolysolve::Options options;
            options.method = method;
            options.goal = certpolysolve::Goal::approximation();
            options.bits = bits;
            options.max_bits = cloud_work_bits;
            verify_planted(certpolysolve::solve(input(p), options), roots, "off-centre tight cloud");
            exact_rescaling(p, -3, method, "off-centre tight cloud");
        }
    }
    std::cout << "equivalence properties hold\n";
    return 0;
}
