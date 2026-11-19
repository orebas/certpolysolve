// The certificate attempt's conjugate pairs (src/certify/attempts/newton_cover.cpp), through the public API.
//
// On a polynomial with real coefficients the attempt certifies one point of each conjugate pair of points and gives
// the other the reflection of the region it proved (tests/test_conjugate_pairs.cpp checks the pairing and the attempt
// from inside). What must hold in an answer, whatever the pairing guessed:
// - every region holds exactly its count of planted roots, and the counts sum to the degree;
// - a reflected region lies strictly on one side of the real axis, is the exact mirror image of a region that was
//   proved (the real ends equal, the imaginary ends exchanged and negated), and says so in its certificate;
// - a polynomial with complex coefficients gets no reflected region;
// - a conjugate pair so close to the axis that its rectangle would meet it, two pairs closer together than the
//   pairing's tolerance, and real roots beside pairs are all solved, each region with exactly its root.
// And what the pairing must find: every pair of the fixtures whose roots are well apart, among them roots with round
// coordinates (i, 1 + i, a lattice), which an earlier draft's grid missed, and lines of roots that share a coordinate.
#include "exact_winding_api.hpp"
#include <certpolysolve/certpolysolve.hpp>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

using namespace reference;

namespace {
// reference::require takes a literal; the messages here carry the fixture's name.
void require(bool condition, const std::string &message) { reference::require(condition, message.c_str()); }

const std::string reflected_mark = "reflected in the real axis";

struct Outcome {
    uint64_t regions = 0, reflected = 0;
};
// Solves, checks every region against the planted roots, and checks each reflected region against its mirror image.
Outcome solve_and_check(const std::vector<Complex> &planted, certpolysolve::Method method, const std::string &name) {
    const auto input = reference::input(from_roots(planted));
    certpolysolve::Options options;
    options.method = method;
    options.goal = certpolysolve::Goal::approximation();
    options.bits = 53;
    options.max_bits = 4096;
    const auto result = certpolysolve::solve(input, options);
    require(result.goal_met() && result.stop() == certpolysolve::Stop::none && result.complete(), name + ": not solved");
    Outcome out;
    uint64_t total = 0;
    const auto roots = result.roots();
    for (const auto &region : roots) {
        const Rational real_lower = endpoint(region.box.real.lo), real_upper = endpoint(region.box.real.hi);
        const Rational imag_lower = endpoint(region.box.imag.lo), imag_upper = endpoint(region.box.imag.hi);
        uint64_t inside = 0;
        for (const auto &root : planted)
            inside += real_lower <= root.real && root.real <= real_upper && imag_lower <= root.imag &&
                      root.imag <= imag_upper;
        require(inside == region.count && region.count == 1, name + ": a region does not hold exactly its planted root");
        total += region.count;
        ++out.regions;
        if (region.certificate.find(reflected_mark) == std::string::npos)
            continue;
        ++out.reflected;
        require(imag_upper < Rational(0) || Rational(0) < imag_lower, name + ": a reflected region meets the real axis");
        // Its mirror image is a proved region on the other side, end for end.
        uint64_t mirrors = 0;
        for (const auto &other : roots)
            mirrors += endpoint(other.box.real.lo) == real_lower && endpoint(other.box.real.hi) == real_upper &&
                       endpoint(other.box.imag.lo) == -imag_upper && endpoint(other.box.imag.hi) == -imag_lower &&
                       other.certificate.find(reflected_mark) == std::string::npos;
        require(mirrors == 1, name + ": a reflected region is not the mirror image of exactly one proved region");
    }
    require(total == planted.size(), name + ": the cover lost roots");
    return out;
}
// The roots x + iy and x - iy for each (x, y) given, and the real roots given.
std::vector<Complex> pairs_and_real(const std::vector<std::pair<Rational, Rational>> &upper,
                                    const std::vector<Rational> &real) {
    std::vector<Complex> planted;
    for (const auto &[x, y] : upper) {
        planted.push_back(Complex(x, y));
        planted.push_back(Complex(x, -y));
    }
    for (const Rational &x : real)
        planted.push_back(Complex(x, Rational(0)));
    return planted;
}
} // namespace

int main() {
    certpolysolve::ThreadScope lifetime;
    try {
        uint64_t reflected = 0;
        for (const auto method : {certpolysolve::Method::automatic, certpolysolve::Method::secular, certpolysolve::Method::classical}) {
            const std::string by = ", method " + std::to_string(static_cast<int>(method));
            // The methods whose attempt sees the whole polynomial find every pair; the classical method's count is
            // reported and bounded.
            const auto require_pairs = [&](const Outcome &outcome, uint64_t pairs, const std::string &name) {
                if (method == certpolysolve::Method::classical)
                    require(outcome.reflected <= pairs, name + by + ": more regions reflected than pairs");
                else
                    require(outcome.reflected == pairs, name + by + ": " + std::to_string(outcome.reflected) +
                                                            " regions were reflected for " + std::to_string(pairs) +
                                                            " pairs");
                reflected += outcome.reflected;
            };
            // 1. Twelve conjugate pairs at rational points and three real roots. The pairs k = 2, 4 and 9 have a
            //    coordinate that is a power of two (x = -2, y = 1, x = 1).
            {
                std::vector<std::pair<Rational, Rational>> upper;
                for (int k = 1; k <= 12; ++k)
                    upper.push_back({Rational(3 * k - 20) / 7, Rational(2 * k + 1) / 9});
                const auto planted = pairs_and_real(upper, {Rational(-3) / 2, Rational(1) / 2, Rational(2)});
                const Outcome outcome = solve_and_check(planted, method, "twelve pairs and three real roots" + by);
                require(outcome.regions == planted.size(), "twelve pairs: a region holds more than one root");
                require_pairs(outcome, 12, "twelve pairs");
            }
            // 2. The same roots turned by (3 + 4i)/5: complex coefficients, nothing to reflect.
            {
                std::vector<Complex> planted;
                const Rational c = Rational(3) / 5, s = Rational(4) / 5;
                for (int k = 1; k <= 12; ++k) {
                    const Rational x = Rational(3 * k - 20) / 7, y = Rational(2 * k + 1) / 9;
                    planted.push_back(Complex(c * x - s * y, s * x + c * y));
                    planted.push_back(Complex(c * x + s * y, s * x - c * y));
                }
                const Outcome outcome = solve_and_check(planted, method, "the pairs turned off the axis" + by);
                require(outcome.reflected == 0, "complex coefficients: a region was reflected");
            }
            // 3. Round coordinates: i, 1 + i, -2 + i, (1 + i)/2, 2 + 3i, and real roots at 3 and -3.
            {
                const auto planted = pairs_and_real({{Rational(0), Rational(1)},
                                                     {Rational(1), Rational(1)},
                                                     {Rational(-2), Rational(1)},
                                                     {Rational(1) / 2, Rational(1) / 2},
                                                     {Rational(2), Rational(3)},
                                                     {Rational(-1), Rational(2)}},
                                                    {Rational(3), Rational(-3)});
                require_pairs(solve_and_check(planted, method, "round coordinates" + by), 6, "round coordinates");
            }
            // 4. A lattice of Gaussian half-integers, a vertical line and a horizontal line of roots: many roots share
            //    each coordinate.
            {
                std::vector<std::pair<Rational, Rational>> lattice, vertical, horizontal;
                for (int a = -2; a <= 2; ++a)
                    for (int b = 1; b <= 3; ++b)
                        lattice.push_back({Rational(a) / 2, Rational(b) / 2});
                for (int k = 1; k <= 10; ++k) {
                    vertical.push_back({Rational(1) / 4, Rational(3 * k) / 8});
                    horizontal.push_back({Rational(3 * k) / 8 - 2, Rational(1)});
                }
                require_pairs(solve_and_check(pairs_and_real(lattice, {Rational(1) / 8}), method, "a lattice" + by), 15,
                              "a lattice");
                require_pairs(solve_and_check(pairs_and_real(vertical, {}), method, "a vertical line" + by), 10,
                              "a vertical line");
                require_pairs(solve_and_check(pairs_and_real(horizontal, {}), method, "a horizontal line" + by), 10,
                              "a horizontal line");
            }
            // 5. A pair 2^-40 from the axis beside ordinary pairs and a real root: whatever is reflected is right.
            {
                const Rational tiny = Rational(1) / Rational("1099511627776");
                std::vector<std::pair<Rational, Rational>> upper = {{Rational(1) / 3, tiny}};
                for (int k = 1; k <= 6; ++k)
                    upper.push_back({Rational(k) / 5, Rational(k + 2) / 7});
                const Outcome outcome =
                    solve_and_check(pairs_and_real(upper, {Rational(-5) / 4}), method, "a pair close to the axis" + by);
                require(outcome.reflected <= 7, "a pair close to the axis: more regions reflected than pairs");
                reflected += outcome.reflected;
            }
            // 6. Two pairs 2^-30 apart, closer than the pairing's tolerance: no pair is guessed for those four
            //    points, and the cover is right.
            {
                const Rational gap = Rational(1) / Rational("1073741824");
                std::vector<std::pair<Rational, Rational>> upper = {{Rational(2) / 3, Rational(3) / 4},
                                                                    {Rational(2) / 3 + gap, Rational(3) / 4}};
                for (int k = 1; k <= 5; ++k)
                    upper.push_back({Rational(-k) / 4, Rational(k) / 3});
                const Outcome outcome =
                    solve_and_check(pairs_and_real(upper, {}), method, "two pairs closer than the tolerance" + by);
                require(outcome.reflected <= 5, "two pairs closer than the tolerance: one of them was reflected");
                reflected += outcome.reflected;
            }
        }
        require(reflected > 0, "no region was reflected in any solve");
        std::cout << "conjugate reflection tests passed (" << reflected << " reflected regions)\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "FAILED: " << error.what() << '\n';
        return 1;
    }
}
