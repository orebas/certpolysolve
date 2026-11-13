// A real-only query whose non-real roots lie closer to the real axis than the complex route's certificates
// resolve at the goal: their boxes meet the axis without being proved real, and the real-only route refines
// them until they are proved one way or the other, within the precision cap (the axis fuzz family's finding of
// 16 October 2026, seed 15841096042299828660: pairs -1 +- i 2^-7 and 9/32 +- i 2^-24 beside the real root
// -2101/8192, a 16-bit goal, a cap of 128 bits; the route stopped with CPS_PRECISION_LIMIT at its first
// accuracy, 96 working bits, whatever the cap). Checked on that input: the answer under a cap of 128 bits (one
// real root, complete) and a session that stops under the smallest cap the API allows (32 bits) and answers
// once the cap is raised; and, with the closer pair at 2^-40 (the complex route needs 192 working bits for it),
// the honest precision stop under the 32-bit cap.
#include "exact_winding_api.hpp"
#include <certpolysolve/certpolysolve.hpp>
#include <iostream>

using namespace reference;

namespace {
certpolysolve::Options query(uint64_t bits, uint64_t cap) {
    certpolysolve::Options options;
    options.method = certpolysolve::Method::automatic;
    options.goal = certpolysolve::Goal::accuracy();
    options.bits = bits;
    options.max_bits = cap;
    return options;
}
certpolysolve::Polynomial input_with_pair_at(int k) {
    const Rational one_128 = Rational(1) / Rational(128), small = Rational(1) / Rational(Integer(1) << k);
    const std::vector<Complex> roots{Complex(Rational(-1), one_128), Complex(Rational(-1), -one_128),
                                     Complex(Rational("9/32"), small), Complex(Rational("9/32"), -small),
                                     Complex(Rational("-2101/8192"), Rational(0))};
    return reference::input(from_roots(roots));
}
void check_answer(const certpolysolve::Result &result, const char *what) {
    require(result.goal_met() && result.complete() && cps_result_termination(result.get()) == CPS_SUCCESS &&
                cps_result_stop_reason(result.get()) == CPS_STOP_NONE,
            what);
    const auto &regions = result.roots();
    require(regions.size() == 1 && regions[0].count == 1 && regions[0].realness == certpolysolve::Realness::real,
            "the real-only answer is not the one real root");
    const Rational root("-2101/8192");
    require(endpoint(regions[0].box.real.lo) <= root && root <= endpoint(regions[0].box.real.hi),
            "the real root's interval misses it");
}
void check_precision_stop(const certpolysolve::Result &result, const char *what) {
    require(!result.goal_met() && !result.complete() && cps_result_termination(result.get()) == CPS_PRECISION_LIMIT &&
                cps_result_stop_reason(result.get()) == CPS_STOP_PRECISION,
            what);
}
} // namespace

int main() {
    certpolysolve::ThreadScope lifetime;
    try {
        const auto finding = input_with_pair_at(24);
        // The finding: a 16-bit goal, a 128-bit cap; the route refines past its first accuracy and answers.
        check_answer(certpolysolve::real_roots(finding, query(16, 128)), "the real-only query under a 128-bit cap did not answer");
        // The smallest cap answers it or stops honestly for the precision.
        {
            const auto first = certpolysolve::real_roots(finding, query(16, 32));
            if (!first.goal_met())
                check_precision_stop(first, "the query under a 32-bit cap neither answered nor stopped for the precision");
        }
        // A pair the smallest cap cannot resolve: an honest precision stop, not an answer.
        check_precision_stop(certpolysolve::real_roots(input_with_pair_at(40), query(16, 32)),
                             "the real-only query on the 2^-40 pair under a 32-bit cap did not stop for the precision");
        std::cout << "real-only realness tests passed\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "FAILED: " << error.what() << '\n';
        return 1;
    }
}
