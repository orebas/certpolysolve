// A solve runs in the default floating-point environment, whatever the calling thread has set, and gives the caller
// its own environment back (docs/api-spec-0.5.md section 8.2, change C28; the thread-safety audit's C3 and E4). Each
// method solves one input under the default environment, then under upward rounding with FE_INEXACT already raised
// and, on x86-64, flush-to-zero and denormals-are-zero on. The two answers must agree region by region, and after the
// call the rounding mode, the raised flags (only the caller's) and the MXCSR bits must be the caller's.
#include "exact_reference_api.hpp"
#include <certpolysolve/certpolysolve.hpp>
#include <cfenv>
#include <iostream>
#include <string>
#include <vector>
#if defined(__x86_64__) || defined(_M_X64)
#include <xmmintrin.h>
#define CERTPOLYSOLVE_TEST_MXCSR 1
#endif

using namespace reference;

namespace {
std::string signature(const certpolysolve::Result &r) {
    std::string s = std::to_string(cps_result_termination(r.get())) + "|";
    for (const auto &e : r.roots())
        s += e.box.real.lo.mantissa() + "@" + e.box.real.lo.exponent() + "," + e.box.real.hi.mantissa() + "@" +
             e.box.real.hi.exponent() + "," + e.box.imag.lo.mantissa() + "@" + e.box.imag.lo.exponent() + "," +
             e.box.imag.hi.mantissa() + "@" + e.box.imag.hi.exponent() + "#" + std::to_string(e.count) + ";";
    return s;
}

// Integer coefficients in a fixed pattern in [-500, 499] (test_threads' input).
certpolysolve::Polynomial patterned(uint64_t degree) {
    certpolysolve::Polynomial p;
    for (uint64_t k = 0; k <= degree; ++k)
        p.set(k, static_cast<double>(static_cast<int>((k * 7919 + 13) % 1000) - 500));
    return p;
}

#ifdef CERTPOLYSOLVE_TEST_MXCSR
// MXCSR's flush-to-zero and denormals-are-zero bits.
constexpr unsigned flush_to_zero = 1u << 15, denormals_are_zero = 1u << 6;
#endif

// What the caller's environment looked like right after a call.
struct Seen {
    int rounding = 0, flags = 0;
    unsigned mxcsr = 0;
};

void method_under_caller_environment(const certpolysolve::Polynomial &p, certpolysolve::Method method) {
    certpolysolve::Options options;
    options.bits = 100;
    options.max_bits = 4096;
    options.method = method;
    const auto expected = signature(certpolysolve::solve(p, options));
    std::fenv_t own;
    require(std::fegetenv(&own) == 0, "fegetenv failed");
    std::fesetround(FE_UPWARD);
    std::feclearexcept(FE_ALL_EXCEPT);
    std::feraiseexcept(FE_INEXACT);
#ifdef CERTPOLYSOLVE_TEST_MXCSR
    _mm_setcsr(_mm_getcsr() | flush_to_zero | denormals_are_zero);
#endif
    const auto result = certpolysolve::solve(p, options);
    Seen seen;
    seen.rounding = std::fegetround();
    seen.flags = std::fetestexcept(FE_ALL_EXCEPT);
#ifdef CERTPOLYSOLVE_TEST_MXCSR
    seen.mxcsr = _mm_getcsr();
#endif
    std::fesetenv(&own);
    const std::string name = "method " + std::to_string(static_cast<int>(method)) + ": ";
    require(signature(result) == expected, (name + "the answer depends on the caller's environment").c_str());
    require(seen.rounding == FE_UPWARD, (name + "the caller's rounding mode was not given back").c_str());
    require(seen.flags == FE_INEXACT, (name + "the caller's exception flags were not given back as they were").c_str());
#ifdef CERTPOLYSOLVE_TEST_MXCSR
    require((seen.mxcsr & (flush_to_zero | denormals_are_zero)) == (flush_to_zero | denormals_are_zero),
            (name + "the caller's flush-to-zero and denormals-are-zero bits were not given back").c_str());
#endif
}
} // namespace

int main() {
    certpolysolve::ThreadScope lifetime;
    try {
        const auto p = patterned(64);
        for (const certpolysolve::Method method : {certpolysolve::Method::automatic, certpolysolve::Method::classical,
                                                certpolysolve::Method::secular, certpolysolve::Method::piecewise})
            method_under_caller_environment(p, method);
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
    std::cout << "floating-point environment tests passed\n";
    return 0;
}
