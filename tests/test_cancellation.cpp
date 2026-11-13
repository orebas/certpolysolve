// The time limit and the cancellation token (the 0.5 API's stage 1, 18 October 2026: added beside the cancellation
// callback, which a later stage removes), the version macros, and cps_last_error, which a successful call clears.
// The slow input is a degree-1,000 polynomial at 1,000 bits of accuracy: its whole solve takes about 2 s on the
// development machine, so a limit of 50 ms, or a request 20 ms after the start, falls inside it; no test here runs
// it to the end. The bounds on elapsed time are generous (a solve stops at its next checkpoint, after the longest
// step that cannot be interrupted).
#include "exact_reference_api.hpp"
#include <certpolysolve/certpolysolve.hpp>
#include <chrono>
#include <cmath>
#include <cstring>
#include <iostream>
#include <limits>
#include <string>
#include <thread>

using namespace reference;

namespace {
using Clock = std::chrono::steady_clock;
double seconds_since(Clock::time_point start) { return std::chrono::duration<double>(Clock::now() - start).count(); }

certpolysolve::Polynomial slow_input() {
    certpolysolve::Polynomial p;
    for (uint64_t k = 0; k <= 1000; ++k)
        p.set(k, static_cast<double>(static_cast<int>((k * 7919 + 13) % 1000) - 500));
    return p;
}
certpolysolve::Options slow_options() {
    certpolysolve::Options options;
    options.bits = 1000;
    options.max_bits = 65536;
    return options;
}

void versions() {
    require(std::strcmp(cps_version(), CPS_VERSION_STRING) == 0, "the library's version is not the headers'");
    require(std::strcmp(certpolysolve::version(), CPS_VERSION_STRING) == 0, "the C++ version is not the C one");
    require(CPS_VERSION_MAJOR == 0 && CPS_VERSION_MINOR == 5, "the version macros");
}

void last_error_is_cleared() {
    cps_options *o = nullptr;
    require(cps_options_create(&o) == CPS_OK, "options");
    require(cps_options_set_time_limit(o, -1.0) == CPS_INVALID_ARGUMENT, "a negative time limit was accepted");
    require(std::strlen(cps_last_error()) > 0, "a failed call left no message");
    require(cps_options_set_time_limit(o, 1.0) == CPS_OK, "a positive time limit was refused");
    require(std::strlen(cps_last_error()) == 0, "a successful call left the previous message");
    require(cps_options_set_time_limit(o, std::nan("")) == CPS_INVALID_ARGUMENT, "a NaN time limit was accepted");
    require(cps_options_set_time_limit(o, std::numeric_limits<double>::infinity()) == CPS_OK &&
                cps_options_set_time_limit(o, 0.0) == CPS_OK,
            "an infinite or zero time limit (none) was refused");
    cps_options_destroy(o);
}

// A limit that the solve does not reach leaves the answer as it is without one.
void limit_not_reached() {
    const Polynomial p = multiply(multiply({-2, 0, 1}, {1, 0, 1}), {-1, 3}); // (x^2 - 2)(x^2 + 1)(3x - 1)
    const auto source = input(p);
    const auto a = certpolysolve::solve(source);
    for (double seconds : {60.0, 1e10, 1e100, std::numeric_limits<double>::max(),
                           std::numeric_limits<double>::infinity()}) {
        certpolysolve::Options limited;
        limited.time_limit = std::chrono::duration<double>(seconds);
        const auto b = certpolysolve::solve(source, limited);
        require(cps_result_termination(b.get()) == CPS_SUCCESS && b.goal_met(),
                "a solve under an unreached limit did not succeed");
        const auto &ra = a.roots(), &rb = b.roots();
        require(ra.size() == rb.size(), "an unreached limit changed the answer");
        for (size_t i = 0; i < ra.size(); ++i)
            require(endpoint(ra[i].box.real.lo) == endpoint(rb[i].box.real.lo) &&
                        endpoint(ra[i].box.real.hi) == endpoint(rb[i].box.real.hi) &&
                        endpoint(ra[i].box.imag.lo) == endpoint(rb[i].box.imag.lo) &&
                        endpoint(ra[i].box.imag.hi) == endpoint(rb[i].box.imag.hi) && ra[i].count == rb[i].count,
                    "an unreached limit changed a region");
    }
}

void limit_reached() {
    auto options = slow_options();
    options.time_limit = std::chrono::duration<double>(0.05);
    const auto start = Clock::now();
    const auto result = certpolysolve::solve(slow_input(), options);
    const double elapsed = seconds_since(start);
    require(cps_result_termination(result.get()) == CPS_CANCELLED &&
                cps_result_stop_reason(result.get()) == CPS_STOP_TIME_LIMIT && !result.goal_met(),
            ("the time limit did not stop the solve: " + result.message()).c_str());
    require(elapsed < 1.5, ("the solve ran " + std::to_string(elapsed) + " s past a 50 ms limit").c_str());
    // The real-only route's solves inherit the call's deadline.
    auto real = slow_options();
    real.time_limit = std::chrono::duration<double>(0.05);
    const auto real_result = certpolysolve::real_roots(slow_input(), real);
    require(cps_result_termination(real_result.get()) == CPS_CANCELLED &&
                cps_result_stop_reason(real_result.get()) == CPS_STOP_TIME_LIMIT,
            ("the time limit did not stop the real-only solve: " + real_result.message()).c_str());
}

void token_requested_before() {
    std::stop_source token;
    token.request_stop();
    require(token.stop_requested(), "a requested token does not say so");
    auto options = slow_options();
    options.stop = token.get_token();
    const auto start = Clock::now();
    const auto result = certpolysolve::solve(slow_input(), options);
    require(cps_result_termination(result.get()) == CPS_CANCELLED &&
                cps_result_stop_reason(result.get()) == CPS_STOP_CANCELLED,
            ("a requested token did not stop the solve: " + result.message()).c_str());
    require(seconds_since(start) < 1.5, "a solve whose token was requested before it started ran on");
}

void token_requested_from_another_thread() {
    std::stop_source token;
    auto options = slow_options();
    options.stop = token.get_token();
    std::thread requester([&] {
        std::this_thread::sleep_for(std::chrono::milliseconds(20));
        token.request_stop();
    });
    const auto start = Clock::now();
    const auto result = certpolysolve::solve(slow_input(), options);
    const double elapsed = seconds_since(start);
    requester.join();
    require(cps_result_termination(result.get()) == CPS_CANCELLED &&
                cps_result_stop_reason(result.get()) == CPS_STOP_CANCELLED,
            ("a token requested from another thread did not stop the solve: " + result.message()).c_str());
    require(elapsed < 1.5, ("the solve ran " + std::to_string(elapsed) + " s after the request").c_str());
}

// Options share the token: destroying the caller's stop source leaves it in force for the options that name it.
void token_outlives_its_handle() {
    auto options = slow_options();
    {
        std::stop_source token;
        options.stop = token.get_token();
        token.request_stop();
    }
    const auto result = certpolysolve::solve(slow_input(), options);
    require(cps_result_termination(result.get()) == CPS_CANCELLED, "options lost the token when its handle was destroyed");
    options.stop = std::stop_token{};
    options.time_limit = std::chrono::duration<double>(0.05);
    const auto cleared = certpolysolve::solve(slow_input(), options);
    require(cps_result_stop_reason(cleared.get()) == CPS_STOP_TIME_LIMIT, "a cleared token still cancelled the solve");
}
} // namespace

int main() {
    certpolysolve::ThreadScope lifetime;
    try {
        versions();
        last_error_is_cleared();
        limit_not_reached();
        limit_reached();
        token_requested_before();
        token_requested_from_another_thread();
        token_outlives_its_handle();
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
    std::cout << "cancellation tests passed\n";
    return 0;
}
