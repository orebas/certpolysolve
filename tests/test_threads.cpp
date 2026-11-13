// Threads (the 0.5 API's thread guarantees, docs/api-spec-0.5.md section 8 and
// its O17): solves run at the same time on many threads; one input and one set of options may be shared read-only
// by all of them; and results (and inputs) may move between threads and be destroyed on any thread, also after the
// thread that made them has run its cleanup (cps_cleanup_thread, here through ThreadScope). Each thread's answer is
// compared with the same solve on the main thread, region by region.
// Wider cases (19 October 2026; the thread-safety audit's F1, docs/thread-safety-audit-20261018.md): every method at
// once on one shared input; a goal above 4,592 bits, where FLINT keeps pi and log 2 in thread-local caches; one
// program shared by every thread; and one cancellation token stopping many solves at once, beside solves stopped by
// their own time limits. The tsan preset runs this file under ThreadSanitizer, against a FLINT built with it.
#include "exact_reference_api.hpp"
#include <certpolysolve/certpolysolve.hpp>
#include <atomic>
#include <chrono>
#include <iostream>
#include <memory>
#include <string>
#include <thread>
#include <vector>

using namespace reference;

namespace {
constexpr int threads = 8;

Polynomial test_input() {
    // (x^2 - 2)(x^2 + 1)(3x - 1)(x - 7/5)(5x^2 + 3x + 9): real and complex roots, distinct
    return multiply(multiply(multiply(multiply({-2, 0, 1}, {1, 0, 1}), {-1, 3}), {Rational(-7) / 5, 1}), {9, 3, 5});
}

std::string signature(const certpolysolve::Result &r) {
    std::string s = std::to_string(cps_result_termination(r.get())) + "|";
    for (const auto &e : r.roots())
        s += e.box.real.lo.mantissa() + "@" + e.box.real.lo.exponent() + "," + e.box.real.hi.mantissa() + "@" +
             e.box.real.hi.exponent() + "," + e.box.imag.lo.mantissa() + "@" + e.box.imag.lo.exponent() + "," +
             e.box.imag.hi.mantissa() + "@" + e.box.imag.hi.exponent() + "#" + std::to_string(e.count) + ";";
    return s;
}

// The complex roots, or the real roots only.
certpolysolve::Result solve_for(const certpolysolve::Polynomial &p, const certpolysolve::Options &options, bool real_only) {
    return real_only ? certpolysolve::real_roots(p, options) : certpolysolve::solve(p, options);
}

// One input and one options object, shared read-only by every thread, each solving it several times at once.
void shared_read_only(const certpolysolve::Polynomial &p, const certpolysolve::Options &options, bool real_only,
                      const std::string &expected) {
    std::vector<std::string> seen(threads);
    std::vector<std::thread> pool;
    for (int t = 0; t < threads; ++t)
        pool.emplace_back([&, t] {
            certpolysolve::ThreadScope scope;
            for (int round = 0; round < 3; ++round) {
                const auto result = solve_for(p, options, real_only);
                const auto s = signature(result);
                if (round == 0) seen[t] = s;
                else if (s != seen[t]) seen[t] = "changed between rounds";
            }
        });
    for (auto &thread : pool) thread.join();
    for (int t = 0; t < threads; ++t)
        require(seen[t] == expected, ("thread " + std::to_string(t) + " saw another answer: " + seen[t]).c_str());
}

// Results made on worker threads, read and destroyed on this one after every worker ran its cleanup; inputs made
// on this thread, destroyed on the workers.
void across_threads(const Polynomial &exact, const certpolysolve::Options &options, bool real_only,
                    const std::string &expected) {
    std::vector<std::unique_ptr<certpolysolve::Result>> results(threads);
    std::vector<std::unique_ptr<certpolysolve::Polynomial>> inputs;
    for (int t = 0; t < threads; ++t) inputs.push_back(std::make_unique<certpolysolve::Polynomial>(input(exact)));
    std::vector<std::thread> pool;
    for (int t = 0; t < threads; ++t)
        pool.emplace_back([&, t] {
            certpolysolve::ThreadScope scope; // its cleanup runs when the worker ends, before the results are read
            results[t] = std::make_unique<certpolysolve::Result>(solve_for(*inputs[t], options, real_only));
            inputs[t].reset(); // an input made on the main thread, destroyed here
        });
    for (auto &thread : pool) thread.join();
    for (int t = 0; t < threads; ++t) {
        require(signature(*results[t]) == expected, ("a result read on another thread differs: " + std::to_string(t)).c_str());
        results[t].reset(); // destroyed on the main thread, after its worker's cleanup
    }
}

bool succeeded(const std::string &signature) { return signature.rfind(std::to_string(CPS_SUCCESS) + "|", 0) == 0; }

// Every thread runs solve_on(t) `rounds` times, each from its own ThreadScope; each round must give expected(t).
template <class Solve, class Expected>
void on_threads(const Solve &solve_on, const Expected &expected, int rounds, const std::string &what) {
    std::vector<std::string> seen(threads);
    std::vector<std::thread> pool;
    for (int t = 0; t < threads; ++t)
        pool.emplace_back([&, t] {
            certpolysolve::ThreadScope scope;
            for (int round = 0; round < rounds; ++round) {
                const auto s = signature(solve_on(t));
                if (round == 0) seen[t] = s;
                else if (s != seen[t]) seen[t] = "changed between rounds";
            }
        });
    for (auto &thread : pool) thread.join();
    for (int t = 0; t < threads; ++t)
        require(seen[t] == expected(t), (what + ": thread " + std::to_string(t) + " saw another answer: " + seen[t]).c_str());
}

// Integer coefficients in a fixed pattern in [-500, 499]; at degree 1,000 this is test_cancellation's slow input.
certpolysolve::Polynomial patterned(uint64_t degree) {
    certpolysolve::Polynomial p;
    for (uint64_t k = 0; k <= degree; ++k)
        p.set(k, static_cast<double>(static_cast<int>((k * 7919 + 13) % 1000) - 500));
    return p;
}

// One degree-64 input solved by the four methods at the same time: two threads a method, each twice, against the
// method's answer on this thread.
void every_method_at_once() {
    const auto p = patterned(64);
    const std::vector<certpolysolve::Method> methods = {certpolysolve::Method::automatic, certpolysolve::Method::classical,
                                                    certpolysolve::Method::secular, certpolysolve::Method::piecewise};
    const int count = static_cast<int>(methods.size());
    std::vector<certpolysolve::Options> options(methods.size());
    std::vector<std::string> expected;
    for (int m = 0; m < count; ++m) {
        options[m].bits = 100;
        options[m].max_bits = 4096;
        options[m].method = methods[m];
        expected.push_back(signature(certpolysolve::solve(p, options[m])));
        require(succeeded(expected.back()), ("method " + std::to_string(static_cast<int>(methods[m])) + " did not succeed").c_str());
    }
    on_threads([&](int t) { return certpolysolve::solve(p, options[t % count]); },
               [&](int t) { return expected[t % count]; }, 2, "every method at once");
}

// A goal of 6,000 bits, above the 4,592 below which FLINT reads pi and log 2 from static tables: above it they come
// from thread-local caches built at the highest precision asked for so far (the audit's B3). Every solve runs on a
// fresh thread, the reference too, so the caches start alike and the answers must agree bit for bit.
void high_precision_on_fresh_threads(const certpolysolve::Polynomial &p) {
    certpolysolve::Options options;
    options.bits = 6000;
    options.max_bits = 65536;
    std::string expected;
    std::thread([&] {
        certpolysolve::ThreadScope scope;
        expected = signature(certpolysolve::solve(p, options));
    }).join();
    require(succeeded(expected), "the 6,000-bit reference solve did not succeed");
    on_threads([&](int) { return certpolysolve::solve(p, options); }, [&](int) { return expected; }, 1, "6,000 bits");
}

// One program shared read-only by every thread. Its output is a product of distinct factors (Mandelbrot's p_4 times
// the Chebyshev polynomial T_8, both by their recurrences), so the product lift runs.
void one_program_on_threads() {
    const auto program = certpolysolve::Program::from([](auto x) {
        auto mandelbrot = x * 0 + 1;
        for (int j = 0; j < 4; ++j) mandelbrot = x * mandelbrot * mandelbrot + 1;
        auto previous = x * 0 + 1, chebyshev = x;
        for (int k = 1; k < 8; ++k) {
            auto next = 2 * x * chebyshev - previous;
            previous = chebyshev;
            chebyshev = next;
        }
        return mandelbrot * chebyshev;
    });
    certpolysolve::Options options;
    options.bits = 100;
    options.max_bits = 4096;
    const auto reference = certpolysolve::solve(program, options);
    const auto expected = signature(reference);
    require(succeeded(expected), "the program's reference solve did not succeed");
    // The degree checks the recording itself (p_4 has degree 15, T_8 degree 8): every thread solves the same program,
    // so their agreement alone would pass a wrong one. Under GCC the recorder once read each node before the call
    // that wrote it and recorded the variable alone (19 October 2026).
    const std::string degree = "the recorded program has degree " + std::to_string(reference.degree()) + ", not 23";
    require(reference.degree() == 15 + 8, degree.c_str());
    on_threads([&](int) { return certpolysolve::solve(program, options); }, [&](int) { return expected; }, 2, "one program");
}

// One cancellation token in the options of half the solves, requested from this thread once every solve has started;
// the other half stop at their own time limit of 50 ms. The whole solve takes about 2 s, so each must stop as asked,
// and an answer that says it is complete must count every root.
void one_token_many_solves() {
    constexpr uint64_t degree = 1000;
    const auto p = patterned(degree);
    std::stop_source token;
    certpolysolve::Options by_token, by_clock;
    by_token.bits = 1000;
    by_token.max_bits = 65536;
    by_token.stop = token.get_token();
    by_clock.bits = 1000;
    by_clock.max_bits = 65536;
    by_clock.time_limit = std::chrono::duration<double>(0.05);
    std::atomic<int> started{0};
    std::vector<std::unique_ptr<certpolysolve::Result>> results(threads);
    std::vector<std::thread> pool;
    for (int t = 0; t < threads; ++t)
        pool.emplace_back([&, t] {
            certpolysolve::ThreadScope scope;
            started.fetch_add(1);
            results[t] = std::make_unique<certpolysolve::Result>(certpolysolve::solve(p, t % 2 ? by_clock : by_token));
        });
    while (started.load() < threads) std::this_thread::yield();
    std::this_thread::sleep_for(std::chrono::milliseconds(20));
    token.request_stop();
    for (auto &thread : pool) thread.join();
    for (int t = 0; t < threads; ++t) {
        const auto &result = *results[t];
        const auto reason = t % 2 ? CPS_STOP_TIME_LIMIT : CPS_STOP_CANCELLED;
        require(cps_result_termination(result.get()) == CPS_CANCELLED && cps_result_stop_reason(result.get()) == reason,
                ("thread " + std::to_string(t) + " did not stop as asked: " + result.message()).c_str());
        if (result.complete()) {
            uint64_t roots = 0;
            for (const auto &e : result.roots()) roots += e.count;
            require(roots == degree, ("thread " + std::to_string(t) + ": a complete answer misses roots").c_str());
        }
    }
}
} // namespace

int main() {
    certpolysolve::ThreadScope lifetime;
    try {
        const auto exact = test_input();
        const auto p = input(exact);
        for (const bool real_only : {false, true}) {
            certpolysolve::Options options;
            options.bits = 200;
            options.max_bits = 4096;
            const auto expected = signature(solve_for(p, options, real_only));
            require(expected.rfind(std::to_string(CPS_SUCCESS) + "|", 0) == 0, "the reference solve did not succeed");
            shared_read_only(p, options, real_only, expected);
            across_threads(exact, options, real_only, expected);
        }
        every_method_at_once();
        high_precision_on_fresh_threads(p);
        one_program_on_threads();
        one_token_many_solves();
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
    std::cout << "thread tests passed\n";
    return 0;
}
