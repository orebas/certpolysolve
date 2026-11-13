// A time limit, and (advanced) stopping a solve from another thread.
#include <certpolysolve/certpolysolve.hpp>
#include <chrono>
#include <print>
#include <stop_token>
#include <thread>

using namespace std::chrono_literals;

certpolysolve::Polynomial slow_input() { // degree 1,000; at 1,000 bits its whole solve takes seconds
    certpolysolve::Polynomial p;
    for (uint64_t k = 0; k <= 1000; ++k) p.set(k, static_cast<int>((k * 7919 + 13) % 1000) - 500);
    return p;
}

void report(const certpolysolve::Polynomial &p, const certpolysolve::Options &options) {
    const auto start = std::chrono::steady_clock::now();
    const certpolysolve::Result result = certpolysolve::solve(p, options);
    const std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - start;

    uint64_t roots_in_boxes = 0;
    for (const certpolysolve::Root &root : result.roots()) roots_in_boxes += root.count;
    std::println("stopped ({}) after {:.3f} s: {} boxes holding {} of {} roots, every root covered: {}", result.stop(),
                 elapsed.count(), result.roots().size(), roots_in_boxes, result.degree(), result.complete() ? "yes" : "no");
}

int main() {
    report(slow_input(), {.bits = 1000, .time_limit = 100ms}); // at most 0.1 s of wall-clock time

    // Advanced: a std::stop_token; the solve stops at its next checkpoint after another thread asks.
    std::stop_source source;
    certpolysolve::Options stoppable;
    stoppable.bits = 1000;
    stoppable.stop = source.get_token();
    std::jthread stopper([&] {
        std::this_thread::sleep_for(50ms);
        source.request_stop();
    });
    report(slow_input(), stoppable);
}
