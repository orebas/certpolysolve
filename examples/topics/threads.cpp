// Many solves at once: the stability sweep of stability.cpp over 400 gains, on four threads.
#include <certpolysolve/certpolysolve.hpp>
#include <atomic>
#include <print>
#include <thread>
#include <vector>

int verdict(double gain) { // 1 proved stable, -1 proved unstable, 0 undecided
    const certpolysolve::Result result = certpolysolve::solve(certpolysolve::Polynomial{1 + gain, 3, 3, 1});
    if (!result.ok()) return 0;
    bool all_left = true;
    for (const certpolysolve::Root &root : result.roots()) {
        if (root.box.real.lo > 0) return -1;
        if (root.box.real.hi >= 0) all_left = false;
    }
    return all_left ? 1 : 0;
}

int main() {
    const int gains = 400, threads = 4;
    const double first = 6, step = 0.01; // K = 6.00, 6.01, ..., 9.99
    std::vector<int> verdicts(gains);
    std::atomic<int> next{0};

    {
        std::vector<std::jthread> pool;
        for (int t = 0; t < threads; ++t)
            pool.emplace_back([&] {
                certpolysolve::ThreadScope scope; // first thing in the thread: frees its FLINT caches when it ends
                for (int i; (i = next++) < gains;) verdicts[i] = verdict(first + step * i);
            });
    } // the jthreads join here

    for (int i = 0; i < gains; ++i)
        if (verdicts[i] != 1) {
            std::println("stable proved up to K = {:.2f}; at K = {:.2f} {}", first + step * (i - 1), first + step * i,
                         verdicts[i] < 0 ? "unstable" : "undecided");
            break;
        }
}
