// W9's first cut, its fork and join (src/platform/parallel.hpp): without a team parallel_for runs its blocks in order on the
// calling thread; with one (TeamScope) every block runs exactly once and the call returns after the last; a
// parallel_for inside a block runs serially and in order; an exception from a block is rethrown after the join, the
// lowest block's, every block having run; a worker's block runs in the default floating-point environment and the
// flags it raises reach the calling thread; a worker starts with every signal blocked and the calling thread keeps
// its mask; a job at a split site counts once, on a team only; a nested TeamScope keeps the outer team; four teams
// on four threads run at once; and a team survives thousands of short jobs (the race and lost-wakeup shake-out, also
// under ThreadSanitizer).
#include "platform/parallel.hpp"
#include <atomic>
#include <cfenv>
#include <iostream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>
#if defined(__unix__) || defined(__APPLE__)
#include <pthread.h>
#include <signal.h>
#endif

namespace d = certpolysolve::detail;

namespace {
void require(bool condition, const std::string &message) {
    if (!condition) throw std::runtime_error(message);
}

void serial_without_a_team() {
    require(d::parallel_width() == 1, "a width without a team");
    std::vector<size_t> order;
    d::parallel_for(5, [&](size_t b) { order.push_back(b); });
    require(order == std::vector<size_t>{0, 1, 2, 3, 4}, "blocks out of order without a team");
}

void each_block_once(unsigned workers, size_t blocks) {
    const d::TeamScope team(workers);
    require(d::parallel_width() == workers, "the team's width");
    std::vector<std::atomic<int>> runs(blocks);
    std::vector<size_t> out(blocks, 0);
    d::parallel_for(blocks, [&](size_t b) {
        runs[b].fetch_add(1, std::memory_order_relaxed);
        out[b] = b * b + 1;
    });
    for (size_t b = 0; b < blocks; ++b)
        require(runs[b].load() == 1 && out[b] == b * b + 1, "block " + std::to_string(b) + " ran " +
                                                                 std::to_string(runs[b].load()) + " times");
}

void nested_runs_serially() {
    const d::TeamScope team(4);
    std::vector<std::vector<size_t>> inner(8);
    std::vector<unsigned> widths(8, 0);
    d::parallel_for(8, [&](size_t b) {
        widths[b] = d::parallel_width();
        d::parallel_for(6, [&](size_t i) { inner[b].push_back(i); });
    });
    for (size_t b = 0; b < 8; ++b)
        require(widths[b] == 1 && inner[b] == std::vector<size_t>{0, 1, 2, 3, 4, 5}, "a nested loop was split");
}

void exceptions_come_home() {
    const d::TeamScope team(4);
    std::atomic<int> ran{0};
    bool caught = false;
    try {
        d::parallel_for(16, [&](size_t b) {
            ran.fetch_add(1);
            if (b == 3 || b == 11) throw std::runtime_error("block " + std::to_string(b));
        });
    } catch (const std::runtime_error &e) {
        caught = std::string(e.what()) == "block 3";
    }
    require(caught && ran.load() == 16, "the lowest block's exception after every block ran");
    // The team still works after an exception.
    std::atomic<int> after{0};
    d::parallel_for(10, [&](size_t) { after.fetch_add(1); });
    require(after.load() == 10, "the team after an exception");
}

// A worker's block runs in the default floating-point environment, the solve's thread's blocks in that thread's own
// (a solve call sets it to the default too); a flag raised in any block is raised on the calling thread after the
// join, as after a serial loop, and a loop that raises nothing leaves the flags clear.
void floating_point_environment() {
    const d::TeamScope team(4);
    const int caller = std::fegetround();
    const auto here = std::this_thread::get_id();
    std::fesetround(FE_UPWARD);
    std::vector<int> modes(48, -1), own(48, 0);
    d::parallel_for(48, [&](size_t b) {
        modes[b] = std::fegetround();
        own[b] = std::this_thread::get_id() == here;
    });
    std::fesetround(caller);
    for (size_t b = 0; b < 48; ++b)
        require(modes[b] == (own[b] ? FE_UPWARD : FE_TONEAREST), "a block ran in the wrong rounding mode");
    std::feclearexcept(FE_ALL_EXCEPT);
    d::parallel_for(16, [&](size_t b) {
        if (b == 5) std::feraiseexcept(FE_DIVBYZERO);
    });
    require(std::fetestexcept(FE_DIVBYZERO) != 0, "a block's flag was lost at the join");
    std::feclearexcept(FE_ALL_EXCEPT);
    d::parallel_for(16, [&](size_t) {});
    require(std::fetestexcept(FE_ALL_EXCEPT) == 0, "a loop that raised nothing left a flag");
}

#if defined(__unix__) || defined(__APPLE__)
// The blocks that run on a worker see SIGINT, SIGTERM, SIGALRM and SIGPROF blocked (signals a host or its runtime
// may use); the calling thread's mask is the same after the team as before.
void workers_block_signals() {
    sigset_t before, after;
    pthread_sigmask(SIG_SETMASK, nullptr, &before);
    const auto here = std::this_thread::get_id();
    std::vector<int> blocked(64, 0), own(64, 0);
    {
        const d::TeamScope team(4);
        d::parallel_for(64, [&](size_t b) {
            sigset_t mask;
            pthread_sigmask(SIG_SETMASK, nullptr, &mask);
            blocked[b] = sigismember(&mask, SIGINT) == 1 && sigismember(&mask, SIGTERM) == 1 &&
                         sigismember(&mask, SIGALRM) == 1 && sigismember(&mask, SIGPROF) == 1;
            own[b] = std::this_thread::get_id() == here;
        });
    }
    for (size_t b = 0; b < 64; ++b)
        require(own[b] || blocked[b], "a worker's block ran with a host signal unblocked");
    pthread_sigmask(SIG_SETMASK, nullptr, &after);
    for (const int number : {SIGINT, SIGTERM, SIGALRM, SIGPROF})
        require(sigismember(&after, number) == sigismember(&before, number), "the calling thread's mask changed");
}
#endif

// A parallel_for at a split site counts one job when it runs on the team: not without a team, not for one block, and
// not inside a block (where it runs serially).
void split_jobs_counted() {
    require(d::split_jobs() == d::SplitJobs{}, "jobs counted without a team");
    d::parallel_for(d::SplitSite::TaylorRow, 8, [](size_t) {});
    require(d::split_jobs() == d::SplitJobs{}, "a job counted without a team");
    const d::TeamScope team(3);
    d::parallel_for(d::SplitSite::TaylorRow, 8, [](size_t) {});
    d::parallel_for(d::SplitSite::TaylorRow, 1, [](size_t) {});
    d::parallel_for(d::SplitSite::FixedTransform, 4,
                    [](size_t) { d::parallel_for(d::SplitSite::SecularPacket, 4, [](size_t) {}); });
    d::SplitJobs expected{};
    expected[static_cast<size_t>(d::SplitSite::TaylorRow)] = 1;
    expected[static_cast<size_t>(d::SplitSite::FixedTransform)] = 1;
    require(d::split_jobs() == expected, "the jobs counted by site");
}

void nested_scope_keeps_the_team() {
    const d::TeamScope outer(6);
    {
        const d::TeamScope inner(3);
        require(d::parallel_width() == 6, "a nested solve started its own team");
    }
    require(d::parallel_width() == 6, "the outer team lost");
}

void four_teams_at_once() {
    std::vector<std::thread> threads;
    std::vector<int> ok(4, 0);
    for (int t = 0; t < 4; ++t)
        threads.emplace_back([&, t] {
            const d::TeamScope team(3);
            long sum = 0;
            for (int round = 0; round < 200; ++round) {
                std::vector<long> part(32, 0);
                d::parallel_for(32, [&](size_t b) { part[b] = static_cast<long>(b) + round; });
                for (const long v : part) sum += v;
            }
            ok[t] = sum == 200L * (31 * 32 / 2) + 32L * (199 * 200 / 2);
        });
    for (auto &thread : threads) thread.join();
    for (const int v : ok) require(v == 1, "a team among four gave a wrong sum");
}

void many_short_jobs() {
    const d::TeamScope team(8);
    std::atomic<long> total{0};
    for (int round = 0; round < 5000; ++round)
        d::parallel_for(3 + round % 7, [&](size_t b) { total.fetch_add(static_cast<long>(b) + 1); });
    long expected = 0;
    for (int round = 0; round < 5000; ++round) {
        const long n = 3 + round % 7;
        expected += n * (n + 1) / 2;
    }
    require(total.load() == expected, "short jobs lost a block");
}
} // namespace

int main() {
    try {
        serial_without_a_team();
        for (const unsigned workers : {2u, 3u, 5u, 8u})
            for (const size_t blocks : {size_t{2}, size_t{7}, size_t{64}, size_t{1000}}) each_block_once(workers, blocks);
        nested_runs_serially();
        exceptions_come_home();
        floating_point_environment();
#if defined(__unix__) || defined(__APPLE__)
        workers_block_signals();
#endif
        split_jobs_counted();
        nested_scope_keeps_the_team();
        four_teams_at_once();
        many_short_jobs();
        serial_without_a_team(); // every scope has ended
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
    std::cout << "parallel tests passed\n";
    return 0;
}
