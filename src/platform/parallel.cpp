// W9's first cut (parallel.hpp): the worker team a solve call starts, and its one construct, fork and join over fixed
// blocks.
//
// A job is the body, its number of blocks, an atomic counter that hands the blocks out, an atomic count of blocks
// left and one exception slot a block. It lives on the solve's stack for the duration of run(). The team's mutex
// guards the job pointer, its generation, the count of workers inside the job and the stop flag; nobody holds it
// while a block runs. A worker that sees a new generation takes the job, counts itself in, claims blocks until none
// is left, and counts itself out under the lock, waking the solve's thread when the last one leaves. The solve's
// thread claims blocks too, then waits until no block is left and no worker is inside, and only then forgets the job:
// so no worker touches it after run() returns, and the lock's release and acquire order every block's writes before
// the solve's thread reads them.
//
// State: thread_local current team and in-block flag (parallel.hpp); the team's count of the jobs each split loop
// ran on it, kept on the solve's thread alone. Entry points: parallel.hpp's.
// Design notes: docs/multithreading-design-20261023.md ("The first cut: low-risk only", "The threads").
#include "platform/parallel.hpp"
#include <atomic>
#include <cfenv>
#include <condition_variable>
#include <cstdint>
#include <exception>
#include <flint/flint.h>
#include <mutex>
#include <new>
#include <system_error>
#include <thread>
#include <vector>
#if defined(__unix__) || defined(__APPLE__)
#include <pthread.h>
#include <signal.h>
#define CERTPOLYSOLVE_POSIX_SIGNAL_MASK 1
#endif

namespace certpolysolve::detail {
namespace {
// The solve call's team, on the solve's thread (TeamScope sets it); never set on a worker.
thread_local WorkerTeam *current = nullptr;
// Inside a block, on whichever thread runs it: a parallel_for there runs serially (no nesting).
thread_local bool in_block = false;

// While it lives, the calling thread blocks every signal; then its own mask comes back. A thread inherits the mask of
// the thread that creates it, so workers created inside one start with every signal blocked: a signal meant for the
// host (an interrupt, a profiler's timer, a runtime's own) is delivered to one of the host's threads, never to a
// library thread. A fault inside a worker still ends the process, as it would on the solve's thread. Without POSIX
// signal masks (Windows) it does nothing.
class SignalsBlocked {
#ifdef CERTPOLYSOLVE_POSIX_SIGNAL_MASK
    sigset_t previous_;
    bool blocked_ = false;

  public:
    SignalsBlocked() {
        sigset_t every;
        sigfillset(&every);
        blocked_ = pthread_sigmask(SIG_SETMASK, &every, &previous_) == 0;
    }
    ~SignalsBlocked() {
        if (blocked_) pthread_sigmask(SIG_SETMASK, &previous_, nullptr);
    }
#else
  public:
    SignalsBlocked() {}
#endif
};
} // namespace

class WorkerTeam {
    struct Job {
        const std::function<void(size_t)> *body = nullptr;
        size_t blocks = 0;
        std::atomic<size_t> next{0}, left{0};
        std::vector<std::exception_ptr> errors;
        std::vector<int> flags; // the floating-point exception flags a worker's block raised
    };
    std::mutex mutex_;
    std::condition_variable start_, done_;
    Job *job_ = nullptr;
    uint64_t generation_ = 0;
    unsigned inside_ = 0;
    bool stopping_ = false;
    std::vector<std::thread> threads_;
    SplitJobs jobs_{}; // the solve's thread's alone (parallel_for at a site)

    // Claims and runs blocks until none is left. On a worker, each block starts with its flags clear and leaves the
    // flags it raised in its slot; the solve's thread runs its blocks in its own environment, where they raise their
    // flags as a serial loop would.
    static void work(Job &job, bool worker) {
        for (;;) {
            const size_t block = job.next.fetch_add(1, std::memory_order_relaxed);
            if (block >= job.blocks) return;
            in_block = true;
            if (worker) std::feclearexcept(FE_ALL_EXCEPT);
            try {
                (*job.body)(block);
            } catch (...) {
                job.errors[block] = std::current_exception();
            }
            if (worker) job.flags[block] = std::fetestexcept(FE_ALL_EXCEPT);
            in_block = false;
            job.left.fetch_sub(1, std::memory_order_acq_rel);
        }
    }
    void loop() {
        // A worker runs in the default floating-point environment: round to nearest, every exception masked.
        std::fesetenv(FE_DFL_ENV);
        uint64_t seen = 0;
        std::unique_lock<std::mutex> lock(mutex_);
        for (;;) {
            start_.wait(lock, [&] { return stopping_ || (job_ && generation_ != seen); });
            if (stopping_) break;
            seen = generation_;
            Job *job = job_;
            ++inside_;
            lock.unlock();
            work(*job, true);
            lock.lock();
            if (--inside_ == 0) done_.notify_all();
        }
        lock.unlock();
        flint_cleanup(); // this thread's FLINT caches, as cps_cleanup_thread frees a caller's
    }

  public:
    explicit WorkerTeam(unsigned workers) {
        threads_.reserve(workers - 1);
        const SignalsBlocked blocked; // the workers start with every signal blocked
        try {
            for (unsigned i = 1; i < workers; ++i) threads_.emplace_back([this] { loop(); });
        } catch (...) {
            stop();
            throw;
        }
    }
    ~WorkerTeam() { stop(); }
    WorkerTeam(const WorkerTeam &) = delete;
    WorkerTeam &operator=(const WorkerTeam &) = delete;

    void stop() noexcept {
        {
            std::lock_guard<std::mutex> lock(mutex_);
            stopping_ = true;
        }
        start_.notify_all();
        for (auto &thread : threads_)
            if (thread.joinable()) thread.join();
        threads_.clear();
    }
    unsigned size() const { return static_cast<unsigned>(threads_.size()) + 1; }
    void count(SplitSite site) { ++jobs_[static_cast<size_t>(site)]; }
    const SplitJobs &jobs() const { return jobs_; }

    void run(size_t blocks, const std::function<void(size_t)> &body) {
        Job job;
        job.body = &body;
        job.blocks = blocks;
        job.left.store(blocks, std::memory_order_relaxed);
        job.errors.resize(blocks);
        job.flags.assign(blocks, 0);
        {
            std::lock_guard<std::mutex> lock(mutex_);
            job_ = &job;
            ++generation_;
        }
        start_.notify_all();
        work(job, false);
        {
            std::unique_lock<std::mutex> lock(mutex_);
            done_.wait(lock, [&] { return job.left.load(std::memory_order_acquire) == 0 && inside_ == 0; });
            job_ = nullptr;
        }
        // The workers' flags raised here, so that a window the solve's thread holds (feholdexcept, fetestexcept) sees
        // every flag the loop raised, as after a serial loop.
        int raised = 0;
        for (const int f : job.flags) raised |= f;
        if (raised) std::feraiseexcept(raised);
        for (auto &error : job.errors)
            if (error) std::rethrow_exception(error);
    }
};

TeamScope::TeamScope(unsigned workers) : previous_(current) {
    if (workers < 2 || current || in_block) return;
    try {
        team_ = std::make_unique<WorkerTeam>(workers);
    } catch (const std::system_error &) {
        return; // no threads to be had: the solve runs serially
    } catch (const std::bad_alloc &) {
        return; // no memory for the team's bookkeeping: the solve runs serially
    }
    current = team_.get();
}
TeamScope::~TeamScope() {
    if (!team_) return;
    current = previous_;
    team_.reset();
}

void parallel_for(size_t blocks, const std::function<void(size_t)> &body) {
    WorkerTeam *team = in_block ? nullptr : current;
    if (!team || blocks < 2) {
        for (size_t block = 0; block < blocks; ++block) body(block);
        return;
    }
    team->run(blocks, body);
}
void parallel_for(SplitSite site, size_t blocks, const std::function<void(size_t)> &body) {
    if (current && !in_block && blocks >= 2) current->count(site);
    parallel_for(blocks, body);
}

unsigned parallel_width() { return current && !in_block ? current->size() : 1; }
SplitJobs split_jobs() { return current ? current->jobs() : SplitJobs{}; }
} // namespace certpolysolve::detail
