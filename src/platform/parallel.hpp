#pragma once
// W9's first cut (docs/multithreading-design-20261023.md, "The first cut: low-risk only"): fork and join over fixed
// blocks, on worker threads that a solve call starts and joins (TeamScope), reached by the loops through the calling
// thread's current team. The rules that make it low-risk:
// - The solve's thread runs blocks beside the workers and goes on when every block is done (parallel_for).
// - A block never waits for another, never takes a lock and never starts a parallel loop of its own: a parallel_for
//   inside a block runs serially. So a loop always finishes, even with no worker free.
// - A block writes only its own output; everything else is read-only while the loop runs. The join orders every
//   block's writes before the solve's thread reads them.
// - The only lock hands out the job and signals the join; nobody holds it while a block runs.
// - An exception in a block is kept and rethrown on the solve's thread after the join (the lowest block's).
// - A worker runs in the default floating-point environment (round to nearest, every exception masked), the solve's
//   thread in its own, which a solve call sets to the same. The flags a worker's blocks raise are raised on the
//   solve's thread after the join, so a window that thread holds sees every flag the loop raised, as after a
//   serial loop.
// - Nothing outlives the call: the scope's destructor stops and joins the workers, and each worker frees its FLINT
//   caches (flint_cleanup) before it ends.
// - A worker starts with every signal blocked, so that a signal meant for the host's threads never runs its handler
//   on a library thread.
// The worker count changes no answer: the callers fix their blocks by the data, never by the count.
//
// State: the calling thread's current team and whether it is inside a block (thread-local, set by TeamScope and by
// the team while it runs a block); the team's count of the jobs each split loop ran on it (for the tests). Entry
// points: TeamScope (the solve calls in c_api.cpp), parallel_for and parallel_width (the split loops), split_jobs
// (c_api.cpp's testing counter). Design notes: docs/multithreading-design-20261023.md.
#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>

namespace certpolysolve::detail {
class WorkerTeam;

// The split loops, for the count of the jobs each ran on a team (certpolysolve_testing.h,
// cps_testing_result_split_jobs): which of them a solve reached with workers, so that a test can show that its input
// exercises the split it means to. The values are the testing API's site numbers.
enum class SplitSite : unsigned {
    SecularPacket = 0,     // the multiprecision secular packet's Jacobi sweep (candidates.cpp)
    MachineEvaluation = 1, // the binary64 engine's evaluations of a sweep (machine_candidates.cpp)
    TaylorRow = 2,         // a Taylor row by residue class (piecewise_models.cpp)
    FixedTransform = 3,    // the fixed-point transform's levels, conversion and enclosure (fixed_transform.cpp)
    BlockedEvaluation = 4, // the finisher's blocked dot products (disk_blocked_polynomial_evaluator.hpp)
    Count
};
inline constexpr size_t split_sites = static_cast<size_t>(SplitSite::Count);
using SplitJobs = std::array<uint64_t, split_sites>;

// For one solve call: with more than one worker, and no team already current on this thread (a nested solve keeps
// the outer call's), starts workers - 1 threads, stopped and joined when the scope ends. Otherwise, or where the
// threads or their bookkeeping cannot be had (std::system_error, std::bad_alloc), nothing: the solve runs serially.
class TeamScope {
    std::unique_ptr<WorkerTeam> team_;
    WorkerTeam *previous_ = nullptr;

  public:
    explicit TeamScope(unsigned workers);
    ~TeamScope();
    TeamScope(const TeamScope &) = delete;
    TeamScope &operator=(const TeamScope &) = delete;
};

// Runs body(0), ..., body(blocks - 1), each once, on the current team and the calling thread, and returns when
// every block has run; serially and in order where there is no team, inside a block, or for fewer than two blocks.
// An exception from a block is rethrown here after the join. body must write only its own block's output.
void parallel_for(size_t blocks, const std::function<void(size_t)> &body);
// parallel_for at one of the split loops: a call that runs on the team counts as one job of its site.
void parallel_for(SplitSite site, size_t blocks, const std::function<void(size_t)> &body);

// The threads a parallel_for here would use: the team's size, or 1 where it would run serially.
unsigned parallel_width();

// The jobs each split loop has run on the calling thread's current team (all zero without a team).
SplitJobs split_jobs();
} // namespace certpolysolve::detail
