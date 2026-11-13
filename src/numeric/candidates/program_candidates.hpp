#pragma once
#include "numeric/evaluate/program_machine.hpp"
#include "core/query.hpp"
#include "core/engine_statistics.hpp"
#include <cstddef>
#include <cstdint>
#include <memory>
#include <span>
#include <vector>

namespace certpolysolve::detail {
// The binary64 stage of the program route's preliminary phase (program_advance, solver_program.cpp): Jacobi Aberth
// sweeps in machine arithmetic from the seeds, every point evaluated at the sweep's start positions in one
// batched call, every correction from those positions, then the update. Proposals only, never certificates.
//
// A point stops Converged when its relative step falls below 2^program_machine_small_step_log2 (MachineEngine's
// SmallStep rule), Floor when its value is at the evaluator's floor and its correction is not small (noise: the
// point cannot be refined further in this arithmetic) or its derivative is at the floor (a cluster or a critical
// point: the ratio is unreadable), and Held when its evaluation is not finite (out of range for this
// arithmetic). Stopped points stay in the others' sums, as fixed roots. A correction that is not finite, or a
// coincident pair, holds the point for the sweep. A point whose evaluation is not finite in binary64 (its
// values out of the range, as a recurrence's at a seed far from its roots) is evaluated again, that sweep, in
// binary64 with a coarse exponent (machine::BlockDouble, block_double.hpp; its own batched call for those points
// only), and
// its correction taken from that: a point pays the scaled arithmetic's price only while its values are out of
// range, and the stage is one. A point whose relative step has not fallen by a bit against
// its smallest so far for program_machine_patience sweeps is at the floor too (stagnation: a straggler whose
// evaluation's noise is above the floor estimate, as Laguerre's largest roots, where the recurrence amplifies
// earlier roundings; Aberth near a simple root halves the step every sweep). The phase's window rule is
// applied as well: over a window of (n + 2) / 2 sweeps the largest relative step must halve against the
// smallest window so far, or every point still active is at the floor of this arithmetic (a stall). The
// stage is finished when no point is active.
//
// The driver hands the stage's points to the ball phase: the converged ones settled (inactive in the ball
// packets, in the sums as fixed roots), the floor and held ones iterated in balls after a read of their
// derivative balls (a raise straight to the precision the root counts needed when they have no correct bits).
class ProgramMachineStage {
  public:
    using C = machine::Complex<double>;
    enum class Stop : uint8_t { Active, Converged, Floor, Held };
    struct Advance {
        uint64_t sweeps = 0, scaled_evaluations = 0; // evaluations taken in the scaled arithmetic, in all
        size_t active = 0, converged = 0, floor = 0, held = 0;
        bool finished = false, stalled = false, cancelled = false;
    };
    // From the seeds' midpoints; with `active` (n flags), only the flagged points are iterated and the others
    // stand as fixed roots (Converged) in the sums. valid() is false when the evaluator refuses (a constant
    // outside binary64's range), a seed is outside binary64's range (non-finite, or a nonzero coordinate
    // rounding to zero), the arrays of a chunk do not fit storage_limit, or the rounding mode is not to nearest; nothing else may be called then. With `zero_order`
    // v > 0 the stage iterates the roots of q = p / x^v (ProgramMachine's zero order): n is q's degree.
    ProgramMachineStage(const Program &program, acb_srcptr seeds, size_t n, size_t chunk, uint64_t storage_limit,
                        const uint8_t *active = nullptr, uint64_t zero_order = 0);
    bool valid() const { return valid_; }
    // Up to `sweeps` sweeps; returns early when finished or cancelled. The counts are the state after the
    // call; `sweeps` the sweeps taken by it.
    Advance advance(unsigned sweeps, const cps_options &options, certpolysolve::detail::CandidateStatistics &statistics);
    void export_points(acb_ptr out) const; // exact
    std::span<const Stop> stops() const { return stops_; }

  private:
    void sweep(certpolysolve::detail::CandidateStatistics &statistics);
    Advance state() const;
    ProgramMachine<double> evaluator_;
    std::unique_ptr<ProgramMachine<machine::BlockDouble>> scaled_; // the range tier, built on first need
    uint64_t scaled_evaluations_ = 0;
    size_t n_ = 0, chunk_ = 0;
    std::vector<C> points_, batch_, sums_;
    std::vector<Stop> stops_;
    std::vector<size_t> active_;
    std::vector<double> best_step_;      // each point's smallest relative step so far (log2)
    std::vector<uint32_t> stagnant_;     // sweeps since that fell by a bit
    std::vector<ProgramMachine<double>::Evaluation> evaluations_;
    uint64_t window_sweeps_ = 0;
    double previous_window_ = 0, current_window_ = 0, last_largest_ = 0;
    bool valid_ = false, finished_ = false, stalled_ = false;
};

} // namespace certpolysolve::detail
