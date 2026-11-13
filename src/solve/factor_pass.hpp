#pragma once

#include <certpolysolve/certpolysolve.h>
#include <cassert>
#include <cstddef>

namespace certpolysolve::detail {

// A run of finished or blocked factors since the last productive step. A full
// pass may stop the solve, or tighten overlapping covers if no factor is blocked.
// Blockages accumulate across that pass; the reported reason must not depend on
// factor order. This scheduling state owns no roots, budgets or numerical work.
class FactorPass {
    size_t stalled_ = 0;
    bool cells_ = false, precision_ = false, input_ = false, iterations_ = false;

public:
    static constexpr bool blocks(cps_stop_reason reason) noexcept {
        return reason == CPS_STOP_CELL_BUDGET || reason == CPS_STOP_PRECISION ||
               reason == CPS_STOP_INPUT || reason == CPS_STOP_ITERATION_BUDGET;
    }

    [[nodiscard]] bool visit_finished(size_t factors) noexcept {
        return ++stalled_ >= factors;
    }

    [[nodiscard]] bool visit_blocked(cps_stop_reason reason, size_t factors) noexcept {
        assert(blocks(reason));
        cells_ = cells_ || reason == CPS_STOP_CELL_BUDGET;
        precision_ = precision_ || reason == CPS_STOP_PRECISION;
        input_ = input_ || reason == CPS_STOP_INPUT;
        iterations_ = iterations_ || reason == CPS_STOP_ITERATION_BUDGET;
        return visit_finished(factors);
    }

    [[nodiscard]] bool has_blockage() const noexcept {
        return iterations_ || cells_ || precision_ || input_;
    }

    [[nodiscard]] cps_stop_reason stop_reason() const noexcept {
        // An unresolved pass without an explicit refusal retains the precision
        // label. An input refusal wins only when no resource blocks progress.
        return iterations_ ? CPS_STOP_ITERATION_BUDGET
               : cells_    ? CPS_STOP_CELL_BUDGET
               : precision_ || !input_ ? CPS_STOP_PRECISION
                                      : CPS_STOP_INPUT;
    }

    // A productive factor step, or a newly tightened separation goal, starts a
    // fresh pass. A factor's persistent iteration_wait is owned by FactorState;
    // visiting it again records its blockage in this new pass.
    void restart() noexcept {
        stalled_ = 0;
        cells_ = precision_ = input_ = iterations_ = false;
    }
};

} // namespace certpolysolve::detail
