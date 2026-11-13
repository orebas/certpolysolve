#include "solve/factor_pass.hpp"
#include <algorithm>
#include <array>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {
void require(bool condition, const char *message) {
    if (!condition) throw std::runtime_error(message);
}
}

int main() {
    using certpolysolve::detail::FactorPass;
    try {
        // All subsets and visit orders: a finished factor cannot erase an
        // earlier blockage, and the last factor cannot determine the label.
        constexpr std::array reasons{CPS_STOP_INPUT, CPS_STOP_PRECISION,
                                     CPS_STOP_CELL_BUDGET, CPS_STOP_ITERATION_BUDGET};
        size_t orders = 0;
        for (unsigned mask = 0; mask < (1u << reasons.size()); ++mask) {
            std::vector<cps_stop_reason> visits{CPS_STOP_NONE, CPS_STOP_NONE};
            auto expected = CPS_STOP_PRECISION;
            for (size_t i = 0; i < reasons.size(); ++i) {
                if (mask & (1u << i)) {
                    visits.push_back(reasons[i]);
                    expected = reasons[i];
                }
            }
            std::sort(visits.begin(), visits.end());
            do {
                FactorPass pass;
                for (size_t i = 0; i < visits.size(); ++i) {
                    const bool exhausted = visits[i] == CPS_STOP_NONE
                        ? pass.visit_finished(visits.size())
                        : pass.visit_blocked(visits[i], visits.size());
                    require(exhausted == (i + 1 == visits.size()),
                            "factor pass exhausted before or after one full pass");
                }
                require(pass.has_blockage() == (mask != 0),
                        "finished factors changed whether the pass was blocked");
                require(pass.stop_reason() == expected,
                        "factor order changed the stop reason");
                ++orders;
            } while (std::next_permutation(visits.begin(), visits.end()));
        }
        // Progress invalidates old pass-wide refusals. A later input-limited
        // pass must not inherit a cell/iteration stop from earlier work.
        FactorPass pass;
        require(!pass.visit_blocked(CPS_STOP_ITERATION_BUDGET, 3), "early iteration stop");
        require(!pass.visit_blocked(CPS_STOP_CELL_BUDGET, 3), "early cell stop");
        pass.restart();
        require(!pass.has_blockage(), "progress retained an old blockage");
        require(!pass.visit_finished(3), "progress retained an old visit count");
        require(!pass.visit_blocked(CPS_STOP_INPUT, 3), "restarted pass ended early");
        require(pass.visit_finished(3) && pass.stop_reason() == CPS_STOP_INPUT,
                "restarted pass inherited the old stop label");
        // A persistent factor wait is reported again after another factor
        // progresses, so resetting the pass cannot hide an iteration stop.
        pass.restart();
        require(!pass.visit_finished(2), "fresh pass ended early");
        require(pass.visit_blocked(CPS_STOP_ITERATION_BUDGET, 2) &&
                    pass.stop_reason() == CPS_STOP_ITERATION_BUDGET,
                "revisited factor wait lost its iteration stop");
        require(!FactorPass::blocks(CPS_STOP_NONE) && !FactorPass::blocks(CPS_STOP_CANCELLED) &&
                    !FactorPass::blocks(CPS_STOP_WORK_BUDGET),
                "a call-level stop was treated as a factor blockage");
        std::cout << "factor-pass contract: " << orders << " visit orders and restart checks passed\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "FAILED: " << error.what() << '\n';
        return 1;
    }
}
