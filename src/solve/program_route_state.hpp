#pragma once
// State and one-shot controller for program input. The driver supplies resolved
// options, policy, result and initial precision; cancellation is already latched.
// Member order retains reverse destruction of owned sources and borrowing passes.
// Created only for program input. Owning that input directly keeps the route
// and its input in one allocation.
#include "numeric/exact/program_preparation.hpp"
#include "certify/attempts/separated_roots.hpp"
#include "numeric/candidates/program_recurrence_seeds.hpp"
#include "numeric/evaluate/program_machine.hpp"
#include "numeric/candidates/program_candidates.hpp"
#include "numeric/candidates/program_seeds.hpp"
#include "numeric/evaluate/program_source.hpp"
#include "methods/secular/secular_refinement.hpp"
#include <cstdint>
#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace certpolysolve::detail {

class ProgramRouteState {
    class Controller;
    // Program input (docs/slp-input-spec-20261015.md): prepared in the advance's first work unit (the options'
    // budgets price the preparation), then solved by the controller on a ProgramSource.
    cps_program input;
    ProgramPreparation preparation;
    std::unique_ptr<ProgramSource> source;
    std::unique_ptr<SecularRefinement> refinement;
    // The program route's preliminary phase: the points between work units, and the windows of
    // sweeps whose progress decides when it ends.
    struct Preliminary {
        std::vector<ComplexBall> points;
        // The start's placement from root counts on circles (program_seeds), for the diagnostics.
        ProgramSeedReport seeds;
        bool done = false;
        // The phase's precision: the initial precision, doubled at a stall whose evaluations have no correct bits
        // in the controller, and the precision the refinement then starts at.
        slong precision = 0;
        uint64_t sweeps = 0, window_sweeps = 0, raises = 0;
        // log2 of the largest movement of a point over a packet, the most over the window before this one and
        // over this one so far.
        double previous_window = std::numeric_limits<double>::infinity();
        double current_window = -std::numeric_limits<double>::infinity();
        // The binary64 stage before the ball packets (program_machine.hpp): its state across work units, then its
        // counts for the diagnostics, and the points it settled (inactive in the ball packets).
        std::unique_ptr<ProgramMachineStage> machine;
        bool machine_done = false, machine_refused = false, machine_stalled = false, polish = false;
        uint64_t machine_sweeps = 0, scaled_evaluations = 0; // the stage's range tier, evaluations in all
        size_t machine_converged = 0, machine_floor = 0, machine_held = 0;
        std::vector<uint8_t> settled;
        RecurrenceSeedReport recurrence; // the structural seeds' levels and counts, when they placed the seeds
    } preliminary;
    // The second version's separated-roots pass (program_certificate.hpp, docs/slp-second-version-note-20261017.md),
    // between the preliminary phase and the refinement, with its own source: the refinement takes `source`
    // when it starts.
    struct Disks {
        std::unique_ptr<ProgramSource> source;
        // For a program with a root of order v at zero (the source is q = p / x^v): p itself, of degree n + v, which
        // the pass's cluster step counts on where a circle lies far from zero (SeparatedRootsPass::advance).
        std::unique_ptr<ProgramSource> undivided;
        std::unique_ptr<SeparatedRootsPass> pass;
    } disks;
    // A power output c q^k (program_structure.hpp; docs/slp-structural-seeding-note-20261018.md): the route
    // solves q, the cone of node `base`, and lifts its cover once, when it publishes it. History:
    // docs/code-history-20261021.md, solver.cpp 1.
    struct Lift {
        uint64_t base = 0, exponent = 1, degree = 0, nodes = 0; // the output's degree, q's node count
    } lift;
    // A root of order v at zero (program_zero_order): the route solves q = p / x^v (the deflated ProgramSource) at
    // degree n - v, and its finish publishes the point 0, count and multiplicity v, beside q's cover, once.
    uint64_t zero_order = 0;
    std::optional<Program> base;
    // A product of distinct factors (output_factors): the diagnostics of the product lift, or of its abandonment.
    std::string product_note;

  public:
    explicit ProgramRouteState(const cps_program &p) : input(p) {}
    // Validate the output and expose its upper degree bound before solving.
    int64_t degree_bound() const;
    // Called once by the one-shot driver. A call-scoped controller borrows the
    // arguments; the state retains no driver or result reference. Owned stages
    // keep their existing source and policy contracts.
    void solve(const cps_options &, const EnginePolicy &, cps_result &, uint64_t initial_precision);
};

} // namespace certpolysolve::detail
