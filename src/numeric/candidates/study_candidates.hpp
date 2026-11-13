#pragma once
#include "numeric/candidates/machine_candidates.hpp"

namespace certpolysolve::detail {
// One fixed iteration family; only arithmetic can change. The factor owner
// supplies transformed coordinates and invalidates preparation on source changes.
class StudyCandidates {
    EnginePolicy policy_;
    CandidateWorkspace adaptive_;
    std::unique_ptr<NumericalEngine> machine_;
    ComplexPolynomial polynomial_;
    slong precision_ = 0;
    bool machine_done_ = false, scaled_ = false;

  public:
    explicit StudyCandidates(const EnginePolicy &policy)
      : policy_(policy)
      , adaptive_(&policy) {}
    void coordinates_changed();
    CandidateProgress advance(acb_ptr,
                              const acb_poly_t,
                              acb_srcptr,
                              slong,
                              unsigned,
                              const cps_options &,
                              certpolysolve::detail::CandidateStatistics &,
                              std::span<const uint8_t>,
                              bool);
};

// Persistent converted secular representation. Node/residue revision is supplied
// by the owning refinement, so changing precision or coordinates cannot reuse
// stale coefficients. A failed attempt never changes the retained proof cover.
class SecularMachineWorkspace {
    struct Impl;
    std::unique_ptr<Impl> impl_;

  public:
    SecularMachineWorkspace();
    ~SecularMachineWorkspace();
    // The machine owns the iteration in its initial stage only: while no
    // active node's iteration precision exceeds this bound for the goal (the
    // goal's bits and the driver's 32 guard bits, at least 64).
    static slong stage_bits(uint64_t goal_bits);
    bool advance(CandidateProgress &,
                 acb_ptr,
                 acb_srcptr,
                 std::span<const ComplexBall>,
                 std::span<const ComplexBall>,
                 std::span<const uint8_t>,
                 std::span<const slong>,
                 uint64_t,
                 unsigned,
                 const cps_options &,
                 certpolysolve::detail::CandidateStatistics &,
                 const EnginePolicy &, ModelStoppingState *model = nullptr, bool keep_unmoved = false,
                 bool floor_steps = false);
};
} // namespace certpolysolve::detail
