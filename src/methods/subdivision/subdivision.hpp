#pragma once
#include "certify/counting.hpp"
#include <functional>

namespace certpolysolve::detail {
// The enclosure callback always starts from the immutable exact/uncertain
// coefficients. Increasing precision never shrinks declared input uncertainty.
using Encloser = std::function<void(acb_poly_t, slong)>;
// Opt in only for an immutable exact square-free polynomial whose arithmetic
// enclosure converges as precision grows. Declared input width is not rounding.
// Legacy remains the numerical routes' default; the exact standalone recipe
// and prepared-factor research controls opt into PreparedExact.
enum class NewtonAccuracy { Legacy, PreparedExact };
// Opt-in for immutable exact square-free input. A trial keeps its
// geometry fixed while requesting precision, and always uses PreparedExact.
// Existing numerical solve paths construct the default Retained policy.
enum class SubdivisionPrecision { Retained, TrialLocal };
struct SubdivisionPrecisionState {
    uint64_t base_bits = 0, trial_bits = 0, newton_bits = 0, count = 0;
    bool splitting = false, blocked = false, finished = false;
    RealBall certificate_radius;
};
struct SubdivisionCompleteEntry {
    bool used_warm_cover = false;
    bool padded_width = false;
    bool exact_point = false;
    uint64_t normalization_bits = 0;
};
class Subdivision {
    struct Impl;
    std::unique_ptr<Impl> impl_;

  public:
    Subdivision(const std::vector<Root>& cover, uint64_t multiplicity, uint64_t precision,
                certpolysolve::detail::SearchStatistics&,
                SubdivisionPrecision = SubdivisionPrecision::Retained);
    ~Subdivision();
    // A complete-route entry for the SAME immutable exact square-free factor
    // as global_bound and warm_cover. global_bound must be an origin-centered
    // certified bound containing every root, with multiplicity included in its
    // count. A certified, disjoint, complete warm cover may supply a smaller
    // all-roots square; otherwise the global bound supplies it. The caller keeps
    // its fine certificates for publication: this entry retains coarse geometry,
    // not an arbitrary warm forest's search history.
    //
    // The source-derived seed_scale_bits must bound log2 of the global square's
    // half-width. It bounds normalization storage and sets the support half-width
    // floor 2^-seed_scale_bits. The initial arithmetic is
    // always 64 bits, TrialLocal, with reset Newton speed. Containment supplies
    // the initial all-root annulus; later certificates are the native ones.
    // Null means invalid global-bound/count/scale arguments, not an inconclusive
    // numerical count. The caller must also enforce its ordinary resource caps.
    static std::unique_ptr<Subdivision> complete_start(
        const Root &global_bound, const std::vector<Root> &warm_cover,
        uint64_t multiplicity, uint64_t seed_scale_bits,
        certpolysolve::detail::SearchStatistics &, SubdivisionCompleteEntry * = nullptr);
    // Only complete_start establishes this provenance. It survives ordinary
    // step/reopen continuation; the generic warm constructor does not establish it.
    bool complete_entry() const;
    SubdivisionPrecision precision_policy() const;
    // One controller action. The previous complete cover remains
    // owned by the caller until complete_cover returns a replacement.
    // `allowance`: the work (in work()'s units) the step may spend on Taylor
    // shifts and Graeffe steps; a disk test that would need more stops before
    // the operation (suspended(), pending_work()) and the next step resumes it.
    // Dovetailing passes its share's allowance; infinite never stops a test.
    // Source enclosure, geometry and Newton work are outside this allowance.
    cps_stop_reason step(const cps_options &, const Encloser &, uint64_t cells_available,
                        double allowance = std::numeric_limits<double>::infinity(),
                        NewtonAccuracy = NewtonAccuracy::Legacy);
    // The last step stopped before an operation of this cost (work()'s units);
    // idle: before computing anything (no Taylor shift, Graeffe or Newton step).
    bool suspended() const;
    bool idle() const;
    double pending_work() const;
    void research_identity_shift_meter(bool enabled);
    void research_service_meter(const EnginePolicy *policy);
    void reopen(const cps_options &);
    bool complete_cover(std::vector<Root> &out) const;
    uint64_t cells() const;
    uint64_t pending() const;
    // Binary64-equivalent operations spent in step() (dovetail.hpp).
    double work() const;
    uint64_t min_bits() const;
    uint64_t max_bits() const;
    void progress(std::vector<cps_component_progress> &out, uint64_t factor) const;
    // Owned internal snapshots for research transition checks, never C API state.
    void precision_state(std::vector<SubdivisionPrecisionState> &out) const;
};
} // namespace certpolysolve::detail
