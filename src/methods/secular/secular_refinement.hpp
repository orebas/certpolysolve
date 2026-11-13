#pragma once
#include "core/cover.hpp"
#include <map>
#include <optional>
#include "numeric/candidates/candidates.hpp"
#include "numeric/evaluate/monic_source.hpp"
#include "numeric/candidates/study_candidates.hpp"
#include "numeric/candidates/local_model_solver.hpp"
#include "certify/secular_disk_count.hpp"

namespace certpolysolve::detail {
class PolynomialSource;
struct PartialNewtonCover; // certify/attempts/newton_cover.hpp
// The source uses y; the requested roots use x = y + offset. Keep the exact
// offset owned and use the same outward-rounded map for stopping and export.
class RootCoordinates {
    GaussianRational offset_;
    // The offset as an exact ball, when it is one (centres are low-denominator dyadics): what
    // original_point and the accuracy targets near the original origin use.
    ComplexBall exact_;
    bool exact_available_ = false;

  public:
    RootCoordinates() = default;
    explicit RootCoordinates(GaussianRational offset);
    bool identity() const { return offset_.zero(); }
    const GaussianRational &offset() const { return offset_; }
    Root original(const Root &, slong precision) const;
    bool meets(const Root &, const cps_options &, slong precision) const;
    // The exact offset, or null for the identity or an offset that is not exact.
    const acb_struct *exact_offset() const { return exact_available_ && !identity() ? exact_.value() : nullptr; }
    // log2 of |z| / |z + offset| when positive: the bits a coordinate near the original origin needs beyond
    // those of its own magnitude for the goals, which are relative to that origin. Zero otherwise.
    double origin_excess_log2(const acb_t z) const;
};

// Exact dyadic interpolation nodes; each residue encloses -p(b_i)/Q'(b_i)
// for the ORIGINAL monic polynomial. Partial updates retain that identity.
// The cost of a regeneration in binary64-equivalent operations (dovetail.hpp):
// every source evaluation and node product at the precision it ran at. Nodes
// regenerate at their own precisions, often far below the largest one, so
// the largest precision is no measure of this work. source_terms is what the
// same evaluations added to statistics.evaluation_terms.
struct RegenerationWork {
    double operations = 0;
    uint64_t source_terms = 0;
    // Fixed-width integer node products do not run at the source's precision.
    // Keep their factor count separate from the existing ball/Horner proxy.
    double service_operations = 0;
    double limb_horner_work = 0;
    uint64_t short_product_factors = 0;
};
bool regenerate_residues(std::vector<ComplexBall> &out, MonicSource &source,
                         std::span<const ComplexBall> nodes, std::span<const ComplexBall> old_nodes,
                         std::span<const ComplexBall> old_residues, std::span<const slong> iteration_bits,
                         const cps_options &options, certpolysolve::detail::CandidateStatistics &statistics, uint64_t &maximum_precision,
                         const EnginePolicy *policy = nullptr, bool partial = true,
                         std::span<slong> source_precision_hints = {}, RegenerationWork *work = nullptr,
                         size_t *budget_limited = nullptr, std::vector<uint8_t> *input_limited = nullptr,
                         const acb_struct *origin_offset = nullptr, std::span<const uint8_t> estimated = {});
// `estimated` (sized to the nodes, or empty) marks the nodes whose hint is an estimate of the precision their
// residue needs and not a precision that served before: a first evaluation there that falls short by a
// measured number of bits takes those bits once, in whole limbs, before the ladder's doublings.
// budget_limited receives the number of residues accepted short of the
// accuracy regeneration demands: that happens only when their evaluation
// reached options.max_work_bits. input_limited (sized to the nodes) marks
// the residues of an inexact source accepted short of that accuracy because
// their error is the input's: a doubled evaluation precision that fails to
// halve a residue's error bound cannot have been the larger part of it (the
// rounding part falls by 2^-p), so no precision shrinks it further. Such a
// residue is the best the input allows; it is not counted as budget-limited.

// The nodes a model cannot steer. The step the model proposes at a node b_k is
//     a_k / D_k,    D_k = sum_{j != k} (a_j + a_k) / (b_k - b_j) - 1,
// and the residues a_j are enclosures. A node is unsteered when that step,
// evaluated on the enclosures, is nonzero and has no correct bit (or is
// unbounded, when D_k's enclosure contains zero): the model then does not say
// where b_k should go, however accurate a_k itself is. A cluster the precision
// budget cannot resolve is the typical case: the neighbours' unresolved
// residues dominate D_k. Entries of inactive nodes are 0.
// The bits a ball evaluation of a polynomial of this degree by Horner's rule loses, at the worst of these points.
// A ball product by a point of argument t re-boxes the error rectangle, whose half-perimeter grows by
// |cos t| + |sin t| against the value's own growth: none on the axes, half a bit on the diagonals, and over the
// degree's products degree log2(|cos t| + |sin t|).
double ball_horner_loss_bits(acb_srcptr points, slong count, slong degree);
std::vector<uint8_t> unsteered_nodes(std::span<const ComplexBall> nodes, std::span<const ComplexBall> residues,
                                     std::span<const uint8_t> inactive);

// Outward binary64 bounds of a rectangle, for deciding cheaply that two rectangles are apart. The cover's merge
// loop tests every pair, and acb_overlaps subtracts the midpoints at their full length, which with long nodes
// can be a large part of a solve, spent mostly on pairs whose rectangles are far apart. A pair these bounds do
// not separate still gets the exact test, so the answer to every test is unchanged: the bounds are rounded
// outward, and arf_get_d rounds in the stated direction beyond binary64's range as well.
// History: docs/code-history-20261021.md, secular_refinement.hpp 1.
struct OutwardBounds {
    double real_lo, real_hi, imag_lo, imag_hi;
    explicit OutwardBounds(const acb_t box) {
        Float lo, hi;
        arb_get_interval_arf(lo.value(), hi.value(), acb_realref(box), std::numeric_limits<double>::digits);
        real_lo = arf_get_d(lo.value(), ARF_RND_FLOOR);
        real_hi = arf_get_d(hi.value(), ARF_RND_CEIL);
        arb_get_interval_arf(lo.value(), hi.value(), acb_imagref(box), std::numeric_limits<double>::digits);
        imag_lo = arf_get_d(lo.value(), ARF_RND_FLOOR);
        imag_hi = arf_get_d(hi.value(), ARF_RND_CEIL);
    }
    bool apart(const OutwardBounds &other) const {
        return real_hi < other.real_lo || other.real_hi < real_lo || imag_hi < other.imag_lo ||
               other.imag_hi < imag_lo;
    }
};
// kept (indexed like the nodes, with kept_roots): nodes that hold a region
// certified by the certificate attempt (PartialNewtonCover). A component whose
// every member is kept is given their regions instead of its own disk, and
// its members stay where they were certified (inactive): each region holds at
// least one root, every region of the cover is disjoint from every other and
// the counts sum to the degree, so each holds exactly its count (the attempt's
// closing test). Secular's own disk D(b_i, n |w_i|) is n times wider than such
// a region and would not meet the goal the region meets. A kept node in a
// component whose uncertified member's own disk reaches it (a node far from
// its root) stays where it was certified (inactive) while the component's
// other members iterate: its region stays valid, the component is one region
// of the cover until it splits, and a whole-component restart, which would
// replace kept nodes too, does not take it (restart_members takes active
// nodes only). Where no such disk reaches it, the kept node iterates with the
// component as before (a cluster with one root kept needs its restart). Two
// nodes far from their roots otherwise merged 29 kept certificates of Hermite
// 64's reduced polynomial into their component and restarted all 31 at 234
// bits (28 October 2026). A kept node in a component whose regions meet another
// region (the cover without the kept regions) is the iteration's as before.
bool secular_inclusion_cover(std::vector<Root> &out, std::vector<uint8_t> &inactive,
                             std::span<const ComplexBall> nodes, std::span<const ComplexBall> residues,
                             bool real, uint64_t multiplicity, slong precision, const cps_options &options,
                             const RootCoordinates &coordinates = {}, std::span<const Root> kept_roots = {},
                             std::span<const uint8_t> kept = {});

// Restart hints for the k members of a certified component: a centre from
// p^(k-1), then the degree-k local Taylor model about it. Where the working
// precision determines that model, a bounded Aberth packet supplies its
// roots; where it does not, the hints go on the circles of the model's
// Newton polygon (scale_candidates). The hints are proposals only; the
// original source certifies. working may grow, but only on a measured
// deficiency of the local model: a model that is not finite or has no
// leading coefficient, and, with `resolution` bits given (the precision the
// members are iterated at), a model whose undetermined hints would lie
// farther from the centre than 2^-resolution of its modulus (the function's
// definition has the reasoning). Zero asks for no resolution. `unresolved`,
// where given, is set when the hints returned are not all roots the model
// determines: scales only, or some placed on a circle of their own spread.
// Borrowed only while component_seeds runs. The actual charge belongs to the
// owning refinement's existing work meter. The raw-term counter prevents its
// later whole-step meter from charging those same terms again.
struct ResearchLocalPacketAccounting {
    double &work;
    uint64_t &prepaid_terms;
    ResearchLocalPacketService service = nullptr;
    void *data = nullptr;
    const EnginePolicy *policy = nullptr;
};
bool component_seeds(acb_ptr out, PolynomialSource &source, const Root &component, slong k, slong &working,
                     const cps_options &options, certpolysolve::detail::CandidateStatistics &statistics, const EnginePolicy &policy,
                     const std::optional<ComplexBall> &center_hint, slong resolution = 0,
                     bool *unresolved = nullptr, slong counted_hint_resolution = 0
                     , ResearchLocalPacketAccounting *accounting = nullptr
                     );

// Candidate repair only. Inactive points remain exact and unchanged;
// regeneration independently rejects any collisions that remain.
void separate_interpolation_nodes(std::span<ComplexBall>, std::span<const uint8_t> inactive,
                                  std::span<const slong> iteration_bits);

class SecularRefinement {
    std::unique_ptr<MonicSource> source_;
    // source_ as a PolynomialSource when it is one, captured once (source_ never changes): the attempt,
    // preliminary, the local models and the restarts take only that source.
    PolynomialSource *polynomial_ = nullptr;
    RootCoordinates coordinates_;
    std::vector<ComplexBall> nodes_, residues_, proposals_;
    std::vector<slong> iteration_bits_;
    // The iteration bits and the precision budget nodes_ and residues_ were regenerated for, and how many of
    // those residues that budget cut short: what a request for the same model keeps (step()).
    std::vector<slong> model_bits_;
    uint64_t model_budget_ = 0;
    size_t budget_limited_ = 0;
    // Precision growths a node missed while certified (step()): a node a validation returns to the active set takes
    // them on its return, so it iterates at the precision it would have reached had it stayed active.
    std::vector<uint8_t> missed_growths_;
    // Precision suggestions only, owned by this immutable source and node
    // ordering. No value, error bound or certificate is reused through them.
    std::vector<slong> source_precision_hints_;
    // The nodes whose hint is a restart's estimate (local_component_restart), for the regeneration that follows.
    std::vector<uint8_t> restart_estimated_;
    // For the members of a restart component whose model had to be raised: the bits that model took beyond
    // k bits a bit of its hints' resolution, which its evaluation cancels (at least one; zero for a node no
    // raised model has had). The same component's next model starts with them, and any other's sizing takes
    // the most of them (local_component_restart). Kept for the solve.
    std::vector<slong> restart_need_;
    // The members of a restart component whose last model did not determine its hints (component_seeds):
    // restarted again, the component starts at its sizing under the accuracy goal (local_component_restart).
    std::vector<uint8_t> restart_unresolved_;
    // Inactivity refers to unchanged candidate indices in the retained cover.
    // A rejected proposal must not replace either its flags or restart scale.
    struct CertifiedState {
        // The certified cover in the refinement's coordinates; the restarts read its components too. (Until
        // 26 October 2026 a second copy, restart_cover, held the same regions: consider set both from one cover and
        // nothing changed either after.)
        std::vector<Root> cover;
        std::vector<uint8_t> inactive;
        bool available = false;
        explicit CertifiedState(size_t size) : inactive(size, 0) {}
        void consider(std::vector<Root>, std::vector<uint8_t>, slong precision);
    } certified_;
    uint64_t multiplicity_, maximum_precision_ = 0;
    double work_ = 0; // binary64-equivalent operations (dovetail.hpp)
    uint64_t research_prepaid_local_terms_ = 0;
    bool research_packet_started_ = false;
    // One certificate re-seed is tried at the precision limit before the
    // solve stops for precision (reset by a goal change).
    bool limit_reseeded_ = false;
    // The binary64 machine owns the iteration while no node's precision is
    // above the initial stage's (SecularMachineWorkspace::advance). Under a
    // precision budget at or below that, the growth cannot raise the
    // precision past it, so the machine is retired when its arithmetic is
    // exhausted at the budget, and the multiprecision engine iterates at the
    // budget's precision before the precision stop (which the growth takes
    // only once that engine, too, is exhausted there). The retirement is for
    // the refinement's life: binary64 arithmetic that reached its floor has
    // nothing to add at a later goal of the same solve.
    bool machine_retired_ = false, last_packet_machine_ = false;
    // Nodes whose last residue was accepted as input-limited (regenerate_
    // residues): the growth raises no precision for them, and a factor whose
    // every active node is one stops for the input, not for precision.
    std::vector<uint8_t> input_limited_;
    EnginePolicy policy_;
    LocalModelSolver local_model_solver_;
    SecularMachineWorkspace machine_;
    ModelStoppingState model_stopping_;
    // The multiprecision group's model state in a mixed round (EnginePolicy::mixed_packets): the binary64
    // packet prepares model_stopping_ for its own nodes, and one preparation resets the other's flags.
    ModelStoppingState model_stopping_high_;
    std::unique_ptr<StudyCandidates> direct_machine_;
    enum class GrowthStage { None, RaisePrecision, RestartCandidates };
    GrowthStage growth_stage_ = GrowthStage::None;
    // nullopt continues this step with the refreshed candidates. A returned
    // reason ends the step, including CPS_STOP_NONE for a productive yield.
    std::optional<cps_stop_reason> advance_growth(const cps_options &, SearchStatistics &, CandidateStatistics &);
    // A stack-local packet record shared with the whole-step meter. The caller
    // retains it across cancellation, budget exits and unwinding; it is not
    // persistent refinement state. The accounting order stays in the phase.
    struct PacketWork {
        std::optional<uint64_t> terms = std::nullopt;
        int64_t bits = 0;
        uint64_t interactions = 0;
        // Mixed packets charge binary64 terms/interactions at binary64's weight;
        // the remaining units use bits (the multiprecision group's maximum).
        uint64_t binary64_units = 0;
        double service_evaluation = 0, service_interaction = 0;
        uint64_t service_scaled_binary64 = 0;
    };
    cps_stop_reason advance_packet(const cps_options &, SearchStatistics &, CandidateStatistics &,
                                  bool interpolated_model, PacketWork &);
    bool need_regeneration_ = true, pending_grow_ = false;
    bool exhausted_before_refresh_ = false; // consumed at the next model comparison
    // Compare after separating interpolation nodes. Exhausted means the
    // unchanged model follows an accepted exhausted packet under that policy.
    enum class ModelMatch { Different, Unchanged, Exhausted };
    ModelMatch model_match(const cps_options &) const;
    // A completed refresh may finish this step (proof, growth or restart),
    // or leave a validated model ready for the packet in the same step.
    enum class ModelRefresh { PacketReady, StepComplete, Cancelled };
    ModelRefresh refresh_model(const cps_options &, SearchStatistics &, CandidateStatistics &, RegenerationWork &);
    void validate_model(const cps_options &, SearchStatistics &, CandidateStatistics &);
    bool preliminary_done_ = false;
    bool newton_attempted_ = false;
    // Certificates a refused attempt held (PartialNewtonCover), by node: the
    // node stays in its certified region and the validation's cover uses the
    // region (secular_inclusion_cover). Kept while they meet the goal.
    std::vector<Root> kept_roots_;
    std::vector<uint8_t> kept_;
    void keep_certificates(PartialNewtonCover &, const cps_options &, slong precision, certpolysolve::detail::CandidateStatistics &);
    // The certificate attempt refused a goal above its ladder's first rung
    // from initial candidates: the iteration meets the library's default goal
    // first (`staged_`), then the attempt is retried from those points at the
    // caller's goal. Once per refinement (`staged_once_`).
    bool staged_ = false, staged_once_ = false;
    // The bits of the goal the caller asked for at the last step: the default goal staged first is a part of
    // the way there, and what a restart prepares for is sized by this one (local_component_restart).
    uint64_t caller_bits_ = 0;
    unsigned packets_ = 0, rounds_ = 0;
    double previous_correction_ = 0;
    // log2 of the last packet's largest correction (the iteration's own
    // progress, for the dovetail schedule's records); +infinity before any.
    double last_correction_log2_ = std::numeric_limits<double>::infinity();
    uint64_t last_goal_ = 0;
    cps_goal last_query_ = CPS_ROOT_ACCURACY;
    Number last_absolute_, last_radius_;
    void preliminary(const cps_options &, certpolysolve::detail::CandidateStatistics &);
    // `attempted` reports whether the attempt ran (it does not on recurrence
    // or real-only queries, other goals, or a source it cannot read).
    bool newton_cover(const cps_options &, certpolysolve::detail::CandidateStatistics &, bool *attempted = nullptr);
    // Members of a restart component whose residues all lack correct bits
    // (residue_unresolved): the model has no information for them at this
    // precision, but the certificate still depends on their geometry.
    bool unresolved_members(const std::vector<size_t> &members) const;
    // Place such members on the circle that minimises their component's
    // certified radius (see the definition). Returns false if no radius is
    // available (a member with an exact residue, or a degenerate geometry).
    bool certificate_reseed(const std::vector<size_t> &members, slong precision);
    // Precision growth re-initializes each certified restart component from
    // these members (empty: it leaves the component alone).
    std::vector<size_t> restart_members(const Root &component, std::span<const ComplexBall> points) const;
    // Growth re-initializes a component of this many members from a local
    // model (component_seeds); otherwise it restarts them on a circle.
    bool local_restart(size_t members) const;
    // That re-initialization for one component at the round's iteration precision: None if the component takes
    // no local model or the model gave no hints, Cancelled if the call was cancelled before the hints were set.
    enum class LocalRestart { None, Initialized, Cancelled };
    LocalRestart local_component_restart(const Root &component, const std::vector<size_t> &members, slong prec,
                                         const cps_options &, certpolysolve::detail::CandidateStatistics &, certpolysolve::detail::SearchStatistics &);
    // A component of more members than a local model takes, restarted by the groups of its nodes that the
    // model's exact counts show about their own roots (the definition has the reasoning). True if such a group
    // was found; the component is then not re-seeded on a circle.
    bool regroup_restart(const std::vector<size_t> &members, slong prec, const cps_options &, certpolysolve::detail::CandidateStatistics &,
                         certpolysolve::detail::SearchStatistics &);
    // The geometric restart (EnginePolicy::geometric_local_restart): after a
    // validation, each group of nodes (node_groups) whose disk's root count K
    // is certified exactly (an isolating plateau of the model's counts) and
    // whose local model at the roots' mean is justified (K dominates it,
    // Pellet, on a circle well inside the other roots) is re-initialized from
    // that model's roots, with nodes released or taken to make K: only a
    // group that crawls (its spread shrinking by a steady factor between one
    // half and one, round after round), once per halving of its spread.
    // Returns true if any was.
    bool geometric_restart(const cps_options &, SearchStatistics &, CandidateStatistics &);
    // Groups of active nodes whose diameter is below their distance to every
    // other active node (single linkage over the midpoints), as candidates
    // for the local model; each is a sorted index list.
    std::vector<std::vector<size_t>> node_groups() const;
    // The local model of a group of nodes (see group_model): the disk's
    // certified root count, whether the truncation is justified and
    // determined, the precisions used, its centre (the roots' mean), and the
    // count roots nearest the centre in original coordinates.
    struct GroupModel {
        slong count = -1;
        bool justified = false, determined = false;
        // The roots' mean was found: `center` holds it, and `outer_log2` is
        // the distance from it to the other roots (log2), whatever the model.
        bool centred = false;
        double outer_log2 = 0;
        slong working = 0, packet_bits = 0, extra = 0;
        ComplexBall center;
        std::vector<ComplexBall> hints;
    };
    GroupModel group_model(const std::vector<size_t> &members, const cps_options &, SearchStatistics &, CandidateStatistics &);
    // Members (sorted indices) -> log2 of the radius at the last attempt.
    std::map<std::vector<size_t>, double> geometric_attempts_;
    // Members -> log2 of the spread at the previous restart check (every
    // group then) and its shrink in bits since the one before: a group is
    // attempted only while it crawls.
    std::map<std::vector<size_t>, std::pair<double, double>> group_spreads_;
    // Clusters whose local model was declined: K roots about their mean, the
    // other roots beyond 2^outer_log2. The decision depends on the roots
    // (the Taylor shift at their mean), not on the nodes, so a group that
    // finds the same count about a centroid inside that disk is not modelled
    // again. Cleared when the goal changes.
    struct DeclinedCluster {
        ComplexBall mean;
        double outer_log2 = 0;
        slong count = 0;
    };
    std::vector<DeclinedCluster> declined_clusters_;
    bool inside_declined(const acb_t point, slong count) const;
    // Metered work (work_'s unit) of the restart's declined attempts: an
    // attempt is made only while it is at most the share the policy gives
    // Subdivision (2^-dovetail_share_log2) of the refinement's work.
    double declined_work_ = 0;
    void snapshot(const char *stage) const;

  public:
    SecularRefinement(std::unique_ptr<MonicSource> source, acb_srcptr seeds, uint64_t multiplicity,
                      slong initial_precision, EnginePolicy policy, LocalModelSolver,
                      RootCoordinates coordinates = {});
    cps_stop_reason step(const cps_options &options, certpolysolve::detail::SearchStatistics &search, certpolysolve::detail::CandidateStatistics &statistics);
    bool cover(std::vector<Root> &out) const;
    uint64_t precision() const { return maximum_precision_; }
    // Binary64-equivalent operations spent in step() (dovetail.hpp).
    double work() const { return work_; }
    double last_correction_log2() const { return last_correction_log2_; }
    size_t active() const;
    unsigned rounds() const { return rounds_; }
};

// The options of the nested solve that places a group's local model's roots as hints (the geometric restart,
// secular_restarts.cpp): the approximation goal at `bits` by the secular method, within the caller's precision cap,
// storage, cancellation and deadline. Its exact budget is zero. The model is a rounded Taylor polynomial of the input
// about the group's mean, so a square-free proof or exact factors of it say nothing about the input, and its answer
// is a hint that the global validation certifies again. A zero budget refuses every exact step (exact_storage_limit():
// the import, Q(i), the modular tests, exact centring), as the compact-exponent rule refused most models before
// 26 October's allowance admitted their write-out.
// Refused when the calling state's work or sweep allowance is spent. Otherwise
// the child gets their remainder, with the existing default local work ceiling.
std::optional<cps_options> local_model_options(const cps_options &options, uint64_t bits,
                                             const SearchStatistics &, const CandidateStatistics &);
} // namespace certpolysolve::detail
