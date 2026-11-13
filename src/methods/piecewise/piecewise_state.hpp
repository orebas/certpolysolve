// PiecewiseRootState::Impl, the piecewise method's state across advance calls (W8's S4: the class out of
// piecewise_roots.cpp, its member functions defined in the files of their parts). Private to them.
// Axis mode is delegated to PiecewiseAxisPass; the source finisher is in piecewise_finisher.cpp,
// the rounds in piecewise_rounds.cpp, the rest (the constructor, advance, configure, snapshot, the guards and the
// observers) in piecewise_roots.cpp; the tiny members two parts call stay defined here. What the class does:
// piecewise_roots.cpp's header.
#pragma once

#include "methods/piecewise/piecewise_roots.hpp"
#include "methods/piecewise/piecewise_proposals.hpp"
#include "methods/piecewise/axis/piecewise_axis.hpp"
#include "methods/piecewise/piecewise_ownership.hpp"
#include "methods/piecewise/root_index.hpp"
#include "core/wide_int.hpp"
#include <algorithm>
#include <optional>
#include <limits>
#include <set>

namespace certpolysolve::detail {
class PiecewiseFinisherSource;
class PiecewiseSourceCertificate;
class PiecewisePointIndex;

struct PiecewiseRootState::Impl {
    DensePolynomial input;
    PiecewiseRootConfig config;
    PiecewiseRootStatistics stats;
    cps_result result;
    cps_options settings;
    size_t length = 0, zeros = 0;
    slong degree = 0, model_bits = piecewise_initial_model_bits, precision = 128;
    uint64_t covered = 0, previous_models = 0, previous_rounding = 0;
    uint64_t previous_candidates = 0, previous_moves = 0, previous_covered = 0;
    uint64_t previous_improvements = 0;
    uint64_t previous_discovery = 0;
    bool proposal_deferred = false;
    slong progress_model_bits = 0, progress_precision = 0;
    std::vector<double> best_corrections;
    bool supplied = false, configured = false, refused = false, pool_ready = false;
    bool live_pool = false;
    bool discovery_snapshot_emitted = false;
    bool reflect = false;
    bool force_reseed = false;
    bool recovery_interrupted = false;
    // A discovery recovery put off for one round (PiecewiseRootConfig::seeds_before_recovery).
    bool recovery_owed = false;
    slong searched_model_bits = 0, searched_precision = 0;
    RootIndex index;
    GuardIndex guards, seed_guards;
    std::unique_ptr<PiecewiseOwnership> ownership;
    std::vector<size_t> candidate_owner, owner_candidate;
    // Scheduling metadata only. Publication still requires the original-source
    // proof, owner identity and final disjointness check.
    std::vector<bool> accurate_owner;
    std::set<size_t> unfinished_owners;
    std::vector<ComplexBall> carried, round_points;
    std::shared_ptr<const DensePolynomial> diagnostic_source;
    std::vector<std::pair<uint64_t, PiecewiseCandidateEvent>> candidate_history;
    std::vector<Root> retained; // Proofs survive tighter goals and fallback.
    // Constructed only when axis mode starts a round. Its private implementation
    // owns all axis charts/brackets and borrows this stable input and statistics.
    std::optional<PiecewiseAxisPass> axis;
    std::unique_ptr<PiecewisePointIndex> point_index;
    std::unique_ptr<PiecewiseModelCursor> models;
    std::unique_ptr<RootVector> pending;
    std::vector<size_t> owners;
    slong pending_count = 0, pending_cursor = 0;
    bool pending_polished = false;
    std::vector<bool> pending_stable;
    struct Finishing {
        // Candidate IDs and proof IDs are stable, but their liveness is checked
        // again before every step. No backend view outlives its owner.
        std::vector<std::pair<size_t, size_t>> owners;
        size_t position = 0;
        // Counted work (counted_work's unit): what the cohort may spend, what
        // it has spent, and the parts of that spent on forming the cohort and
        // on preparing the source.
        double spent = 0, budget = 0, preparation = 0, admission = 0;
        double last_correction = std::numeric_limits<double>::infinity();
    };
    std::unique_ptr<PiecewiseFinisherSource> finisher_source;
    // Deflated discovery: the nonzero part's source (the zero roots removed),
    // counted work at the round's start, and the rung it last ran at.
    std::unique_ptr<PiecewiseFinisherSource> deflation_source;
    double round_start_work = 0;
    slong deflated_model_bits = 0, deflated_precision = 0;
    // Pool indices of the points deflated discovery stopped at a round's resolution.
    std::vector<size_t> limited_seeds;
    std::optional<Finishing> finishing;
    slong finisher_precision = 0, failed_finisher_precision = 0;
    uint64_t failed_finisher_covered = 0, failed_finisher_improvements = 0;
    // The first certified root of a precision is the pilot: its counted work,
    // and the part of it that every further root will cost.
    double pilot_work = 0, marginal_pilot_work = 0;
    static constexpr size_t new_owner = std::numeric_limits<size_t>::max();

    Impl(const DensePolynomial &p, PiecewiseRootConfig c, const std::vector<ComplexBall> *points);
    ~Impl();
    // Sources and cursors retain references to this state's members; its address must stay stable.
    Impl(const Impl &) = delete;
    Impl &operator=(const Impl &) = delete;
    Impl(Impl &&) = delete;
    Impl &operator=(Impl &&) = delete;
    bool valid(const cps_options &o) const;
    bool finished() const { return config.axis_only ? (axis && axis->complete()) : (length && covered == length - 1); }
    // The ladder's last rung is the cap: a rung that would pass it is cut to
    // it when the present one is below it, so a 256-bit cap is tried at 256
    // bits, not refused at 196 (queue item 3.2). A rung at the cap that still
    // cannot finish is the precision stop.
    slong clamped_rung(slong next, const cps_options &o) const;
    // A complete cover at a stop that is not the goal (audit fix 3, queue item
    // 3.2): complete_cover_at_stop (certification.cpp) on the certified regions
    // and the input. Such a cover is honest, not a rescue: its coarse region
    // meets no goal, and the state keeps its certified regions for the next
    // call. A cancelled call returns the certified regions as they are: the
    // caller asked for no further work.
    std::vector<Root> cover_at_stop(const cps_options &o) const;
    bool discovered() const {
        if (config.axis_only) return axis && axis->complete();
        // A counted cluster may need splitting after a stricter goal. Preserve
        // the ordinary regional recovery until its individual owners exist.
        return ownership && ownership->singletons() && ownership->count() == static_cast<uint64_t>(degree);
    }
    size_t owner_of(size_t candidate) const {
        return candidate < candidate_owner.size() && ownership->active(candidate_owner[candidate]) ?
            candidate_owner[candidate] : PiecewiseOwnership::none;
    }
    void remember_ownership(const Root &proof, size_t candidate = new_owner, bool mirror = true);
    bool needed_candidate(size_t i);
    bool owner_finished(size_t id) const;
    void snapshot(bool terminal, const cps_options &o, bool force = false);
    void mark_finished_owner(const Root &proof);
    // The source finisher may spend what the suspended model batch would
    // cost: next_work, that batch's counted construction work. Everything it
    // is compared with is counted work too, so the decision does not depend
    // on the machine's load.
    void begin_finishing(double next_work, const cps_options &o);
    // PiecewiseRootConfig::deflated_discovery, after a round that owned no new
    // root or left only discovery to do. For a free point z, Aberth's correction
    //     z - N / (1 - N S),  N = q(z) / q'(z),  S = sum_j w_j / (z - r_j),
    // with r_j every owned root (w_j its proof's count) and every other free
    // point (weight one), is Newton's step on q / prod_j (z - r_j)^w_j, whose
    // zeros are the roots not yet owned; q is the nonzero part. The free
    // points start on Fujiwara's bound on the roots' moduli and iterate at the
    // round's precision until a step is below the model's resolution
    // (2^-model_bits of |z|: converged), or the enclosure of q(z) holds zero
    // first (the round's resolution: stopped, unseeded), or
    // initial_sweep_allowance sweeps; together they may spend what the stalled
    // round spent in counted work. Differences z - r_j are taken on two-double
    // coordinates, exact to about 2^-106 of |z|: the held roots of a cluster
    // stay apart. Converged points that no owner or final guard contains join
    // the carried pool (mirrored into the upper half-plane when charts there
    // stand for both); the next round polishes and certifies them as any seed.
    void seed_missing_roots(const cps_options &o);
    void observe_finisher(const Finishing &task, const char *phase,
                          double step = 0, bool certificate = false) const;
    // Only a source-created certificate crosses this admission boundary;
    // owner identity, disjointness and remaining degree are checked here.
    bool admit_source_certificate(PiecewiseSourceCertificate certificate, size_t owner);
    void advance_finisher(const cps_options &o);
    void configure(const cps_options &o);
    // The outer precision a call reports (cps_result_work_bits: the last outer precision, not a backend's peak).
    void report_precision() {
        stats.model_bits = static_cast<uint64_t>(model_bits);
        stats.arithmetic_bits = static_cast<uint64_t>(precision);
        result.work_bits = stats.arithmetic_bits;
        result.progress.max_outer_bits = std::max(result.progress.max_outer_bits, result.work_bits);
    }
    void start_round(const cps_options &o);
    enum class RoundStep { Continue, Stop };
    // Observe and close a completed round before scheduling another one.
    // Axis limits stop here; interrupted discovery is never marked searched.
    RoundStep finish_round(const cps_options &o);
    bool seed_round_discovery(const cps_options &o);
    void schedule_next_round(const cps_options &o, bool seeded);
    void observe_population(bool terminal);
    void observe_selection(const RootModel &model, size_t i);
    void prepare_model(const RootModel &model, const cps_options &o);
    // Whether a fresh proposal (chart coordinates, not yet polished) lies where process_candidate retires a new
    // proposal without looking further: inside a certified root's guard, or inside a region another candidate
    // owns with one root. Such a point is a neighbouring chart's copy of a root already held; it is counted as
    // process_candidate counts it and dropped before polishing. The iteration's point and its polished position
    // differ by the iteration's residual, far below a guard's width, and a point kept here meets the same two
    // tests again after polishing; either way it was only a proposal.
    bool held_copy(const RootModel &model, const acb_t local);
    bool certify_guard(Root &root, const RootModel &model, const acb_t point,
                       const cps_options &o, ChartFailure *failure = nullptr);
    // The jet of the candidate in work, shared by its certificates (certify_root_model's shared_jet): one
    // candidate, one model, one precision from process_candidate's entry to its return.
    CertificateJet candidate_jet;
    void process_candidate(const RootModel &model, const cps_options &o, slong i);
    // Exact real coefficients: a nonreal certificate proves its conjugate
    // too. Lower-half-plane models are skipped under that symmetry, so every
    // path that publishes a nonreal root must also publish and finish its
    // mirror; otherwise the mirror's owner can never be finished.
    void publish_reflection();
    void process_pending(const RootModel &model, const cps_options &o);
    // See PiecewiseRootConfig::predicted_refinement for the break-even.
    bool skip_doubled_rung(slong refused_bits, const cps_options &o) const;
    // The regime exit (see the config). With the selector's estimate (Auto's
    // route, or the caller's piecewise where the selector would choose it):
    // counted work against the estimate at its overhead plus the global
    // budget, in binary64 operations. The caller's piecewise without one:
    // against the budget in binary64 operations. Else (a direct caller of the
    // state; the axis attempt has no regime exit, configure()): against the
    // budget at the state's own work per construction unit
    // (construction_work / construction_units, both counts).
    // An axis-mode attempt has taken its rounds: its last allowed round has ended, and it starts no other.
    bool axis_spent() const {
        return config.axis_only && config.axis_rounds && stats.rounds >= config.axis_rounds && !models;
    }
    bool regime_exit_due() const;
    cps_result advance(const cps_options &o);
};

} // namespace certpolysolve::detail
