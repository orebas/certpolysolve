#pragma once
// The engines' switches: CandidateEngine (the candidate stage's engines), TraceEvent (the trace hook's record), the
// method enums, RegenerationPolicy (the secular regeneration's switches), StudyStatistics (internal counters for
// tests and benchmark workers), NumericalSnapshot (owned diagnostic state), EnginePolicy (every switch of the global
// driver). The maintained recipes are built separately in solve/policies.hpp and solve/policies.cpp.
#include "core/engine_config.hpp"
#include "core/engine_statistics.hpp"
#include "core/number.hpp"
#include "core/heuristic_budget.hpp"
#include <cstdint>
#include <functional>
#include <memory>
#include <string>
#include <vector>

namespace certpolysolve::detail {
// The candidate stage's engines (EnginePolicy::candidates). Every recipe leaves the field at Adaptive, which
// factor_step resolves to Centered or Aberth; the methods pass Aberth and Secular to the packets themselves, and
// Aberth, Centered and Secular stay policy values that tests use as levers. The others are the research build's:
// its test engines (since 22 October 2026) and the Portfolio.
enum class CandidateEngine {
#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
    Legacy,     // Aberth packets whose next packet moves every centre (no inactive marks)
    DoubleSeed, // Legacy from a binary64 start (binary64_seeds)
#endif
    Aberth,
    Secular,
    Centered,
#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
    Hybrid, // binomials in closed form at every degree; else Aberth from degree 128, Centered or Legacy below
    // Five that send every factor of degree 2 or more to the structured route (test_structured_engine):
    Persistent, // with Jacobi sweeps on native secular input too
    Tree,       // as Persistent, with the interaction tree at every degree
    Seeded,     // no more (its seeds were the Mandelbrot recurrence's, deleted in W8's S2)
    Series,     // on the local Taylor cache (EvaluationMethod::TaylorCache)
    Direct,     // evaluating the source directly (PolynomialSourceOptions::candidate_direct)
#endif
    Adaptive,
#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
    DurandKerner, // Weierstrass steps in candidate_packet (also the Portfolio's)
    Portfolio, // the research build's engine portfolio (research/solvers/portfolio.cpp)
#endif
};
// Whether e is one of the research build's engines above; each is constant false in the default build, which has
// none of them, so the conditions that exclude them read the same in both builds.
#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
inline constexpr bool portfolio_engine(CandidateEngine e) { return e == CandidateEngine::Portfolio; }
inline constexpr bool legacy_centers_engine(CandidateEngine e) {
    return e == CandidateEngine::Legacy || e == CandidateEngine::DoubleSeed;
}
inline constexpr bool hybrid_engine(CandidateEngine e) { return e == CandidateEngine::Hybrid; }
// The engines that send every factor to the structured route (factor_step).
inline constexpr bool test_structured_engine(CandidateEngine e) {
    return e == CandidateEngine::Persistent || e == CandidateEngine::Tree || e == CandidateEngine::Seeded ||
           e == CandidateEngine::Series || e == CandidateEngine::Direct;
}
inline constexpr bool persistent_engine(CandidateEngine e) { return e == CandidateEngine::Persistent; }
inline constexpr bool tree_engine(CandidateEngine e) { return e == CandidateEngine::Tree; }
inline constexpr bool direct_engine(CandidateEngine e) { return e == CandidateEngine::Direct; }
inline constexpr bool durand_kerner_engine(CandidateEngine e) { return e == CandidateEngine::DurandKerner; }
#else
inline constexpr bool portfolio_engine(CandidateEngine) { return false; }
inline constexpr bool legacy_centers_engine(CandidateEngine) { return false; }
inline constexpr bool hybrid_engine(CandidateEngine) { return false; }
inline constexpr bool test_structured_engine(CandidateEngine) { return false; }
inline constexpr bool persistent_engine(CandidateEngine) { return false; }
inline constexpr bool tree_engine(CandidateEngine) { return false; }
inline constexpr bool direct_engine(CandidateEngine) { return false; }
inline constexpr bool durand_kerner_engine(CandidateEngine) { return false; }
#endif
struct ResearchSchedulingMetrics;
struct TraceEvent {
    const char *stage, *phase;
    uint64_t precision, active;
    const certpolysolve::detail::CandidateStatistics *candidates;
    const certpolysolve::detail::SearchStatistics *search;
    const char *engine = "", *tier = "", *reason = "";
    uint64_t preparation_bytes = 0, certified = 0;
    bool evaluation_reliable = false;
    uint64_t arithmetic_bits = 0, representation_revision = 0;
    uint64_t model_settled = 0, model_evaluations = 0, model_error_stops = 0, model_small_stops = 0;
    uint64_t algorithm_cutoff = 0;
    uint64_t exact_budget = 0, required_bits_lower_bound = 0;
    uint64_t exact_input_bits = 0, exact_output_bits = 0, preparation_ns = 0;
    bool preparation_accepted = false;
    uint64_t factor = 0, factor_power = 1;
    bool centered = false;
    double metered_work = 0; // diagnostic charge, never a timing-based decision
};
enum class RefinementMethod { Automatic, Aberth, Secular };
enum class StudyMethod { None, Aberth, Secular };
enum class SecularModelPolicy { Legacy, Continuation };
// The secular regeneration's switches, each independent of the others: off by default (D5), all five on in the
// maintained secular recipe (maintained_method_policy).
struct RegenerationPolicy {
    bool value_only = false;
    bool reuse_products = false;
    bool remember_precision = false;
    // While every part of every node is at most two limbs long (a binary64 number in every regeneration but the
    // last, in practice), a node's product of differences, and the factors of an unmoved node's update, are
    // formed in integer mantissas with a bound known in advance (node_products.hpp) where that bound serves the
    // residue; ball arithmetic otherwise.
    bool integer_node_products = false;
    // A dense source's value at a node whose parts are at most two limbs long (a binary64 node, or a last
    // regeneration's of 115 bits) is taken by Horner's rule on integer limbs with a short multiplier and a
    // disk error (limb_horner_evaluator.hpp), at every precision; the block dot products or the source's
    // own evaluation otherwise.
    bool limb_horner = false;
};
// Internal counters of the engines, for the tests and benchmark workers that pass EnginePolicy::study_statistics;
// never part of the installed statistics ABI.
struct StudyStatistics {
    uint64_t machine_packets = 0, scaled_packets = 0, adaptive_packets = 0;
    uint64_t preparations = 0, reuses = 0, continuations = 0, ingested_factors = 0;
    uint64_t same_precision_handoffs = 0, double_double_packets = 0;
    uint64_t scaled_recoveries = 0, secular_preparations = 0, secular_reuses = 0;
    uint64_t rectangular_retries = 0, rectangular_accepted = 0;
#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
    uint64_t exact_retries = 0, exact_accepted = 0; // exact_residues_at_budget's (research build)
#endif
    uint64_t source_values = 0, source_joints = 0, product_builds = 0, product_reuses = 0;
    // Of product_builds, the ones formed by ShortNodeProducts.
    uint64_t integer_node_products = 0;
    // Updates of an unmoved node's residue whose factors ShortNodeProducts formed.
    uint64_t integer_node_updates = 0;
    // Of the regeneration's evaluations, the ones by LimbHornerEvaluator.
    uint64_t limb_horner_evaluations = 0;
    uint64_t source_precision_retries = 0, product_precision_retries = 0;
    uint64_t source_peak_bits = 0, product_peak_bits = 0;
    uint64_t source_hint_starts = 0, deferred_validations = 0;
    uint64_t machine_cover_precision_retries = 0, machine_cover_peak_bits = 0;
    // The certificate attempt's Newton steps and its block-evaluated points.
    uint64_t machine_cover_newton_steps = 0, machine_cover_block_evaluations = 0;
    // Points the attempt certified at their Newton point from the evaluation at the point itself.
    uint64_t machine_cover_stepped_certificates = 0;
    // Its evaluations on the Newton ladder, below the certificate precision.
    uint64_t machine_cover_rung_evaluations = 0;
    // Packets the attempt made on points its first pass refused.
    uint64_t machine_cover_rounds = 0;
    // Roots the attempt certified by an inclusion disk and the closing count.
    uint64_t machine_cover_disks = 0;
    // Roots the attempt certified by its compensated tier (compensated_certificate.hpp).
    uint64_t machine_cover_compensated_certificates = 0;
    // The attempt's double-double repairs of points the compensated tier refused (EnginePolicy::double_double_repair).
    uint64_t machine_cover_double_double_repairs = 0;
    // Groups of refused points the attempt's cluster start certified, every member (EnginePolicy::cluster_start).
    uint64_t machine_cover_cluster_starts = 0;
    uint64_t candidate_restarts = 0, candidates_replaced = 0;
    uint64_t preliminary_packets = 0, preliminary_completions = 0, ordinary_polishes = 0;
    // Preliminary validations not run because ball Horner loses the whole precision at a candidate.
    uint64_t preliminary_validations_skipped = 0;
};
// Owned diagnostic state. Model/candidate coordinates are proposals, never
// certificates. The coordinate relation is x^power = y + offset for this
// exact factor; multiplicity belongs to the factor, not numerical clustering.
struct NumericalSnapshot {
    std::string stage;
    uint64_t factor = 0, power = 1, multiplicity = 1, revision = 0, source_precision = 0;
    uint64_t zero_count = 0;
    GaussianRational offset;
    std::vector<ComplexBall> points, nodes, residues;
    std::vector<slong> coordinate_precision;
    std::vector<uint8_t> inactive;
    std::vector<uint64_t> counts;
};
// Switches for the global driver (SolveState). Field order is historical:
// replay patches and archived campaign builds match this text, so fields are
// classified here rather than regrouped.
//
// Maintained recipes. production_policy() in solve/policies.hpp is D5, which public Auto uses
// for the queries the explicit recipes do not support:
//   study_machine, prepared_source_refinement, early_power_composition,
//   dovetail_subdivision (every recipe below inherits these).
// maintained_method_policy() in solve/policies.cpp adds, for both explicit recipes:
//   complex_centering, centering_refinement, scale_invariant_centering,
//   centred_origin_basis, pellet_real_shift; for explicit classical:
//   study = refinement = Aberth, classical_error_stopping,
//   defer_classical_validation, newton_point_refinement,
//   candidate_allow_underflow, coefficient_newton_bound;
// and for explicit secular:
//   study = refinement = Secular, model_stopping = Continuation,
//   regeneration.{value_only, reuse_products, remember_precision, integer_node_products,
//   limb_horner},
//   model_driven_refresh, secular_preliminary, newton_point_refinement,
//   candidate_allow_underflow, coefficient_newton_bound, newton_point_handover,
//   geometric_local_restart, mixed_packets, machine_keeps_unmoved, packet_lane_sums.
// The standalone_subdivision_policy() recipe disables backend and enables
//   prepared_subdivision for exact complex isolation and private refinement.
// The numerical recipes use backend (true), candidates = Adaptive and
//   dovetail_share_log2.
//
// Instrumentation, which observes but never changes results: trace, trace_data,
//   snapshot*, study_statistics.
//
// Test and private controls: backend = false without prepared_subdivision (the
// older Subdivision control), candidates = Aberth, Centered or Secular (levers)
// and test_frozen_secular.
// The research build (CERTPOLYSOLVE_BUILD_RESEARCH) also has exact_residues_at_budget
// (kept as research), the test engines (the other
// candidates values) and reads iteration, arithmetic, real_method, compensated,
// jit and candidates = Portfolio.
//
// Switches that no maintained recipe set, or that were always on, have been
// removed. History: docs/code-history-20261021.md, internal.hpp 5.
// A public method may use one only through its recipe, after qualification.
// A synchronous service boundary in a determined local Aberth packet. Counts
// are cumulative within that packet; work_units is filled by the owning
// SolveState, including its in-flight outer action. These are observations,
// never root certificates. The callback may advance a different, already-owned
// fallback, but may not advance/finalize/fork this numerical state.
struct ResearchLocalPacketProgress {
    uint64_t evaluation_terms = 0, root_interactions = 0, candidate_sweeps = 0;
    uint64_t precision = 0, work_units = 0;
};
enum class ResearchServiceAction { Continue, StopNumerical };
using ResearchLocalPacketService = ResearchServiceAction (*)(void *, const ResearchLocalPacketProgress &);
struct EnginePolicy {
    bool backend = true;
    // Exact standalone route: require prepared square-free factors, start from
    // their native Cauchy covers, and use guarded Newton with trial-local
    // precision. Off in every existing accuracy/approximation recipe.
    bool prepared_subdivision = false;
    // M4's bounded heuristic phase and exclusive complete service afterwards.
    // Off until the Auto integration gate. Overrides support small boundary tests.
    bool complete_fallback = false;
    uint64_t heuristic_packets = HeuristicBudget::default_packets;
    uint64_t heuristic_precision_multiplier = HeuristicBudget::default_precision_multiplier;
    // The certificate attempt's compensated tier (compensated_certificate.hpp) tries each point before the ball
    // certificate. On in every recipe, and for an attempt called without a policy; off only in the tests of the ball
    // certificate's own steps, ladder and evaluations, which must still reach them.
    bool compensated_certificates = true;
    // Proposal 1 (dev/docs/cluster-double-double-design-20261027.md): before any ball certificate, the attempt tries
    // the compensated tier alone on every point, repairs up to heuristics::double_double_refine_points refused points
    // in double-double (CompensatedCertificate::refine) and tries the tier on them again. Needs the compensated tier.
    bool double_double_repair = true;
    // The cluster start (research/d8-clusters/README.md): before the ball certificate's pass, a group of points the
    // tiers refuse whose Newton disks all meet is moved by an Aberth packet from a circle about the cluster's centre.
    // Needs the compensated tier.
    bool cluster_start = true;
    CandidateEngine candidates = CandidateEngine::Adaptive;
    void (*trace)(void *, const TraceEvent &) = nullptr;
    void *trace_data = nullptr;
    std::function<void(NumericalSnapshot)> snapshot;
    uint64_t snapshot_factor = 0, snapshot_power = 1;
    GaussianRational snapshot_offset;
    // D5's choice of update source for secular packets (prepared_source_policy,
    // secular_refinement.cpp): a sparse coefficient source is evaluated
    // directly instead of through a degree-sized secular sum, and an exact
    // direct source's updates move to the interpolating representation. On in
    // production_policy(); inert under the explicit recipes, which set their
    // study and refinement.
    bool prepared_source_refinement = false;
    // Secular's exact centering may also be admitted and accepted by heights
    // at the polynomial's own root scale (centering.hpp), invariant under
    // x -> 2^s x, and skips shifts that no route can accept.
    bool scale_invariant_centering = false;
    // Neither centering route accepts a translate that carries two or more roots out of the binary64 stage's
    // reach of the new origin when no certified group at the origin holds them (unseparated_origin_roots,
    // centering.hpp): the translated iteration would peel them a scale at a time in multiprecision where the
    // untranslated one finds them in binary64. On in the maintained secular recipe, the only one whose
    // centering this is (prepare_secular_center). History: docs/code-history-20261021.md, internal.hpp 7.
    bool centering_keeps_origin_scales = false;
    // A centred secular source evaluates the uncentred polynomial where its expansion has the smaller terms
    // (PolynomialSourceOptions::origin): a root group near the original origin keeps the accuracy that the
    // goals, relative to that origin, ask of it, at the precision its own spacing needs.
    bool centred_origin_basis = false;
    // A secular round whose active nodes iterate at different tiers (the coordinate rule near the original
    // origin, restarted members) runs a binary64 packet on the nodes within the initial stage and a
    // multiprecision packet on the others, each group holding the other's points, instead of sending every
    // node to the multiprecision engine because one needs it. On in the maintained secular recipe
    // (maintained_method_policy). History: docs/code-history-20261021.md, internal.hpp 8.
    bool mixed_packets = false;
    // A node a binary64 packet does not move keeps its entry point, not that point's binary64 rounding: in a mixed
    // round's binary64 sweeps, and in an ordinary binary64 packet while the goal is the default one (binary64's own).
    // Nodes placed by a local restart or a multiprecision round can hold far more than binary64's bits, and rounding
    // them widens their inclusion disks until the disks merge. At a tighter goal an ordinary packet still writes back
    // the rounding. That is scaffolding, not a design (no rule may rest on a lucky
    // outcome): its gain at large degree and high goals is rounding that lands on roots which are binary64 numbers
    // (clustered-769, 4.7 times at 3,322 bits); with the roots 2^-56 off the binary64 grid it costs 1.9 times.
    // That 80-input screen found no effect at degrees 16-64; this is not a degree limit. Private degree-3 packets
    // reach the rule; a stage-based replacement remains experimental (dev/docs/tabled.md).
    // On in the maintained secular recipe. History: docs/code-history-20261021.md, internal.hpp 9.
    bool machine_keeps_unmoved = false;
    // A binary64 secular packet takes its sums over the model's terms and over the other points with one
    // reciprocal a term in four lanes (PacketSums, study_candidates.cpp) where every part is in the lanes' range
    // and no distance is zero; the one-accumulator loop with Complex's division otherwise. The points differ by
    // rounding alone. On in the maintained secular recipe.
    bool packet_lane_sums = false;
    // At a fresh model's floor, the binary64 secular packet's sub-floor steps are taken node-relative (NodeOffsets,
    // study_candidates.cpp), at the default goal. The packet rounds the points to binary64, so none leaves it more
    // accurate than about 53 bits, and a solve whose certificate needs a little more (the default goal's) climbs to a
    // multiprecision packet with model rebuilds around it. When every node of a fresh model settles without its step,
    // those steps are taken once with each point as node + binary64 offset and written back in multiprecision: a node
    // 2^-53|b| off its root comes out about 2^-106|b| off, and the validation there usually certifies. Every other
    // packet is unchanged, and so is a higher goal's staged default solve (SecularRefinement::step: its
    // multiprecision packet's points start the higher goal's certificate attempt, which floor-stepped points made
    // dearer). History: docs/code-history-20261021.md, internal.hpp 10.
    bool machine_node_offsets = false;
    // Pellet cluster counts shift by the real |c| after rotating the
    // polynomial by c/|c|, so off-axis centres do not re-box the error at
    // every step of the Taylor shift (pellet_count in counting.cpp).
    bool pellet_real_shift = false;
    RefinementMethod refinement = RefinementMethod::Automatic;
    // The input's nonzero part (the exact zero roots set apart) written as q(x^k) before the exact import
    // (prepare_power_composition, in preprocess in solver.cpp). Every recipe sets it (production_policy); a default policy, as tests build one, imports the input
    // as given and reduces each factor's power after the factorization.
    bool early_power_composition = false;
#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
    // Private scheduling experiment: surface a legacy handoff as a stop,
    // rather than silently starting a third engine. Never set in Auto.
    bool research_no_subdivision = false;
#endif
    // Integrated scheduling. Enabled only by its qualified controller recipe.
    // The outer share is separate from the local-restart admission fraction.
    enum class DovetailStart { Ready, PrecisionCluster, Never };
    bool research_integrated_dovetail = false;
    DovetailStart research_dovetail_start = DovetailStart::PrecisionCluster;
    unsigned research_dovetail_share_log2 = 4;
    double research_dovetail_growth = 1;
    ResearchSchedulingMetrics *research_scheduling_metrics = nullptr;
    bool research_warm_subdivision = false;
    // Private accounting repair. Changes service admission, never enabled by
    // maintained recipes: a zero-centered Taylor shift is a coefficient copy.
    bool research_identity_shift_meter = false;
    // Operation-priced native accounting; private and off in maintained recipes.
    bool research_service_meter = false;
    // Private packet-service experiment. Account at each local packet's
    // actual precision in both its no-service control and serviced arm.
    bool research_local_packet_meter = false;
    // Installed by SolveState's private driver. Synchronous borrow only;
    // warm forks explicitly discard it, as do ordinary maintained recipes.
    ResearchLocalPacketService research_local_packet_service = nullptr;
    void *research_local_packet_data = nullptr;
#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
    // The research build's Portfolio engine's choices (research/solvers/portfolio.cpp): its iteration and
    // arithmetic, its real-only method, and the machine engine's compensated evaluation and JIT.
    IterationMethod iteration = IterationMethod::Automatic;
    ArithmeticTier arithmetic = ArithmeticTier::Automatic;
    RealMethod real_method = RealMethod::Automatic;
    bool compensated = false, jit = false;
#endif
    StudyMethod study = StudyMethod::None;
    bool study_machine = false;
    SecularModelPolicy model_stopping = SecularModelPolicy::Legacy;
    RegenerationPolicy regeneration;
    // Dovetail a structured (secular) factor with Subdivision on a metered
    // work share (dovetail.hpp): once the factor has a certified cover,
    // Subdivision gets 2^-dovetail_share_log2 of the work.
    // production_policy() switches it on. The share of 2^-4 is a design
    // choice: when secular finishes alone, s = 1/16 of the metered work
    // costs s / (1 - s) = 1/15 of secular's own work at most.
    // History: docs/code-history-20261021.md, internal.hpp 11.
    bool dovetail_subdivision = false;
    unsigned dovetail_share_log2 = 4;
    // Test-only: once secular has a certified cover, its packets still do
    // their work but never move a node, a nonterminating iteration that
    // exercises dovetailing's guarantee.
    bool test_frozen_secular = false;
#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
    // Research control, off. At the precision budget, a residue the plain
    // evaluation leaves without its demanded accuracy is computed from the
    // factor's exact value at the node (regenerate_residues;
    // structured_source.hpp, evaluate_exact), within the exact budgets
    // (exact_fits). On, a solve under a tight precision budget returns
    // answers it otherwise stops short of, at a cost the default exact
    // budgets do not bound: about degree^2 x node bits^2 bit operations per
    // node per regeneration. Kept as research (since 22 October 2026), in the
    // research build only (research/solvers/exact_residues.cpp).
    // History: docs/code-history-20261021.md, internal.hpp 12.
    bool exact_residues_at_budget = false;
#endif
    bool secular_preliminary = false;
    bool model_driven_refresh = false;
    // Classical's scheduling, on in the maintained classical recipe
    // (maintained_method_policy): a cancellation floor pauses proposals at
    // this precision; only original-input checks certify roots.
    bool classical_error_stopping = false;
    bool defer_classical_validation = false;
    // Independent original-source Newton certificates tried before a global
    // inclusion representation is built (from the PWPoly campaign); on in the
    // maintained classical and secular recipes.
    bool newton_point_refinement = false;
    bool candidate_allow_underflow = false;
    bool coefficient_newton_bound = false;
    // Variant C, on in the maintained secular recipe: when that attempt fails
    // after its machine run converged, its machine points seed secular's
    // ordinary preliminary phase instead of being discarded. They are hints
    // only.
    bool newton_point_handover = false;
    // Re-initialize a group of nodes that crawls onto a cluster from its
    // local model (SecularRefinement::geometric_restart): the disk's root
    // count certified exactly, the model at the roots' mean, instead of
    // waiting for a precision growth that a moving iteration never triggers
    // (the crawl: the nodes converge only linearly onto the cluster, over
    // many regenerations). On in the maintained secular recipe
    // (maintained_method_policy); a private control otherwise.
    // History: docs/code-history-20261021.md, internal.hpp 13.
    bool geometric_local_restart = false;
    // The maintained classical and secular recipes (Auto's routes to them
    // included) admit exact Gaussian centering by monic coefficient height;
    // D5 does not.
    bool complex_centering = false;
    // Maintained explicit methods refine an accepted centre while that halves
    // the geometric-mean root distance (refine_center in factor_state.cpp).
    bool centering_refinement = false;
    std::shared_ptr<StudyStatistics> study_statistics;
    bool explicit_source_policy() const {
        return refinement != RefinementMethod::Automatic;
    }
};
} // namespace certpolysolve::detail
