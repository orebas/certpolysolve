// Two machine-arithmetic packet engines, the file named for the first. StudyCandidates runs a factor's packets of one
// fixed iteration on a coefficient polynomial: binary64, scaled binary64, then multiprecision. SecularMachineWorkspace
// runs the secular refinement's binary64 packet on its model of nodes and residues. Neither certifies a root.
//
// With study_machine and StudyMethod::Aberth, StudyCandidates::advance seeds from the caller or initial_candidates and
// tries at most twice with prepare_machine_engine's Aberth engine, in binary64, then in scaled binary64 if preparation
// or arithmetic fails, each within cps_options::sweeps_left. An engine that stops returns its points marked converged,
// for validation first, and ends the machine stage for good (machine_done_), as a higher precision does; other packets
// run the multiprecision candidate_packet. coordinates_changed drops the engine, not machine_done_.
//
// SecularMachineWorkspace::advance runs while no active node's bits exceed stage_bits(goal) = max(64, goal + 32). Per
// revision SecularPoints<double>::prepare (else SecularPoints<machine::ScaledDouble>) converts the nodes b_j and
// residues w_j of p(z) = prod_j (z - b_j) (1 - sum_j w_j / (z - b_j)) and zeroes a residue below 2^-2b |b_j|, with
// b = binary64_bits. SecularPoints::advance makes Jacobi sweeps of the Ehrlich-Aberth iteration: ratio gives Newton's
// q = p/p' (at a node, its removable limit), the step is q / (1 - q A), A = sum_{j != i} 1/(z_i - z_j), or -1/A where q
// is infinite (NewtonPole), and a packet ends once ModelStoppingState settles every node (SecularModelPolicy::Legacy:
// at a largest relative step below 2^(8-b)). Under packet_lane_sums the binary64 sums run in fixed-order lanes
// (PacketSums). A zero divisor, a nonfinite step or FE_INVALID, FE_DIVBYZERO, FE_OVERFLOW or (binary64) FE_UNDERFLOW
// rejects the packet; that revision is not tried again. At a fresh model's floor in binary64, with floor_steps and at
// most the default goal, NodeOffsets::refine takes the steps binary64 could not, node-relative and exactly.
//
// State: StudyCandidates' members (study_candidates.hpp), among them machine_ (the prepared engine) and machine_done_;
// SecularMachineWorkspace::Impl, two SecularPoints and a NodeOffsets with the revisions they were prepared for.
// Entry points: SecularMachineWorkspace::advance from SecularRefinement::step (secular_refinement.cpp, machine_), and
// stage_bits from step, geometric_restart and prepare_secular_center (factor_state.cpp); StudyCandidates from
// isolate_complex (classical_isolation.cpp, CandidateState::study), SecularRefinement::preliminary (a local one) and
// step (direct_machine_), with coordinates_changed from factor_step (solver.cpp) and refine_center (factor_state.cpp).
// Not here: the coefficient engine (MachineEngine, machine_candidates.cpp), multiprecision packets (candidates.cpp),
// the step's and the sums' terms' operations (aberth_step.hpp), stopping rules (model_stopping.hpp), validation
// (secular_refinement.cpp, classical_isolation.cpp).
// Design notes: docs/campaign-d1-algorithms.md (the two engines), docs/mathematics.md (Aberth's correction),
// docs/secular-scaling-findings-20260928.md (the 2^-2b rule), docs/final-report-20261011/README.md (item 7).
#include "core/query.hpp"
#include "numeric/candidates/study_candidates.hpp"
#include "numeric/candidates/aberth_step.hpp"
#include "numeric/candidates/machine_adapter.hpp"
#include "core/telemetry.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <type_traits>
#include <numeric>
#include <cfenv>
#include <optional>
#include <tuple>

namespace certpolysolve::detail {
namespace {
void
event(const EnginePolicy &p,
      const certpolysolve::detail::CandidateStatistics &s,
      const char *stage,
      const char *tier,
      const char *reason,
      slong bits,
      size_t active,
      uint64_t arithmetic_bits = 0,
      uint64_t revision = 0) {
    if (!p.trace) return;
    TraceEvent e{ stage, "packet", static_cast<uint64_t>(bits), active, &s, nullptr };
    e.engine = p.study == StudyMethod::Secular ? "secular" : "aberth";
    e.tier = tier;
    e.reason = reason;
    e.arithmetic_bits = arithmetic_bits;
    e.representation_revision = revision;
    p.trace(p.trace_data, e);
}
}
void
StudyCandidates::coordinates_changed() {
    machine_.reset();
    scaled_ = false;
    // The owner has translated retained proposals. Discard prepared coordinates,
    // but keep the arithmetic stage and requested precision: a translation
    // must not restart binary64 after this factor has asked for more accuracy.
}
CandidateProgress
StudyCandidates::advance(acb_ptr out,
                         const acb_poly_t p,
                         acb_srcptr initial,
                         slong precision,
                         unsigned sweeps,
                         const cps_options &o,
                         certpolysolve::detail::CandidateStatistics &s,
                         std::span<const uint8_t> inactive,
                         bool tree) {
    const auto n = acb_poly_degree(p);
    auto *counts = policy_.study_statistics.get();
    const uint64_t machine_bits = machine::binary64_bits;
    // The call's iteration budget (cps_options::max_iterations) bounds every attempt: the sweeps binary64 made
    // before its arithmetic failed are counted sweeps (statistics.candidate_sweeps, which sweeps_left reads),
    // so the scaled arithmetic's attempt, and the adaptive packet below, get what the budget leaves of the
    // packet, and none once it is spent: an attempt given the packet's whole length after a failure would
    // overrun the budget. Without a budget nothing changes: sweeps_left is then unbounded.
    // History: docs/code-history-20261021.md, study_candidates.cpp 1.
    CandidateProgress spent; // a failed attempt's progress, returned when the budget allows no other
    const auto allowed = [&] { return static_cast<unsigned>(std::min<uint64_t>(sweeps, o.sweeps_left(s))); };
    if (!machine_done_ && precision_ > 0 && precision > precision_) {
        // The owning solver requested more accurate arithmetic. Rebuilding a
        // binary64 workspace at each higher precision cannot satisfy that
        // request, even when a few machine candidates remain active.
        machine_done_ = true;
        machine_.reset();
        if (counts) ++counts->continuations;
    }
    if (policy_.study_machine && policy_.study == StudyMethod::Aberth && !machine_done_ &&
        o.max_work_bits >= machine_bits) {
        RootVector seeds(n);
        if (initial)
            _acb_vec_set(seeds.data(), initial, n);
        else if (!initial_candidates(seeds.data(), p, precision))
            machine_done_ = true;
        if (precision_ != precision || !acb_poly_equal(polynomial_.value(), p)) {
            machine_.reset();
            precision_ = precision;
            acb_poly_set(polynomial_.value(), p);
        }
        for (unsigned attempt = 0; !machine_done_ && attempt < 2; ++attempt) {
            if (!allowed()) return spent;
            if (!machine_) {
                TraceScope trace(
                  &policy_, "study_preparation", static_cast<uint64_t>(precision), static_cast<uint64_t>(n), &s);
                machine_ = prepare_machine_engine(IterationMethod::Aberth,
                                                  scaled_ ? ArithmeticTier::ScaledBinary64 : ArithmeticTier::Binary64,
                                                  p,
                                                  seeds.view(),
                                                  false, false, false, o.temporary_bytes);
                if (counts) ++counts->preparations;
            } else if (counts)
                ++counts->reuses;
            if (!machine_) {
                if (!scaled_) {
                    scaled_ = true;
                    continue;
                }
                machine_done_ = true;
                break;
            }
            const auto step = machine_->advance(allowed(), o, s);
            if (counts) {
                ++counts->machine_packets;
                if (scaled_) ++counts->scaled_packets;
            }
            event(policy_,
                  s,
                  "study_candidates",
                  scaled_ ? "scaled-binary64" : "binary64",
                  step.arithmetic_failure ? "arithmetic-failure" :
                  step.arithmetic_limited ? "accuracy-floor" :
                  step.exhausted ? "iteration-exhausted" :
                  step.candidates.converged ? "small-step" : "progress",
                  precision,
                  step.active,
                  machine_bits);
            if (step.candidates.cancelled) return step.candidates;
            if (step.arithmetic_failure) {
                spent = step.candidates;
                machine_.reset();
                if (!scaled_) {
                    scaled_ = true;
                    continue;
                }
                machine_done_ = true;
                break;
            }
            machine_->export_candidates(out);
            auto progress = step.candidates;
            if (step.arithmetic_limited || step.exhausted || progress.converged) {
                // Return proposals for original-input validation first. If they
                // are insufficient, the next call continues at adaptive precision.
                machine_done_ = true;
                machine_.reset();
                if (counts) ++counts->continuations;
                progress.converged = true; // Frozen D2 control behavior.
            }
            return progress;
        }
    }
    if (counts) ++counts->adaptive_packets;
    const auto method = policy_.study == StudyMethod::Secular ? CandidateEngine::Secular : CandidateEngine::Aberth;
    if (policy_.trace)
        event(policy_, s, "study_candidates", "adaptive", "fixed-method", precision,
              inactive.size() == static_cast<size_t>(n)
                ? static_cast<size_t>(std::count(inactive.begin(), inactive.end(), uint8_t(0)))
                : static_cast<size_t>(n), static_cast<uint64_t>(precision));
    if (!allowed()) return spent;
    return candidate_packet(out, p, initial, precision, allowed(), method, o, s, inactive, tree, &adaptive_);
}

namespace {
// The binary64 packet's sums over the model's terms and over the other points.
//
// A sweep takes, for every point z, sum_j w_j / (z - b_j) with its square-weighted and unweighted companions, and
// the Aberth sum of 1 / (z - z_j): two complex inverses for every pair of a point and a node, n^2 pairs a sweep.
// Written with Complex's division each inverse is Smith's algorithm, three real divisions behind a branch on
// which part of the divisor is the larger, which the processor cannot predict; with the node products and the
// evaluations in integers these sums are a large part of a hard secular solve. An inverse is also
// conj(d) / |d|^2: one division, no branch. Taken over separate arrays of real and imaginary parts in
// packet_lanes independent partial sums, the processor overlaps the divisions of neighbouring terms and the
// compiler packs them. History: docs/code-history-20261021.md, study_candidates.cpp 2.
//
// The partial sums are added in one fixed order, ((0 + 1) + (2 + 3)), so the result does not depend on how the
// compiler packs the lanes. It differs from the one-accumulator sum by rounding alone, and a packet's points are
// proposals: the model and the certificates are computed from them elsewhere.
//
// A packet runs with the floating-point exceptions watched (SecularMachineWorkspace::advance rejects it on an
// invalid operation, a division by zero, an overflow or an underflow), so the lanes must raise none that the
// general loop does not. |d|^2 needs room that d does not, and a point on a node gives d = 0, the model's pole.
// So the lanes run only while every part of every node and point is zero or between 2^-packet_part_exponent
// and 2^packet_part_exponent in size: a difference of two such parts is zero or at least a unit in the last
// place of the lower bound, 2^-(packet_part_exponent + 52), and a squared distance is zero or within
// [2^-904, 2^803], where neither it nor its reciprocal leaves binary64's range. A zero distance is recorded and
// its term is divided by the least normal number instead (the term is then zero): the sum it is in is discarded
// and the caller's general loop runs, as it did before. The guard is a maximum and not a test: a compiler
// moves a division into the arm of a test and computes it for every term, zero divisor included, and the
// division-by-zero flag it raises would reject every packet with a point on a node.
// History: docs/code-history-20261021.md, study_candidates.cpp 3.
inline constexpr size_t packet_lanes = 4;
inline constexpr int packet_part_exponent = 400;
static_assert(2 * (packet_part_exponent + std::numeric_limits<double>::digits - 1) < -std::numeric_limits<double>::min_exponent &&
                  2 * (packet_part_exponent + 1) + 1 < std::numeric_limits<double>::max_exponent,
              "a squared distance of parts in the lanes' range is a normal binary64 number");
inline bool
packet_parts_in_range(const std::vector<double> &parts) {
    const double low = std::ldexp(1.0, -packet_part_exponent), high = std::ldexp(1.0, packet_part_exponent);
    bool every = true;
    for (const double x : parts) {
        const double size = std::abs(x);
        every = every && (size == 0 || (size >= low && size <= high));
    }
    return every;
}
struct PacketSums {
    // sum w / d, sum w / d^2, sum 1 / d, sum of the sizes |Re| + |Im| of the first's terms; and the least
    // |d|^2 met.
    double sum_re[packet_lanes] = {}, sum_im[packet_lanes] = {}, slope_re[packet_lanes] = {},
           slope_im[packet_lanes] = {}, pole_re[packet_lanes] = {}, pole_im[packet_lanes] = {},
           size[packet_lanes] = {}, least[packet_lanes];
    PacketSums() {
        for (size_t l = 0; l < packet_lanes; ++l) least[l] = std::numeric_limits<double>::infinity();
    }
    static double total(const double *lanes) { return (lanes[0] + lanes[1]) + (lanes[2] + lanes[3]); }
    // Whether no term had a zero distance.
    bool apart() const {
        bool every = true;
        for (size_t l = 0; l < packet_lanes; ++l) every = every && least[l] > 0;
        return every;
    }
};
static_assert(packet_lanes == 4, "PacketSums::total adds four lanes");
inline void
packet_term(PacketSums &s, size_t lane, double dr, double di, double wr, double wi) {
    const double square = dr * dr + di * di;
    s.least[lane] = square < s.least[lane] ? square : s.least[lane];
    const double reciprocal = 1.0 / std::fmax(square, std::numeric_limits<double>::min());
    const double ir = dr * reciprocal, ii = -(di * reciprocal);
    const double tr = wr * ir - wi * ii, ti = wr * ii + wi * ir;
    s.size[lane] += std::abs(tr) + std::abs(ti);
    s.sum_re[lane] += tr;
    s.sum_im[lane] += ti;
    s.slope_re[lane] += tr * ir - ti * ii;
    s.slope_im[lane] += tr * ii + ti * ir;
    s.pole_re[lane] += ir;
    s.pole_im[lane] += ii;
}
inline void
packet_pole(PacketSums &s, size_t lane, double dr, double di) {
    const double square = dr * dr + di * di;
    const double reciprocal = 1.0 / std::fmax(square, std::numeric_limits<double>::min());
    s.least[lane] = square < s.least[lane] ? square : s.least[lane];
    s.pole_re[lane] += dr * reciprocal;
    s.pole_im[lane] -= di * reciprocal;
}
inline void
packet_quotient(PacketSums &s, size_t lane, double dr, double di, double wr, double wi) {
    const double square = dr * dr + di * di;
    const double reciprocal = 1.0 / std::fmax(square, std::numeric_limits<double>::min());
    s.least[lane] = square < s.least[lane] ? square : s.least[lane];
    const double ir = dr * reciprocal, ii = -(di * reciprocal);
    s.sum_re[lane] += wr * ir - wi * ii;
    s.sum_im[lane] += wr * ii + wi * ir;
}
// The model's three sums at z over n terms. The accumulators are local to the loop (returned by value): held
// through a reference they could alias the arrays, and every term would be stored and reloaded.
PacketSums
packet_model_sums(double zr, double zi, const double *node_re, const double *node_im, const double *residue_re,
                  const double *residue_im, size_t n) {
    PacketSums s;
    size_t j = 0;
    for (; j + packet_lanes <= n; j += packet_lanes)
        for (size_t l = 0; l < packet_lanes; ++l)
            packet_term(s, l, zr - node_re[j + l], zi - node_im[j + l], residue_re[j + l], residue_im[j + l]);
    for (size_t l = 0; j < n; ++j, ++l) packet_term(s, l, zr - node_re[j], zi - node_im[j], residue_re[j], residue_im[j]);
    return s;
}
// sum (w_j + w) / (z - b_j) over the terms [begin, end): a node's own limit, whose terms carry its residue w.
PacketSums
packet_limit_sums(double zr, double zi, const double *node_re, const double *node_im, const double *residue_re,
                  const double *residue_im, size_t begin, size_t end, double wr, double wi) {
    PacketSums s;
    size_t j = begin;
    for (; j + packet_lanes <= end; j += packet_lanes)
        for (size_t l = 0; l < packet_lanes; ++l)
            packet_quotient(s, l, zr - node_re[j + l], zi - node_im[j + l], residue_re[j + l] + wr, residue_im[j + l] + wi);
    for (size_t l = 0; j < end; ++j, ++l)
        packet_quotient(s, l, zr - node_re[j], zi - node_im[j], residue_re[j] + wr, residue_im[j] + wi);
    return s;
}
// sum 1 / (z - z_j) over the points [begin, end).
PacketSums
packet_pole_sums(double zr, double zi, const double *point_re, const double *point_im, size_t begin, size_t end) {
    PacketSums s;
    size_t j = begin;
    for (; j + packet_lanes <= end; j += packet_lanes)
        for (size_t l = 0; l < packet_lanes; ++l) packet_pole(s, l, zr - point_re[j + l], zi - point_im[j + l]);
    for (size_t l = 0; j < end; ++j, ++l) packet_pole(s, l, zr - point_re[j], zi - point_im[j]);
    return s;
}

template<class R>
class SecularPoints {
    using C = machine::Complex<R>;
    // Whether this arithmetic is binary64 itself, whose sums run in lanes (PacketSums).
    static constexpr bool binary64 = std::is_same_v<R, double>;
    std::vector<C> nodes_, residues_, points_, next_;
    // The nodes, the residues and a sweep's points as separate arrays of parts, for the lanes (binary64 only),
    // and whether the nodes' parts are in the lanes' range.
    std::vector<double> node_re_, node_im_, residue_re_, residue_im_, point_re_, point_im_;
    bool node_lanes_ = false;
    static void split(std::vector<double> &re, std::vector<double> &im, const std::vector<C> &z) {
        re.resize(z.size());
        im.resize(z.size());
        if constexpr (binary64)
            for (size_t i = 0; i < z.size(); ++i) {
                re[i] = z[i].re;
                im[i] = z[i].im;
            }
    }
    // The correction from the model's sums at a point that is no node: q = value / (derivative + value poles),
    // value = 1 - sum.
    bool quotient(C &q, C z, C sum, C derivative, C poles, R absolute_sum, ModelStop *stop) const {
        const C value = C(R(1)) - sum;
        const C denominator = derivative + value * poles;
        if (denominator == C{}) return false;
        q = value / denominator;
        if (stop && machine::finite(q))
            *stop = model_stop(machine::logabs(value), machine::logabs(absolute_sum),
                               machine::logabs(q), machine::logabs(z), nodes_.size(), machine::binary64_bits);
        return machine::finite(q);
    }
    std::vector<C> entry_; // the points as the packet received them, rounded to this arithmetic
    std::vector<uint8_t> unresolved_; // see residue_unresolved
    // The Newton correction q = p/p' at z under the model. When it is infinite
    // (see NewtonPole) q is unset, *pole says so, and advance() takes the
    // Aberth limit; that needs model stopping, and without it (the D5 recipe)
    // this branch keeps its older rules.
    // `own` names the node a packet's point started on: where it still is, the scan for a node under the point
    // ends at once.
    bool ratio(C &q, C z, ModelStop *stop = nullptr, NewtonPole *pole = nullptr, bool lanes = false,
               size_t own = std::numeric_limits<size_t>::max()) const {
        // p(z) = product(z-b_j) * (1 - sum(a_j/(z-b_j))).
        // Compute p/p' without expanding p. At a pole b_k, cancellation
        // gives the removable limit in the first branch below.
        const bool on_own = own < nodes_.size() && z == nodes_[own];
        if constexpr (binary64) {
            if (lanes && !on_own) {
                // In lanes, unless the point is on a node.
                const PacketSums sums = packet_model_sums(z.re, z.im, node_re_.data(), node_im_.data(),
                                                          residue_re_.data(), residue_im_.data(), nodes_.size());
                if (sums.apart())
                    return quotient(q, z, C(PacketSums::total(sums.sum_re), PacketSums::total(sums.sum_im)),
                                    C(PacketSums::total(sums.slope_re), PacketSums::total(sums.slope_im)),
                                    C(PacketSums::total(sums.pole_re), PacketSums::total(sums.pole_im)),
                                    R(1) + PacketSums::total(sums.size), stop);
            }
        }
        C sum, derivative, poles;
        size_t k = on_own ? own : 0;
        while (k < nodes_.size() && !(z == nodes_[k])) ++k;
        if (k < nodes_.size()) {
            // The limit q = w_k / D, D = sum_{j != k}(w_j + w_k)/(b_k - b_j) - 1.
            bool summed = false;
            if constexpr (binary64) {
                if (lanes) {
                    // The terms before and after the node in lanes; no other node is under the point (prepare
                    // refuses equal nodes), and a sum that met one all the same is left to the loop below.
                    const PacketSums before = packet_limit_sums(z.re, z.im, node_re_.data(), node_im_.data(),
                                                                residue_re_.data(), residue_im_.data(), 0, k,
                                                                residues_[k].re, residues_[k].im);
                    const PacketSums after = packet_limit_sums(z.re, z.im, node_re_.data(), node_im_.data(),
                                                               residue_re_.data(), residue_im_.data(), k + 1,
                                                               nodes_.size(), residues_[k].re, residues_[k].im);
                    if (before.apart() && after.apart()) {
                        sum = C(PacketSums::total(before.sum_re) + PacketSums::total(after.sum_re),
                                PacketSums::total(before.sum_im) + PacketSums::total(after.sum_im));
                        summed = true;
                    }
                }
            }
            if (!summed)
                for (size_t j = 0; j < nodes_.size(); ++j)
                    if (j != k) {
                        const C difference = z - nodes_[j];
                        if (difference == C{}) return false;
                        sum += (residues_[j] + residues_[k]) / difference;
                    }
            const C denominator = sum - C(R(1));
            if (denominator == C{}) {
                if (stop && pole) {
                    // D's rounding error is needed only here: a second
                    // pass over the same terms sums their magnitudes.
                    R absolute_sum(1);
                    for (size_t j = 0; j < nodes_.size(); ++j)
                        if (j != k) {
                            const C term = (residues_[j] + residues_[k]) / (z - nodes_[j]);
                            absolute_sum = absolute_sum + machine::abs(term.re) + machine::abs(term.im);
                        }
                    pole->at = true;
                    pole->unresolved = unresolved_[k] != 0;
                    pole->margin = machine::logabs(residues_[k]) - machine::logabs(absolute_sum) -
                                   evaluation_roundoff_log2(nodes_.size(), machine::binary64_bits);
                    return true;
                }
                // Without model stopping no node settles alone. A residue
                // below 2^-2b |z| (b = binary64_bits) is taken to mean the
                // node already is a root in this arithmetic: no step. (A rule
                // of thumb: it assumes the true D is at least the rounding
                // level.) Otherwise the packet is rejected, as it always was.
                if (!(machine::logabs(residues_[k]) <= machine::logabs(z) - 2 * machine::binary64_bits)) return false;
                q = C{};
            } else
                q = residues_[k] / denominator;
            if (stop && machine::finite(q))
                *stop = unresolved_[k] && q != C{}
                    ? ModelStop::EvaluationFloor
                    : model_stop(NAN, NAN, machine::logabs(q), machine::logabs(z), nodes_.size(), machine::binary64_bits);
            return machine::finite(q);
        }
        R absolute_sum(1);
        for (size_t j = 0; j < nodes_.size(); ++j) {
            const C inverse = C(R(1)) / (z - nodes_[j]);
            const C term = residues_[j] * inverse;
            if (stop) absolute_sum = absolute_sum + machine::abs(term.re) + machine::abs(term.im);
            sum += term;
            derivative += term * inverse;
            poles += inverse;
        }
        return quotient(q, z, sum, derivative, poles, absolute_sum, stop);
    }

  public:
    bool prepare(std::span<const ComplexBall> nodes, std::span<const ComplexBall> residues,
                 const ModelStoppingState *model, uint64_t revision) {
        nodes_.resize(nodes.size());
        residues_.resize(nodes.size());
        points_.resize(nodes.size());
        next_.resize(nodes.size());
        unresolved_.resize(nodes.size());
        for (size_t i = 0; i < nodes.size(); ++i) {
            if (!from_ball(nodes_[i], nodes[i].value())) return false;
            unresolved_[i] = residue_unresolved(residues[i].value()) || (model && model->unsteered(revision, i));
            // A residue whose enclosure lies below 2^-2b of its node (b =
            // binary64_bits) is zero in this arithmetic, resolved or not, in
            // either representation: the node is a root here, and its term is
            // below the rounding of the model's constant at every point
            // binary64 can tell from the node (ratio() takes such a node as a
            // root too). Kept, its products leave binary64's range: a step
            // underflows, and every packet of the round goes to the
            // multiprecision engine.
            // History: docs/code-history-20261021.md, study_candidates.cpp 4.
            Magnitude bound;
            acb_get_mag(bound.value(), residues[i].value());
            if (mag_get_d_log2_approx(bound.value()) < machine::logabs(nodes_[i]) - 2 * machine::binary64_bits)
                residues_[i] = C{};
            else if (!from_ball(residues_[i], residues[i].value()))
                return false;
        }
        // Two equal nodes refuse the model. Found by sorting (O(n log n)):
        // comparing every pair would cost each preparation O(n^2),
        // 3.4 x 10^10 comparisons at degree 262,144.
        // History: docs/code-history-20261021.md, study_candidates.cpp 5.
        std::vector<size_t> order(nodes.size());
        std::iota(order.begin(), order.end(), size_t{0});
        const auto key = [&](size_t i) { return std::pair{nodes_[i].re, nodes_[i].im}; };
        std::sort(order.begin(), order.end(), [&](size_t a, size_t b) { return key(a) < key(b); });
        for (size_t k = 1; k < order.size(); ++k)
            if (nodes_[order[k]] == nodes_[order[k - 1]]) return false;
        split(node_re_, node_im_, nodes_);
        split(residue_re_, residue_im_, residues_);
        node_lanes_ = binary64 && packet_parts_in_range(node_re_) && packet_parts_in_range(node_im_);
        return true;
    }
    bool advance(CandidateProgress &progress,
                 acb_ptr out,
                 acb_srcptr initial,
                 std::span<const uint8_t> inactive,
                 unsigned sweeps,
                 const cps_options &o,
                 certpolysolve::detail::CandidateStatistics &s,
                 ModelStoppingState *model = nullptr, bool keep_unmoved = false, bool lane_sums = false) {
        const auto n = nodes_.size();
        for (size_t i = 0; i < n; ++i)
            if (!from_ball(points_[i], initial + static_cast<slong>(i)))
                return false;
        entry_ = points_;
        for (unsigned sweep = 0; sweep < sweeps; ++sweep) {
            double largest = -std::numeric_limits<double>::infinity();
            next_ = points_;
            split(point_re_, point_im_, points_);
            // This sweep's sums run in lanes where the recipe asks and every part is in their range.
            const bool lanes = lane_sums && node_lanes_ && packet_parts_in_range(point_re_) &&
                               packet_parts_in_range(point_im_);
            for (size_t i = 0; i < n; ++i) {
                if (o.cancelled()) {
                    progress.cancelled = true;
                    return true;
                }
                if (inactive[i] || (model && model->settled(i))) continue;
                C q;
                ModelStop reason = ModelStop::Active;
                NewtonPole pole;
                {
                    StageTimer timer(s.evaluation_ns);
                    ++s.evaluations;
                    s.evaluation_terms += n;
                    if (model) model->evaluated();
                    if (!ratio(q, points_[i], model ? &reason : nullptr, &pole, lanes, i)) return false;
                }
                if (model && reason != ModelStop::Active) {
                    model->stop(i, reason);
                    continue;
                }
                C sum;
                {
                    StageTimer timer(s.interaction_ns);
                    bool summed = false;
                    if constexpr (binary64) {
                        if (lanes) {
                            // The other points in lanes, before and after this one, unless two coincide.
                            const PacketSums before = packet_pole_sums(points_[i].re, points_[i].im, point_re_.data(),
                                                                       point_im_.data(), 0, i);
                            const PacketSums after = packet_pole_sums(points_[i].re, points_[i].im, point_re_.data(),
                                                                      point_im_.data(), i + 1, n);
                            if (before.apart() && after.apart()) {
                                sum = C(PacketSums::total(before.pole_re) + PacketSums::total(after.pole_re),
                                        PacketSums::total(before.pole_im) + PacketSums::total(after.pole_im));
                                s.root_interactions += n - 1;
                                summed = true;
                            }
                        }
                    }
                    if (!summed)
                        for (size_t j = 0; j < n; ++j)
                            if (j != i) {
                                const C difference = points_[i] - points_[j];
                                if (difference == C{}) return false;
                                sum += machine::interaction_term(difference);
                                ++s.root_interactions;
                            }
                }
                C delta;
                if (pole.at) {
                    // Newton's correction is infinite here; the Aberth step is
                    // its limit -1/A when that is determined (NewtonPole), and
                    // otherwise this node's step is not: it settles at the
                    // evaluation floor while the other roots' packet stands.
                    if (!pole.resolved(machine::logabs(sum))) {
                        model->stop(i, ModelStop::EvaluationFloor);
                        continue;
                    }
                    if (sum == C{}) return false;
                    delta = machine::aberth_pole_limit(sum);
                } else {
                    const C denominator = machine::aberth_denominator(q, sum);
                    if (denominator == C{}) return false;
                    delta = q / denominator;
                }
                next_[i] = points_[i] - delta;
                if (!machine::finite(delta) || !machine::finite(next_[i])) return false;
                const double radius = machine::logabs(points_[i]);
                largest = std::max(largest, machine::logabs(delta) - (std::isfinite(radius) ? radius : 0));
            }
            points_.swap(next_);
            ++s.candidate_sweeps;
            progress.correction_log2 = largest;
            progress.model_exhausted = model && model->exhausted(inactive);
            progress.model_immediate = model && model->immediate();
            progress.converged = model ? progress.model_exhausted : largest < 8 - machine::binary64_bits;
            if (progress.converged) break;
        }
        // With keep_unmoved, a node the packet did not move keeps its entry point, not that point's rounding to
        // this arithmetic (EnginePolicy::machine_keeps_unmoved, in a mixed round).
        for (size_t i = 0; i < n; ++i) {
            if (inactive[i] || (keep_unmoved && points_[i] == entry_[i]))
                acb_set(out + static_cast<slong>(i), initial + static_cast<slong>(i));
            else
                to_ball(out + static_cast<slong>(i), points_[i]);
        }
        progress.finite = true;
        return true;
    }
};

// The steps a binary64 secular packet cannot take (EnginePolicy::machine_node_offsets). SecularPoints rounds the
// points to binary64, so a node about 2^-53|b| off its root settles without its step (the packet's floor); when every
// node of a fresh model settles so (ModelStoppingState::immediate), the solve goes on to validation and up the
// precision ladder. NodeOffsets then takes those steps once, node-relative: point i is node b_i plus a binary64
// offset, the node differences b_i - b_j come from the nodes' double-double values, (hi_i - hi_j) + (lo_i - lo_j),
// with binary64's relative accuracy unless two nodes agree to about machine::double_double_bits, and the point goes
// back as node + offset in multiprecision. A node 2^-53|b| off its root comes out about 2^-106|b| off, and the
// validation that follows meets nodes 53 bits closer. Only steps below binary64's floor for the point and above
// the stage's (2^(2 - stage)|z|; a closer node already holds the accuracy the stage works to) are taken. Every
// other packet stays SecularPoints' own: taking these steps in every packet, or writing every moved point back as
// node + offset, would make points of two limbs that buy nothing until the floor and make each model rebuild
// dearer. History: docs/code-history-20261021.md, study_candidates.cpp 6.
class NodeOffsets {
    using C = machine::Complex<double>;
    std::vector<C> high_, low_, residues_, offsets_, steps_;
    std::vector<uint8_t> unresolved_, refined_; // see residue_unresolved; a step taken here

    C node_difference(size_t i, size_t j) const { return (high_[i] - high_[j]) + (low_[i] - low_[j]); }

  public:
    bool prepare(std::span<const ComplexBall> nodes, std::span<const ComplexBall> residues,
                 const ModelStoppingState *model, uint64_t revision) {
        const size_t n = nodes.size();
        high_.resize(n);
        low_.resize(n);
        residues_.resize(n);
        unresolved_.resize(n);
        for (size_t i = 0; i < n; ++i) {
            machine::Complex<machine::DoubleDouble> node;
            if (!from_ball(node, nodes[i].value())) return false;
            high_[i] = {node.re.hi, node.im.hi};
            low_[i] = {node.re.lo, node.im.lo};
            unresolved_[i] = residue_unresolved(residues[i].value()) || (model && model->unsteered(revision, i));
            // As in SecularPoints::prepare: a residue below 2^-2b of its node (b = binary64_bits) is zero here.
            Magnitude bound;
            acb_get_mag(bound.value(), residues[i].value());
            if (mag_get_d_log2_approx(bound.value()) < machine::logabs(high_[i]) - 2 * machine::binary64_bits)
                residues_[i] = C{};
            else if (!from_ball(residues_[i], residues[i].value()))
                return false;
        }
        // Two equal nodes refuse the model (sorted, O(n log n), as in SecularPoints::prepare).
        std::vector<size_t> order(n);
        std::iota(order.begin(), order.end(), size_t{0});
        const auto key = [&](size_t i) { return std::tuple{high_[i].re, high_[i].im, low_[i].re, low_[i].im}; };
        std::sort(order.begin(), order.end(), [&](size_t a, size_t b) { return key(a) < key(b); });
        for (size_t k = 1; k < order.size(); ++k)
            if (key(order[k]) == key(order[k - 1])) return false;
        return true;
    }
    // The sub-floor steps of the active points that sit at their own nodes (a Jacobi sweep from the packet's
    // output), written into `out`. Returns the number of steps taken, or nothing if this arithmetic refused (the
    // packet's output then stands).
    std::optional<size_t> refine(acb_ptr out, std::span<const ComplexBall> nodes, std::span<const uint8_t> inactive,
                                 slong stage, certpolysolve::detail::CandidateStatistics &s) {
        const size_t n = high_.size();
        if (nodes.size() != n || inactive.size() != n) return std::nullopt;
        offsets_.resize(n);
        steps_.resize(n);
        refined_.assign(n, 0);
        Float difference;
        for (size_t i = 0; i < n; ++i) {
            // d_i = z_i - b_i rounded once to binary64 (arf_sub rounds the exact difference).
            const auto k = static_cast<slong>(i);
            arf_sub(difference.value(), arb_midref(acb_realref(out + k)), arb_midref(acb_realref(nodes[i].value())),
                    machine::binary64_bits, ARF_RND_NEAR);
            if (!from_arf(offsets_[i].re, difference.value())) return std::nullopt;
            arf_sub(difference.value(), arb_midref(acb_imagref(out + k)), arb_midref(acb_imagref(nodes[i].value())),
                    machine::binary64_bits, ARF_RND_NEAR);
            if (!from_arf(offsets_[i].im, difference.value())) return std::nullopt;
        }
        size_t taken = 0;
        for (size_t i = 0; i < n; ++i) {
            if (inactive[i] || offsets_[i] != C{} || unresolved_[i]) continue;
            // Newton's correction at the node, q = w_i / D_i with D_i = sum_{j != i}(w_j + w_i)/(b_i - b_j) - 1
            // (SecularPoints::ratio's first branch).
            C sum;
            ++s.evaluations;
            s.evaluation_terms += n;
            for (size_t j = 0; j < n; ++j)
                if (j != i) sum += (residues_[j] + residues_[i]) / node_difference(i, j);
            const C denominator = sum - C(1.0);
            if (denominator == C{}) continue; // the Newton limit is infinite here (NewtonPole): no sub-floor step
            const C q = residues_[i] / denominator;
            const double relative = machine::logabs(q) - machine::logabs(high_[i]);
            if (!machine::finite(q) || !(relative < 2.0 - machine::binary64_bits) ||
                relative < 2.0 - static_cast<double>(stage))
                continue;
            // Aberth's step from the node, with the others where the packet left them.
            C interaction;
            for (size_t j = 0; j < n; ++j)
                if (j != i) {
                    const C between = node_difference(i, j) - offsets_[j];
                    if (between == C{}) return std::nullopt;
                    interaction += machine::interaction_term(between);
                    ++s.root_interactions;
                }
            const C step = machine::aberth_correction(q, interaction);
            if (!machine::finite(step)) return std::nullopt;
            steps_[i] = step;
            refined_[i] = 1;
            ++taken;
        }
        Float offset;
        for (size_t i = 0; i < n; ++i) {
            if (!refined_[i]) continue;
            // node - step, exactly: the step's bits lie below the node's, and no rounding may drop them.
            const auto k = static_cast<slong>(i);
            arf_set_d(offset.value(), -steps_[i].re);
            arf_add(arb_midref(acb_realref(out + k)), arb_midref(acb_realref(nodes[i].value())), offset.value(),
                    ARF_PREC_EXACT, ARF_RND_DOWN);
            mag_zero(arb_radref(acb_realref(out + k)));
            arf_set_d(offset.value(), -steps_[i].im);
            arf_add(arb_midref(acb_imagref(out + k)), arb_midref(acb_imagref(nodes[i].value())), offset.value(),
                    ARF_PREC_EXACT, ARF_RND_DOWN);
            mag_zero(arb_radref(acb_imagref(out + k)));
        }
        return taken;
    }
};
}
struct SecularMachineWorkspace::Impl {
    SecularPoints<double> ordinary;
    SecularPoints<machine::ScaledDouble> scaled;
    NodeOffsets offsets;
    uint64_t revision = UINT64_MAX, offsets_revision = UINT64_MAX;
    bool ordinary_ready = false, scaled_ready = false, offsets_ready = false;
};
SecularMachineWorkspace::SecularMachineWorkspace() = default;
slong SecularMachineWorkspace::stage_bits(uint64_t goal_bits) {
    return std::max<slong>(64, resources::precision_sum(goal_bits, 32));
}
SecularMachineWorkspace::~SecularMachineWorkspace() = default;
bool
SecularMachineWorkspace::advance(CandidateProgress &progress,
                                 acb_ptr out,
                                 acb_srcptr initial,
                                 std::span<const ComplexBall> nodes,
                                 std::span<const ComplexBall> residues,
                                 std::span<const uint8_t> inactive,
                                 std::span<const slong> bits,
                                 uint64_t revision,
                                 unsigned sweeps,
                                 const cps_options &o,
                                 certpolysolve::detail::CandidateStatistics &s,
                                 const EnginePolicy &policy, ModelStoppingState *model, bool keep_unmoved,
                                 bool floor_steps) {
    progress.service_scaled_binary64 = 0;
    const auto n = nodes.size();
    if (!n || residues.size() != n || inactive.size() != n || bits.size() != n || o.max_work_bits < machine::binary64_bits ||
        std::fegetround() != FE_TONEAREST)
        return false;
    // Try machine iteration during the initial refinement stage. Once the
    // existing per-root scheduler grows precision, retain its adaptive choice.
    const auto initial_bits = stage_bits(o.bits);
    for (size_t i = 0; i < n; ++i)
        if (!inactive[i] && bits[i] > initial_bits) return false;
    fenv_t environment;
    if (std::feholdexcept(&environment) != 0) return false;
    struct Restore {
        fenv_t *p;
        ~Restore() { std::fesetenv(p); }
    } restore{ &environment };
    if (!impl_) impl_ = std::make_unique<Impl>();
    auto &cache = *impl_;
    auto *counts = policy.study_statistics.get();
    StageTimer timer(s.candidate_ns);
    if (cache.revision != revision) {
        TraceScope trace(&policy, "study_secular_preparation", machine::binary64_bits, n, &s);
        cache.revision = revision;
        cache.ordinary_ready = cache.ordinary.prepare(nodes, residues, model, revision);
        cache.scaled_ready = false;
        if (!cache.ordinary_ready) cache.scaled_ready = cache.scaled.prepare(nodes, residues, model, revision);
        if (counts) ++counts->preparations;
    } else if (counts)
        ++counts->reuses;
    if (!cache.ordinary_ready && !cache.scaled_ready) return false;
    std::feclearexcept(FE_ALL_EXCEPT);
    const bool ordinary = cache.ordinary_ready;
    // The exponent-managed scalar kernel is different from the ordinary
    // binary64 kernel. Retain its actual terms even on refusal/cancellation.
    struct MachineMeter {
        CandidateProgress &progress;
        const CandidateStatistics &statistics;
        bool enabled;
        uint64_t before;
        ~MachineMeter() {
            if (enabled)
                progress.service_scaled_binary64 += statistics.evaluation_terms + statistics.root_interactions - before;
        }
    } meter{progress, s, policy.research_service_meter && !ordinary,
            s.evaluation_terms + s.root_interactions};
    if (policy.model_stopping == SecularModelPolicy::Legacy) model = nullptr;
    std::vector<slong> machine_bits;
    if (model) {
        machine_bits.assign(n, machine::binary64_bits);
        model->prepare(revision, ordinary ? 1u : 2u, machine_bits, inactive);
    }
    bool accepted = ordinary ? cache.ordinary.advance(progress, out, initial, inactive, sweeps, o, s, model, keep_unmoved,
                                                      policy.packet_lane_sums)
                             : cache.scaled.advance(progress, out, initial, inactive, sweeps, o, s, model, keep_unmoved);
    const int exceptions = std::fetestexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | (ordinary ? FE_UNDERFLOW : 0));
    if (progress.cancelled) return true;
    accepted = accepted && !exceptions;
    if (counts) {
        ++counts->machine_packets;
        if (!ordinary) ++counts->scaled_packets;
        if (!accepted) ++counts->continuations;
    }
    if (policy.trace) event(policy,
          s,
          "study_secular_candidates",
          ordinary ? "binary64" : "scaled-binary64",
          accepted ? "regenerated-representation" : "adaptive-continuation",
          *std::max_element(bits.begin(), bits.end()),
          static_cast<size_t>(std::count(inactive.begin(), inactive.end(), uint8_t(0))),
          machine::binary64_bits, revision);
    if (!accepted) {
        // Do not retry an unchanged failing representation on every packet.
        cache.ordinary_ready = false;
        cache.scaled_ready = false;
        progress = {};
        if (model) model->invalidate();
    }
    if (accepted && model) model->trace(policy, s, inactive, ordinary ? "binary64" : "scaled-binary64", machine::binary64_bits);
    // At a fresh model's floor, the steps the packet could not take (NodeOffsets), when the caller asks
    // (`floor_steps`: EnginePolicy::machine_node_offsets, outside a higher goal's staged default solve) and the goal
    // is at most the default. A refusal keeps the packet's output.
    if (accepted && ordinary && model && progress.model_exhausted && progress.model_immediate && floor_steps &&
        o.bits <= cps_options{}.bits) {
        if (cache.offsets_revision != revision) {
            cache.offsets_revision = revision;
            cache.offsets_ready = cache.offsets.prepare(nodes, residues, model, revision);
        }
        if (cache.offsets_ready) {
            RootVector refined(static_cast<slong>(n));
            _acb_vec_set(refined.data(), out, static_cast<slong>(n));
            std::feclearexcept(FE_ALL_EXCEPT);
            const auto taken = cache.offsets.refine(refined.data(), nodes, inactive, initial_bits, s);
            if (taken && !std::fetestexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW)) {
                _acb_vec_set(out, refined.data(), static_cast<slong>(n));
                if (policy.trace)
                    event(policy, s, "study_secular_candidates", "offset-binary64", "floor-steps",
                          *std::max_element(bits.begin(), bits.end()), *taken, machine::binary64_bits, revision);
            } else
                cache.offsets_ready = false;
        }
    }
    return accepted;
}
} // namespace certpolysolve::detail
