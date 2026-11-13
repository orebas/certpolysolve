// Program candidate iteration: batched evaluator readings, Jacobi Aberth
// corrections and the binary64 stage's bounded stopping rules. Every published
// root still needs a later original-source certificate.
#include "numeric/candidates/program_candidates.hpp"
#include "numeric/candidates/aberth_step.hpp"
#include "numeric/candidates/interaction_tree.hpp"
#include "numeric/arithmetic/machine_conversion.hpp"
#include "core/heuristics.hpp"
#include "core/resource_limits.hpp"
#include "core/telemetry.hpp"
#include <algorithm>
#include <cfenv>
#include <cmath>
#include <complex>
#include <limits>

namespace certpolysolve::detail {
namespace {
double relative_step_log2(machine::Complex<double> delta, machine::Complex<double> z) {
    const double d = machine::logabs(delta), m = machine::logabs(z);
    if (d == -std::numeric_limits<double>::infinity())
        return d;
    return d - m; // -inf - (-inf) cannot arise: a zero z with a nonzero delta gives +inf, a large step
}
} // namespace

ProgramMachineStage::ProgramMachineStage(const Program &program, acb_srcptr seeds, size_t n, size_t chunk,
                                          uint64_t storage_limit, const uint8_t *active, uint64_t zero_order)
    : evaluator_(program, zero_order), n_(n), chunk_(std::max<size_t>(1, std::min(chunk, std::max<size_t>(n, 1)))) {
    if (!evaluator_.valid() || n == 0 || std::fegetround() != FE_TONEAREST)
        return;
    if (!resources::storage_fits(resources::add(evaluator_.storage_bytes(chunk_),
                                                resources::multiply(static_cast<uint64_t>(n), 4 * sizeof(C) + 16)),
                                 storage_limit))
        return;
    points_.resize(n);
    for (size_t i = 0; i < n; ++i)
        if (!from_ball(points_[i], seeds + static_cast<slong>(i)))
            return;
    stops_.assign(n, Stop::Active);
    if (active)
        for (size_t i = 0; i < n; ++i)
            stops_[i] = active[i] ? Stop::Active : Stop::Converged;
    best_step_.assign(n, std::numeric_limits<double>::infinity());
    stagnant_.assign(n, 0);
    sums_.resize(n);
    previous_window_ = std::numeric_limits<double>::infinity();
    current_window_ = -std::numeric_limits<double>::infinity();
    valid_ = true;
}

ProgramMachineStage::Advance ProgramMachineStage::state() const {
    Advance a;
    for (const auto stop : stops_) {
        a.active += stop == Stop::Active;
        a.converged += stop == Stop::Converged;
        a.floor += stop == Stop::Floor;
        a.held += stop == Stop::Held;
    }
    a.finished = finished_;
    a.stalled = stalled_;
    a.scaled_evaluations = scaled_evaluations_;
    return a;
}

ProgramMachineStage::Advance
ProgramMachineStage::advance(unsigned sweeps, const cps_options &options, certpolysolve::detail::CandidateStatistics &statistics) {
    uint64_t taken = 0;
    for (unsigned s = 0; s < sweeps && !finished_; ++s) {
        if (options.cancelled()) {
            auto a = state();
            a.sweeps = taken;
            a.cancelled = true;
            return a;
        }
        sweep(statistics);
        ++taken;
    }
    auto a = state();
    a.sweeps = taken;
    return a;
}

void ProgramMachineStage::sweep(certpolysolve::detail::CandidateStatistics &statistics) {
    StageTimer timer(statistics.candidate_ns);
    const double floor_factor = static_cast<double>(heuristics::program_machine_floor_factor);
    active_.clear();
    for (size_t i = 0; i < n_; ++i)
        if (stops_[i] == Stop::Active)
            active_.push_back(i);
    if (active_.empty()) {
        finished_ = true;
        return;
    }
    // Every active point evaluated at the sweep's start positions, in one batched call.
    batch_.resize(active_.size());
    evaluations_.resize(active_.size());
    for (size_t k = 0; k < active_.size(); ++k)
        batch_[k] = points_[active_[k]];
    // Per point, from either tier: finite, at the floor, and the Newton correction.
    struct Reading {
        bool finite = false, value_floor = false, derivative_floor = false;
        C newton;
    };
    std::vector<Reading> readings(active_.size());
    {
        StageTimer evaluation(statistics.evaluation_ns);
        evaluator_.evaluate(batch_, evaluations_, chunk_);
        statistics.evaluations += active_.size();
        statistics.evaluation_terms += static_cast<uint64_t>(active_.size()) * evaluator_.terms();
        std::vector<size_t> overflowing;
        for (size_t k = 0; k < active_.size(); ++k) {
            const auto &e = evaluations_[k];
            auto &r = readings[k];
            if (!e.finite) {
                overflowing.push_back(k);
                continue;
            }
            r.finite = true;
            r.derivative_floor = e.derivative == C{} || machine::norm(e.derivative) < floor_factor * e.noise_derivative;
            r.value_floor = e.value == C{} || machine::norm(e.value) < floor_factor * e.noise;
            if (!r.derivative_floor)
                r.newton = e.value / e.derivative;
        }
        // The range tier for the points whose values left binary64, this sweep only.
        if (!overflowing.empty()) {
            using S = machine::BlockDouble;
            using SC = machine::Complex<S>;
            if (!scaled_)
                scaled_ = std::make_unique<ProgramMachine<S>>(evaluator_);
            if (scaled_->valid()) {
                std::vector<SC> points(overflowing.size());
                for (size_t m = 0; m < overflowing.size(); ++m) {
                    const C z = batch_[overflowing[m]];
                    points[m] = {S(z.re), S(z.im)};
                }
                std::vector<ProgramMachine<S>::Evaluation> scaled(overflowing.size());
                scaled_->evaluate(points, scaled, chunk_);
                scaled_evaluations_ += overflowing.size();
                statistics.evaluations += overflowing.size();
                statistics.evaluation_terms += static_cast<uint64_t>(overflowing.size()) * scaled_->terms();
                const S factor = S(heuristics::program_machine_floor_factor);
                for (size_t m = 0; m < overflowing.size(); ++m) {
                    const auto &e = scaled[m];
                    auto &r = readings[overflowing[m]];
                    if (!e.finite)
                        continue; // held: out of the range of both tiers
                    r.derivative_floor = e.derivative == SC{} || machine::norm(e.derivative) < factor * e.noise_derivative;
                    r.value_floor = e.value == SC{} || machine::norm(e.value) < factor * e.noise;
                    if (!r.derivative_floor) {
                        const SC newton = e.value / e.derivative;
                        r.newton = {static_cast<double>(newton.re), static_cast<double>(newton.im)};
                        if (!machine::finite(r.newton))
                            continue; // a correction beyond binary64: held
                    }
                    r.finite = true;
                }
            }
        }
    }
    // The interaction sums of the active points over every point's start position, in binary64, the stage's arithmetic
    // (only values leave the range, which is what the range tier serves). From machine_interaction_tree_degree on,
    // while program_machine_tree_active points or more are active (fewer cost less summed directly than a tree built
    // over all n), the interaction tree built once a sweep from the start positions (every sum reads them: Jacobi needs
    // no move) gives each active point its far field and the near points to sum directly; below it, the direct sum
    // transposed: for each j, every active i accumulates its term, so each sum is its own fixed-order accumulator
    // (elementwise across i). A coincident pair holds the point for the sweep.
    std::vector<uint8_t> coincident(active_.size(), 0);
    {
        StageTimer interaction(statistics.interaction_ns);
        using D = machine::Complex<double>;
        // Corrections update points_ only after every interaction sum is formed.
        // This borrow lasts just for that read-only phase; no snapshot copy is needed.
        const auto &positions = points_;
        std::fill_n(sums_.begin(), active_.size(), C{});
        std::vector<uint8_t> summed(active_.size(), 0);
        if (n_ >= heuristics::machine_interaction_tree_degree &&
            active_.size() >= heuristics::program_machine_tree_active) {
            std::vector<std::complex<double>> standard(n_);
            for (size_t j = 0; j < n_; ++j)
                standard[j] = {positions[j].re, positions[j].im};
            const InteractionTreeOf<heuristics::program_tree_order, heuristics::ProgramFarField> tree(standard);
            statistics.root_interactions += tree.build_terms();
            std::vector<size_t> near;
            for (size_t k = 0; k < active_.size(); ++k) {
                const size_t i = active_[k];
                std::complex<double> far;
                double tail = 0;
                size_t far_fields = 0;
                if (!tree.evaluate(i, far, near, tail, far_fields))
                    continue; // the direct sum below
                D sum{far.real(), far.imag()};
                for (const size_t j : near) {
                    const D d = positions[i] - positions[j];
                    const double q = d.re * d.re + d.im * d.im;
                    if (q == 0) {
                        coincident[k] = 1;
                        continue;
                    }
                    sum += D{d.re / q, -d.im / q};
                }
                sums_[k] = sum;
                summed[k] = 1;
                statistics.root_interactions += near.size() + far_fields * heuristics::program_tree_order;
            }
        }
        // Only the active points the tree did not serve: with the tree serving them all, the loop over every pair would
        // cost n^2 empty steps a sweep (a billion at degree 32,767). Each sum still accumulates its terms in j order,
        // so its bits are the same.
        std::vector<size_t> unsummed;
        for (size_t k = 0; k < active_.size(); ++k)
            if (!summed[k])
                unsummed.push_back(k);
        if (!unsummed.empty())
            for (size_t j = 0; j < n_; ++j) {
                const D zj = positions[j];
                for (const size_t k : unsummed) {
                    const size_t i = active_[k];
                    if (i == j)
                        continue;
                    const D d = positions[i] - zj;
                    const double q = d.re * d.re + d.im * d.im;
                    if (q == 0) {
                        coincident[k] = 1;
                        continue;
                    }
                    sums_[k] += D{d.re / q, -d.im / q}; // 1 / d
                }
            }
        const uint64_t direct = unsummed.size();
        statistics.root_interactions += direct * (n_ - 1);
    }
    // The corrections, from the start positions; the stops; then the update. A point's patience is the
    // window: during Aberth's sorting-out of the angles a point's step can hold for many sweeps before it
    // falls, and a shorter patience floors points that are still sorting out.
    // History: docs/slp-seed-placement-note-20261016.md, section 9.
    const uint32_t patience = std::max<uint32_t>(heuristics::program_machine_patience, static_cast<uint32_t>((n_ + 2) / 2));
    double largest = -std::numeric_limits<double>::infinity();
    for (size_t k = 0; k < active_.size(); ++k) {
        const size_t i = active_[k];
        const auto &r = readings[k];
        if (!r.finite) {
            stops_[i] = Stop::Held;
            continue;
        }
        const C z = points_[i];
        if (r.derivative_floor) {
            stops_[i] = Stop::Floor;
            continue;
        }
        if (coincident[k]) {
            // Held this sweep; held for as long as its patience, it is at the floor (a pair that never parts would
            // otherwise keep the stage running while any other point moves).
            if (++stagnant_[i] >= patience)
                stops_[i] = Stop::Floor;
            continue;
        }
        const bool value_floor = r.value_floor;
        const C newton = r.newton;
        const C delta = machine::aberth_correction(newton, sums_[k]);
        const C next = z - delta;
        if (!machine::finite(delta) || !machine::finite(next)) {
            if (++stagnant_[i] >= patience)
                stops_[i] = Stop::Floor;
            continue; // held this sweep
        }
        const double step = relative_step_log2(delta, z);
        const bool small = step < static_cast<double>(heuristics::program_machine_small_step_log2);
        if (value_floor) {
            stops_[i] = small ? Stop::Converged : Stop::Floor;
            if (small)
                points_[i] = next;
            continue;
        }
        points_[i] = next;
        if (small) {
            stops_[i] = Stop::Converged;
            continue;
        }
        if (step < best_step_[i] - 1) {
            best_step_[i] = step;
            stagnant_[i] = 0;
        } else if (++stagnant_[i] >= patience) {
            stops_[i] = Stop::Floor;
            continue;
        }
        largest = std::max(largest, step);
    }
    ++statistics.candidate_sweeps;
    last_largest_ = largest;
    // The window rule (program_advance's, solver_program.cpp): a window whose largest step does not halve against the
    // smallest window so far is a stall, and so is a window in which no point moved.
    const uint64_t window = (n_ + 2) / 2;
    current_window_ = std::max(current_window_, largest);
    ++window_sweeps_;
    if (window_sweeps_ >= window) {
        const bool moved = current_window_ > -std::numeric_limits<double>::infinity();
        const bool stalled = !moved || current_window_ > previous_window_ - 1;
        previous_window_ = std::min(previous_window_, current_window_);
        current_window_ = -std::numeric_limits<double>::infinity();
        window_sweeps_ = 0;
        if (stalled) {
            stalled_ = true;
            for (auto &stop : stops_)
                if (stop == Stop::Active)
                    stop = Stop::Floor;
        }
    }
    bool any_active = false;
    for (const auto stop : stops_)
        any_active = any_active || stop == Stop::Active;
    finished_ = !any_active;
}

void ProgramMachineStage::export_points(acb_ptr out) const {
    for (size_t i = 0; i < n_; ++i)
        to_ball(out + static_cast<slong>(i), points_[i]);
}

} // namespace certpolysolve::detail
