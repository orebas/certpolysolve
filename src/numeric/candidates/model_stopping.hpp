#pragma once
#include "core/engine_policy.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <span>

namespace certpolysolve::detail {
enum class ModelStop { Active, EvaluationFloor, SmallCorrection };

// log2 of the rounding error of a sum of `terms` model terms at `bits` bits,
// relative to the sum of their absolute values. Rounded-model scheduling
// information, never an original-input certificate: a serial operation-count
// scale rather than MPSolve's pairwise-sum constant.
inline double evaluation_roundoff_log2(size_t terms, slong bits) {
    return std::log2(32.0 * static_cast<double>(terms) + 16.0) + 1.0 - static_cast<double>(bits);
}

// Original-input inclusion predicates remain the acceptance boundary.
inline ModelStop model_stop(double value_log, double terms_log, double correction_log, double point_log,
                            size_t terms, slong bits) {
    if (std::isfinite(terms_log) && !std::isnan(value_log)) {
        if (terms_log + evaluation_roundoff_log2(terms, bits) >= value_log)
            return ModelStop::EvaluationFloor;
    }
    if (correction_log == -std::numeric_limits<double>::infinity() ||
        (std::isfinite(point_log) && correction_log - point_log < 2.0 - static_cast<double>(bits)))
        return ModelStop::SmallCorrection;
    return ModelStop::Active;
}

// At a node b_k the model's value, normalised by the other poles, is the
// residue itself: p(b_k) = -w_k prod_{j!=k}(b_k - b_j). A residue enclosure
// with no correct bits therefore leaves p(b_k) indistinguishable from zero:
// b_k is a root at the precision that produced w_k, and no step the model
// proposes for it is informed. This is the evaluation floor at a pole.
// Regeneration accepts such a residue only at its precision budget (adequate()
// otherwise demands positive relative accuracy), so this marks the nodes the
// budget blocks. An exact zero residue is an exact root, not a floor.
inline bool residue_unresolved(const acb_t residue) {
    return !acb_is_zero(residue) && acb_rel_accuracy_bits(residue) <= 0;
}

// Where Newton's correction p/p' is infinite in the model's arithmetic: a
// node b_k whose limit denominator D_k = sum_{j != k}(w_j + w_k)/(b_k - b_j) - 1
// rounds to exactly zero. The step a sweep takes is Aberth's,
//     w_k / (D_k - w_k A),   A = sum_{j != k} 1/(z_k - z_j),
// and its limit for D_k -> 0 is -1/A, the repulsion from the other points
// alone, which is how Aberth's method steps at a critical point. (A tiny
// nonzero D_k gives the same step through a huge Newton correction; only the
// exact zero needs the limit.) The limit is the step when w_k A stands above
// the rounding error of D_k: the computed D_k is then noise beside the term
// that decides D_k - w_k A. Otherwise that denominator is below its rounding
// error and this arithmetic determines no step. `margin` is log2|w_k| less
// log2 of D_k's rounding error; the caller supplies log2|A|. A residue with
// no correct bit (residue_unresolved) determines nothing either way.
struct NewtonPole {
    bool at = false;
    bool unresolved = false;
    double margin = 0;
    bool resolved(double interaction_log) const { return at && !unresolved && margin + interaction_log > 0; }
};

class ModelStoppingState {
    std::vector<ModelStop> roots_;
    std::vector<slong> bits_;
    uint64_t revision_ = 0;
    unsigned tier_ = 0;
    bool ready_ = false;
    size_t initial_active_ = 0;
    uint64_t evaluations_ = 0;
    std::vector<uint8_t> unsteered_;
    uint64_t unsteered_revision_ = 0;

  public:
    void invalidate() { ready_ = false; }
    // The owning refinement's verdict on the model of this revision: the nodes
    // whose step that model leaves undetermined (unsteered_nodes). The packet
    // engines settle such a node at its first evaluation, as they do a node
    // whose own residue is unresolved.
    void steering(uint64_t revision, std::vector<uint8_t> unsteered) {
        unsteered_revision_ = revision;
        unsteered_ = std::move(unsteered);
    }
    bool unsteered(uint64_t revision, size_t i) const {
        return unsteered_revision_ == revision && i < unsteered_.size() && unsteered_[i];
    }
    void prepare(uint64_t revision, unsigned tier, std::span<const slong> bits,
                 std::span<const uint8_t> certified) {
        if (!ready_ || revision_ != revision || tier_ != tier || bits_.size() != bits.size() ||
            !std::equal(bits_.begin(), bits_.end(), bits.begin())) {
            roots_.assign(bits.size(), ModelStop::Active);
            bits_.assign(bits.begin(), bits.end());
            revision_ = revision;
            tier_ = tier;
            initial_active_ = static_cast<size_t>(std::count(certified.begin(), certified.end(), uint8_t(0)));
            evaluations_ = 0;
            ready_ = true;
        }
    }
    bool settled(size_t i) const { return ready_ && roots_.at(i) != ModelStop::Active; }
    ModelStop reason(size_t i) const { return ready_ ? roots_.at(i) : ModelStop::Active; }
    void evaluated() { ++evaluations_; }
    void evaluated(uint64_t count) { evaluations_ += count; } // a sweep's evaluations at once (W9's split packet)
    void stop(size_t i, ModelStop reason) { roots_.at(i) = reason; }
    size_t count(ModelStop reason) const {
        return static_cast<size_t>(std::count(roots_.begin(), roots_.end(), reason));
    }
    uint64_t evaluations() const { return evaluations_; }
    bool exhausted(std::span<const uint8_t> certified) const {
        if (!ready_ || certified.size() != roots_.size())
            return false;
        for (size_t i = 0; i < roots_.size(); ++i)
            if (!certified[i] && !settled(i))
                return false;
        return true;
    }
    bool immediate() const { return ready_ && evaluations_ <= initial_active_; }
    void trace(const EnginePolicy &policy, const certpolysolve::detail::CandidateStatistics &statistics,
               std::span<const uint8_t> certified, const char *tier, slong bits) const {
        if (!policy.trace)
            return;
        TraceEvent e{"secular_model_stopping",
                     "packet",
                     static_cast<uint64_t>(bits),
                     static_cast<uint64_t>(std::count(certified.begin(), certified.end(), uint8_t(0))),
                     &statistics,
                     nullptr};
        e.engine = "secular";
        e.tier = tier;
        e.reason =
            exhausted(certified) ? (immediate() ? "fresh-model-floor" : "model-exhausted") : "model-progress";
        e.representation_revision = revision_;
        e.arithmetic_bits = static_cast<uint64_t>(bits);
        e.model_evaluations = evaluations_;
        e.model_error_stops = count(ModelStop::EvaluationFloor);
        e.model_small_stops = count(ModelStop::SmallCorrection);
        e.model_settled = e.model_error_stops + e.model_small_stops;
        policy.trace(policy.trace_data, e);
    }
};
} // namespace certpolysolve::detail
