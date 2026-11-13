#pragma once
// Cached coefficient evaluators used by secular residue regeneration. A scoped
// selection binds the evaluator, point/origin, precision and value/joint mode.
#include "numeric/evaluate/structured_source.hpp"
#include "numeric/evaluate/disk_blocked_polynomial_evaluator.hpp"
#include "numeric/evaluate/disk_sparse_polynomial_evaluator.hpp"
#include "numeric/evaluate/limb_horner_evaluator.hpp"
#include <cmath>
#include <map>
#include <memory>
#include <optional>

namespace certpolysolve::detail {
// The source's value at a node with disk errors, for a coefficient source:
// block dot products for a dense one (DiskBlockedPolynomialEvaluator), the
// terms' powers for a sparse one (DiskSparsePolynomialEvaluator), one
// evaluator for each precision the ladder visits; exact translations and
// recurrences keep the source's own evaluation. That evaluation is ball
// Horner, which re-boxes its error rectangle at every term: up to half a bit
// each (a factor |cos t| + |sin t| at argument t), the degree's half in all,
// and the residues' ladder then climbs far beyond the precision that suffices.
//
// The blocks serve while the precision is at most the degree, where that loss
// can be half the precision or more and the blocks are also the faster. Above
// it the loss is a smaller share of the precision, and the blocks' long
// products cost more where the nodes are short: the degree is the crossover
// measured against ball Horner. History:
// docs/code-history-20261021.md, secular_refinement.cpp 1.
class DiskValues {
    PolynomialSource *source_;
    // Whether a node whose parts are short (a binary64 node, or a last regeneration's of 115 bits) takes the
    // evaluation written for one (limb_horner_evaluator.hpp): Horner's rule on limbs with a short multiplier
    // and a disk error, at every precision. It loses nothing to the rectangles either, so it also serves
    // above the degree, where the source's ball Horner still loses up to half a bit a term. It follows the
    // source's own rule for a centred polynomial (PolynomialSource::origin_point): a node nearer the
    // original origin is evaluated on the uncentred coefficients at its distance from that origin, where it
    // is short too. History: docs/code-history-20261021.md, secular_refinement.cpp 2.
    bool limb_horner_;
    struct Evaluators {
        std::unique_ptr<DiskBlockedPolynomialEvaluator> dense;
        std::unique_ptr<DiskSparsePolynomialEvaluator> sparse;
        std::unique_ptr<LimbHornerEvaluator> limb_horner, origin_limb_horner;
        ComplexBall leading, origin_leading;
        bool is_dense = false;
    };
    std::map<slong, Evaluators> cache_;
    ComplexBall scratch_, from_origin_;

  public:
    // Borrows its DiskValues owner and input node. Both must remain alive and
    // the node unchanged through evaluation; the next prepare() invalidates a
    // selection (origin coordinates and evaluation scratch are reused). Keep it
    // in the caller's immediate evaluation scope. Preparation selects a schedule;
    // evaluation can still refuse a point, such as zero for a sparse joint call.
    class Prepared {
        friend class DiskValues;
        DiskValues &owner_;
        Evaluators &entry_;
        const ComplexBall &point_, &leading_;
        slong precision_;
        bool value_only_;
        LimbHornerEvaluator *limb_;

        Prepared(DiskValues &owner, Evaluators &entry, const ComplexBall &point,
                 const ComplexBall &leading, slong precision, bool value_only,
                 LimbHornerEvaluator *limb = nullptr)
            : owner_(owner), entry_(entry), point_(point), leading_(leading),
              precision_(precision), value_only_(value_only), limb_(limb) {}

      public:
        Prepared(const Prepared &) = delete;
        Prepared(Prepared &&) noexcept = default;
        Prepared &operator=(const Prepared &) = delete;
        Prepared &operator=(Prepared &&) = delete;
        bool uses_limb_horner() const { return limb_ != nullptr; }

        // The derivative is untouched for a value-only selection. The point,
        // precision and mode cannot be replaced between preparation and use.
        double research_limb_work(size_t terms) const {
            // One coefficient pass plus short-multiplier limb passes over
            // the evaluator's actual window, at the actual mapped point.
            const double multiplier = std::ceil(static_cast<double>(acb_bits(point_.value())) / FLINT_BITS);
            return static_cast<double>(terms) * LimbHornerEvaluator::limbs_for(precision_) * (1 + multiplier);
        }
        bool evaluate(acb_t value, acb_t derivative) const {
            if (value_only_) {
                const bool done = limb_ ? limb_->evaluate_value(owner_.scratch_.value(), point_.value())
                                  : entry_.dense ? entry_.dense->evaluate_value(owner_.scratch_.value(), point_.value())
                                                 : entry_.sparse->evaluate_value(owner_.scratch_.value(), point_.value());
                if (!done)
                    return false;
                acb_div(value, owner_.scratch_.value(), leading_.value(), precision_);
                return acb_is_finite(value);
            }
            const bool done = entry_.dense ? entry_.dense->evaluate(value, derivative, point_.value())
                                          : entry_.sparse->evaluate(value, derivative, point_.value());
            if (!done)
                return false;
            acb_div(value, value, leading_.value(), precision_);
            acb_div(derivative, derivative, leading_.value(), precision_);
            return acb_is_finite(value) && acb_is_finite(derivative);
        }
    };

  private:
    // The source's coefficients at a precision, for an evaluator to own. PolynomialSource evaluates p / a_n.
    // The evaluators get p itself and the value is divided once, as the source does: dividing the coefficients
    // first would make short exact coefficients full-precision numbers.
    ComplexPolynomial coefficients(slong precision) {
        ComplexPolynomial copy;
        acb_poly_set(copy.value(), source_->coefficient_enclosure(precision)->polynomial.value());
        return copy;
    }
    // What serves a precision, or nothing where the source has no coefficient schedule there.
    Evaluators *evaluators(slong precision) {
        if (const auto found = cache_.find(precision); found != cache_.end())
            return &found->second;
        const bool dense = source_->rectangular_available(precision);
        if (!dense && !source_->sparse_available(precision))
            return nullptr;
        const auto enclosed = source_->coefficient_enclosure(precision);
        const auto p = enclosed->polynomial.value();
        const slong n = acb_poly_degree(p);
        if (n < 1 || n != source_->degree() || acb_contains_zero(p->coeffs + n))
            return nullptr;
        Evaluators entry;
        entry.is_dense = dense;
        acb_set(entry.leading.value(), p->coeffs + n);
        return &cache_.emplace(precision, std::move(entry)).first->second;
    }
    // The limb evaluation of a short node, on the coefficients the source's rule gives it.
    std::optional<Prepared> choose_limb_horner(Evaluators &entry, slong precision, const ComplexBall &node) {
        const bool origin = source_->origin_point(from_origin_.value(), node.value());
        if (!LimbHornerEvaluator::serves(origin ? from_origin_.value() : node.value()))
            return std::nullopt;
        auto &evaluator = origin ? entry.origin_limb_horner : entry.limb_horner;
        if (!evaluator) {
            ComplexPolynomial copy;
            if (origin) {
                const auto p = source_->origin_enclosure(precision)->polynomial.value();
                const slong n = acb_poly_degree(p);
                // An uncentred enclosure that does not show the degree leaves an unusable evaluator behind.
                if (n == source_->degree() && !acb_contains_zero(p->coeffs + n)) {
                    acb_poly_set(copy.value(), p);
                    acb_set(entry.origin_leading.value(), p->coeffs + n);
                }
            } else
                copy = coefficients(precision);
            evaluator = std::make_unique<LimbHornerEvaluator>(std::move(copy), precision);
        }
        if (!evaluator->usable())
            return std::nullopt;
        return Prepared(*this, entry, origin ? from_origin_ : node,
                        origin ? entry.origin_leading : entry.leading, precision, true, evaluator.get());
    }

  public:
    DiskValues(MonicSource &source, bool limb_horner)
        : source_(dynamic_cast<PolynomialSource *>(&source)), limb_horner_(limb_horner) {}
    // Selections borrow the owner's scratch and cache: its address stays fixed.
    DiskValues(const DiskValues &) = delete;
    DiskValues(DiskValues &&) = delete;
    DiskValues &operator=(const DiskValues &) = delete;
    DiskValues &operator=(DiskValues &&) = delete;

    // Prepare a disk evaluator for this precision, node and mode (built on first
    // use). No allocation per selection; reusable evaluators/scratch stay here.
    std::optional<Prepared> prepare(slong precision, const ComplexBall &node, bool value_only) & {
        if (!source_)
            return std::nullopt;
        const bool below_degree = precision <= source_->degree();
        const bool limbs = limb_horner_ && value_only;
        if (!below_degree && !limbs)
            return std::nullopt;
        Evaluators *entry = evaluators(precision);
        if (!entry)
            return std::nullopt;
        if (limbs && entry->is_dense)
            if (auto selected = choose_limb_horner(*entry, precision, node))
                return selected;
        if (!below_degree)
            return std::nullopt;
        if (entry->is_dense && !entry->dense)
            entry->dense = std::make_unique<DiskBlockedPolynomialEvaluator>(coefficients(precision), precision);
        if (!entry->is_dense && !entry->sparse)
            entry->sparse = std::make_unique<DiskSparsePolynomialEvaluator>(coefficients(precision), precision);
        return Prepared(*this, *entry, node, entry->leading, precision, value_only);
    }
    std::optional<Prepared> prepare(slong, const ComplexBall &, bool) && = delete;
    std::optional<Prepared> prepare(slong, ComplexBall &&, bool) & = delete;
    std::optional<Prepared> prepare(slong, const ComplexBall &&, bool) & = delete;
};
} // namespace certpolysolve::detail
