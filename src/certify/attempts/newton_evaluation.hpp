#pragma once
// Evaluation storage for one Newton-cover attempt. The source is borrowed for
// that call; normalized polynomials and all prepared evaluators are owned here.
// Work budgets, precision decisions, certificates and points belong to the caller.
#include "numeric/evaluate/disk_blocked_polynomial_evaluator.hpp"
#include "numeric/evaluate/disk_sparse_polynomial_evaluator.hpp"
#include "numeric/evaluate/polynomial_remainder.hpp"
#include "numeric/evaluate/structured_source.hpp"
#include <algorithm>
#include <memory>
#include <optional>
#include <utility>
#include <vector>

namespace certpolysolve::detail {

class NewtonCoverEvaluation {
    PolynomialSource &source_;
    const slong degree_;
    slong precision_ = 0;
    // Remainders and evaluations describe the same normalized polynomial p/a_n.
    std::optional<PolynomialLinearRemainder> bounds_;
    std::unique_ptr<DiskBlockedPolynomialEvaluator> blocks_;
    std::unique_ptr<DiskSparsePolynomialEvaluator> sparse_;
    struct Rung {
        slong bits = 0;
        std::unique_ptr<DiskBlockedPolynomialEvaluator> blocks;
        std::unique_ptr<DiskSparsePolynomialEvaluator> sparse;
    };
    std::vector<Rung> rungs_;

  public:
    NewtonCoverEvaluation(PolynomialSource &source, slong degree)
        : source_(source), degree_(degree) {}

    // The caller invokes this on each raise and rollback, at the original
    // accounting points. A false result requires refusal or another prepare;
    // it never authorizes evaluation at the requested precision. Invalid
    // degree/leading-coefficient checks precede changes to the old cache.
    bool prepare(slong bits) {
        // PolynomialSource evaluates p/a_n. The majorant must use those same
        // units, including small or nonreal leading coefficients.
        ComplexPolynomial monic;
        const auto enclosed = source_.coefficient_enclosure(bits);
        const auto p = enclosed->polynomial.value();
        if (acb_poly_degree(p) != degree_)
            return false;
        if (acb_contains_zero(p->coeffs + degree_))
            return false;
        acb_poly_scalar_div(monic.value(), p, p->coeffs + degree_, bits);
        acb_poly_set_coeff_si(monic.value(), degree_, 1);
        // Readiness is committed only after every required owner is prepared.
        // Invalid source checks above preserve the old cache; an allocation
        // exception during rebuilding cannot leave it advertised as ready.
        precision_ = 0;
        blocks_.reset();
        sparse_.reset();
        const bool dense = source_.rectangular_available(bits);
        const bool sparse = !dense && source_.sparse_available(bits);
        // Sparse evaluation owns its own remainder bound. Other paths need
        // only these magnitudes, so the dense evaluator can take the polynomial.
        if (sparse)
            bounds_.reset();
        else
            bounds_.emplace(monic.coefficients());
        if (dense) {
            blocks_ = std::make_unique<DiskBlockedPolynomialEvaluator>(std::move(monic), bits);
        } else if (sparse) {
            sparse_ = std::make_unique<DiskSparsePolynomialEvaluator>(monic, bits);
        }
        precision_ = bits;
        return true;
    }

    bool available() const { return precision_ != 0; }
    bool sparse() const { return static_cast<bool>(sparse_); }
    bool has_disk_evaluator() const { return blocks_ || sparse_; }

    // After a successful prepare, dense sources use block dot products and
    // sparse sources use disk-bounded powers. Sources with exact translations
    // retain their own evaluation. All paths enclose the same polynomial.
    bool evaluate(acb_t value, acb_t derivative, const acb_t center) {
        return blocks_   ? blocks_->evaluate(value, derivative, center)
               : sparse_ ? sparse_->evaluate(value, derivative, center)
                         : source_.evaluate(value, derivative, center, precision_);
    }

    void linear_remainder(arb_t out, const acb_t center, const arb_t radius) const {
        if (sparse_)
            sparse_->linear_remainder(out, center, radius);
        else
            bounds_->linear_remainder(out, center, radius);
    }

    // Lower Newton-ladder rungs are prepared lazily and retained for the whole
    // attempt, including certificate-precision raises and rollbacks. These
    // evaluations guide Newton steps; they do not certify roots on a rung.
    bool rung_evaluate(slong bits, acb_t value, acb_t derivative, const acb_t center) {
        auto it = std::find_if(rungs_.begin(), rungs_.end(), [&](const Rung &r) { return r.bits == bits; });
        if (it == rungs_.end()) {
            Rung rung;
            rung.bits = bits;
            const auto enclosed = source_.coefficient_enclosure(bits);
            const auto p = enclosed->polynomial.value();
            if (acb_poly_degree(p) != degree_ || acb_contains_zero(p->coeffs + degree_))
                return false;
            ComplexPolynomial monic;
            acb_poly_scalar_div(monic.value(), p, p->coeffs + degree_, bits);
            acb_poly_set_coeff_si(monic.value(), degree_, 1);
            if (source_.rectangular_available(bits)) {
                rung.blocks = std::make_unique<DiskBlockedPolynomialEvaluator>(std::move(monic), bits);
            } else if (source_.sparse_available(bits)) {
                rung.sparse = std::make_unique<DiskSparsePolynomialEvaluator>(monic, bits);
            }
            rungs_.push_back(std::move(rung));
            it = rungs_.end() - 1;
        }
        return it->blocks   ? it->blocks->evaluate(value, derivative, center)
               : it->sparse ? it->sparse->evaluate(value, derivative, center)
                            : source_.evaluate(value, derivative, center, bits);
    }
};

} // namespace certpolysolve::detail
