#pragma once
#include "core/number.hpp"
#include <memory>
#include <stdexcept>
#include <vector>
#include <span>

namespace certpolysolve::detail {
// Internal enclosure contract: values contain p(z), p'(z) for EVERY point in z
// and every input allowed by the source coefficient enclosures. false means
// unresolved (including a secular pole), never a root-count conclusion.
// Outputs must not alias each other or z. Retained sources own their storage;
// backend pointers supplied to a call are borrowed only for that call.
class Evaluator {
  public:
    virtual ~Evaluator() = default;
    virtual slong degree() const = 0;
    virtual bool evaluate(acb_t value, acb_t derivative, const acb_t z, slong precision) const = 0;
    virtual bool evaluate_value(acb_t value, const acb_t z, slong precision) const {
        ComplexBall derivative;
        return evaluate(value, derivative.value(), z, precision);
    }
};

// A read-only coefficient snapshot shared by its prepared source and evaluators.
// Raw-input construction copies; owned-input construction transfers the value.
// Callers cannot mutate the storage or lend a handle from a temporary snapshot.
class PolynomialSnapshot {
    std::shared_ptr<const ComplexPolynomial> polynomial_;
    const ComplexPolynomial &owned_polynomial() const {
        if (!polynomial_) throw std::logic_error("access to a moved polynomial snapshot");
        return *polynomial_;
    }

  public:
    explicit PolynomialSnapshot(const acb_poly_t polynomial) {
        if (!polynomial) throw std::invalid_argument("null polynomial snapshot");
        auto owned = std::make_shared<ComplexPolynomial>();
        acb_poly_set(owned->value(), polynomial);
        polynomial_ = std::move(owned);
    }
    explicit PolynomialSnapshot(ComplexPolynomial polynomial)
        : polynomial_(std::make_shared<ComplexPolynomial>(std::move(polynomial))) {}
    const acb_poly_struct *value() const & {
        return owned_polynomial().value();
    }
    const acb_poly_struct *value() const && = delete;
    // A read-only borrow, valid while a snapshot retains this immutable
    // storage. The view itself does not increase the shared-owner count.
    ConstComplexBallView coefficients() const & { return owned_polynomial().coefficients(); }
    ConstComplexBallView coefficients() const && = delete;
};

class CoefficientEvaluator final : public Evaluator {
    PolynomialSnapshot polynomial_;
    std::vector<slong> indices_;
    bool sparse_;

  public:
    explicit CoefficientEvaluator(PolynomialSnapshot polynomial, bool sparse = true);
    explicit CoefficientEvaluator(const acb_poly_t polynomial, bool sparse = true)
        : CoefficientEvaluator(PolynomialSnapshot(polynomial), sparse) {}
    slong degree() const override { return acb_poly_degree(polynomial_.value()); }
    bool sparse() const { return sparse_; }
    bool evaluate(acb_t value, acb_t derivative, const acb_t z, slong precision) const override;
    bool evaluate_value(acb_t value, const acb_t z, slong precision) const override;
};

#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
// Reusable local Taylor polynomials with rigorous discarded-term bounds for
// BOTH p and p'. The source coefficients and ball arithmetic remain authoritative.
// This is an experiment, independent of the numerical root iteration policy: the
// research build's, for its Series test engine
// (research/solvers/local_polynomial_evaluator.cpp).
class LocalPolynomialEvaluator final : public Evaluator {
    PolynomialSnapshot polynomial_;
    CoefficientEvaluator fallback_;
    struct Local {
        ComplexBall center;
        ComplexPolynomial polynomial;
        slong precision = 0;
    };
    mutable std::vector<std::unique_ptr<Local>> cache_;
    mutable uint64_t retained_ = 0, discarded_ = 0, cache_hits_ = 0;

  public:
    explicit LocalPolynomialEvaluator(PolynomialSnapshot p)
        : polynomial_(std::move(p)), fallback_(polynomial_) {}
    explicit LocalPolynomialEvaluator(const acb_poly_t p)
        : LocalPolynomialEvaluator(PolynomialSnapshot(p)) {}
    slong degree() const override { return acb_poly_degree(polynomial_.value()); }
    bool evaluate(acb_t value, acb_t derivative, const acb_t z, slong precision) const override;
    uint64_t discarded() const { return discarded_; }
    uint64_t retained() const { return retained_; }
    uint64_t cache_hits() const { return cache_hits_; }
};
#endif

} // namespace certpolysolve::detail
