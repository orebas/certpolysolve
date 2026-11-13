#pragma once
// Original polynomial evaluation contract, independent of coefficient caches
// and concrete coefficient, secular or program representations.
#include "core/raii.hpp"
#include <cstdint>

namespace certpolysolve::detail {
enum class InputSemantics { Exact, CoefficientEnclosures, Unknown };
enum class SourceRepresentation { Coefficients, SecularNumerator, Program };
enum class EvaluationStatus { Enclosed, Inconclusive, DenominatorContainsZero };
struct SourceMetadata {
    slong degree;
    bool real, monic;
    SourceRepresentation representation;
    InputSemantics input;
    bool coefficients_available;
};
struct ValueEvaluation {
    ComplexBall value;
    EvaluationStatus status = EvaluationStatus::Inconclusive;
    slong precision = 0;
};

// A monic POLYNOMIAL evaluator with known degree. A rational secular function
// cannot satisfy this contract merely by reporting its number of poles.
class MonicSource {
  public:
    virtual ~MonicSource() = default;
    virtual slong degree() const = 0;
    virtual bool real() const = 0;
    virtual SourceMetadata metadata() const = 0;
    virtual uint64_t terms() const { return static_cast<uint64_t>(degree()); }
    virtual bool evaluate(acb_t value, acb_t derivative, const acb_t z, slong precision) = 0;
    virtual bool evaluate_value(acb_t value, const acb_t z, slong precision) {
        ComplexBall derivative;
        return evaluate(value, derivative.value(), z, precision);
    }
    // An owning result cannot alias the input. Enclosed means a uniform ball
    // bound over the whole query and source uncertainty, not a promised
    // relative accuracy. Secular sources describe the polynomial numerator,
    // including limits at original poles.
    ValueEvaluation value(const ComplexBall &z, slong precision);
    // Optional full local coefficient model for p(center+x), preserving the
    // same monic normalization and all input uncertainty. false is fallback.
    virtual bool taylor(acb_poly_t, const acb_t, slong) { return false; }
};

} // namespace certpolysolve::detail
