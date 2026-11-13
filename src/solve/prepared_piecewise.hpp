#pragma once
#include "methods/piecewise/piecewise_config.hpp"
#include "methods/piecewise/piecewise_statistics.hpp"
#include "core/cover.hpp"
#include "core/query.hpp"
#include "input/coefficients/polynomial_input.hpp"
#include "numeric/exact/prepared_structure.hpp"
#include <memory>
#include <optional>

namespace certpolysolve::detail {
// End-to-end exact-input variant. Each nonzero square-free factor is solved by
// the same piecewise engine; no global numerical fallback is used.
class PreparedPiecewiseState {
    struct Impl;
    std::unique_ptr<Impl> impl_;
  public:
    explicit PreparedPiecewiseState(DensePolynomial,
                                    PiecewiseRootConfig = maintained_piecewise_config());
    ~PreparedPiecewiseState();
    cps_result advance(const cps_options &);
    const PiecewiseRootStatistics &statistics() const;
    std::optional<PreparedStructure> prepared_structure() const;
};
cps_result solve_prepared_piecewise_roots(const DensePolynomial &, const cps_options &,
                                        PiecewiseRootStatistics &,
                                        PiecewiseRootConfig config = maintained_piecewise_config());
} // namespace certpolysolve::detail
