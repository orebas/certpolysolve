#pragma once
#include <cstddef>
#include <flint/flint.h>
#include <memory>
#include <vector>

struct cps_options;

namespace certpolysolve::detail {
struct DensePolynomial;
struct PiecewiseRootStatistics;
struct Root;
struct RootModel;

// One query's axis pass. The outer driver owns model construction, budgets and
// precision rungs. This component owns chart classification and root refinement.
// Input and statistics are borrowed from the stable outer state; models are copied
// on arrival. No caller result, model view or callback is retained.
class PiecewiseAxisPass {
    struct Impl;
    std::unique_ptr<Impl> impl_;

  public:
    PiecewiseAxisPass(const DensePolynomial &, size_t zeros, slong degree, PiecewiseRootStatistics &);
    ~PiecewiseAxisPass();
    PiecewiseAxisPass(const PiecewiseAxisPass &) = delete;
    PiecewiseAxisPass &operator=(const PiecewiseAxisPass &) = delete;
    PiecewiseAxisPass(PiecewiseAxisPass &&) = delete;
    PiecewiseAxisPass &operator=(PiecewiseAxisPass &&) = delete;

    // Call before the first arrival and before retrying a failed/unfinished round.
    // Discard its charts and pending roots; the source enclosure survives and is
    // rebuilt only for a different precision. A completed pass ends this query.
    void start_round(slong model_bits, slong precision);
    // Models must describe this input's nonzero part at the current round's
    // precisions. The outer cursor establishes their source/domain provenance.
    void arrive(const RootModel &, const cps_options &);
    // Called after the model stream ends or failed() becomes true. Append roots
    // only after the whole pass succeeds. Existing destination roots (including
    // the outer driver's exact zero roots) are preserved.
    void finish_round(const cps_options &, std::vector<Root> &destination);
    bool failed() const;
    bool complete() const;
};

} // namespace certpolysolve::detail
