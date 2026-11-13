#pragma once
#include "methods/piecewise/piecewise_config.hpp"
#include "methods/piecewise/piecewise_statistics.hpp"
#include <span>
#include <algorithm>
#include <array>
#include "core/cover.hpp"
#include "input/coefficients/polynomial_input.hpp"
#include "core/query.hpp"
#include "methods/piecewise/piecewise_work.hpp"
#include "methods/piecewise/piecewise_estimate.hpp"
#include "numeric/models/root_model.hpp"
#include "certify/chart_certificate.hpp"
#include <functional>
#include <limits>
#include <optional>

namespace certpolysolve::detail {
struct PiecewiseUnfinishedOwner {
    uint64_t id = 0, candidate_id = UINT64_MAX, attempts = 0;
    Root proof;
    std::optional<ComplexBall> candidate;
    std::optional<PiecewiseCandidateEvent> last_attempt;
};
struct PiecewiseRootSnapshot {
    // Source is the exact factor used by this state, before any local chart.
    // Sharing this immutable owned copy avoids duplicating giant inputs.
    std::shared_ptr<const DensePolynomial> source;
    uint64_t source_id = 0, multiplicity = 1, bits = 0;
    PiecewiseRootStatistics statistics;
    std::vector<PiecewiseUnfinishedOwner> unfinished;
    bool terminal = false, singleton_ownership = false;
};
struct PiecewisePopulationPoint {
    ComplexBall point;
    std::optional<ComplexBall> indexed_point;
    uint64_t owner = UINT64_MAX, containing_owner = UINT64_MAX;
    bool inside_final_guard = false;
};
struct PiecewisePopulationSnapshot {
    uint64_t source = 0, round = 0, model_bits = 0, arithmetic_bits = 0;
    uint64_t discovered = 0, accurate = 0;
    bool terminal = false;
    std::vector<PiecewisePopulationPoint> points;
};
struct PiecewiseTransformCapture {
    PiecewiseTransformEvent event;
    std::vector<ComplexBall> input, output;
};

// Outward coverage proof for every radius in [inner,outer] and every angle.
// center_uncertainty bounds the distance from the exact center to every point
// in its stored ball. False requests the established angular cover.
bool piecewise_annulus_covered(const RealBall &inner, const RealBall &outer,
    const RealBall &gamma, const RealBall &radius, size_t sectors,
    const RealBall &center_uncertainty, slong precision);

// A ring's sector count before bit_ceil, the ceiling of 4 pi gamma / R, for the ring from 2^start to 2^next: at
// 96 bits from its edges (ring_sectors, what the construction uses), and in binary64 from its width in octaves where
// that provably gives the same count after bit_ceil (ring_sectors_binary64; nullopt near a power of two, or for a
// ring too narrow for its edges' accuracy). The admission's estimate uses the second, with the first as its fallback;
// tests compare them.
double ring_sectors(double start, double next);
std::optional<double> ring_sectors_binary64(double start, double next);

enum class PiecewiseModelStep { Work, BatchBoundary, Model, Done, Refused };
// Owns a bounded construction batch. model() is valid until the next advance;
// the cursor owns all backend storage and never retains a caller callback.
class PiecewiseModelCursor {
    struct Impl;
    std::unique_ptr<Impl> impl_;
  public:
    PiecewiseModelCursor(const acb_poly_t, slong, slong, PiecewiseRootStatistics &,
                         std::function<bool(const RootModel &)> = {}, PiecewiseRootConfig = {});
    ~PiecewiseModelCursor();
    PiecewiseModelStep advance(const cps_options &);
    const RootModel &model() const;
    // The suspended batch's construction work if every Taylor row is built
    // (construction_work's unit). Zero unless the cursor is at a batch boundary.
    double next_batch_work() const;
};

class PiecewiseRootState {
    struct Impl;
    std::unique_ptr<Impl> impl_;
  public:
    explicit PiecewiseRootState(const DensePolynomial &, PiecewiseRootConfig = {},
                                const std::vector<ComplexBall> *supplied = nullptr);
    ~PiecewiseRootState();
    cps_result advance(const cps_options &);
    const PiecewiseRootStatistics &statistics() const;
    // The regions the state has certified so far (what a stop short of the
    // goal completes with a counted disk; a cancelled call returns them as
    // they are).
    size_t certified_regions() const;
};

// Private experimental method. A callback streams one ring's owned models,
// avoiding a degree-by-model cache. It never calls the general root solver.
bool piecewise_root_models(const acb_poly_t polynomial, slong model_bits, slong precision, const cps_options &,
                           PiecewiseRootStatistics &, const std::function<bool(const RootModel &)> &,
                           const std::function<bool(const RootModel &)> &interest = {},
                           PiecewiseRootConfig config = {});
// Private replay entry: proposal work only, using the ordinary local polisher.
// The caller must separately certify every proposed result against its source.
void polish_piecewise_candidates(RootVector &, const RootModel &, const cps_options &,
                                 slong, PiecewiseRootStatistics &, PiecewiseRootConfig);
// Diagnostic replay of a fresh chart: ordinary machine proposals followed by
// the same two polishing entries and stability flags used in prepare/process.
// No carried candidates, ownership pruning, or certificate is supplied here.
// The caller owns a vector with exactly the model degree and certifies outputs.
slong propose_piecewise_candidates(RootVector &, const RootModel &, const cps_options &,
                                  slong, PiecewiseRootStatistics &, PiecewiseRootConfig);
cps_result solve_piecewise_roots(const DensePolynomial &, const cps_options &, PiecewiseRootStatistics &,
                                PiecewiseRootConfig config = {});
cps_result certify_piecewise_points(const DensePolynomial &, const cps_options &, PiecewiseRootStatistics &,
                                   const std::vector<ComplexBall> &points, PiecewiseRootConfig config = {});
} // namespace certpolysolve::detail
