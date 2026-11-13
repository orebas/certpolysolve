// Piecewise's local models: Taylor models of p with rigorous error bounds on charts covering the rings where its roots
// can lie, streamed by PiecewiseModelCursor; the admission's count of their work, unbuilt (estimate_piecewise_models);
// model-error bounds and restriction are shared operations in numeric/models/root_model.cpp.
//
// It follows Imbach and Moroz's piecewise approximation. CoefficientEnvelope is p's Newton polygon, the
// concave upper hull of integer upper bounds on log2 |a_j|. for_each_ring cuts the window of root moduli (root_window)
// into rings of at most an octave, each with the band l..u of terms within m = 16 + model bits of the largest, and
// none where one term dominates (no root lies there); estimate_on_envelope counts the walk's models and work, each
// ring's sector count in binary64 where that provably equals the 96-bit count (ring_sectors_binary64). A ring
// of mid-radius gamma has S charts of radius R = 3/4 of its width, centred at gamma times the S-th roots of unity
// (angular_grid), S the least power of two >= 4 pi gamma / R. Geometry is planned in binary64, bounds in outward
// balls. model_stream, a coroutine (ModelTask) yielding after each unit of work, builds in batches the charts its
// interest predicate accepts: row t of a batch is one length-S DFT (ModelTransform; axis_only: two sums) of the
// weights 2^-s a_(l+j) gamma^j C(j, t) rho^t, rho = R / gamma, summed by (j - t) mod S, 2^s bounding the band at
// gamma. So p(z) / (z^l 2^s) = F(w) + E(w) on z = c + R w, |w| <= 1, F the band's Taylor polynomial at c, |E| at most
// the envelope's bound() on the omitted terms (below l at gamma - R, above u at gamma + R) plus the tail's: Cauchy's
// M 2^-L (M the band's bound on |w| = 2, L = min(u - l, 4 m)) or, under adaptive_taylor, the geometric series of the
// first omitted row once it falls below 2^-m times the band at gamma - R, which ends the rows. On a rectangle,
// root_model_error takes, if smaller, the omitted terms at its radii (local_omitted_band_bound) and the tail T as
// T r^k where |w| <= r <= 1 on it, k the tail's first order.
//
// State: PiecewiseModelCursor::Impl copies p and the config (model_storage_bytes capped per advance by temporary_bytes)
// and owns the coroutine, whose frame holds a ring's grid, transform, scratch and a batch that model() views until the
// next advance. Under local_band_error the models share the envelope, kept past the cursor (RootModel::error_envelope).
// Entry points: PiecewiseModelCursor, made in start_round and stepped in advance by PiecewiseRootState
// (piecewise_roots.cpp), which hands a BatchBoundary's next_batch_work to begin_finishing; estimate_piecewise_models
// (admit_impl, admit_terms: piecewise_admission.cpp); root_model_error (certify_root_model and the axis mode in
// piecewise_roots.cpp, polish in piecewise_proposals.cpp); restricted_root_model (proposals, piecewise_proposals.cpp).
// Tests alone call piecewise_root_models; this file and tests call piecewise_annulus_covered, local_omitted_band_bound.
// Not here: model-error bounds and restriction (numeric/models/root_model.cpp), rounds, ownership, certificates
// (piecewise_roots.cpp), proposals and polishing (piecewise_proposals.cpp),
// admission rules (piecewise_admission.cpp), transform kernels (fixed_transform.cpp).
// Design notes: docs/piecewise-root-certification.md (models and their error bounds), docs/piecewise-selection.md
// (the estimate), docs/piecewise-local-band-derivation-20260922.md, docs/c1-axis-charts-note-20261015.md (axis mode).
#include "numeric/models/coefficient_envelope.hpp"
#include "methods/piecewise/piecewise_roots.hpp"
#include "core/dovetail.hpp"
#include "numeric/transforms/fixed_transform.hpp"
#include "platform/parallel.hpp"
#include "core/wide_int.hpp"
#include <span>
#include "core/telemetry.hpp"
#include <algorithm>
#include <bit>
#include <cmath>
#include <coroutine>
#include <exception>
#include <flint/acb_dft.h>
#include <limits>
#include <numbers>
#include <numeric>

namespace certpolysolve::detail {
namespace {
// W9's split of a Taylor row (model_stream, parallel.hpp): rows of at least row_split_min_terms terms, on a solve
// with workers, in row_split_blocks blocks of residue classes at most. Both constants fix the blocks by the row
// alone; shorter rows cost too little to share.
constexpr slong row_split_min_terms = 4096;
constexpr size_t row_split_blocks = 64;
} // namespace

bool piecewise_annulus_covered(const RealBall &inner, const RealBall &outer,
    const RealBall &gamma, const RealBall &radius, size_t sectors,
    const RealBall &center_uncertainty, slong precision) {
    if (!sectors || !arb_is_positive(inner.value()) || !arb_ge(outer.value(),inner.value()) ||
        !arb_is_positive(gamma.value()) || !arb_is_positive(radius.value()) ||
        !arb_is_nonnegative(center_uncertainty.value())) return false;
    RealBall angle, cosine, angular, radial, distance;
    arb_one(angle.value());
    arb_div_ui(angle.value(),angle.value(),static_cast<ulong>(sectors),precision);
    arb_cos_pi(cosine.value(),angle.value(),precision);
    arb_one(angular.value());
    arb_sub(angular.value(),angular.value(),cosine.value(),precision);
    arb_mul(angular.value(),angular.value(),gamma.value(),precision);
    arb_mul_2exp_si(angular.value(),angular.value(),1);
    // Distance squared is convex in the radial coordinate, so its maximum
    // over the entire annulus occurs at one of these two radial endpoints.
    for (const auto *endpoint : {&inner,&outer}) {
        arb_sub(radial.value(),endpoint->value(),gamma.value(),precision);
        arb_sqr(radial.value(),radial.value(),precision);
        arb_mul(distance.value(),angular.value(),endpoint->value(),precision);
        arb_add(distance.value(),distance.value(),radial.value(),precision);
        arb_sqrt(distance.value(),distance.value(),precision);
        arb_add(distance.value(),distance.value(),center_uncertainty.value(),precision);
        if (!arb_lt(distance.value(),radius.value())) return false;
    }
    return true;
}

namespace {
using Envelope = CoefficientEnvelope;

struct Ring {
    slong lower, upper;
    RealBall gamma, radius, inner, outer, proposal_radius;
    size_t sectors;
};
std::unique_ptr<RootVector> angular_grid(Ring &ring, slong precision,
    const PiecewiseRootConfig &config) {
    auto roots = std::make_unique<RootVector>(static_cast<slong>(ring.sectors));
    _acb_vec_unit_roots(roots->data(), roots->size(), roots->size(), precision);
    if (config.interior_proposals) {
        // The stored center ball may contain points besides its exact center.
        // Twice the sum of component radii bounds that whole ball's diameter.
        Float maximum, bound;
        RealBall real, imag, uncertainty;
        for (slong j = 0; j < roots->size(); ++j) {
            ComplexBall center;
            acb_mul_arb(center.value(), (*roots)[j], ring.gamma.value(), precision);
            arb_get_rad_arb(real.value(), acb_realref(center.value()));
            arb_get_rad_arb(imag.value(), acb_imagref(center.value()));
            arb_add(real.value(), real.value(), imag.value(), precision);
            arb_mul_2exp_si(real.value(), real.value(), 1);
            arb_get_ubound_arf(bound.value(), real.value(), precision);
            if (arf_cmp(bound.value(), maximum.value()) > 0) arf_set(maximum.value(), bound.value());
        }
        arb_set_arf(uncertainty.value(), maximum.value());
        // Prove the smallest binary64 radius fraction covering the annulus.
        // The stopping criterion is representability, not a tuned tolerance.
        const auto covered = [&](double fraction) {
            RealBall proposed;
            arb_set_d(proposed.value(), fraction);
            arb_mul(proposed.value(), proposed.value(), ring.radius.value(), precision);
            return piecewise_annulus_covered(ring.inner, ring.outer, ring.gamma, proposed,
                                             ring.sectors, uncertainty, precision);
        };
        if (covered(1)) {
            double lower = 0, upper = 1;
            for (;;) {
                const double middle = std::midpoint(lower, upper);
                if (middle == lower || middle == upper) break;
                if (covered(middle)) upper = middle;
                else lower = middle;
            }
            if (upper < 1 && covered(upper)) {
                // Split the available margin between proposal displacement
                // and the outer proof domain. Approximation roots may still
                // escape, so the controller retains full-chart recovery.
                arb_set_d(ring.proposal_radius.value(), upper);
                arb_add_ui(ring.proposal_radius.value(), ring.proposal_radius.value(), 1, precision);
                arb_mul_2exp_si(ring.proposal_radius.value(), ring.proposal_radius.value(), -1);
            }
        }
    }
    return roots;
}
void power_two(arb_t out, double exponent) {
    RealBall log_two;
    arb_set_d(out, exponent);
    arb_const_log2(log_two.value(), 96);
    arb_mul(out, out, log_two.value(), 96);
    arb_exp(out, out, 96);
    arb_get_mid_arb(out, out);
}
// A ring's edges 2^start and 2^next, mid-radius gamma and chart radius R = 3/4 of its width at 96 bits, and its
// sector count before bit_ceil, the ceiling of 4 pi gamma / R in binary64.
struct RingShape {
    RealBall inner, outer, gamma, radius;
    double sectors = 0;
};
RingShape ring_shape(double start, double next) {
    RingShape s;
    RealBall ratio;
    power_two(s.inner.value(), start);
    power_two(s.outer.value(), next);
    arb_add(s.gamma.value(), s.inner.value(), s.outer.value(), 96);
    arb_mul_2exp_si(s.gamma.value(), s.gamma.value(), -1);
    arb_get_mid_arb(s.gamma.value(), s.gamma.value());
    arb_sub(s.radius.value(), s.outer.value(), s.inner.value(), 96);
    arb_mul_ui(s.radius.value(), s.radius.value(), 3, 96);
    arb_mul_2exp_si(s.radius.value(), s.radius.value(), -2);
    arb_get_mid_arb(s.radius.value(), s.radius.value());
    arb_div(ratio.value(), s.gamma.value(), s.radius.value(), 96);
    s.sectors = std::ceil(4 * std::numbers::pi * arf_get_d(arb_midref(ratio.value()), ARF_RND_CEIL));
    return s;
}
// Walks the rings from the innermost outward. visit(Ring &&) returns false to
// stop early. The result is false when the geometry is unavailable.
//
// A ring spans at most one octave of radius. Where the model would hold a
// single term, that term is a vertex of the coefficient envelope and every
// other term is below it by the model's margin (2^-(bits + 16) each), so the
// polynomial has no root there and no ring is made: the rings follow the
// places where terms meet, not the width of the root window. (With a ring
// for every octave, x^3 + 2^e x - 1 would have about 1.5 e of them, all but
// a few with nothing to find, and at e = 2^30 too many to allocate.)
// History: docs/code-history-20261021.md, piecewise_models.cpp 1.
// With Geometry false (the admission's estimate, which reads a ring's band and sector count only), a ring's edges,
// gamma and R are left unset, and its sector count is ring_sectors_binary64's wherever that is proved equal.
template <bool Geometry, class Visit> bool for_each_ring(const Envelope &envelope, slong n, slong bits, Visit &&visit) {
    const auto h = envelope.negative_heights(n);
    const auto vertices = envelope.indices();
    const double m = static_cast<double>(bits) + 16;
    std::vector<double> low(static_cast<size_t>(n + 1), -std::numeric_limits<double>::infinity());
    std::vector<double> high(static_cast<size_t>(n + 1), std::numeric_limits<double>::infinity());
    for (slong j = 0; j <= n; ++j) {
        const size_t split =
            static_cast<size_t>(std::lower_bound(vertices.begin(), vertices.end(), j) - vertices.begin());
        if (split) {
            auto slope = [&](size_t k) {
                const auto i = vertices[k];
                return (h[static_cast<size_t>(j)] - m - h[static_cast<size_t>(i)]) /
                       static_cast<double>(j - i);
            };
            size_t a = 0, b = split - 1;
            while (a < b) {
                const size_t k = a + (b - a) / 2;
                if (slope(k) < slope(k + 1))
                    a = k + 1;
                else
                    b = k;
            }
            low[static_cast<size_t>(j)] = slope(a);
        }
        const size_t right =
            static_cast<size_t>(std::upper_bound(vertices.begin(), vertices.end(), j) - vertices.begin());
        if (right < vertices.size()) {
            auto slope = [&](size_t k) {
                const auto i = vertices[k];
                return (h[static_cast<size_t>(i)] - h[static_cast<size_t>(j)] + m) /
                       static_cast<double>(i - j);
            };
            size_t a = right, b = vertices.size() - 1;
            while (a < b) {
                const size_t k = a + (b - a) / 2;
                if (slope(k) > slope(k + 1))
                    a = k + 1;
                else
                    b = k;
            }
            high[static_cast<size_t>(j)] = slope(a);
        }
    }
    double start = low[1], end = high[static_cast<size_t>(n - 1)];
    const auto [inner, outer] = envelope.root_window(h);
    start = std::max(start, inner);
    end = std::min(end, outer);
    slong lower = 0, upper = 0;
    while (start < end) {
        while (lower < n && high[static_cast<size_t>(lower)] <= start)
            ++lower;
        upper = std::max(lower, upper);
        while (upper < n &&
               static_cast<double>(upper + 1 - lower) * (low[static_cast<size_t>(upper + 1)] - start) < m)
            ++upper;
        if (upper == lower && std::binary_search(vertices.begin(), vertices.end(), lower)) {
            // A stretch with a single term: go on, in the same one-octave
            // steps, to the first at which the next term is within reach (the
            // test above, which is strict), or past the window's end.
            const double reach = upper < n ? low[static_cast<size_t>(upper + 1)] - m : end;
            const double after = start + std::max(1.0, std::floor(std::min(reach, end) - start) + 1.0);
            if (!(after > start) || !std::isfinite(after))
                return false;
            start = after;
            continue;
        }
        const double next =
            std::min({end, start + 1.0,
                      upper > lower ? start + m / static_cast<double>(upper - lower + 1)
                                    : (upper < n ? low[static_cast<size_t>(upper + 1)] : end)});
        if (!(next > start) || !std::isfinite(next))
            return false;
        Ring ring;
        ring.lower = lower;
        ring.upper = upper;
        std::optional<double> estimated;
        if constexpr (!Geometry)
            estimated = ring_sectors_binary64(start, next);
        double sectors;
        if (estimated) {
            sectors = *estimated;
        } else {
            auto shape = ring_shape(start, next);
            ring.inner = std::move(shape.inner);
            ring.outer = std::move(shape.outer);
            ring.gamma = std::move(shape.gamma);
            ring.radius = std::move(shape.radius);
            sectors = shape.sectors;
        }
        if (!std::isfinite(sectors) || sectors > static_cast<double>(std::bit_floor(static_cast<size_t>(PTRDIFF_MAX) / sizeof(acb_struct))))
            return false;
        ring.sectors = std::bit_ceil(static_cast<size_t>(std::max(1.0, sectors)));
        if (!visit(std::move(ring)))
            return true;
        start = next;
    }
    return true;
}
// Every ring, or none when the geometry is unavailable.
std::vector<Ring> rings(const Envelope &envelope, slong n, slong bits) {
    std::vector<Ring> result;
    if (!for_each_ring<true>(envelope, n, bits, [&](Ring &&ring) {
            result.push_back(std::move(ring));
            return true;
        }))
        result.clear();
    return result;
}
} // namespace

namespace {
// The estimate on an envelope of degree n, the ends' balls tested (`ends_contain_zero`): shared by the coefficient
// polynomial and the terms.
PiecewiseModelEstimate estimate_on_envelope(const Envelope &envelope, slong n, slong bits,
                                            std::optional<uint64_t> maximum_envelope_bend,
                                            const std::function<bool(uint64_t)> &beyond);
} // namespace

double ring_sectors(double start, double next) { return ring_shape(start, next).sectors; }

std::optional<double> ring_sectors_binary64(double start, double next) {
    // The count after bit_ceil, bit_ceil(ceil(x)) for x = 4 pi gamma / R, changes only where x crosses a power of two.
    // x = (8 pi / 3) (2^w + 1) / (2^w - 1) for the width w = next - start in octaves. Both computations stay within
    // 2^-47 of x, relative.
    // - ring_shape's: an edge 2^e at 96 bits is within about 2^-94 (1 + |e|) of its value, relative (the exponent
    //   e ln 2 is formed at 96 bits), so outer - inner is within 2^-93 (1 + max |e|) / (w ln 2) of the width: at most
    //   2^-52.5 where w >= sector_margin (1 + max |e|). Then one rounding up to binary64, pi in binary64, one product.
    // - This one: the difference next - start within an ulp of w, expm1 within a few ulps, and five roundings.
    // So where x is farther than sector_margin from every power of two, both lie on the same side of each, and the
    // counts agree. Elsewhere, and for narrower rings, nullopt: the caller takes ring_shape's count.
    constexpr double sector_margin = 0x1p-40;
    const double width = next - start;
    if (!std::isfinite(width) || !(width >= sector_margin * (1 + std::max(std::abs(start), std::abs(next)))))
        return std::nullopt;
    const double x = 4 * std::numbers::pi * (2.0 / 3.0) * (1 + 2 / std::expm1(width * std::numbers::ln2));
    if (!std::isfinite(x))
        return std::nullopt;
    const double power = std::exp2(std::round(std::log2(x)));
    if (std::abs(x - power) <= sector_margin * power)
        return std::nullopt;
    return std::ceil(x);
}

PiecewiseModelEstimate estimate_piecewise_models(const acb_poly_t p, slong bits,
                                                std::optional<uint64_t> maximum_envelope_bend,
                                                const std::function<bool(uint64_t)> &beyond) {
    const slong n = acb_poly_degree(p);
    if (n < 1 || bits < 2 || static_cast<uint64_t>(bits) > resources::precision_max ||
        acb_contains_zero(p->coeffs) || acb_contains_zero(p->coeffs + n))
        return {};
    return estimate_on_envelope(Envelope(p), n, bits, maximum_envelope_bend, beyond);
}

PiecewiseModelEstimate estimate_piecewise_models(std::span<const std::pair<slong, ComplexBall>> terms, slong bits,
                                                std::optional<uint64_t> maximum_envelope_bend,
                                                const std::function<bool(uint64_t)> &beyond) {
    if (terms.empty()) return {};
    const slong n = terms.back().first;
    if (n < 1 || bits < 2 || static_cast<uint64_t>(bits) > resources::precision_max ||
        terms.front().first != 0 || acb_contains_zero(terms.front().second.value()) ||
        acb_contains_zero(terms.back().second.value()))
        return {};
    return estimate_on_envelope(Envelope(terms), n, bits, maximum_envelope_bend, beyond);
}

namespace {
PiecewiseModelEstimate estimate_on_envelope(const Envelope &envelope, slong n, slong bits,
                                            std::optional<uint64_t> maximum_envelope_bend,
                                            const std::function<bool(uint64_t)> &beyond) {
    PiecewiseModelEstimate result;
    const auto bend = envelope.bend();
    result.envelope_bend = static_cast<double>(bend.excess) / static_cast<double>(bend.span);
    // Admission must reject unsuitable geometry before allocating its rings.
    // An enormous radial gap can otherwise make even a cubic's probe costly.
    // The cap is compared exactly: excess / span > cap, that is excess > cap x span, below 2^127 (a cap below 2^64,
    // a span below 2^63).
    if (maximum_envelope_bend && bend.excess > static_cast<int128>(*maximum_envelope_bend) * bend.span) {
        result.bend_exceeded = true;
        return result;
    }
    // The rings are counted, not stored: their number follows the radial
    // width of the root window, not the degree.
    bool overflow = false;
    const bool walked = for_each_ring<false>(envelope, n, bits, [&](Ring &&ring) {
        const auto delta = ring.upper - ring.lower;
        const auto degree = std::min(delta, resources::precision_product(resources::precision_sum(bits, 16), 4));
        const uint64_t terms = static_cast<uint64_t>(degree + 1);
        if (ring.sectors > (std::numeric_limits<uint64_t>::max() - result.coefficients) / terms ||
            ring.sectors > std::numeric_limits<uint64_t>::max() - result.models) {
            overflow = true;
            return false;
        }
        ++result.rings;
        result.models += ring.sectors;
        result.coefficients += ring.sectors * terms;
        result.maximum_degree = std::max(result.maximum_degree, static_cast<uint64_t>(degree));
        // Counts, exact; saturating, as a caller's bound stops the walk far below 2^64.
        const uint64_t sectors = ring.sectors;
        const uint64_t levels = sectors ? static_cast<uint64_t>(std::bit_width(sectors) - 1) : 0;
        result.construction_work = resources::add(
            result.construction_work,
            resources::multiply(terms, resources::add(resources::nonnegative(delta + 1), resources::multiply(sectors, levels))));
        result.local_work = resources::add(result.local_work, resources::multiply(resources::multiply(sectors, terms), terms));
        // The work only grows from here, so a caller whose decision is settled
        // once it passes a bound need not see the rest.
        if (beyond && beyond(resources::add(result.construction_work, result.local_work))) {
            result.truncated = true;
            return false;
        }
        return true;
    });
    if (overflow)
        return {};
    if (!walked) {
        PiecewiseModelEstimate unavailable;
        unavailable.envelope_bend = result.envelope_bend;
        return unavailable;
    }
    result.available = result.rings != 0;
    return result;
}
} // namespace

namespace {
struct ModelYield {
    PiecewiseModelStep step = PiecewiseModelStep::Work;
    const RootModel *model = nullptr; // Checked view into the suspended frame.
    double next_batch_work = 0;
};
class ModelTransform {
    struct Plan {
        acb_dft_rad2_t value;
        Plan(int exponent, slong precision) { acb_dft_rad2_init(value, exponent, precision); }
        ~Plan() { acb_dft_rad2_clear(value); }
        Plan(const Plan &) = delete;
    } plan_;
    slong precision_;
    std::unique_ptr<FixedTransform> fixed_;
  public:
    ModelTransform(int exponent, slong precision, bool fixed)
      : plan_(exponent, precision), precision_(precision) {
        if (fixed && FixedTransform::limbs_for(precision, exponent)) fixed_ = std::make_unique<FixedTransform>(exponent);
    }
    ~ModelTransform() = default;
    ModelTransform(const ModelTransform &) = delete;
    ModelTransform &operator=(const ModelTransform &) = delete;
    // `wanted`: the outputs the caller will read (all when empty). The fixed-point transform publishes those
    // alone; FLINT's ball transform fills every output.
    void apply(RootVector &out, const RootVector &in, PiecewiseRootStatistics &stats,
               std::span<const slong> wanted = {}) {
        // A radix-2 transform of length S has (S / 2) log2 S butterflies.
        const double points = static_cast<double>(out.size());
        const double levels = static_cast<double>(std::bit_width(static_cast<uint64_t>(out.size())) - 1);
        // The ledger charges a transform at the working precision the same whichever arithmetic delivers it: the
        // regime exit and the finisher's admission were qualified with this charge, and pricing the fixed-point
        // transform at its own cost is a separate change to those two rules.
        stats.work.butterflies += points / 2 * levels * operation_weight(precision_);
        if (fixed_) {
            FixedTransformStatistics s;
            const bool transformed = fixed_->apply_into(out, in, precision_, wanted, s);
            stats.fft_fixed_accepted += s.accepted;
            stats.fft_fixed_declined += s.declined;
            stats.fft_fixed_ns += s.conversion_ns + s.execution_ns + s.enclosure_ns;
            if (transformed) return;
        }
        acb_dft_rad2_precomp(out.data(), in.data(), plan_.value, precision_);
    }
};
// The frame owns rings, FFT scratch and one batch. Destruction releases every
// RAII value even when a caller cancels or abandons a partially built model.
class ModelTask {
  public:
    struct promise_type {
        ModelYield current;
        std::exception_ptr exception;
        ModelTask get_return_object() { return ModelTask(Handle::from_promise(*this)); }
        std::suspend_always initial_suspend() noexcept { return {}; }
        std::suspend_always final_suspend() noexcept { return {}; }
        std::suspend_always yield_value(ModelYield value) noexcept { current = value; return {}; }
        void return_void() noexcept {}
        void unhandled_exception() { exception = std::current_exception(); }
    };
    using Handle = std::coroutine_handle<promise_type>;
    explicit ModelTask(Handle handle) : handle_(handle) {}
    ModelTask(const ModelTask &) = delete;
    ~ModelTask() { if (handle_) handle_.destroy(); }
    ModelYield advance() {
        if (!handle_.done()) handle_.resume();
        if (handle_.promise().exception) std::rethrow_exception(handle_.promise().exception);
        return handle_.done() ? ModelYield{PiecewiseModelStep::Done, nullptr} : handle_.promise().current;
    }
  private:
    Handle handle_;
};
ModelTask model_stream(const acb_poly_t p, slong bits, slong precision,
                       PiecewiseRootStatistics &stats,
                       const std::function<bool(const RootModel &)> &interest,
                       const PiecewiseRootConfig &config) {
    const slong n = acb_poly_degree(p);
    if (n < 1 || bits < 2 || static_cast<uint64_t>(bits) > resources::precision_max ||
        acb_contains_zero(p->coeffs) || acb_contains_zero(p->coeffs + n)) {
        co_yield {PiecewiseModelStep::Refused};
        co_return;
    }
    std::optional<Envelope> built;
    try {
        StageTimer timer(stats.setup_ns);
        built.emplace(p);
    } catch (const std::length_error &) {
        // Exponents beyond the envelope's working range: refuse the input, as
        // for the other unsupportable inputs above. It is not an internal
        // error, and the state's owner can hand the solve to another method.
        ++stats.range_refusals;
    }
    if (!built) {
        co_yield {PiecewiseModelStep::Refused};
        co_return;
    }
    const Envelope &envelope = *built;
    const auto shared_envelope = config.local_band_error ? std::make_shared<const Envelope>(envelope) : nullptr;
    std::vector<Ring> domains;
    { StageTimer timer(stats.setup_ns); domains = rings(envelope, n, bits); }
    co_yield {PiecewiseModelStep::Work};
    for (auto ring : domains) {
        const slong delta = ring.upper - ring.lower, L = std::min(delta, resources::precision_product(resources::precision_sum(bits, 16), 4));
        const size_t bytes = sizeof(acb_struct) +
            4 * static_cast<size_t>(resources::precision_sum(precision, FLINT_BITS - 1) / FLINT_BITS) * sizeof(ulong);
        // Check sizes before allocating even the root-of-unity grid. Include
        // FFT buffers, weights, centers, metadata and the polynomial matrix.
        // Exact and saturating, as per_model below.
        const uint64_t scratch_estimate = resources::multiply(
            bytes, resources::add(resources::multiply(4, ring.sectors), resources::nonnegative(delta + 1)));
        const uint64_t per_model = resources::add(resources::multiply(bytes, resources::add(resources::nonnegative(L), 5)),
                                                       sizeof(RootModel) + sizeof(size_t));
        if (scratch_estimate >= static_cast<uint64_t>(PTRDIFF_MAX) ||
            !resources::storage_fits(per_model, config.model_storage_bytes) ||
            scratch_estimate > config.model_storage_bytes ||
            config.model_storage_bytes - scratch_estimate < per_model) {
            ++stats.storage_refusals;
            co_yield {PiecewiseModelStep::Refused};
            co_return;
        }
        std::unique_ptr<RootVector> roots;
        {
            StageTimer timer(stats.construction_ns);
            OptionalStageTimer detail(stats.fft_setup_ns,config.detailed_accounting);
            roots=angular_grid(ring,precision,config);
            stats.work.grid += static_cast<double>(ring.sectors) * operation_weight(precision);
        }
        const size_t capacity = std::min(ring.sectors,
            (config.model_storage_bytes - static_cast<size_t>(scratch_estimate)) / per_model);
        const slong count = static_cast<slong>(ring.sectors);
        RootVector weights(delta + 1), groups(count), values(count);
        std::unique_ptr<ModelTransform> transform;
        RealBall rho, power, norm, lower_radius, upper_radius, outer, error, tail, part, band_error, tail_allowance;
        slong scale;
        {
            StageTimer timer(stats.construction_ns);
            {
                OptionalStageTimer detail(stats.fft_setup_ns, config.detailed_accounting);
                transform = std::make_unique<ModelTransform>(static_cast<int>(std::bit_width(ring.sectors) - 1),
                                                            precision, config.fixed_fft);
                stats.work.grid += static_cast<double>(ring.sectors) / 2 * operation_weight(precision);
            }
            OptionalStageTimer detail(stats.coefficient_preparation_ns, config.detailed_accounting);
            arb_div(rho.value(), ring.radius.value(), ring.gamma.value(), precision);
            envelope.bound(norm.value(), ring.lower, ring.upper, ring.gamma.value(), ring.lower, 0);
            Float magnitude;
            arb_get_ubound_arf(magnitude.value(), norm.value(), 64);
            scale = arf_abs_bound_lt_2exp_si(magnitude.value());
            arb_sub(lower_radius.value(), ring.gamma.value(), ring.radius.value(), precision);
            arb_add(upper_radius.value(), ring.gamma.value(), ring.radius.value(), precision);
            envelope.bound(error.value(), 0, ring.lower - 1, lower_radius.value(), ring.lower, scale);
            envelope.bound(part.value(), ring.upper + 1, n, upper_radius.value(), ring.lower, scale);
            arb_add(error.value(), error.value(), part.value(), precision);
            band_error = error;
            if (L < delta) {
                arb_mul_2exp_si(outer.value(), ring.radius.value(), 1);
                arb_add(outer.value(), outer.value(), ring.gamma.value(), precision);
                envelope.bound(tail.value(), ring.lower, ring.upper, outer.value(), ring.lower, scale);
                arb_mul_2exp_si(tail.value(), tail.value(), -L);
                arb_add(error.value(), error.value(), tail.value(), precision);
            }
            tail_allowance = tail;
            if (config.adaptive_taylor) {
                envelope.bound(tail_allowance.value(), ring.lower, ring.upper, lower_radius.value(), ring.lower, scale);
                arb_mul_2exp_si(tail_allowance.value(), tail_allowance.value(),
                               -resources::precision_sum(bits, 16));
            }
            ++stats.rings;
        }
        co_yield {PiecewiseModelStep::Work};
        size_t sector = 0;
        while (sector < ring.sectors) {
            std::vector<RootModel> models;
            std::vector<size_t> active;
            models.reserve(capacity);
            active.reserve(capacity);
            const size_t first_sector = sector;
            {
                StageTimer timer(stats.construction_ns);
                {
                OptionalStageTimer detail(stats.model_assembly_ns, config.detailed_accounting);
                while (sector < ring.sectors && models.size() < capacity) {
                    RootModel model;
                    model.proposal_radius = ring.proposal_radius;
                    acb_mul_arb(model.center.value(), (*roots)[static_cast<slong>(sector)], ring.gamma.value(), precision);
                    model.radius = ring.radius;
                    if (!interest || interest(model)) {
                        model.error = error;
                        if (config.error_components || config.local_band_error || config.adaptive_taylor) {
                            model.omitted_band_error = band_error;
                            model.taylor_model_error = tail;
                            model.error_components = config.error_components;
                        }
                        if (config.adaptive_taylor) model.taylor_error_order = L + 1;
                        model.error_envelope = shared_envelope;
                        model.lower_index = ring.lower;
                        model.upper_index = ring.upper;
                        model.scale = scale;
                        // FLINT otherwise doubles capacity as Taylor rows are
                        // appended, exceeding the matrix storage estimate.
                        acb_poly_fit_length(model.polynomial.value(), L + 1);
                        active.push_back(sector);
                        models.push_back(std::move(model));
                    } else ++stats.skipped_models;
                    ++sector;
                }
                }
                if (!models.empty()) {
                    OptionalStageTimer detail(stats.coefficient_preparation_ns, config.detailed_accounting);
                    arb_one(power.value());
                    arb_mul_2exp_si(power.value(), power.value(), -scale);
                    for (slong j = 0; j <= delta; ++j) {
                        acb_mul_arb(weights[j], p->coeffs + ring.lower + j, power.value(), precision);
                        arb_mul(power.value(), power.value(), ring.gamma.value(), precision);
                    }
                    ++stats.batches;
                    stats.work.weights += 2 * static_cast<double>(delta + 1) * operation_weight(precision);
                }
                stats.work.weights += static_cast<double>(sector - first_sector) * operation_weight(precision);
            }
            co_yield {PiecewiseModelStep::Work};
            if (models.empty()) continue;
            // Charged recurrence terms, FFT butterflies and coefficient writes.
            // Calibrate their time from completed construction work; this is
            // scheduling metadata and never a certification estimate.
            const uint64_t row_units = resources::add(static_cast<uint64_t>(delta) + 1,
                resources::add(resources::multiply(ring.sectors, resources::nonnegative(std::bit_width(ring.sectors))), models.size()));
            if (config.source_finisher) {
                // The rows t = 0..L in counted work: sum_t (4 (delta - t) + 1)
                // recurrence terms and one transform each.
                const double rows = static_cast<double>(L) + 1;
                const double recurrence_terms =
                    rows * (4 * static_cast<double>(delta) + 1) - 2 * static_cast<double>(L) * rows;
                const double points = static_cast<double>(ring.sectors);
                const double levels = static_cast<double>(std::bit_width(ring.sectors) - 1);
                const auto &f = piecewise_work_factors;
                const double transform = f.butterflies * points / 2 * levels * operation_weight(precision);
                co_yield {PiecewiseModelStep::BatchBoundary, nullptr,
                          f.recurrence * recurrence_terms * operation_weight(precision) + rows * transform};
            }
            // The outputs the batch's models read: one a sector of interest (a transform capture reads them all).
            std::vector<slong> wanted;
            if (!config.capture_transform) {
                wanted.reserve(active.size());
                for (const size_t sector : active) wanted.push_back((count - static_cast<slong>(sector)) % count);
            }
            // C1: only the axis charts are active, whose rows need outputs 0 and S/2 alone.
            const bool axis_outputs = config.axis_only && !config.capture_transform &&
                std::all_of(active.begin(), active.end(), [&](size_t sector) {
                    return sector == 0 || static_cast<slong>(sector) == count / 2;
                });
            slong last_row = L;
            std::vector<Magnitude> row_magnitudes; // a split row's term magnitudes, one a j (below)
            for (slong t = 0; t <= L; ++t) {
                const uint64_t fft_before = stats.fft_execution_ns;
                bool tail_proved = false;
                {
                    StageTimer timer(stats.construction_ns);
                    {
                    OptionalStageTimer detail(stats.taylor_rows_ns, config.detailed_accounting);
                    RealBall ratio, denominator, proved_tail;
                    Magnitude first_omitted, term;
                    bool geometric = false;
                    if (config.adaptive_taylor && t < L) {
                        // After this row, weights[j] encloses Taylor term t+1.
                        // Every subsequent term ratio is at most
                        // rho * (delta-t-1)/(t+2), and decreases. Here rho is
                        // already the chart radius divided by gamma.
                        arb_mul_si(ratio.value(), rho.value(), delta - t - 1, precision);
                        arb_div_si(ratio.value(), ratio.value(), t + 2, precision);
                        arb_one(denominator.value());
                        arb_sub(denominator.value(), denominator.value(), ratio.value(), precision);
                        geometric = arb_is_positive(denominator.value());
                    }
                    _acb_vec_zero(groups.data(), count);
                    // Term t + 1 from term t: weights[j] rho (j - t) / (t + 1). The real factor is formed first,
                    // rho / (t + 1) once a row and a small-integer product a term, so a weight takes one real
                    // scaling. History: docs/code-history-20261021.md, piecewise_models.cpp 2.
                    RealBall row_step, factor;
                    arb_div_si(row_step.value(), rho.value(), t + 1, precision);
                    if (parallel_width() < 2 || delta - t + 1 < row_split_min_terms || count < 2) {
                        for (slong j = t; j <= delta; ++j) {
                            acb_add(groups[(j - t) % count], groups[(j - t) % count], weights[j], precision);
                            if (j > t) {
                                arb_mul_si(factor.value(), row_step.value(), j - t, precision);
                                acb_mul_arb(weights[j], weights[j], factor.value(), precision);
                                if (geometric) {
                                    acb_get_mag(term.value(), weights[j]);
                                    mag_add(first_omitted.value(), first_omitted.value(), term.value());
                                }
                            }
                        }
                    } else {
                        // W9's split (parallel.hpp): by residue class, in fixed blocks of classes, each class's terms
                        // in increasing j, so groups[c] receives its weights as the loop above adds them. A term's
                        // magnitude goes to its own slot, and the tail sum adds them in increasing j after the join,
                        // as the running mag_add above does: the same bits at every worker count.
                        // The blocks take copies: the data of groups and weights (made with count and delta + 1
                        // entries) and the row's scalars. Captured by reference, model_stream's locals would have
                        // their addresses taken, and the serial loop above would reload them after every FLINT
                        // call: at one worker, 206 million more instructions in model_stream on
                        // massive-dense_uniform-32768 (the W9 gate, 24 October 2026).
                        if (groups.size() < count || weights.size() <= delta)
                            throw std::logic_error("model_stream: a Taylor row's split outgrows its vectors");
                        if (geometric) row_magnitudes.resize(static_cast<size_t>(delta) + 1);
                        const auto classes = static_cast<size_t>(count);
                        const size_t blocks = std::min(classes, row_split_blocks);
                        parallel_for(SplitSite::TaylorRow, blocks,
                                     [group = groups.data(), weight = weights.data(), step_source = row_step.value(),
                                      magnitude = row_magnitudes.data(), t, delta, count, precision, geometric,
                                      classes, blocks](size_t b) {
                            RealBall step;
                            for (size_t c = b * classes / blocks; c < (b + 1) * classes / blocks; ++c)
                                for (slong j = t + static_cast<slong>(c); j <= delta; j += count) {
                                    acb_add(group + c, group + c, weight + j, precision);
                                    if (j == t) continue;
                                    arb_mul_si(step.value(), step_source, j - t, precision);
                                    acb_mul_arb(weight + j, weight + j, step.value(), precision);
                                    if (geometric) acb_get_mag(magnitude[static_cast<size_t>(j)].value(), weight + j);
                                }
                        });
                        if (geometric)
                            for (slong j = t + 1; j <= delta; ++j)
                                mag_add(first_omitted.value(), first_omitted.value(),
                                        row_magnitudes[static_cast<size_t>(j)].value());
                    }
                    if (geometric) {
                        Float bound;
                        arf_set_mag(bound.value(), first_omitted.value());
                        arb_set_arf(proved_tail.value(), bound.value());
                        arb_div(proved_tail.value(), proved_tail.value(), denominator.value(), precision);
                        arb_get_ubound_arf(bound.value(), proved_tail.value(), precision);
                        arb_set_arf(proved_tail.value(), bound.value());
                        tail_proved = arb_le(proved_tail.value(), tail_allowance.value());
                        if (tail_proved) {
                            last_row = t;
                            for (auto &model : models) {
                                model.taylor_model_error = proved_tail;
                                model.taylor_error_order = t + 1;
                                arb_add(model.error.value(), band_error.value(), proved_tail.value(), precision);
                            }
                        }
                    }
                    }
                    {
                        OptionalStageTimer detail(stats.fft_execution_ns, config.detailed_accounting);
                        if (axis_outputs) {
                            // The transform's outputs 0 and S/2 (X_k = sum_j x_j w^(jk), w = e^(-2 pi i/S),
                            // src/fixed_transform.hpp): the sum and the alternating sum, in ball arithmetic.
                            acb_zero(values[0]);
                            acb_zero(values[count / 2]);
                            for (slong j = 0; j < count; ++j) {
                                acb_add(values[0], values[0], groups[j], precision);
                                if (j & 1) acb_sub(values[count / 2], values[count / 2], groups[j], precision);
                                else acb_add(values[count / 2], values[count / 2], groups[j], precision);
                            }
                        } else
                            transform->apply(values, groups, stats, wanted);
                    }
                    {
                    OptionalStageTimer detail(stats.model_assembly_ns, config.detailed_accounting);
                    // Row t of the transform is coefficient t of the models, one output a sector. An output goes
                    // to one model only and the next row's transform overwrites them all, so it is moved into its
                    // slot (reserved above, never written before: zero) where acb_poly_set_coeff_acb would copy
                    // it: the same polynomial, without two allocations a coefficient. A transform capture reads
                    // the outputs afterwards and keeps the copy.
                    const bool moved = !config.capture_transform;
                    for (size_t k = 0; k < active.size(); ++k) {
                        const acb_ptr output = values[(count - static_cast<slong>(active[k])) % count];
                        acb_poly_struct *polynomial = models[k].polynomial.value();
                        if (!moved) {
                            acb_poly_set_coeff_acb(polynomial, t, output);
                            continue;
                        }
                        acb_swap(polynomial->coeffs + t, output);
                        _acb_poly_set_length(polynomial, t + 1);
                        _acb_poly_normalise(polynomial);
                    }
                    }
                    ++stats.fft_transforms;
                    stats.construction_units = resources::add(stats.construction_units, row_units);
                    // One addition per term and three scalings per term after the first: the charge the ledger's
                    // decisions were qualified with (a term now takes one scaling and a small product).
                    stats.work.recurrence += (4 * static_cast<double>(delta - t) + 1) * operation_weight(precision);
                }
                if (config.observe_transform || config.capture_transform) {
                    StageTimer diagnostic(stats.diagnostic_ns);
                    const PiecewiseTransformEvent event{static_cast<uint64_t>(count),
                        static_cast<uint64_t>(bits), static_cast<uint64_t>(precision),
                        static_cast<uint64_t>(t), stats.fft_execution_ns - fft_before};
                    if (config.observe_transform) config.observe_transform(event);
                    if (config.capture_transform && config.want_transform && config.want_transform(event)) {
                        PiecewiseTransformCapture saved;
                        saved.event = event;
                        saved.input.resize(static_cast<size_t>(count));
                        saved.output.resize(static_cast<size_t>(count));
                        for (slong j = 0; j < count; ++j) {
                            acb_set(saved.input[static_cast<size_t>(j)].value(), groups[j]);
                            acb_set(saved.output[static_cast<size_t>(j)].value(), values[j]);
                        }
                        config.capture_transform(std::move(saved));
                    }
                }
                co_yield {PiecewiseModelStep::Work};
                if (tail_proved) break;
            }
            stats.models += models.size();
            stats.coefficients += models.size() * static_cast<size_t>(last_row + 1);
            stats.maximum_degree = std::max(stats.maximum_degree, static_cast<uint64_t>(last_row));
            for (const auto &model : models) co_yield {PiecewiseModelStep::Model, &model};
        }
    }
}
} // namespace
struct PiecewiseModelCursor::Impl {
    ComplexPolynomial polynomial;
    PiecewiseRootStatistics &statistics;
    std::function<bool(const RootModel &)> interest;
    PiecewiseRootConfig config;
    ModelTask task;
    ModelYield current;
    Impl(const acb_poly_t p, slong bits, slong precision, PiecewiseRootStatistics &s,
         std::function<bool(const RootModel &)> wanted, PiecewiseRootConfig config)
        : statistics(s), interest(std::move(wanted)), config(config),
          task(model_stream(polynomial.value(), bits, precision, statistics, interest, this->config)) {
        acb_poly_set(polynomial.value(), p);
    }
};
PiecewiseModelCursor::PiecewiseModelCursor(const acb_poly_t p, slong bits, slong precision,
    PiecewiseRootStatistics &s, std::function<bool(const RootModel &)> interest, PiecewiseRootConfig config)
    : impl_(std::make_unique<Impl>(p, bits, precision, s, std::move(interest), config)) {}
PiecewiseModelCursor::~PiecewiseModelCursor() = default;
PiecewiseModelStep PiecewiseModelCursor::advance(const cps_options &options) {
    if (options.cancelled()) return PiecewiseModelStep::Work;
    impl_->config.model_storage_bytes = static_cast<size_t>(std::min<uint64_t>(
        impl_->config.model_storage_bytes, options.temporary_bytes));
    impl_->current = impl_->task.advance();
    return impl_->current.step;
}
const RootModel &PiecewiseModelCursor::model() const {
    if (impl_->current.step != PiecewiseModelStep::Model || !impl_->current.model)
        throw std::logic_error("piecewise cursor has no current model");
    return *impl_->current.model;
}
double PiecewiseModelCursor::next_batch_work() const {
    return impl_->current.step == PiecewiseModelStep::BatchBoundary ? impl_->current.next_batch_work : 0;
}
bool piecewise_root_models(const acb_poly_t p, slong bits, slong precision, const cps_options &options,
                           PiecewiseRootStatistics &stats,
                           const std::function<bool(const RootModel &)> &visit,
                           const std::function<bool(const RootModel &)> &interest,
                           PiecewiseRootConfig config) {
    PiecewiseModelCursor cursor(p, bits, precision, stats, interest, config);
    while (!options.cancelled()) {
        switch (cursor.advance(options)) {
        case PiecewiseModelStep::Work:
        case PiecewiseModelStep::BatchBoundary: break;
        case PiecewiseModelStep::Model: if (!visit(cursor.model())) return false; break;
        case PiecewiseModelStep::Done: return true;
        case PiecewiseModelStep::Refused: return false;
        }
    }
    return false;
}
} // namespace certpolysolve::detail
