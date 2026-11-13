// The machine-arithmetic candidate engine: n root proposals in binary64 or scaled binary64 (machine::ScaledDouble,
// a binary64 mantissa with its own exponent), moved in packets of Aberth-Ehrlich sweeps (MachineEngine<R>, a
// NumericalEngine), and its factory prepare_machine_engine.
//
// prepare_machine_engine requires round-to-nearest and a storage estimate within storage_limit before allocating; the
// constructor scales coefficients and seeds into R's range (prepare_machine, machine_adapter.hpp) and prepares an
// evaluator (PreparedNewton in binary64, PreparedJets in scaled binary64). A sweep of advance evaluates p and p' at
// the active points, then moves each in turn, its S = sum_{j != i} 1/(z_i - z_j) seeing the points already moved
// (Gauss-Seidel). Aberth-Ehrlich's correction is N/(1 - N S), N = p/p', or 1/(p'/p - S) where p' = 0 or, in binary64,
// the exponent fields (field_exponent) show p' subnormal or N or N S out of range. From
// heuristics::machine_interaction_tree_degree on, binary64 with interaction_tree takes S from a
// MachineInteractionTree, built for a sweep when its counted terms undercut the direct sums. A point stops when its
// relative step is below 2^(8 - bits_) (Stop::SmallStep) or |p| below its Horner roundoff estimate (Stop::Roundoff).
// An arithmetic failure (an invalid, divide-by-zero or overflow flag, a disallowed underflow, a non-finite or
// coincident point) restores the packet's starting points (backup_) and clears valid_; export_candidates writes the
// points unscaled, as zero-radius balls.
//
// State: MachineEngine<R>, between advance calls: coefficients_, scale_, points_, stops_, backup_, newton_evaluator_ or
// evaluator_, valid_, and the tree's last costs (tree_build_terms_, tree_point_terms_).
// Entry points: prepare_machine_engine, advance and export_candidates, called by machine_newton_cover
// (certificate_attempt.cpp: binary64 with non-dyadic scaling and the tree, then scaled binary64) and
// StudyCandidates::advance (study_candidates.cpp: binary64 then scaled binary64). The research build hands every other
// request (the other methods, tiers and evaluators: Jenkins-Traub, Laguerre and Durand-Kerner, double-double and the
// experimental tiers, the JIT, compensated evaluation) to its copy of this engine as it was before W8's S10
// (prepare_research_machine_engine, research/solvers/machine_methods.cpp).
// Not here: the certificates (the callers' routes), the evaluators (candidate_jets.hpp), the tree
// (interaction_tree.hpp), the correction's and the sums' terms' operations (aberth_step.hpp), the multiprecision
// packets (candidates.cpp), the companion engines (research/solvers/).
// Design notes: docs/campaign-c3-algorithms.md (methods, tiers), docs/fast-methods-assessment-20261008.md (the tree).
#include "core/query.hpp"
#include "numeric/candidates/machine_candidates.hpp"
#include "numeric/candidates/aberth_step.hpp"
#include "numeric/candidates/candidate_jets.hpp"
#include "core/heuristics.hpp"
#include "numeric/candidates/interaction_tree.hpp"
#include "numeric/candidates/machine_adapter.hpp"
#include "platform/parallel.hpp"
#include "core/telemetry.hpp"
#include <algorithm>
#include <bit>
#include <cfenv>
#include <cstdint>
#include <limits>
#include <type_traits>

namespace certpolysolve::detail {
namespace {
// Range checks for Aberth's correction (MachineEngine::advance), in binary64's exponent range; the scaled tier carries
// its own exponents.
template <class R> constexpr bool binary64_exponent = std::is_same_v<R, double>;
// floor(log2 |x|) of a normal binary64, read from its exponent field: no arithmetic, so no floating-point
// exception. Zero and subnormal x read as min_exponent - 2, below every normal exponent.
int field_exponent(double x) {
    constexpr int fraction_bits = std::numeric_limits<double>::digits - 1;
    constexpr int bias = std::numeric_limits<double>::max_exponent - 1;
    constexpr uint64_t field = (uint64_t{1} << (std::numeric_limits<uint64_t>::digits - 1 - fraction_bits)) - 1;
    return static_cast<int>((std::bit_cast<uint64_t>(x) >> fraction_bits) & field) - bias;
}
// A complex number's, from its larger component (exponent e): |z| is at least 2^e when that component is
// normal, and below 2^(e + 3/2) (the component is below 2^(e + 1), and |z| at most sqrt(2) times it); the
// bound's 3/2, rounded up to a whole exponent, is complex_bound_log2.
template <class R> int field_exponent(machine::Complex<R> z) {
    return std::max(field_exponent(z.re), field_exponent(z.im));
}
constexpr int complex_bound_log2 = 2;
// N and N S below 2^representable_exponent are finite, and so is 1 - N S: one binade below binary64's overflow.
constexpr int representable_exponent = std::numeric_limits<double>::max_exponent - 1;
template <class R> class MachineEngine final : public NumericalEngine {
    using C = machine::Complex<R>;
    // These states describe candidate arithmetic, never certified roots.
    enum class Stop : uint8_t { Active, SmallStep, Roundoff };
    MachineScaling scale_;
    std::vector<C> coefficients_, points_, backup_;
    std::vector<Stop> stops_;
    std::vector<machine::Jet<R>> jets_;
    std::unique_ptr<machine::PreparedJets<R, 1>> evaluator_;
    std::unique_ptr<machine::PreparedNewton<R>> newton_evaluator_;
    std::vector<machine::NewtonEvaluation<R>> newton_values_;
    // Aberth's interaction sums from the interaction tree at large degree
    // (binary64 only), equal to the direct sum's to within its rounding: the
    // far fields carry as many terms as binary64 warrants, and the tree moves
    // each point as the sweep does, so far and near points alike are seen
    // where they are now (Gauss-Seidel, as the direct sum). Far fields left at
    // each sweep's starting positions would move the proposals enough to
    // change the certificate attempt's outcome. The proposals change, never a
    // certificate.
    // History: docs/code-history-20261021.md, machine_candidates.cpp 1.
    bool tree_ = false;
    std::vector<std::complex<double>> positions_;
    std::vector<size_t> near_;
    // What the tree costs, in the root interactions the packet counts: its
    // last construction's terms, and its terms a moving point in its last
    // sweep (the evaluation and the move). A tree term (a complex
    // multiply-add) and a direct pair (a complex division) cost about the
    // same.
    // History: docs/code-history-20261021.md, machine_candidates.cpp 2.
    uint64_t tree_build_terms_ = 0;
    double tree_point_terms_ = 0;
    bool valid_ = false;
    bool allow_underflow_ = false;
    // The significand bits of both arithmetics: binary64's, which is also ScaledDouble's mantissa.
    static constexpr double bits_ = machine::binary64_bits;

  public:
    MachineEngine(const acb_poly_t p, ConstRootView initial, bool allow_underflow, bool nondyadic_scaling,
                  bool interaction_tree)
        : tree_(interaction_tree && std::is_same_v<R, double>), allow_underflow_(allow_underflow) {
        valid_ = std::fegetround() == FE_TONEAREST &&
                 prepare_machine(coefficients_, points_, scale_, p, initial, nondyadic_scaling);
        if (!valid_)
            return;
        backup_.resize(points_.size());
        stops_.assign(points_.size(), Stop::Active);
        if constexpr (std::is_same_v<R, double>) {
            newton_evaluator_ = std::make_unique<machine::PreparedNewton<R>>(coefficients_);
            newton_values_.resize(points_.size());
        }
        if (!newton_evaluator_) {
            jets_.resize(points_.size());
            evaluator_ = std::make_unique<machine::PreparedJets<R, 1>>(coefficients_);
        }
    }
    bool valid() const { return valid_; }
#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
    // The preparation's bytes and work, which the research build's Portfolio trace reads
    // (research/solvers/portfolio.cpp).
    EngineCapabilities capabilities() const override {
        EngineCapabilities c;
        c.preparation_bytes = static_cast<uint64_t>(sizeof(C)) *
                              (coefficients_.capacity() + points_.capacity() + backup_.capacity());
        c.preparation_bytes += jets_.capacity() * sizeof(machine::Jet<R>);
        if (evaluator_)
            c.preparation_bytes += evaluator_->storage_bytes();
        if (newton_evaluator_)
            c.preparation_bytes += newton_evaluator_->storage_bytes() +
                                   newton_values_.capacity() * sizeof(machine::NewtonEvaluation<R>);
        c.preparation_bytes += stops_.capacity() * sizeof(Stop);
        c.preparation_work = coefficients_.size() * (newton_evaluator_ ? 2 : 1);
        return c;
    }
#endif
    EngineAdvance advance(unsigned sweeps, const cps_options &options,
                          certpolysolve::detail::CandidateStatistics &statistics) override {
        StageTimer timer(statistics.candidate_ns);
        EngineAdvance out;
        out.active = points_.size();
        if (!valid_ || std::fegetround() != FE_TONEAREST) {
            out.arithmetic_failure = true;
            return out;
        }
        // Restore the caller's flags; do not leave rounding-mode/global state
        // changes behind. FE_INEXACT is expected for candidate arithmetic.
        fenv_t environment;
        if (std::feholdexcept(&environment) != 0) {
            out.arithmetic_failure = true;
            return out;
        }
        struct Restore {
            fenv_t *value;
            ~Restore() { std::fesetenv(value); }
        } restore{&environment};
        backup_ = points_;
        for (unsigned sweep = 0; sweep < sweeps; ++sweep) {
            double largest = -std::numeric_limits<double>::infinity();
            std::unique_ptr<MachineInteractionTree> tree;
            double tree_terms = 0, tree_points = 0;
            if constexpr (std::is_same_v<R, double>) {
                if (tree_ && points_.size() >= heuristics::machine_interaction_tree_degree) {
                    // A sweep sums with the tree when its construction over every point and its terms at the
                    // moving ones count fewer than the direct sums there (n - 1 pairs a point). With every
                    // point moving that holds from machine_interaction_tree_degree on; with a few hundred
                    // left, as in a long packet's last quanta, the direct sums cost less.
                    // History: docs/code-history-20261021.md, machine_candidates.cpp 3.
                    const auto moving =
                        static_cast<double>(std::count(stops_.begin(), stops_.end(), Stop::Active));
                    if (!tree_build_terms_ || static_cast<double>(tree_build_terms_) + moving * tree_point_terms_ <
                                                  moving * static_cast<double>(points_.size() - 1)) {
                        StageTimer interaction_timer(statistics.interaction_ns);
                        positions_.resize(points_.size());
                        for (size_t j = 0; j < points_.size(); ++j)
                            positions_[j] = {points_[j].re, points_[j].im};
                        tree = std::make_unique<MachineInteractionTree>(positions_);
                        tree_build_terms_ = tree->build_terms();
                        statistics.root_interactions += tree_build_terms_;
                    }
                }
            }
            {
                StageTimer evaluation_timer(statistics.evaluation_ns);
                // The sweep's evaluations read only points_ and the evaluator, which evaluate leaves as it is, and
                // each writes its own slot. W9's split (parallel.hpp): in blocks of machine_split_block points from
                // machine_split_points points, on a solve with workers; the flags the workers raise reach this
                // thread's window at the join, so the packet's verdict below is today's.
                const size_t count = points_.size();
                const auto evaluate_range = [&](size_t first, size_t last) {
                    if (newton_evaluator_) {
                        for (size_t j = first; j < last; ++j)
                            if (stops_[j] == Stop::Active) newton_values_[j] = newton_evaluator_->evaluate(points_[j]);
                    } else
                        for (size_t j = first; j < last; ++j)
                            if (stops_[j] == Stop::Active) jets_[j] = evaluator_->evaluate(points_[j]);
                };
                if (count >= heuristics::machine_split_points && parallel_width() > 1) {
                    const size_t block = heuristics::machine_split_block;
                    parallel_for(SplitSite::MachineEvaluation, (count + block - 1) / block,
                                 [&](size_t b) { evaluate_range(b * block, std::min(count, (b + 1) * block)); });
                } else {
                    evaluate_range(0, count);
                }
                if (newton_evaluator_)
                    for (size_t j = 0; j < count; ++j)
                        if (stops_[j] == Stop::Active) statistics.evaluation_terms += newton_values_[j].terms;
                const auto active =
                    static_cast<uint64_t>(std::count(stops_.begin(), stops_.end(), Stop::Active));
                statistics.evaluations += active;
                if (!newton_evaluator_)
                    statistics.evaluation_terms += active * evaluator_->terms();
            }
            for (size_t i = 0; i < points_.size(); ++i) {
                if (options.cancelled()) {
                    out.candidates.cancelled = true;
                    break;
                }
                if (stops_[i] != Stop::Active)
                    continue;
                const C z = points_[i];
                const C value = newton_evaluator_ ? newton_values_[i].value : jets_[i].value;
                const C derivative = newton_evaluator_ ? newton_values_[i].denominator : jets_[i].derivative;
                const R absolute_sum = newton_evaluator_ ? newton_values_[i].absolute_sum : jets_[i].absolute_sum;
                // Complex Horner accumulates O(n*u) absolute roundoff.
                // This midpoint condition estimate is deliberately NOT a
                // certificate: it only asks for a more accurate tier before
                // a cancellation-dominated quotient damages the proposals.
                const R noise = R(std::exp2(-bits_)) * R(32 * static_cast<double>(coefficients_.size())) * absolute_sum;
                if (value == C{} || machine::norm(value) < noise) {
                    // This proposal cannot improve at the current tier.
                    // Other Aberth proposals can still make useful progress
                    // before the whole factor continues adaptively.
                    stops_[i] = Stop::Roundoff;
                    continue;
                }
                // Aberth's correction is N / (1 - N S), with N = p/p' and S the sum of 1/(z - z_j) over the
                // other points, or equally 1 / (p'/p - S). The first form is used wherever the arithmetic
                // represents N and N S; the second where it does not. N is out of range exactly where p'/p is
                // tiny: at a polynomial whose lowest nonconstant term has a high degree k, a point that a
                // sweep throws inside the unit circle sees p' ~ k c_k z^(k-1) underflow, to zero once
                // |z|^(k-1) < 2^-1074. That zero is no critical point, and Aberth's correction, -1/S there,
                // needs none: the second form takes it, as it does at a true critical point, and neither the
                // zero nor an N out of range ends the packet.
                // History: docs/code-history-20261021.md, machine_candidates.cpp 4.
                bool logarithmic = false;
                int newton_exponent = 0; // log2 |N| < newton_exponent (binary64 exponents only)
                if (derivative == C{})
                    logarithmic = true;
                else if constexpr (binary64_exponent<R>) {
                    newton_exponent = field_exponent(value) + complex_bound_log2 - field_exponent(derivative);
                    if (newton_evaluator_)
                        newton_exponent += field_exponent(newton_values_[i].coordinate_factor) + complex_bound_log2;
                    logarithmic = field_exponent(derivative) < std::numeric_limits<double>::min_exponent - 1 ||
                                  newton_exponent > representable_exponent;
                }
                // p'/p, the logarithmic derivative (the reversed evaluator's quotient is z p'/p).
                const auto logarithmic_derivative = [&] {
                    if (derivative == C{})
                        return C{};
                    const C quotient = derivative / value;
                    return newton_evaluator_ ? quotient / newton_values_[i].coordinate_factor : quotient;
                };
                // Newton's correction N for the first form, the logarithmic derivative p'/p for the second.
                C newton, log_derivative;
                if (logarithmic)
                    log_derivative = logarithmic_derivative();
                else if (newton_evaluator_) {
                    const auto ratio = newton_values_[i].ratio();
                    if (!ratio) {
                        out.exhausted = true;
                        continue;
                    }
                    newton = *ratio;
                } else
                    newton = value / derivative;
                C sum;
                {
                    StageTimer interaction_timer(statistics.interaction_ns);
                    bool summed = false;
                    if constexpr (std::is_same_v<R, double>) {
                        std::complex<double> far;
                        double tail = 0;
                        size_t far_fields = 0;
                        if (tree && tree->evaluate(i, far, near_, tail, far_fields)) {
                            // A far field counts as the terms of its series.
                            sum = C(far.real(), far.imag());
                            for (const size_t j : near_) {
                                const C difference = z - points_[j];
                                if (difference == C{}) {
                                    out.arithmetic_failure = true;
                                    break;
                                }
                                sum += machine::interaction_term(difference);
                            }
                            const auto terms = near_.size() + far_fields * MachineInteractionTree::order;
                            statistics.root_interactions += terms;
                            tree_terms += static_cast<double>(terms);
                            ++tree_points;
                            summed = true;
                        }
                    }
                    for (size_t j = 0; !summed && j < points_.size(); ++j)
                        if (j != i) {
                            const C difference = z - points_[j];
                            if (difference == C{}) {
                                out.arithmetic_failure = true;
                                break;
                            }
                            sum += machine::interaction_term(difference);
                            ++statistics.root_interactions;
                        }
                    if constexpr (binary64_exponent<R>)
                        if (!logarithmic &&
                            newton_exponent + field_exponent(sum) + complex_bound_log2 > representable_exponent) {
                            logarithmic = true; // N S is out of range: the second form
                            log_derivative = logarithmic_derivative();
                        }
                }
                const C delta = logarithmic ? machine::logarithmic_aberth_correction(log_derivative, sum)
                                            : machine::aberth_correction(newton, sum);
                if (out.arithmetic_failure || !machine::finite(delta) || !machine::finite(z - delta)) {
                    out.arithmetic_failure = true;
                    break;
                }
                const double size = machine::logabs(z);
                largest = std::max(largest, machine::logabs(delta) - (std::isfinite(size) ? size : 0));
                points_[i] = z - delta;
                if constexpr (std::is_same_v<R, double>) {
                    if (tree) {
                        const auto terms = tree->move(i, {points_[i].re, points_[i].im});
                        statistics.root_interactions += terms;
                        tree_terms += static_cast<double>(terms);
                    }
                }
                if (machine::logabs(delta) - (std::isfinite(size) ? size : 0) < 8 - bits_)
                    stops_[i] = Stop::SmallStep;
                ++out.work;
            }
            if (tree_points > 0)
                tree_point_terms_ = tree_terms / tree_points;
            ++statistics.candidate_sweeps;
            out.candidates.correction_log2 = largest;
            out.active = static_cast<uint64_t>(std::count(stops_.begin(), stops_.end(), Stop::Active));
            out.candidates.converged = out.active == 0;
            out.arithmetic_limited =
                out.active == 0 && std::find(stops_.begin(), stops_.end(), Stop::Roundoff) != stops_.end();
            if (out.arithmetic_failure || out.arithmetic_limited || out.exhausted || out.candidates.cancelled ||
                out.candidates.converged)
                break;
        }
        const auto exceptions = std::fetestexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);
        // Exponent-alignment underflow is intentional in ScaledDouble, but its
        // mantissa multiplication remains normalized. Binary64 rejects it unless the caller allows it.
        out.arithmetic_failure = out.arithmetic_failure ||
                                 (exceptions & (FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW)) ||
                                 (!allow_underflow_ && !std::is_same_v<R, machine::ScaledDouble> && (exceptions & FE_UNDERFLOW));
        for (auto z : points_)
            out.arithmetic_failure = out.arithmetic_failure || !machine::finite(z);
        if (out.arithmetic_failure) {
            points_ = backup_;
            valid_ = false;
        }
        out.candidates.finite = !out.arithmetic_failure;
        return out;
    }
    void export_candidates(acb_ptr output) const override {
        for (size_t i = 0; i < points_.size(); ++i) {
            auto z = points_[i];
            if constexpr (std::is_same_v<R, double>)
                if (scale_.fraction != 0) {
                    const double factor = std::exp2(scale_.fraction);
                    z = {z.re * factor, z.im * factor};
                }
            to_ball(output + i, z, scale_.variable);
        }
    }
};
template <class R>
std::unique_ptr<NumericalEngine> make(const acb_poly_t p, ConstRootView initial, bool allow_underflow = false,
                                      bool nondyadic_scaling = false, bool interaction_tree = false) {
    auto engine = std::make_unique<MachineEngine<R>>(p, initial, allow_underflow, nondyadic_scaling, interaction_tree);
    if (!engine->valid())
        return {};
    return engine;
}
} // namespace
std::unique_ptr<NumericalEngine> prepare_machine_engine(IterationMethod method, ArithmeticTier tier,
                                                        const acb_poly_t p, ConstRootView initial,
                                                        bool compensated, bool jit, bool allow_underflow, uint64_t storage_limit,
                                                        bool nondyadic_scaling, bool interaction_tree) {
    // This engine serves Aberth in binary64 and scaled binary64, without compensated evaluation or the JIT: the only
    // engine the recipes ask for. The research build hands every other request to its copy of the engine as it was
    // before W8's S10; the default build has no other engine.
    const bool served = method == IterationMethod::Aberth && !compensated && !jit &&
                        (tier == ArithmeticTier::Binary64 || tier == ArithmeticTier::ScaledBinary64);
    if (!served) {
#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
        return prepare_research_machine_engine(method, tier, p, initial, compensated, jit, allow_underflow,
                                               storage_limit, nondyadic_scaling, interaction_tree);
#else
        return {};
#endif
    }
    if (std::fegetround() != FE_TONEAREST)
        return {};
    fenv_t environment;
    if (std::feholdexcept(&environment) != 0)
        return {};
    struct Restore {
        fenv_t *value;
        ~Restore() { std::fesetenv(value); }
    } restore{&environment};
    const auto n = acb_poly_degree(p);
    // Admission before allocation: 16 * 64 bytes a coefficient, the estimate of the engine that also had the
    // research methods' scratch.
    if (n < 1 || !resources::storage_fits(resources::multiply(static_cast<uint64_t>(n) + 1, 16 * 64), storage_limit))
        return {};
    if (tier == ArithmeticTier::Binary64)
        return make<double>(p, initial, allow_underflow, nondyadic_scaling, interaction_tree);
    return make<machine::ScaledDouble>(p, initial);
}
} // namespace certpolysolve::detail
