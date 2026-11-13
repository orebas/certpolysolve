#pragma once
#include <certpolysolve/certpolysolve.h>
#include "numeric/evaluate/monic_source.hpp"
#include "numeric/evaluate/evaluator.hpp"
#include "input/secular/secular.hpp"
#include "core/service_meter.hpp"
#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
#include "core/engine_config.hpp"
#include "core/number.hpp"
#endif
#include <functional>
#include <memory>
#include <optional>
#include <map>
#include <vector>

namespace certpolysolve::detail {
using PolynomialEncloser = std::function<void(acb_poly_t, slong)>;
using ExactTaylorEncloser = std::function<bool(acb_poly_t, const acb_t, slong)>;

struct PolynomialSourceOptions {
#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
    EvaluationMethod evaluation = EvaluationMethod::Coefficients; // TaylorCache: the Series test engine's
#endif
    bool candidate_direct = false;
    bool sparse_coefficients = false;
    InputSemantics input = InputSemantics::Unknown;
    ExactTaylorEncloser exact_taylor = {};
    // Reads the owned exact cache's cumulative raw service work. Present
    // only in the private service-meter experiment; never borrows a caller.
    std::function<double()> exact_taylor_work;
    // The precision from which a source with exact translations takes them for every point: the solve's
    // precision budget, where no more precision can be asked for. Below it the value from the coefficients is
    // taken wherever it has a correct bit (PolynomialSource::evaluate_impl). Zero: translations at every
    // precision.
    slong exact_precision = 0;
#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
    // The exact polynomial this source evaluates (one of the two), shared
    // with its encloser, for one use: the exact value at a node whose
    // residue the precision budget leaves inadequate (evaluate_exact; the
    // research build's EnginePolicy::exact_residues_at_budget).
    std::shared_ptr<const IntegerPolynomial> exact_integer;
    std::shared_ptr<const std::vector<GaussianRational>> exact_gaussian;
#endif
    // A centred source, whose polynomial is p(z + c) up to a constant factor: the uncentred p and c (exact).
    // The accuracy goals measure each root against its distance from p's origin. Where p's own expansion
    // has smaller terms at z + c than the centred one at z (a group of roots near p's origin, which the
    // centring carried away from the new one), values and local models come from p (evaluate_impl, taylor):
    // the same polynomial, with the rounding error of an expansion about the nearer origin. The term sizes
    // are log2 |coefficient| (-inf for a zero one) of the centred and the uncentred polynomial, from their
    // exact coefficients. Empty without a centre.
    PolynomialEncloser origin;
    ComplexBall origin_offset;
    std::vector<double> centred_log2, origin_log2;
};

class PolynomialSource final : public MonicSource {
    PolynomialEncloser enclose_;
    slong degree_;
    bool real_;
    PolynomialSourceOptions options_;

  public:
    struct Prepared {
        PolynomialSnapshot polynomial;
        std::unique_ptr<Evaluator> evaluator;
        slong precision;
        Prepared(ComplexPolynomial coefficients, slong bits)
            : polynomial(std::move(coefficients)), precision(bits) {}
    };
    using Coefficients = std::shared_ptr<const Prepared>;

  private:
    // The largest term of a monic expansion at |x| = 2^t, log2, from the upper hull of (k, log2 |a_k|).
    struct TermHull {
        std::vector<std::pair<double, double>> vertices; // (k, log2 |a_k|) in increasing k
        double leading = 0;
        bool available = false;
        void build(const std::vector<double> &log2_coefficients);
        double largest_log2(double t) const;
    };
    std::map<slong, Coefficients> prepared_, origin_prepared_;
    uint64_t taylor_shifts_ = 0, taylor_pairs_ = 0, taylor_heads_ = 0;
    bool service_meter_ = false;
    ServiceTaylorQuote last_taylor_quote_;
    TermHull centred_terms_, origin_terms_;
    std::vector<std::pair<ComplexBall, Coefficients>> translated_;
    bool evaluate_impl(acb_t value, acb_t derivative, const acb_t z, slong precision, bool joint);
    std::optional<bool> coefficient_sparsity(slong precision);
    enum class Translated { Unavailable, Evaluated, Failed };
    Translated evaluate_translated(const ExactTaylorEncloser &encloser, acb_t value, acb_t derivative,
                                   const acb_t z, slong precision, bool joint);

  public:
    PolynomialSource(PolynomialEncloser enclose, slong degree, bool real,
                     PolynomialSourceOptions options = {})
        : enclose_(std::move(enclose)), degree_(degree), real_(real), options_(options) {
        if (options_.origin) {
            centred_terms_.build(options_.centred_log2);
            origin_terms_.build(options_.origin_log2);
        }
    }
    slong degree() const override { return degree_; }
    bool real() const override { return real_; }
    SourceMetadata metadata() const override {
        return {degree_,
                real_,
                true,
                SourceRepresentation::Coefficients,
                options_.input,
                true};
    }
    // Whether the source evaluates through the research build's local Taylor cache (the Series test engine;
    // LocalPolynomialEvaluator); false in the default build, which evaluates the coefficients.
#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
    bool local() const { return options_.evaluation == EvaluationMethod::TaylorCache; }
#else
    bool local() const { return false; }
#endif
    bool direct() const { return options_.candidate_direct; }
    bool sparse_coefficients() const { return options_.sparse_coefficients; }
    bool exact_local_models() const { return static_cast<bool>(options_.exact_taylor); }
    bool midpoint_evaluation() const {
        return direct() && !local() && !options_.exact_taylor;
    }
    uint64_t terms() const override {
        return static_cast<uint64_t>(degree_ + 1);
    }
    bool evaluate(acb_t value, acb_t derivative, const acb_t z, slong precision) override;
    bool evaluate_value(acb_t value, const acb_t z, slong precision) override;
    // Whether z is evaluated through the uncentred polynomial (a centred source, at a point where that
    // expansion's terms are the smaller); w receives z + c, exactly. An evaluator that takes the source's
    // coefficients for itself follows the same rule: origin_enclosure's at w, coefficient_enclosure's at z.
    bool origin_point(acb_t w, const acb_t z);
    Coefficients origin_enclosure(slong precision);
    // A handle pins polynomial AND evaluator storage across cache eviction and
    // source destruction. Candidates may use midpoints; proof paths use balls.
    Coefficients coefficient_enclosure(slong precision);
    bool rectangular_available(slong precision);
    // The same coefficient source, with few nonzero coefficients
    // (CoefficientEvaluator's own classification).
    bool sparse_available(slong precision);
#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
    // The research build's exact residues (research/solvers/exact_residues.cpp).
    // The value (and derivative) of the monic polynomial at an exact z,
    // computed exactly in Gaussian-integer arithmetic and rounded once to
    // the precision: no cancellation beyond that rounding, at a cost of
    // about degree x (coefficient bits + degree x z's bits) bit operations.
    // For a residue the precision budget leaves inadequate
    // (regenerate_residues). False when the source holds no exact
    // polynomial or z is not exact; nothing is written then.
    bool exact_available() const {
        return options_.exact_integer || options_.exact_gaussian;
    }
    bool evaluate_exact(acb_t value, acb_t derivative, const acb_t z, slong precision, bool joint);
    // Whether evaluate_exact at z fits the exact budgets: its Horner
    // accumulator holds about the coefficient bits plus degree x (z's
    // mantissa bits + |z's exponent|), each of the degree steps multiplies
    // it by z's mantissa. Refused beyond exact_operations or the exact
    // storage limit, and, without exponent expansion, beyond the owned
    // polynomial's storage (as the exact centering is). operations is the
    // count of word operations, for the meter: a count and not its log2, whose
    // round trip through the C library's exp2 can differ in its last bits
    // between platforms.
    // History: docs/code-history-20261021.md, structured_source.hpp 1.
    bool exact_fits(const acb_t z, const cps_options &options, double &operations) const;
#endif
    void research_service_meter(bool enabled) { service_meter_ = enabled; }
    double research_last_taylor_work() const { return last_taylor_quote_.work; }
    ServiceOperation research_last_taylor_operation() const { return last_taylor_quote_.operation; }
    double research_exact_taylor_work() const {
        return options_.exact_taylor_work ? options_.exact_taylor_work() : 0;
    }
    bool taylor(acb_poly_t out, const acb_t center, slong precision) override;
    // Coefficients k - 1 and k of taylor()'s expansion (the same polynomial and normalisation), without the
    // expansion: with q = p^(k-1) / (k-1)!, they are q(center) and q'(center) / k, one joint evaluation of a
    // polynomial of degree n - k + 1. A Newton step on p^(k-1), which is what a cluster's centre moves by, reads
    // exactly these two. An exact local model (small inputs) is expanded whole.
    bool taylor_pair(acb_t lower, acb_t upper, const acb_t center, slong precision, slong k);
    // The first `length` coefficients of taylor()'s expansion (the same polynomial and normalisation) by Horner's
    // rule on `length` accumulators: over the coefficients a_i from the top, T_j <- T_j x + T_(j-1) for j from
    // length - 1 down to 1 and T_0 <- T_0 x + a_i leave the expansion's coefficient j in T_j. The accumulators
    // are exact numbers and each one's error is one disk radius, so nothing is lost to the rectangles of ball
    // arithmetic at a complex centre (half a bit a term). `length` passes of n multiplications at the precision,
    // where FLINT's whole shift costs many such passes at a restart's lengths and precisions. Length one is the
    // value at the centre. Refused for an exact local model (small inputs: taylor() is their expansion).
    // History: docs/code-history-20261021.md, structured_source.hpp 2.
    //
    // Where the polynomial has few terms the head is summed term by term instead: coefficient j is the sum over
    // the terms a_i x^i of C(i, j) a_i center^(i - j), a power by repeated squaring and `length` products a term
    // (termwise_head_multiplications), taken where that is fewer multiplications than the passes'.
    bool taylor_head(acb_poly_t out, const acb_t center, slong precision, slong length);
    // What taylor_head(length) costs, in passes of degree + 1 multiplications: `length` by Horner's rule, less
    // term by term for a polynomial of few terms. Nothing for a source that is not evaluated from its
    // coefficients (taylor_head's own refusals, and experimental evaluators).
    std::optional<double> taylor_head_passes(slong precision, slong length);
    // An upper bound on sum over j >= length of |e_j| |offset|^j, where e_j are taylor()'s coefficients at the
    // exact `center`: the terms an expansion's head drops, at center + offset. From magnitudes only, term by
    // term (see the definition), so no evaluation's cancellation enters it; it ignores cancellation between
    // terms too, and is tight where one term carries the tail. False where taylor_head refuses.
    bool taylor_tail_bound(mag_t out, const acb_t center, const acb_t offset, slong precision, slong length);
    // The expansions computed so far (diagnostics and tests): whole shifts, pairs and heads.
    uint64_t taylor_shifts() const { return taylor_shifts_; }
    uint64_t taylor_pairs() const { return taylor_pairs_; }
    uint64_t taylor_heads() const { return taylor_heads_; }
};

// Evaluates Q*(1-sum a/(x-b)), without materializing its coefficients. At an
// original pole use the polynomial limit, including its derivative.
class SecularSource final : public MonicSource {
    cps_secular input_;
    slong precision_ = 0;
    std::vector<ComplexBall> residues_, poles_;
    bool real_ = true;

  public:
    explicit SecularSource(const cps_secular &input);
    slong degree() const override { return static_cast<slong>(input_.terms.size()); }
    bool real() const override { return real_; }
    SourceMetadata metadata() const override {
        return {degree(), real_, true, SourceRepresentation::SecularNumerator, InputSemantics::Exact, false};
    }
    void enclose_terms(slong precision);
    const std::vector<ComplexBall> &residues() const { return residues_; }
    const std::vector<ComplexBall> &poles() const { return poles_; }
    bool evaluate(acb_t value, acb_t derivative, const acb_t z, slong precision) override;
};
// Candidate-only eigenvalues of strongly coupled 2-by-2 principal blocks of
// diag(b)+a*1^T. No coefficient expansion or global eigenvalue claim is made.
bool secular_block_seeds(acb_ptr out, SecularSource &source, slong precision, const cps_options &options);
} // namespace certpolysolve::detail
