#pragma once
// Owns one exact coefficient source and its precision-keyed evaluator cache.
// The full-input finisher and zero-deflated discovery use separate instances.
// Statistics are borrowed from the stable outer state. Calls are exclusive;
// no backend view, caller point or options reference is retained.
#include "certify/chart_certificate.hpp"
#include "numeric/evaluate/disk_blocked_polynomial_evaluator.hpp"
#include "numeric/evaluate/polynomial_remainder.hpp"
#include "core/dovetail.hpp"
#include "numeric/arithmetic/modulus_bound.hpp"
#include "methods/piecewise/piecewise_statistics.hpp"
#include "methods/piecewise/piecewise_source_certificate.hpp"
#include "core/query.hpp"
#include "core/telemetry.hpp"
#include <memory>
#include <optional>
#include <utility>
#include <vector>

namespace certpolysolve::detail {
class PiecewiseFinisherSource final {
    // Newton ratios and the certificate are homogeneous. Keep the original
    // coefficients for both evaluation and the remainder; this object has
    // no monic-source interface or unused normalization callback.
    const std::vector<Coefficient> coefficients;
    PiecewiseRootStatistics &statistics;
    // A precision label, remainder bound and evaluator are one cache. Build a
    // replacement completely before committing it, so an exception cannot
    // leave an old label paired with values from different precisions.
    struct Cache {
        const slong precision;
        PolynomialLinearRemainder remainder;
        SharedDiskBlockedPolynomialEvaluator blocked;
        Cache(PolynomialSnapshot snapshot, slong bits)
            : precision(bits), remainder(snapshot.coefficients()), blocked(std::move(snapshot), bits) {}
    };
    std::unique_ptr<Cache> cache;
    PolynomialSnapshot original_snapshot(slong bits) const {
        ComplexPolynomial polynomial;
        ComplexBall coefficient;
        for (size_t i = 0; i < coefficients.size(); ++i) {
            coefficients[i].enclose(coefficient.value(),bits);
            acb_poly_set_coeff_acb(polynomial.value(),static_cast<slong>(i),coefficient.value());
        }
        return PolynomialSnapshot(std::move(polynomial));
    }
    Cache &prepare_cache(slong bits) {
        if (!cache || cache->precision != bits) {
            StageTimer preparation(statistics.finisher_preparation_ns);
            // Preserve the established conservative work charge and admission
            // estimate while measuring the eliminated coefficient copy.
            statistics.work.finisher_preparation +=
                3 * static_cast<double>(coefficients.size()) * operation_weight(bits);
            auto replacement = std::make_unique<Cache>(original_snapshot(bits), bits);
            cache = std::move(replacement);
        }
        return *cache;
    }
  public:
    PiecewiseFinisherSource(std::vector<Coefficient> owned_coefficients, PiecewiseRootStatistics &s)
        : coefficients(std::move(owned_coefficients)), statistics(s) {}
    PiecewiseFinisherSource(const PiecewiseFinisherSource &) = delete;
    PiecewiseFinisherSource &operator=(const PiecewiseFinisherSource &) = delete;
    PiecewiseFinisherSource(PiecewiseFinisherSource &&) = delete;
    PiecewiseFinisherSource &operator=(PiecewiseFinisherSource &&) = delete;
    bool evaluate(acb_t value, acb_t derivative, const acb_t point, slong bits) {
        ++statistics.finisher_evaluations;
        statistics.finisher_terms += coefficients.size();
        statistics.work.finisher += 2 * static_cast<double>(coefficients.size()) * operation_weight(bits);
        auto &prepared = prepare_cache(bits);
        ++statistics.finisher_block_evaluations;
        return prepared.blocked.evaluate(value,derivative,point);
    }
    std::optional<PiecewiseSourceCertificate> certify(const acb_t point,const cps_options &options,slong bits) {
        if (!acb_is_finite(point) || acb_is_zero(point) || options.cancelled()) return std::nullopt;
        Root root;
        RootModel model;
        acb_get_mid(model.center.value(),point);
        upper(model.radius.value(),model.center.value());
        arb_mul_2exp_si(model.radius.value(),model.radius.value(),-static_cast<slong>(options.bits)+4);
        ComplexBall value,derivative,linear,zero;
        if (!evaluate(value.value(),derivative.value(),model.center.value(),bits)) return std::nullopt;
        prepare_cache(bits).remainder.linear_remainder(model.error.value(),model.center.value(),model.radius.value());
        acb_poly_set_coeff_acb(model.polynomial.value(),0,value.value());
        acb_mul_arb(linear.value(),derivative.value(),model.radius.value(),bits);
        acb_poly_set_coeff_acb(model.polynomial.value(),1,linear.value());
        if (!certify_root_model(root,model,zero.value(),options,bits,nullptr,nullptr,nullptr,nullptr,0,false,
                                &statistics.work.certificate)) return std::nullopt;
        root.certificate="original-coefficient disk evaluation and second-derivative Taylor remainder";
        return PiecewiseSourceCertificate(*this, std::move(root));
    }
};

} // namespace certpolysolve::detail
