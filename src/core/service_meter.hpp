#pragma once
// Private operation prices for the scheduling experiment. No clock enters a
// charge. One owner adds each returned price to its native work ledger.
#include "core/dovetail.hpp"
#include "core/engine_policy.hpp"
#include <array>

namespace certpolysolve::detail {
enum class ServiceOperation : size_t {
    Binary64, ScaledBinary64, Evaluation, Interaction, Regeneration, Taylor, Graeffe, Auxiliary, ExactTaylor,
    ShortProduct, LimbHorner, GraeffeWords, TaylorDivideConquer, Count
};
inline constexpr std::array<const char *, static_cast<size_t>(ServiceOperation::Count)> service_names{
    "binary64", "scaled_binary64", "evaluation", "interaction", "regeneration", "taylor", "graeffe", "auxiliary", "exact_taylor",
    "short_product", "limb_horner", "graeffe_words", "taylor_divide_conquer"};
// Frozen calibration: subdivision-service-meter-inventory4-20261110/prices-frozen-4.json.
// SHA256 4f284e0fc12c3f26d34124c6c02d74d7f694dcd5f6efbb5538e779d287b0b1f1
// The first twelve prices retain the version8 inventory (the exact-cache
// basis retains its earlier supported price). D&C Taylor's separate affine
// word calibration is documented at service_taylor_quote below.
inline constexpr std::array<double, service_names.size()> service_prices{
    2.518769401821626, 124.29011657372389, 2.4061813353521755, 9.009513124565773, 3.8581367328511074, 0.03720070013556002, 0.5407947955112974, 0.10307652283197645, 1.3632375668513286, 107.15751875272649, 1.5235059966736468, 8.40451284599325, 1.422632600479756};
inline double service_price(ServiceOperation operation, double work) {
    return work * service_prices[static_cast<size_t>(operation)];
}
struct ServicePolynomialWork {
    double scalar = 0, transform = 0;
    ServicePolynomialWork &operator+=(const ServicePolynomialWork &other) {
        scalar += other.scalar; transform += other.transform; return *this;
    }
    double total() const { return scalar + transform; }
};
inline double service_graeffe_price(ServicePolynomialWork work, slong bits) {
    // Small block pairs retain scalar multiplication cost. Large polynomial
    // transforms pay word work per butterfly, not one full acb multiplication.
    // The crossover is the existing min(schoolbook, transform) proxy.
    return service_price(ServiceOperation::Graeffe, work.scalar * operation_weight(bits)) +
        service_price(ServiceOperation::GraeffeWords, work.transform * std::ceil(static_cast<double>(bits) / FLINT_BITS));
}
inline void service_observe(const EnginePolicy *policy, ServiceOperation operation,
                            double work, int64_t bits = 0) {
    if (!policy || !policy->research_service_meter || !policy->trace || work == 0) return;
    TraceEvent event{"service_charge", "work", static_cast<uint64_t>(std::max<int64_t>(0, bits)),
                     0, nullptr, nullptr};
    event.reason = service_names[static_cast<size_t>(operation)];
    event.metered_work = work;
    event.arithmetic_bits = event.precision;
    policy->trace(policy->trace_data, event);
}
inline double service_charge(const EnginePolicy *policy, ServiceOperation operation,
                             double work, int64_t bits = 0) {
    if (!policy || !policy->research_service_meter) return work;
    const double priced = service_price(operation, work);
    service_observe(policy, operation, priced, bits);
    return priced;
}
enum class ServiceTaylorAlgorithm { Identity, Horner, Convolution, DivideConquer };
inline ServiceTaylorAlgorithm service_taylor_algorithm(const acb_poly_t source, const acb_t center, slong bits) {
    const double n = static_cast<double>(acb_poly_length(source));
    if (acb_is_zero(center)) return ServiceTaylorAlgorithm::Identity;
    // FLINT 3.6 acb_poly/taylor_shift.c also chooses Horner for short
    // centers at high precision. These are that kernel's dispatch cutoffs,
    // not solver thresholds. Square the positive inequality to avoid sqrt.
    const double center_bits = static_cast<double>(acb_bits(center));
    const bool horner = n <= 30 ||
        (n <= 500 && center_bits == 1 && (n - 30) * (n - 30) < 9 * static_cast<double>(bits)) ||
        (n <= 100 && 100 * center_bits < static_cast<double>(bits));
    if (horner) return ServiceTaylorAlgorithm::Horner;
    return static_cast<double>(bits) > 2 * n ? ServiceTaylorAlgorithm::Convolution
                                           : ServiceTaylorAlgorithm::DivideConquer;
}
inline double service_taylor_base(const acb_poly_t source, const acb_t center, slong bits) {
    const auto algorithm = service_taylor_algorithm(source, center, bits);
    const double n = static_cast<double>(acb_poly_length(source));
    if (algorithm == ServiceTaylorAlgorithm::Identity) return n * operation_weight(bits);
    const double operations = algorithm == ServiceTaylorAlgorithm::Horner ? n * n :
        block_taylor_shift_operations(source, bits);
    return operations * operation_weight(bits);
}
struct ServiceTaylorQuote {
    ServiceOperation operation = ServiceOperation::Taylor;
    double work = 0;
    double price() const { return service_price(operation, work); }
    void observe(const EnginePolicy *policy, slong bits) const {
        service_observe(policy, operation, price(), bits);
    }
};
inline ServiceTaylorQuote service_taylor_quote(const acb_poly_t source, const acb_t center, slong bits,
                                               bool rotated_source = false) {
    if (service_taylor_algorithm(source, center, bits) != ServiceTaylorAlgorithm::DivideConquer)
        return {ServiceOperation::Taylor, service_taylor_base(source, center, bits)};
    // Frozen primitive calibration: subdivision-m4-integration-20261111/
    // taylor-price-2/price-frozen.json. Fixed overhead plus limb work, with
    // one/two real products for a real shift of real/complex coefficients,
    // and four for complex intermediate products. This is a scheduling
    // proxy; exact short centers can be substantially cheaper.
    constexpr double overhead_words = 11.586346649051508;
    const double real_products = !arb_is_zero(acb_imagref(center)) ? 4 :
        (!rotated_source && acb_poly_is_real(source) ? 1 : 2);
    const double words = std::ceil(static_cast<double>(bits) / FLINT_BITS);
    return {ServiceOperation::TaylorDivideConquer,
            block_taylor_shift_operations(source, bits) * real_products * (overhead_words + words)};
}

// A Graeffe transform squares the even and odd polynomials separately. Each
// complex square performs a², b² and ab, with exactly zero components cheap.
// Inspect their actual exponent profiles, retaining the existing block-pair
// proxy. No coefficient copy, backend call or clock is needed for this price.
inline ServicePolynomialWork service_graeffe_work(const acb_poly_t source, slong bits) {
    struct Profile {
        std::vector<std::pair<slong, double>> exponents;
        slong length = 0;
    };
    const auto profile = [&](slong parity, bool imaginary) {
        Profile p;
        for (slong j = parity; j < acb_poly_length(source); j += 2) {
            const arb_struct *value = imaginary ? acb_imagref(source->coeffs + j) : acb_realref(source->coeffs + j);
            if (!arb_is_zero(value)) p.length = j / 2 + 1;
            const arf_struct *midpoint = arb_midref(value);
            if (!arf_is_zero(midpoint) && arf_is_finite(midpoint) && fmpz_fits_si(ARF_EXPREF(midpoint)))
                p.exponents.emplace_back(j / 2, static_cast<double>(fmpz_get_si(ARF_EXPREF(midpoint))));
        }
        return p;
    };
    const auto blocks = [&](const Profile &p, double slope) {
        std::vector<double> lengths;
        if (!p.length) return lengths;
        if (p.exponents.size() < 2) return std::vector<double>{static_cast<double>(p.length)};
        const double height = 3 * static_cast<double>(bits) + 512; // existing block proxy
        double top = 0, bottom = 0;
        slong start = p.exponents.front().first;
        bool open = false;
        for (const auto &[index, exponent] : p.exponents) {
            const double value = exponent - slope * static_cast<double>(index);
            if (!open) { top = bottom = value; open = true; continue; }
            const double next_top = std::max(top, value), next_bottom = std::min(bottom, value);
            if (next_top - next_bottom < height) { top = next_top; bottom = next_bottom; }
            else { lengths.push_back(static_cast<double>(index - start)); start = index; top = bottom = value; }
        }
        lengths.push_back(static_cast<double>(p.length - start));
        return lengths;
    };
    const auto product = [&](const Profile &a, const Profile &b) {
        if (!a.length || !b.length) return ServicePolynomialWork{};
        double difference = 0, span = 0;
        for (const auto *p : {&a, &b})
            if (!p->exponents.empty()) {
                difference += p->exponents.back().second - p->exponents.front().second;
                span += static_cast<double>(p->exponents.back().first - p->exponents.front().first);
            }
        // The mixed product uses one common slope, as the backend does.
        const double slope = span != 0 ? std::floor(difference / span + .5) : 0;
        const auto x = blocks(a, slope), y = blocks(b, slope);
        ServicePolynomialWork work;
        for (double u : x) for (double v : y) {
            const double scalar = u*v, transform = 3*(u+v)*portable::log2(2*(u+v));
            if (scalar <= transform) work.scalar += scalar;
            else work.transform += transform;
        }
        return work;
    };
    ServicePolynomialWork work;
    for (slong parity : {0, 1}) {
        const auto a = profile(parity, false), b = profile(parity, true);
        work += product(a, a); work += product(b, b); work += product(a, b);
    }
    // Two decimation/copy passes and two assembly/sign passes.
    work.transform += 4 * static_cast<double>(acb_poly_length(source));
    return work;
}
} // namespace certpolysolve::detail
