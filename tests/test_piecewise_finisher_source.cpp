// Cache replacement retains a usable previous precision after ordinary C++
// allocation failures. FLINT's allocator-failure policy is independent.
#include "methods/piecewise/finishing/piecewise_finisher_source.hpp"
#include <certpolysolve/certpolysolve.hpp>
#include <cstdlib>
#include <iostream>
#include <new>
#include <stdexcept>

namespace d = certpolysolve::detail;
static thread_local long remaining = -1;
void *operator new(std::size_t size) {
    if (remaining == 0) { remaining = -1; throw std::bad_alloc(); }
    if (remaining > 0) --remaining;
    if (auto *p = std::malloc(size ? size : 1)) return p;
    throw std::bad_alloc();
}
void *operator new[](std::size_t size) { return ::operator new(size); }
void operator delete(void *p) noexcept { std::free(p); }
void operator delete[](void *p) noexcept { std::free(p); }
void operator delete(void *p, std::size_t) noexcept { std::free(p); }
void operator delete[](void *p, std::size_t) noexcept { std::free(p); }

static void require(bool value, const char *message) {
    if (!value) throw std::runtime_error(message);
}

static_assert(!std::is_default_constructible_v<d::PiecewiseSourceCertificate>);
static_assert(!std::is_constructible_v<d::PiecewiseSourceCertificate, const d::PiecewiseFinisherSource &, d::Root &&>);
static_assert(!std::is_copy_constructible_v<d::PiecewiseSourceCertificate>);
static_assert(!std::is_move_assignable_v<d::PiecewiseSourceCertificate>);
static_assert(std::is_nothrow_move_constructible_v<d::PiecewiseSourceCertificate>);

static void require_consumed(d::PiecewiseSourceCertificate &certificate) {
    bool view_refused = false, release_refused = false;
    try { (void)certificate.proof(); } catch (const std::logic_error &) { view_refused = true; }
    try { (void)std::move(certificate).release(); } catch (const std::logic_error &) { release_refused = true; }
    require(view_refused && release_refused, "consumed source certificate remained usable");
}

static void source_certificates() {
    std::vector<d::Coefficient> coefficients(2);
    coefficients[0].real = d::Number("-1");
    coefficients[1].real = d::Number("1");
    d::PiecewiseRootStatistics statistics, other_statistics;
    d::PiecewiseFinisherSource source(coefficients, statistics), other(coefficients, other_statistics);
    cps_options options;
    options.goal = CPS_ROOT_APPROXIMATION;
    options.bits = 53;
    d::ComplexBall point;
    acb_one(point.value());
    auto certificate = source.certify(point.value(), options, 128);
    require(certificate.has_value(), "linear original-source certificate failed");
    require(certificate->belongs_to(source) && !certificate->belongs_to(other),
            "certificate lost its originating source identity");
    const d::Root &proof = certificate->proof();
    // Independent exact planted root: the published rectangle must contain 1.
    require(proof.count == 1 && arf_cmp_si(proof.endpoints[0].value(), 1) <= 0 &&
                arf_cmp_si(proof.endpoints[1].value(), 1) >= 0 &&
                arf_sgn(proof.endpoints[2].value()) <= 0 && arf_sgn(proof.endpoints[3].value()) >= 0,
            "source certificate misses the exact linear root");
    d::PiecewiseSourceCertificate moved(std::move(*certificate));
    require(!certificate->belongs_to(source) && moved.belongs_to(source),
            "moving a certificate left two source capabilities");
    require_consumed(*certificate);
    d::Root released = std::move(moved).release();
    require(released.count == 1 && !moved.belongs_to(source), "certificate release lost its proof or remained live");
    require_consumed(moved);

    cps_options cancelled = options;
    cancelled.cancel = [](void *) { return 1; };
    const auto before = statistics.finisher_evaluations;
    require(!source.certify(point.value(), cancelled, 256) && statistics.finisher_evaluations == before,
            "cancelled certification prepared or evaluated its source");
    acb_zero(point.value());
    require(!source.certify(point.value(), options, 128), "finisher accepted the separately handled exact zero point");
    acb_indeterminate(point.value());
    require(!source.certify(point.value(), options, 128), "nonfinite point produced a source certificate");
    acb_set_si(point.value(), 2);
    require(!source.certify(point.value(), options, 128), "root-free region produced a source certificate");

    // A double root cannot produce the singleton linear-model proof used here.
    coefficients.resize(3);
    coefficients[0].real = d::Number("1");
    coefficients[1].real = d::Number("-2");
    coefficients[2].real = d::Number("1");
    d::PiecewiseFinisherSource repeated(coefficients, other_statistics);
    acb_one(point.value());
    require(!repeated.certify(point.value(), options, 128), "multiple root produced a singleton source certificate");
}

int main() {
    certpolysolve::ThreadScope lifetime;
    try {
        source_certificates();
        std::vector<d::Coefficient> coefficients(17);
        for (size_t j = 0; j < coefficients.size(); ++j) {
            coefficients[j].real = d::Number(std::to_string(static_cast<int>(j % 7) - 3) + "/" +
                                             std::to_string(2 * j + 5));
            coefficients[j].imag = d::Number(std::to_string(static_cast<int>(j % 5) - 2) + "/" +
                                             std::to_string(2 * j + 7));
        }
        d::ComplexBall point, reference128, derivative128, reference256, derivative256;
        acb_set_d_d(point.value(), .5, .25);
        d::PiecewiseRootStatistics reference_statistics;
        d::PiecewiseFinisherSource reference(coefficients, reference_statistics);
        require(reference.evaluate(reference128.value(), derivative128.value(), point.value(), 128),
                "reference evaluation at 128 bits failed");
        require(reference.evaluate(reference256.value(), derivative256.value(), point.value(), 256),
                "reference evaluation at 256 bits failed");
        unsigned failures = 0, completed = 0;
        for (long fault = 0; fault < 16; ++fault) {
            d::PiecewiseRootStatistics statistics;
            d::PiecewiseFinisherSource source(coefficients, statistics);
            d::ComplexBall value, derivative;
            require(source.evaluate(value.value(), derivative.value(), point.value(), 128),
                    "initial cache preparation failed");
            bool failed = false;
            try {
                remaining = fault;
                const bool finite = source.evaluate(value.value(), derivative.value(), point.value(), 256);
                remaining = -1;
                require(finite, "successful replacement refused a finite point");
                ++completed;
            } catch (const std::bad_alloc &) {
                remaining = -1;
                failed = true;
                ++failures;
            }
            const double preparation_before_reuse = statistics.work.finisher_preparation;
            require(source.evaluate(value.value(), derivative.value(), point.value(), 128),
                    "old precision could not be reused after replacement");
            require(acb_equal(value.value(), reference128.value()) &&
                        acb_equal(derivative.value(), derivative128.value()),
                    "old-precision cache changed its enclosure");
            if (failed)
                require(statistics.work.finisher_preparation == preparation_before_reuse,
                        "failed replacement discarded the retained precision's cache");
            require(source.evaluate(value.value(), derivative.value(), point.value(), 256),
                    "replacement could not be retried");
            require(acb_equal(value.value(), reference256.value()) &&
                        acb_equal(derivative.value(), derivative256.value()),
                    "retried cache changed its enclosure");
        }
        require(failures >= 3 && completed > 0, "fault injection missed cache allocations or completion");
        std::cout << failures << " refused cache replacements; " << completed
                  << " replacements beyond the injected fault; old/new precision enclosures preserved\n";
    } catch (const std::exception &error) {
        remaining = -1;
        std::cerr << error.what() << '\n';
        return 1;
    }
}
