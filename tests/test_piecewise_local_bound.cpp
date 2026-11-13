#include "numeric/arithmetic/polynomial_arithmetic.hpp"
#include "exact_winding_api.hpp"
#include "methods/piecewise/piecewise_roots.hpp"
#include <iostream>

namespace d = certpolysolve::detail;
using namespace reference;

static d::RealBall direct_bound(const d::ComplexPolynomial &p, const d::RootModel &m,
                               const d::ComplexBall &rectangle) {
    d::Float lo, hi, coefficient;
    acb_get_abs_lbound_arf(lo.value(), rectangle.value(), 256);
    acb_get_abs_ubound_arf(hi.value(), rectangle.value(), 256);
    d::RealBall result, radius, term, power;
    for (slong j = 0; j <= acb_poly_degree(p.value()); ++j) {
        if (j >= m.lower_index && j <= m.upper_index) continue;
        acb_get_abs_ubound_arf(coefficient.value(), p.value()->coeffs + j, 256);
        arb_set_arf(term.value(), coefficient.value());
        const slong exponent = j - m.lower_index;
        arb_set_arf(radius.value(), exponent < 0 ? lo.value() : hi.value());
        arb_pow_ui(power.value(), radius.value(), static_cast<ulong>(std::abs(exponent)), 256);
        if (exponent < 0) arb_div(term.value(), term.value(), power.value(), 256);
        else arb_mul(term.value(), term.value(), power.value(), 256);
        arb_add(result.value(), result.value(), term.value(), 256);
    }
    arb_mul_2exp_si(result.value(), result.value(), -m.scale);
    return result;
}

static void whole_rectangle_and_lifetime() {
    d::ComplexPolynomial original;
    d::ComplexBall coefficient;
    for (slong j = 0; j <= 64; ++j) {
        acb_set_si_si(coefficient.value(), j % 5 - 2, j % 7 - 3);
        acb_mul_2exp_si(coefficient.value(), coefficient.value(), (j % 9) - 4);
        if (j == 0 || j == 64) acb_one(coefficient.value());
        else {
            // The envelope must include uncertain coefficients, including
            // balls that contain zero, not only their midpoints.
            arb_add_error_2exp_si(acb_realref(coefficient.value()), -20);
            arb_add_error_2exp_si(acb_imagref(coefficient.value()), -19);
        }
        acb_poly_set_coeff_acb(original.value(), j, coefficient.value());
    }
    d::RootModel retained;
    {
        d::ComplexPolynomial input;
        acb_poly_set(input.value(), original.value());
        d::PiecewiseRootStatistics statistics;
        d::PiecewiseRootConfig config;
        config.local_band_error = true;
        d::PiecewiseModelCursor cursor(input.value(), 8, 192, statistics, {}, config);
        cps_options query;
        for (;;) {
            const auto step = cursor.advance(query);
            require(step != d::PiecewiseModelStep::Refused && step != d::PiecewiseModelStep::Done,
                    "test polynomial supplied no model");
            if (step == d::PiecewiseModelStep::Model) break;
        }
        retained = cursor.model();
        require(bool(retained.error_envelope), "model lacks an owning envelope");
        acb_poly_zero(input.value());
    }
    // The copied model remains usable after both input and cursor destruction.
    // Test both tails, each individual tail, and exact empty omitted bands.
    for (const auto [lower, upper] : {std::pair<slong, slong>{7, 49}, {0, 49}, {7, 64}, {0, 64}}) {
        retained.lower_index = lower;
        retained.upper_index = upper;
        retained.scale = 11;
        for (const auto [real, imag] : {std::pair<double, double>{1, 0}, {0, 2}, {-1, 2}, {0.125, -0.25}}) {
            d::ComplexBall rectangle;
            acb_set_d_d(rectangle.value(), real, imag);
            arb_add_error_2exp_si(acb_realref(rectangle.value()), -9);
            arb_add_error_2exp_si(acb_imagref(rectangle.value()), -8);
            const auto bound = d::local_omitted_band_bound(retained, rectangle.value());
            require(bool(bound), "finite origin-excluding rectangle refused");
            const auto direct = direct_bound(original, retained, rectangle);
            // Shared stopping/proof arithmetic retains the unchanged Taylor
            // error and uses the global bound unless provably smaller.
            retained.error = *bound;
            arb_mul_2exp_si(retained.error.value(), retained.error.value(), 2);
            arb_set_d(retained.taylor_model_error.value(), 0.125);
            arb_add(retained.error.value(), retained.error.value(), retained.taylor_model_error.value(), 256);
            const auto combined = d::root_model_error(retained, rectangle.value(), 256);
            d::RealBall expected = direct;
            arb_add(expected.value(), expected.value(), retained.taylor_model_error.value(), 256);
            require(arb_ge(combined.value(), expected.value()), "combined error lost a remainder component");
            require((arb_is_zero(bound->value()) && arb_is_zero(direct.value())) ||
                    arb_ge(bound->value(), direct.value()), "local envelope underestimates direct whole-box sum");
        }
    }
    d::ComplexBall origin;
    arb_add_error_2exp_si(acb_realref(origin.value()), -100);
    arb_add_error_2exp_si(acb_imagref(origin.value()), -100);
    retained.lower_index = 1;
    require(!d::local_omitted_band_bound(retained, origin.value()), "negative powers admitted the origin");
    retained.lower_index = 0;
    retained.upper_index = 64;
    const auto empty = d::local_omitted_band_bound(retained, origin.value());
    require(empty && arb_is_zero(empty->value()), "empty bands acquired a nonzero error");
    retained.error_envelope.reset();
    require(!d::local_omitted_band_bound(retained, origin.value()), "missing snapshot invented a bound");
}

static void original_polynomial_certificates() {
    const std::vector<Complex> roots{{1}, {2}, {-1, 2}, {2, -3}, {-2, -1}, {1, 4}};
    for (const auto &scale : {Complex(1), Complex(Rational(1)/8), Complex(0, 2)})
        for (const auto &shift : {Complex(0), Complex(2, -1)}) {
            std::vector<Complex> transformed;
            for (const auto &z : roots) transformed.push_back(z * scale + shift);
            const auto polynomial = from_roots(transformed);
            auto input = reference::input(polynomial);
            d::PiecewiseRootConfig config;
            config.local_band_error = true;
            config.aligned_polish = config.reuse_polish = true;
            config.error_components = true;
            config.durable_ownership = config.reuse_guards = config.recover_discovery = true;
            d::PiecewiseRootStatistics statistics;
            cps_options query;
            query.goal = CPS_ROOT_APPROXIMATION;
            const auto result = d::solve_piecewise_roots(certpolysolve::detail::dense(*input.get()), query, statistics, config);
            require(result.goal_met && d::disjoint(result.roots), "local-bound solve lost coverage");
            uint64_t count = 0;
            for (const auto &root : result.roots) {
                auto endpoint = [](const d::Endpoint &e) {
                    return reference::endpoint(certpolysolve::Dyadic{e.mantissa_text(), e.exponent_text()});
                };
                require(rectangle_count(polynomial, endpoint(root.endpoints[0]), endpoint(root.endpoints[1]),
                          endpoint(root.endpoints[2]), endpoint(root.endpoints[3])) == root.count,
                        "local-bound certificate failed independent exact winding");
                count += root.count;
            }
            require(count == transformed.size() && statistics.certificate.model_error.local_band_attempts > 0,
                    "local-bound test missed the certificate path");
        }
}

int main() {
    certpolysolve::ThreadScope lifetime;
    try {
        whole_rectangle_and_lifetime();
        original_polynomial_certificates();
        std::cout << "local band bounds, uncertainty, lifetime and original-input counts passed\n";
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
