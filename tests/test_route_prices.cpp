// The selector's prices of the certificate attempt by kernel variant (core/heuristics.hpp, AttemptPrices;
// platform/kernel_variant.hpp, price_variant), and the pass limit of the attempt that a variant's prices put before
// the piecewise method (machine_newton_cover's pass_limit; dev/docs/instruction-set-prices-20261117.md):
//   - the reference table's price is four times the price before the tables, exactly, over degrees, supports and
//     sparse products on both sides of every case of the formula, and a variant's price is never above it;
//   - the variant whose prices are read is the kernels' variant, and the pin holds it at the reference;
//   - the pass limit: an attempt whose compensated tier passes on more points than the limit ends there, before a
//     repair, a cluster start, a ball certificate or a round; with the limit at or above what the tier passes on,
//     it does what it does without a limit;
//   - the packet's charge: with an allowance that the packet's price in the reference code uses up, the attempt
//     has nothing left for its double-double repair; charged at the AVX2 table's price it has, and repairs.
// The order rule itself, under each table and under the pin, is in test_auto_selection.cpp.
#include "certify/attempts/newton_cover.hpp"
#include "core/engine_policy.hpp"
#include "core/engine_statistics.hpp"
#include "core/heuristics.hpp"
#include "core/query.hpp" // cps_options
#include "core/raii.hpp"
#include "numeric/evaluate/structured_source.hpp"
#include "platform/kernel_variant.hpp"
#include <certpolysolve/certpolysolve.hpp>
#include <flint/acb_poly.h>
#include <complex>
#include <cstdint>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace d = certpolysolve::detail;
namespace h = certpolysolve::detail::heuristics;

namespace {
void require(bool condition, const std::string &message) {
    if (!condition)
        throw std::runtime_error(message);
}

void table_tests() {
    require(h::attempt_price_quarters == 4, "the tables are not in quarters of an operation");
    uint64_t dense = 0, sparse = 0, direct = 0, tree = 0, cheaper = 0;
    const std::vector<uint64_t> degrees = {1,    2,    7,    12,    100,   1000,  2047,
                                           2048, 4096, 8192, 16384, 32768, 65536, 1000000};
    for (const uint64_t degree : degrees)
        for (const uint64_t support : std::vector<uint64_t>{1, 2, 3, degree / 8 + 1, degree / 4, degree / 4 + 1,
                                                            degree / 2 + 1, degree + 1})
            for (const uint64_t products : std::vector<uint64_t>{0, 5, 44}) {
                if (support > degree + 1 || support == 0)
                    continue;
                const uint64_t before = h::certificate_attempt_point_operations(degree, support, products);
                const uint64_t reference =
                    h::certificate_attempt_point_quarters(h::reference_attempt_prices, degree, support, products);
                const uint64_t variant =
                    h::certificate_attempt_point_quarters(h::avx2_fma_attempt_prices, degree, support, products);
                require(reference == h::attempt_price_quarters * before,
                        "the reference table's price is not four times the price before the tables");
                require(variant <= reference, "a variant's price is above the reference's");
                (4 * support <= degree ? sparse : dense) += 1;
                (degree >= h::machine_interaction_tree_degree ? tree : direct) += 1;
                cheaper += variant < reference;
            }
    require(dense && sparse && direct && tree && cheaper, "the grid misses a case of the price's formula");
    // A dense source at a tree degree: the evaluation at a quarter, the tree as it is, the certificate at 2 for 14.
    const uint64_t n = 16384, levels = 15;
    require(h::certificate_attempt_point_quarters(h::avx2_fma_attempt_prices, n, n + 1) ==
                h::machine_cover_sweeps * ((n + 1) + 4 * h::machine_tree_level_weight * levels) + 8 * (n + 1),
            "the AVX2 table's price of a dense source is not its three weights");
    // The pass limit: the points whose ball certificates cost what the tier's pass over every point cost.
    require(h::attempt_pass_limit(h::avx2_fma_attempt_prices, 16384) == 16384 / 34 &&
                h::attempt_pass_limit(h::avx2_fma_attempt_prices, 32768) == 32768 / 34,
            "the AVX2 table's pass limit is not a thirty-fourth of the degree");
}

void pin_tests() {
    const d::KernelVariant variant_before = d::kernel_variant();
    const bool pin_before = d::set_reference_prices(false);
    require(d::price_variant() == d::kernel_variant(), "free prices are not the kernels' variant's");
    d::set_kernel_variant(d::KernelVariant::Reference);
    require(d::price_variant() == d::KernelVariant::Reference, "the reference kernels with another variant's prices");
    if (d::kernel_variant_available(d::KernelVariant::Avx2Fma)) {
        d::set_kernel_variant(d::KernelVariant::Avx2Fma);
        require(d::price_variant() == d::KernelVariant::Avx2Fma, "the AVX2 kernels without their prices");
        require(!d::set_reference_prices(true), "set_reference_prices does not return the pin before");
        require(d::price_variant() == d::KernelVariant::Reference && d::kernel_variant() == d::KernelVariant::Avx2Fma,
                "the pin does not hold the prices at the reference's, or it moved the kernels");
        require(d::set_reference_prices(false), "set_reference_prices does not return the pin before");
        require(d::price_variant() == d::KernelVariant::Avx2Fma, "the prices stay pinned after the pin is freed");
    }
    d::set_kernel_variant(variant_before);
    d::set_reference_prices(pin_before);
}

// The monic polynomial with the given roots, exactly (dyadic coordinates, a precision far above their products').
d::ComplexPolynomial from_roots(const std::vector<std::complex<double>> &roots) {
    const slong n = static_cast<slong>(roots.size());
    d::RootVector points(n);
    for (slong i = 0; i < n; ++i)
        acb_set_d_d(points[i], roots[static_cast<size_t>(i)].real(), roots[static_cast<size_t>(i)].imag());
    d::ComplexPolynomial p;
    acb_poly_product_roots(p.value(), points.data(), n, 16384);
    for (slong k = 0; k <= n; ++k)
        require(acb_is_exact(acb_poly_get_coeff_ptr(p.value(), k)), "a test polynomial's coefficient is not exact");
    return p;
}

struct Outcome {
    bool certified = false;
    std::vector<d::Root> cover;
    d::PartialNewtonCover partial;
    d::StudyStatistics counts;
};

// The attempt as the piecewise recipe calls it, on a source without real structure (no conjugate pairs: the points
// the tier passes on are then the points it refuses), with the given pass limit.
Outcome attempt(const d::ComplexPolynomial &polynomial, uint64_t pass_limit, uint64_t fallback_operations = 0,
                const h::AttemptPrices *packet_prices = nullptr) {
    const slong degree = acb_poly_degree(polynomial.value());
    d::PolynomialSourceOptions source_options;
    source_options.input = d::InputSemantics::Exact;
    d::PolynomialSource source(
        [&](acb_poly_t out, slong precision) {
            d::ComplexBall a;
            acb_poly_zero(out);
            for (slong k = 0; k <= degree; ++k) {
                acb_set_round(a.value(), acb_poly_get_coeff_ptr(polynomial.value(), k), precision);
                acb_poly_set_coeff_acb(out, k, a.value());
            }
        },
        degree, false, source_options);
    cps_options options;
    options.goal = CPS_ROOT_APPROXIMATION;
    options.bits = 53;
    d::EnginePolicy policy;
    policy.study_statistics = std::make_shared<d::StudyStatistics>();
    Outcome out;
    d::CandidateStatistics work{};
    out.certified = d::machine_newton_cover(out.cover, source, nullptr, options, 85, true, work, true, &policy, nullptr,
                                            nullptr, nullptr, true, &out.partial, {}, fallback_operations, pass_limit,
                                            packet_prices);
    out.counts = *policy.study_statistics;
    return out;
}

uint64_t held(const Outcome &outcome) {
    uint64_t count = 0;
    for (const uint8_t flag : outcome.partial.certified)
        count += flag;
    return count;
}

// What an attempt did after its tier's first pass.
uint64_t later_work(const d::StudyStatistics &counts) {
    return counts.machine_cover_double_double_repairs + counts.machine_cover_cluster_starts + counts.machine_cover_rounds +
           counts.machine_cover_newton_steps + counts.machine_cover_precision_retries +
           counts.machine_cover_rung_evaluations + counts.machine_cover_compensated_later_tries;
}

void pass_limit_tests() {
    // Twelve well-separated roots and three double roots: the tier proves the twelve and passes the six on.
    std::vector<std::complex<double>> roots = {{0.0, 1.0}, {1.0, 1.0},  {-2.0, 1.0}, {0.5, 0.5},  {2.0, 3.0},  {-1.0, -2.0},
                                               {3.0, 0.5}, {-3.0, 1.5}, {1.5, -2.5}, {-0.5, -3.0}, {4.0, -1.0}, {-4.0, 2.0}};
    const size_t separated = roots.size();
    for (const std::complex<double> z : {std::complex<double>{0.75, 0.25}, {-1.25, 0.5}, {2.5, -1.75}}) {
        roots.push_back(z);
        roots.push_back(z);
    }
    const d::ComplexPolynomial p = from_roots(roots);
    const size_t n = roots.size(), passed_on = n - separated;
    const Outcome free = attempt(p, 0);
    require(!free.certified, "the fixture with double roots was certified: it tests nothing");
    require(later_work(free.counts) > 0, "the attempt without a limit did nothing after its tier's pass");
    // Under the limit the attempt ends at the tier's pass: no later work, and the tier's certificates are held.
    const Outcome limited = attempt(p, passed_on - 1);
    require(!limited.certified, "an attempt over its pass limit returned a cover");
    require(later_work(limited.counts) == 0, "an attempt over its pass limit went on after the tier's pass");
    require(held(limited) == separated, "an attempt over its pass limit does not hold the tier's certificates");
    require(limited.counts.machine_cover_compensated_certificates == separated,
            "the tier did not prove the separated roots before the limit ended the attempt");
    // At the limit, and far above it, the attempt is the attempt without a limit.
    for (const uint64_t limit : {static_cast<uint64_t>(passed_on), static_cast<uint64_t>(n), uint64_t{1} << 40}) {
        const Outcome same = attempt(p, limit);
        require(same.certified == free.certified && held(same) == held(free) &&
                    same.partial.certified == free.partial.certified &&
                    same.partial.certificate_bits == free.partial.certificate_bits &&
                    later_work(same.counts) == later_work(free.counts) &&
                    same.counts.machine_cover_compensated_certificates == free.counts.machine_cover_compensated_certificates,
                "an attempt within its pass limit differs from the attempt without a limit");
    }
    // An input the attempt finishes is finished under any limit that its tier's pass stays within.
    roots.resize(separated);
    const d::ComplexPolynomial easy = from_roots(roots);
    const Outcome whole = attempt(easy, 0), bounded = attempt(easy, 1);
    require(whole.certified && bounded.certified && whole.cover.size() == separated && bounded.cover.size() == separated,
            "the separated roots are not certified, with a pass limit or without");
}
// The packet's charge. Sixty-two roots on a grid and one pair 2^-45 apart: the tier passes the pair on, and the
// double-double repair is what the attempt has for it. The allowance is the packet's one quantum at the reference's
// price, so that charge leaves nothing; the AVX2 table's charge leaves the difference of the two prices.
void packet_charge_tests() {
    require(h::machine_packet_point_quarters(h::reference_attempt_prices, 64, 65) ==
                h::attempt_price_quarters * h::machine_packet_point_operations(64, 65),
            "the reference table's packet is not the packet's price before the tables");
    std::vector<std::complex<double>> roots;
    for (int a = -4; a <= 3; ++a)
        for (int b = -4; b <= 3; ++b)
            if (roots.size() < 62)
                roots.emplace_back(a + 0.5, b + 0.25);
    roots.emplace_back(0.125, 0.0625);
    roots.emplace_back(0.125 + 0x1p-45, 0.0625);
    const d::ComplexPolynomial p = from_roots(roots);
    const uint64_t n = roots.size();
    const uint64_t quantum = h::machine_packet_point_operations(n, n + 1) * n;
    require(h::machine_packet_point_quarters(h::avx2_fma_attempt_prices, n, n + 1) * n / h::attempt_price_quarters < quantum,
            "the AVX2 table's packet is not below the reference's");
    const Outcome reference = attempt(p, 0, quantum), variant = attempt(p, 0, quantum, &h::avx2_fma_attempt_prices);
    require(reference.counts.machine_cover_compensated_certificates >= 62 &&
                variant.counts.machine_cover_compensated_certificates >= 62,
            "the tier did not prove the roots on the grid");
    require(reference.counts.machine_cover_double_double_repairs == 0,
            "an attempt whose packet used its allowance up still repaired");
    require(variant.counts.machine_cover_double_double_repairs >= 1,
            "an attempt whose packet is charged at the AVX2 table's price did not repair with what that leaves");
}
} // namespace

int main() {
    certpolysolve::ThreadScope lifetime;
    try {
        table_tests();
        pin_tests();
        pass_limit_tests();
        packet_charge_tests();
        std::cout << "route prices tests passed\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "FAILED: " << error.what() << '\n';
        return 1;
    }
}
