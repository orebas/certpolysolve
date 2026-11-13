#include "numeric/arithmetic/polynomial_arithmetic.hpp"
#include "numeric/exact/modular_squarefree.hpp"
#include "exact_winding_api.hpp"
#include "methods/piecewise/piecewise_admission.hpp"
#include <flint/ulong_extras.h>
#include <iostream>
#include <limits>

namespace d = certpolysolve::detail;
using namespace reference;

static void check() {
    require(n_is_prime(998244353), "modular admission modulus is not prime");
    cps_options options;
    options.max_work_bits = 4096;
    auto simple = reference::input(from_roots({Complex(1, 2), Complex(3), Complex(-2, 1)}));
    auto repeated = reference::input(from_roots({Complex(1, 2), Complex(1, 2), Complex(-2, 1)}));
    require(d::nonzero_part_squarefree_modular(certpolysolve::detail::dense(*simple.get()), options),
            "simple Q(i) polynomial not proved squarefree");
    require(!d::nonzero_part_squarefree_modular(certpolysolve::detail::dense(*repeated.get()), options),
            "modular proof accepted a repeated Q(i) root");
    // A proof refused by the allowance is no proof: four coefficients need four units.
    cps_options below;
    below.max_work_bits = 4096;
    below.exact_operations = 3;
    require(!d::nonzero_part_squarefree_modular(certpolysolve::detail::dense(*simple.get()), below),
            "a modular proof refused by exact_operations was reported as a proof");
    certpolysolve::Polynomial bad_denominator, vanished_leading, characteristic;
    bad_denominator.set(0, certpolysolve::Number("1/998244353"));
    bad_denominator.set(1, certpolysolve::Number("1"));
    vanished_leading.set(0, certpolysolve::Number("1"));
    vanished_leading.set(1, certpolysolve::Number("998244353"));
    characteristic.set(0, certpolysolve::Number("998244353"));
    characteristic.set(2, certpolysolve::Number("1"));
    require(!d::nonzero_part_squarefree_modular(certpolysolve::detail::dense(*bad_denominator.get()), options) &&
                !d::nonzero_part_squarefree_modular(certpolysolve::detail::dense(*vanished_leading.get()), options) &&
                !d::nonzero_part_squarefree_modular(certpolysolve::detail::dense(*characteristic.get()), options),
            "bad reductions must be inconclusive, including squarefree characteristic-zero inputs");
    certpolysolve::Polynomial compressed;
    compressed.set(0, certpolysolve::Number("2@-2000000000000000000000000000000"));
    compressed.set(1, certpolysolve::Number("-3@-1000000000000000000000000000000"));
    compressed.set(2, certpolysolve::Number("1"));
    require(d::nonzero_part_squarefree_modular(certpolysolve::detail::dense(*compressed.get()), options),
            "compressed exponents were not reduced modulo p-1");

    // Same structure under exact common and variable power-of-two scaling.
    // Admission needs no benchmark identifiers and makes no root-count claim.
    certpolysolve::Polynomial dense, scaled, bowed, composition;
    constexpr size_t n = 16384;
    for (size_t j = 0; j <= n; ++j) {
        dense.set(j, certpolysolve::Number("1"));
        scaled.set(j, certpolysolve::Number("1@" + std::to_string(1024 * static_cast<int64_t>(j) - 9000)));
        bowed.set(j, certpolysolve::Number(j == n / 2 ? "1@32" : "1"));
        if (!(j % 2))
            composition.set(j, certpolysolve::Number("1"));
    }
    const auto ordinary = d::admit_piecewise(certpolysolve::detail::dense(*dense.get()), options);
    const auto transformed = d::admit_piecewise(certpolysolve::detail::dense(*scaled.get()), options);
    require(ordinary.selected && transformed.selected && ordinary.squarefree_proved &&
                transformed.squarefree_proved,
            "cheap balanced dense case was not admitted under exact scaling");
    require(ordinary.envelope_bend == transformed.envelope_bend && ordinary.support == transformed.support &&
                ordinary.model_bits == transformed.model_bits,
            "structural features changed under exact common/variable scaling");
    require(!d::admit_piecewise(certpolysolve::detail::dense(*bowed.get()), options).selected, "broad radial geometry was admitted");
    require(!d::admit_piecewise(certpolysolve::detail::dense(*composition.get()), options).selected,
            "selector intercepted an existing power-composition route");
    // A cheap rejection must not allocate one ring per unit of an enormous
    // exponent gap. This input would require infeasible work in that order.
    certpolysolve::Polynomial enormous_gap;
    enormous_gap.set(0, certpolysolve::Number("1"));
    enormous_gap.set(1, certpolysolve::Number("1"));
    enormous_gap.set(2, certpolysolve::Number("1@1000000000000"));
    enormous_gap.set(3, certpolysolve::Number("1"));
    const auto cheap_refusal = d::admit_piecewise(certpolysolve::detail::dense(*enormous_gap.get()), options);
    require(!cheap_refusal.selected && std::string(cheap_refusal.reason) == "broad_radial_geometry" &&
                cheap_refusal.estimated_models == 0 && cheap_refusal.envelope_bend > 1000000000,
            "wide-envelope guard did expensive chart work before declining");
    certpolysolve::Polynomial low;
    for (size_t j = 0; j <= 256; ++j)
        low.set(j, certpolysolve::Number("1"));
    require(!d::admit_piecewise(certpolysolve::detail::dense(*low.get()), options).selected,
            "small input admitted despite predicted overhead");
    cps_options real = options;
    real.real_only = true;
    require(!d::admit_piecewise(certpolysolve::detail::dense(*dense.get()), real).selected, "unsupported real-only query intercepted");

    {
        d::PiecewiseAdmissionPolicy admission;
        admission.minimum_degree = 2049;
        admission.require_attempt_fits_call = true;
        auto tiny = options;
        tiny.work_units = 1;
        const auto refusal = d::admit_piecewise(certpolysolve::detail::dense(*dense.get()), tiny, admission);
        require(!refusal.selected && std::string_view(refusal.reason) == "remaining_call_work" &&
                    refusal.estimated_attempt_work > tiny.work_units && !refusal.squarefree_proved,
                "admission ignored the remaining caller allowance");
        for (size_t degree : {size_t(2048), size_t(2049)}) {
            certpolysolve::Polynomial below;
            // Deflation lowers the effective degree to two. The huge exponent
            // must never reach enclosure/model estimation in this probe.
            below.set(degree - 2, certpolysolve::Number("1"));
            below.set(degree - 1, certpolysolve::Number("1@1000000000000"));
            below.set(degree, certpolysolve::Number("1"));
            const auto reduced = d::admit_piecewise(certpolysolve::detail::dense(*below.get()), options, admission);
            require(!reduced.selected && std::string_view(reduced.reason) == "small_degree" &&
                        reduced.degree == 2 && reduced.estimated_models == 0,
                    "zero deflation did not protect small factors from model preparation");
        }
    }
}

int main() {
    certpolysolve::ThreadScope lifetime;
    try {
        check();
        std::cout << "structural selection and modular proof passed\n";
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
