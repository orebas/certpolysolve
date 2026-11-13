#include <certpolysolve/certpolysolve.hpp>
#include "numeric/candidates/candidates.hpp"
#include "core/raii.hpp"
#include <iostream>
#include <random>

// Exact variable rescaling x -> 2^s x must rescale the initial candidates
// exactly. Otherwise equivalent inputs start different numerical trajectories;
// near-multiple degree-769 then took 32 s on q(x) but 572 s on q(8x).
namespace d = certpolysolve::detail;

namespace {
bool rescaled_seeds_match(const d::ComplexPolynomial &polynomial, slong shift, slong precision) {
    const slong n = acb_poly_degree(polynomial.value());
    d::ComplexPolynomial scaled;
    acb_poly_set(scaled.value(), polynomial.value());
    for (slong k = 0; k <= n; ++k)
        acb_mul_2exp_si(scaled.value()->coeffs + k, scaled.value()->coeffs + k, shift * k);
    d::RootVector original(n), rescaled(n);
    if (!d::initial_candidates(original.data(), polynomial.value(), precision) ||
        !d::initial_candidates(rescaled.data(), scaled.value(), precision))
        return false;
    for (slong i = 0; i < n; ++i) {
        acb_mul_2exp_si(rescaled.data() + i, rescaled.data() + i, shift);
        if (!acb_equal(original.data() + i, rescaled.data() + i)) return false;
    }
    return true;
}
} // namespace

int main() {
    certpolysolve::ThreadScope lifetime;
    std::mt19937_64 random(0x5ca1e);
    // Random magnitudes give Newton-polygon edges whose lengths do not divide
    // their height drops, where rounding drop/count used to break the property.
    for (unsigned trial = 0; trial < 40; ++trial) {
        const slong n = 20 + static_cast<slong>(random() % 180);
        d::ComplexPolynomial polynomial;
        acb_poly_fit_length(polynomial.value(), n + 1);
        for (slong k = 0; k <= n; ++k) {
            const auto mantissa = static_cast<slong>(random() % 2001) - 1000;
            acb_set_si(polynomial.value()->coeffs + k, mantissa ? mantissa : 1);
            acb_mul_2exp_si(polynomial.value()->coeffs + k, polynomial.value()->coeffs + k,
                            static_cast<slong>(random() % 400) - 200);
        }
        _acb_poly_set_length(polynomial.value(), n + 1);
        for (const slong shift : {1, 3, -2, 7, -11}) {
            if (!rescaled_seeds_match(polynomial, shift, 128)) {
                std::cerr << "initial candidates are not exactly equivariant under x -> 2^" << shift
                          << " x (trial " << trial << ", degree " << n << ")\n";
                return 1;
            }
        }
    }
    return 0;
}
