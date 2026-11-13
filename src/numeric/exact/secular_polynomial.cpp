// Native secular input as a polynomial (secular_polynomial): for residues a_j and distinct poles b_j, the exact
// numerator P = sum_j a_j Q / (x - b_j) - Q, Q = prod_j (x - b_j), whose roots solve sum_j a_j / (x - b_j) = 1, as a
// DensePolynomial of exact Gaussian rationals, of degree the number of terms and leading coefficient -1.
//
// Past the effort gate (cps_options::exact_effort_fits on 12 (n + 1)^2 Gaussian operations), one pass over the terms
// multiplies Q and the numerator by x - b_j (linear_product) and adds a_j times the old Q; Q is subtracted at the end.
// A second construction divides Q by each x - b_j exactly (synthetic division, remainder checked) and sums
// a_j Q / (x - b_j) - Q again; the two must agree in every coefficient (equal), else std::logic_error. Products and
// polynomials are held to cps_options::exact_storage_limit bits (multiply, budget); a refusal or a cancellation
// returns false with its reason in the diagnostic.
//
// State: none.
// Entry points: secular_polynomial, called by solver.cpp's SolveState::Impl::advance where native_advance does not run
// (a real-only query, EnginePolicy::backend off) or has fallen back. secular.hpp also declares SecularTerm and
// cps_secular, built by cps_secular_create (c_api.cpp, which requires nonzero residues and distinct poles).
// Not here: the native secular route (native_advance in solver.cpp, secular_refinement.cpp).
// Design notes: docs/mathematics.md ("Secular inputs and evaluation groundwork"), docs/api-spec-0.5.md (section 3.2),
// docs/resource-policy.md (the conversion's operation estimate).
#include "numeric/exact/secular_polynomial.hpp"
#include "numeric/exact/exact_complex.hpp"
#include <stdexcept>

namespace certpolysolve::detail {
namespace {
struct ExpansionLimit {};
uint64_t bits(const GaussianRational &x) {
    return resources::add(resources::add(fmpz_bits(fmpq_numref(x.real.value())), fmpz_bits(fmpq_denref(x.real.value()))),
                          resources::add(fmpz_bits(fmpq_numref(x.imag.value())), fmpz_bits(fmpq_denref(x.imag.value()))));
}
void budget(const GaussianPolynomial &p, uint64_t maximum) {
    uint64_t used = 0;
    for (const auto &x : p) {
        const auto cost = bits(x);
        if (cost > maximum - used)
            throw ExpansionLimit{};
        used += cost;
    }
}
GaussianRational add(const GaussianRational &a, const GaussianRational &b, bool subtract = false) {
    GaussianRational r;
    if (subtract) {
        fmpq_sub(r.real.value(), a.real.value(), b.real.value());
        fmpq_sub(r.imag.value(), a.imag.value(), b.imag.value());
    } else {
        fmpq_add(r.real.value(), a.real.value(), b.real.value());
        fmpq_add(r.imag.value(), a.imag.value(), b.imag.value());
    }
    return r;
}
GaussianRational multiply(const GaussianRational &a, const GaussianRational &b, uint64_t maximum) {
    if (bits(a) > maximum || bits(b) > maximum || resources::add(bits(a), bits(b)) == UINT64_MAX || resources::add(bits(a), bits(b)) > maximum)
        throw ExpansionLimit{};
    GaussianRational r;
    Rational t;
    fmpq_mul(r.real.value(), a.real.value(), b.real.value());
    fmpq_mul(t.value(), a.imag.value(), b.imag.value());
    fmpq_sub(r.real.value(), r.real.value(), t.value());
    fmpq_mul(r.imag.value(), a.real.value(), b.imag.value());
    fmpq_mul(t.value(), a.imag.value(), b.real.value());
    fmpq_add(r.imag.value(), r.imag.value(), t.value());
    if (bits(r) > maximum)
        throw ExpansionLimit{};
    return r;
}
GaussianPolynomial linear_product(const GaussianPolynomial &p, const GaussianRational &pole,
                                  uint64_t maximum) {
    GaussianPolynomial out(p.size() + 1);
    for (size_t j = 0; j < p.size(); ++j) {
        out[j] = add(out[j], multiply(p[j], pole, maximum), true);
        out[j + 1] = add(out[j + 1], p[j]);
    }
    budget(out, maximum);
    return out;
}
bool equal(const GaussianRational &a, const GaussianRational &b) {
    return fmpq_equal(a.real.value(), b.real.value()) && fmpq_equal(a.imag.value(), b.imag.value());
}
} // namespace

bool secular_polynomial(DensePolynomial &out, const cps_secular &source, const cps_options &options,
                        std::string &diagnostic) {
    try {
        const uint64_t n = source.terms.size();
        // Two streaming products plus the independent synthetic reconstruction:
        // conservative Gaussian add/multiply count, not a bit-complexity bound.
        if (!options.exact_effort_fits(resources::multiply(12, resources::multiply(n + 1, n + 1)))) {
            diagnostic = "secular exact numerator effort estimate exceeds caller allowance";
            return false;
        }
        if (options.exact_storage_limit() < 8 || source.terms.size() >= GaussianPolynomial().max_size())
            throw ExpansionLimit{};
        GaussianPolynomial product(1), numerator(1);
        fmpq_one(product[0].real.value());
        for (const auto &term : source.terms) {
            if (options.cancelled()) {
                diagnostic = "secular conversion cancelled";
                return false;
            }
            auto next = linear_product(numerator, term.pole, options.exact_storage_limit());
            for (size_t j = 0; j < product.size(); ++j)
                next[j] = add(next[j], multiply(term.residue, product[j], options.exact_storage_limit()));
            numerator = std::move(next);
            product = linear_product(product, term.pole, options.exact_storage_limit());
            budget(numerator, options.exact_storage_limit());
        }
        for (size_t j = 0; j < product.size(); ++j)
            numerator[j] = add(numerator[j], product[j], true);
        budget(numerator, options.exact_storage_limit());

        // Independently reconstruct sum a_i * Q/(x-b_i) - Q by exact synthetic
        // division of the finished denominator. This differs from the streaming
        // numerator recurrence above and catches sign/index/copy errors.
        GaussianPolynomial reconstructed(product.size());
        for (size_t j = 0; j < product.size(); ++j)
            reconstructed[j] = add(reconstructed[j], product[j], true);
        for (const auto &term : source.terms) {
            if (options.cancelled()) {
                diagnostic = "secular reconstruction cancelled";
                return false;
            }
            GaussianPolynomial quotient(product.size() - 1);
            quotient.back() = product.back();
            for (size_t j = quotient.size() - 1; j > 0; --j)
                quotient[j - 1] = add(product[j], multiply(term.pole, quotient[j], options.exact_storage_limit()));
            if (!add(product[0], multiply(term.pole, quotient[0], options.exact_storage_limit())).zero())
                throw std::logic_error("secular denominator division was not exact");
            for (size_t j = 0; j < quotient.size(); ++j)
                reconstructed[j] =
                    add(reconstructed[j], multiply(term.residue, quotient[j], options.exact_storage_limit()));
            budget(reconstructed, options.exact_storage_limit());
        }
        for (size_t j = 0; j < numerator.size(); ++j)
            if (!equal(numerator[j], reconstructed[j]))
                throw std::logic_error("secular numerator identity failed");
        DensePolynomial proposed;
        proposed.coefficients.resize(numerator.size());
        for (size_t j = 0; j < numerator.size(); ++j) {
            fmpq_set(proposed.coefficients[j].real.significand.value(), numerator[j].real.value());
            fmpq_set(proposed.coefficients[j].imag.significand.value(), numerator[j].imag.value());
        }
        out = std::move(proposed);
        diagnostic = "exact secular numerator; reconstructed identity verified";
        return true;
    } catch (ExpansionLimit &) {
        diagnostic = "secular exact numerator exceeds the supported expansion budget";
        return false;
    }
}
} // namespace certpolysolve::detail
