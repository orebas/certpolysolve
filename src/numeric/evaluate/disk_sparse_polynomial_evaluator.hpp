#pragma once

#include "core/raii.hpp"
#include "core/resource_limits.hpp"
#include <bit>
#include <span>
#include <stdexcept>
#include <vector>

namespace certpolysolve::detail {
// Value and slope of a polynomial with few nonzero coefficients, with disk
// errors: the sparse sibling of DiskBlockedPolynomialEvaluator.
//
// p(z) = sum_j c_j z^(e_j) over the support e_0 < e_1 < ... The powers
// W_j = z^(e_j) are built along the support, W_j = W_(j-1) z^(e_j - e_(j-1)),
// and a gap's power comes from a table of z^(d 2^(w i)) (the gap's digits d of
// w bits): one multiplication for each nonzero digit. Then two FLINT dot
// products, p(z) = c . W and z p'(z) = (e c) . W.
//
// Every product is kept as a midpoint, a disk radius and a bound of the
// midpoint's size,
//     |a b - a^ b^| <= E_a (|b^| + E_b) + |a^| E_b,   plus the product's rounding,
// so no error is re-boxed into a rotated rectangle. Ball Horner over the terms
// does that at every term and loses up to half a bit each time (a factor
// |cos t| + |sin t| at argument t).
// History: docs/code-history-20261021.md, disk_sparse_polynomial_evaluator.hpp 1.
//
// A power z^e built from products carries e rounding errors whatever chain
// builds it (a squaring doubles the relative error it is given). The powers
// are therefore computed with bit_width(degree) guard bits: their error is
// then that of one rounding at the caller's precision, as the dot products'
// own is.
//
// This object owns its coefficients and scratch. Calls on the same object
// must be exclusive.
class DiskSparsePolynomialEvaluator {
  public:
    // How the powers are built: digits of `width` bits, `windows` of them, and
    // the complex multiplications one evaluation makes (the table's and the
    // chain's). The width is the one with the fewest multiplications for this
    // support; any width is correct.
    struct Plan {
        unsigned width = 1;
        slong windows = 0;
        uint64_t multiplications = 0;
    };
    // `exponents` ascending and nonnegative: the support.
    static Plan plan(std::span<const slong> exponents) {
        slong largest = 0, previous = 0;
        for (const slong e : exponents) {
            if (e < previous)
                throw std::invalid_argument("sparse support must be ascending and nonnegative");
            largest = std::max(largest, e - previous);
            previous = e;
        }
        Plan best;
        const auto bits = static_cast<unsigned>(std::bit_width(static_cast<ulong>(largest)));
        if (!bits)
            return best;
        bool found = false;
        for (unsigned width = 1; width <= bits; ++width) {
            const uint64_t windows = (bits + width - 1) / width;
            // A window holds z^(d 2^(w i)) for d = 1 .. 2^w - 1: its first entry
            // is the base it was handed, every further entry one product, and
            // one more product hands the next window its base.
            const uint64_t entries = (uint64_t(1) << width) - 2;
            uint64_t count = resources::add(resources::multiply(windows, entries), windows - 1);
            if (found && count >= best.multiplications)
                break; // wider tables only cost more
            const ulong mask = (ulong(1) << width) - 1;
            previous = 0;
            for (const slong e : exponents) {
                for (auto gap = static_cast<ulong>(e - previous); gap; gap >>= width)
                    count += (gap & mask) != 0;
                previous = e;
            }
            if (!found || count < best.multiplications) {
                best = {width, static_cast<slong>(windows), count};
                found = true;
            }
        }
        return best;
    }

    DiskSparsePolynomialEvaluator(const ComplexPolynomial &polynomial, slong precision)
        : precision_(precision), exponents_(support(polynomial)), plan_(plan(exponents_)),
          coefficients_(terms()), slopes_(terms()), powers_(terms()), table_(table_size()) {
        if (precision < 2)
            throw std::invalid_argument("evaluation precision must be at least two bits");
        const slong degree = exponents_.empty() ? 0 : exponents_.back();
        power_precision_ =
            resources::precision_sum(precision, std::bit_width(static_cast<ulong>(degree)));
        second_.resize(exponents_.size());
        table_error_.resize(static_cast<size_t>(table_.size()));
        table_bound_.resize(static_cast<size_t>(table_.size()));
        const auto coefficients = polynomial.coefficients();
        for (slong j = 0; j < terms(); ++j) {
            const slong e = exponents_[static_cast<size_t>(j)];
            const acb_srcptr c = coefficients[e];
            finite_ = finite_ && acb_is_finite(c);
            acb_set(coefficients_[j], c);
            acb_mul_si(slopes_[j], c, e, ARF_PREC_EXACT);
            // e (e - 1) |c|: the term's share of the second derivative's majorant.
            auto &m = second_[static_cast<size_t>(j)];
            acb_get_mag(m.value(), c);
            mag_mul_ui(m.value(), m.value(), static_cast<ulong>(e));
            mag_mul_ui(m.value(), m.value(), static_cast<ulong>(e > 0 ? e - 1 : 0));
        }
    }
    DiskSparsePolynomialEvaluator(const DiskSparsePolynomialEvaluator &) = delete;
    DiskSparsePolynomialEvaluator &operator=(const DiskSparsePolynomialEvaluator &) = delete;

    slong terms() const { return static_cast<slong>(exponents_.size()); }
    const Plan &plan() const { return plan_; }
    // The complex multiplications made so far: plan().multiplications for each
    // evaluation (the selector's cost bound counts them; a test holds it to that).
    uint64_t multiplications() const { return multiplications_; }

    // Enclosures of p and p' at a finite point whose ball excludes zero (the
    // slope is recovered from z p'(z)).
    bool evaluate(acb_t value, acb_t derivative, const acb_t point) {
        // Invalid aliasing is refused before touching any caller-owned value.
        if (value == derivative || value == point || derivative == point)
            return false;
        if (!finite_ || !acb_is_finite(point) || acb_contains_zero(point))
            return false;
        build_powers(point);
        acb_dot(value, nullptr, 0, coefficients_.data(), 1, powers_.data(), 1, terms(), precision_);
        acb_dot(derivative, nullptr, 0, slopes_.data(), 1, powers_.data(), 1, terms(), precision_);
        acb_div(derivative, derivative, point, precision_);
        return acb_is_finite(value) && acb_is_finite(derivative);
    }
    // The value alone, at any finite point.
    bool evaluate_value(acb_t value, const acb_t point) {
        if (value == point || !finite_ || !acb_is_finite(point))
            return false;
        build_powers(point);
        acb_dot(value, nullptr, 0, coefficients_.data(), 1, powers_.data(), 1, terms(), precision_);
        return acb_is_finite(value);
    }


    // DiskPolynomialEvaluator::linear_remainder over the support: for
    // |z - c| <= r, Taylor's integral remainder is at most
    // (r^2 / 2) sum_{e_j >= 2} e_j (e_j - 1) |c_j| (|c| + r)^(e_j - 2).
    void linear_remainder(arb_t out, const acb_t center, const arb_t radius) const {
        Magnitude norm, r, bound, power;
        acb_get_mag(norm.value(), center);
        arb_get_mag(r.value(), radius);
        mag_add(norm.value(), norm.value(), r.value());
        // Horner over the support, from the highest term down to exponent 2.
        slong previous = 2;
        for (slong j = terms(); j-- > 0 && exponents_[static_cast<size_t>(j)] >= 2;) {
            const slong e = exponents_[static_cast<size_t>(j)];
            if (previous > e) {
                mag_pow_ui(power.value(), norm.value(), static_cast<ulong>(previous - e));
                mag_mul(bound.value(), bound.value(), power.value());
            }
            mag_add(bound.value(), bound.value(), second_[static_cast<size_t>(j)].value());
            previous = e;
        }
        mag_pow_ui(power.value(), norm.value(), static_cast<ulong>(previous - 2));
        mag_mul(bound.value(), bound.value(), power.value());
        mag_mul(bound.value(), bound.value(), r.value());
        mag_mul(bound.value(), bound.value(), r.value());
        mag_mul_2exp_si(bound.value(), bound.value(), -1);
        arb_zero(out);
        arf_set_mag(arb_midref(out), bound.value());
    }

  private:
    slong precision_, power_precision_ = 0;
    std::vector<slong> exponents_;
    Plan plan_;
    RootVector coefficients_, slopes_, powers_, table_;
    // A table entry F: its midpoint, its disk radius E_F, and |F^| + E_F.
    std::vector<Magnitude> second_, table_error_, table_bound_;
    Magnitude rounding_;
    bool finite_ = true;
    uint64_t multiplications_ = 0;

    static std::vector<slong> support(const ComplexPolynomial &polynomial) {
        std::vector<slong> result;
        const auto coefficients = polynomial.coefficients();
        for (slong k = 0; k < coefficients.size(); ++k)
            if (!acb_is_zero(coefficients[k]))
                result.push_back(k);
        return result;
    }
    slong table_size() const {
        return plan_.windows * static_cast<slong>((ulong(1) << plan_.width) - 1);
    }
    // (a^, E_a, N_a) <- (a^, E_a, N_a) times the table's entry F, where N_a
    // bounds |a^|: with r the l1 radius of the rounded product's rectangle
    // (it dominates the Euclidean one),
    //     E <- E_a (|F^| + E_F) + N_a E_F + r,     N <- N_a (|F^| + E_F) + r.
    void multiply(acb_t midpoint, mag_t error, mag_t norm, slong entry) {
        ++multiplications_;
        const auto index = static_cast<size_t>(entry);
        mag_mul(error, error, table_bound_[index].value());
        mag_addmul(error, norm, table_error_[index].value());
        acb_mul(midpoint, midpoint, table_[entry], power_precision_);
        mag_add(rounding_.value(), arb_radref(acb_realref(midpoint)), arb_radref(acb_imagref(midpoint)));
        acb_get_mid(midpoint, midpoint);
        mag_add(error, error, rounding_.value());
        mag_mul(norm, norm, table_bound_[index].value());
        mag_add(norm, norm, rounding_.value());
    }
    // W_j = z^(e_j) along the support, each as its midpoint and disk radius.
    void build_powers(const acb_t point) {
        build_table(point);
        ComplexBall power;
        Magnitude error, norm;
        acb_one(power.value());
        mag_one(norm.value());
        const ulong mask = (ulong(1) << plan_.width) - 1;
        slong previous = 0;
        for (slong j = 0; j < terms(); ++j) {
            const slong e = exponents_[static_cast<size_t>(j)];
            slong window = 0;
            for (auto gap = static_cast<ulong>(e - previous); gap; gap >>= plan_.width, ++window)
                if (const ulong digit = gap & mask)
                    multiply(power.value(), error.value(), norm.value(),
                             window * static_cast<slong>(mask) + static_cast<slong>(digit) - 1);
            previous = e;
            // The rectangle of half-widths E holds the disk of radius E.
            acb_set(powers_[j], power.value());
            acb_add_error_mag(powers_[j], error.value());
        }
    }
    void store(slong entry, const acb_t midpoint, const mag_t error, const mag_t norm) {
        const auto index = static_cast<size_t>(entry);
        acb_set(table_[entry], midpoint);
        mag_set(table_error_[index].value(), error);
        mag_add(table_bound_[index].value(), norm, error);
    }
    void build_table(const acb_t point) {
        if (!plan_.windows)
            return;
        const auto entries = static_cast<slong>((ulong(1) << plan_.width) - 1);
        // The window's base z^(2^(w i)), then its multiples d = 2 .. 2^w - 1;
        // one more product hands the next window its base.
        ComplexBall power;
        Magnitude error, norm;
        acb_get_mid(power.value(), point);
        mag_add(error.value(), arb_radref(acb_realref(point)), arb_radref(acb_imagref(point)));
        acb_get_mag(norm.value(), power.value());
        for (slong window = 0; window < plan_.windows; ++window) {
            const slong first = window * entries;
            store(first, power.value(), error.value(), norm.value());
            for (slong digit = 1; digit < entries; ++digit) {
                multiply(power.value(), error.value(), norm.value(), first);
                store(first + digit, power.value(), error.value(), norm.value());
            }
            if (window + 1 < plan_.windows)
                multiply(power.value(), error.value(), norm.value(), first);
        }
    }
};
} // namespace certpolysolve::detail
