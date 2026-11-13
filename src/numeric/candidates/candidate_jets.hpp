#pragma once
#include <array>
#include "numeric/arithmetic/machine_arithmetic.hpp"
#include <bit>
#include <optional>
#include <span>
#include <vector>

namespace certpolysolve::detail::machine {
template <class R> struct Jet {
    Complex<R> value, derivative, second_derivative;
    R absolute_sum{}; // A condition estimate only, not an outward error bound.
};
template <class R, unsigned Order = 2> Jet<R> jet(std::span<const Complex<R>> coefficients, Complex<R> z) {
    Jet<R> q;
    if (coefficients.empty())
        return q;
    q.value = coefficients.back();
    q.absolute_sum = norm(q.value);
    const auto radius = norm(z);
    for (size_t j = coefficients.size() - 1; j-- > 0;) {
        if constexpr (Order >= 2)
            q.second_derivative = q.second_derivative * z + q.derivative * Complex<R>(R(2));
        q.derivative = q.derivative * z + q.value;
        q.value = q.value * z + coefficients[j];
        q.absolute_sum = q.absolute_sum * radius + norm(coefficients[j]);
    }
    return q;
}
// Immutable support and magnitudes belong to one candidate engine. A sparse Horner step propagates the entire
// derivative jet; evaluating at zero never divides by z. This is a midpoint evaluator only. evaluate writes nothing
// of the object (the sparse path's power tables are its own locals, at most 64 entries: a gap below 2^64), so calls
// may run at once on several threads (W9's split engine sweep, machine_candidates.cpp).
template <class R, unsigned Order = 2> class PreparedJets {
    using C = Complex<R>;
    static constexpr size_t power_slots = 64;
    std::span<const C> coefficients_;
    std::vector<size_t> indices_;
    std::vector<R> magnitudes_;
    size_t power_count_ = 0;
    static C power(const std::array<C, power_slots> &powers, size_t exponent) {
        C value(R(1));
        size_t bit = 0;
        while (exponent) {
            if (exponent & 1)
                value *= powers[bit];
            exponent >>= 1;
            ++bit;
        }
        return value;
    }
    static R radius_power(const std::array<R, power_slots> &radius_powers, size_t exponent) {
        R value(1);
        size_t bit = 0;
        while (exponent) {
            if (exponent & 1)
                value = value * radius_powers[bit];
            exponent >>= 1;
            ++bit;
        }
        return value;
    }

  public:
    explicit PreparedJets(std::span<const C> coefficients) : coefficients_(coefficients) {
        for (auto c : coefficients)
            magnitudes_.push_back(norm(c));
        for (size_t k = coefficients.size(); k-- > 0;)
            if (coefficients[k] != C{} || k == 0)
                indices_.push_back(k);
        size_t maximum_gap = 0;
        for (size_t j = 1; j < indices_.size(); ++j)
            maximum_gap = std::max(maximum_gap, indices_[j - 1] - indices_[j]);
        power_count_ = static_cast<size_t>(std::bit_width(maximum_gap));
    }
    size_t terms() const {
        return indices_.size() * 4 < coefficients_.size() ? indices_.size() : coefficients_.size();
    }
    size_t storage_bytes() const {
        // The power tables counted as the vectors they were (one entry a bit of the largest gap), so the
        // storage the resource checks see is unchanged.
        return indices_.capacity() * sizeof(size_t) + (magnitudes_.capacity() + power_count_) * sizeof(R) +
               power_count_ * sizeof(C);
    }
    Jet<R> evaluate(C z) const {
        Jet<R> q;
        if (coefficients_.empty())
            return q;
        q.value = coefficients_.back();
        q.absolute_sum = magnitudes_.back();
        const R radius = norm(z);
        if (terms() == coefficients_.size()) {
            for (size_t j = coefficients_.size() - 1; j-- > 0;) {
                if constexpr (Order >= 2)
                    q.second_derivative = q.second_derivative * z + q.derivative * C(R(2));
                q.derivative = q.derivative * z + q.value;
                q.value = q.value * z + coefficients_[j];
                q.absolute_sum = q.absolute_sum * radius + magnitudes_[j];
            }
            return q;
        }
        std::array<C, power_slots> powers;
        std::array<R, power_slots> radius_powers;
        if (power_count_) {
            powers[0] = z;
            radius_powers[0] = radius;
        }
        for (size_t j = 1; j < power_count_; ++j) {
            powers[j] = powers[j - 1] * powers[j - 1];
            radius_powers[j] = radius_powers[j - 1] * radius_powers[j - 1];
        }
        for (size_t j = 1; j < indices_.size(); ++j) {
            const size_t gap = indices_[j - 1] - indices_[j];
            const C factor = power(powers, gap),
                    derivative = power(powers, gap - 1) * C(R(static_cast<double>(gap)));
            if constexpr (Order >= 2) {
                const C second = gap < 2 ? C{} : power(powers, gap - 2) *
                    C(R(static_cast<double>(gap)) * R(static_cast<double>(gap - 1)));
                q.second_derivative = q.second_derivative * factor +
                    q.derivative * derivative * C(R(2)) + q.value * second;
            }
            q.derivative = q.derivative * factor + q.value * derivative;
            q.value = q.value * factor + coefficients_[indices_[j]];
            q.absolute_sum = q.absolute_sum * radius_power(radius_powers, gap) + magnitudes_[indices_[j]];
        }
        return q;
    }
};
// A factored Newton ratio with a residual condition estimate for the represented
// polynomial. These fields are not promised to be the original p and p'.
template <class R> struct NewtonEvaluation {
    using C = Complex<R>;
    C value, denominator, coordinate_factor;
    R absolute_sum{};
    size_t terms = 0;
    std::optional<C> ratio() const {
        if (denominator == C{})
            return std::nullopt;
        const C result = coordinate_factor * (value / denominator);
        return finite(result) ? std::optional<C>(result) : std::nullopt;
    }
};

// The factor owner keeps the forward coefficient storage alive. Reverse
// coefficients and both sets of support/power scratch live in this evaluator.
template <class R> class PreparedNewton {
    using C = Complex<R>;
    PreparedJets<R, 1> forward_;
    std::vector<C> reversed_coefficients_;
    PreparedJets<R, 1> backward_;
    size_t degree_;
    static std::vector<C> reverse(std::span<const C> coefficients) {
        std::vector<C> result(coefficients.rbegin(), coefficients.rend());
        // p can contain an exact zero factor. q(w)=w^n p(1/w) then has
        // degree below n; remove only its leading zero padding.
        while (result.size() > 1 && result.back() == C{})
            result.pop_back();
        return result;
    }

  public:
    explicit PreparedNewton(std::span<const C> coefficients)
        : forward_(coefficients), reversed_coefficients_(reverse(coefficients)),
          backward_(reversed_coefficients_), degree_(coefficients.empty() ? 0 : coefficients.size() - 1) {}
    // backward_ borrows our reverse buffer. Keep this workspace at a stable
    // address; callers transfer its unique owner instead of copying views.
    PreparedNewton(const PreparedNewton &) = delete;
    PreparedNewton &operator=(const PreparedNewton &) = delete;
    PreparedNewton(PreparedNewton &&) = delete;
    PreparedNewton &operator=(PreparedNewton &&) = delete;
    NewtonEvaluation<R> evaluate(C z) const {
        if (R(1) < norm(z)) {
            const C w = C(R(1)) / z;
            const auto q = backward_.evaluate(w);
            // p/p' = z*q(w) / (n*q(w)-w*q'(w)), w=1/z.
            return {q.value, C(R(static_cast<double>(degree_))) * q.value - w * q.derivative,
                    z, q.absolute_sum, backward_.terms()};
        }
        const auto p = forward_.evaluate(z);
        return {p.value, p.derivative, C(R(1)), p.absolute_sum, forward_.terms()};
    }
    size_t storage_bytes() const {
        return forward_.storage_bytes() + backward_.storage_bytes() +
               reversed_coefficients_.capacity() * sizeof(C);
    }
};

template <class R> Complex<R> laguerre_step(const Jet<R> &q, size_t degree) {
    if (q.value == Complex<R>{})
        return {};
    const auto g = q.derivative / q.value;
    const auto h = g * g - q.second_derivative / q.value;
    const auto radical = sqrt(Complex<R>(R(static_cast<double>(degree - 1))) *
                              (Complex<R>(R(static_cast<double>(degree))) * h - g * g));
    const auto plus = g + radical, minus = g - radical;
    return Complex<R>(R(static_cast<double>(degree))) / (norm(plus) < norm(minus) ? minus : plus);
}
// Skowron--Gould's local switch uses F=p*p''/(p')^2. Away from the
// Newton basin use Laguerre; small F permits the second-order correction.
template <class R> Complex<R> laguerre_newton_step(const Jet<R> &q, size_t degree) {
    if (q.value == Complex<R>{})
        return {};
    if (q.derivative != Complex<R>{}) {
        const auto newton = q.value / q.derivative;
        const auto f = newton * (q.second_derivative / q.derivative);
        if (norm(f) < R(0.05))
            return newton;
        if (norm(f) < R(0.5))
            return newton * (Complex<R>(R(1)) + f * Complex<R>(R(0.5)));
    }
    return laguerre_step(q, degree);
}
} // namespace certpolysolve::detail::machine
