#pragma once
// A real polynomial with few terms, p(x) = sum_j c_j x^(e_j) with exact real coefficients and ascending exponents,
// evaluated term by term: enclosures of p and p' at an exact real point at a precision, with the powers built
// along the gaps of the support by repeated squaring (about 2 log2(gap) multiplications a term), never through
// the dense coefficient list. For B6 (docs/b6-sparse-real-roots-note-20261015.md): the evaluator of
// real_sign and real_refine (real_sign_refinement.hpp) over the input and its fractional derivatives, where the
// dense evaluation would cost the degree an evaluation (10^6 at fewnomial-1000000-30).
//
// p'(x) = sum_j e_j c_j x^(e_j - 1) is formed from the powers x^(e_j - 1) directly, and p from x times them (plus
// the constant term), so x = 0 needs no division. The coefficients' enclosures are cached per precision, as the
// axis charts cache their dense source (piecewise_roots.cpp, axis_source).
#include "core/number.hpp"
#include "core/raii.hpp"
#include <flint/arb.h>
#include <bit>
#include <cstdint>
#include <stdexcept>
#include <utility>
#include <vector>

namespace certpolysolve::detail {

struct SparseRealTerm {
    slong exponent = 0;
    Number coefficient; // exact
};

class SparseRealEvaluator {
    std::vector<SparseRealTerm> terms_;
    mutable std::vector<RealBall> enclosed_; // the coefficients at enclosed_bits_
    mutable slong enclosed_bits_ = 0;
    mutable uint64_t multiplications_ = 0;

    const std::vector<RealBall> &enclosed(slong bits) const {
        if (enclosed_bits_ != bits) {
            enclosed_.resize(terms_.size());
            for (size_t j = 0; j < terms_.size(); ++j) terms_[j].coefficient.enclose(enclosed_[j].value(), bits);
            enclosed_bits_ = bits;
        }
        return enclosed_;
    }

  public:
    // Terms with ascending, distinct, nonnegative exponents (a Number is exact by construction).
    explicit SparseRealEvaluator(std::vector<SparseRealTerm> terms) : terms_(std::move(terms)) {
        slong previous = -1;
        for (const auto &term : terms_) {
            if (term.exponent <= previous) throw std::invalid_argument("sparse terms must have ascending exponents");
            previous = term.exponent;
        }
    }
    // The dense polynomial's nonzero coefficients (for tests and for the input's own polynomial).
    static SparseRealEvaluator from_dense(const std::vector<Number> &coefficients) {
        std::vector<SparseRealTerm> terms;
        for (size_t j = 0; j < coefficients.size(); ++j)
            if (!coefficients[j].zero()) terms.push_back({static_cast<slong>(j), coefficients[j]});
        return SparseRealEvaluator(std::move(terms));
    }
    size_t terms() const { return terms_.size(); }
    slong degree() const { return terms_.empty() ? -1 : terms_.back().exponent; }
    // The real multiplications made by every evaluation so far (the work meter's count).
    uint64_t multiplications() const { return multiplications_; }

    // Enclosures of p(x) and p'(x) at `bits`: the interface of real_sign and real_refine.
    void evaluate2(arb_t value, arb_t slope, const arb_t x, slong bits) const {
        const auto &c = enclosed(bits);
        arb_zero(value);
        arb_zero(slope);
        if (terms_.empty()) return;
        // sum_{e_j >= 1} c_j x^(e_j - 1) and sum_{e_j >= 1} e_j c_j x^(e_j - 1), the powers along the gaps.
        RealBall power, gap_power, term, tail, tail_slope;
        bool have_power = false;
        slong at = 0; // the exponent `power` holds, when have_power
        for (size_t j = 0; j < terms_.size(); ++j) {
            const slong e = terms_[j].exponent;
            if (e == 0) {
                arb_add(value, value, c[j].value(), bits);
                continue;
            }
            const slong want = e - 1;
            if (!have_power) {
                arb_pow_ui(power.value(), x, static_cast<ulong>(want), bits);
                multiplications_ += want ? static_cast<uint64_t>(2 * std::bit_width(static_cast<ulong>(want))) : 0;
                have_power = true;
            } else if (want > at) {
                arb_pow_ui(gap_power.value(), x, static_cast<ulong>(want - at), bits);
                arb_mul(power.value(), power.value(), gap_power.value(), bits);
                multiplications_ += static_cast<uint64_t>(2 * std::bit_width(static_cast<ulong>(want - at))) + 1;
            }
            at = want;
            arb_mul(term.value(), c[j].value(), power.value(), bits);
            arb_add(tail.value(), tail.value(), term.value(), bits);
            arb_mul_si(term.value(), term.value(), e, bits);
            arb_add(tail_slope.value(), tail_slope.value(), term.value(), bits);
            multiplications_ += 2;
        }
        arb_addmul(value, tail.value(), x, bits);
        arb_set(slope, tail_slope.value());
        ++multiplications_;
    }
};

} // namespace certpolysolve::detail
