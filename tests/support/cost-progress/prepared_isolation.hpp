#pragma once
// Private prepared-factor research service. Inputs are Gaussian integers;
// square-freeness is an explicit fixture-manifest premise, not inferred here.
#include "../scaling_io.hpp"
#include "core/query.hpp"
#include "methods/subdivision/subdivision.hpp"
#include <algorithm>

namespace isolation_research {
namespace d = certpolysolve::detail;

class IntegerFactor {
    struct Coefficient { d::Integer real, imag; };
    std::vector<Coefficient> coefficients_;

  public:
    explicit IntegerFactor(const std::string &path) {
        const auto rows = scaling_io::read(path);
        if (rows.size() < 3) throw std::invalid_argument("prepared isolation requires degree at least two");
        coefficients_.resize(rows.size());
        for (size_t j = 0; j < rows.size(); ++j) {
            if (!d::Number(rows[j][2]).zero() || !d::Number(rows[j][3]).zero())
                throw std::invalid_argument("prepared factor must have exact coefficients");
            if (fmpz_set_str(coefficients_[j].real.value(), rows[j][0].c_str(), 10) ||
                fmpz_set_str(coefficients_[j].imag.value(), rows[j][1].c_str(), 10))
                throw std::invalid_argument("prepared factor requires decimal Gaussian-integer coefficients");
        }
        const auto &leading = coefficients_.back();
        if (fmpz_is_zero(leading.real.value()) && fmpz_is_zero(leading.imag.value()))
            throw std::invalid_argument("prepared factor must have fixed nonzero leading coefficient");
    }
    size_t degree() const { return coefficients_.size() - 1; }
    void enclose(acb_poly_t out, slong precision) const {
        acb_poly_zero(out);
        d::ComplexBall coefficient;
        for (size_t j = 0; j < coefficients_.size(); ++j) {
            arb_set_round_fmpz(acb_realref(coefficient.value()), coefficients_[j].real.value(), precision);
            arb_set_round_fmpz(acb_imagref(coefficient.value()), coefficients_[j].imag.value(), precision);
            acb_poly_set_coeff_acb(out, static_cast<slong>(j), coefficient.value());
        }
    }
};

inline cps_options isolation_options(uint64_t cap) {
    cps_options options;
    options.goal = CPS_ROOT_ISOLATION;
    options.bits = 0; // no requested-bits preparation/refinement on the explicit cut
    options.max_work_bits = cap;
    return options;
}

class Run {
    IntegerFactor input_;

  public:
    cps_options options;
    d::SearchStatistics statistics{};
    uint64_t turns = 0, enclosure_calls = 0, enclosure_bits_sum = 0, peak_enclosure_bits = 0;
    std::vector<d::Root> cover;
    bool goal_met = false;

  private:
    // The engine borrows statistics. Moving Run would invalidate that borrow.
    std::unique_ptr<d::Subdivision> engine_;
    void enclose(acb_poly_t p, slong precision) {
        ++enclosure_calls;
        enclosure_bits_sum = d::resources::add(enclosure_bits_sum, static_cast<uint64_t>(precision));
        peak_enclosure_bits = std::max(peak_enclosure_bits, static_cast<uint64_t>(precision));
        input_.enclose(p, precision);
    }

  public:
    Run(IntegerFactor input, uint64_t initial_bits, cps_options configured, d::SubdivisionPrecision policy)
        : input_(std::move(input)), options(std::move(configured)) {
        if (initial_bits < 64 || initial_bits > options.max_work_bits ||
            options.goal != CPS_ROOT_ISOLATION || options.bits != 0)
            throw std::invalid_argument("prepared isolation requires P0>=64, cap>=P0 and the explicit isolation cut");
        d::ComplexPolynomial p;
        enclose(p.value(), static_cast<slong>(initial_bits));
        d::Root global;
        if (!d::cauchy_cluster(global, p.value(), static_cast<slong>(initial_bits)))
            throw std::runtime_error("native Cauchy initialization was inconclusive");
        cover.push_back(std::move(global));
        engine_ = std::make_unique<d::Subdivision>(cover, 1, initial_bits, statistics, policy);
    }
    Run(const Run &) = delete;
    Run &operator=(const Run &) = delete;
    Run(Run &&) = delete;
    Run &operator=(Run &&) = delete;
    size_t degree() const { return input_.degree(); }
    d::Subdivision &engine() { return *engine_; }
    std::vector<d::SubdivisionPrecisionState> state() const {
        std::vector<d::SubdivisionPrecisionState> out;
        engine_->precision_state(out);
        return out;
    }
    // Continue the SAME prepared factor and search. This private driver has
    // a deliberate first-isolation boundary; the public solve is one-shot.
    void refine(uint64_t bits) {
        if (!goal_met || (options.goal != CPS_ROOT_ISOLATION && options.goal != CPS_ROOT_ACCURACY) ||
            !bits || bits > d::resources::precision_max)
            throw std::invalid_argument("refinement requires a completed isolation and positive representable bits");
        options.goal = CPS_ROOT_ACCURACY;
        options.bits = bits;
        options.absolute_tolerance = d::Number();
        engine_->reopen(options);
        goal_met = std::all_of(cover.begin(), cover.end(), [&](const d::Root &r) {
            return d::meets_goal(r, options, static_cast<slong>(peak_enclosure_bits));
        });
    }
    cps_stop_reason tick(double allowance = std::numeric_limits<double>::infinity()) {
        if (goal_met || !engine_->pending()) return CPS_STOP_NONE;
        ++turns;
        ++statistics.work_units;
        const auto used = engine_->cells();
        const auto reason = engine_->step(options, [&](acb_poly_t p, slong bits) { enclose(p, bits); },
                                          options.live_cells > used ? options.live_cells - used : 0,
                                          allowance, d::NewtonAccuracy::PreparedExact);
        std::vector<d::Root> next;
        if (engine_->complete_cover(next)) {
            const bool attained = std::all_of(next.begin(), next.end(), [&](const d::Root &r) {
                return d::meets_goal(r, options, static_cast<slong>(peak_enclosure_bits));
            });
            if (attained || d::prefer_cover(next, cover, static_cast<slong>(peak_enclosure_bits)))
                cover = std::move(next);
        }
        goal_met = std::all_of(cover.begin(), cover.end(), [&](const d::Root &r) {
            return d::meets_goal(r, options, static_cast<slong>(peak_enclosure_bits));
        });
        if (!engine_->pending() && !goal_met)
            throw std::logic_error("finished prepared search has not attained its goal");
        return reason;
    }
};
} // namespace isolation_research
