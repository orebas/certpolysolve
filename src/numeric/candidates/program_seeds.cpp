// Starting points for the program route at the roots' radii (program_seeds, declared in program_seeds.hpp with
// circle_count and circle_sample): a program has no coefficients for a Newton polygon, so the numbers of roots inside
// circles about the origin, read from the program's own evaluations, place the seeds.
//
// circle_count takes one count: the Cauchy sum (1/q) sum z p'(z) / p(z) over q equally spaced samples (circle_sample)
// is about the number of roots inside the circle, with an error from roots near it that its imaginary part and its ball
// radius measure (Pan 2022, Theorem 8). With a machine evaluator the samples go through binary64 first, in one batched
// call (circle_count_machine), then, for a circle whose values leave binary64's range, through the range tier
// (machine::BlockDouble); a sample at the evaluator's floor, or a sum whose error estimate exceeds
// program_count_tolerance, sends the circle to the balls at the current precision, which say why a count fails: p and
// p' without correct bits (cancellation) or a root near a sample (near_root).
//
// program_seeds, within a budget of evaluations: the coarse ladder, the radius halving from the root bound (whose count
// is the degree) until no root is inside; then the annuli between rungs, the fullest bisected in log radius until each
// holds at most the capacity or reaches the resolution. Each count is clamped between its neighbours' (the counts are
// monotone in the radius) and has its remedies (take): no correct bits, the precision doubles; a root at a sample, the
// samples turn; an ambiguous count, the samples double; the precision and the samples are never lowered. Then each
// annulus's roots go on its middle circle at the placement angles (placement_angle, with a circle index of its own),
// and the seeds no count placed on the root bound's circle. The counts are proposals: a wrong one costs sweeps, and
// nothing here certifies.
//
// State: none. The report (ProgramSeedReport) records the circles, the remedies, the precision reached and the annuli.
// Entry points: program_seeds, from solver_program.cpp (program_advance: the preliminary phase's start, with the
// binary64 evaluator for the counts) and program_recurrence_seeds.cpp (a sub-program no structural rule seeds);
// circle_count and circle_sample, from tests/test_program.cpp.
// Not here: the evaluators (program_machine.cpp; ProgramSource for the balls), the coefficient route's start from the
// Newton polygon (initialize in candidates.cpp), the structural seeds (program_recurrence_seeds.cpp).
// Design notes: docs/slp-seed-placement-note-20261016.md (section 3, method A; section 6, the decision; section 8, the
// counts in binary64).
#include "numeric/candidates/program_seeds.hpp"
#include "numeric/evaluate/monic_source.hpp"
#include "numeric/candidates/candidates.hpp" // placement_angle
#include "core/heuristics.hpp"
#include "numeric/candidates/machine_adapter.hpp"
#include "core/resource_limits.hpp"
#include <algorithm>
#include <cmath>
#include <numbers>
#include <memory>
#include <optional>
#include <queue>

namespace certpolysolve::detail {
namespace {
// placement_angle's offset: the samples and the seeds share no symmetry with the source.
constexpr double sample_offset = 0.2718281828459045;
// The radius 2^log2_radius as an exact floating-point number: 2^fraction in binary64 times an exact power of
// two, so a bound beyond binary64's range still has its circle (initialize in candidates.cpp splits its radii
// the same way).
void exact_radius(arf_t out, double log2_radius) {
    const double whole = std::floor(log2_radius);
    arf_set_d(out, std::exp2(log2_radius - whole));
    arf_mul_2exp_si(out, out, static_cast<slong>(whole));
}
// The point at `angle` on the circle of radius 2^log2_radius, as an exact midpoint.
void circle_point(acb_t out, double log2_radius, double angle) {
    Float radius, direction;
    exact_radius(radius.value(), log2_radius);
    acb_zero(out);
    arf_set_d(direction.value(), std::cos(angle));
    arf_mul(arb_midref(acb_realref(out)), radius.value(), direction.value(), 53, ARF_RND_NEAR);
    arf_set_d(direction.value(), std::sin(angle));
    arf_mul(arb_midref(acb_imagref(out)), radius.value(), direction.value(), 53, ARF_RND_NEAR);
}
} // namespace

void circle_sample(acb_t out, double log2_radius, unsigned samples, double rotation, unsigned g) {
    circle_point(out, log2_radius,
                 2 * std::numbers::pi * (static_cast<double>(g) + rotation) / static_cast<double>(samples) + sample_offset);
}

namespace {
// The count from a machine evaluator (binary64, or the range tier machine::BlockDouble), or nullopt when
// this circle is for another tier: `out_of_range` set when a sample or its evaluation leaves the arithmetic's
// range (the range tier may read it), clear when a sample is at the floor (a root at a sample included: its
// value is unreadable in this arithmetic, and the balls then say which) or the error estimate exceeds the
// tolerance (the balls). The terms z p'(z) / p(z) are read in binary64 whatever the arithmetic: they are of the
// order of the degree on any circle a count is taken on.
template <class R>
std::optional<CircleCount> circle_count_machine(ProgramMachine<R> &machine, double log2_radius, unsigned samples,
                                                double rotation, certpolysolve::detail::CandidateStatistics &statistics, bool &out_of_range) {
    using C = machine::Complex<R>;
    using D = machine::Complex<double>;
    out_of_range = false;
    std::vector<C> points(samples);
    ComplexBall z;
    for (unsigned g = 0; g < samples; ++g) {
        circle_sample(z.value(), log2_radius, samples, rotation, g);
        if (!from_ball(points[g], z.value())) {
            out_of_range = true;
            return std::nullopt;
        }
    }
    std::vector<typename ProgramMachine<R>::Evaluation> evaluations(samples);
    machine.evaluate(points, evaluations, heuristics::program_machine_chunk);
    statistics.evaluations += samples;
    statistics.evaluation_terms += static_cast<uint64_t>(samples) * machine.terms();
    CircleCount count;
    count.evaluations = samples;
    count.binary64 = true;
    const R floor_factor = R(heuristics::program_machine_floor_factor);
    const double u = 0x1p-53;
    D sum;
    double error = 0;
    for (unsigned g = 0; g < samples; ++g) {
        const auto &e = evaluations[g];
        if (!e.finite) {
            out_of_range = true;
            return std::nullopt;
        }
        const R p = machine::norm(e.value), d = machine::norm(e.derivative);
        if (e.value == C{} || p < floor_factor * e.noise || e.derivative == C{} || d < floor_factor * e.noise_derivative)
            return std::nullopt;
        const C term = points[g] * e.derivative / e.value;
        const D t{static_cast<double>(term.re), static_cast<double>(term.im)};
        const double magnitude = machine::norm(t);
        if (!std::isfinite(magnitude))
            return std::nullopt;
        sum += t;
        error += magnitude * (static_cast<double>(e.noise / p) + static_cast<double>(e.noise_derivative / d)) +
                 4 * u * magnitude;
    }
    const double q = static_cast<double>(samples);
    count.real = sum.re / q;
    count.imag = sum.im / q;
    count.radius = error / q;
    count.finite = std::isfinite(count.real) && std::isfinite(count.imag) && std::isfinite(count.radius);
    if (!count.finite || count.radius > heuristics::program_count_tolerance)
        return std::nullopt;
    return count;
}
} // namespace

CircleCount circle_count(MonicSource &source, double log2_radius, unsigned samples, double rotation, slong precision,
                         certpolysolve::detail::CandidateStatistics &statistics, ProgramMachine<double> *machine,
                         ProgramMachine<machine::BlockDouble> *block) {
    uint64_t machine_evaluations = 0; // tried and not read: counted toward the budget by the caller
    if (machine) {
        bool out_of_range = false;
        if (auto count = circle_count_machine(*machine, log2_radius, samples, rotation, statistics, out_of_range))
            return *count;
        machine_evaluations += samples;
        if (out_of_range && block) {
            if (auto count = circle_count_machine(*block, log2_radius, samples, rotation, statistics, out_of_range)) {
                count->evaluations += machine_evaluations;
                return *count;
            }
            machine_evaluations += samples;
        }
    }
    CircleCount count;
    count.evaluations = machine_evaluations;
    ComplexBall z, value, derivative, term, sum;
    for (unsigned g = 0; g < samples; ++g) {
        circle_sample(z.value(), log2_radius, samples, rotation, g);
        ++statistics.evaluations;
        ++count.evaluations;
        if (!source.evaluate(value.value(), derivative.value(), z.value(), precision)) {
            count.cancellation = true;
            return count;
        }
        if (acb_contains_zero(value.value())) {
            // p without correct bits: with p' as bad, the precision does not evaluate the program here (a
            // recurrence cancelling at its large roots); with p' fine, a root near the sample.
            if (acb_contains_zero(derivative.value()))
                count.cancellation = true;
            else
                count.near_root = true;
            return count;
        }
        acb_div(term.value(), derivative.value(), value.value(), precision);
        acb_mul(term.value(), term.value(), z.value(), precision);
        acb_add(sum.value(), sum.value(), term.value(), precision);
    }
    acb_div_ui(sum.value(), sum.value(), samples, precision);
    count.real = arf_get_d(arb_midref(acb_realref(sum.value())), ARF_RND_NEAR);
    count.imag = arf_get_d(arb_midref(acb_imagref(sum.value())), ARF_RND_NEAR);
    count.radius = mag_get_d(arb_radref(acb_realref(sum.value()))) + mag_get_d(arb_radref(acb_imagref(sum.value())));
    count.finite = std::isfinite(count.real) && std::isfinite(count.imag) && std::isfinite(count.radius);
    return count;
}

bool program_seeds(acb_ptr out, MonicSource &source, slong degree, const mag_t root_bound, slong precision,
                   const cps_options &options, certpolysolve::detail::CandidateStatistics &statistics, const EnginePolicy *policy,
                   ProgramSeedReport *report_out, ProgramMachine<double> *machine) {
    ProgramSeedReport report;
    const auto n = static_cast<uint64_t>(std::max<slong>(degree, 0));
    const double top = mag_get_d_log2_approx(root_bound);
    const uint64_t budget = std::max<uint64_t>(heuristics::program_seed_budget_floor,
                                               std::min<uint64_t>(n, uint64_t{1} << 40) * heuristics::program_seed_budget_sweeps);
    const slong max_bits = std::max<slong>(
        precision, static_cast<slong>(std::min<uint64_t>(options.max_work_bits, resources::precision_max)));
    slong bits = std::min(precision, max_bits);
    unsigned samples = heuristics::program_count_samples;
    // The range tier for circles whose samples leave binary64's range (the compositions' and the recurrences'
    // outer circles), from the binary64 evaluator.
    std::unique_ptr<ProgramMachine<machine::BlockDouble>> block;
    if (machine && machine->valid())
        block = std::make_unique<ProgramMachine<machine::BlockDouble>>(*machine);
    const unsigned max_samples = heuristics::program_count_samples_max;

    // The count on one circle, clamped into [lower, upper] (the counts are monotone in the radius and the
    // bound's count is proved). Each failure has its remedy: no correct bits, the precision doubles; a root at
    // a sample, the samples turn; an ambiguous value (not within the tolerance of an integer with an imaginary
    // part as small), the samples double. The precision and the samples are never lowered again: what one
    // circle needed, the next is likely to need. nullopt: cancelled, the budget spent, or a circle whose count
    // could not be taken (left out of the ladder).
    const auto take = [&](double log2_radius, uint64_t lower, uint64_t upper) -> std::optional<uint64_t> {
        double rotation = 0;
        unsigned rotations = 0;
        for (;;) {
            if (options.cancelled()) {
                report.cancelled = true;
                return std::nullopt;
            }
            if (report.evaluations + samples > budget) {
                report.budget_spent = true;
                return std::nullopt;
            }
            const auto count = circle_count(source, log2_radius, samples, rotation, bits, statistics, machine, block.get());
            report.evaluations += count.evaluations;
            const bool unreadable = count.cancellation || (count.finite && count.radius > heuristics::program_count_tolerance);
            if (unreadable) {
                if (bits >= max_bits)
                    return std::nullopt;
                bits = std::min(max_bits, resources::precision_product(bits, 2));
                ++report.raises;
                continue;
            }
            if (count.near_root) {
                if (rotations >= 2)
                    return std::nullopt;
                rotation = rotations ? 0.25 : 0.5;
                ++rotations;
                ++report.rotations;
                continue;
            }
            if (!count.finite)
                return std::nullopt;
            const double rounded = std::round(count.real);
            const bool sharp = std::abs(count.real - rounded) <= heuristics::program_count_tolerance &&
                               std::abs(count.imag) <= heuristics::program_count_tolerance;
            if (!sharp && samples < max_samples) {
                samples = std::min(max_samples, samples * 2);
                ++report.sample_doublings;
                continue;
            }
            ++report.circles;
            if (!sharp)
                ++report.unsharp;
            if (count.binary64)
                ++report.binary64_circles;
            const double clamped = std::clamp(rounded, static_cast<double>(lower), static_cast<double>(upper));
            return static_cast<uint64_t>(clamped);
        }
    };

    // The coarse ladder: the radius halves from the bound, whose count is the degree, until no root is inside.
    struct Rung {
        double log2_radius;
        uint64_t inside;
    };
    std::vector<Rung> rungs{{top, n}};
    bool stopped = false; // the budget or a cancellation
    for (unsigned i = 1; rungs.back().inside > 0 && i <= heuristics::program_ladder_max_rungs && !stopped; ++i) {
        const double log2_radius = top - static_cast<double>(i) * heuristics::program_ladder_step_log2;
        const auto inside = take(log2_radius, 0, rungs.back().inside);
        if (!inside) {
            stopped = report.budget_spent || report.cancelled;
            continue; // a circle without a count: the next rung's annulus is wider
        }
        rungs.push_back({log2_radius, *inside});
    }
    if (report.cancelled) {
        if (report_out)
            *report_out = std::move(report);
        return false;
    }
    // The annuli between consecutive rungs, the fullest bisected first (in log radius) until each holds at most
    // the capacity or is at the resolution.
    struct Annulus {
        double outer, inner;
        uint64_t outer_inside, inner_inside;
        uint64_t count() const { return outer_inside - inner_inside; }
    };
    std::vector<Annulus> annuli;
    std::priority_queue<std::pair<uint64_t, size_t>> fullest;
    for (size_t i = 1; i < rungs.size(); ++i) {
        annuli.push_back({rungs[i - 1].log2_radius, rungs[i].log2_radius, rungs[i - 1].inside, rungs[i].inside});
        fullest.emplace(annuli.back().count(), annuli.size() - 1);
    }
    const uint64_t capacity =
        std::max<uint64_t>(heuristics::program_annulus_capacity_min, n / heuristics::program_annulus_capacity_fraction);
    while (!fullest.empty() && !stopped) {
        const auto [count, index] = fullest.top();
        fullest.pop();
        const Annulus annulus = annuli[index];
        if (count <= capacity || annulus.outer - annulus.inner <= heuristics::program_ladder_resolution_log2 * 1.5)
            continue;
        const double middle = (annulus.outer + annulus.inner) / 2;
        const auto inside = take(middle, annulus.inner_inside, annulus.outer_inside);
        if (!inside) {
            stopped = report.budget_spent || report.cancelled;
            continue;
        }
        annuli[index] = {annulus.outer, middle, annulus.outer_inside, *inside};
        annuli.push_back({middle, annulus.inner, *inside, annulus.inner_inside});
        fullest.emplace(annuli[index].count(), index);
        fullest.emplace(annuli.back().count(), annuli.size() - 1);
    }
    if (report.cancelled) {
        if (report_out)
            *report_out = std::move(report);
        return false;
    }
    // What lies inside the lowest rung (the ladder stopped short of an empty circle) is an annulus a step below it.
    // Without a single count the ladder says nothing, and every seed goes on the bound's circle below.
    if (rungs.size() > 1 && rungs.back().inside > 0)
        annuli.push_back({rungs.back().log2_radius, rungs.back().log2_radius - heuristics::program_ladder_step_log2,
                          rungs.back().inside, 0});
    std::sort(annuli.begin(), annuli.end(), [](const Annulus &a, const Annulus &b) { return a.outer > b.outer; });

    // The placement: each annulus's roots on its middle circle, at the placement angles with a circle index of
    // its own (a shared index would leave two circles' placements symmetric about the same axis).
    uint64_t position = 0;
    slong circle = 1;
    for (const auto &annulus : annuli) {
        const uint64_t count = annulus.count();
        if (!count)
            continue;
        const double middle = (annulus.outer + annulus.inner) / 2;
        for (uint64_t j = 0; j < count; ++j, ++position)
            circle_point(out + static_cast<slong>(position), middle,
                         placement_angle(static_cast<slong>(j), static_cast<slong>(count), circle));
        report.annuli.push_back({annulus.outer, annulus.inner, count});
        ++circle;
    }
    report.placed = position;
    // The seeds no count placed, on the circle of the bound as before.
    const uint64_t remaining = n - position;
    for (uint64_t j = 0; j < remaining; ++j, ++position)
        circle_point(out + static_cast<slong>(position), top,
                     placement_angle(static_cast<slong>(j), static_cast<slong>(remaining), 0));
    report.precision = bits;
    if (policy && policy->trace) {
        // The seeds placed, the program evaluations the counts took, the precision they ended at. The counts
        // themselves are in the diagnostics (ProgramPreliminary::describe): nothing here is certified.
        TraceEvent event{"program_seeds", "placement", static_cast<uint64_t>(bits), report.placed, &statistics, nullptr};
        event.model_evaluations = report.evaluations;
        event.reason = report.budget_spent ? "budget-spent" : "";
        policy->trace(policy->trace_data, event);
    }
    if (report_out)
        *report_out = std::move(report);
    return true;
}

} // namespace certpolysolve::detail
