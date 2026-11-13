#include "solve/policies.hpp"
// An input whose coefficients are balls bounds what any precision can
// certify. A residue whose error is the input's is accepted as the best the
// input allows, the secular iteration stops for the input instead of climbing
// to the precision budget, and the factor goes to the complete method at the
// precision that stop was taken at (regenerate_residues, the growth in
// SecularRefinement::step, structured_step's handoff).
//
// The explicit methods refuse inexact input through the public session, so
// these checks drive the maintained secular recipe through the internal
// driver, as test_method_handoff does.
#include "numeric/arithmetic/polynomial_arithmetic.hpp"
#include "exact_winding_api.hpp"
#include "internal.hpp"
#include <iostream>
#include <memory>
#include <string_view>

namespace d = certpolysolve::detail;
using namespace reference;

namespace {
// Twelve conjugate pairs with three-decimal coordinates: moduli from 0.27 to
// 2.6, the closest two roots 0.19 apart.
std::vector<Complex> planted() {
    const std::pair<const char *, const char *> upper[] = {
        {"-1847/1000", "177/1000"}, {"-1763/1000", "233/125"},   {"-857/500", "493/1000"},
        {"-213/125", "841/500"},    {"-1629/1000", "1129/1000"}, {"-1383/1000", "809/1000"},
        {"-337/500", "971/500"},    {"-503/1000", "597/500"},    {"-28/125", "857/1000"},
        {"39/500", "11/25"},        {"97/500", "193/1000"},      {"333/500", "99/1000"}};
    std::vector<Complex> roots;
    for (const auto &[real, imag] : upper) {
        roots.emplace_back(Rational(real), Rational(imag));
        roots.emplace_back(Rational(real), -Rational(imag));
    }
    return roots;
}
// The polynomial with every coefficient a ball of relative radius 2^-exponent.
certpolysolve::Polynomial with_radii(const ComplexPolynomial &p, int exponent) {
    certpolysolve::Polynomial result;
    for (size_t j = 0; j < p.size(); ++j) {
        require(p[j].imag == 0, "the fixture has real coefficients");
        const Rational radius = absolute(p[j].real) * power2(-exponent);
        result.set(j, certpolysolve::ball(certpolysolve::Number(p[j].real.str()), certpolysolve::Number(radius.str())));
    }
    return result;
}
struct Handoffs {
    unsigned input = 0, other = 0;
    uint64_t bits = 0; // the factor's precision at the input-floor handoff
};
struct Solve {
    cps_result result;
    Handoffs handoffs;
    uint64_t source_peak_bits = 0;
};
// `repair`: the certificate attempt's double-double repair (EnginePolicy::double_double_repair, on in every recipe).
// Off for the checks of the secular method's input stop and handoff, which the repair would otherwise pre-empt where
// it lets the binary64 phase prove the roots.
Solve solve(const ComplexPolynomial &polynomial, int radius_exponent, uint64_t cap, bool repair = false) {
    const auto input = with_radii(polynomial, radius_exponent);
    auto policy = d::maintained_method_policy(CPS_METHOD_SECULAR);
    policy.double_double_repair = repair;
    policy.study_statistics = std::make_shared<d::StudyStatistics>();
    Solve out;
    policy.trace_data = &out.handoffs;
    policy.trace = [](void *p, const d::TraceEvent &event) {
        if (std::string_view(event.stage) != "secular_subdivision_handoff")
            return;
        auto &seen = *static_cast<Handoffs *>(p);
        if (std::string_view(event.reason) == "input-floor") {
            ++seen.input;
            seen.bits = event.precision;
        } else
            ++seen.other;
    };
    cps_options options;
    options.goal = CPS_ROOT_APPROXIMATION;
    options.bits = 53;
    options.max_work_bits = cap;
    options.work_units = 100000;
    options.live_cells = 65536;
    d::SolveState solver(certpolysolve::detail::dense(*input.get()), policy);
    out.result = solver.advance(&options);
    out.source_peak_bits = policy.study_statistics->source_peak_bits;
    return out;
}
// Every exported region against the roots of the midpoint polynomial, which
// is one of the polynomials the balls admit: a region's count is the number
// of planted roots in its closed rectangle (exact rational comparisons).
void check_cover(const std::vector<Complex> &roots, cps_result result, const char *message) {
    const auto published = d::make_result(CPS_OK, d::published_result(std::move(result)));
    uint64_t total = 0;
    for (const auto &region : published.roots()) {
        const Rational real_lower = endpoint(region.box.real.lo), real_upper = endpoint(region.box.real.hi);
        const Rational imag_lower = endpoint(region.box.imag.lo), imag_upper = endpoint(region.box.imag.hi);
        uint64_t inside = 0;
        for (const auto &root : roots)
            inside += real_lower <= root.real && root.real <= real_upper && imag_lower <= root.imag &&
                      root.imag <= imag_upper;
        require(inside == region.count, message);
        total += region.count;
    }
    require(total == roots.size(), "a cover of an inexact input lost roots");
}
} // namespace

int main() {
    certpolysolve::ThreadScope lifetime;
    try {
        const auto roots = planted();
        const auto polynomial = from_roots(roots);
        // The precision budget of these checks. A ladder that climbs runs to
        // it; the stops asserted below are taken two rungs or more beneath.
        constexpr uint64_t cap = 4096;
        constexpr uint64_t rungs_beneath = 4; // a factor of four: two doublings
        // Relative radii of 2^-20 forbid a 53-bit approximation. The promoted
        // library regenerated every residue up the whole ladder (an error
        // bound that is the input's never becomes adequate) and stopped for
        // precision at the budget.
        {
            auto run = solve(polynomial, 20, cap);
            require(run.result.complete && !run.result.goal_met,
                    "a goal the coefficient radii forbid was reported met, or the cover is incomplete");
            require(run.handoffs.other == 0,
                    "secular stopped for precision, not for the input, on residues the input limits");
            require(run.handoffs.input == 1, "secular's input stop did not hand the factor to the complete method");
            require(run.source_peak_bits * rungs_beneath <= cap,
                    "regeneration climbed the precision ladder on residues the input limits");
            require(run.handoffs.bits * rungs_beneath <= cap,
                    "the iteration's precision climbed before the input stop");
            check_cover(roots, std::move(run.result),
                        "a region of the input-limited cover disagrees with the midpoint polynomial's roots");
        }
        // Radii of 2^-80 allow the goal, but not through secular's disks,
        // which are a factor of the degree wider than the input's floor: the
        // input stop hands the factor to the complete method, which meets it.
        {
            auto run = solve(polynomial, 80, cap);
            require(run.result.complete && run.result.goal_met && run.result.termination == CPS_SUCCESS &&
                        d::disjoint(run.result.roots),
                    "a goal within the input's radii was not met after secular's input stop");
            require(run.handoffs.input == 1 && run.handoffs.other == 0,
                    "the reachable goal did not go through the input-floor handoff");
            require(run.source_peak_bits * rungs_beneath <= cap,
                    "regeneration climbed the precision ladder on residues the input limits");
            check_cover(roots, std::move(run.result),
                        "a region of the reachable cover disagrees with the midpoint polynomial's roots");
            // With the repair (the recipes' default) the binary64 phase proves these roots itself: the goal is met
            // with no handoff at all.
            auto repaired = solve(polynomial, 80, cap, true);
            require(repaired.result.complete && repaired.result.goal_met && repaired.result.termination == CPS_SUCCESS &&
                        d::disjoint(repaired.result.roots),
                    "a goal within the input's radii was not met with the double-double repair");
            require(repaired.handoffs.input == 0 && repaired.handoffs.other == 0,
                    "the repair's solve still went through a handoff");
            check_cover(roots, std::move(repaired.result),
                        "a region of the repaired cover disagrees with the midpoint polynomial's roots");
        }
        // Exact input is untouched: the same polynomial with zero radii
        // meets the goal with no input stop.
        {
            const auto input = reference::input(polynomial);
            Handoffs handoffs;
            auto policy = d::maintained_method_policy(CPS_METHOD_SECULAR);
            policy.trace_data = &handoffs;
            policy.trace = [](void *p, const d::TraceEvent &event) {
                if (std::string_view(event.stage) == "secular_subdivision_handoff")
                    ++static_cast<Handoffs *>(p)->other;
            };
            cps_options options;
            options.goal = CPS_ROOT_APPROXIMATION;
            options.bits = 53;
            options.max_work_bits = cap;
            d::SolveState solver(certpolysolve::detail::dense(*input.get()), policy);
            auto result = solver.advance(&options);
            require(result.complete && result.goal_met && result.termination == CPS_SUCCESS && handoffs.other == 0,
                    "the exact fixture did not finish in secular alone");
            check_cover(roots, std::move(result), "a region of the exact cover has the wrong count");
        }
        std::cout << "inexact input tests passed\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "FAILED: " << error.what() << '\n';
        return 1;
    }
}
