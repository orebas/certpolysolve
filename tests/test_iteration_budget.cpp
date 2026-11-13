// The iteration budget of a solve (Options::iteration_budget): at most N
// sweeps of the simultaneous iterations, counted by
// Result::statistics().sweeps; a spent budget returns the last certified
// cover with CPS_WORK_LIMIT and CPS_STOP_ITERATION_BUDGET; zero sets no limit and
// a budget the solve never reaches changes nothing. (Until 18 October 2026 a
// session's next call resumed with a fresh budget; the 0.5 API has one-shot
// solves only.) Checked through the public API on every method that
// sweeps the input (classical, secular, Auto), and on the piecewise method,
// whose own work the budget does not count.
#include "api_test_support.hpp"
#include "exact_winding_api.hpp"
#include <certpolysolve/certpolysolve.hpp>
#include <iostream>
#include <string>

using namespace reference;

namespace {
// Twelve conjugate pairs with three-decimal coordinates: moduli from 0.27 to
// 2.6, the closest two roots 0.19 apart. Binary64 needs a few dozen sweeps
// from its start on this polynomial.
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
// Every region holds exactly its count of the planted roots; a complete
// result's counts sum to the degree.
void check_cover(const certpolysolve::Result &result, const std::vector<Complex> &roots, const char *message) {
    uint64_t total = 0;
    for (const auto &region : result.roots()) {
        const Rational real_lower = endpoint(region.box.real.lo), real_upper = endpoint(region.box.real.hi);
        const Rational imag_lower = endpoint(region.box.imag.lo), imag_upper = endpoint(region.box.imag.hi);
        uint64_t inside = 0;
        for (const auto &root : roots)
            inside += real_lower <= root.real && root.real <= real_upper && imag_lower <= root.imag &&
                      root.imag <= imag_upper;
        require(inside == region.count, message);
        total += region.count;
    }
    require(!result.complete() || total == roots.size(), "a complete cover under an iteration budget lost roots");
}
std::string signature(const certpolysolve::Result &result) {
    std::string out;
    for (const auto &region : result.roots())
        out += region.box.real.lo.mantissa() + "@" + region.box.real.lo.exponent() + " " +
               region.box.real.hi.mantissa() + "@" + region.box.real.hi.exponent() + " " +
               region.box.imag.lo.mantissa() + "@" + region.box.imag.lo.exponent() + " " +
               region.box.imag.hi.mantissa() + "@" + region.box.imag.hi.exponent() + " " +
               std::to_string(region.count) + ";";
    return out;
}
certpolysolve::Options query(certpolysolve::Method method, uint64_t sweeps) {
    certpolysolve::Options options;
    options.method = method;
    options.goal = certpolysolve::Goal::approximation();
    options.bits = 53;
    options.max_bits = 4096;
    options.iteration_budget = sweeps;
    return options;
}
} // namespace

int main() {
    certpolysolve::ThreadScope lifetime;
    try {
        const auto roots = planted();
        const auto input = reference::input(from_roots(roots));
        // The sweeps of a call that spends a small budget, on every method
        // that sweeps the input. Since the compensated tier's Newton retry
        // (28 October 2026) Auto's certificate attempt proves this fixture's
        // roots after five sweeps; three leave it unfinished.
        constexpr uint64_t budget = 3;
        for (const auto method : {certpolysolve::Method::classical, certpolysolve::Method::secular, certpolysolve::Method::automatic}) {
            const auto unlimited = certpolysolve::solve(input, query(method, 0));
            require(unlimited.goal_met() && cps_result_termination(unlimited.get()) == CPS_SUCCESS,
                    "the unlimited solve of the fixture did not finish");
            const uint64_t needed = unlimited.statistics().sweeps;
            require(needed > 2 * budget, "the fixture finishes within the test's budget: the test tests nothing");

            // One call: exactly the budget's sweeps, the new stop, a sound cover.
            const auto stopped = certpolysolve::solve(input, query(method, budget));
            require(!stopped.goal_met() && cps_result_termination(stopped.get()) == CPS_WORK_LIMIT &&
                        cps_result_stop_reason(stopped.get()) == CPS_STOP_ITERATION_BUDGET,
                    "a spent iteration budget was not reported as CPS_WORK_LIMIT / CPS_STOP_ITERATION_BUDGET");
            require(stopped.statistics().sweeps == budget,
                    "a call under an iteration budget did not make exactly the sweeps the budget allows");
            check_cover(stopped, roots, "a region returned at the iteration budget has a wrong count");

            // A budget the solve never reaches changes nothing.
            const auto generous = certpolysolve::solve(input, query(method, needed + 1));
            require(generous.goal_met() && generous.statistics().sweeps == needed &&
                        signature(generous) == signature(unlimited),
                    "an iteration budget above the solve's own sweeps changed the solve");
            // One sweep short: the call makes no more than that (it stops
            // there, unless the solve's last sweep was one it did not need).
            const auto short_by_one = certpolysolve::solve(input, query(method, needed - 1));
            require(short_by_one.statistics().sweeps <= needed - 1 &&
                        (short_by_one.goal_met() ||
                         (cps_result_stop_reason(short_by_one.get()) == CPS_STOP_ITERATION_BUDGET &&
                          short_by_one.statistics().sweeps == needed - 1)),
                    "a budget one sweep short of the solve was exceeded or misreported");
        }
        // The piecewise method's own work is not counted: a forced piecewise
        // solve under a small budget still answers honestly (it finishes, or
        // it stops for a budget that is named).
        {
            const auto result = certpolysolve::solve(input, query(certpolysolve::Method::piecewise, budget));
            check_cover(result, roots, "a region of the piecewise solve under an iteration budget has a wrong count");
            const auto stop = cps_result_stop_reason(result.get());
            require(result.goal_met() ? cps_result_termination(result.get()) == CPS_SUCCESS && stop == CPS_STOP_NONE
                                      : cps_result_termination(result.get()) == CPS_WORK_LIMIT && stop == CPS_STOP_ITERATION_BUDGET,
                    "a piecewise solve under an iteration budget ended with another stop");
        }
        // Native secular input runs its own loop under the same budget: a solve under each budget up to the
        // sweeps it needs makes no more sweeps than allowed and stops honestly.
        {
            const certpolysolve::SecularEquation secular{
                {.residue = certpolysolve::Number("36/5"), .pole = certpolysolve::Number("-2")},
                {.residue = certpolysolve::Number("-4/3"), .pole = certpolysolve::Number("0")},
                {.residue = certpolysolve::Number("2/15"), .pole = certpolysolve::Number("3")}};
            require(secular.ok(), "the native secular fixture was refused");
            const std::vector<Complex> known{Complex(Rational(1)), Complex(Rational(2)), Complex(Rational(4))};
            certpolysolve::Options unlimited;
            unlimited.bits = 53;
            unlimited.max_bits = 1024;
            const auto reference_result = certpolysolve::solve(secular, unlimited);
            require(reference_result.goal_met(), "the native secular fixture did not finish");
            const uint64_t needed = reference_result.statistics().sweeps;
            require(needed > 2, "the native secular fixture needs no sweeps: the test tests nothing");
            for (uint64_t budget = 1; budget <= needed; ++budget) {
                certpolysolve::Options limited;
                limited.bits = 53;
                limited.max_bits = 1024;
                limited.iteration_budget = budget;
                const auto result = certpolysolve::solve(secular, limited);
                require(result.statistics().sweeps <= budget,
                        "a native secular solve exceeded its iteration budget");
                check_cover(result, known, "a region of a native secular solve under a budget has a wrong count");
                if (!result.goal_met())
                    require(cps_result_termination(result.get()) == CPS_WORK_LIMIT &&
                                cps_result_stop_reason(result.get()) == CPS_STOP_ITERATION_BUDGET,
                            "an unfinished native secular solve under an iteration budget reported another stop");
            }
        }
        // A solve that meets the goal says so and nothing else: termination
        // CPS_SUCCESS with stop reason CPS_STOP_NONE, even when the step that met
        // it spent the budget (the D1 fuzz family, 16 October 2026, seed
        // 9894797150310748061: a real-only session under a one-sweep budget on
        // an input without real roots ended goal met with the budget's stop
        // reason). Here (x^2 + 1)(x^2 + 4)(x^2 + 9), real-only, under every
        // budget up to the sweeps the solve needs: the last one meets the goal
        // in the sweep that spends it.
        {
            Polynomial p{1};
            for (const Rational a : {Rational(1), Rational(4), Rational(9)}) p = multiply(p, {a, 0, 1});
            certpolysolve::Options unlimited;
            unlimited.bits = 53;
            unlimited.max_bits = 1024;
            const auto reference_result = certpolysolve::real_roots(reference::input(p), unlimited);
            require(reference_result.goal_met(), "the real-only fixture without real roots did not finish");
            const uint64_t needed = reference_result.statistics().sweeps;
            require(needed >= 1, "the real-only fixture needs no sweeps: the test tests nothing");
            for (uint64_t budget = 1; budget <= needed; ++budget) {
                certpolysolve::Options limited;
                limited.bits = 53;
                limited.max_bits = 1024;
                limited.iteration_budget = budget;
                const auto result = certpolysolve::real_roots(reference::input(p), limited);
                if (result.goal_met()) {
                    require(result.complete() && result.roots().empty(),
                            "a real-only query on an input without real roots answered with regions");
                    require(cps_result_termination(result.get()) == CPS_SUCCESS && cps_result_stop_reason(result.get()) == CPS_STOP_NONE,
                            "a solve that met the goal reported a stop reason");
                } else {
                    require(cps_result_termination(result.get()) == CPS_WORK_LIMIT &&
                                cps_result_stop_reason(result.get()) == CPS_STOP_ITERATION_BUDGET,
                            "an unfinished real-only solve under an iteration budget reported another stop");
                }
            }
        }
        // A real-only query on exact real input is answered by the complex
        // route (since 12 October 2026), whose sweeps are under the budget: a
        // solve under each budget up to the sweeps it needs makes no more
        // sweeps than allowed, stops honestly, and the last one answers.
        {
            const std::vector<Complex> known{Complex(Rational(-3)), Complex(Rational(-1)), Complex(Rational(2)),
                                             Complex(Rational(4)),  Complex(Rational(7)),  Complex(Rational(11))};
            const auto real_input = reference::input(from_roots(known));
            certpolysolve::Options unlimited;
            unlimited.bits = 53;
            unlimited.max_bits = 1024;
            const auto reference_result = certpolysolve::real_roots(real_input, unlimited);
            require(reference_result.goal_met(), "the real-only fixture did not finish");
            const uint64_t needed = reference_result.statistics().sweeps;
            for (uint64_t budget = 1; budget <= needed; ++budget) {
                certpolysolve::Options limited;
                limited.bits = 53;
                limited.max_bits = 1024;
                limited.iteration_budget = budget;
                const auto result = certpolysolve::real_roots(real_input, limited);
                require(result.statistics().sweeps <= budget,
                        "a real-only solve made more sweeps than its iteration budget");
                if (result.goal_met())
                    check_cover(result, known,
                                "a region of the real-only solve under an iteration budget has a wrong count");
                else
                    require(cps_result_termination(result.get()) == CPS_WORK_LIMIT &&
                                cps_result_stop_reason(result.get()) == CPS_STOP_ITERATION_BUDGET,
                            "an unfinished real-only solve under an iteration budget reported another stop");
                if (budget == needed)
                    require(result.goal_met(), "a real-only solve under the sweeps it needs did not answer");
            }
        }
        // The axis fuzz family's first finding (16 October 2026, seed 12355245785107802074, minimized): roots
        // -i 2^-530 and 262/143 as the family wrote them, with a Gaussian coefficient, under the cluster goal with
        // radius 2^-8. Binary64 cannot hold the small root: the study's binary64 attempt fails after its sweeps,
        // and the scaled arithmetic's retry, like the adaptive packet after it, must get what the call's budget
        // leaves, not the packet's length again (two sweeps under a budget of one before the fix). Auto's route
        // (the secular route's preliminary packet), one solve under budgets of one and four sweeps.
        // The family's classical findings (seeds 10270399729761490471, 13605017879532207035) go through the
        // same study; the forced classical recipe refuses this query's cluster goal, so it is not repeated here.
        {
            certpolysolve::Polynomial finding;
            finding.set(0, certpolysolve::Number::complex(certpolysolve::Number("0"), certpolysolve::Number("-262")));
            finding.set(1, certpolysolve::Number::complex(certpolysolve::Number("-92087141732056050960653210079197441173739112890242356083781090626898550479369517345592864043689015"
                           "9822874436395298808220531691958233989467593930451083338226597888"), certpolysolve::Number("143")));
            finding.set(2, certpolysolve::Number::complex(certpolysolve::Number("50261302548412272089211484890554328579559897493529224885422503662772872971564278551220532665066905"
                           "6697217726734838662502045923473387253793381420055362280024440832"), certpolysolve::Number("0")));
            require(finding.ok(), "the fuzz finding was refused as input");
            const Rational tiny = Rational(1) / Rational(Integer(1) << 530);
            const std::vector<Complex> roots{Complex(Rational(0), -tiny), Complex(Rational("262/143"), Rational(0))};
            const auto budgeted = [](certpolysolve::Method method, uint64_t sweeps) {
                certpolysolve::Options options;
                options.method = method;
                options.bits = 53;
                options.max_bits = 4096;
                options.work_budget = 100000;
                options.cell_budget = 65536;
                options.goal = certpolysolve::Goal::clusters(certpolysolve::Number("1@-8"));
                options.iteration_budget = sweeps;
                return options;
            };
            for (const auto method : {certpolysolve::Method::automatic})
                for (const uint64_t budget : {uint64_t(1), uint64_t(4)}) {
                    const auto result = certpolysolve::solve(finding, budgeted(method, budget));
                    const uint64_t sweeps = result.statistics().sweeps;
                    require(sweeps <= budget, "a call whose binary64 arithmetic failed exceeded its iteration budget");
                    check_cover(result, roots, "a region returned after an arithmetic failure has a wrong count");
                    if (!result.goal_met())
                        require(cps_result_termination(result.get()) == CPS_WORK_LIMIT &&
                                    cps_result_stop_reason(result.get()) == CPS_STOP_ITERATION_BUDGET && sweeps == budget,
                                ("a call stopped after an arithmetic failure misreported its iteration budget (method " +
                                 std::to_string(static_cast<int>(method)) + ", budget " + std::to_string(budget) + ": termination " +
                                 std::to_string(cps_result_termination(result.get())) + ", stop " +
                                 std::to_string(cps_result_stop_reason(result.get())) + ", sweeps " + std::to_string(sweeps) +
                                 ", complete " + std::to_string(result.complete()) + ")")
                                    .c_str());
                }
        }
        std::cout << "iteration budget tests passed\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "FAILED: " << error.what() << '\n';
        return 1;
    }
}
