// The Auto selector: which maintained recipe answers an Auto query, decided from the query and the input alone,
// before any root work (select_auto_route; AutoRoute and AutoSelection are declared in auto_selection.hpp).
//
// The decisions, in their order: a query outside the explicit methods (real-only, or a goal other than approximation,
// accuracy and the cluster radius) and inexact coefficients keep the D5 route (selection_prefix); every other query
// takes the maintained secular recipe ("global"), unless piecewise is cheaper. A precision cap below piecewise's first
// working precision decides for secular before any estimate. Otherwise piecewise's own admission (admit_piecewise,
// piecewise_admission.hpp) runs with its development caps lifted (lifted_policy), so only its structural refusals and
// its work estimate remain, judged at the library's default goal or above (judged_query): an estimate below the global
// candidate stage's minimum sweep budget routes to piecewise (route_by_admission), which then begins with the
// certificate attempt when its estimate is at least the golden ratio times the attempt's cost. Every comparison is of
// exact integer counts or of binary64 ratios, so the route is the same on every platform.
//
// State: none. Each call is a function of its arguments; its own time goes to AutoSelection::elapsed_ns.
// Entry points: select_auto_route on a coefficient list (DensePolynomial), and on the input object's terms
// (cps_polynomial), which answers only where the decision needs no coefficient list (C40's second step). Callers:
// session.cpp (Auto's dispatch, and forced piecewise's regime exit, which asks where Auto would go), and
// real_only_route.cpp (the complex route of a real-only query, without the square-free proof).
// Not here: running a route, the admission's estimate (piecewise_admission.cpp), the D5 recipe.
// Design notes: docs/selector-campaign-20260928.md (the selector), docs/same-routes-every-platform-note-20261021.md
// (its decisions in integers and binary64), docs/d5-route-note-20260930.md (the cluster-radius goal).
#include "solve/auto_selection.hpp"
#include "numeric/evaluate/disk_sparse_polynomial_evaluator.hpp"
#include "core/heuristics.hpp"
#include "methods/piecewise/piecewise_admission.hpp"
#include "methods/piecewise/piecewise_config.hpp"
#include "core/telemetry.hpp"
#include <algorithm>
#include <limits>
#include <numbers>

namespace certpolysolve::detail {
namespace {
// The selector's decisions before piecewise's admission, in their order, from the query and whether every coefficient
// is exact (`all_exact`, called only when the query gets that far): shared by the coefficient list and the input
// object's terms (C40's second step). True when one of them decided; otherwise the route is Secular ("global") and
// the admission comes next.
template <class AllExact> bool selection_prefix(AutoSelection &result, const cps_options &options, AllExact &&all_exact) {
    // The cluster-radius goal on the maintained recipes (docs/d5-route-note-20260930.md, section 2): a cluster-radius
    // query on exact input enters the maintained secular route. Piecewise's
    // admission still refuses this goal; real-only queries and inexact coefficients keep D5.
    // History: docs/code-history-20261021.md, auto_selection.cpp 1.
    const bool explicit_goal = engine_goal(options) == CPS_ROOT_ACCURACY || engine_goal(options) == CPS_ROOT_APPROXIMATION ||
                               options.goal == CPS_CLUSTER_RADIUS;
    if (options.real_only || !explicit_goal) {
        result.reason = "query-outside-explicit-methods";
        return true;
    }
    if (!all_exact()) {
        result.reason = "inexact-coefficients";
        return true;
    }
    result.route = AutoRoute::Secular;
    result.reason = "global";
    // Piecewise refuses a precision cap below its first working precision
    // before any work; the route would only hand the query back.
    if (options.max_work_bits < static_cast<uint64_t>(piecewise_start_precision(options.bits))) {
        result.reason = "precision-cap-below-piecewise-start";
        return true;
    }
    return false;
}
// The query the admission judges and its policy (see select_auto_route).
cps_options judged_query(const cps_options &options) {
    cps_options judged = options;
    judged.bits = std::max(options.bits, cps_options{}.bits);
    return judged;
}
PiecewiseAdmissionPolicy lifted_policy(bool prove_squarefree) {
    PiecewiseAdmissionPolicy lifted;
    lifted.maximum_envelope_bend = std::nullopt;
    lifted.maximum_work_ratio = std::nullopt;
    lifted.maximum_global_sweeps = heuristics::initial_sweep_allowance;
    lifted.prove_squarefree = prove_squarefree;
    return lifted;
}
} // namespace

namespace {
// The route from piecewise's admission (the selector's last decision): its estimate against the global sweep
// budget, and the certificate attempt's order. `exponents(zero_order)` lists the input's nonzero exponents less its
// zero order, from the coefficient list or the terms.
template <class Exponents> void route_by_admission(AutoSelection &result, const PiecewiseAdmission &admission,
                                                   Exponents &&exponents) {
    if (admission.selected) {
        result.route = AutoRoute::Piecewise;
        result.reason = "piecewise-estimate-below-global-sweep-budget";
        // The global recipes begin with a certificate attempt (a binary64
        // packet and two evaluation passes: machine_newton_cover) that
        // finishes a well-conditioned input in binary64's range outright,
        // and that can be tried before piecewise too. Whether it will
        // succeed is not known beforehand, so the order is chosen by cost.
        // W is the piecewise estimate and A the attempt allowance, priced by
        // heuristics::certificate_attempt_point_operations, both in global
        // sweeps. The empirical threshold W/A >= phi balances the proxies
        // 1 + A/W and W/A. It is not a competitive-ratio guarantee: W is an
        // estimate, and a successful attempt can cost a << A, so choosing
        // piecewise first can cost W/a, far above W/A. Historical motivation:
        // docs/code-history-20261021.md, auto_selection.cpp 2.
        //
        // In counts (docs/same-routes-every-platform-note-20261021.md): a global sweep is degree x (degree +
        // support) operations, the admission's work is W in sweeps times that, and the attempt's operations at a
        // point (heuristics::certificate_attempt_point_operations) are A times degree + support; so W >= phi A is
        // work >= phi x (operations at a point) x degree, a ratio against an irrational, compared in binary64.
        result.piecewise_operations = admission.work;
        // The attempt evaluates a source with at most a quarter of its
        // coefficients nonzero term by term, building the terms' powers
        // (prepared_piecewise.cpp; the evaluator's plan counts the products).
        uint64_t products = 0;
        if (admission.support > 2 && admission.support * 4 <= admission.degree) {
            products = DiskSparsePolynomialEvaluator::plan(exponents(admission.zero_order)).multiplications;
        }
        const double attempt = static_cast<double>(heuristics::certificate_attempt_point_operations(
                                   admission.degree, admission.support, products)) *
                               static_cast<double>(admission.degree);
        if (static_cast<double>(admission.work) >= std::numbers::phi * attempt) {
            result.attempt_first = true;
            result.reason = "certificate-attempt-before-piecewise";
        }
    }
}
} // namespace

std::optional<AutoSelection> select_auto_route(const cps_polynomial &input, const cps_options &options,
                                               bool prove_squarefree) {
    AutoSelection result;
    bool decided = false;
    {
        StageTimer timer(result.elapsed_ns);
        decided = selection_prefix(result, options, [&] {
            return std::all_of(input.terms.begin(), input.terms.end(), [](const auto &t) { return t.second.exact(); });
        });
        // Without the square-free proof the admission runs on the terms (admit_piecewise's terms form); with it,
        // only an admission its structure refuses is decided here, which leaves the route Secular ("global") as for
        // the coefficient list.
        if (!decided && !prove_squarefree) {
            route_by_admission(result, admit_piecewise(input, judged_query(options), lifted_policy(false)),
                               [&](uint64_t zero_order) {
                                   std::vector<slong> exponents;
                                   for (const auto &term : input.terms)
                                       exponents.push_back(static_cast<slong>(term.first - zero_order));
                                   return exponents;
                               });
            decided = true;
        } else if (!decided) {
            decided = admission_by_structure(input_structure_of_terms(input.terms), judged_query(options),
                                             lifted_policy(prove_squarefree))
                          .has_value();
        }
    }
    if (!decided)
        return std::nullopt;
    return result;
}

AutoSelection select_auto_route(const DensePolynomial &input, const cps_options &options, bool prove_squarefree) {
    AutoSelection result;
    // The timer's scope closes before the return, as in the other overload: a timer that ends at `return result;`
    // adds its time after the return value is formed, which reaches the caller only where the compiler builds that
    // value in the caller's storage (GCC and Clang do; MSVC copies, and the time read 0). History:
    // docs/code-history-20261021.md, auto_selection.cpp 3.
    {
        StageTimer timer(result.elapsed_ns);
        if (!selection_prefix(result, options, [&] {
                return std::all_of(input.coefficients.begin(), input.coefficients.end(),
                                   [](const Coefficient &c) { return c.exact(); });
            })) {
            // Piecewise's own admission, with its development caps lifted so that
            // only its structural refusals remain (power composites, square-freeness
            // not proved, unavailable geometry). Its work estimate counts complex
            // coefficient operations: model construction plus local solves over its
            // precision ladder, as work_ratio x degree x support.
            //
            // One global Aberth sweep costs degree evaluations of `support` terms plus
            // degree^2 root interactions. The global candidate stage budgets at least
            // initial_sweep_allowance such sweeps at a precision before escalating
            // (classical_isolation.cpp), so piecewise is preferred when its whole estimate is
            // below that minimum budget: work_ratio * degree * support <
            // initial_sweep_allowance * degree * (degree + support). The admission
            // applies that bound itself, and stops estimating once it is passed: the
            // estimate walks one ring per octave of the root window, which for
            // x^3 + 2^(2^30) x - 1 is 1.6e9 rings.
            //
            // The estimate, the sweep budget and the order rule below were calibrated
            // together at the library's default goal (5.3). Below it neither route
            // gets cheaper the way the estimate does: the estimate sums piecewise's
            // model rungs up to twice the goal's bits, so it drops rungs as the goal
            // falls, while piecewise's discovery at its first rung (16-bit models, at
            // least 128 bits of working precision) is most of its work at any goal
            // and the global recipes begin in binary64 at any goal, so judged at its
            // own goal the estimate would send queries to piecewise where the global
            // route is faster (the measurement: docs/code-history-20261021.md,
            // auto_selection.cpp 4). A goal below the default is judged at the default goal.
            route_by_admission(result, admit_piecewise(input, judged_query(options), lifted_policy(prove_squarefree)),
                               [&](uint64_t zero_order) {
                                   std::vector<slong> exponents;
                                   for (size_t j = zero_order; j < input.coefficients.size(); ++j)
                                       if (!input.coefficients[j].zero())
                                           exponents.push_back(static_cast<slong>(j - zero_order));
                                   return exponents;
                               });
        }
    }
    return result;
}
} // namespace certpolysolve::detail
