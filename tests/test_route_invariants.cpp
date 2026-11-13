#include "solve/prepared_piecewise.hpp"
#include "solve/policies.hpp"
#include "numeric/arithmetic/polynomial_arithmetic.hpp"
#include "internal.hpp"
#include "exact_winding_api.hpp"
#include "methods/piecewise/piecewise_roots.hpp"
#include <iostream>

namespace d = certpolysolve::detail;
using namespace reference;
static Rational endpoint(const d::Endpoint &e) {
    return reference::endpoint(certpolysolve::Dyadic(e.mantissa_text(), e.exponent_text()));
}
static void verify(const cps_result &result, const std::vector<Complex> &known) {
    require(result.complete && result.goal_met && d::disjoint(result.roots), "route lost a disjoint complete cover");
    uint64_t total = 0;
    for (const auto &root : result.roots) {
        std::vector<Complex> inside;
        for (const auto &z : known)
            if (endpoint(root.endpoints[0]) <= z.real && z.real <= endpoint(root.endpoints[1]) &&
                endpoint(root.endpoints[2]) <= z.imag && z.imag <= endpoint(root.endpoints[3])) inside.push_back(z);
        require(inside.size() == root.count && !inside.empty(), "route certificate disagrees with exact planted roots");
        const bool same = std::all_of(inside.begin(), inside.end(), [&](const Complex &z) {
            return z.real == inside[0].real && z.imag == inside[0].imag;
        });
        if (root.multiplicity) require(same && root.multiplicity == inside.size(), "route fabricated multiplicity");
        if (root.realness == CPS_REAL)
            require(std::all_of(inside.begin(), inside.end(), [](const Complex &z) { return z.imag == 0; }),
                    "route fabricated realness");
        if (root.realness == CPS_NONREAL)
            require(std::none_of(inside.begin(), inside.end(), [](const Complex &z) { return z.imag == 0; }),
                    "route fabricated nonrealness");
        certpolysolve::Root box;
        box.box = {{certpolysolve::Dyadic(root.endpoints[0].mantissa_text(), root.endpoints[0].exponent_text()),
                    certpolysolve::Dyadic(root.endpoints[1].mantissa_text(), root.endpoints[1].exponent_text())},
                   {certpolysolve::Dyadic(root.endpoints[2].mantissa_text(), root.endpoints[2].exponent_text()),
                    certpolysolve::Dyadic(root.endpoints[3].mantissa_text(), root.endpoints[3].exponent_text())}};
        check_accuracy(box, 53);
        total += root.count;
    }
    require(total == known.size(), "route omitted or duplicated exact roots");
}
int main() {
    certpolysolve::ThreadScope lifetime;
    try {
        d::Root finite_zero;
        finite_zero.publish(64);
        require(finite_zero.endpoints[0].mantissa_text() == "0", "zero endpoint conversion failed");
        for (bool infinity : {false, true}) {
            d::Root bad;
            bad.publish(64);
            if (infinity) arb_pos_inf(acb_realref(bad.box.value()));
            else arb_indeterminate(acb_imagref(bad.box.value()));
            require(!d::disjoint(std::vector<d::Root>{bad}), "nonfinite singleton accepted as disjoint cover");
            bool refused = false;
            try { bad.publish(64); } catch (const std::logic_error &) { refused = true; }
            require(refused && bad.endpoints[0].mantissa_text() == "0", "nonfinite publication was not rejected before export");
        }
        const std::vector<Complex> base{{0}, {1}, {1}, {-2}, {0, 1}, {0, -1},
                                       {Rational(1)/3, Rational(2)/5}, {3, -2}};
        for (const auto &scale : {Complex(1), Complex(0, 2), Complex(Rational(3)/8, Rational(1)/4)}) {
            std::vector<Complex> known;
            for (const auto &z : base) known.push_back(z * scale + Complex(2, -1));
            const auto input = reference::input(from_roots(known));
            cps_options options;
            options.goal = CPS_ROOT_APPROXIMATION;
            options.max_work_bits = 2048;
            for (const auto lane : {d::StudyMethod::None, d::StudyMethod::Aberth, d::StudyMethod::Secular}) {
                auto policy = d::production_policy();
                policy.study = lane;
                // A one-unit solve stops honestly, a cancelled one says so, and the whole solve is right: three
                // fresh solves (one advance a solve state since 18 October 2026; until then a cancelled later call
                // had to keep the certified output of the call before it).
                options.work_units = 1;
                d::SolveState unit(certpolysolve::detail::dense(*input.get()), policy);
                const auto partial = unit.advance(&options);
                require(partial.goal_met || partial.termination == CPS_WORK_LIMIT, "a one-unit route solve stopped otherwise");
                options.cancel = [](void *) { return 1; };
                d::SolveState cancelled_state(certpolysolve::detail::dense(*input.get()), policy);
                require(cancelled_state.advance(&options).termination == CPS_CANCELLED, "a route ignored cancellation");
                options.cancel = nullptr;
                options.work_units = 100000;
                d::SolveState whole(certpolysolve::detail::dense(*input.get()), policy);
                verify(whole.advance(&options), known);
            }
            d::PiecewiseRootConfig config;
            config.reuse_guards = config.reuse_polish = config.blocked_polish = true;
            config.conjugate_symmetry = config.recover_discovery = config.durable_ownership = true;
            config.local_band_error = config.aligned_polish = true;
            d::PiecewiseRootStatistics statistics;
            verify(d::solve_prepared_piecewise_roots(certpolysolve::detail::dense(*input.get()), options, statistics, config), known);
        }
        std::cout << "exact planted-root route, rotation, scale, multiplicity, realness and resume checks passed\n";
    } catch (const std::exception &e) { std::cerr << e.what() << '\n'; return 1; }
}
