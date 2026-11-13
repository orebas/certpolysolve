// Public isolation is a count/multiplicity contract without an accuracy target.
// Check the exported rectangles using exact planted roots or independent
// Sturm--Cauchy counts; never use a second numerical solver as the oracle.
#include "api_test_support.hpp"
#include "exact_winding_api.hpp"
#include <algorithm>
#include <chrono>
#include <iostream>
#include <stop_token>
#include <vector>

using reference::Complex;
using reference::ComplexPolynomial;
using reference::Rational;
using reference::require;

namespace {
struct KnownRoot {
    Complex value;
    uint64_t multiplicity = 1;
};

certpolysolve::Options isolation_options() {
    certpolysolve::Options options;
    options.goal = certpolysolve::Goal::isolation();
    options.max_bits = 2048;
    options.time_limit = std::chrono::seconds(10);
    return options;
}

ComplexPolynomial coefficients(const std::vector<KnownRoot> &known) {
    std::vector<Complex> roots;
    for (const auto &root : known)
        for (uint64_t j = 0; j < root.multiplicity; ++j) roots.push_back(root.value);
    return reference::from_roots(roots);
}

bool contains(const reference::Rectangle &box, const Complex &point) {
    return box.real_lower <= point.real && point.real <= box.real_upper &&
           box.imag_lower <= point.imag && point.imag <= box.imag_upper;
}

void goal_and_values(const certpolysolve::Result &result) {
    require(result.get() && cps_result_goal(result.get()) == CPS_ROOT_ISOLATION,
            "public result lost the isolation goal");
    require(result.values().empty(), "isolation published root approximations");
    for (size_t j = 0; j < result.roots().size(); ++j) {
        require(!result.roots()[j].value, "isolation region has a C++ approximation");
        double real = 17, imag = 19;
        require(cps_result_get_approximation(result.get(), j, &real, &imag) == CPS_UNRESOLVED &&
                    real == 17 && imag == 19,
                "C approximation getter accepted isolation or changed refused outputs");
    }
    require(result.statistics().sweeps == 0, "complex isolation ran simultaneous iteration");
}

void planted_cover(const certpolysolve::Result &result, const std::vector<KnownRoot> &known, bool attained,
                   bool reported_complete = true) {
    goal_and_values(result);
    require(result.complete() == reported_complete, "unexpected complete flag");
    require(result.goal_met() == attained, "unexpected isolation attainment");
    reference::check_disjoint(result.roots());
    std::vector<unsigned> owners(known.size());
    uint64_t degree = 0;
    for (const auto &root : known) degree += root.multiplicity;
    require(result.degree() == static_cast<int64_t>(degree), "input degree changed");
    for (const auto &root : result.roots()) {
        const auto box = reference::rectangle(root);
        uint64_t count = 0, multiplicity = 0;
        size_t distinct = 0;
        for (size_t j = 0; j < known.size(); ++j)
            if (contains(box, known[j].value)) {
                ++owners[j];
                ++distinct;
                count += known[j].multiplicity;
                multiplicity = known[j].multiplicity;
            }
        require(count > 0 && count == root.count, "planted-root rectangle count is wrong");
        if (root.multiplicity)
            require(distinct == 1 && *root.multiplicity == multiplicity,
                    "a cluster was given a root multiplicity");
        if (attained)
            require(distinct == 1 && root.multiplicity.value_or(0) == root.count,
                    "successful isolation retained a cluster");
    }
    require(std::all_of(owners.begin(), owners.end(), [](unsigned n) { return n == 1; }),
            "a planted root was lost or counted twice");
}

void winding_cover(const ComplexPolynomial &polynomial, const certpolysolve::Result &result,
                   std::vector<uint64_t> multiplicities) {
    goal_and_values(result);
    require(result.ok() && result.complete(), "irrational-root isolation did not finish");
    reference::check_disjoint(result.roots());
    uint64_t total = 0;
    std::vector<uint64_t> observed;
    for (const auto &root : result.roots()) {
        const auto counted = reference::rectangle_count(polynomial, root);
        require(counted && *counted > 0 && static_cast<uint64_t>(*counted) == root.count,
                "independent exact winding count differs");
        require(root.multiplicity.value_or(0) == root.count, "isolated multiplicity is unknown");
        total += root.count;
        observed.push_back(root.count);
    }
    std::sort(observed.begin(), observed.end());
    std::sort(multiplicities.begin(), multiplicities.end());
    require(observed == multiplicities && total + 1 == polynomial.size(),
            "isolated roots do not have the known multiplicities");
}

void planted_roots() {
    const Rational delta = reference::power2(-96);
    const std::vector<std::vector<KnownRoot>> fixtures{
        {{Complex(Rational(13) / 7), 1}},
        {{Complex(Rational(1) / 3), 1}, {Complex(Rational(-7) / 5), 1}, {Complex(2), 1}},
        {{Complex(0), 3}, {Complex(Rational(1) / 3, Rational(1) / 5), 2},
         {Complex(-2, 3), 1}, {Complex(2, -1), 1}},
        // Different multiplicities force separate exact factors. Their first
        // individually isolating boxes need not be disjoint from each other.
        {{Complex(Rational(1) / 3), 1}, {Complex(Rational(1) / 3 + delta), 2}},
        {{Complex(0), 5}}
    };
    for (const auto &known : fixtures) {
        const auto result = certpolysolve::solve(reference::input(coefficients(known)), isolation_options());
        if (!result.ok()) throw std::runtime_error("planted isolation stopped: " + result.message());
        planted_cover(result, known, true);
    }
}

void irrational_roots_and_no_width() {
    for (const ComplexPolynomial polynomial : {ComplexPolynomial{{-2}, {0}, {1}},
                                              ComplexPolynomial{{0, -1}, {0}, {1}}}) {
        const auto input = reference::input(polynomial);
        auto options = isolation_options();
        options.bits = 1;
        const auto low = certpolysolve::solve(input, options);
        winding_cover(polynomial, low, {1, 1});
        options.bits = 8192; // exceeds max_bits: isolation must ignore the accuracy setting
        options.absolute_tolerance = certpolysolve::Number("1e100");
        const auto high = certpolysolve::solve(input, options);
        winding_cover(polynomial, high, {1, 1});
        require(low.roots().size() == high.roots().size(), "ignored accuracy settings changed isolation");
        bool wide = false;
        for (size_t j = 0; j < low.roots().size(); ++j) {
            require(low.roots()[j].box == high.roots()[j].box &&
                        low.roots()[j].count == high.roots()[j].count &&
                        low.roots()[j].multiplicity == high.roots()[j].multiplicity,
                    "requested bits or absolute tolerance changed isolated regions");
            const auto box = reference::rectangle(low.roots()[j]);
            const Rational width = box.real_upper - box.real_lower + box.imag_upper - box.imag_lower;
            const Rational scale = reference::absolute(box.real_lower) + reference::absolute(box.real_upper) +
                                   reference::absolute(box.imag_lower) + reference::absolute(box.imag_upper);
            wide = wide || width > reference::power2(-CPS_DEFAULT_BITS) * scale;
        }
        require(wide, "isolation unnecessarily refined every box to the default accuracy");
        const auto a = low.statistics(), b = high.statistics();
        const auto ca = low.engine_counters(), cb = high.engine_counters();
        require(a.work_units == b.work_units && a.highest_work_bits == b.highest_work_bits &&
                    a.precision_increases == b.precision_increases && ca.disk_tests == cb.disk_tests &&
                    ca.taylor_shifts == cb.taylor_shifts && ca.graeffe_steps == cb.graeffe_steps &&
                    ca.boxes_split == cb.boxes_split,
                "ignored accuracy settings changed isolation work");
    }
    // The nonzero quadratic has two well-separated roots, but one is very
    // close to the separately extracted double zero. Final assembly must
    // refine against that zero as well as against other nonzero factors.
    const Rational delta = reference::power2(-80);
    const ComplexPolynomial near_zero{{0}, {0}, {-delta}, {-2}, {1}};
    const auto result = certpolysolve::solve(reference::input(near_zero), isolation_options());
    winding_cover(near_zero, result, {1, 1, 2});
    size_t zeros = 0;
    for (const auto &root : result.roots())
        if (contains(reference::rectangle(root), Complex(0))) {
            require(root.count == 2 && root.multiplicity == 2, "extracted zero multiplicity changed");
            ++zeros;
        }
    require(zeros == 1, "the extracted zero is not isolated exactly once");
}

void stops() {
    const std::vector<KnownRoot> known{{Complex(-1), 1}, {Complex(1), 1}, {Complex(0, 1), 1}};
    const auto input = reference::input(coefficients(known));
    auto options = isolation_options();
    options.work_budget = 1;
    auto result = certpolysolve::solve(input, options);
    require(result.stop() == certpolysolve::Stop::work_budget, "one work unit did not stop honestly");
    planted_cover(result, known, false);

    options = isolation_options();
    options.cell_budget = 1;
    result = certpolysolve::solve(input, options);
    require(result.stop() == certpolysolve::Stop::cell_budget, "one cell did not stop honestly");
    planted_cover(result, known, false);

    const Rational delta = reference::power2(-96);
    const std::vector<KnownRoot> close{{Complex(1 - delta), 1}, {Complex(1 + delta), 1}};
    options = isolation_options();
    options.max_bits = CPS_MIN_MAX_WORK_BITS;
    result = certpolysolve::solve(reference::input(coefficients(close)), options);
    require(result.stop() == certpolysolve::Stop::precision_cap, "tight precision cap did not stop honestly");
    planted_cover(result, close, false);

    options = isolation_options();
    std::stop_source cancellation;
    cancellation.request_stop();
    options.stop = cancellation.get_token();
    result = certpolysolve::solve(input, options);
    require(result.stop() == certpolysolve::Stop::cancelled && !result.goal_met(),
            "a pre-cancelled query claimed isolation");
    goal_and_values(result);

    options = isolation_options();
    options.resources.emplace();
    options.resources->exact_materialization_bits = 0;
    const std::vector<KnownRoot> repeated{{Complex(1), 2}, {Complex(2), 1}};
    result = certpolysolve::solve(reference::input(coefficients(repeated)), options);
    require(result.stop() == certpolysolve::Stop::unsupported && !result.complete() &&
                !result.goal_met() && result.preparation_limited(),
            "exact preparation refusal violates the result contract");
    // I4 marks an unsupported query incomplete, while keeping the regions
    // certified before refusal. The original counted cover must survive.
    planted_cover(result, repeated, false, false);
    const std::vector<KnownRoot> zero{{Complex(0), 5}};
    result = certpolysolve::solve(reference::input(coefficients(zero)), options);
    require(result.ok(), "pure zero roots required exact factor preparation");
    planted_cover(result, zero, true);
}

void unsupported_and_empty_inputs() {
    const auto options = isolation_options();
    const auto usage = [](const certpolysolve::Result &r) {
        require(r.stop() == certpolysolve::Stop::error && !r.get(), "unsupported option combination was not a usage error");
    };
    const certpolysolve::Polynomial exact{-2, 0, 1};
    for (const auto method : {certpolysolve::Method::classical, certpolysolve::Method::secular, certpolysolve::Method::piecewise}) {
        auto forced = options;
        forced.method = method;
        usage(certpolysolve::solve(exact, forced));
    }
    certpolysolve::Polynomial uncertain;
    uncertain.set(0, -2).set(2, certpolysolve::ball(1, "1/1000"));
    auto result = certpolysolve::solve(uncertain, options);
    require(result.stop() == certpolysolve::Stop::unsupported && !result.complete() &&
                result.statistics().work_units == 0,
            "ball coefficients were admitted to exact complex isolation");
    const certpolysolve::Polynomial gaussian{certpolysolve::Number::complex(0, -1), 0, 1};
    result = certpolysolve::real_roots(gaussian, options);
    require(result.stop() == certpolysolve::Stop::unsupported, "real isolation admitted nonreal coefficients");
    const certpolysolve::SecularEquation secular{{.residue = 1, .pole = 1}, {.residue = 1, .pole = 2}};
    usage(certpolysolve::solve(secular, options));
    const auto program = certpolysolve::Program::from([](auto z) { return z * z - 2; });
    usage(certpolysolve::solve(program, options));
    certpolysolve_test::COptions region(options);
    require(region.ok() && cps_options_set_region(region.get(), "0", "1", nullptr, nullptr) == CPS_OK,
            "could not configure a region");
    usage(certpolysolve_test::solve(exact, region.get()));

    result = certpolysolve::solve(certpolysolve::Polynomial(3), options);
    require(result.ok() && result.complete() && result.roots().empty(), "constant input has a root");
    goal_and_values(result);
    result = certpolysolve::solve(certpolysolve::Polynomial(), options);
    require(result.stop() == certpolysolve::Stop::zero_polynomial && !result.goal_met(),
            "zero polynomial was given a finite root enumeration");
    goal_and_values(result);
}
} // namespace

int main() {
    certpolysolve::ThreadScope lifetime;
    try {
        planted_roots();
        irrational_roots_and_no_width();
        stops();
        unsupported_and_empty_inputs();
        std::cout << "public complex isolation exact-root and contract checks passed\n";
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
