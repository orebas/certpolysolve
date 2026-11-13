#include "exact_reference_api.hpp"
#include "api_test_support.hpp"
#include <future>
#include <iostream>
#include <random>

using namespace reference;

static void real_certificate(const Polynomial &p, int bits = 128, bool real_only = false) {
    certpolysolve::Options options;
    options.bits = static_cast<uint64_t>(bits);
    options.max_bits = std::max<uint64_t>(4096, static_cast<uint64_t>(bits) + 128);
    const auto result = real_only ? certpolysolve::real_roots(input(p), options) : certpolysolve::solve(input(p), options);
    if (!result.complete() || !result.goal_met()) {
        std::string message = "real test did not complete: bits=" + std::to_string(bits) +
                              " real_only=" + std::to_string(real_only) + " coefficients=";
        for (const auto &c : p)
            message += c.str() + ",";
        throw std::runtime_error(message + "\n" + result.message());
    }
    const auto &enclosures = result.roots();
    check_disjoint(enclosures);
    uint64_t count = 0;
    for (const auto &root : enclosures) {
        if (root.realness != certpolysolve::Realness::real) {
            std::string message =
                "realness was not proved, degree=" + std::to_string(degree(p)) + " coefficients=";
            for (const auto &c : p)
                message += c.str() + ",";
            throw std::runtime_error(message + "\n" + result.message());
        }
        require(endpoint(root.box.imag.lo) == 0 && endpoint(root.box.imag.hi) == 0,
                "real root has nonzero imaginary enclosure");
        require(roots_closed(p, endpoint(root.box.real.lo), endpoint(root.box.real.hi)) ==
                    static_cast<int>(root.count),
                "independent exact Sturm count disagrees with certificate");
        require(root.count == root.multiplicity.value_or(0), "isolated root has unknown multiplicity");
        check_accuracy(root, bits);
        count += root.count;
    }
    require(count == static_cast<uint64_t>(degree(p)), "real construction lost a root");
}

static void numerical_families() {
    for (int bits : {8, 53, 128, 512, 4096}) {
        real_certificate({-2, 0, 1}, bits);
        real_certificate({-2, 0, 1}, bits, true);
    }
    real_certificate({Rational(1) / 3, Rational(-2) / 3, Rational(1) / 3});
    real_certificate(multiply(multiply({0, 1}, {0, 1}), multiply({-2, 0, 1}, {-2, 0, 1})));
    real_certificate(multiply(multiply({0, 1}, {0, 1}), multiply({-2, 0, 1}, {-2, 0, 1})), 128, true);
    Polynomial wilkinson{1};
    for (int i = 1; i <= 15; ++i)
        wilkinson = multiply(wilkinson, {-i, 1});
    real_certificate(wilkinson);
    std::reverse(wilkinson.begin(), wilkinson.end());
    real_certificate(wilkinson);
    const Rational gap = power2(-200);
    const auto close = multiply({-1, 1}, {-1 - gap, 1});
    real_certificate(close);
    auto adaptive = certpolysolve::solve(input(close));
    require(adaptive.goal_met() && adaptive.work_bits() > 85,
            "hard case did not demonstrate precision escalation");
    std::mt19937 generator(20260906);
    for (int sample = 0; sample < 32; ++sample) {
        Polynomial p{1};
        const int n = 2 + static_cast<int>(generator() % 9);
        for (int i = 0; i < n; ++i) {
            const Rational root = Rational(static_cast<int>(generator() % 25) - 12) / (1 + generator() % 5);
            p = multiply(p, {-root, 1});
        }
        if (sample % 2)
            for (auto &c : p)
                c *= -17; // Scaling must preserve counts.
        real_certificate(p, 80, sample % 3 == 0);
    }
}

static void complex_certificates() {
    for (int n : {3, 5, 7, 12}) {
        Polynomial p(static_cast<size_t>(n) + 1);
        p[0] = -2;
        p.back() = 1;
        certpolysolve::Options options;
        options.bits = 128;
        auto result = certpolysolve::solve(input(p), options);
        require(result.goal_met() && result.roots().size() == static_cast<size_t>(n),
                "complex family did not isolate");
        std::vector<Complex> cp;
        for (const auto &c : p)
            cp.emplace_back(c);
        const auto &enclosures = result.roots();
        check_disjoint(enclosures);
        for (const auto &root : enclosures) {
            check_accuracy(root, 128);
            if (root.realness == certpolysolve::Realness::real) {
                require(roots_closed(p, endpoint(root.box.real.lo), endpoint(root.box.real.hi)) == 1,
                        "independent real root count failed");
            } else {
                // The closed forms (degree 3 here) publish tight boxes whose root may sit off the centre.
                require(rouche_inside_rectangle(cp, root) || rouche_inside_rectangle_after_newton(cp, root),
                        "independent exact Rouche inclusion failed");
            }
        }
    }
    // Near-real conjugates must not be rounded onto the real axis.
    const Rational epsilon = power2(-200);
    auto near = certpolysolve::solve(input({1 + epsilon * epsilon, -2, 1}));
    require(near.goal_met() && near.roots().size() == 2, "near-real pair did not isolate");
    for (const auto &root : near.roots())
        require(root.realness == certpolysolve::Realness::not_real, "false real-root classification");
    auto real = certpolysolve::real_roots(input({1 + epsilon * epsilon, -2, 1}));
    require(real.goal_met() && real.roots().size() == 0, "dedicated real solver invented near-real roots");

    certpolysolve::Polynomial complex_linear;
    complex_linear.set(0, certpolysolve::Number::complex(certpolysolve::Number("-1/3"), certpolysolve::Number("-2/7"))).set(1, "1");
    auto linear = certpolysolve::solve(complex_linear);
    require(linear.goal_met(), "complex rational linear polynomial did not solve");
    const auto z = linear.roots()[0];
    require(endpoint(z.box.real.lo) <= Rational(1) / 3 && endpoint(z.box.real.hi) >= Rational(1) / 3,
            "real component lost");
    require(endpoint(z.box.imag.lo) <= Rational(2) / 7 && endpoint(z.box.imag.hi) >= Rational(2) / 7,
            "imag component lost");
}

static void uncertainty_and_limits() {
    certpolysolve::Polynomial interval;
    interval.set(0, certpolysolve::ball(certpolysolve::Number("-1"), certpolysolve::Number("1e-8"))).set(1, "1");
    certpolysolve::Options options;
    options.bits = 53;
    options.max_bits = 256;
    auto unresolved = certpolysolve::solve(interval, options);
    require(unresolved.complete() && !unresolved.goal_met(),
            "invented precision from uncertain coefficients");
    options.absolute_tolerance = certpolysolve::Number("1e-7");
    auto resolved = certpolysolve::solve(interval, options);
    require(resolved.complete() && resolved.goal_met(), "absolute tolerance not honored");
    const auto root = resolved.roots()[0];
    require(endpoint(root.box.real.lo) <= Rational("99999999/100000000"),
            "uncertainty lower endpoint excluded");
    require(endpoint(root.box.real.hi) >= Rational("100000001/100000000"),
            "uncertainty upper endpoint excluded");

    for (const auto &exponent : {"10000000", "-10000000", "100000000000000000000"}) {
        certpolysolve::Polynomial p;
        p.set(0, std::string("-1@") + exponent).set(1, "1");
        auto result = certpolysolve::solve(p);
        require(result.goal_met(), "compressed-exponent linear solve failed");
        auto r = result.roots()[0];
        require(r.box.real.lo.mantissa() == "1" && r.box.real.lo.exponent() == exponent,
                "compressed exponent was lost");
    }
    auto decimal = certpolysolve::solve(certpolysolve::Polynomial(std::vector<std::string>{"-1e10000", "1"}));
    require(decimal.goal_met(), "large decimal exponent failed");
    auto decimal_fraction =
        certpolysolve::solve(certpolysolve::Polynomial(std::vector<std::string>{"0.1", "-0.2", "0.1"}));
    require(decimal_fraction.goal_met() && decimal_fraction.roots()[0].multiplicity.value_or(0) == 2,
            "decimal exactness lost");

    // No 0.5 Options field for an exact allowance alone: C options, the answer as a certpolysolve::Result.
    const auto no_exact = certpolysolve_test::raw_options();
    require(certpolysolve_test::set_exact_allowance(no_exact.get(), 0) == CPS_OK, "exact allowance refused");
    require(cps_options_set_precision(no_exact.get(), 128) == CPS_OK && cps_options_set_max_work_bits(no_exact.get(), 128) == CPS_OK,
            "precision refused");
    auto cluster = certpolysolve_test::solve(certpolysolve::Polynomial{1, -2, 1}, no_exact.get());
    require(cluster.complete() && !cluster.goal_met(), "uncertified repeated root accepted");
    auto enclosure = cluster.roots()[0];
    require(enclosure.count == 2 && enclosure.multiplicity.value_or(0) == 0, "cluster count confused with multiplicity");
    require(endpoint(enclosure.box.real.lo) < 1 && endpoint(enclosure.box.real.hi) > 1,
            "global bound excludes actual root");

    certpolysolve::Polynomial possibly_more_zeros;
    possibly_more_zeros.set(0, "0").set(1, certpolysolve::ball(certpolysolve::Number("0"), certpolysolve::Number("1"))).set(2, "1");
    auto zero_family = certpolysolve_test::solve(possibly_more_zeros, no_exact.get());
    require(zero_family.complete() && !zero_family.goal_met(), "uncertain zero multiplicity was guessed");
    require(zero_family.roots().size() == 1 && zero_family.roots()[0].count == 2 &&
                zero_family.roots()[0].multiplicity.value_or(0) == 0,
            "zero deflation produced overlapping/double-counted certificates");

    // Dynamic range in a NONLINEAR problem, with no expanded integer polynomial.
    for (const char *exponent : {"2000000", "-2000000"}) {
        certpolysolve::Polynomial p;
        p.set(0, std::string("-1@") + exponent).set(2, "1");
        auto result = certpolysolve::solve(p);
        require(result.goal_met() && result.roots().size() == 2, "scaled quadratic did not isolate");
        for (const auto &r : result.roots()) {
            require(r.box.real.lo.mantissa().size() < 200,
                    "scaled quadratic expanded an exponent into its mantissa");
            require(r.realness == certpolysolve::Realness::real, "scaled quadratic realness unknown");
        }
    }

    // Two reciprocal roots with enormous separation in scale. Exact center +
    // radius is compact even when expanding their sum would be expensive.
    certpolysolve::Polynomial mixed_scale;
    mixed_scale.set(0, "1").set(1, certpolysolve::ball(certpolysolve::Number("-1@4096"), certpolysolve::Number("1@-4096"))).set(2, "1");
    certpolysolve::Options mixed_options;
    mixed_options.bits = 128;
    auto mixed = certpolysolve::solve(mixed_scale, mixed_options);
    require(mixed.goal_met() && mixed.roots().size() == 2, "mixed-scale uncertain polynomial did not isolate");
    check_disjoint(mixed.roots());
    for (const auto &r : mixed.roots()) {
        check_accuracy(r, 128);
        for (int sign : {-1, 1}) {
            const Polynomial endpoint_polynomial{1, -power2(4096) + sign * power2(-4096), 1};
            require(roots_closed(endpoint_polynomial, endpoint(r.box.real.lo), endpoint(r.box.real.hi)) == 1,
                    "mixed-scale coefficient endpoint excluded");
        }
    }
}

static void concurrent_contexts() {
    std::vector<std::future<void>> jobs;
    for (int i = 0; i < 4; ++i)
        jobs.push_back(std::async(std::launch::async, [i] {
            certpolysolve::ThreadScope thread_lifetime;
            real_certificate({-2 - i, 0, 1}, 64 + i * 32);
        }));
    for (auto &job : jobs)
        job.get();
}

static void local_cluster_certificates() {
    const auto limited = certpolysolve_test::raw_options();
    require(certpolysolve_test::set_exact_allowance(limited.get(), 0) == CPS_OK, "exact allowance refused");
    require(cps_options_set_precision(limited.get(), 128) == CPS_OK && cps_options_set_max_work_bits(limited.get(), 128) == CPS_OK,
            "precision refused");
    const auto p = multiply(multiply(multiply({-1, 1}, {-1, 1}), multiply({2, 1}, {2, 1})), {-4, 1});
    auto result = certpolysolve_test::solve(input(p), limited.get());
    require(result.complete() && !result.goal_met(), "unresolved multiplicity misreported");
    if (result.roots().size() != 3)
        throw std::runtime_error("did not recover three local clusters: " + result.message());
    const auto &enclosures = result.roots();
    check_disjoint(enclosures);
    for (const auto &root : enclosures) {
        // The count's certificate family depends on the route Auto takes (D5's
        // Pellet counts, or an explicit method's cover); the independent Sturm
        // count below checks the claim itself. test_auto_selection keeps D5's
        // Pellet path covered through an explicit D5 session.
        require(!root.certificate.empty(), "local cluster lacks a count certificate");
        require(roots_closed(p, endpoint(root.box.real.lo), endpoint(root.box.real.hi)) ==
                    static_cast<int>(root.count),
                "independent Sturm count disagrees with local cluster");
        require(endpoint(root.box.real.hi) - endpoint(root.box.real.lo) < Rational(1) / 8,
                "local cluster is needlessly global");
        require(root.multiplicity.value_or(0) == (root.count == 1 ? 1u : 0u), "unknown cluster multiplicity was guessed");
    }
    const auto gap = power2(-200);
    const auto close = multiply(multiply({-1, 1}, {-1 - gap, 1}), {-4, 1});
    auto pair = certpolysolve_test::solve(input(close), limited.get());
    require(pair.complete() && !pair.goal_met() && pair.roots().size() == 2,
            "tight pair not recovered separately from distant root");
    for (const auto &root : pair.roots())
        require(roots_closed(close, endpoint(root.box.real.lo), endpoint(root.box.real.hi)) ==
                    static_cast<int>(root.count),
                "distinct-root cluster count failed independent check");

    // Every member of this (independent-component) coefficient family must be
    // covered, even though the central double root can split along either axis.
    certpolysolve::Polynomial family;
    family.set(0, certpolysolve::ball(certpolysolve::Number("0"), certpolysolve::Number("1@-38")))
        .set(1, certpolysolve::ball(certpolysolve::Number("0"), certpolysolve::Number("1@-40")))
        .set(2, "-4")
        .set(3, "1");
    auto uniform = certpolysolve_test::solve(family, limited.get());
    if (!uniform.complete() || uniform.goal_met() || uniform.roots().size() != 2)
        throw std::runtime_error("uncertain collision did not yield two local counts: regions=" +
                                 std::to_string(uniform.roots().size()) + "; " + uniform.message());
    uint64_t count = 0;
    for (const auto &root : uniform.roots()) {
        count += root.count;
        require(root.certificate.find("pellet") == 0, "uncertain cluster count not certified locally");
    }
    require(count == 3, "uncertain family coverage lost roots");
}

static void explicit_cluster_goal() {
    const auto options = certpolysolve_test::raw_options();
    require(certpolysolve_test::set_exact_allowance(options.get(), 0) == CPS_OK, "exact allowance refused");
    require(cps_options_set_cluster_radius(options.get(), "1e-8") == CPS_OK, "cluster radius refused");
    const auto p = multiply(multiply({-1, 1}, {-1, 1}), {-4, 1});
    auto result = certpolysolve_test::solve(input(p), options.get());
    require(cps_result_goal(result.get()) == CPS_CLUSTER_RADIUS && result.goal_met() && result.complete(),
            "explicit cluster-radius goal failed");
    require(result.roots().size() == 2 && result.work_bits() < 4096,
            "clustering goal needlessly tried full separation");
    for (const auto &root : result.roots()) {
        const Rational error = (endpoint(root.box.real.hi) - endpoint(root.box.real.lo) +
                                endpoint(root.box.imag.hi) - endpoint(root.box.imag.lo)) /
                               2;
        require(error <= Rational("1/100000000"), "exported cluster radius exceeds goal");
        require(roots_closed(p, endpoint(root.box.real.lo), endpoint(root.box.real.hi)) ==
                    static_cast<int>(root.count),
                "cluster-goal independent count failed");
    }
    // The old approximations() threw when the C getter refused a root; here the getter's code is read.
    bool refused = false;
    for (size_t i = 0; i < result.roots().size(); ++i) {
        double re = 0, im = 0;
        if (cps_result_get_approximation(result.get(), i, &re, &im) != CPS_OK) refused = true;
    }
    require(refused, "cluster goal exposed unknown multiplicities as individual roots");
    auto loose = certpolysolve::solve(certpolysolve::Polynomial{-2, 0, 1}, {.goal = certpolysolve::Goal::clusters(certpolysolve::Number("100"))});
    require(loose.goal_met() && loose.roots().size() == 1 && loose.roots()[0].count == 2,
            "loose global cluster goal failed");
    require(cps_options_set_cluster_radius(options.get(), "0") == CPS_OK, "cluster radius refused");
    require(cps_options_set_precision(options.get(), 128) == CPS_OK && cps_options_set_max_work_bits(options.get(), 128) == CPS_OK,
            "precision refused");
    auto strict = certpolysolve_test::solve(input(p), options.get());
    require(cps_result_goal(strict.get()) == CPS_ROOT_ACCURACY && !strict.goal_met(),
            "zero radius did not restore root-accuracy semantics");
    require(cps_options_set_cluster_radius(options.get(), "1e-8") == CPS_OK, "cluster radius refused");
    require(cps_options_set_real_only(options.get(), 1) == CPS_OK, "real-only refused");
    auto unsupported = certpolysolve_test::solve(input(p), options.get());
    require(unsupported.stop() == certpolysolve::Stop::error && !unsupported.get() && !unsupported.complete(),
            "real-only clustering was not a usage error");
}

int main() {
    certpolysolve::ThreadScope main_lifetime;
    try {
        numerical_families();
        complex_certificates();
        uncertainty_and_limits();
        local_cluster_certificates();
        explicit_cluster_goal();
        concurrent_contexts();
        std::cout << "Independent exact Sturm/Rouche, adversarial and concurrency tests passed\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
