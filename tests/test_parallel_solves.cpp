// W9's first cut through the public C API: the same answer at every worker count (docs/multithreading-design-
// 20261023.md, decision 1). Each input is solved at 1, 3 and 8 workers (cps_options_set_workers), and the answers
// must agree in everything but time: termination, completeness, goal, work bits, square-free verdict, the
// diagnostics text, every region's count, multiplicity, realness, certificate and four exact endpoints, the
// statistics and the engine counters (their times left out). Each solve must meet its goal, and with workers each
// case must reach the split loop it is there for (cps_testing_result_split_jobs): forced piecewise at degree 5,000
// (a Taylor row, the fixed-point transform, the finisher's blocked evaluation), forced secular on 277 roots in four
// clusters 2^-48 apart (the multiprecision packet's Jacobi sweep, which binary64 leaves to it) and forced classical
// at degree 600 (the binary64 engine's evaluations). Two
// more are the input kinds the gate corpus lacks, a native secular equation of 300 terms and a program of degree 600,
// for the same answer alone. Under ThreadSanitizer (the tsan preset) the same solves look for races.
#include <certpolysolve/certpolysolve.h>
#include <certpolysolve/certpolysolve_testing.h>
#include <array>
#include <cstdint>
#include <functional>
#include <initializer_list>
#include <iostream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace {
void check(cps_error e) {
    if (e != CPS_OK) throw std::runtime_error(cps_last_error());
}
std::string text(cps_dyadic_view d) { return std::string(d.mantissa) + "@" + d.exponent; }

// The split sites of certpolysolve_testing.h.
enum Site : unsigned { secular_packet, machine_evaluation, taylor_row, fixed_transform, blocked_evaluation };
using Splits = std::array<uint64_t, CPS_TESTING_SPLIT_SITES>;

struct Outcome {
    std::string answer; // every field but the times
    Splits splits{};
    cps_termination termination{};
};

// The outcome of a result, which is destroyed.
Outcome describe(cps_result *r) {
    Outcome o;
    o.termination = cps_result_termination(r);
    std::ostringstream out;
    out << cps_result_termination(r) << ' ' << cps_result_complete(r) << ' ' << cps_result_goal_met(r) << ' '
        << cps_result_work_bits(r) << ' ' << cps_result_squarefree(r) << '\n'
        << cps_result_diagnostics(r) << '\n';
    for (size_t i = 0; i < cps_result_size(r); ++i) {
        cps_root_view v{};
        check(cps_result_get_root(r, i, &v));
        out << v.count << ' ' << v.multiplicity << ' ' << v.realness << ' ' << (v.certificate ? v.certificate : "")
            << ' ' << text(v.real_lower) << ' ' << text(v.real_upper) << ' ' << text(v.imag_lower) << ' '
            << text(v.imag_upper) << '\n';
    }
    cps_statistics s{};
    s.size = sizeof s;
    check(cps_result_get_statistics(r, &s));
    out << s.work_units << ' ' << s.sweeps << ' ' << s.precision_increases << ' ' << s.final_work_bits << ' '
        << s.highest_work_bits << ' ' << s.regions << ' ' << s.roots << ' ' << s.resolved_regions << ' '
        << s.exact_factors << '\n';
    cps_engine_counters e{};
    e.size = sizeof e;
    check(cps_result_get_engine_counters(r, &e));
    for (const uint64_t c : {e.precision_attempts, e.isolation_passes, e.disk_tests, e.graeffe_steps, e.taylor_shifts,
                             e.taylor_cache_hits, e.boxes_split, e.boxes_excluded, e.newton_attempts,
                             e.newton_successes, e.backend_calls, e.evaluations, e.evaluation_terms,
                             e.root_interactions, e.validation_calls, e.secular_regenerations, e.double_seed_calls,
                             e.candidate_stagnations, e.real_refinement_steps, e.real_bisections,
                             e.pending_components, e.live_cells, e.min_pending_bits, e.max_pending_bits,
                             e.input_uncertainty})
        out << c << ' ';
    o.answer = out.str();
    for (unsigned site = 0; site < CPS_TESTING_SPLIT_SITES; ++site) o.splits[site] = cps_testing_result_split_jobs(r, site);
    cps_result_destroy(r);
    return o;
}

cps_options *options_for(cps_method method, unsigned workers) {
    cps_options *o = nullptr;
    check(cps_options_create(&o));
    check(cps_options_set_method(o, method));
    check(cps_options_set_goal(o, CPS_ROOT_APPROXIMATION));
    check(cps_options_set_workers(o, workers));
    return o;
}

Outcome polynomial_solve(const cps_polynomial *p, cps_method method, unsigned workers) {
    cps_options *o = options_for(method, workers);
    cps_result *r = nullptr;
    check(cps_solve(p, o, &r));
    cps_options_destroy(o);
    return describe(r);
}

Outcome coefficient_solve(const std::vector<std::string> &coefficients, cps_method method, unsigned workers) {
    cps_polynomial *p = nullptr;
    check(cps_polynomial_create(&p));
    for (size_t i = 0; i < coefficients.size(); ++i)
        check(cps_polynomial_set_ball_coefficient(p, i, coefficients[i].c_str(), "0", "0", "0"));
    auto outcome = polynomial_solve(p, method, workers);
    cps_polynomial_destroy(p);
    return outcome;
}

// Real roots in four clusters, as the gate corpus's selector-fresh-clustered family makes them
// (benchmarks/selector/fresh.py, jittered_cluster): root j is (j mod 4) - 2 plus (16 (j + 1) + jitter) / 2^(s + 4)
// with jitter below 16, so consecutive roots of a cluster lie 49 to 79 units of 2^-(s + 4) apart, all of them less
// than degree times 2^-s from the centre. At s = 48 binary64 cannot tell them apart, and the secular method refines
// them in multiprecision packets over every point.
cps_polynomial *clustered(size_t degree, int separation_bits, uint64_t seed) {
    std::mt19937_64 random(seed);
    const int64_t denominator_log2 = separation_bits + 4;
    std::vector<std::string> roots(degree);
    std::vector<const char *> views(degree);
    for (size_t j = 0; j < degree; ++j) {
        const int64_t centre = static_cast<int64_t>(j % 4) - 2, jitter = static_cast<int64_t>(random() % 16);
        const int64_t numerator = centre * (int64_t{1} << denominator_log2) + 16 * static_cast<int64_t>(j + 1) + jitter;
        roots[j] = std::to_string(numerator) + "@-" + std::to_string(denominator_log2);
        views[j] = roots[j].c_str();
    }
    cps_polynomial *p = nullptr;
    check(cps_polynomial_from_roots(views.data(), nullptr, nullptr, degree, &p));
    return p;
}

// The native secular equation sum r_i / (x - i) - 1 over i = 1, ..., n, its residues r_i nonzero integers.
Outcome secular_solve(const std::vector<std::string> &residues, unsigned workers) {
    std::vector<std::string> poles(residues.size());
    std::vector<cps_secular_term_view> terms(residues.size());
    for (size_t i = 0; i < residues.size(); ++i) {
        poles[i] = std::to_string(i + 1);
        terms[i] = {residues[i].c_str(), nullptr, poles[i].c_str(), nullptr};
    }
    cps_secular *s = nullptr;
    check(cps_secular_create(terms.data(), terms.size(), nullptr, &s));
    cps_options *o = options_for(CPS_METHOD_SECULAR, workers);
    cps_result *r = nullptr;
    check(cps_solve_secular(s, o, &r));
    cps_options_destroy(o);
    cps_secular_destroy(s);
    return describe(r);
}

// The program x^600 + 3 x^7 - 1, under Auto.
Outcome program_solve(unsigned workers) {
    cps_program *p = nullptr;
    check(cps_program_create(&p));
    uint64_t x, high, low, three, scaled, sum, one, output;
    check(cps_program_variable(p, &x));
    check(cps_program_power(p, x, 600, &high));
    check(cps_program_power(p, x, 7, &low));
    check(cps_program_constant(p, "3", nullptr, &three));
    check(cps_program_multiply(p, three, low, &scaled));
    check(cps_program_add(p, high, scaled, &sum));
    check(cps_program_constant(p, "1", nullptr, &one));
    check(cps_program_subtract(p, sum, one, &output));
    check(cps_program_set_output(p, output));
    cps_options *o = options_for(CPS_METHOD_AUTO, workers);
    cps_result *r = nullptr;
    check(cps_solve_program(p, o, &r));
    cps_options_destroy(o);
    cps_program_destroy(p);
    return describe(r);
}

std::vector<std::string> random_coefficients(size_t degree, uint64_t seed) {
    std::mt19937_64 random(seed);
    std::vector<std::string> c(degree + 1);
    for (auto &value : c) value = std::to_string(static_cast<long>(random() % 2001) - 1000);
    if (c.back() == "0") c.back() = "1";
    if (c.front() == "0") c.front() = "7";
    return c;
}

std::vector<std::string> random_residues(size_t terms, uint64_t seed) {
    std::mt19937_64 random(seed);
    std::vector<std::string> r(terms);
    for (auto &value : r) {
        const long v = static_cast<long>(random() % 19) - 9;
        value = std::to_string(v ? v : 1);
    }
    return r;
}

// The option: 1 by default (NULL options too), 0 refused with the options unchanged, any other count kept and copied.
void worker_option() {
    cps_options *o = nullptr;
    check(cps_options_create(&o));
    if (cps_options_workers(o) != 1 || cps_options_workers(nullptr) != 1) throw std::runtime_error("the default count");
    if (cps_options_set_workers(o, 0) != CPS_INVALID_ARGUMENT || cps_options_workers(o) != 1)
        throw std::runtime_error("a count of 0 accepted");
    if (cps_options_set_workers(nullptr, 2) != CPS_INVALID_ARGUMENT) throw std::runtime_error("null options accepted");
    check(cps_options_set_workers(o, 5));
    cps_options *copy = nullptr;
    check(cps_options_copy(o, &copy));
    if (cps_options_workers(o) != 5 || cps_options_workers(copy) != 5) throw std::runtime_error("a count kept and copied");
    cps_options_destroy(copy);
    cps_options_destroy(o);
}

// The same outcome at 1, 3 and 8 workers, the goal met; with workers (where the machine has more than one hardware
// thread, which bounds the count) every site of `reached` ran a job on them.
void same_at_every_count(const std::function<Outcome(unsigned)> &solve, std::initializer_list<Site> reached,
                         const std::string &name) {
    const auto one = solve(1);
    if (one.termination != CPS_SUCCESS) throw std::runtime_error(name + ": the goal was not met");
    for (const unsigned workers : {3u, 8u}) {
        const auto other = solve(workers);
        if (other.answer != one.answer)
            throw std::runtime_error(name + ": the answer at " + std::to_string(workers) + " workers differs from one");
        if (std::thread::hardware_concurrency() < 2) continue;
        for (const Site site : reached)
            if (!other.splits[site])
                throw std::runtime_error(name + ": split site " + std::to_string(site) + " ran no job at " +
                                         std::to_string(workers) + " workers");
    }
}
} // namespace

int main() {
    try {
        worker_option();
        const auto piecewise = random_coefficients(5000, 20261023), classical = random_coefficients(600, 20261025);
        const auto residues = random_residues(300, 20261026);
        cps_polynomial *clusters = clustered(277, 48, 20261024);
        same_at_every_count([&](unsigned w) { return coefficient_solve(piecewise, CPS_METHOD_PIECEWISE, w); },
                            {taylor_row, fixed_transform, blocked_evaluation}, "piecewise, degree 5,000");
        same_at_every_count([&](unsigned w) { return polynomial_solve(clusters, CPS_METHOD_SECULAR, w); },
                            {secular_packet}, "secular, 277 clustered roots");
        cps_polynomial_destroy(clusters);
        same_at_every_count([&](unsigned w) { return coefficient_solve(classical, CPS_METHOD_CLASSICAL, w); },
                            {machine_evaluation}, "classical, degree 600");
        same_at_every_count([&](unsigned w) { return secular_solve(residues, w); }, {}, "native secular, 300 terms");
        same_at_every_count(program_solve, {}, "program, degree 600");
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        cps_cleanup_thread();
        return 1;
    }
    cps_cleanup_thread();
    std::cout << "parallel solve tests passed\n";
    return 0;
}
