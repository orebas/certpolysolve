#pragma once
// The fuzz case model shared by certpolysolve_fuzz and the trace probe: seeded
// generation, exact polynomial construction and the replayable case file.
// Generation must stay identical for a given harness_version.
//
// The inputs are the 0.5 C++ header's values (certpolysolve.hpp), which remember their first bad value instead of throwing;
// the harness's own rule for a bad input (it cannot arise from these exact texts) is the exception it always had.
#include "exact_winding_api.hpp"
#include <certpolysolve/certpolysolve.hpp>

#include <algorithm>
#include <array>
#include <cstdio>
#include <fstream>
#include <optional>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fuzz {
using namespace reference;

// Bump when generation or checks change: a seed names a case only together
// with this version. Failure records also carry exact coefficient strings.
// fuzz-3: the inexact family (generate_inexact); the exact cases and their
// checks are those of fuzz-2.
// fuzz-4: one case in four also runs under an iteration budget
// (cps_options_set_max_iterations), in one call and as a resumed session; the
// cases and the other solves are those of fuzz-3.
// fuzz-5 (18 October 2026, the 0.5 API with one-shot solves only): the
// iteration budget runs in one solve, the D1 family's "session" variant is two
// fresh solves (isolation, then the tightened goal), and the axis family's
// raised cap is a fresh solve; generation is fuzz-4's, case for case.
// fuzz-6 (18 October 2026, the 0.5 API's time limit): the soft deadline is the
// solve's time limit (cps_options_set_time_limit) instead of a cancellation
// callback, and a solve it stops reports stop reason TIME_LIMIT; generation is
// fuzz-5's (so fuzz-4's), case for case.
// fuzz-7 (19 October 2026): the program family's solve with a planned cancellation that the soft deadline stops
// first (stop TIME_LIMIT) is a deadline cut (slow_cancelled), not a stop_mismatch, and its checkpoints are not
// compared with the base's; generation is fuzz-6's, case for case, and the other families' checks are unchanged.
// fuzz-7 on the 0.5 inputs (20 October 2026, stage 3b-2): the harness no longer includes the pre-0.5 C++ header. Its
// inputs are the 0.5 header's values and its solves, options and result getters the C API (the poll counter of
// certpolysolve_testing.h in place of a counting cancellation callback). Its records, messages and case files are
// unchanged, so the version stays fuzz-7.
// fuzz-8 (20 October 2026, stage 3b-3): the harness follows the C API of stage 3b-3 (docs/api-spec-0.5.md, C7 to C20).
// Generation is fuzz-7's, case for case (the case files differ in the harness name of their header line only). A
// query that no input of the entry point's kind can answer is a usage error, CPS_INVALID_ARGUMENT from the solve and no
// result (section 5.2), and the harness checks each cell of that table: the solve of a usage-error cell is recorded
// with "usage_error" (term USAGE_ERROR) where fuzz-7 recorded a refusal, a refusal answered or an answer refused as a
// usage error is a failure (usage_error_missing, unexpected_usage_error), and a reversed or inexpandable region is such
// a cell at the solve, not at the setter (C8). The work a solve spent and the highest precision it reached are
// cps_statistics' work_units and highest_work_bits (the pre-0.5 progress block is gone), whose values can differ from the
// old numbers; the public work budget has no zero (0 selects the automatic budget), so a plan of 0 work units runs
// as a one-unit budget and its records say so ("work_planned":0). It also checks, on every result, the goal the result
// reports (C13: goal_default_reported, goal_not_applied), the regions' exact midpoint order (C26: region_order) and the
// result invariants of section 6.2 (C27): I1 to I4 on every termination, including a refusal and a result the soft deadline
// stopped, and I9 (real_region_off_axis); the exception for the planned cancellation that C25 once left
// (known_c25_goal_met_cancelled) is gone. Everything else in the records is fuzz-7's.
// fuzz-9 (23 October 2026): one program-family check changes meaning. The library proves a program's region real by
// reflection when every constant of the program is real (step (1) of "Real roots on program input"),
// so program_real_label fires only on a region labelled real of a program with a non-real constant (the proof needs
// real coefficients); whether a real label is right is the planted roots' check (false_realness) and, on a recurrence,
// Sturm's counts. Generation is fuzz-8's, case for case.
// fuzz-10 (27 October 2026): a new family, --high: the standing case planted at degree 49 to 99 (the degree picker and
// the case's budgets for clusters and zeros), where the small-degree fast paths act (the compensated binary64 tier and
// the double-double cluster repair) and no family reached (every family planted at most 48). Every other family's
// generation, checks and records are fuzz-9's, case for case, and keep its name; only the high family's records and
// case files name fuzz-10.
inline bool high_degree_family = false;
// fuzz-11 (27 October 2026): three new families, --edge53 (generate_edge53, the review's case T08: binary64-exact roots
// beside partners binary64 cannot hold apart, at the goals around binary64's 53 bits), --origin (generate_origin, case
// T12: a tiny root beside a far cloud), --algebraic (generate_algebraic, case T09: irreducible clusters) and --composed
// (generate_composed, case T11: almost structured inputs), the last two judged by the new Factored oracle. Every other
// family's generation, checks and records are fuzz-10's, case for case, and keep their names.
inline bool edge53_family = false, origin_family = false, algebraic_family = false, composed_family = false;
// fuzz-12: the exact-coefficient Auto isolation query is available for complex
// roots. Update its query checks without changing any older generator, and add
// a separately drawn small planted-root isolation family.
inline bool complex_isolation_family = false;
// fuzz-13: a flagged exact complex-isolation preparation refusal may retain
// certified regions; their mathematical checks still run. Generation is fuzz-12's.
// fuzz-14: explicitly enforce I4's coverage flags on that retained-cover path.
inline const char *harness_version() { return "fuzz-14"; }

// ---------------------------------------------------------------- utilities
inline uint64_t splitmix(uint64_t x) {
    x += UINT64_C(0x9e3779b97f4a7c15);
    x = (x ^ (x >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    x = (x ^ (x >> 27)) * UINT64_C(0x94d049bb133111eb);
    return x ^ (x >> 31);
}
inline uint64_t case_seed(uint64_t campaign, uint64_t index) { return splitmix(campaign ^ splitmix(index)); }

class Rng {
    std::mt19937_64 engine_;
  public:
    explicit Rng(uint64_t seed) : engine_(splitmix(seed)) {}
    uint64_t below(uint64_t n) { return engine_() % n; }
    int64_t range(int64_t lo, int64_t hi) { return lo + static_cast<int64_t>(below(static_cast<uint64_t>(hi - lo + 1))); }
    bool percent(unsigned p) { return below(100) < p; }
    template <class T> T weighted(std::initializer_list<std::pair<T, unsigned>> table) {
        unsigned total = 0;
        for (const auto &entry : table) total += entry.second;
        auto x = below(total);
        for (const auto &entry : table) {
            if (x < entry.second) return entry.first;
            x -= entry.second;
        }
        return table.begin()->first;
    }
};

inline std::string json_string(const std::string &s) {
    std::string out = "\"";
    for (unsigned char c : s) {
        switch (c) {
        case '"': out += "\\\""; break;
        case '\\': out += "\\\\"; break;
        case '\n': out += "\\n"; break;
        case '\r': out += "\\r"; break;
        case '\t': out += "\\t"; break;
        default:
            if (c < 0x20 || c >= 0x7f) {
                char buffer[8];
                std::snprintf(buffer, sizeof buffer, "\\u%04x", c);
                out += buffer;
            } else out += static_cast<char>(c);
        }
    }
    return out + "\"";
}

inline const char *method_name(cps_method m) {
    switch (m) {
    case CPS_METHOD_AUTO: return "auto";
    case CPS_METHOD_CLASSICAL: return "classical";
    case CPS_METHOD_SECULAR: return "secular";
    case CPS_METHOD_PIECEWISE: return "piecewise";
    }
    return "?";
}
inline const char *goal_name(cps_goal g) {
    switch (g) {
    case CPS_ROOT_ACCURACY: return "accuracy";
    case CPS_CLUSTER_RADIUS: return "cluster";
    case CPS_ROOT_APPROXIMATION: return "approximation";
    case CPS_ROOT_ISOLATION: return "isolation";
    case CPS_GOAL_DEFAULT: return "default"; // no case draws it: the options' initial goal
    }
    return "?";
}
inline const char *termination_name(cps_termination t) {
    switch (t) {
    case CPS_SUCCESS: return "SUCCESS";
    case CPS_PRECISION_LIMIT: return "PRECISION_LIMIT";
    case CPS_CANCELLED: return "CANCELLED";
    case CPS_ZERO_POLYNOMIAL: return "ZERO_POLYNOMIAL";
    case CPS_INDETERMINATE_DEGREE: return "INDETERMINATE_DEGREE";
    case CPS_UNSUPPORTED_QUERY: return "UNSUPPORTED_QUERY";
    case CPS_WORK_LIMIT: return "WORK_LIMIT";
    }
    return "?";
}
inline const char *stop_name(cps_stop_reason s) {
    switch (s) {
    case CPS_STOP_NONE: return "NONE";
    case CPS_STOP_WORK_BUDGET: return "WORK_BUDGET";
    case CPS_STOP_CELL_BUDGET: return "CELL_BUDGET";
    case CPS_STOP_PRECISION: return "PRECISION";
    case CPS_STOP_CANCELLED: return "CANCELLED";
    case CPS_STOP_INPUT: return "INPUT";
    case CPS_STOP_ITERATION_BUDGET: return "ITERATION_BUDGET";
    case CPS_STOP_TIME_LIMIT: return "TIME_LIMIT";
    }
    return "?";
}
inline std::optional<cps_method> parse_method(const std::string &s) {
    for (auto m : {CPS_METHOD_AUTO, CPS_METHOD_CLASSICAL, CPS_METHOD_SECULAR, CPS_METHOD_PIECEWISE})
        if (s == method_name(m)) return m;
    return std::nullopt;
}
inline std::optional<cps_goal> parse_goal(const std::string &s) {
    for (auto g : {CPS_ROOT_ACCURACY, CPS_CLUSTER_RADIUS, CPS_ROOT_APPROXIMATION, CPS_ROOT_ISOLATION})
        if (s == goal_name(g)) return g;
    return std::nullopt;
}
inline Rational parse_rational(const std::string &text) {
    const auto slash = text.find('/');
    if (slash == std::string::npos) return Rational(Integer(text));
    return Rational(Integer(text.substr(0, slash))) / Rational(Integer(text.substr(slash + 1)));
}
inline Complex conjugate(const Complex &z) { return {z.real, -z.imag}; }
inline std::string complex_text(const Complex &z) { return z.real.str() + " " + z.imag.str(); }

// ------------------------------------------------------------------ the case
// Expansion: the program family's recurrences, whose roots are not planted: the counts are decided on the
// program's exact expansion (expand_program). Factored (harness fuzz-11, the algebraic family): the polynomial is a
// product of exact factors with multiplicities (Spec::factors), distinct and pairwise coprime, so their product is
// the square-free part; the counts and multiplicities are decided by exact Sturm-Cauchy counts on each factor.
enum class Oracle { Planted, Mignotte, Expansion, Factored };
inline const char *oracle_name(Oracle o) {
    return o == Oracle::Planted ? "planted" : o == Oracle::Mignotte ? "mignotte" : o == Oracle::Expansion ? "expansion"
                                                                                                        : "factored";
}
struct Spec {
    uint64_t seed = 0;
    std::string family;
    Oracle oracle = Oracle::Planted;
    std::vector<Complex> roots; // planted roots, repeated by multiplicity
    // Oracle::Factored: the exact factors (coefficients lowest degree first) and their multiplicities.
    std::vector<std::pair<ComplexPolynomial, int>> factors;
    Complex multiplier{1};
    int mignotte_n = 0, mignotte_k = 0;
    cps_method method = CPS_METHOD_AUTO;
    cps_goal goal = CPS_ROOT_APPROXIMATION;
    uint64_t bits = 53, max_work_bits = 4096, work_units = 100000, live_cells = 65536;
    bool real_only = false, monotone = true, determinism = false;
    int cluster_exponent = 20; // CPS_CLUSTER_RADIUS radius 2^-cluster_exponent
    // Inexact input (the --inexact family): 0 for exact coefficients. Otherwise
    // the input's coefficient balls are the smallest boxes that hold two
    // polynomials with planted roots: the one with `roots`, and its partner,
    // whose root `partner[i]` is `roots[i]` moved by a few 2^-inexact of its
    // own modulus. Every certified statement must hold for both.
    int inexact = 0;
    std::vector<Complex> partner;
    // The D1 family (--d1, generate_d1): 0 for none, else a D1Variant; the region's ends in the coefficient
    // grammar ("" for an unbounded side); and, for the refusal variant, the refusal probed (D1Refusal).
    int d1 = 0, d1_refusal = 0;
    std::string region_lower, region_upper, region_imag_lower, region_imag_upper;
    // The axis family (--axis, generate_axis): 0 for none, else 1; the largest exponent k of its pairs' 2^-k
    // (the resolution session raises the cap past it); its transform (0 none, 1 a rotation, 2 a non-real
    // scalar on the real polynomial).
    int axis = 0, axis_k = 0, axis_transform = 0;
    // The secular family (--secular, generate_secular): 0 for none, else 1; the poles b_k, one a planted root
    // (the residues follow: secular_residues); for the refusal variant, the refusal probed (SecularRefusal).
    int secular = 0, secular_refusal = 0;
    std::vector<Complex> poles;
    // The program family (--program, generate_program): 0 for none, else its shape (ProgramShape); the
    // program's degree; its text after the "PSPROGRAM 1" header, one operation a line (x, c, add, sub, mul,
    // neg, pow, out). A planted shape keeps its roots in `roots` (the program computes a nonzero multiple of
    // their product); a recurrence has none (Oracle::Expansion). The stops, drawn from the family's stream:
    // a solve cancelled at checkpoint 1 + program_cancel mod the base solve's count of checkpoints (0: no
    // such solve); the cap-then-cancellation solve, at the cap program_cap_cancel, cancelled at checkpoint
    // 1 + program_cap_cancel_draw mod that cap's count (0: none); a solve under a time limit of
    // program_time_limit_us microseconds (0: none).
    int program = 0, program_degree = 0;
    std::vector<std::string> program_lines;
    uint64_t program_cancel = 0, program_cap_cancel = 0, program_cap_cancel_draw = 0, program_time_limit_us = 0;
    bool region_set() const {
        return !region_lower.empty() || !region_upper.empty() || !region_imag_lower.empty() || !region_imag_upper.empty();
    }
};
// D1 (docs/real-interval-contract-note-20261015.md): the isolation goal and the closed region of a real-only query.
enum class D1Variant : int { None = 0, Isolation = 1, Region = 2, Session = 3, Refusal = 4 };
enum class D1Refusal : int {
    None = 0,
    RegionWithoutRealOnly = 1,
    IsolationWithoutRealOnly = 2,
    ExplicitMethod = 3,
    ClusterGoal = 4,
    NonRealInput = 5,
    InexactInput = 6,
    InexpandableEnd = 7,
    SetterReversed = 8,
    SetterMalformed = 9
};
// The secular family's refusals: what the constructor (two equal poles, a zero residue) or the solve (the empty
// equation answered as anything but complete and empty; an explicit classical or piecewise method; real_only on
// a numerator that is not real) must refuse.
enum class SecularRefusal : int { None = 0, DuplicatePole = 1, ZeroResidue = 2, Empty = 3, ExplicitMethod = 4, RealOnlyNonReal = 5 };

inline ComplexPolynomial build(const Spec &s) {
    ComplexPolynomial p;
    if (s.oracle == Oracle::Factored) {
        p = {Complex(1)};
        for (const auto &[factor, multiplicity] : s.factors)
            for (int m = 0; m < multiplicity; ++m) p = multiply(p, factor);
    } else if (s.oracle == Oracle::Mignotte) {
        const int n = s.mignotte_n;
        const Rational a = power2(s.mignotte_k);
        p.assign(static_cast<size_t>(n) + 1, Complex(0));
        p[0] = Complex(Rational(-2));
        p[1] = Complex(Rational(4) * a);
        p[2] = Complex(Rational(-2) * a * a);
        p[static_cast<size_t>(n)] = p[static_cast<size_t>(n)] + Complex(Rational(1));
    } else {
        // Clear each root's denominators: (d x - d z) has Gaussian-integer
        // coefficients, so the product needs no rational normalization.
        using GaussianInteger = std::pair<Integer, Integer>;
        std::vector<GaussianInteger> q{{Integer(1), Integer(0)}};
        for (const auto &z : s.roots) {
            const Integer dr = denominator(z.real), di = denominator(z.imag);
            const Integer d = boost::multiprecision::lcm(dr, di);
            const Integer a = -numerator(z.real) * (d / dr), b = -numerator(z.imag) * (d / di);
            std::vector<GaussianInteger> next(q.size() + 1, {Integer(0), Integer(0)});
            for (size_t j = 0; j < q.size(); ++j) {
                next[j].first += q[j].first * a - q[j].second * b;
                next[j].second += q[j].first * b + q[j].second * a;
                next[j + 1].first += q[j].first * d;
                next[j + 1].second += q[j].second * d;
            }
            q = std::move(next);
        }
        for (auto &c : q) p.emplace_back(Rational(c.first), Rational(c.second));
    }
    for (auto &c : p) c = c * s.multiplier;
    return p;
}
inline bool real_coefficients(const ComplexPolynomial &p) {
    return std::all_of(p.begin(), p.end(), [](const Complex &c) { return c.imag == 0; });
}
// The exact input polynomial: each coefficient's exact text, never a double. The pre-0.5 setter threw on a bad value;
// the 0.5 polynomial remembers the first one, so the refusal is raised here, where it always was.
inline certpolysolve::Polynomial accepted(certpolysolve::Polynomial input) {
    if (!input.ok()) throw std::runtime_error(input.status());
    return input;
}
inline certpolysolve::Polynomial to_input(const ComplexPolynomial &p) { return accepted(reference::input(p)); }
// The two planted polynomials of an inexact case on one scale: monic products
// times the multiplier (build() clears each root's denominators, which would
// scale the two differently).
inline ComplexPolynomial build_monic(const std::vector<Complex> &roots, const Complex &multiplier) {
    ComplexPolynomial p{Complex(1)};
    for (const auto &z : roots) {
        ComplexPolynomial next(p.size() + 1, Complex(0));
        for (size_t j = 0; j < p.size(); ++j) {
            next[j] = next[j] - p[j] * z;
            next[j + 1] = next[j + 1] + p[j];
        }
        p = std::move(next);
    }
    for (auto &c : p) c = c * multiplier;
    return p;
}
// Coefficient j is the box with corners p[j] and q[j]: midpoint (p + q)/2,
// radii |p - q|/2 in the real and the imaginary part. A coefficient the two
// share stays exact.
inline certpolysolve::Polynomial to_input(const ComplexPolynomial &p, const ComplexPolynomial &q) {
    require(p.size() == q.size(), "the planted polynomials of an inexact case differ in degree");
    certpolysolve::Polynomial input;
    const auto half_distance = [](const Rational &a, const Rational &b) { return (a < b ? b - a : a - b) / 2; };
    for (size_t j = 0; j < p.size(); ++j) {
        const Rational real = (p[j].real + q[j].real) / 2, imag = (p[j].imag + q[j].imag) / 2;
        const Rational real_radius = half_distance(p[j].real, q[j].real),
                       imag_radius = half_distance(p[j].imag, q[j].imag);
        const auto centre = certpolysolve::Number::complex(certpolysolve::Number(real.str()), certpolysolve::Number(imag.str()));
        if (real_radius == 0 && imag_radius == 0) input.set(j, centre);
        else
            input.set(j, certpolysolve::ball(centre, certpolysolve::Number(real_radius.str()), certpolysolve::Number(imag_radius.str())));
    }
    return accepted(std::move(input));
}

// A program's text (the PSPROGRAM 1 grammar of cps_program_parse) and its exact expansion, computed by the harness
// from the same lines (exact Gaussian rationals, no FLINT): the program family's oracle polynomial.
inline std::string program_text(const Spec &s) {
    std::string text = "PSPROGRAM 1\n";
    for (const auto &line : s.program_lines) text += line + "\n";
    return text;
}
inline ComplexPolynomial trimmed(ComplexPolynomial p) {
    while (!p.empty() && zero(p.back())) p.pop_back();
    return p;
}
inline ComplexPolynomial add(const ComplexPolynomial &p, const ComplexPolynomial &q, bool subtract = false) {
    ComplexPolynomial sum(std::max(p.size(), q.size()), Complex(0));
    for (size_t j = 0; j < p.size(); ++j) sum[j] = sum[j] + p[j];
    for (size_t j = 0; j < q.size(); ++j) sum[j] = subtract ? sum[j] - q[j] : sum[j] + q[j];
    return trimmed(std::move(sum));
}
inline ComplexPolynomial expand_program(const std::vector<std::string> &lines) {
    std::vector<ComplexPolynomial> value;
    std::optional<size_t> output;
    const auto operand = [&](const std::string &text) {
        const auto index = static_cast<size_t>(std::stoull(text));
        require(index < value.size(), "a program operand is not an earlier node");
        return index;
    };
    for (const auto &line : lines) {
        std::istringstream fields(line);
        std::string op, a, b;
        fields >> op >> a >> b;
        if (op == "x") value.push_back({Complex(0), Complex(1)});
        else if (op == "c") value.push_back(trimmed({Complex(parse_rational(a), b.empty() ? Rational(0) : parse_rational(b))}));
        else if (op == "add" || op == "sub") value.push_back(add(value[operand(a)], value[operand(b)], op == "sub"));
        else if (op == "mul") value.push_back(multiply(value[operand(a)], value[operand(b)]));
        else if (op == "neg") value.push_back(add({}, value[operand(a)], true));
        else if (op == "pow") {
            const ComplexPolynomial base = value[operand(a)];
            ComplexPolynomial power{Complex(1)};
            for (uint64_t e = std::stoull(b); e > 0; --e) power = multiply(power, base);
            value.push_back(std::move(power));
        } else if (op == "out") output = operand(a);
        else throw std::runtime_error("unknown program operation: " + op);
    }
    require(output.has_value(), "the program has no output");
    return value[*output];
}
// Whether a constant of the program is non-real ("c re im" with im not 0). The library proves a region real by
// reflection only on a program whose constants are all real, whose coefficients are then real.
inline bool program_has_nonreal_constant(const std::vector<std::string> &lines) {
    for (const auto &line : lines) {
        std::istringstream fields(line);
        std::string op, a, b;
        fields >> op >> a >> b;
        if (op == "c" && !b.empty() && parse_rational(b) != 0) return true;
    }
    return false;
}

inline std::string describe(const Spec &s, bool coefficients = true) {
    std::ostringstream out;
    out << "# certpolysolve fuzz case (harness " << harness_version() << ")\n"
        << "seed " << s.seed << "\nfamily " << s.family << "\noracle "
        << oracle_name(s.oracle) << "\nmethod " << method_name(s.method)
        << "\ngoal " << goal_name(s.goal) << "\nbits " << s.bits << "\nmax_work_bits " << s.max_work_bits
        << "\nwork_units " << s.work_units << "\nlive_cells " << s.live_cells << "\nreal_only " << s.real_only
        << "\ncluster_exponent " << s.cluster_exponent << "\nmonotone " << s.monotone << "\ndeterminism "
        << s.determinism << "\nmultiplier " << complex_text(s.multiplier) << '\n';
    if (s.oracle == Oracle::Mignotte) out << "mignotte " << s.mignotte_n << ' ' << s.mignotte_k << '\n';
    for (const auto &z : s.roots) out << "root " << complex_text(z) << '\n';
    for (const auto &[factor, multiplicity] : s.factors) {
        out << "factor " << multiplicity;
        for (const auto &c : factor) out << ' ' << complex_text(c);
        out << '\n';
    }
    if (s.inexact) {
        out << "inexact " << s.inexact << '\n';
        for (const auto &z : s.partner) out << "partner " << complex_text(z) << '\n';
    }
    // Every region outside the program family belongs to a D1 query (draw_d1_query sets s.d1), so those case files
    // are unchanged; a program case's region is its refusal probe.
    if (s.d1) out << "d1 " << s.d1 << ' ' << s.d1_refusal << '\n';
    if (s.d1 || s.region_set()) {
        const auto end = [](const std::string &text) { return text.empty() ? std::string("-") : text; };
        out << "region " << end(s.region_lower) << ' ' << end(s.region_upper) << ' ' << end(s.region_imag_lower) << ' '
            << end(s.region_imag_upper) << '\n';
    }
    if (s.axis) out << "axis " << s.axis << ' ' << s.axis_k << ' ' << s.axis_transform << '\n';
    if (s.secular) {
        out << "secular " << s.secular << ' ' << s.secular_refusal << '\n';
        for (const auto &b : s.poles) out << "pole " << complex_text(b) << '\n';
    }
    if (s.program) {
        out << "program " << s.program << ' ' << s.program_degree << '\n'
            << "program_stops " << s.program_cancel << ' ' << s.program_cap_cancel << ' ' << s.program_cap_cancel_draw
            << ' ' << s.program_time_limit_us << '\n';
        for (const auto &line : s.program_lines) out << "op " << line << '\n';
    }
    if (coefficients && s.program) {
        const auto p = expand_program(s.program_lines);
        out << "# the program's exact expansion, ascending (its coefficient twin)\n";
        for (size_t j = 0; j < p.size(); ++j) out << "coefficient " << j << ' ' << complex_text(p[j]) << '\n';
    } else if (coefficients && s.inexact) {
        const auto p = build_monic(s.roots, s.multiplier), q = build_monic(s.partner, s.multiplier);
        out << "# ascending coefficient boxes: the two corners (inputs to cps_polynomial_set_ball_coefficient\n"
            << "# are their midpoint and half their distance in each part)\n";
        for (size_t j = 0; j < p.size(); ++j)
            out << "coefficient " << j << ' ' << complex_text(p[j]) << " | " << complex_text(q[j]) << '\n';
    } else if (coefficients) {
        const auto p = build(s);
        out << "# ascending exact coefficients (input to cps_polynomial_set_coefficient)\n";
        for (size_t j = 0; j < p.size(); ++j) out << "coefficient " << j << ' ' << complex_text(p[j]) << '\n';
    }
    return out.str();
}

inline Spec load_spec(const std::string &path) {
    std::ifstream in(path);
    require(static_cast<bool>(in), "cannot open case file");
    Spec s;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream fields(line);
        std::string key;
        fields >> key;
        std::string a, b;
        if (key == "seed") fields >> s.seed;
        else if (key == "family") { std::getline(fields, a); s.family = a.empty() ? a : a.substr(1); }
        else if (key == "oracle") {
            fields >> a;
            s.oracle = a == "mignotte" ? Oracle::Mignotte : a == "expansion" ? Oracle::Expansion
                       : a == "factored" ? Oracle::Factored : Oracle::Planted;
        }
        else if (key == "method") { fields >> a; s.method = parse_method(a).value(); }
        else if (key == "goal") { fields >> a; s.goal = parse_goal(a).value(); }
        else if (key == "bits") fields >> s.bits;
        else if (key == "max_work_bits") fields >> s.max_work_bits;
        else if (key == "work_units") fields >> s.work_units;
        else if (key == "live_cells") fields >> s.live_cells;
        else if (key == "real_only") fields >> s.real_only;
        else if (key == "cluster_exponent") fields >> s.cluster_exponent;
        else if (key == "monotone") fields >> s.monotone;
        else if (key == "determinism") fields >> s.determinism;
        else if (key == "multiplier") { fields >> a >> b; s.multiplier = {parse_rational(a), parse_rational(b)}; }
        else if (key == "mignotte") fields >> s.mignotte_n >> s.mignotte_k;
        else if (key == "root") { fields >> a >> b; s.roots.emplace_back(parse_rational(a), parse_rational(b)); }
        else if (key == "factor") {
            int multiplicity = 0;
            fields >> multiplicity;
            ComplexPolynomial factor;
            while (fields >> a >> b) factor.emplace_back(parse_rational(a), parse_rational(b));
            s.factors.emplace_back(std::move(factor), multiplicity);
        }
        else if (key == "inexact") fields >> s.inexact;
        else if (key == "partner") { fields >> a >> b; s.partner.emplace_back(parse_rational(a), parse_rational(b)); }
        else if (key == "d1") fields >> s.d1 >> s.d1_refusal;
        else if (key == "axis") fields >> s.axis >> s.axis_k >> s.axis_transform;
        else if (key == "secular") fields >> s.secular >> s.secular_refusal;
        else if (key == "pole") { fields >> a >> b; s.poles.emplace_back(parse_rational(a), parse_rational(b)); }
        else if (key == "program") fields >> s.program >> s.program_degree;
        else if (key == "program_stops")
            fields >> s.program_cancel >> s.program_cap_cancel >> s.program_cap_cancel_draw >> s.program_time_limit_us;
        else if (key == "op") { std::getline(fields, a); s.program_lines.push_back(a.empty() ? a : a.substr(1)); }
        else if (key == "region") {
            std::string c, d;
            fields >> a >> b >> c >> d;
            const auto end = [](const std::string &text) { return text == "-" ? std::string() : text; };
            s.region_lower = end(a);
            s.region_upper = end(b);
            s.region_imag_lower = end(c);
            s.region_imag_upper = end(d);
        }
        else if (key == "coefficient") continue; // derived from the roots; informational
        else throw std::runtime_error("unknown case-file key: " + key);
    }
    return s;
}

// ---------------------------------------------------------------- generation
inline Rational coordinate(Rng &r, int exponent) {
    // A random dyadic or non-dyadic rational of magnitude about 2^exponent.
    if (r.below(4) == 0) {
        const int64_t q = 1 + static_cast<int64_t>(r.below(999));
        const int64_t p = r.range(-2 * q, 2 * q);
        return Rational(Integer(p)) / Rational(Integer(q)) * power2(exponent);
    }
    const int bits = 1 + static_cast<int>(r.below(24));
    const int64_t m = r.range(-(int64_t(1) << bits), int64_t(1) << bits);
    return Rational(Integer(m)) * power2(exponent - bits);
}

struct Planter {
    Rng &r;
    bool symmetric; // keep the multiset closed under conjugation
    std::vector<Complex> roots;
    int size() const { return static_cast<int>(roots.size()); }
    int cost(const Complex &z, int copies = 1) const { return copies * (symmetric && z.imag != 0 ? 2 : 1); }
    void add(const Complex &z, int copies = 1) {
        for (int c = 0; c < copies; ++c) roots.push_back(z);
        if (symmetric && z.imag != 0)
            for (int c = 0; c < copies; ++c) roots.push_back(conjugate(z));
    }
};

// Primitive Pythagorean triples: (a + b i)/c lies on the unit circle.
inline constexpr std::array<std::array<int, 3>, 12> triples{{{3, 4, 5}, {5, 12, 13}, {8, 15, 17}, {7, 24, 25},
                                                      {20, 21, 29}, {12, 35, 37}, {9, 40, 41}, {28, 45, 53},
                                                      {11, 60, 61}, {33, 56, 65}, {16, 63, 65}, {48, 55, 73}}};

// The degree a case plants at most: the standing families' 48 (fuzz-1 to fuzz-9), and the high family's 99, the top of
// the direction set on 27 October 2026 (easy and near-easy inputs below degree 100).
inline constexpr int standing_max_degree = 48;
inline constexpr int high_max_degree = 99;
inline int max_planted_degree() { return high_degree_family ? high_max_degree : standing_max_degree; }
inline int pick_degree(Rng &r) {
    if (high_degree_family) return static_cast<int>(r.range(standing_max_degree + 1, high_max_degree));
    const auto band = r.below(10);
    if (band < 4) return static_cast<int>(r.range(1, 8));
    if (band < 7) return static_cast<int>(r.range(9, 20));
    return static_cast<int>(r.range(21, standing_max_degree));
}

inline Spec generate(uint64_t seed) {
    Rng r(seed);
    const int cap = max_planted_degree();
    Spec s;
    s.seed = seed;
    const int n = pick_degree(r);
    std::ostringstream family;
    const auto kind = r.below(100);
    if (kind < 8) {
        s.oracle = Oracle::Mignotte;
        s.mignotte_n = std::max(3, n);
        s.mignotte_k = static_cast<int>(r.range(1, 24));
        family << "mignotte(n=" << s.mignotte_n << ",k=" << s.mignotte_k << ")";
        if (r.percent(30)) {
            const int k = static_cast<int>(r.range(-60, 60));
            s.multiplier = Complex(power2(k));
            family << "*2^" << k;
        }
    } else {
        // Reserve room for the modifiers, then fill the base set.
        const bool cluster = r.percent(35), repeats = r.percent(30), zeros = r.percent(20);
        int reserved = 0;
        const int cluster_size = cluster ? static_cast<int>(r.range(2, 8)) : 0;
        const int zero_count = zeros ? static_cast<int>(r.range(1, 5)) : 0;
        reserved = std::min(n - 1, cluster_size + zero_count);
        const int base_size = std::max(1, n - reserved);
        enum Base { Uniform, Real, Conjugate, Circle, Wilkinson, Grid, Mixed };
        const Base base = r.weighted<Base>({{Uniform, 24}, {Real, 16}, {Conjugate, 16}, {Circle, 10},
                                            {Wilkinson, 8}, {Grid, 10}, {Mixed, 16}});
        const bool symmetric = base == Real || base == Conjugate || base == Wilkinson ||
                               ((base == Circle || base == Grid || base == Mixed) && r.percent(50));
        Planter p{r, symmetric, {}};
        switch (base) {
        case Uniform:
            family << "uniform";
            while (p.size() < base_size) p.add({coordinate(r, 0), coordinate(r, 0)});
            break;
        case Real:
            family << "real";
            while (p.size() < base_size) p.add({coordinate(r, 0)});
            break;
        case Conjugate:
            family << "conjugate";
            while (p.size() < base_size) {
                Complex z{coordinate(r, 0), r.percent(70) ? coordinate(r, 0) : Rational(0)};
                if (p.size() + p.cost(z) > base_size) z.imag = 0;
                p.add(z);
            }
            break;
        case Circle: {
            const int radius = r.percent(30) ? static_cast<int>(r.range(-20, 20)) : 0;
            family << "circle(r=2^" << radius << ")";
            while (p.size() < base_size) {
                const auto &t = triples[r.below(triples.size())];
                Rational a(t[0], t[2]), b(t[1], t[2]);
                if (r.percent(50)) std::swap(a, b);
                if (r.percent(50)) a = -a;
                if (r.percent(50)) b = -b;
                Complex z{a * power2(radius), b * power2(radius)};
                if (p.size() + p.cost(z) > base_size) z = Complex(r.percent(50) ? power2(radius) : Rational(-power2(radius)));
                p.add(z);
            }
            break;
        }
        case Wilkinson: {
            const int shift = r.percent(50) ? static_cast<int>(r.range(-8, 0)) : 0;
            family << "wilkinson(2^" << shift << ")";
            for (int k = 1; p.size() < base_size; ++k) p.add(Complex(Rational(k) * power2(shift)));
            break;
        }
        case Grid:
            family << "grid";
            while (p.size() < base_size) {
                Complex z{Rational(r.range(-3, 3)), Rational(r.range(-3, 3))};
                if (p.size() + p.cost(z) > base_size) z.imag = 0;
                p.add(z);
            }
            break;
        case Mixed:
            family << "mixed";
            while (p.size() < base_size) {
                const int e = static_cast<int>(r.range(-60, 60));
                Complex z{coordinate(r, e), r.percent(75) ? coordinate(r, e) : Rational(0)};
                if (p.size() + p.cost(z) > base_size) z.imag = 0;
                p.add(z);
            }
            break;
        }
        if (cluster && p.size() + 2 <= cap) {
            // Distinct small Gaussian-integer offsets at scale 2^-depth about a
            // center, optionally with a nested sub-cluster.
            const int depth = static_cast<int>(r.range(10, 120));
            // Braced initializers fix the evaluation order of the random draws,
            // so GCC and Clang builds generate the same case from a seed.
            Complex center = r.percent(50) && !p.roots.empty() ? p.roots[r.below(p.roots.size())]
                                                               : Complex{coordinate(r, 0), p.symmetric ? Rational(0) : coordinate(r, 0)};
            if (p.symmetric && r.percent(60)) center.imag = 0;
            const bool real_offsets = p.symmetric && center.imag == 0;
            std::set<std::pair<int, int>> used;
            int added = 0;
            const int budget = std::min(cluster_size, cap - p.size());
            family << "+cluster(k=" << budget << ",2^-" << depth << ")";
            while (added < budget && used.size() < (real_offsets ? 8u : 24u)) {
                std::pair<int, int> d{static_cast<int>(r.range(-4, 4)), real_offsets ? 0 : static_cast<int>(r.range(-2, 2))};
                if (!used.insert(d).second) continue;
                const Complex z{center.real + d.first * power2(-depth), center.imag + d.second * power2(-depth)};
                if (p.size() + p.cost(z) > cap || added + p.cost(z) > budget + 1) continue;
                p.add(z);
                added += p.cost(z);
            }
            if (r.percent(20) && p.size() + 2 <= cap) {
                const int deeper = depth + static_cast<int>(r.range(10, 60));
                family << "+nested(2^-" << deeper << ")";
                const Complex z = p.roots.back();
                p.add({z.real + power2(-deeper), z.imag});
                if (p.size() + 1 <= cap && !p.symmetric) p.add({z.real, z.imag + power2(-deeper)});
            }
        }
        if (repeats && !p.roots.empty()) {
            const int picks = static_cast<int>(r.range(1, 3));
            family << "+repeat";
            for (int j = 0; j < picks; ++j) {
                const Complex z = p.roots[r.below(p.roots.size())];
                const int extra = static_cast<int>(r.range(1, 5)); // multiplicity 2..6
                const int fit = std::min(extra, (cap - p.size()) / p.cost(z));
                if (fit > 0) {
                    p.add(z, fit);
                    family << "(m+" << fit << ")";
                }
            }
        }
        if (zeros && p.size() < cap) {
            const int k = std::min(zero_count, cap - p.size());
            p.add(Complex(0), k);
            family << "+x^" << k;
        }
        // Transformations of the whole root set.
        Complex rotation{1};
        if (r.percent(15)) { rotation = {0, 1}; family << "*i"; }
        else if (r.percent(12)) { rotation = {Rational(3, 5), Rational(4, 5)}; family << "*(3+4i)/5"; }
        if (r.percent(30)) {
            const int k = static_cast<int>(r.range(-60, 60));
            rotation = rotation * Complex(power2(k));
            family << "*2^" << k;
        }
        Complex shift{0};
        if (r.percent(8)) {
            shift = {coordinate(r, 0), p.symmetric && rotation.imag == 0 ? Rational(0) : coordinate(r, 0)};
            family << "+shift";
        }
        for (auto &z : p.roots) z = z * rotation + shift;
        s.roots = std::move(p.roots);
        const bool real_polynomial_possible = p.symmetric && rotation.imag == 0 && shift.imag == 0;
        if (r.percent(40)) {
            const auto m = r.below(3);
            if (m == 0) {
                const int k = static_cast<int>(r.range(-100, 100));
                s.multiplier = Complex(power2(k));
                family << ";c=2^" << k;
            } else if (m == 1 || real_polynomial_possible) {
                const int64_t a = r.range(1, 999), b = r.range(1, 999);
                s.multiplier = Complex(Rational(Integer(r.percent(50) ? a : -a)) / Rational(Integer(b)));
                family << ";c=rational";
            } else {
                s.multiplier = {Rational(r.range(-9, 9)) / 7, Rational(r.range(1, 9)) / 11};
                family << ";c=gaussian";
            }
        }
    }
    s.family = family.str();
    const bool real_polynomial = real_coefficients(build(s));
    s.method = static_cast<cps_method>(r.below(4));
    const bool explicit_method = s.method != CPS_METHOD_AUTO;
    const auto g = r.below(100);
    s.goal = explicit_method ? (g < 3 ? CPS_CLUSTER_RADIUS : g < 51 ? CPS_ROOT_ACCURACY : CPS_ROOT_APPROXIMATION)
                             : (g < 25 ? CPS_CLUSTER_RADIUS : g < 60 ? CPS_ROOT_ACCURACY : CPS_ROOT_APPROXIMATION);
    s.real_only = real_polynomial ? r.percent(explicit_method ? 5 : 30) : r.percent(3);
    if (s.real_only && s.goal == CPS_CLUSTER_RADIUS && !r.percent(10))
        s.goal = r.percent(50) ? CPS_ROOT_ACCURACY : CPS_ROOT_APPROXIMATION;
    s.cluster_exponent = static_cast<int>(r.range(1, 80));
    s.bits = r.weighted<uint64_t>({{16, 1}, {53, 1}, {100, 1}});
    s.max_work_bits = r.weighted<uint64_t>({{32, 2}, {48, 3}, {64, 14}, {72, 10}, {80, 10}, {96, 10}, {128, 12},
                                            {192, 8}, {256, 8}, {512, 8}, {4096, 15}});
    s.work_units = r.weighted<uint64_t>({{1, 8}, {10, 12}, {100, 18}, {1000, 22}, {100000, 40}});
    s.live_cells = r.weighted<uint64_t>({{16, 15}, {256, 20}, {65536, 65}});
    s.determinism = r.percent(5);
    return s;
}

// Separate root stream; older families' seeds and draws remain unchanged.
// Degree at most eight keeps every case within the small exact-count regime.
inline constexpr uint64_t complex_isolation_stream = UINT64_C(0x63706c782d69736f); // "cplx-iso"
inline Spec generate_complex_isolation(uint64_t seed) {
    Spec s = generate(seed); // reuse the standing resource/budget distribution
    Rng r(seed ^ complex_isolation_stream);
    s.oracle = Oracle::Planted;
    s.mignotte_n = s.mignotte_k = 0;
    s.roots.clear();
    const bool real_roots = r.percent(50);
    const auto degree = static_cast<size_t>(r.range(1, 8));
    for (size_t j = 0; j < degree; ++j) {
        Complex point{coordinate(r, 0), real_roots ? Rational(0) : coordinate(r, 0)};
        if (r.percent(15)) point = Complex(0);
        else if (j && r.percent(30)) point = s.roots[r.below(j)];
        else if (j && r.percent(30)) {
            point = s.roots[r.below(j)];
            point.real += power2(-static_cast<int>(r.range(8, 96)));
        }
        s.roots.push_back(std::move(point));
    }
    s.method = CPS_METHOD_AUTO;
    s.real_only = false;
    s.goal = CPS_ROOT_ISOLATION;
    s.family = "complex-isolation(n=" + std::to_string(degree) + ")";
    return s;
}

// The edge53 family (--edge53, harness fuzz-11): the review's case T08 (codex-notes/review-20261014/adversarial-cases.md;
// dev/docs/adversarial-corpus-note-20261016.md). A binary64-exact root beside a partner that binary64 cannot hold
// apart, solved at the goals around binary64's 53 bits: a certificate written at a rounded point of the handover from
// binary64 to multiprecision (the review's H01, H02: the binary64 packet rounds unmoved points back above the default
// goal) would show here as a wrong or inaccurate answer, and lost digits in the records' work_bits and sweeps. From
// the family's own stream, the standing case of the seed is not used:
// - one anchor (70%) or two: an odd 53-bit integer times 2^(e - 52), e in -40..20 (binary64-exact), real or with a
//   binary64-exact imaginary part; its partner at the anchor plus 2^(e - k), k in 53..200, along 1, i or 1 + i (70%),
//   or plus 2^(e - k)/3, k in 50..120, a non-dyadic offset (30%);
// - one to six other planted roots (coordinate draws of magnitude 2^-4 to 2^4), so that the solve works a packet;
//   degree at most 32;
// - a real polynomial in half the seeds (real anchors and offsets; the other roots in conjugate pairs or real), and a
//   real-only query in a fifth of those;
// - the query: Auto (70%) or an explicit method; the accuracy goal (60%) or approximation (40%); goal bits 52, 53
//   (twice as likely), 54, 55, 100 or 333; the precision cap 4,096 bits (60%) or the library's default; a work
//   budget of 100,000 units and 65,536 cells, so that the cap and the budgets are not the story.
inline constexpr uint64_t edge53_stream = UINT64_C(0x6564676535332d66); // "edge53-f"
inline constexpr int edge53_max_degree = 32;
inline Spec generate_edge53(uint64_t seed) {
    Rng r(seed ^ edge53_stream);
    Spec s;
    s.seed = seed;
    s.oracle = Oracle::Planted;
    const bool real = r.percent(50);
    Planter p{r, real, {}};
    std::ostringstream family;
    family << "edge53";
    constexpr int binary64_bits = 53;
    const int anchors = r.percent(70) ? 1 : 2;
    for (int a = 0; a < anchors; ++a) {
        const int e = static_cast<int>(r.range(-40, 20));
        auto exact = [&] { // an odd 53-bit integer times 2^(e - 52): binary64-exact, its last bit set
            const int64_t m = static_cast<int64_t>(r.below(UINT64_C(1) << (binary64_bits - 1))) |
                              (int64_t(1) << (binary64_bits - 1)) | 1;
            return Rational(Integer(r.percent(50) ? m : -m)) * power2(e - (binary64_bits - 1));
        };
        const Rational re = exact();
        const Complex z{re, real || r.percent(60) ? Rational(0) : exact()};
        Complex offset;
        if (r.percent(70)) {
            const int k = static_cast<int>(r.range(binary64_bits, 200));
            const Rational d = power2(e - k);
            const auto direction = real ? 0 : r.below(3);
            offset = direction == 0 ? Complex{d, 0} : direction == 1 ? Complex{0, d} : Complex{d, d};
            family << "+pair(2^-" << k << (direction == 0 ? "" : direction == 1 ? " i" : " (1+i)") << ")";
        } else {
            const int k = static_cast<int>(r.range(50, 120));
            offset = Complex{power2(e - k) / 3, 0};
            family << "+pair(2^-" << k << "/3)";
        }
        p.add(z);
        p.add(z + offset);
    }
    const int others = static_cast<int>(r.range(1, 6));
    for (int j = 0; j < others; ++j) {
        const int e = static_cast<int>(r.range(-4, 4));
        const Complex w{coordinate(r, e), real && r.percent(40) ? Rational(0) : coordinate(r, e)};
        if (p.size() + p.cost(w) > edge53_max_degree) break;
        p.add(w);
    }
    s.roots = p.roots;
    family << "+" << others << " others" << (real ? ";real" : "");
    s.family = family.str();
    s.method = r.percent(70) ? CPS_METHOD_AUTO : static_cast<cps_method>(1 + r.below(3));
    s.goal = r.percent(60) ? CPS_ROOT_ACCURACY : CPS_ROOT_APPROXIMATION;
    s.real_only = real && s.method == CPS_METHOD_AUTO && r.percent(20);
    s.bits = r.weighted<uint64_t>({{52, 1}, {53, 2}, {54, 1}, {55, 1}, {100, 1}, {333, 1}});
    s.max_work_bits = r.percent(60) ? 4096 : CPS_DEFAULT_MAX_WORK_BITS;
    s.work_units = 100000;
    s.live_cells = 65536;
    s.determinism = r.percent(5);
    return s;
}

// The origin family (--origin, harness fuzz-11): the review's case T12 (dev/docs/adversarial-corpus-note-20261016.md).
// A tiny root beside a far cloud: a box accurate relative to the centred cloud and useless at the tiny root (centering
// or power lifting publishing in the wrong frame) fails the accuracy contract, which the checker judges at each true
// root in the input's coordinates. From the family's own stream:
// - a root 2^-t/3, t one of 52, 64, 128 and 1,024 (real, or along i in a third of the non-real seeds); an exact zero
//   root beside it in 30%;
// - a cloud of 2 to 8 roots C + j 2^-u/3 (j = 1, 2, ...), C = 1,000 or 2^20, u in 20..80;
// - one or two ordinary planted roots (coordinate draws) in half the seeds;
// - every root turned by i or (3 + 4i)/5 in 30% (the polynomial no longer real); the coefficients multiplied by a
//   rational scalar or by 2^k, k in -200..200, in half the seeds (the roots unchanged);
// - the query: Auto (70%) or an explicit method; the accuracy goal (45%), approximation (35%) or, under Auto, the
//   cluster radius 2^-N, N in 10..120 (20%); goal bits 53 or 100; the precision cap 4,096 bits (60%) or the
//   library's default; 100,000 work units and 65,536 cells; a real-only query in a fifth of the real Auto seeds.
inline constexpr uint64_t origin_stream = UINT64_C(0x6f726967696e2d66); // "origin-f"
inline Spec generate_origin(uint64_t seed) {
    Rng r(seed ^ origin_stream);
    Spec s;
    s.seed = seed;
    s.oracle = Oracle::Planted;
    std::ostringstream family;
    const int t = r.weighted<int>({{52, 1}, {64, 1}, {128, 1}, {1024, 1}});
    const Rational tiny = power2(-t) / 3;
    const auto turn = r.below(10);
    const bool rotated = turn < 3;
    const Complex rotation = turn == 0 ? Complex{0, 1} : Complex{Rational(3) / 5, Rational(4) / 5};
    std::vector<Complex> roots;
    const bool tiny_along_i = rotated ? false : r.percent(15);
    roots.push_back(tiny_along_i ? Complex{0, tiny} : Complex{tiny, 0});
    if (tiny_along_i) roots.push_back(Complex{0, -tiny});
    family << "origin(2^-" << t << "/3" << (tiny_along_i ? " i" : "") << ")";
    if (r.percent(30)) {
        roots.push_back(Complex{0, 0});
        family << "+zero";
    }
    const Rational centre = r.percent(50) ? Rational(1000) : power2(20);
    const int u = static_cast<int>(r.range(20, 80));
    const int cloud = static_cast<int>(r.range(2, 8));
    for (int j = 1; j <= cloud; ++j) roots.push_back(Complex{centre + Rational(j) * power2(-u) / 3, 0});
    family << "+cloud(" << cloud << " at " << (centre == Rational(1000) ? "1000" : "2^20") << ", 2^-" << u << "/3)";
    if (r.percent(50)) {
        const int ordinary = static_cast<int>(r.range(1, 2));
        for (int j = 0; j < ordinary; ++j) {
            const int e = static_cast<int>(r.range(-4, 4));
            const Complex w{coordinate(r, e), coordinate(r, e)};
            roots.push_back(w);
            if (!rotated) roots.push_back(conjugate(w)); // keeps the unrotated polynomial real
        }
        family << "+" << ordinary << " ordinary";
    }
    if (rotated) {
        for (auto &z : roots) z = z * rotation;
        family << (turn == 0 ? ";turned i" : ";turned (3+4i)/5");
    }
    if (r.percent(50)) {
        if (r.percent(50)) {
            const int k = static_cast<int>(r.range(-200, 200));
            s.multiplier = Complex(power2(k));
            family << ";c=2^" << k;
        } else {
            const int64_t a = r.range(1, 999), b = r.range(1, 999);
            s.multiplier = Complex(Rational(Integer(r.percent(50) ? a : -a)) / Rational(Integer(b)));
            family << ";c=rational";
        }
    }
    s.roots = roots;
    s.family = family.str();
    const bool real = !rotated;
    s.method = r.percent(70) ? CPS_METHOD_AUTO : static_cast<cps_method>(1 + r.below(3));
    const auto g = r.below(100);
    s.goal = s.method == CPS_METHOD_AUTO && g < 20 ? CPS_CLUSTER_RADIUS : g < 65 ? CPS_ROOT_ACCURACY : CPS_ROOT_APPROXIMATION;
    s.cluster_exponent = static_cast<int>(r.range(10, 120));
    s.real_only = real && s.method == CPS_METHOD_AUTO && s.goal != CPS_CLUSTER_RADIUS && r.percent(20);
    s.bits = r.percent(50) ? 53 : 100;
    s.max_work_bits = r.percent(60) ? 4096 : CPS_DEFAULT_MAX_WORK_BITS;
    s.work_units = 100000;
    s.live_cells = 65536;
    s.determinism = r.percent(5);
    return s;
}

// The algebraic family (--algebraic, harness fuzz-11): the review's case T09 (dev/docs/adversarial-corpus-note-20261016.md).
// Every planted root of the other families is Gaussian-rational, so binary64 can land on it exactly and a region can
// close by an exact coincidence; these roots have no finite representation. Factors irreducible over Q(i):
//   q2(z) = (z - a)^2 - 2 c^2 (roots a +- sqrt(2) c) and q4(z) = (z - a)^4 + 2 c^4 (roots a + c 2^(1/4) e^(i pi (2k+1)/4)),
// with c = d w: a a Gaussian-rational centre (coordinate draws, magnitude 2^-3 to 2^3), d = 2^-k (k in 0..12, or in
// 20..120 for a close pair or a tight cluster in half the factors), w a rotation (1, i, (3+4i)/5 or (5+12i)/13). The
// polynomial (Oracle::Factored) is the product of:
// - two to five such factors, each repeated (multiplicity 2 or 3) in 30%;
// - up to four linear factors z - b (b a coordinate draw), repeated in 15%;
// with distinct factors (a factor equal to one drawn already is skipped), degree at most 32; a rational or 2^k
// multiplier in a third of the seeds. The query: Auto (70%) or an explicit method; the accuracy goal (40%),
// approximation (40%) or, under Auto, the cluster radius 2^-N, N in 10..120 (20%); goal bits 16, 53 or 100; the
// precision cap 64, 128, 256, 1,024, 4,096 bits or the library's default; 100,000 work units, 65,536 cells; no
// real-only query (its expected real count would need another oracle).
inline constexpr uint64_t algebraic_stream = UINT64_C(0x616c67656272612d); // "algebra-"
inline constexpr int algebraic_max_degree = 32;
inline Spec generate_algebraic(uint64_t seed) {
    Rng r(seed ^ algebraic_stream);
    Spec s;
    s.seed = seed;
    s.oracle = Oracle::Factored;
    std::ostringstream family;
    family << "algebraic";
    int degree = 0;
    const auto add = [&](ComplexPolynomial factor, int multiplicity, const std::string &label) {
        const int d = static_cast<int>(factor.size()) - 1;
        if (degree + d * multiplicity > algebraic_max_degree) return;
        const auto same = [](const ComplexPolynomial &u, const ComplexPolynomial &v) {
            if (u.size() != v.size()) return false;
            for (size_t i = 0; i < u.size(); ++i)
                if (u[i].real != v[i].real || u[i].imag != v[i].imag) return false;
            return true;
        };
        for (const auto &[old, m] : s.factors)
            if (same(old, factor)) return;
        s.factors.emplace_back(std::move(factor), multiplicity);
        degree += d * multiplicity;
        family << "+" << label << (multiplicity > 1 ? "^" + std::to_string(multiplicity) : "");
    };
    const std::array<Complex, 4> rotations{{Complex{1}, Complex{0, 1}, Complex{Rational(3) / 5, Rational(4) / 5},
                                            Complex{Rational(5) / 13, Rational(12) / 13}}};
    const int irreducible = static_cast<int>(r.range(2, 5));
    for (int j = 0; j < irreducible; ++j) {
        const int e = static_cast<int>(r.range(-3, 3));
        const Complex a{coordinate(r, e), r.percent(50) ? Rational(0) : coordinate(r, e)};
        const bool tight = r.percent(50);
        const int k = tight ? static_cast<int>(r.range(20, 120)) : static_cast<int>(r.range(0, 12));
        const Complex c = rotations[r.below(rotations.size())] * Complex(power2(-k));
        const Complex c2 = c * c, a2 = a * a;
        const int multiplicity = r.percent(30) ? static_cast<int>(r.range(2, 3)) : 1;
        if (r.percent(60)) { // (z - a)^2 - 2 c^2 = z^2 - 2 a z + (a^2 - 2 c^2)
            add({a2 + Complex(Rational(-2)) * c2, Complex(Rational(-2)) * a, Complex(1)}, multiplicity,
                "q2(2^-" + std::to_string(k) + ")");
        } else { // (z - a)^4 + 2 c^4
            const Complex a3 = a2 * a, a4 = a2 * a2, c4 = c2 * c2;
            add({a4 + Complex(Rational(2)) * c4, Complex(Rational(-4)) * a3, Complex(Rational(6)) * a2,
                 Complex(Rational(-4)) * a, Complex(1)},
                multiplicity, "q4(2^-" + std::to_string(k) + ")");
        }
    }
    const int linear = static_cast<int>(r.range(0, 4));
    for (int j = 0; j < linear; ++j) {
        const int e = static_cast<int>(r.range(-3, 3));
        const Complex b{coordinate(r, e), r.percent(50) ? Rational(0) : coordinate(r, e)};
        add({Complex{-b.real, -b.imag}, Complex(1)}, r.percent(15) ? 2 : 1, "lin");
    }
    if (s.factors.empty()) add({Complex(Rational(-2)), Complex(0), Complex(1)}, 1, "q2(1)"); // z^2 - 2
    if (r.percent(33)) {
        if (r.percent(50)) {
            const int k = static_cast<int>(r.range(-60, 60));
            s.multiplier = Complex(power2(k));
            family << ";c=2^" << k;
        } else {
            const int64_t m = r.range(1, 999), n = r.range(1, 999);
            s.multiplier = Complex(Rational(Integer(r.percent(50) ? m : -m)) / Rational(Integer(n)));
            family << ";c=rational";
        }
    }
    s.family = family.str();
    s.method = r.percent(70) ? CPS_METHOD_AUTO : static_cast<cps_method>(1 + r.below(3));
    const auto g = r.below(100);
    s.goal = s.method == CPS_METHOD_AUTO && g < 20 ? CPS_CLUSTER_RADIUS : g < 60 ? CPS_ROOT_ACCURACY : CPS_ROOT_APPROXIMATION;
    s.cluster_exponent = static_cast<int>(r.range(10, 120));
    s.bits = r.weighted<uint64_t>({{16, 1}, {53, 2}, {100, 1}});
    s.max_work_bits = r.weighted<uint64_t>({{64, 1}, {128, 1}, {256, 1}, {1024, 1}, {4096, 2}, {CPS_DEFAULT_MAX_WORK_BITS, 2}});
    s.work_units = 100000;
    s.live_cells = 65536;
    s.determinism = r.percent(5);
    return s;
}

// Exact arithmetic over Q(i) for the composed family's square-free proof: an inverse, the derivative, the remainder
// of a division and the degree of a gcd (Euclid's algorithm on exact rationals; degree 32 at most here).
inline Complex complex_inverse(const Complex &z) {
    const Rational norm = z.real * z.real + z.imag * z.imag;
    return {z.real / norm, -z.imag / norm};
}
inline ComplexPolynomial derivative_of(const ComplexPolynomial &p) {
    ComplexPolynomial d;
    for (size_t j = 1; j < p.size(); ++j) d.push_back(Complex(Rational(static_cast<int64_t>(j))) * p[j]);
    return d;
}
inline void trim(ComplexPolynomial &p) {
    while (!p.empty() && zero(p.back())) p.pop_back();
}
inline ComplexPolynomial remainder_of(ComplexPolynomial a, const ComplexPolynomial &b) {
    trim(a);
    const Complex lead = complex_inverse(b.back());
    while (a.size() >= b.size()) {
        const Complex q = a.back() * lead;
        const size_t shift = a.size() - b.size();
        for (size_t j = 0; j < b.size(); ++j) a[shift + j] = a[shift + j] + Complex{-q.real, -q.imag} * b[j];
        a.pop_back(); // its coefficient is now zero
        trim(a);
    }
    return a;
}
inline int gcd_degree(ComplexPolynomial a, ComplexPolynomial b) {
    trim(a);
    trim(b);
    while (!b.empty()) {
        ComplexPolynomial r = remainder_of(a, b);
        a = std::move(b);
        b = std::move(r);
    }
    return static_cast<int>(a.size()) - 1;
}

// The composed family (--composed, harness fuzz-11): the review's case T11 (dev/docs/adversarial-corpus-note-20261016.md),
// almost structured inputs p = z^s q(z^k) + e z^j, which an exact reduction (power_composition) must take only when the
// support's gcd allows it: q of degree 1 to 6 with planted roots (coordinate draws, nonzero; one repeated, multiplicity 2
// or 3, in 40%; one near the negative real axis, -|r| + 2^-30 i, in 30%), s in 0..4, k in {2, 3, 4, 5, 6, 8}; e = 0
// (40%), a coarse rational (30%), 2^-53 (15%) or 2^-200 (15%), j in 0..degree - 1; degree at most 32 for e = 0 and at
// most 12 otherwise (the check's exact count on the one dense factor g grows past any campaign's time above that: a
// trial of 27 October spent minutes on single degree-20 checks). Its oracle
// (Oracle::Factored):
// - e = 0: the factors z (multiplicity s) and z^k - r for each distinct root r of q, with r's multiplicity: each is
//   square-free (r is not zero) and they are pairwise coprime, so their counts decide the roots and multiplicities;
// - e != 0: p = z^m g with m = min(s, j) and g(0) != 0, proved square-free here (the gcd of g and g' has degree 0);
//   the factors z (multiplicity m) and g. A g that is not square-free (rare) is replaced by the e = 0 case.
// The query: Auto (70%) or an explicit method; accuracy (40%), approximation (40%) or, under Auto, the cluster radius
// 2^-N, N in 10..120 (20%); goal bits 16, 53 or 100; the precision cap 128, 1,024, 4,096 bits or the library's
// default; 100,000 work units, 65,536 cells; no real-only query.
inline constexpr uint64_t composed_stream = UINT64_C(0x636f6d706f73652d); // "compose-"
inline constexpr int composed_max_degree = 32, composed_perturbed_max_degree = 12;
inline Spec generate_composed(uint64_t seed) {
    Rng r(seed ^ composed_stream);
    Spec s;
    s.seed = seed;
    s.oracle = Oracle::Factored;
    const std::array<int, 6> powers{{2, 3, 4, 5, 6, 8}};
    const int k = powers[r.below(powers.size())];
    const int shift = static_cast<int>(r.range(0, 4));
    const auto kind = r.below(100);
    const int cap = kind < 40 ? composed_max_degree : composed_perturbed_max_degree;
    const int room = (cap - shift) / k;
    const int qdegree = static_cast<int>(r.range(1, std::max(1, std::min(6, room))));
    // q's distinct roots with multiplicities (their total is qdegree).
    std::vector<std::pair<Complex, int>> q;
    int planted = 0;
    const bool repeat = qdegree >= 2 && r.percent(40);
    while (planted < qdegree) {
        const int e = static_cast<int>(r.range(-2, 2));
        Complex z{coordinate(r, e), r.percent(50) ? Rational(0) : coordinate(r, e)};
        if (q.empty() && r.percent(30)) z = Complex{-absolute(z.real) - Rational(1, 2), power2(-30)}; // near the negative axis
        if (zero(z)) continue;
        bool duplicate = false;
        for (const auto &[w, m] : q) duplicate = duplicate || (w.real == z.real && w.imag == z.imag);
        if (duplicate) continue;
        const int m = repeat && q.empty() ? std::min(qdegree, static_cast<int>(r.range(2, 3))) : 1;
        q.emplace_back(z, std::min(m, qdegree - planted));
        planted += q.back().second;
    }
    const int degree = shift + k * qdegree;
    // p0 = z^shift q(z^k), with q monic.
    ComplexPolynomial qpoly{Complex(1)};
    for (const auto &[z, m] : q)
        for (int c = 0; c < m; ++c) qpoly = multiply(qpoly, {Complex{-z.real, -z.imag}, Complex(1)});
    ComplexPolynomial p0(static_cast<size_t>(degree) + 1, Complex(0));
    for (size_t j = 0; j < qpoly.size(); ++j) p0[static_cast<size_t>(shift) + j * static_cast<size_t>(k)] = qpoly[j];
    Rational epsilon = kind < 40 ? Rational(0) : kind < 70 ? Rational(Integer(r.range(1, 9))) / Rational(Integer(r.range(2, 11)))
                     : kind < 85 ? power2(-53) : power2(-200);
    const int j = static_cast<int>(r.range(0, degree - 1));
    std::ostringstream family;
    family << "composed(z^" << shift << " q(z^" << k << "), deg q " << qdegree << (repeat ? ", repeated" : "") << ")";
    if (epsilon != 0) {
        ComplexPolynomial p = p0;
        p[static_cast<size_t>(j)] = p[static_cast<size_t>(j)] + Complex(epsilon);
        const int m = std::min(shift, j);
        const ComplexPolynomial g(p.begin() + m, p.end());
        if (gcd_degree(g, derivative_of(g)) == 0) {
            if (m) s.factors.emplace_back(ComplexPolynomial{Complex(0), Complex(1)}, m);
            s.factors.emplace_back(g, 1);
            family << "+eps z^" << j << (kind < 70 ? " coarse" : kind < 85 ? " 2^-53" : " 2^-200");
        } else
            epsilon = 0;
    }
    if (epsilon == 0) {
        if (shift) s.factors.emplace_back(ComplexPolynomial{Complex(0), Complex(1)}, shift);
        for (const auto &[z, m] : q) {
            ComplexPolynomial f(static_cast<size_t>(k) + 1, Complex(0));
            f[0] = Complex{-z.real, -z.imag};
            f[static_cast<size_t>(k)] = Complex(1);
            s.factors.emplace_back(std::move(f), m);
        }
    }
    s.family = family.str();
    s.method = r.percent(70) ? CPS_METHOD_AUTO : static_cast<cps_method>(1 + r.below(3));
    const auto g = r.below(100);
    s.goal = s.method == CPS_METHOD_AUTO && g < 20 ? CPS_CLUSTER_RADIUS : g < 60 ? CPS_ROOT_ACCURACY : CPS_ROOT_APPROXIMATION;
    s.cluster_exponent = static_cast<int>(r.range(10, 120));
    s.bits = r.weighted<uint64_t>({{16, 1}, {53, 2}, {100, 1}});
    s.max_work_bits = r.weighted<uint64_t>({{128, 1}, {1024, 1}, {4096, 2}, {CPS_DEFAULT_MAX_WORK_BITS, 2}});
    s.work_units = 100000;
    s.live_cells = 65536;
    s.determinism = r.percent(5);
    return s;
}

// The inexact family: the exact case of the same seed with a partner
// polynomial planted beside it (Spec::inexact). The exact generation above is
// untouched, and the family draws from its own stream, so a seed names the
// same exact case with or without it.
//
// - Each nonzero root moves by (a + b i) 2^-k of itself, a and b in -3..3
//   (30% of the roots stay; at least one moves). A real polynomial keeps a
//   real partner: real roots move along the axis and conjugates move
//   together. Copies of a repeated root move independently, so the partner
//   has distinct roots where the case has a multiple one.
// - k is drawn so that the radii forbid, allow, or sit at the goal.
// - The explicit methods and real-only queries refuse inexact input; one case
//   in ten of each keeps the query, to check the refusal.
// - A quarter of the cases run at the library's default precision cap
//   (CPS_DEFAULT_MAX_WORK_BITS), where a goal the radii forbid costs the most.
inline constexpr uint64_t inexact_stream = UINT64_C(0x696e6578616374); // "inexact"
inline Spec generate_inexact(uint64_t seed) {
    Spec s = generate(seed);
    if (s.oracle != Oracle::Planted)
        return s; // Mignotte cases have no planted roots to move: they stay exact
    Rng r(seed ^ inexact_stream);
    const int k = r.weighted<int>({{static_cast<int>(r.range(8, 40)), 40},
                                   {static_cast<int>(r.range(41, 120)), 40},
                                   {static_cast<int>(r.range(121, 400)), 20}});
    const bool real_polynomial = real_coefficients(build(s));
    const Rational unit = power2(-k);
    std::vector<Complex> partner(s.roots.size());
    std::vector<char> done(s.roots.size(), 0);
    bool moved = false;
    for (size_t i = 0; i < s.roots.size(); ++i) {
        if (done[i]) continue;
        const Complex &z = s.roots[i];
        Complex step{0};
        const bool zero_root = z.real == 0 && z.imag == 0;
        if (!zero_root && (r.percent(70) || (!moved && i + 1 == s.roots.size()))) {
            const bool axis = real_polynomial && z.imag == 0;
            Rational a(r.range(-3, 3)), b(axis ? 0 : r.range(-3, 3));
            if (a == 0 && b == 0) a = 1;
            step = Complex(a * unit, b * unit);
            moved = true;
        }
        partner[i] = z + z * step;
        done[i] = 1;
        if (real_polynomial && z.imag != 0)
            for (size_t j = i + 1; j < s.roots.size(); ++j)
                if (!done[j] && s.roots[j].real == z.real && s.roots[j].imag == -z.imag) {
                    partner[j] = conjugate(partner[i]);
                    done[j] = 1;
                    break;
                }
    }
    if (!moved)
        return s; // every root is zero: nothing to move, the case stays exact
    s.inexact = k;
    s.partner = std::move(partner);
    s.family += ";inexact(2^-" + std::to_string(k) + ")";
    if (s.method != CPS_METHOD_AUTO && !r.percent(10)) s.method = CPS_METHOD_AUTO;
    if (s.real_only && !r.percent(10)) s.real_only = false;
    if (r.percent(25)) s.max_work_bits = CPS_DEFAULT_MAX_WORK_BITS;
    return s;
}

// The D1 family (--d1): the exact case of the same seed, made a real polynomial where it is not one (every
// unpaired non-real root gets its conjugate, pairs dropped from the end past degree 48, a non-real multiplier
// replaced by a real one: the planted real roots and their multiplicities are unchanged, and they are the
// oracle), queried on the real-only path with the isolation goal (D1 Part 1) or a closed region (Part 2). The
// exact generation above is untouched, and the family draws from its own stream, so a seed names the same
// standing case with or without it.
//
// - isolation (30%): goal CPS_ROOT_ISOLATION, no region.
// - region (45%): goal isolation, accuracy or approximation; each real end independently unbounded (20%), a
//   planted real root exactly (25%: dyadic or p/q as planted), a planted root +- 2^-k (25%: k from 1 to the
//   case's precision cap, above it in a tenth of them: the stop at the cap is exercised there), or a plain
//   rational (30%); the ends swapped if reversed, equal ends allowed. The imaginary sides unbounded or a band
//   holding 0, and in 15% a rectangle off the axis (an empty, complete answer with no work).
// - session (15%; the name kept so that the cases are drawn as before): isolation first, then the accuracy or
//   approximation goal, each a fresh solve, holding the same roots. (Until 18 October 2026 one session, which
//   then had to refuse a changed region and the goal set back.)
// - refusal (10%): the query a D1 setting must refuse (D1Refusal), or a malformed region at the setter.
inline constexpr uint64_t d1_stream = UINT64_C(0x64312d66757a7a); // "d1-fuzz"
inline std::string end_text(const Rational &x) { return x.str(); } // "p/q" or "p": the coefficient grammar
inline void make_real(Spec &s) {
    std::vector<Complex> roots;
    std::vector<char> paired(s.roots.size(), 0);
    for (size_t i = 0; i < s.roots.size(); ++i) {
        const Complex &z = s.roots[i];
        if (z.imag == 0) {
            roots.push_back(z);
            continue;
        }
        if (paired[i]) continue;
        // Its conjugate, if planted and not yet paired, pairs with it; otherwise the conjugate is added.
        bool found = false;
        for (size_t j = i + 1; j < s.roots.size() && !found; ++j)
            if (!paired[j] && s.roots[j].real == z.real && s.roots[j].imag == -z.imag) {
                paired[j] = 1;
                found = true;
            }
        if (roots.size() + 2 > 48) continue; // past the fuzzer's degree: the pair is dropped
        roots.push_back(z);
        roots.push_back(conjugate(z));
    }
    s.roots = std::move(roots);
    if (s.multiplier.imag != 0) s.multiplier = Complex(s.multiplier.real != 0 ? s.multiplier.real : Rational(1));
}
// The D1 query drawn on a Spec: the D1 family's own (generate_d1, below) and the axis family's D1 variant
// (generate_axis). Verbatim from generate_d1 of 16 October 2026: the draws and their order are the D1
// family's, so its cases are unchanged (checked by --describe over its seeds).
inline void draw_d1_query(Spec &s, Rng &r, D1Variant variant, D1Refusal refusal) {
    s.family += ";d1";
    s.method = CPS_METHOD_AUTO;
    s.real_only = true;
    s.d1 = static_cast<int>(variant);
    s.d1_refusal = static_cast<int>(refusal);
    // The planted real roots, distinct, for the ends.
    std::vector<Rational> real_roots;
    if (s.oracle == Oracle::Planted)
        for (const auto &z : s.roots)
            if (z.imag == 0 && std::find(real_roots.begin(), real_roots.end(), z.real) == real_roots.end())
                real_roots.push_back(z.real);
    static const char *plain[] = {"1/3", "7/5", "-2/7", "22/7", "0", "1", "-1", "3/2", "-5/4", "9"};
    const auto draw_end = [&](std::optional<Rational> &out) {
        switch (r.weighted<int>({{0, 20}, {1, 25}, {2, 25}, {3, 30}})) {
        case 0: out.reset(); return;
        case 1:
            if (!real_roots.empty()) {
                out = real_roots[r.below(real_roots.size())];
                return;
            }
            [[fallthrough]];
        case 2:
            if (!real_roots.empty()) {
                const Rational x = real_roots[r.below(real_roots.size())];
                const auto cap = static_cast<int64_t>(s.max_work_bits);
                const int k = static_cast<int>(r.percent(10) ? cap + r.range(1, 64) : r.range(1, cap));
                out = x + (r.percent(50) ? Rational(1) : Rational(-1)) * power2(-k);
                return;
            }
            [[fallthrough]];
        default: out = parse_rational(plain[r.below(sizeof plain / sizeof plain[0])]); return;
        }
    };
    // `required`: a region there must be (the refusal probes), so four unbounded sides are drawn again.
    const auto draw_region = [&](bool off_axis_allowed, bool required = false) {
      for (;;) {
        std::optional<Rational> lower, upper;
        draw_end(lower);
        draw_end(upper);
        if (lower && upper && *lower > *upper) std::swap(lower, upper);
        s.region_lower = lower ? end_text(*lower) : "";
        s.region_upper = upper ? end_text(*upper) : "";
        const auto band = r.weighted<int>({{0, 60}, {1, 15}, {2, 5}, {3, 5}, {4, off_axis_allowed ? 15u : 0u}});
        s.region_imag_lower = band == 1 ? "-1" : band == 3 ? "0" : band == 4 ? (r.percent(50) ? "1" : "-3") : "";
        s.region_imag_upper = band == 1 ? "1" : band == 2 ? "0" : band == 4 ? (s.region_imag_lower == "1" ? "2" : "-1/32") : "";
        if (required && !s.region_set()) continue;
        s.family += "(region [" + (s.region_lower.empty() ? "-inf" : s.region_lower) + ", " +
                    (s.region_upper.empty() ? "inf" : s.region_upper) + "]" +
                    (band == 4 ? " off the axis" : "") + ")";
        return;
      }
    };
    const auto pick_goal = [&] {
        return r.weighted<cps_goal>({{CPS_ROOT_ISOLATION, 30}, {CPS_ROOT_ACCURACY, 35}, {CPS_ROOT_APPROXIMATION, 35}});
    };
    switch (variant) {
    case D1Variant::Isolation:
        s.goal = CPS_ROOT_ISOLATION;
        s.family += "(isolation)";
        break;
    case D1Variant::Region:
        s.goal = pick_goal();
        draw_region(true);
        s.family += goal_name(s.goal);
        break;
    case D1Variant::Session:
        // The first call's goal is isolation; s.goal is the tightened call's.
        s.goal = r.percent(50) ? CPS_ROOT_ACCURACY : CPS_ROOT_APPROXIMATION;
        if (r.percent(50)) draw_region(false);
        s.family += "(session to " + std::string(goal_name(s.goal)) + ")";
        break;
    case D1Variant::Refusal:
        switch (refusal) {
        case D1Refusal::RegionWithoutRealOnly:
            s.real_only = false;
            s.goal = r.percent(50) ? CPS_ROOT_ACCURACY : CPS_ROOT_APPROXIMATION;
            draw_region(false, true);
            break;
        case D1Refusal::IsolationWithoutRealOnly:
            s.real_only = false;
            s.goal = CPS_ROOT_ISOLATION;
            break;
        case D1Refusal::ExplicitMethod:
            s.method = static_cast<cps_method>(r.range(1, 3));
            if (r.percent(50)) s.goal = CPS_ROOT_ISOLATION;
            else {
                s.goal = CPS_ROOT_ACCURACY;
                draw_region(false, true);
            }
            break;
        case D1Refusal::ClusterGoal:
            s.goal = CPS_CLUSTER_RADIUS;
            draw_region(false, true);
            break;
        case D1Refusal::NonRealInput:
            s.goal = r.percent(50) ? CPS_ROOT_ISOLATION : CPS_ROOT_ACCURACY;
            if (s.goal != CPS_ROOT_ISOLATION) draw_region(false, true);
            break;
        case D1Refusal::InexactInput:
            s.goal = CPS_ROOT_ISOLATION;
            break;
        case D1Refusal::InexpandableEnd:
            s.goal = CPS_ROOT_ACCURACY;
            s.region_lower = "1@100000";
            break;
        case D1Refusal::SetterReversed:
        case D1Refusal::SetterMalformed:
            s.goal = CPS_ROOT_ACCURACY; // the solve after the refused setter: the whole axis
            break;
        case D1Refusal::None: break;
        }
        s.family += "(refusal " + std::to_string(s.d1_refusal) + ")";
        break;
    case D1Variant::None: break;
    }
}
inline Spec generate_d1(uint64_t seed) {
    Spec s = generate(seed);
    Rng r(seed ^ d1_stream);
    const auto variant = r.weighted<D1Variant>(
        {{D1Variant::Isolation, 30}, {D1Variant::Region, 45}, {D1Variant::Session, 15}, {D1Variant::Refusal, 10}});
    const bool real_input = real_coefficients(build(s));
    auto refusal = D1Refusal::None;
    if (variant == D1Variant::Refusal) {
        refusal = r.weighted<D1Refusal>({{D1Refusal::RegionWithoutRealOnly, 1}, {D1Refusal::IsolationWithoutRealOnly, 1},
                                         {D1Refusal::ExplicitMethod, 1}, {D1Refusal::ClusterGoal, 1},
                                         {D1Refusal::NonRealInput, 1}, {D1Refusal::InexactInput, 1},
                                         {D1Refusal::InexpandableEnd, 1}, {D1Refusal::SetterReversed, 1},
                                         {D1Refusal::SetterMalformed, 1}});
        if (refusal == D1Refusal::NonRealInput && (real_input || s.oracle != Oracle::Planted))
            refusal = D1Refusal::RegionWithoutRealOnly; // nothing non-real to refuse
        if (refusal == D1Refusal::InexactInput) {
            const Spec inexact = generate_inexact(seed);
            if (inexact.inexact && real_input) {
                s.inexact = inexact.inexact;
                s.partner = inexact.partner;
            } else
                refusal = D1Refusal::ExplicitMethod;
        }
    }
    if (s.oracle == Oracle::Planted && !real_input && refusal != D1Refusal::NonRealInput) make_real(s);
    draw_d1_query(s, r, variant, refusal);
    return s;
}

// The axis family (--axis): the review's case T13 (codex-notes/review-20261014/adversarial-cases.md;
// docs/adversarial-corpus-note-20261016.md), near-real conjugate pairs against true real roots. A real polynomial
// built from the family's own stream; the standing case of the seed lends its precision cap, budgets, goal bits
// and determinism, nothing else, so the caps and budgets have the standing distribution.
//
// - one to three pairs a +- i 2^-k: k up to the case's precision cap (60%), in the cap's upper half (25%), or
//   past the cap by 1 to 64 (15%); a at a chart boundary (0, +-1/2, +-1, +-2^m with m in -8..8) in half the
//   pairs, else a planted-style rational (coordinate);
// - up to six true real roots: at a pair's own a (30%), at a +- 2^-m with m up to that pair's k (30%), or
//   anywhere (40%); a repeated real factor, multiplicity 2 to 4, in 30% (one of the real roots, or a pair's a);
// - ordinary pairs away from the axis in 30% (one to four);
// - a transform: none (70%), a rotation by i or (3+4i)/5 (15%: the pairs beside the rotated axis, the
//   polynomial no longer real), or a non-real scalar on the real polynomial (15%: the roots unchanged, the
//   coefficients not real); then the standing power-of-two scaling in 30%;
// - the query: Auto (70%) or an explicit method; the goal accuracy (40%), approximation (40%) or the cluster
//   radius (20%, Auto only); on exact real input a real-only query in 40% of the Auto seeds, half of them
//   through the D1 drawing (draw_d1_query: isolation, a region or the session, checked as D1 is); on the rotated
//   or scalar variant a real-only query in a tenth of the Auto seeds, to check the refusal.
inline constexpr uint64_t axis_stream = UINT64_C(0x617869732d66757a); // "axis-fuz"
inline Spec generate_axis(uint64_t seed) {
    const Spec standing = generate(seed);
    Rng r(seed ^ axis_stream);
    Spec s;
    s.seed = seed;
    s.bits = standing.bits;
    s.max_work_bits = standing.max_work_bits;
    s.work_units = standing.work_units;
    s.live_cells = standing.live_cells;
    s.determinism = standing.determinism;
    s.axis = 1;
    const auto cap = static_cast<int64_t>(s.max_work_bits);
    std::ostringstream family;
    family << "axis(";
    // The pairs. (Braced initializers fix the order of the draws, as in generate.)
    const int pairs = r.weighted<int>({{1, 50}, {2, 30}, {3, 20}});
    std::vector<std::pair<Rational, int>> near; // (a, k)
    for (int j = 0; j < pairs; ++j) {
        const int k = r.weighted<int>({{static_cast<int>(r.range(1, cap)), 60},
                                       {static_cast<int>(r.range(cap / 2 + 1, cap)), 25},
                                       {static_cast<int>(cap + r.range(1, 64)), 15}});
        Rational a;
        if (r.percent(50)) {
            const Rational sign(r.percent(50) ? 1 : -1);
            switch (r.below(4)) {
            case 0: a = Rational(0); break;
            case 1: a = sign / Rational(2); break;
            case 2: a = sign; break;
            default: a = sign * power2(static_cast<int>(r.range(-8, 8))); break;
            }
        } else a = coordinate(r, 0);
        near.emplace_back(a, k);
        s.roots.emplace_back(a, power2(-k));
        s.roots.emplace_back(a, -power2(-k));
        s.axis_k = std::max(s.axis_k, k);
        family << (j ? "," : "") << "a=" << a.str() << ",k=" << k;
    }
    family << ")";
    // The true real roots.
    std::vector<Rational> reals;
    const int count = static_cast<int>(r.range(0, 6));
    for (int j = 0; j < count; ++j) {
        const auto &[a, k] = near[r.below(near.size())];
        const auto kind = r.weighted<int>({{0, 30}, {1, 30}, {2, 40}});
        Rational x;
        if (kind == 0) x = a;
        else if (kind == 1) {
            // One draw a statement: the operands of a product are unsequenced.
            const Rational sign(r.percent(50) ? 1 : -1);
            const int m = static_cast<int>(r.range(1, k));
            x = a + sign * power2(-m);
        } else x = coordinate(r, 0);
        reals.push_back(x);
        s.roots.emplace_back(x, Rational(0));
    }
    if (count) family << "+real" << count;
    // A repeated real factor.
    if (r.percent(30)) {
        const int extra = static_cast<int>(r.range(1, 3)); // multiplicity 2 to 4
        const Rational x = !reals.empty() && r.percent(70) ? reals[r.below(reals.size())] : near[r.below(near.size())].first;
        if (std::find(reals.begin(), reals.end(), x) == reals.end()) s.roots.emplace_back(x, Rational(0));
        for (int c = 0; c < extra; ++c) s.roots.emplace_back(x, Rational(0));
        family << "+repeat(m+" << extra << ")";
    }
    // Ordinary pairs away from the axis.
    if (r.percent(30)) {
        const int ordinary = static_cast<int>(r.range(1, 4));
        for (int j = 0; j < ordinary && s.roots.size() + 2 <= 48; ++j) {
            Complex z{coordinate(r, 0), coordinate(r, 0)};
            if (z.imag == 0) z.imag = Rational(1);
            s.roots.push_back(z);
            s.roots.push_back(conjugate(z));
        }
        family << "+pairs" << ordinary;
    }
    // The transform, then the standing scaling.
    s.axis_transform = r.weighted<int>({{0, 70}, {1, 15}, {2, 15}});
    if (s.axis_transform == 1) {
        const Complex rotation = r.percent(50) ? Complex{Rational(0), Rational(1)} : Complex{Rational(3, 5), Rational(4, 5)};
        for (auto &z : s.roots) z = z * rotation;
        family << (rotation.real == 0 ? "*i" : "*(3+4i)/5");
    } else if (s.axis_transform == 2) {
        s.multiplier = {Rational(r.range(-9, 9)) / Rational(7), Rational(r.range(1, 9)) / Rational(11)};
        family << ";c=gaussian";
    }
    if (r.percent(30)) {
        const int e = static_cast<int>(r.range(-60, 60));
        for (auto &z : s.roots) z = z * Complex(power2(e));
        family << "*2^" << e;
    }
    s.family = family.str();
    // The query.
    const bool real_input = real_coefficients(build(s));
    s.method = r.percent(70) ? CPS_METHOD_AUTO : static_cast<cps_method>(r.range(1, 3));
    const bool explicit_method = s.method != CPS_METHOD_AUTO;
    s.goal = r.weighted<cps_goal>({{CPS_ROOT_ACCURACY, 40}, {CPS_ROOT_APPROXIMATION, 40}, {CPS_CLUSTER_RADIUS, explicit_method ? 0u : 20u}});
    s.cluster_exponent = static_cast<int>(r.range(1, 80));
    if (!explicit_method && real_input && r.percent(40)) {
        s.real_only = true;
        if (s.goal == CPS_CLUSTER_RADIUS) s.goal = r.percent(50) ? CPS_ROOT_ACCURACY : CPS_ROOT_APPROXIMATION;
        if (r.percent(50)) {
            const auto variant = r.weighted<D1Variant>({{D1Variant::Isolation, 40}, {D1Variant::Region, 40}, {D1Variant::Session, 20}});
            draw_d1_query(s, r, variant, D1Refusal::None);
        }
    } else if (!explicit_method && !real_input && r.percent(10)) {
        s.real_only = true; // the input is not real: the query must be refused (axis_refusal_missing otherwise)
        if (s.goal == CPS_CLUSTER_RADIUS) s.goal = CPS_ROOT_ACCURACY;
    }
    return s;
}

// The secular family (--secular): the review's case T17 (codex-notes/review-20261014/adversarial-cases.md;
// docs/adversarial-corpus-note-20261016.md), native secular input against its exact numerator. The equation
// S(x) = sum_k a_k/(x - b_k) - 1 with the planted roots r_i and the poles b_k has the numerator
// sum_k a_k prod_(j != k)(x - b_j) - prod_k (x - b_k) = -P(x), P = prod_i (x - r_i), when
// a_k = -P(b_k)/Q'(b_k), Q = prod_k (x - b_k): the partial fractions of (Q - P)/Q, whose numerator's degree is
// below the poles' count (both monic). So the planted roots are the equation's roots and the standing checker
// applies as it is. A pole at a root makes its residue zero and two equal poles make Q' vanish: the constructor
// must refuse both (the refusal variant). The self-test checks S(r_i) = 0 exactly.
inline Complex complex_divide(const Complex &a, const Complex &b) {
    const Rational norm = b.real * b.real + b.imag * b.imag;
    const Complex product = a * conjugate(b);
    return {product.real / norm, product.imag / norm};
}
inline Complex evaluate_at(const ComplexPolynomial &p, const Complex &z) {
    Complex value;
    for (auto i = p.rbegin(); i != p.rend(); ++i) value = value * z + *i;
    return value;
}
inline std::vector<Complex> secular_residues(const std::vector<Complex> &roots, const std::vector<Complex> &poles) {
    const ComplexPolynomial P = build_monic(roots, Complex(Rational(1)));
    std::vector<Complex> residues;
    for (size_t k = 0; k < poles.size(); ++k) {
        Complex derivative{Rational(1)};
        for (size_t j = 0; j < poles.size(); ++j)
            if (j != k) derivative = derivative * (poles[k] - poles[j]);
        if (zero(derivative)) { // two equal poles: no partial fraction (the constructor refuses them)
            residues.emplace_back(Rational(0));
            continue;
        }
        const Complex value = evaluate_at(P, poles[k]);
        residues.push_back(complex_divide(Complex(-value.real, -value.imag), derivative));
    }
    return residues;
}
inline std::vector<certpolysolve::SecularEquation::Term> secular_terms(const Spec &s) {
    std::vector<certpolysolve::SecularEquation::Term> terms;
    const auto residues = secular_residues(s.roots, s.poles);
    const auto number = [](const Complex &z) {
        return certpolysolve::Number::complex(certpolysolve::Number(z.real.str()), certpolysolve::Number(z.imag.str()));
    };
    for (size_t k = 0; k < s.poles.size(); ++k) terms.push_back({number(residues[k]), number(s.poles[k])});
    return terms;
}
// S(z) exactly (z not a pole): zero at every planted root.
inline Complex secular_value(const Spec &s, const Complex &z) {
    const auto residues = secular_residues(s.roots, s.poles);
    Complex sum{Rational(-1)};
    for (size_t k = 0; k < s.poles.size(); ++k) sum = sum + complex_divide(residues[k], z - s.poles[k]);
    return sum;
}
// - The roots are the standing case's planted roots (a Mignotte seed, which has none, gets a few planted in
//   their place). The native route's crossover at degree 256 (heuristics::native_secular_fallback_degree) is left to a
//   list of its own: a grid of 257 roots at a 128-bit cap ran to the 30 s deadline in the first trial.
// - The poles, one a root: planted-style coordinates (60%), or beside a root at a distance 2^-m, m up to 200
//   (20%); a real numerator takes a real-only query in 10% (the exact numerator fallback); a refusal in 10%:
//   two equal poles, a pole at a root (a zero residue), the empty equation, an explicit classical or piecewise
//   method, or real_only on a numerator that is not real.
// - The query otherwise is the standing case's, with classical and piecewise mapped to the secular method
//   (native input supports Auto and Secular only) and the cluster goal kept for Auto. The numerator has no
//   scalar (the multiplier is one). The standing case's roots, budgets and caps are what the seed would have
//   had, so its records tell the two apart by `family` and the `secular` fields only.
inline constexpr uint64_t secular_stream = UINT64_C(0x7365632d66757a7a); // "sec-fuzz"
inline Spec generate_secular(uint64_t seed) {
    Spec s = generate(seed);
    Rng r(seed ^ secular_stream);
    if (s.oracle != Oracle::Planted) {
        s.oracle = Oracle::Planted;
        const int count = std::max(1, std::min(s.mignotte_n, 12));
        s.mignotte_n = s.mignotte_k = 0;
        for (int j = 0; j < count; ++j) {
            // One draw a statement: the order of a call's arguments is unspecified.
            const Rational re = coordinate(r, 0);
            const Rational im = r.percent(50) ? coordinate(r, 0) : Rational(0);
            s.roots.emplace_back(re, im);
        }
        s.family = "planted" + std::to_string(count);
    }
    s.multiplier = Complex(Rational(1));
    s.secular = 1;
    const size_t n = s.roots.size();
    const auto variant = r.weighted<int>({{0, 60}, {1, 20}, {2, 10}, {3, 10}});
    const auto taken = [&](const Complex &b) {
        const auto same = [&](const Complex &w) { return w.real == b.real && w.imag == b.imag; };
        return std::any_of(s.roots.begin(), s.roots.end(), same) || std::any_of(s.poles.begin(), s.poles.end(), same);
    };
    for (size_t k = 0; k < n; ++k) {
        Complex b;
        do {
            if (variant == 1 && r.percent(50)) {
                const Complex &root = s.roots[r.below(n)];
                const int m = static_cast<int>(r.range(4, 200));
                const Rational re = Rational(r.range(1, 3)) * power2(-m); // one draw a statement (the order of a
                const Rational im = Rational(r.range(-2, 2)) * power2(-m); // call's arguments is unspecified)
                b = root + Complex(re, im);
            } else b = Complex{coordinate(r, 0), r.percent(70) ? coordinate(r, 0) : Rational(0)};
        } while (taken(b));
        s.poles.push_back(b);
    }
    const bool real_numerator = real_coefficients(build(s));
    if (s.method == CPS_METHOD_CLASSICAL || s.method == CPS_METHOD_PIECEWISE) s.method = CPS_METHOD_SECULAR;
    if (s.method == CPS_METHOD_SECULAR && s.goal == CPS_CLUSTER_RADIUS) s.goal = CPS_ROOT_ACCURACY;
    s.real_only = variant == 2 && real_numerator && s.method == CPS_METHOD_AUTO;
    if (s.real_only && s.goal == CPS_CLUSTER_RADIUS) s.goal = r.percent(50) ? CPS_ROOT_ACCURACY : CPS_ROOT_APPROXIMATION;
    std::string name = variant == 1 ? "poles-near" : s.real_only ? "real-only" : "poles";
    auto refusal = SecularRefusal::None;
    if (variant == 3) {
        refusal = r.weighted<SecularRefusal>({{SecularRefusal::DuplicatePole, 1}, {SecularRefusal::ZeroResidue, 1},
                                              {SecularRefusal::Empty, 1}, {SecularRefusal::ExplicitMethod, 1},
                                              {SecularRefusal::RealOnlyNonReal, 1}});
        if (refusal == SecularRefusal::DuplicatePole && n < 2) refusal = SecularRefusal::ZeroResidue;
        if (refusal == SecularRefusal::RealOnlyNonReal && real_numerator) refusal = SecularRefusal::ExplicitMethod;
        switch (refusal) {
        case SecularRefusal::DuplicatePole: s.poles[1] = s.poles[0]; break;
        case SecularRefusal::ZeroResidue: s.poles[0] = s.roots[0]; break;
        case SecularRefusal::Empty:
            s.roots.clear();
            s.poles.clear();
            s.real_only = false;
            break;
        case SecularRefusal::ExplicitMethod:
            s.method = r.percent(50) ? CPS_METHOD_CLASSICAL : CPS_METHOD_PIECEWISE;
            if (s.goal == CPS_CLUSTER_RADIUS) s.goal = CPS_ROOT_ACCURACY;
            s.real_only = false;
            break;
        case SecularRefusal::RealOnlyNonReal:
            s.method = CPS_METHOD_AUTO;
            s.real_only = true;
            if (s.goal == CPS_CLUSTER_RADIUS) s.goal = CPS_ROOT_ACCURACY;
            break;
        case SecularRefusal::None: break;
        }
        name = "refusal " + std::to_string(static_cast<int>(refusal));
    }
    s.secular_refusal = static_cast<int>(refusal);
    s.family += ";secular(" + name + ")";
    return s;
}

// The program family (--program): polynomials given as programs (cps_solve_program; docs/slp-input-spec-20261015.md)
// under the standing case's caps, small work budgets in a share of the cases, and the stops of the family's own
// draws: a cancellation at a drawn checkpoint, a lower cap followed by a cancellation, a time limit of 1 us to 10 ms.
// The shapes (ProgramShape):
// - Product: two or three factors multiplied, each a planted polynomial in Horner form, a product of linear
//   factors, or a power of a small planted factor (the route's product and power lifts); in a share of the factors
//   a root equal to, or 2^-k beside, a root of an earlier factor (shared roots send the output to a whole solve).
// - Generic: one planted polynomial in Horner form (the route without structure), with clusters and repeats.
// - Composition: f(g(x)) with g quadratic, optionally after a linear map; f's roots are g's values at drawn points
//   s, so every root is planted (s and -b/a - s); f in Horner form over the cut node g, or a product of (g - r).
// - Recurrence: Chebyshev T and U, Hermite, Laguerre and Legendre by their three-term recurrences, degree 2 to 32,
//   optionally at a scaled variable: all roots real and simple, so a region's count is Sturm's count on its real
//   interval when it meets the axis and 0 otherwise, decided on the exact expansion (Oracle::Expansion).
// Any shape may be wrapped as (P + Q) - Q with deg Q < deg P (the output is no longer a product, the degree stays
// exact). The query: Auto (70%), classical or secular; the approximation, accuracy or cluster goal; in 8% a query
// the program contract refuses (piecewise, real_only, the isolation goal, a region).
enum class ProgramShape : int { None = 0, Product = 1, Generic = 2, Composition = 3, Recurrence = 4 };
// p divided by its leading coefficient: a recurrence under a Gaussian scalar is real again.
inline ComplexPolynomial monic_of(ComplexPolynomial p) {
    if (p.empty()) return p;
    const Complex leading = p.back();
    for (auto &c : p) c = complex_divide(c, leading);
    return p;
}
inline const char *program_shape_name(int shape) {
    switch (static_cast<ProgramShape>(shape)) {
    case ProgramShape::Product: return "product";
    case ProgramShape::Generic: return "generic";
    case ProgramShape::Composition: return "composition";
    case ProgramShape::Recurrence: return "recurrence";
    case ProgramShape::None: break;
    }
    return "none";
}
enum class ProgramRefusal : int { None = 0, Piecewise = 1, RealOnly = 2, Isolation = 3, Region = 4 };

// Appends program lines and returns node indices. Each call is one statement in the generator, so the lines come
// in the same order from every compiler (the order of a call's arguments is unspecified).
class ProgramWriter {
  public:
    std::vector<std::string> lines;
    uint64_t x() { return emit("x"); }
    uint64_t constant(const Complex &z) { return emit("c " + z.real.str() + " " + z.imag.str()); }
    uint64_t add(uint64_t a, uint64_t b) { return emit("add " + std::to_string(a) + " " + std::to_string(b)); }
    uint64_t sub(uint64_t a, uint64_t b) { return emit("sub " + std::to_string(a) + " " + std::to_string(b)); }
    uint64_t mul(uint64_t a, uint64_t b) { return emit("mul " + std::to_string(a) + " " + std::to_string(b)); }
    uint64_t pow(uint64_t a, uint64_t e) { return emit("pow " + std::to_string(a) + " " + std::to_string(e)); }
    void out(uint64_t node) { lines.push_back("out " + std::to_string(node)); }
    // Horner's rule for the ascending coefficients p at the node `at`.
    uint64_t horner(const ComplexPolynomial &p, uint64_t at) {
        uint64_t value = constant(p.back());
        for (size_t j = p.size() - 1; j-- > 0;) {
            const uint64_t product = mul(value, at);
            const uint64_t coefficient = constant(p[j]);
            value = add(product, coefficient);
        }
        return value;
    }
    // The product of (at - z) over the roots, multiplied in order.
    uint64_t linear_product(const std::vector<Complex> &roots, uint64_t at) {
        uint64_t value = 0;
        for (size_t j = 0; j < roots.size(); ++j) {
            const uint64_t root = constant(roots[j]);
            const uint64_t factor = sub(at, root);
            value = j ? mul(value, factor) : factor;
        }
        return value;
    }

  private:
    uint64_t nodes_ = 0;
    uint64_t emit(std::string line) {
        lines.push_back(std::move(line));
        return nodes_++;
    }
};

// Planted roots of magnitude about 2^scale: 30% real; a symmetric set keeps its non-real roots in conjugate pairs.
inline std::vector<Complex> program_roots(Rng &r, int count, int scale, bool symmetric) {
    std::vector<Complex> roots;
    while (static_cast<int>(roots.size()) < count) {
        const Rational re = coordinate(r, scale);
        const bool real = r.percent(30);
        const Rational im = real ? Rational(0) : coordinate(r, scale);
        if (symmetric && im != 0 && static_cast<int>(roots.size()) + 2 > count) {
            roots.emplace_back(re, Rational(0));
            continue;
        }
        roots.emplace_back(re, im);
        if (symmetric && im != 0) roots.push_back(conjugate(roots.back()));
    }
    // A cluster about the first root (20%): offsets of a few 2^-depth.
    if (count >= 2 && r.percent(20)) {
        const int depth = static_cast<int>(r.range(8, 60));
        const int size = std::min(count, static_cast<int>(r.range(2, 4)));
        for (int j = 1; j < size; ++j) {
            const Rational dr = Rational(r.range(-3, 3)) * power2(-depth);
            const Rational di = Rational(r.range(-3, 3)) * power2(-depth);
            roots[static_cast<size_t>(j)] = roots[0] + Complex(dr, di);
        }
    }
    // A repeated root (15%).
    if (count >= 2 && r.percent(15)) roots[static_cast<size_t>(count) - 1] = roots[0];
    return roots;
}
// A nonzero scalar: a power of two, a rational, or a Gaussian rational.
inline Complex program_scalar(Rng &r) {
    switch (r.below(3)) {
    case 0: return Complex(power2(static_cast<int>(r.range(-40, 40))));
    case 1: {
        const int64_t a = r.range(1, 999);
        const int64_t b = r.range(1, 999);
        const bool negative = r.percent(50);
        return Complex(Rational(Integer(negative ? -a : a)) / Rational(Integer(b)));
    }
    default: {
        const Rational re = Rational(r.range(-9, 9)) / Rational(7);
        const Rational im = Rational(r.range(1, 9)) / Rational(11);
        return {re, im};
    }
    }
}
inline Complex nonzero_coordinate(Rng &r, int scale) {
    const Rational re = coordinate(r, scale);
    const Rational im = r.percent(50) ? coordinate(r, scale) : Rational(0);
    return re == 0 && im == 0 ? Complex(Rational(1)) : Complex(re, im);
}

inline constexpr uint64_t program_stream = UINT64_C(0x70726f672d66757a); // "prog-fuz"
inline Spec generate_program(uint64_t seed) {
    const Spec standing = generate(seed);
    Rng r(seed ^ program_stream);
    Spec s;
    s.seed = seed;
    s.bits = standing.bits;
    s.max_work_bits = standing.max_work_bits;
    s.work_units = standing.work_units;
    s.live_cells = standing.live_cells;
    s.determinism = standing.determinism;
    s.oracle = Oracle::Planted;
    ProgramWriter w;
    const uint64_t x = w.x();
    std::ostringstream family;
    uint64_t output = 0;
    const auto shape = r.weighted<ProgramShape>({{ProgramShape::Product, 45}, {ProgramShape::Generic, 15},
                                                 {ProgramShape::Composition, 20}, {ProgramShape::Recurrence, 20}});
    s.program = static_cast<int>(shape);
    family << "program(" << program_shape_name(s.program);
    switch (shape) {
    case ProgramShape::Product: {
        const int factors = r.weighted<int>({{2, 60}, {3, 40}});
        std::vector<uint64_t> nodes;
        for (int f = 0; f < factors && s.roots.size() < 48; ++f) {
            const int room = 48 - static_cast<int>(s.roots.size());
            const int kind = r.weighted<int>({{0, 50}, {1, 20}, {2, 30}}); // Horner, linear, power
            const int scale = r.percent(30) ? static_cast<int>(r.range(-12, 12)) : 0;
            const bool symmetric = r.percent(40);
            const int exponent = kind == 2 ? static_cast<int>(r.range(2, 5)) : 1;
            const int drawn = kind == 2 ? static_cast<int>(r.range(1, 3)) : static_cast<int>(r.range(1, kind == 1 ? 6 : 12));
            const int count = std::max(1, std::min(drawn, room / exponent));
            if (count * exponent > room) break;
            auto roots = program_roots(r, count, scale, symmetric);
            // A root of an earlier factor, or one 2^-k beside it (35%).
            std::string tag;
            if (!s.roots.empty() && r.percent(35)) {
                const Complex z = s.roots[r.below(s.roots.size())];
                if (r.percent(70)) {
                    roots[0] = z;
                    tag = "+shared";
                } else {
                    const int k = static_cast<int>(r.range(8, 80));
                    roots[0] = z + Complex(power2(-k));
                    tag = "+beside(2^-" + std::to_string(k) + ")";
                }
            }
            const bool scaled = r.percent(30);
            const Complex scalar = scaled ? program_scalar(r) : Complex(Rational(1));
            uint64_t node = 0;
            if (kind == 1) node = w.linear_product(roots, x);
            else {
                const bool horner = kind == 0 || r.percent(50);
                const uint64_t base = horner ? w.horner(build_monic(roots, scalar), x) : w.linear_product(roots, x);
                node = kind == 2 ? w.pow(base, static_cast<uint64_t>(exponent)) : base;
            }
            family << (f ? "," : ":") << (kind == 1 ? "linear" : kind == 2 ? "power" : "horner") << count;
            if (kind == 2) family << "^" << exponent;
            family << tag;
            for (int e = 0; e < exponent; ++e) s.roots.insert(s.roots.end(), roots.begin(), roots.end());
            nodes.push_back(node);
        }
        output = nodes[0];
        for (size_t j = 1; j < nodes.size(); ++j) output = w.mul(output, nodes[j]);
        if (r.percent(25)) {
            const Complex scalar = program_scalar(r);
            const uint64_t c = w.constant(scalar);
            output = w.mul(c, output);
            family << ";scalar";
        }
        break;
    }
    case ProgramShape::Generic: {
        const int count = static_cast<int>(r.range(1, 32));
        const int scale = r.percent(30) ? static_cast<int>(r.range(-12, 12)) : 0;
        const bool symmetric = r.percent(40);
        s.roots = program_roots(r, count, scale, symmetric);
        const bool scaled = r.percent(40);
        const Complex scalar = scaled ? program_scalar(r) : Complex(Rational(1));
        output = w.horner(build_monic(s.roots, scalar), x);
        family << ":" << count;
        break;
    }
    case ProgramShape::Composition: {
        // g(y) = a y^2 + b y + c, after the linear map y = p x + q in 30%.
        const Complex a = nonzero_coordinate(r, 0);
        const Complex b = r.percent(30) ? Complex(Rational(0)) : nonzero_coordinate(r, 0);
        const Complex c = nonzero_coordinate(r, 0);
        const bool mapped = r.percent(30);
        const Complex p = mapped ? nonzero_coordinate(r, 0) : Complex(Rational(1));
        const Complex q = mapped ? nonzero_coordinate(r, 0) : Complex(Rational(0));
        const int count = static_cast<int>(r.range(1, 12));
        const auto seeds = program_roots(r, count, 0, false);
        std::vector<Complex> values; // f's roots: g at the seeds
        const Complex mirror = complex_divide(Complex(-b.real, -b.imag), a); // the two preimages add up to -b/a
        for (const auto &u : seeds) {
            values.push_back((a * u + b) * u + c);
            for (const auto &preimage : {u, mirror - u}) s.roots.push_back(complex_divide(preimage - q, p));
        }
        uint64_t y = x;
        if (mapped) {
            const uint64_t cp = w.constant(p);
            const uint64_t px = w.mul(cp, x);
            const uint64_t cq = w.constant(q);
            y = w.add(px, cq);
        }
        const uint64_t g = w.horner({c, b, a}, y);
        const bool horner = r.percent(60);
        const bool scaled = r.percent(30);
        const Complex scalar = scaled ? program_scalar(r) : Complex(Rational(1));
        output = horner ? w.horner(build_monic(values, scalar), g) : w.linear_product(values, g);
        family << ":" << count << (mapped ? ",mapped" : "") << (horner ? ",horner" : ",linear");
        break;
    }
    case ProgramShape::Recurrence: {
        s.oracle = Oracle::Expansion;
        const int kind = static_cast<int>(r.below(5)); // T, U, Hermite, Laguerre, Legendre
        const int n = static_cast<int>(r.range(2, 32));
        const bool scaled = r.percent(30);
        const int e = scaled ? static_cast<int>(r.range(-6, 6)) : 0;
        uint64_t X = x;
        if (scaled) {
            const uint64_t factor = w.constant(Complex(power2(e)));
            X = w.mul(factor, x);
        }
        const uint64_t one = w.constant(Complex(Rational(1)));
        const uint64_t two = w.constant(Complex(Rational(2)));
        const uint64_t two_x = w.mul(two, X);
        uint64_t previous = one, current = 0;
        switch (kind) {
        case 0: current = X; break;
        case 1: case 2: current = two_x; break;
        case 3: current = w.sub(one, X); break;
        default: current = X; break;
        }
        for (int k = 1; k < n; ++k) {
            uint64_t next = 0;
            if (kind == 0 || kind == 1) {
                const uint64_t t = w.mul(two_x, current);
                next = w.sub(t, previous);
            } else if (kind == 2) {
                const uint64_t t = w.mul(two_x, current);
                const uint64_t ck = w.constant(Complex(Rational(2 * k)));
                const uint64_t u = w.mul(ck, previous);
                next = w.sub(t, u);
            } else {
                // Laguerre ((2k+1 - X) L_k - k L_(k-1))/(k+1); Legendre ((2k+1) X P_k - k P_(k-1))/(k+1).
                const uint64_t odd = w.constant(Complex(Rational(2 * k + 1)));
                const uint64_t a = kind == 3 ? w.sub(odd, X) : w.mul(odd, X);
                const uint64_t t = w.mul(a, current);
                const uint64_t ck = w.constant(Complex(Rational(k)));
                const uint64_t u = w.mul(ck, previous);
                const uint64_t difference = w.sub(t, u);
                const uint64_t inverse = w.constant(Complex(Rational(1, k + 1)));
                next = w.mul(inverse, difference);
            }
            previous = current;
            current = next;
        }
        output = current;
        s.program_degree = n;
        static constexpr std::array<const char *, 5> names{{"chebyshev-t", "chebyshev-u", "hermite", "laguerre", "legendre"}};
        family << ":" << names[static_cast<size_t>(kind)] << n;
        if (scaled) family << ",x*2^" << e;
        if (r.percent(30)) {
            const Complex scalar = program_scalar(r);
            const uint64_t c = w.constant(scalar);
            output = w.mul(output, c);
            family << ";scalar";
        }
        break;
    }
    case ProgramShape::None: break;
    }
    if (s.oracle == Oracle::Planted) s.program_degree = static_cast<int>(s.roots.size());
    // (P + Q) - Q with deg Q < deg P (15%).
    if (s.program_degree >= 1 && r.percent(15)) {
        const int degree = static_cast<int>(r.range(0, s.program_degree - 1));
        ComplexPolynomial q;
        for (int j = 0; j <= degree; ++j) {
            // One draw a statement: the operands of a quotient are unsequenced.
            const Rational re_numerator(r.range(-9, 9));
            const Rational re = re_numerator / Rational(r.range(1, 9));
            Rational im(0);
            if (r.percent(50)) {
                const Rational im_numerator(r.range(-9, 9));
                im = im_numerator / Rational(r.range(1, 9));
            }
            q.emplace_back(re, im);
        }
        if (zero(q.back())) q.back() = Complex(Rational(1));
        const uint64_t noise = w.horner(q, x);
        const uint64_t sum = w.add(output, noise);
        output = w.sub(sum, noise);
        family << ";+q-q(" << degree << ")";
    }
    w.out(output);
    s.program_lines = std::move(w.lines);
    family << ")";
    // The query.
    s.method = r.weighted<cps_method>({{CPS_METHOD_AUTO, 70}, {CPS_METHOD_CLASSICAL, 15}, {CPS_METHOD_SECULAR, 15}});
    s.goal = r.weighted<cps_goal>({{CPS_ROOT_APPROXIMATION, 35}, {CPS_ROOT_ACCURACY, 35}, {CPS_CLUSTER_RADIUS, 30}});
    s.cluster_exponent = static_cast<int>(r.range(1, 80));
    const auto refusal = r.weighted<ProgramRefusal>({{ProgramRefusal::None, 92}, {ProgramRefusal::Piecewise, 2},
                                                    {ProgramRefusal::RealOnly, 2}, {ProgramRefusal::Isolation, 2},
                                                    {ProgramRefusal::Region, 2}});
    switch (refusal) {
    case ProgramRefusal::Piecewise: s.method = CPS_METHOD_PIECEWISE; break;
    case ProgramRefusal::RealOnly: s.real_only = true; break;
    case ProgramRefusal::Isolation: s.goal = CPS_ROOT_ISOLATION; break;
    case ProgramRefusal::Region:
        s.region_lower = "-1";
        s.region_upper = "1";
        break;
    case ProgramRefusal::None: break;
    }
    if (refusal != ProgramRefusal::None) family << ";refusal " << static_cast<int>(refusal);
    // The stops.
    if (r.percent(45)) s.work_units = static_cast<uint64_t>(r.range(1, 40));
    if (r.percent(60)) s.program_cancel = 1 + r.below(UINT64_C(1) << 62);
    if (r.percent(30)) {
        s.program_cap_cancel = r.weighted<uint64_t>({{32, 1}, {48, 1}, {64, 1}, {96, 1}});
        s.program_cap_cancel_draw = 1 + r.below(UINT64_C(1) << 62);
    }
    if (r.percent(25)) s.program_time_limit_us = r.weighted<uint64_t>({{1, 1}, {10, 1}, {100, 1}, {1000, 1}, {10000, 1}});
    s.family = family.str();
    return s;
}
} // namespace fuzz
