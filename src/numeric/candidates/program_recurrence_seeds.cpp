// Structural seeds for a program whose output is a composition f(g(x)) through a cut node or a squaring recurrence
// F(x, q) over an earlier node q (recurrence_seeds, declared in program_recurrence_seeds.hpp): proposals for every root
// from the roots of the smaller programs below, solved roughly first.
//
// recurrence_seeds (seeds_at) declines below a degree threshold, past 64 levels, or where no rule applies; the caller
// then seeds as before. A composition is tried first (squaring_recurrence with no x outside q, then above_cut for f and
// cone for g): composition_seeds solves f and g roughly (rough_roots), starts each g(x) = s, for s a root of f, by a
// Newton step from g's roots (g and g' there from one ProgramMachine), polishes it with the binary64 stage (polish),
// and takes the deg f times deg g roots as the seeds: the output's exact preimages, up to the rough solves' accuracy.
// Otherwise a squaring recurrence (squaring_recurrence): q's roots r from a rough solve of q's cone, q'(r) from the
// binary64 evaluator, and for each r the coefficients of F(r, y) in y (recurrence_coefficients: the nodes below q at
// x = r in binary64, the nodes above it polynomials in y truncated to F's degree d), whose roots y (small_roots:
// directly to degree 2, Aberth's iteration from Cauchy's circle above) give d seeds z = r + y / q'(r): near a root r of
// q, q(z) is about q'(r) (z - r), so the output vanishes where that equals a root of F(r, .). The seeds q's roots leave
// unplaced go on the root bound's circle, where Aberth's sums push them toward the roots no seed holds.
//
// A rough solve (rough_roots) seeds its sub-program by the same rules, or by the circle counts (program_seeds) where
// they do not apply, then runs the binary64 stage (ProgramMachineStage) to its end. Seeds that coincide exactly are
// spread on a small circle about their value (spread_duplicates). The arithmetic is IEEE operations and FLINT's sine
// and cosine only (square_root, unit), for the same bits on every CPU. Proposals only: certification stays with the
// separated-roots pass on the output.
//
// State: none. The report (RecurrenceSeedReport) counts the levels, the compositions, the binary64 sweeps below, and
// the output's seeds placed and left over.
// Entry points: recurrence_seeds, from solver_program.cpp (program_advance: the preliminary phase's start for an output
// of these forms).
// Not here: the detection of the forms (program_structure.cpp), the generic seeds (program_seeds.cpp), the
// certification of the output's roots (SeparatedRootsPass, program_certificate.cpp), the operations of small_roots'
// Aberth step (aberth_step.hpp).
// Design notes: docs/slp-structural-seeding-note-20261018.md (section 3, the seeds; section 9, as built).
#include "numeric/candidates/program_recurrence_seeds.hpp"
#include "numeric/candidates/aberth_step.hpp"
#include "core/heuristics.hpp"
#include "numeric/arithmetic/machine_arithmetic.hpp"
#include "numeric/evaluate/program_machine.hpp"
#include "numeric/candidates/program_candidates.hpp"
#include "numeric/exact/program_preparation.hpp"
#include "numeric/candidates/program_seeds.hpp"
#include "numeric/evaluate/program_source.hpp"
#include "numeric/exact/program_structure.hpp"
#include <algorithm>
#include <cmath>
#include <vector>

namespace certpolysolve::detail {

namespace {
using C = machine::Complex<double>;

double norm(C a) { return a.re * a.re + a.im * a.im; }
bool finite(C a) { return std::isfinite(a.re) && std::isfinite(a.im); }

// The principal square root, from IEEE square roots and basic operations only (the same bits on every CPU).
C square_root(C z) {
    if (z.re == 0 && z.im == 0)
        return {0, 0};
    const double scale = std::max(std::abs(z.re), std::abs(z.im));
    const double re = z.re / scale, im = z.im / scale;
    const double modulus = scale * std::sqrt(re * re + im * im);
    const double t = std::sqrt((modulus + std::abs(z.re)) / 2);
    if (z.re >= 0)
        return {t, z.im / (2 * t)};
    return {std::abs(z.im) / (2 * t), z.im < 0 ? -t : t};
}

// A unit vector at the angle 2 pi (j + 1/2) / m, from FLINT's sin and cos (the same bits on every CPU).
C unit(uint64_t j, uint64_t m) {
    RealBall t, s, c;
    arb_set_ui(t.value(), 2 * j + 1);
    arb_div_ui(t.value(), t.value(), m, 64);
    arb_sin_cos_pi(s.value(), c.value(), t.value(), 64);
    return {arf_get_d(arb_midref(c.value()), ARF_RND_NEAR), arf_get_d(arb_midref(s.value()), ARF_RND_NEAR)};
}

// The roots of c_0 + c_1 y + ... + c_d y^d (d >= 1, c_d nonzero), proposals in binary64: the stable quadratic
// formula for d = 2, Aberth's iteration from a circle otherwise. Empty when the iteration does not settle.
std::vector<C> small_roots(std::vector<C> c) {
    while (c.size() > 1 && c.back().re == 0 && c.back().im == 0)
        c.pop_back();
    const size_t d = c.size() - 1;
    if (d == 0)
        return {};
    if (d == 1)
        return {-c[0] / c[1]};
    if (d == 2) {
        C s = square_root(c[1] * c[1] - C{4, 0} * c[2] * c[0]);
        if (c[1].re * s.re + c[1].im * s.im < 0)
            s = -s;
        const C q = -(c[1] + s) / C{2, 0};
        if (q.re == 0 && q.im == 0)
            return {C{0, 0}, C{0, 0}};
        return {q / c[2], c[0] / q};
    }
    // Aberth from Cauchy's circle, 1 + max |c_j / c_d|, which holds every root (no libm power: the same bits on
    // every CPU).
    double radius = 0;
    for (size_t j = 0; j < d; ++j)
        radius = std::max(radius, std::sqrt(norm(c[j] / c[d])));
    radius += 1;
    std::vector<C> y(d);
    for (size_t j = 0; j < d; ++j)
        y[j] = unit(j, d) * C{radius, 0};
    for (int sweep = 0; sweep < 200; ++sweep) {
        double largest = 0;
        for (size_t i = 0; i < d; ++i) {
            C value = c[d], derivative{0, 0};
            for (size_t j = d; j-- > 0;) {
                derivative = derivative * y[i] + value;
                value = value * y[i] + c[j];
            }
            if (norm(derivative) == 0)
                continue;
            const C ratio = value / derivative;
            C sum{0, 0};
            for (size_t j = 0; j < d; ++j)
                if (j != i && norm(y[i] - y[j]) > 0)
                    sum = sum + machine::interaction_term(y[i] - y[j]);
            const C step = machine::aberth_correction(ratio, sum);
            if (!finite(step))
                continue;
            y[i] = y[i] - step;
            largest = std::max(largest, norm(step) / std::max(norm(y[i]), 1e-300));
        }
        if (largest < 0x1p-100)
            return y;
    }
    return y; // proposals still: the seeds are never trusted
}

// Seeds that coincide exactly (a multiple root of F(r, .), as y = 0 twice for q^2 (x + 3), or a multiple root of f in
// a composition) would hold each other in Aberth's sums: each such group is spread on a circle of radius
// 2^-20 (1 + |z|) about its value (IEEE scaling and FLINT's sin and cos only).
void spread_duplicates(std::vector<C> &seeds) {
    std::vector<size_t> order(seeds.size());
    for (size_t i = 0; i < order.size(); ++i)
        order[i] = i;
    std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
        return seeds[a].re != seeds[b].re ? seeds[a].re < seeds[b].re : seeds[a].im < seeds[b].im;
    });
    for (size_t a = 0; a < order.size();) {
        size_t b = a + 1;
        while (b < order.size() && seeds[order[b]].re == seeds[order[a]].re && seeds[order[b]].im == seeds[order[a]].im)
            ++b;
        if (b - a > 1) {
            const C z = seeds[order[a]];
            const C radius{std::ldexp(1 + std::sqrt(norm(z)), -20), 0};
            for (size_t j = 0; j < b - a; ++j)
                seeds[order[a + j]] = z + unit(j, b - a) * radius;
        }
        a = b;
    }
}

// The exact dyadic a double is.
Number exact(double x) {
    Float a;
    arf_set_d(a.value(), x);
    Rational q;
    arf_get_fmpq(q.value(), a.value());
    char *text = fmpq_get_str(nullptr, 10, q.value());
    Number n{std::string_view(text)};
    flint_free(text);
    return n;
}

// The output's cone above q, as F(x, q): which nodes read q, directly or not.
std::vector<uint8_t> above_base(const Program &program, uint64_t base) {
    const uint64_t output = *program.output;
    std::vector<uint8_t> above(output + 1, 0);
    above[base] = 1;
    for (uint64_t i = base + 1; i <= output; ++i) {
        const auto &n = program.nodes[i];
        switch (n.op) {
        case ProgramOp::Variable:
        case ProgramOp::Constant: break;
        case ProgramOp::Negate:
        case ProgramOp::Power: above[i] = above[n.a]; break;
        case ProgramOp::Add:
        case ProgramOp::Subtract:
        case ProgramOp::Multiply: above[i] = above[n.a] || above[n.b]; break;
        }
    }
    return above;
}

// F(r, .)'s coefficients: the nodes below q at x = r in binary64, q itself the variable y, the nodes above q
// polynomials in y of degree at most d. Empty when a value is not finite.
std::vector<C> recurrence_coefficients(const Program &program, const std::vector<C> &constants,
                                       const std::vector<uint8_t> &above, uint64_t base, size_t d, C r) {
    const uint64_t output = *program.output;
    std::vector<C> value(output + 1);
    std::vector<std::vector<C>> poly(output + 1);
    const auto truncated = [&](std::vector<C> p) {
        if (p.size() > d + 1)
            p.resize(d + 1);
        return p;
    };
    const auto as_poly = [&](uint64_t i) { return above[i] ? poly[i] : std::vector<C>{value[i]}; };
    const auto add = [&](const std::vector<C> &a, const std::vector<C> &b, bool subtract) {
        std::vector<C> out(std::max(a.size(), b.size()), C{0, 0});
        for (size_t j = 0; j < a.size(); ++j)
            out[j] = a[j];
        for (size_t j = 0; j < b.size(); ++j)
            out[j] = subtract ? out[j] - b[j] : out[j] + b[j];
        return out;
    };
    const auto multiply = [&](const std::vector<C> &a, const std::vector<C> &b) {
        std::vector<C> out(std::min(a.size() + b.size() - 1, d + 1), C{0, 0});
        for (size_t i = 0; i < a.size(); ++i)
            for (size_t j = 0; j < b.size() && i + j < out.size(); ++j)
                out[i + j] = out[i + j] + a[i] * b[j];
        return out;
    };
    for (uint64_t i = 0; i <= output; ++i) {
        const auto &n = program.nodes[i];
        if (i == base) {
            poly[i] = {C{0, 0}, C{1, 0}};
            continue;
        }
        if (!above[i]) {
            switch (n.op) {
            case ProgramOp::Variable: value[i] = r; break;
            case ProgramOp::Constant: value[i] = constants[n.value]; break;
            case ProgramOp::Negate: value[i] = -value[n.a]; break;
            case ProgramOp::Add: value[i] = value[n.a] + value[n.b]; break;
            case ProgramOp::Subtract: value[i] = value[n.a] - value[n.b]; break;
            case ProgramOp::Multiply: value[i] = value[n.a] * value[n.b]; break;
            case ProgramOp::Power: {
                C result{1, 0}, square = value[n.a];
                for (uint64_t e = n.value; e; e >>= 1) {
                    if (e & 1)
                        result = result * square;
                    if (e > 1)
                        square = square * square;
                }
                value[i] = result;
                break;
            }
            }
            continue;
        }
        switch (n.op) {
        case ProgramOp::Variable:
        case ProgramOp::Constant: break;
        case ProgramOp::Negate: {
            auto p = poly[n.a];
            for (auto &t : p)
                t = -t;
            poly[i] = std::move(p);
            break;
        }
        case ProgramOp::Add: poly[i] = truncated(add(as_poly(n.a), as_poly(n.b), false)); break;
        case ProgramOp::Subtract: poly[i] = truncated(add(as_poly(n.a), as_poly(n.b), true)); break;
        case ProgramOp::Multiply: poly[i] = multiply(as_poly(n.a), as_poly(n.b)); break;
        case ProgramOp::Power: {
            std::vector<C> result{C{1, 0}}, square = poly[n.a];
            for (uint64_t e = n.value; e; e >>= 1) {
                if (e & 1)
                    result = multiply(result, square);
                if (e > 1)
                    square = multiply(square, square);
            }
            poly[i] = std::move(result);
            break;
        }
        }
    }
    auto out = poly[output];
    for (const auto &t : out)
        if (!finite(t))
            return {};
    return out;
}

bool seeds_at(acb_ptr out, const Program &program, uint64_t degree, const mag_t root_bound, slong precision,
              const cps_options &options, certpolysolve::detail::CandidateStatistics &statistics, RecurrenceSeedReport &report,
              uint64_t threshold, unsigned depth);

// The binary64 stage from the given seeds to its end: proposals for the program's `degree` roots. False when the
// stage refuses, or on cancellation.
bool polish(std::vector<C> &roots, const Program &program, uint64_t degree, acb_srcptr seeds, const cps_options &options,
            certpolysolve::detail::CandidateStatistics &statistics, RecurrenceSeedReport &report) {
    ProgramMachineStage stage(program, seeds, degree, heuristics::program_machine_chunk, options.temporary_bytes);
    if (!stage.valid())
        return false;
    for (unsigned packet = 0; packet < 4096; ++packet) {
        const auto step = stage.advance(heuristics::program_machine_packet_sweeps, options, statistics);
        report.sub_sweeps += step.sweeps;
        if (step.cancelled)
            return false;
        if (step.finished)
            break;
    }
    const auto n = static_cast<slong>(degree);
    RootVector points(n);
    stage.export_points(points.data());
    roots.resize(degree);
    for (uint64_t i = 0; i < degree; ++i) {
        const acb_srcptr z = points.data() + static_cast<slong>(i);
        roots[i] = {arf_get_d(arb_midref(acb_realref(z)), ARF_RND_NEAR), arf_get_d(arb_midref(acb_imagref(z)), ARF_RND_NEAR)};
    }
    return true;
}

// Proposals for the roots of a prepared sub-program: the rule again where it applies, the circle counts otherwise,
// then the binary64 stage to its end. False when the stage refuses, or on cancellation.
bool rough_roots(std::vector<C> &roots, const Program &program, const ProgramPreparation &preparation,
                 slong precision, const cps_options &options, certpolysolve::detail::CandidateStatistics &statistics,
                 RecurrenceSeedReport &report, uint64_t threshold, unsigned depth) {
    const uint64_t degree = preparation.degree;
    if (options.cancelled())
        return false;
    const auto n = static_cast<slong>(degree);
    RootVector seeds(n);
    if (!seeds_at(seeds.data(), program, degree, preparation.root_bound.value(), precision, options, statistics, report,
                  threshold, depth)) {
        ProgramSource source(program, n, preparation.leading);
        if (!program_seeds(seeds.data(), source, n, preparation.root_bound.value(), precision, options, statistics,
                           nullptr))
            return false;
    }
    return polish(roots, program, degree, seeds.data(), options, statistics, report);
}

// A composition f(g(x)) through the cut node `base`: f's roots s, rough, then the roots of g(x) - s for each, rough;
// together the output's degree of seeds, its exact preimages up to the rough solves' accuracy.
bool composition_seeds(acb_ptr out, const Program &program, uint64_t degree, uint64_t base, slong precision,
                       const cps_options &options, certpolysolve::detail::CandidateStatistics &statistics, RecurrenceSeedReport &report,
                       uint64_t threshold, unsigned depth) {
    const auto f = above_cut(program, base);
    if (!f)
        return false;
    const auto g = cone(program, base);
    const auto f_preparation = prepare_program(*f, options);
    const auto g_preparation = prepare_program(g, options);
    if (!f_preparation.ready() || !g_preparation.ready() || f_preparation.degree == 0 || g_preparation.degree == 0 ||
        f_preparation.degree * g_preparation.degree != degree)
        return false;
    std::vector<C> values;
    if (!rough_roots(values, *f, f_preparation, precision, options, statistics, report, threshold, depth + 1))
        return false;
    // g's roots and g' there, once; each g - s is seeded from them by a Newton step toward g = s, r + (s - g(r)) /
    // g'(r), and polished by the binary64 stage (the circle counts for each of the deg f solves cost far more).
    // History: docs/slp-structural-seeding-note-20261018.md, section 12.
    const uint64_t m = g_preparation.degree;
    std::vector<C> g_roots;
    if (!rough_roots(g_roots, g, g_preparation, precision, options, statistics, report, threshold, depth + 1))
        return false;
    ProgramMachine<double> g_machine(g);
    if (!g_machine.valid())
        return false;
    std::vector<C> g_values(m), g_slopes(m);
    for (uint64_t i = 0; i < m; ++i) {
        const auto e = g_machine.evaluate_one(g_roots[i]);
        ++statistics.evaluations;
        statistics.evaluation_terms += g.nodes.size();
        if (!e.finite || norm(e.derivative) == 0)
            return false;
        g_values[i] = e.value;
        g_slopes[i] = e.derivative;
    }
    std::vector<C> seeds;
    seeds.reserve(degree);
    RootVector start(static_cast<slong>(m));
    for (const auto &value : values) {
        if (!finite(value) || options.cancelled())
            return false;
        Program shifted = g;
        Coefficient constant;
        constant.real = exact(value.re);
        constant.imag = exact(value.im);
        ProgramNode subtract;
        subtract.op = ProgramOp::Subtract;
        subtract.a = *shifted.output;
        subtract.b = shifted.append_constant(std::move(constant));
        shifted.output = shifted.append(subtract);
        for (uint64_t i = 0; i < m; ++i) {
            const C z = g_roots[i] + (value - g_values[i]) / g_slopes[i];
            if (!finite(z))
                return false;
            acb_set_d_d(start.data() + static_cast<slong>(i), z.re, z.im);
        }
        std::vector<C> roots;
        if (!polish(roots, shifted, m, start.data(), options, statistics, report))
            return false;
        seeds.insert(seeds.end(), roots.begin(), roots.end());
    }
    if (seeds.size() != degree)
        return false;
    spread_duplicates(seeds);
    for (uint64_t i = 0; i < degree; ++i)
        acb_set_d_d(out + static_cast<slong>(i), seeds[i].re, seeds[i].im);
    ++report.levels;
    ++report.compositions;
    if (depth == 0) {
        report.base_degree = g_preparation.degree;
        report.seeded = degree;
        report.leftover = 0;
    }
    return true;
}

bool seeds_at(acb_ptr out, const Program &program, uint64_t degree, const mag_t root_bound, slong precision,
              const cps_options &options, certpolysolve::detail::CandidateStatistics &statistics, RecurrenceSeedReport &report,
              uint64_t threshold, unsigned depth) {
    if (degree <= threshold || depth >= 64 || !program.output)
        return false;
    // A composition first: no x outside the base (x degree 0), any degree in it.
    if (const auto composition = squaring_recurrence(program, 0, UINT64_MAX)) {
        if (composition_seeds(out, program, degree, composition->base, precision, options, statistics, report,
                              threshold, depth))
            return true;
        if (options.cancelled())
            return false;
    }
    const auto recurrence = squaring_recurrence(program);
    if (!recurrence)
        return false;
    const uint64_t d = recurrence->degree_in_base;
    const auto q = cone(program, recurrence->base);
    // q's degree proved (the bounds can run far above it: the recorder's x 0 + 1 has bound 1).
    const auto q_preparation = prepare_program(q, options);
    if (!q_preparation.ready() || q_preparation.degree == 0 || d * q_preparation.degree > degree)
        return false;
    const uint64_t q_degree = q_preparation.degree;
    std::vector<C> roots;
    if (!rough_roots(roots, q, q_preparation, precision, options, statistics, report, threshold, depth + 1))
        return false;
    // q' at q's roots, and the constants in binary64.
    ProgramMachine<double> q_machine(q);
    if (!q_machine.valid())
        return false;
    std::vector<C> constants(program.constants.size());
    ComplexBall ball;
    for (size_t i = 0; i < constants.size(); ++i) {
        program.constants[i].enclose(ball.value(), 64);
        constants[i] = {arf_get_d(arb_midref(acb_realref(ball.value())), ARF_RND_NEAR),
                        arf_get_d(arb_midref(acb_imagref(ball.value())), ARF_RND_NEAR)};
    }
    const auto above = above_base(program, recurrence->base);
    std::vector<C> seeds;
    seeds.reserve(degree);
    for (const auto &r : roots) {
        const auto e = q_machine.evaluate_one(r);
        ++statistics.evaluations;
        statistics.evaluation_terms += q.nodes.size();
        if (!e.finite || norm(e.derivative) == 0)
            continue;
        const auto f = recurrence_coefficients(program, constants, above, recurrence->base, d, r);
        if (f.empty())
            continue;
        for (const auto &y : small_roots(f)) {
            const C z = r + y / e.derivative;
            if (finite(z) && seeds.size() < degree)
                seeds.push_back(z);
        }
    }
    const uint64_t seeded = seeds.size();
    // The leftover seeds on the root bound's circle.
    const double bound = mag_get_d(root_bound);
    const uint64_t leftover = degree - seeded;
    for (uint64_t j = 0; j < leftover; ++j)
        seeds.push_back(unit(j, leftover) * C{bound, 0});
    spread_duplicates(seeds);
    for (uint64_t i = 0; i < degree; ++i)
        acb_set_d_d(out + static_cast<slong>(i), seeds[i].re, seeds[i].im);
    ++report.levels;
    if (depth == 0) {
        report.base_degree = q_degree;
        report.seeded = seeded;
        report.leftover = leftover;
    }
    return true;
}
} // namespace

bool recurrence_seeds(acb_ptr out, const Program &program, uint64_t degree, const mag_t root_bound, slong precision,
                      const cps_options &options, certpolysolve::detail::CandidateStatistics &statistics, RecurrenceSeedReport &report,
                      uint64_t threshold) {
    return seeds_at(out, program, degree, root_bound, precision, options, statistics, report, threshold, 0);
}

} // namespace certpolysolve::detail
