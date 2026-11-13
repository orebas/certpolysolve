#include "numeric/evaluate/structured_source.hpp"
#include "solve/policies.hpp"
#include "solve/local_model_solver.hpp"
#include "certify/power_cover.hpp"
#include "numeric/arithmetic/polynomial_arithmetic.hpp"
#include "internal.hpp"
#include "exact_winding_api.hpp"
#include "api_test_support.hpp"
#include "methods/secular/secular_refinement.hpp"
#include "numeric/candidates/interaction_tree.hpp"
#include "numeric/exact/power_composition.hpp"
#include "numeric/exact/centering.hpp"
#include <numbers>
#include <iostream>

using namespace reference;
namespace d = certpolysolve::detail;
// A solver snapshot (the internal result) as the header's Result.
static certpolysolve::Result wrap(cps_result snapshot) {
    return d::make_result(CPS_OK, d::published_result(std::move(snapshot)));
}
static Complex quotient(const Complex &a, const Complex &b) {
    const Rational norm = b.real * b.real + b.imag * b.imag;
    return {(a.real * b.real + a.imag * b.imag) / norm, (a.imag * b.real - a.real * b.imag) / norm};
}
static void enclose(acb_t out, const Complex &z, slong bits) {
    d::Number(z.real.str()).enclose(acb_realref(out), bits);
    d::Number(z.imag.str()).enclose(acb_imagref(out), bits);
}
static bool contains(const acb_t value, const Complex &z) {
    d::Rational re, im;
    fmpq_set_str(re.value(), z.real.str().c_str(), 10);
    fmpq_set_str(im.value(), z.imag.str().c_str(), 10);
    return arb_contains_fmpq(acb_realref(value), re.value()) && arb_contains_fmpq(acb_imagref(value), im.value());
}
static void enclosure(acb_poly_t out, const ComplexPolynomial &p, slong bits) {
    d::ComplexBall value;
    acb_poly_zero(out);
    for (size_t i = 0; i < p.size(); ++i) {
        enclose(value.value(), p[i], bits);
        acb_poly_set_coeff_acb(out, static_cast<slong>(i), value.value());
    }
}
static std::vector<Complex> weights(const ComplexPolynomial &p, const std::vector<Complex> &nodes) {
    std::vector<Complex> result;
    for (size_t i = 0; i < nodes.size(); ++i) {
        Complex product(1);
        for (size_t j = 0; j < nodes.size(); ++j)
            if (i != j)
                product = product * (nodes[i] - nodes[j]);
        result.push_back(quotient(Complex(0) - evaluate(p, nodes[i]), product * p.back()));
    }
    return result;
}
static void known(const certpolysolve::Result &result, const std::vector<Complex> &roots, int bits) {
    if (!result.goal_met() || !result.complete())
        throw std::runtime_error("structured known-root solve incomplete: " + result.message());
    check_disjoint(result.roots());
    size_t total = 0;
    for (const auto &root : result.roots()) {
        size_t count = 0;
        for (const auto &z : roots)
            count += endpoint(root.box.real.lo) <= z.real && z.real <= endpoint(root.box.real.hi) &&
                     endpoint(root.box.imag.lo) <= z.imag && z.imag <= endpoint(root.box.imag.hi);
        require(count == root.count && count == root.multiplicity.value_or(0), "structured count disagreement");
        total += count;
        check_accuracy(root, bits);
    }
    require(total == roots.size(), "structured cover lost roots");
}
static void conjugate_cover_realness() {
    // These are exact polynomials and independently computed residues, not
    // synthetic disks. A wide singleton crosses the real axis, but contains a
    // nonreal root whose conjugate is in another member of the complete cover.
    for (bool centered : {false, true})
        for (int exponent : {0, -80})
            for (bool perturb : {false, true})
                for (bool reverse : {false, true})
                    for (uint64_t multiplicity : {1u, 3u}) {
                        const Rational scale = power2(exponent);
                        std::vector<Complex> roots{Complex(10, -scale), Complex(10, scale)};
                        std::vector<Complex> exact_nodes{
                            roots[0] + Complex(perturb ? scale * power2(-20) : Rational(0)),
                            Complex(Rational(10) + scale / 2, 7 * scale / 8)};
                        if (centered) {
                            roots.emplace_back(-20);
                            exact_nodes.emplace_back(-20);
                        }
                        if (reverse)
                            std::reverse(exact_nodes.begin(), exact_nodes.end());
                        const auto exact_weights = weights(from_roots(roots), exact_nodes);
                        std::vector<d::ComplexBall> nodes(roots.size()), residues(roots.size());
                        constexpr slong precision = 256;
                        for (size_t i = 0; i < nodes.size(); ++i) {
                            enclose(nodes[i].value(), exact_nodes[i], precision);
                            enclose(residues[i].value(), exact_weights[i], precision);
                        }
                        const auto saved_nodes = nodes, saved_residues = residues;
                        cps_options options;
                        options.goal = CPS_ROOT_APPROXIMATION;
                        std::vector<d::Root> cover;
                        std::vector<uint8_t> inactive;
                        require(d::secular_inclusion_cover(cover, inactive, nodes, residues,
                                                         true, multiplicity, precision, options),
                                "exact conjugate-cover identity failed");
                        require(cover.size() == roots.size(), "conjugate singleton cover was merged");
                        for (size_t i = 0; i < nodes.size(); ++i)
                            require(acb_equal(nodes[i].value(), saved_nodes[i].value()) &&
                                        acb_equal(residues[i].value(), saved_residues[i].value()),
                                    "conjugation proof modified its input");
                        cps_result snapshot;
                        snapshot.roots = cover;
                        const auto exported = wrap(std::move(snapshot)).roots();
                        check_disjoint(exported);
                        size_t total = 0;
                        for (size_t i = 0; i < exported.size(); ++i) {
                            const auto &box = exported[i];
                            size_t count = 0;
                            for (const auto &z : roots)
                                if (endpoint(box.box.real.lo) <= z.real && z.real <= endpoint(box.box.real.hi) &&
                                    endpoint(box.box.imag.lo) <= z.imag && z.imag <= endpoint(box.box.imag.hi)) {
                                    ++count;
                                    require(contains(cover[i].box.value(), z), "internal cover lost an exact root");
                                    require(box.realness != certpolysolve::Realness::real || z.imag == 0,
                                            "conjugation proof called a nonreal root real");
                                }
                            require(count == 1 && box.count == multiplicity && box.multiplicity.value_or(0) == multiplicity,
                                    "secular conjugation cover lost a known root");
                            total += count;
                        }
                        require(total == roots.size(), "conjugate-cover degree accounting failed");
                    }
}
static void centered_zero_initialization() {
    for (unsigned zeros : {1u, 3u, 5u}) {
        std::vector<Complex> roots(zeros, Complex(0));
        if (zeros != 5) {
            roots.emplace_back(2);
            roots.emplace_back(-3, 1);
        }
        d::ComplexPolynomial polynomial;
        enclosure(polynomial.value(), from_roots(roots), 128);
        d::RootVector seeds(static_cast<slong>(roots.size()));
        require(d::initial_candidates(seeds.data(), polynomial.value(), 128),
                "zero root prevented coefficient-scale initialization");
        unsigned zero_seeds = 0;
        for (slong j = 0; j < seeds.size(); ++j) {
            require(acb_is_finite(seeds[j]), "nonfinite zero-factor seed");
            zero_seeds += acb_is_zero(seeds[j]);
        }
        require(zero_seeds == zeros, "initialization lost zero-factor degree");
    }
    // The eight-bit center of these 64 exact roots is 160, itself a root.
    // The persistent source must still initialize and certify all ORIGINAL roots.
    std::vector<Complex> roots;
    for (int j = 129; j <= 192; ++j) roots.emplace_back(j);
    auto polynomial = input(from_roots(roots));
    auto policy = d::production_policy();
    policy.study = d::StudyMethod::Secular;
    policy.refinement = d::RefinementMethod::Secular;
    policy.model_stopping = d::SecularModelPolicy::Continuation;
    policy.regeneration = {.value_only = true, .reuse_products = true};
    policy.secular_preliminary = true;
    policy.model_driven_refresh = true;
    policy.study_statistics = std::make_shared<d::StudyStatistics>();
    d::SolveState state(certpolysolve::detail::dense(*polynomial.get()), policy);
    cps_options options;
    options.bits = 53;
    options.goal = CPS_ROOT_APPROXIMATION;
    options.max_work_bits = 4096;
    const auto result = state.advance(&options);
    require(result.statistics_v2.taylor_shifts != 0, "zero-centered regression did not center");
    require(policy.study_statistics->preliminary_packets != 0,
            "centered zero rejected the persistent secular driver");
    known(wrap(result), roots, 53);
}
// The preliminary stage's validation is FLINT's, by ball Horner: at a candidate of argument t it loses
// degree log2(|cos t| + |sin t|) bits (ball_horner_loss_bits). Where that is the whole working precision the
// validation is not run; where it is not, it runs as before. 1 + x + ... + x^240 has its roots at every argument
// of the unit circle (120 bits lost of 85); the same of degree 8 loses at most four.
static void hopeless_preliminary_validation() {
    // The loss itself: none on the axes, half a bit a term on a diagonal, the worst point decides.
    d::RootVector points(3);
    acb_set_d_d(points[0], 3, 0);
    acb_set_d_d(points[1], 0, -0.25);
    acb_set_d_d(points[2], 1e300, 1e-300);
    require(d::ball_horner_loss_bits(points.data(), 3, 1000) < 1e-6, "points on the axes lose bits");
    acb_set_d_d(points[1], -0x1p-700, 0x1p-700);
    const double diagonal = d::ball_horner_loss_bits(points.data(), 3, 1000);
    require(diagonal > 499.9 && diagonal < 500.1, "a point on a diagonal does not lose half a bit a term");
    for (const slong degree : {240, 8}) {
        const ComplexPolynomial coefficients(static_cast<size_t>(degree) + 1, Complex(1));
        auto polynomial = input(coefficients);
        auto policy = d::production_policy();
        policy.study = d::StudyMethod::Secular;
        policy.refinement = d::RefinementMethod::Secular;
        policy.model_stopping = d::SecularModelPolicy::Continuation;
        policy.secular_preliminary = true;
        policy.regeneration = {.value_only = true, .reuse_products = true};
        policy.study_statistics = std::make_shared<d::StudyStatistics>();
        d::SolveState state(certpolysolve::detail::dense(*polynomial.get()), policy);
        cps_options options;
        options.goal = CPS_ROOT_APPROXIMATION;
        options.max_work_bits = 4096;
        const auto result = state.advance(&options);
        require(result.complete && result.goal_met, "the all-ones polynomial was not solved");
        require(policy.study_statistics->preliminary_packets == 1, "the all-ones polynomial skipped its preliminary phase");
        require(policy.study_statistics->preliminary_validations_skipped == (degree == 240 ? 1u : 0u),
                degree == 240 ? "a validation that loses 120 of 85 bits was run"
                              : "a validation that loses four bits was not run");
    }
}
static void power_composition() {
    for (const auto &base :
         {Complex(1), Complex(1, 1), Complex(1, 2), Complex(Rational(1) / Integer(1 << 20))}) {
        const Complex i(0, 1);
        std::vector<Complex> roots;
        Complex root = base;
        for (int j = 0; j < 4; ++j) {
            roots.push_back(root);
            root = root * i;
        }
        // Distinct squarefree factors, known repeated roots and an exact zero
        // exercise the transformed count ledger and cross-factor cover union.
        auto repeated = roots;
        roots.insert(roots.end(), repeated.begin(), repeated.end());
        roots.emplace_back(0);
        roots.emplace_back(3);
        auto polynomial = input(from_roots(roots));
        // One solve with the 2,000 units the one-unit calls of a session used to share, then a fresh one at 160
        // bits (one advance a solve state since 18 October 2026).
        d::EnginePolicy policy;
        d::SolveState state(certpolysolve::detail::dense(*polynomial.get()), policy);
        cps_options options;
        options.bits = 53;
        options.absolute_tolerance = d::Number("1e-30");
        options.max_work_bits = 4096;
        options.work_units = 2000;
        cps_result result = state.advance(&options);
        if (!result.goal_met)
            std::cerr << "power initial: " << base.real << "," << base.imag << "\n";
        known(wrap(result), roots, 53);
        require(std::any_of(result.roots.begin(), result.roots.end(),
                            [](const d::Root &r) {
                                return r.certificate.find("exact-power-composition") != std::string::npos;
                            }),
                "power identity was not used");
        options.bits = 160;
        options.absolute_tolerance = d::Number("1e-80");
        options.work_units = 4000;
        d::SolveState tighter(certpolysolve::detail::dense(*polynomial.get()), policy);
        result = tighter.advance(&options);
        if (!result.goal_met) {
            std::cerr << "power tightened: " << base.real << "," << base.imag << " bits " << result.work_bits
                      << " roots " << result.roots.size() << "\n";
            for (const auto &r : result.roots)
                std::cerr << r.count << " " << r.certificate << " acc " << acb_rel_accuracy_bits(r.box.value())
                          << "\n";
        }
        known(wrap(result), roots, 160);
    }
    // A contour crossing the branch point is inconclusive, not a count proof.
    std::vector<d::Root> invalid(1), out;
    arb_add_error_2exp_si(acb_realref(invalid[0].box.value()), -10);
    require(!d::lift_power_cover(out, invalid, 4, 128), "ramified box accepted");
}
// The value and derivative at z through the source's evaluate, enclosed or not (what MonicSource::joint gave
// until W8's S2 removed it: no library code called it).
struct Joint {
    d::ComplexBall value, derivative;
    bool enclosed = false;
};
static Joint joint(d::MonicSource &source, const d::ComplexBall &z, slong precision) {
    Joint result;
    result.enclosed = source.evaluate(result.value.value(), result.derivative.value(), z.value(), precision);
    return result;
}
static void exact_local_evaluation() {
    const auto exact = from_roots({Complex(1), Complex(Rational(1) + power2(-90)), Complex(-2, 3)});
    d::GaussianPolynomial gaussian(exact.size());
    for (size_t j = 0; j < exact.size(); ++j) {
        fmpq_set_str(gaussian[j].real.value(), exact[j].real.str().c_str(), 10);
        fmpq_set_str(gaussian[j].imag.value(), exact[j].imag.str().c_str(), 10);
    }
    d::PolynomialSourceOptions options;
    options.input = d::InputSemantics::Exact;
    unsigned translations = 0;
    options.exact_taylor = [gaussian, &translations](acb_poly_t out, const acb_t center, slong bits) {
        d::GaussianRational shift;
        arf_get_fmpq(shift.real.value(), arb_midref(acb_realref(center)));
        arf_get_fmpq(shift.imag.value(), arb_midref(acb_imagref(center)));
        d::GaussianPolynomial shifted;
        if (!d::center_gaussian(shifted, gaussian, shift, 8000000))
            return false;
        d::enclose_gaussian(out, shifted, bits);
        ++translations;
        return true;
    };
    d::PolynomialSource source([exact](acb_poly_t out, slong bits) { enclosure(out, exact, bits); }, 3, false,
                               options);
    ComplexPolynomial derivative(3);
    for (size_t j = 1; j < exact.size(); ++j)
        derivative[j - 1] = exact[j] * Complex(Rational(j));
    for (int j = 0; j < 12; ++j) {
        d::ComplexBall z;
        const auto center = Complex(Rational(1) + power2(-95) * j, power2(-96));
        enclose(z.value(), center, 256);
        arb_add_error_2exp_si(acb_realref(z.value()), -100);
        arb_add_error_2exp_si(acb_imagref(z.value()), -100);
        const auto result = joint(source, z, 192);
        const auto value = source.value(z, 192);
        require(result.enclosed, "exact local model failed");
        require(value.status == d::EvaluationStatus::Enclosed, "value-only exact translation failed");
        for (int x : {-1, 0, 1})
            for (int y : {-1, 0, 1}) {
                const auto sample = center + Complex(power2(-100) * x, power2(-100) * y);
                require(contains(result.value.value(), evaluate(exact, sample)) &&
                            contains(result.derivative.value(), evaluate(derivative, sample)),
                        "exact local value/derivative lost a whole-ball sample");
                require(contains(value.value.value(), evaluate(exact, sample)),
                        "value-only exact translation lost a whole-ball sample");
            }
    }
    require(translations == 1, "nearby queries did not reuse the exact translation");
}
// Below its precision budget a source with exact translations takes the value from its coefficients wherever
// that has a correct bit, and a caller that needs more of the value's bits asks again with more precision. The
// translation (an exact Taylor shift for every new centre) is for a value with no correct bit and for the
// budget, where no more precision can be asked for.
static void exact_translations_wait_for_the_budget() {
    const auto exact = from_roots({Complex(1), Complex(Rational(1) + power2(-90)), Complex(-2, 3)});
    d::GaussianPolynomial gaussian(exact.size());
    for (size_t j = 0; j < exact.size(); ++j) {
        fmpq_set_str(gaussian[j].real.value(), exact[j].real.str().c_str(), 10);
        fmpq_set_str(gaussian[j].imag.value(), exact[j].imag.str().c_str(), 10);
    }
    d::PolynomialSourceOptions options;
    options.input = d::InputSemantics::Exact;
    constexpr slong budget = 1024, asked = 192;
    options.exact_precision = budget;
    unsigned translations = 0;
    options.exact_taylor = [gaussian, &translations](acb_poly_t out, const acb_t center, slong bits) {
        d::GaussianRational shift;
        arf_get_fmpq(shift.real.value(), arb_midref(acb_realref(center)));
        arf_get_fmpq(shift.imag.value(), arb_midref(acb_imagref(center)));
        d::GaussianPolynomial shifted;
        if (!d::center_gaussian(shifted, gaussian, shift, 8000000))
            return false;
        d::enclose_gaussian(out, shifted, bits);
        ++translations;
        return true;
    };
    d::PolynomialSource source([exact](acb_poly_t out, slong bits) { enclosure(out, exact, bits); }, 3, false,
                               options);
    ComplexPolynomial derivative(3);
    for (size_t j = 1; j < exact.size(); ++j)
        derivative[j - 1] = exact[j] * Complex(Rational(j));
    const auto at = [&](const Complex &point, slong bits) {
        d::ComplexBall z;
        enclose(z.value(), point, 2 * budget);
        const auto result = joint(source, z, bits);
        require(result.enclosed, "a source with exact translations gave no value");
        require(contains(result.value.value(), evaluate(exact, point)) &&
                    contains(result.derivative.value(), evaluate(derivative, point)),
                "a source with exact translations lost its value or derivative");
        return result;
    };
    // Away from the roots, and 2^-40 from two of them: the coefficients' value has bits.
    at(Complex(3, Rational(1) / 2), asked);
    at(Complex(Rational(1) + power2(-40)), asked);
    require(translations == 0, "an exact translation was taken for a value its coefficients give");
    // On a root the coefficients give a ball about zero, the translation the exact zero.
    const auto on_root = at(Complex(1), asked);
    require(translations == 1 && acb_is_zero(on_root.value.value()),
            "a value with no correct bit was not taken from the exact translation");
    // At the budget, every point.
    at(Complex(3, Rational(1) / 2), budget);
    require(translations == 2, "the precision budget's value was not taken from the exact translation");
}
static void regeneration(const d::EnginePolicy *policy = nullptr) {
    std::vector<Complex> roots, exact_nodes;
    for (int i = 0; i < 8; ++i) {
        roots.emplace_back(i + 1);
        exact_nodes.emplace_back(-i - 1);
    }
    const auto p = from_roots(roots);
    d::PolynomialSource source([&](acb_poly_t out, slong bits) { enclosure(out, p, bits); }, 8, true);
    std::vector<d::ComplexBall> nodes(8), a;
    for (size_t i = 0; i < nodes.size(); ++i)
        enclose(nodes[i].value(), exact_nodes[i], 128);
    std::vector<slong> bits(8, 80);
    cps_options o;
    o.bits = 80;
    certpolysolve::detail::CandidateStatistics stats{};
    uint64_t maximum = 0;
    require(d::regenerate_residues(a, source, nodes, {}, {}, bits, o, stats, maximum, policy),
            "full regeneration failed");
    const auto expected = weights(p, exact_nodes);
    for (size_t i = 0; i < a.size(); ++i)
        require(contains(a[i].value(), expected[i]), "wrong full residue");
    const auto old_nodes = nodes, old_a = a;
    exact_nodes[2] = Complex(Rational(11) / 8, Rational(1) / 16);
    enclose(nodes[2].value(), exact_nodes[2], 128);
    const auto before = stats.evaluations;
    require(d::regenerate_residues(a, source, nodes, old_nodes, old_a, bits, o, stats, maximum, policy),
            "partial regeneration failed");
    const auto updated = weights(p, exact_nodes);
    for (size_t i = 0; i < a.size(); ++i)
        require(contains(a[i].value(), updated[i]), "partial identity disagreement");
    require(stats.evaluations - before < 8, "partial regeneration refreshed every node");
    // Regenerated nodes ARE allowed to be exact roots; zero residues must retain them.
    for (size_t i = 0; i < nodes.size(); ++i)
        enclose(nodes[i].value(), roots[i], 128);
    require(d::regenerate_residues(a, source, nodes, {}, {}, bits, o, stats, maximum, policy),
            "zero residue rejected");
    for (const auto &r : a)
        require(acb_is_zero(r.value()), "exact interpolation root is not zero");
    std::vector<d::Root> cover;
    std::vector<uint8_t> inactive;
    require(d::secular_inclusion_cover(cover, inactive, nodes, a, true, 1, 128, o),
            "zero-weight cover failed");
    require(cover.size() == 8 && std::count(inactive.begin(), inactive.end(), uint8_t(1)) == 8,
            "zero weights lost roots");
    for (size_t i = 0; i < cover.size(); ++i)
        require(contains(cover[i].box.value(), roots[i]), "exact node lost from cover");
    nodes[1] = nodes[0];
    require(!d::regenerate_residues(a, source, nodes, {}, {}, bits, o, stats, maximum, policy),
            "duplicate nodes accepted");
}
// The rebuild's evaluations keep their bits at complex nodes. Ball Horner
// re-boxes its error rectangle at every term and loses up to half a bit each;
// at degree 600 that is up to 300 bits, and the residues' ladder (starting at
// the nodes' precision plus 32, doubling) climbed past them. With disk errors
// the residues at nodes within 2^-85 of their roots are adequate at the
// ladder's second rung.
static void regeneration_keeps_bits() {
    constexpr slong degree = 600;
    d::PolynomialSource source(
        [](acb_poly_t out, slong) {
            acb_poly_zero(out);
            for (slong k = 0; k <= degree; ++k)
                acb_poly_set_coeff_si(out, k, 1);
        },
        degree, true);
    // The roots of 1 + x + ... + x^600: the 601st roots of unity but 1,
    // rounded to 85 bits.
    std::vector<d::ComplexBall> nodes(degree);
    for (slong k = 1; k <= degree; ++k) {
        d::ComplexBall &z = nodes[static_cast<size_t>(k - 1)];
        d::RealBall angle;
        arb_set_si(angle.value(), 2 * k);
        arb_div_si(angle.value(), angle.value(), degree + 1, 256);
        arb_sin_cos_pi(acb_imagref(z.value()), acb_realref(z.value()), angle.value(), 256);
        acb_get_mid(z.value(), z.value());
        acb_set_round(z.value(), z.value(), 85);
        acb_get_mid(z.value(), z.value());
    }
    std::vector<slong> bits(degree, 85);
    cps_options o;
    o.bits = 53;
    auto policy = d::maintained_method_policy(CPS_METHOD_SECULAR);
    for (const bool value_only : {true, false}) {
        policy.regeneration.value_only = value_only;
        certpolysolve::detail::CandidateStatistics stats{};
        uint64_t maximum = 0;
        std::vector<d::ComplexBall> residues;
        require(d::regenerate_residues(residues, source, nodes, {}, {}, bits, o, stats, maximum, &policy),
                "the rebuild at roots of unity failed");
        require(residues.size() == static_cast<size_t>(degree), "the rebuild lost residues");
        for (const auto &w : residues)
            require(acb_is_finite(w.value()) && acb_rel_accuracy_bits(w.value()) >= 85, "a residue is short");
        require(maximum <= 2 * (85 + 32), "the rebuild climbed past the ladder's second rung");
    }
}
static void value_only_sources() {
    for (int degree : {3, 32}) {
        ComplexPolynomial exact(static_cast<size_t>(degree + 1));
        exact[0] = Complex(Rational(2) / 3, Rational(-3) / 5);
        exact[1] = Complex(4, -1);
        exact.back() = Complex(2, 1);
        d::PolynomialSource source([exact](acb_poly_t p, slong bits) { enclosure(p, exact, bits); },
                                   degree, false);
        d::ComplexBall query;
        const Complex center(Rational(3) / 4, Rational(1) / 8);
        enclose(query.value(), center, 192);
        arb_add_error_2exp_si(acb_realref(query.value()), -20);
        arb_add_error_2exp_si(acb_imagref(query.value()), -20);
        const auto value = source.value(query, 128);
        require(value.status == d::EvaluationStatus::Enclosed, "value-only coefficient evaluation failed");
        for (int x : {-1, 0, 1}) for (int y : {-1, 0, 1}) {
            const auto point = center + Complex(power2(-20) * x, power2(-20) * y);
            require(contains(value.value.value(), quotient(evaluate(exact, point), exact.back())),
                    "dense/sparse monic value-only enclosure lost exact sample");
        }
        acb_zero(query.value());
        require(contains(value.value.value(), quotient(evaluate(exact, center), exact.back())),
                "owned value result aliases its query");
    }
}
static void independent_regeneration_precision() {
    std::vector<Complex> roots, exact_nodes;
    for (int j = 1; j <= 8; ++j) {
        roots.emplace_back(Rational(1) + j * power2(-80));
        exact_nodes.emplace_back(Rational(1) - j * power2(-80));
    }
    const auto exact = from_roots(roots);
    bool cancel = false, cancel_in_source = false;
    d::PolynomialSource source([&](acb_poly_t p, slong bits) {
        enclosure(p, exact, bits);
        if (cancel_in_source) cancel = true;
    }, 8, true);
    d::EnginePolicy policy;
    policy.model_stopping = d::SecularModelPolicy::Continuation;
    policy.regeneration = {.value_only = true, .reuse_products = true};
    policy.study_statistics = std::make_shared<d::StudyStatistics>();
    std::vector<d::ComplexBall> nodes(8), residues;
    std::vector<slong> bits(8, 80);
    for (size_t i = 0; i < nodes.size(); ++i) enclose(nodes[i].value(), exact_nodes[i], 192);
    cps_options options;
    options.cancel = [](void *flag) -> int { return *static_cast<bool *>(flag); };
    options.cancel_data = &cancel;
    certpolysolve::detail::CandidateStatistics stats{};
    uint64_t peak = 0;
    require(d::regenerate_residues(residues, source, nodes, {}, {}, bits, options, stats, peak, &policy),
            "split-precision cancellation-heavy regeneration failed");
    auto expected = weights(exact, exact_nodes);
    for (size_t i = 0; i < nodes.size(); ++i)
        require(contains(residues[i].value(), expected[i]), "split-precision regeneration lost exact residue");
    const auto &work = *policy.study_statistics;
    require(work.product_builds == 8 && work.product_reuses > 0 &&
                work.source_peak_bits > work.product_peak_bits && work.source_joints == 0,
            "source cancellation rebuilt products or computed derivatives");
    // A new model must rebuild geometry-dependent products, even at the same precision.
    exact_nodes[0] = Complex(Rational(1) - 9 * power2(-80));
    enclose(nodes[0].value(), exact_nodes[0], 192);
    require(d::regenerate_residues(residues, source, nodes, {}, {}, bits, options, stats, peak, &policy),
            "changed-node regeneration failed");
    expected = weights(exact, exact_nodes);
    for (size_t i = 0; i < nodes.size(); ++i)
        require(contains(residues[i].value(), expected[i]), "stale product survived a geometry change");
    require(work.product_builds == 16, "model-local product cache escaped its geometry");
    const auto retained = residues;
    cancel_in_source = true;
    bits.assign(8, 81); // Force an uncached coefficient enclosure.
    require(!d::regenerate_residues(residues, source, nodes, {}, {}, bits, options, stats, peak, &policy),
            "cancelled regeneration was published");
    for (size_t i = 0; i < residues.size(); ++i)
        require(acb_equal(residues[i].value(), retained[i].value()), "partial model replaced the retained model");
}
static void preliminary_contracts() {
    std::vector<Complex> roots;
    for (int j = 1; j <= 8; ++j)
        roots.emplace_back(Rational(j) / 3, Rational(j % 3 - 1) / 5);
    auto polynomial = input(from_roots(roots));
    for (const auto goal : {CPS_ROOT_APPROXIMATION, CPS_ROOT_ACCURACY}) {
        auto policy = d::production_policy();
        policy.study = d::StudyMethod::Secular;
        policy.refinement = d::RefinementMethod::Secular;
        policy.model_stopping = d::SecularModelPolicy::Continuation;
        policy.secular_preliminary = true;
        policy.regeneration = {.value_only = true, .reuse_products = true};
        policy.study_statistics = std::make_shared<d::StudyStatistics>();
        size_t snapshots = 0;
        policy.snapshot = [&](d::NumericalSnapshot event) {
            ++snapshots;
            // The observer owns its copy; corrupting it cannot change the solve.
            for (auto &point : event.points) acb_indeterminate(point.value());
        };
        d::SolveState state(certpolysolve::detail::dense(*polynomial.get()), policy);
        cps_options options;
        options.goal = goal;
        auto result = state.advance(&options);
        known(wrap(result), roots, 53);
        require(policy.study_statistics->preliminary_packets == 1,
                "secular driver skipped or repeated its preliminary phase");
        // 100 bits: a fresh solve, which runs its own preliminary phase once (one advance a solve state since
        // 18 October 2026).
        options.bits = 100;
        policy.study_statistics = std::make_shared<d::StudyStatistics>();
        d::SolveState tighter(certpolysolve::detail::dense(*polynomial.get()), policy);
        result = tighter.advance(&options);
        known(wrap(result), roots, 100);
        require(policy.study_statistics->preliminary_packets == 1,
                "a 100-bit solve skipped or repeated its preliminary phase");
        require(snapshots > 0, "candidate snapshot observer was never called");
    }
    // A huge exact translation can merge individually accurate local boxes
    // after outward rounding. Such a cover must not freeze all candidates.
    const auto exact = from_roots({Complex(-1), Complex(1)});
    auto source = std::make_unique<d::PolynomialSource>(
        [exact](acb_poly_t p, slong bits) { enclosure(p, exact, bits); }, 2, true);
    d::RootVector seeds(2);
    acb_set_si(seeds[0], -1);
    acb_set_si(seeds[1], 1);
    d::GaussianRational offset;
    fmpq_one(offset.real.value());
    fmpq_mul_2exp(offset.real.value(), offset.real.value(), 500);
    auto policy = d::production_policy();
    policy.secular_preliminary = true;
    policy.model_stopping = d::SecularModelPolicy::Continuation;
    d::SecularRefinement refinement(std::move(source), seeds.data(), 1, 85, policy, d::solve_local_model,
                                    d::RootCoordinates(std::move(offset)));
    cps_options options;
    options.goal = CPS_ROOT_APPROXIMATION;
    options.max_work_bits = 128;
    certpolysolve::detail::SearchStatistics search{};
    certpolysolve::detail::CandidateStatistics statistics{};
    refinement.step(options, search, statistics);
    std::vector<d::Root> exported;
    require(refinement.active() > 0 || refinement.cover(exported),
            "preliminary phase froze an unpublishable translated cover");
}
static void centered_source_identity() {
    d::IntegerPolynomial exact;
    fmpz_poly_one(exact.value());
    for (unsigned level = 0; level < 8; ++level) {
        fmpz_poly_sqr(exact.value(), exact.value());
        fmpz_poly_shift_left(exact.value(), exact.value(), 1);
        fmpz_poly_set_coeff_si(exact.value(), 0, 1);
    }
    certpolysolve::Polynomial polynomial;
    const slong degree = fmpz_poly_degree(exact.value());
    for (slong j = 0; j <= degree; ++j) {
        char *coefficient = fmpz_get_str(nullptr, 10, exact.value()->coeffs + j);
        polynomial.set(static_cast<uint64_t>(j), certpolysolve::Number(coefficient));
        flint_free(coefficient);
    }
    auto policy = d::production_policy();
    policy.study = d::StudyMethod::Secular;
    policy.refinement = d::RefinementMethod::Secular;
    policy.model_stopping = d::SecularModelPolicy::Continuation;
    policy.secular_preliminary = policy.model_driven_refresh = true;
    policy.regeneration = {.value_only = true, .reuse_products = true};
    d::SolveState state(certpolysolve::detail::dense(*polynomial.get()), policy);
    cps_options options;
    options.goal = CPS_ROOT_APPROXIMATION;
    const auto result = state.advance(&options);
    require(result.goal_met && result.statistics_v2.taylor_shifts > 0,
            "centered source identity regression did not exercise a completed translation");
    d::ComplexPolynomial original;
    acb_poly_set_fmpz_poly(original.value(), exact.value(), 1024);
    d::ComplexBall value, sum, weighted;
    for (const auto &root : result.roots) {
        acb_poly_evaluate(value.value(), original.value(), root.box.value(), 1024);
        require(acb_contains_zero(value.value()), "exported cover excludes a root of the ORIGINAL polynomial");
        acb_mul_ui(weighted.value(), root.box.value(), root.count, 1024);
        acb_add(sum.value(), sum.value(), weighted.value(), 1024);
    }
    acb_add(sum.value(), sum.value(), original.value()->coeffs + degree - 1, 1024);
    require(acb_contains_zero(sum.value()), "exported root multiset violates the original monic Vieta sum");
}
static void source_contracts() {
    d::PolynomialSource::Coefficients pinned;
    {
        const auto exact = from_roots({Complex(1), Complex(-1)});
        d::PolynomialSource source([exact](acb_poly_t p, slong bits) { enclosure(p, exact, bits); }, 2, true,
                                   {.input = d::InputSemantics::Exact});
        require(source.metadata().monic && source.metadata().input == d::InputSemantics::Exact,
                "source metadata lost exact semantics");
        pinned = source.coefficient_enclosure(80);
        for (slong precision = 81; precision < 100; ++precision)
            source.coefficient_enclosure(precision);
    }
    d::ComplexBall z, value, derivative;
    acb_set_si(z.value(), 2);
    require(pinned->evaluator->evaluate(value.value(), derivative.value(), z.value(), 128) &&
                contains(value.value(), Complex(3)),
            "prepared storage expired with source/cache");
    std::unique_ptr<d::SecularSource> secular;
    {
        certpolysolve::SecularEquation input{{.residue = certpolysolve::Number("1"), .pole = certpolysolve::Number("0")}};
        secular = std::make_unique<d::SecularSource>(*input.get());
    }
    acb_zero(z.value());
    // p'/p at the original pole 0 (p = x - 1 up to its scale): -1.
    const auto at_pole = joint(*secular, z, 128);
    d::ComplexBall quotient;
    acb_div(quotient.value(), at_pole.derivative.value(), at_pole.value.value(), 128);
    require(at_pole.enclosed && !acb_contains_zero(at_pole.value.value()) && contains(quotient.value(), Complex(-1)),
            "owned secular polynomial ratio failed at original pole");
    bool rejected = false;
    try {
        d::RootVector invalid(-1);
    } catch (const std::length_error &) {
        rejected = true;
    }
    require(rejected, "negative vector allocation accepted");
    d::RootVector vector(2);
    rejected = false;
    try {
        vector.view()[2];
    } catch (const std::out_of_range &) {
        rejected = true;
    }
    require(rejected, "out of bounds numerical view accepted");
}
static void interactions() {
    for (bool clustered : {false, true}) {
        std::vector<std::complex<double>> points;
        for (unsigned i = 0; i < 384; ++i) {
            const double angle = 2 * std::numbers::pi * static_cast<double>(i) / 384;
            const double scale = clustered && i < 192 ? 1.e-12 : 1;
            points.emplace_back(scale * std::cos(angle), scale * std::sin(angle));
        }
        d::InteractionTree tree(points);
        size_t near_total = 0;
        for (size_t i = 0; i < points.size(); ++i) {
            std::complex<double> far, direct = 0;
            std::vector<size_t> near;
            double tail;
            require(tree.evaluate(i, far, near, tail), "tree nonfinite");
            for (size_t j : near) {
                require(j != i, "tree included self");
                far += 1.0 / (points[i] - points[j]);
            }
            double absolute_sum = 0;
            for (size_t j = 0; j < points.size(); ++j)
                if (j != i) {
                    const auto term = 1.0 / (points[i] - points[j]);
                    direct += term;
                    absolute_sum += std::abs(term);
                }
            require(std::abs(far - direct) <= tail + 2.e-14 * absolute_sum, "tree/direct disagreement");
            near_total += near.size();
        }
        require(near_total < points.size() * points.size() / 2, "tree never used a far field");
    }
    // Coincident binary64 coordinates remain direct neighbours, so a caller
    // can use the original arbitrary-precision points to resolve the distance.
    const std::vector<std::complex<double>> duplicates(20, {1, 0});
    d::InteractionTree tree(duplicates);
    std::complex<double> far;
    std::vector<size_t> near;
    double tail;
    require(tree.evaluate(3, far, near, tail) && near.size() == 19 && far == std::complex<double>(0),
            "coincident tree coordinates were approximated");
    // Copies one unit apart that stand for points farther apart than that: a
    // target one unit of its imaginary part (2^-72) from twenty copies at one
    // position, every copy within 2^-53 of its point. Without the rounding
    // errors the twenty are one point mass at 2^-72, taken as far fields (the
    // packets' tree on the spiral toward 0, centred: no correct digit).
    std::vector<std::complex<double>> rounded(21, {1, 0x1p-20});
    rounded[20] = {1, 0x1p-20 + 0x1p-72};
    const std::vector<double> errors(rounded.size(), 0x1p-53);
    d::InteractionTree exact_copies(rounded), inexact_copies(rounded, errors);
    require(exact_copies.evaluate(20, far, near, tail) && near.size() < 20,
            "the fixture no longer makes the copies a far field");
    require(inexact_copies.evaluate(20, far, near, tail) && near.size() == 20 && far == std::complex<double>(0),
            "copies the rounding cannot place were taken as a far field");
    // Points of two doubles: twenty whose first parts coincide at 1 and whose
    // second parts spread them over 2^-66, and a target 2^-60 away. With their
    // second parts the twenty are a far field, as accurate as the series; with
    // first parts alone they were neighbours (one position, distance zero),
    // which nested_6_256's clusters paid for in exact interactions.
    std::vector<std::complex<double>> highs(21, {1, 0}), lows(21);
    for (size_t k = 0; k < 20; ++k)
        lows[k] = {std::ldexp(static_cast<double>(k), -70), std::ldexp(static_cast<double>(k % 5), -70)};
    lows[20] = {std::ldexp(1.0, -60), 0};
    const std::vector<double> exact(highs.size(), 0);
    d::InteractionTree two_part(highs, exact, lows);
    require(two_part.evaluate(20, far, near, tail) && near.size() < 20,
            "points told apart by their second parts were not taken as a far field");
    std::complex<double> direct = 0;
    double magnitudes = 0;
    for (size_t k = 0; k < 20; ++k) {
        const auto term = 1.0 / (lows[20] - lows[k]);
        direct += term;
        magnitudes += std::abs(term);
    }
    for (size_t j : near)
        far += 1.0 / (lows[20] - lows[j]);
    require(std::abs(far - direct) <= tail + 2.e-14 * magnitudes, "the two-part far field disagrees with the direct sum");
}
static void native() {
    for (int degree : {8, 16})
        for (bool complex : {false, true}) {
            std::vector<Complex> roots, poles;
            for (int i = 0; i < degree; ++i) {
                roots.emplace_back(i + 1, complex ? Rational((i % 3) - 1) / 8 : Rational(0));
                poles.emplace_back(Rational(2 * i + 1) / 2);
            }
            const auto p = from_roots(roots), a = weights(p, poles);
            std::vector<certpolysolve::SecularEquation::Term> terms;
            for (size_t i = 0; i < a.size(); ++i)
                terms.push_back({.residue = certpolysolve::Number::complex(certpolysolve::Number(a[i].real.str()),
                                                                       certpolysolve::Number(a[i].imag.str())),
                                 .pole = certpolysolve::Number::complex(certpolysolve::Number(poles[i].real.str()),
                                                                    certpolysolve::Number(poles[i].imag.str()))});
            certpolysolve::SecularEquation s(terms);
            d::SecularSource source(*s.get());
            ComplexPolynomial dp(p.size() - 1);
            for (size_t i = 1; i < p.size(); ++i)
                dp[i - 1] = p[i] * Complex(Rational(i));
            for (const auto &z : poles) { // At original poles the monic polynomial limit is finite.
                d::ComplexBall point, value, dvalue;
                enclose(point.value(), z, 256);
                require(source.evaluate(value.value(), dvalue.value(), point.value(), 256), "pole limit nonfinite");
                require(contains(value.value(), evaluate(p, z)) && contains(dvalue.value(), evaluate(dp, z)),
                        "pole-limit identity failed");
            }
            // The native route must certify converged proposals before restarting;
            // a modest precision ceiling catches the former restart/growth cycle.
            // (Until 18 October 2026 a session advanced a unit at a time, at most 151 calls; one solve since.)
            // (An exact allowance alone has no field in the header's Options: C options, solved through the test
            // support.)
            auto o = certpolysolve_test::raw_options();
            require(certpolysolve_test::set_exact_allowance(o.get(), 0) == CPS_OK, "setting the exact allowance failed");
            require(cps_options_set_precision(o.get(), 80) == CPS_OK && cps_options_set_max_work_bits(o.get(), 512) == CPS_OK,
                    "setting the precision failed");
            require(cps_options_set_work_budget(o.get(), 151) == CPS_OK && cps_options_set_live_cells(o.get(), 65536) == CPS_OK,
                    "setting the search budget failed");
            const auto r = certpolysolve_test::solve(s, o.get());
            known(r, roots, 80);
            require(r.statistics().preparation_ns == 0, "native solve expanded coefficients");
            for (const auto &component : r.components())
                require(component.factor_index == 0, "native component has a nonexistent factor index");
            require(cps_options_set_precision(o.get(), 160) == CPS_OK && cps_options_set_max_work_bits(o.get(), 512) == CPS_OK,
                    "setting the precision failed");
            require(cps_options_set_work_budget(o.get(), 1000) == CPS_OK && cps_options_set_live_cells(o.get(), 65536) == CPS_OK,
                    "setting the search budget failed");
            known(certpolysolve_test::solve(s, o.get()), roots, 160);
        }
}
#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
// The local Taylor cache, the research build's Series test engine's evaluator.
static void local_evaluator() {
    // Dense high-degree tail with exact dyadic coefficients. At a small
    // argument both value and derivative tails are tiny but nonzero.
    ComplexPolynomial p(97);
    for (size_t j = 0; j < p.size(); ++j)
        p[j] = Complex(Rational(static_cast<int>(j % 7) + 1) / Integer(1 << static_cast<int>(j % 9)));
    ComplexPolynomial dp(p.size() - 1);
    for (size_t j = 1; j < p.size(); ++j)
        dp[j - 1] = p[j] * Complex(Rational(j));
    d::ComplexPolynomial ball;
    enclosure(ball.value(), p, 256);
    d::LocalPolynomialEvaluator evaluator(ball.value());
    for (const Complex z :
         {Complex(Rational(1) / 1024), Complex(Rational(1) / 2048, Rational(1) / 4096), Complex(0)}) {
        d::ComplexBall argument, v, dv;
        enclose(argument.value(), z, 256);
        require(evaluator.evaluate(v.value(), dv.value(), argument.value(), 128), "local evaluator failed");
        require(contains(v.value(), evaluate(p, z)) && contains(dv.value(), evaluate(dp, z)),
                "local tail enclosure failed");
        // The same cache must enclose an entire uncertain input region.
        arb_add_error_2exp_si(acb_realref(argument.value()), -50);
        require(evaluator.evaluate(v.value(), dv.value(), argument.value(), 128), "local ball evaluator failed");
        const Complex endpoint = z + Complex(Rational(1) / (Integer(1) << 50));
        require(contains(v.value(), evaluate(p, endpoint)) && contains(dv.value(), evaluate(dp, endpoint)),
                "local argument uncertainty was discarded");
    }
    require(evaluator.discarded() > 0 && evaluator.cache_hits() > 0,
            "local test did not exercise band/cache");
}
#endif
static void structured_solves() {
    // Persistent, Tree, Series and Direct are test engines of the research build.
    for (auto engine : {
#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
             d::CandidateEngine::Persistent, d::CandidateEngine::Tree, d::CandidateEngine::Series,
             d::CandidateEngine::Direct,
#endif
             d::CandidateEngine::Adaptive}) {
        for (int separation : {40, 300}) {
            std::vector<Complex> roots{{-3}, {2}, {0, 1}, {0, -1}};
            for (int j = 0; j < 4; ++j)
                roots.emplace_back(Rational(1) + Rational(j) / (Integer(1) << separation),
                                   Rational(j % 2) / (Integer(1) << (separation + 2)));
            auto p = input(from_roots(roots));
            d::SolveState state(certpolysolve::detail::dense(*p.get()), {.backend = true, .candidates = engine});
            cps_options o;
            o.bits = 80;
            o.max_work_bits = 8192;
            o.work_units = 1000;
            const auto result = wrap(state.advance(&o));
            known(result, roots, 80);
        }
    }
#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
    // The research build's Seeded test engine (the structured route; its seeds were the Mandelbrot recurrence's,
    // deleted in W8's S2) on the expanded Mandelbrot polynomials, with independent Sturm--Cauchy root counts.
    ComplexPolynomial p{{1}};
    for (unsigned level = 1; level <= 3; ++level) {
        p = multiply(p, p);
        p.insert(p.begin(), Complex(1));
        auto polynomial = input(p);
        d::SolveState state(certpolysolve::detail::dense(*polynomial.get()),
                            {.backend = true, .candidates = d::CandidateEngine::Seeded});
        cps_options o;
        o.bits = 53;
        const auto result = wrap(state.advance(&o));
        require(result.goal_met(), "seeded recurrence solve incomplete");
        for (const auto &root : result.roots())
            require(rectangle_count(p, root) == root.count, "recurrence independent root-count disagreement");
    }
#endif
}
int main() {
    certpolysolve::ThreadScope lifetime;
    try {
        conjugate_cover_realness();
        regeneration();
        regeneration_keeps_bits();
        for (unsigned mask = 1; mask < 4; ++mask) {
            d::EnginePolicy policy;
            policy.regeneration = {.value_only = bool(mask & 1), .reuse_products = bool(mask & 2)};
            regeneration(&policy);
        }
        value_only_sources();
        independent_regeneration_precision();
        preliminary_contracts();
        centered_source_identity();
        centered_zero_initialization();
        hopeless_preliminary_validation();
        power_composition();
        exact_local_evaluation();
        exact_translations_wait_for_the_budget();
        source_contracts();
        interactions();
#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
        local_evaluator();
#endif
        structured_solves();
        native();
        std::cout << "structured identities and roots passed\n";
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
