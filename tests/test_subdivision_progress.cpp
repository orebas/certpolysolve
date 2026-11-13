#include "reference_math.hpp"
#include "exact_winding_api.hpp"
#include "internal.hpp"
#include "methods/subdivision/subdivision.hpp"
#include <certpolysolve/certpolysolve.hpp>
#include <algorithm>
#include <array>
#include <iostream>
#include <limits>

namespace d = certpolysolve::detail;
using reference::Rational;
using reference::require;

namespace {
void declared_uncertainty_preserves_counted_contraction() {
    // Supplied state, not a native global-start trajectory. Every family member
    // F_eta=(z-1)^2+eta is square-free, eta in [r^2/32,15r^2/32].
    constexpr slong radius_exponent = -54;
    const Rational r = reference::power2(radius_exponent);
    const Rational eta_lower = r*r/32, eta_upper = 15*r*r/32;
    const Rational midpoint = 1+r*r/4, uncertainty = 7*r*r/32;
    require(eta_lower > 0 && eta_upper < r*r, "invalid uniform Pellet margin");
    // Exact family Newton values at x=1+6r: 1-eta/(6r).
    const Rational y_lower = 1-eta_upper/(6*r), y_upper = 1-eta_lower/(6*r);
    require((y_upper-y_lower)/2 == 7*r/192 && (y_upper+y_lower)/2 == 1-r/24,
            "incorrect independent Newton image of uncertain family");
    require(7*r/192 > r/32 && r/24 < r/16, "fixture does not obstruct accuracy before a valid snap");
    for (const slong precision : {256, 512}) {
        const d::Encloser enclose = [&](acb_poly_t out, slong bits) {
            acb_poly_zero(out);
            d::ComplexBall constant;
            d::RealBall error;
            d::Number(midpoint.str()).enclose(acb_realref(constant.value()), bits);
            d::Number(uncertainty.str()).enclose(error.value(), bits);
            require(arb_is_exact(error.value()), "declared coefficient radius rounded");
            arb_add_error_arf(acb_realref(constant.value()), arb_midref(error.value()));
            acb_poly_set_coeff_acb(out, 0, constant.value());
            acb_poly_set_coeff_si(out, 1, -2);
            acb_poly_set_coeff_si(out, 2, 1);
        };
        d::ComplexPolynomial polynomial;
        enclose(polynomial.value(), precision);
        d::ComplexBall point, value, derivative, proposal, center;
        d::Number(Rational(1+6*r).str()).enclose(acb_realref(point.value()), precision);
        acb_poly_evaluate2(value.value(), derivative.value(), polynomial.value(), point.value(), precision);
        require(!acb_contains_zero(derivative.value()), "uncertain-family annular derivative is unresolved");
        acb_div(proposal.value(), value.value(), derivative.value(), precision);
        acb_mul_2exp_si(proposal.value(), proposal.value(), 1);
        acb_sub(proposal.value(), point.value(), proposal.value(), precision);
        // An independent exact endpoint image forces this width at every P.
        d::Float propagated;
        arf_set_mag(propagated.value(), arb_radref(acb_realref(proposal.value())));
        arf_mul_ui(propagated.value(), propagated.value(), 192, ARF_PREC_EXACT, ARF_RND_NEAR);
        d::RealBall minimum;
        d::Number(Rational(7*r).str()).enclose(minimum.value(), precision);
        require(arf_cmp(propagated.value(), arb_midref(minimum.value())) >= 0,
                "Newton enclosure misses the exact family image");
        acb_one(center.value());
        d::RealBall seed_radius;
        d::Number(Rational(2*r).str()).enclose(seed_radius.value(), precision);
        d::Root seed;
        require(d::certify_cluster_box(seed, polynomial.value(), center.value(), arb_midref(seed_radius.value()),
                                      2, 1, precision), "uncertain-family seed lacks uniform count two");
        // Supply the original exact square proved by those disk counts.
        // Publication may inflate its radius; that is not fixture geometry.
        acb_set(seed.box.value(), center.value());
        // Adding an arf error can round its mag upward too. Start from the
        // exact magnitude one and scale by 2^(radius_exponent+1) = 2r.
        for (auto coordinate : {acb_realref(seed.box.value()), acb_imagref(seed.box.value())}) {
            mag_one(arb_radref(coordinate));
            mag_mul_2exp_si(arb_radref(coordinate), arb_radref(coordinate), radius_exponent+1);
        }
        d::Float supplied_half;
        arf_set_mag(supplied_half.value(), arb_radref(acb_realref(seed.box.value())));
        require(arf_equal(supplied_half.value(), arb_midref(seed_radius.value())),
                "supplied-state square radius is not exactly 2r");
        cps_options options;
        options.goal = CPS_ROOT_APPROXIMATION;
        options.bits = 53;
        options.max_work_bits = precision;
        require(!d::meets_goal(seed, options, precision), "uncertain-family seed already meets the goal");
        d::SearchStatistics statistics{};
        d::Subdivision search({seed}, 1, precision, statistics);
        // Default generic policy must preserve an immediately certifiable wide
        // midpoint. Constructor doubles the seed half-width, giving R=4r,e=2.
        require(search.step(options, enclose, 0) == CPS_STOP_NONE,
                "uncertain-family contraction stopped before counting");
        std::vector<d::Root> cover;
        require(search.complete_cover(cover) && cover.size() == 1 && cover[0].count == 2 &&
                    d::meets_goal(cover[0], options, precision) && statistics.newton_successes == 1 &&
                    statistics.precision_increases == 0 && statistics.boxes_split == 0,
                "declared uncertainty prevented a valid legacy counted contraction");
        d::Float half;
        arf_set_mag(half.value(), arb_radref(acb_realref(cover[0].box.value())));
        d::Number(r.str()).enclose(minimum.value(), precision);
        // The returned published box contains the nominal r-square and stays
        // below the next dyadic scale; its outward radius need not equal r.
        d::Float next_scale;
        arf_mul_2exp_si(next_scale.value(), arb_midref(minimum.value()), 1);
        require(arf_cmp(half.value(), arb_midref(minimum.value())) >= 0 &&
                    arf_cmp(half.value(), next_scale.value()) < 0,
                "unexpected supplied-state trial radius");
        std::cout << "uncertain-family precision=" << precision << " uniform_count=2 contractions="
                  << statistics.newton_successes << " goal=1\n";
    }
}

void exact_points_release_support() {
    for (int center_value : {0, 1}) {
        d::ComplexBall center;
        acb_set_si(center.value(), center_value);
        const d::Encloser enclose = [&](acb_poly_t out, slong) {
            acb_poly_zero(out);
            acb_poly_set_coeff_si(out, 0, -center_value);
            acb_poly_set_coeff_si(out, 1, 1);
        };
        d::ComplexPolynomial polynomial;
        enclose(polynomial.value(), 128);
        d::Float radius;
        arf_one(radius.value());
        d::Root point;
        require(d::certify_cluster_box(point, polynomial.value(), center.value(), radius.value(), 1, 2, 128),
                "exact-point fixture was not certified");
        d::SearchStatistics statistics{};
        d::Subdivision search({point}, 2, 128, statistics);
        cps_options options;
        options.goal = CPS_ROOT_APPROXIMATION;
        options.bits = 53;
        require(search.step(options, enclose, 0) == CPS_STOP_NONE && !search.pending(),
                "exact point did not finish without further cells");
        require(search.cells() == 0, "finished exact points still consume the search cell budget");
        options.bits = 4096;
        options.max_work_bits = 8192;
        options.cluster_radius = d::Number("1@-4096");
        for (auto goal : {CPS_ROOT_ACCURACY, CPS_ROOT_APPROXIMATION, CPS_CLUSTER_RADIUS}) {
            options.goal = goal;
            search.reopen(options);
            std::vector<d::Root> cover;
            require(!search.pending() && search.complete_cover(cover) && cover.size() == 1,
                    "a finished exact point reopened for a stronger query");
            require(acb_equal(cover.front().box.value(), center.value()) &&
                        cover.front().count == 2 && cover.front().multiplicity == 2,
                    "exact-point completion changed its location or multiplicity");
        }
    }
}

void small_cell_budget() {
    // Axis fuzz seed 13398634035015802212, its doubled-precision followup.
    // Three nearby exact points must not retain their old four-cell supports.
    const std::vector<reference::Complex> known{
        {1, reference::power2(-57)}, {1, -reference::power2(-57)},
        {Rational("-1603699/2097152")}, {Rational("-122885/262144")},
        {Rational("511/512")}, {Rational("454709/1048576")},
        {Rational("59/133")}, {1}};
    auto polynomial = reference::from_roots(known);
    reference::Polynomial real;
    for (const auto &coefficient : polynomial) {
        require(coefficient.imag == 0, "conjugate fixture has nonreal coefficients");
        real.push_back(coefficient.real);
    }
    reference::make_primitive(real);
    certpolysolve::Options options;
    options.method = certpolysolve::Method::piecewise;
    options.goal = certpolysolve::Goal::approximation();
    options.bits = 100;
    options.max_bits = 192;
    options.work_budget = 1000;
    options.cell_budget = 16;
    const auto result = certpolysolve::solve(reference::input(real), options);
    require(result.ok() && result.complete() && result.roots().size() == known.size(),
            "finished exact-point supports prevent the remaining roots from meeting their goal");
    std::vector<unsigned> visits(known.size());
    for (const auto &root : result.roots()) {
        const auto box = reference::rectangle(root);
        unsigned contained = 0;
        for (size_t j = 0; j < known.size(); ++j) {
            const auto &z = known[j];
            if (box.real_lower <= z.real && z.real <= box.real_upper &&
                box.imag_lower <= z.imag && z.imag <= box.imag_upper) {
                ++contained;
                ++visits[j];
            }
        }
        require(contained == 1 && root.count == 1 && root.multiplicity == 1,
                "small-budget output disagrees with the exact planted roots");
        reference::check_accuracy(root, 100);
    }
    require(std::all_of(visits.begin(), visits.end(), [](unsigned count) { return count == 1; }),
            "small-budget output loses or repeats a planted root");
    reference::check_disjoint(result.roots());
}

void finite_cap_split() {
    // Reduced from algebraic fuzz seed 4328692777287013451. The third
    // quadratic has the close nonreal pair 25/16 +/- i sqrt(2) 2^-26.
    // A precision request for some coefficient must not by itself abandon
    // a split at the cap: geometry still separates this pair at 64 bits.
    const std::vector<reference::ComplexPolynomial> factors{
        {{Rational("3400831/524288")}, {Rational("-163/32")}, {1}},
        {{Rational("-32198007823/68719476736"), Rational("-426591/33554432")},
         {Rational("-47399/131072"), Rational("9/128")}, {1}},
        {{Rational("5497558138880001/2251799813685248")}, {Rational("-25/8")}, {1}},
        {{Rational("8601/16384")}, {1}},
        {{Rational("-396691/65536"), Rational("754/59")}, {1}}};
    reference::ComplexPolynomial polynomial{{1}};
    for (const auto &factor : factors)
        polynomial = reference::multiply(polynomial, factor);
    certpolysolve::Options options;
    options.method = certpolysolve::Method::classical;
    options.goal = certpolysolve::Goal::accuracy();
    options.bits = 16;
    options.max_bits = 64;
    options.work_budget = 100000;
    options.cell_budget = 65536;
    const auto result = certpolysolve::solve(reference::input(polynomial), options);
    require(result.ok() && result.complete() && result.roots().size() == 8 && result.work_bits() <= 64,
            "the soft precision request abandoned a useful split at the caller's cap");
    for (const auto &root : result.roots()) {
        int total = 0;
        for (const auto &factor : factors) {
            const auto count = reference::rectangle_count(factor, root);
            require(count.has_value(), "the finite-cap output touches a root on its boundary");
            total += *count;
        }
        require(total == 1 && root.count == 1 && root.multiplicity == 1,
                "the finite-cap singleton count disagrees with exact Sturm-Cauchy counts");
        reference::check_accuracy(root, 16);
    }
    reference::check_disjoint(result.roots());
}

// Two exact stationary samples, neither a root. The small pair is +/-1/8;
// the other roots have modulus >2. The last a makes sum 1/(1-a) exactly zero.
const std::array<Rational, 5> squares{Rational("1/64"), Rational("17/4"), Rational("9/2"),
                                      Rational(5), Rational("3841/565")};

bool contains_positive_sqrt(const Rational &lo, const Rational &hi, const Rational &a) {
    return hi >= 0 && hi * hi >= a && (lo <= 0 || lo * lo <= a);
}

reference::Rectangle rectangle(const d::Root &root) {
    std::array<Rational, 4> ends;
    for (size_t j = 0; j < ends.size(); ++j)
        ends[j] = reference::dyadic(root.endpoints[j].mantissa_text(), root.endpoints[j].exponent_text());
    return {ends[0], ends[1], ends[2], ends[3], root.count};
}

// Check the actual root contract independently of the search's certificates.
void check_sample_pair(const std::vector<d::Root> &cover, const std::array<Rational, 2> &roots,
                       uint64_t bits) {
    std::array<unsigned, 2> visits{};
    std::vector<reference::Rectangle> boxes;
    for (const auto &root : cover) {
        const auto box = rectangle(root);
        unsigned count = 0;
        for (size_t j = 0; j < roots.size(); ++j)
            if (box.real_lower <= roots[j] && roots[j] <= box.real_upper &&
                box.imag_lower <= 0 && 0 <= box.imag_upper) { ++count; ++visits[j]; }
        require(count == 1 && root.count == 1 && root.multiplicity == 1,
                "sample pair has a wrong count or multiplicity");
        const Rational scale = box.real_lower > 0 ? box.real_lower
            : (box.real_upper < 0 ? Rational(-box.real_upper) : Rational(0));
        const Rational real_width = box.real_upper - box.real_lower;
        const Rational imag_width = box.imag_upper - box.imag_lower;
        require(real_width * real_width + imag_width * imag_width <=
                    4 * scale * scale * reference::power2(-2 * static_cast<int>(bits)),
                "sample pair misses individual relative accuracy");
        boxes.push_back(box);
    }
    require(visits[0] == 1 && visits[1] == 1, "sample pair loses or duplicates a root");
    reference::check_disjoint(boxes);
}

void ring_support_preserves_counts() {
    // The quarter root lies on closed cell boundaries inside an initially
    // connected ring. Scanning old cells must preserve all certified covers
    // throughout partial splits, including exact axis roots and empty supports.
    constexpr slong ring_degree = 16;
    constexpr uint64_t goal_bits = 53;
    constexpr slong initial_bits = goal_bits + 32;
    constexpr unsigned action_limit = 2048; // test guard, not a solver threshold
    reference::ComplexPolynomial exact(static_cast<size_t>(ring_degree + 2));
    exact[0] = {1};
    exact[1] = {-4};
    exact[static_cast<size_t>(ring_degree)] = {-1};
    exact[static_cast<size_t>(ring_degree + 1)] = {4};
    const d::Encloser enclose = [](acb_poly_t out, slong) {
        acb_poly_zero(out);
        acb_poly_set_coeff_si(out, 0, 1);
        acb_poly_set_coeff_si(out, 1, -4);
        acb_poly_set_coeff_si(out, ring_degree, -1);
        acb_poly_set_coeff_si(out, ring_degree + 1, 4);
    };
    d::ComplexPolynomial polynomial;
    enclose(polynomial.value(), initial_bits);
    d::ComplexBall center;
    d::Float radius;
    arf_set_si(radius.value(), 2);
    d::Root seed;
    require(d::certify_cluster_box(seed, polynomial.value(), center.value(), radius.value(),
                                  ring_degree + 1, 1, initial_bits),
            "ring fixture has no initial global certificate");
    d::SearchStatistics statistics{};
    d::Subdivision search({seed}, 1, initial_bits, statistics);
    cps_options options;
    options.goal = CPS_ROOT_ACCURACY;
    options.bits = goal_bits;
    options.max_work_bits = 2048;
    for (unsigned action = 0; action < action_limit && search.pending(); ++action) {
        require(search.step(options, enclose, 10000) == CPS_STOP_NONE, "ring search stopped early");
        std::vector<d::Root> intermediate;
        if (search.complete_cover(intermediate)) {
            uint64_t count = 0;
            for (const auto &root : intermediate) count += root.count;
            require(count == ring_degree + 1 && d::disjoint(intermediate),
                    "partial ring split lost coverage or overlapped certificates");
        }
    }
    std::vector<d::Root> cover;
    require(!search.pending() && search.complete_cover(cover) && cover.size() == ring_degree + 1,
            "ring support did not finish with individual roots");
    unsigned quarter_visits = 0;
    for (const auto &root : cover) {
        const auto box = rectangle(root);
        require(root.count == 1 && root.multiplicity == 1 &&
                    reference::rectangle_count(exact, box.real_lower, box.real_upper,
                                                box.imag_lower, box.imag_upper) == 1,
                "ring output disagrees with independent exact winding counts");
        if (box.real_lower <= Rational(1)/4 && Rational(1)/4 <= box.real_upper &&
            box.imag_lower <= 0 && 0 <= box.imag_upper) ++quarter_visits;
    }
    require(quarter_visits == 1 && d::disjoint(cover), "ring output loses or repeats the quarter root");
}

void resolved_outside_sample_pair(d::NewtonAccuracy accuracy = d::NewtonAccuracy::Legacy) {
    // Roots a and -a+a^3; p(0),p'(0) are resolved, but their modified Newton
    // point is 2a-2/a. The old center-only refusal takes 3,096 actions here.
    constexpr int height = 128;
    constexpr size_t action_limit = 1024; // regression budget, not an algorithm threshold
    constexpr uint64_t goal_bits = 53;
    constexpr slong initial_bits = goal_bits + 32;
    const Rational a = reference::power2(-height);
    const std::array<Rational, 2> roots{a, -a + a * a * a};
    const reference::Polynomial polynomial{-a * a + a * a * a * a, -a * a * a, 1};
    for (const auto &z : roots) require(reference::evaluate(polynomial, z) == 0, "invalid sample pair fixture");
    std::vector<d::Number> coefficients;
    for (const auto &x : polynomial) coefficients.emplace_back(x.str());
    const d::Encloser enclose = [&](acb_poly_t out, slong bits) {
        acb_poly_zero(out);
        d::ComplexBall coefficient;
        for (size_t j = 0; j < coefficients.size(); ++j) {
            coefficients[j].enclose(acb_realref(coefficient.value()), bits);
            acb_poly_set_coeff_acb(out, static_cast<slong>(j), coefficient.value());
        }
    };
    d::ComplexPolynomial p;
    enclose(p.value(), initial_bits);
    d::ComplexBall center;
    d::Float radius;
    arf_one(radius.value());
    d::Root seed;
    require(d::certify_cluster_box(seed, p.value(), center.value(), radius.value(), 2, 1, initial_bits),
            "sample pair has no initial global count");
    d::SearchStatistics statistics{};
    d::Subdivision search({seed}, 1, initial_bits, statistics);
    cps_options options;
    options.goal = CPS_ROOT_ACCURACY;
    options.bits = goal_bits;
    options.max_work_bits = 4096;
    for (size_t action = 0; action < action_limit && search.pending(); ++action)
        require(search.step(options, enclose, 10000, std::numeric_limits<double>::infinity(), accuracy) == CPS_STOP_NONE, "sample pair stopped before its goal");
    std::vector<d::Root> cover;
    require(!search.pending() && search.complete_cover(cover),
            "resolved outside sample keeps bisecting instead of contracting");
    check_sample_pair(cover, roots, goal_bits);
}

void annular_coordinate_precision(d::NewtonAccuracy accuracy = d::NewtonAccuracy::Legacy) {
    // A valid warm cover can be much narrower than an ulp at the new working
    // precision. The intended annulus point rounds to the center at 64 bits.
    constexpr int coordinate_bits = 256;
    constexpr slong seed_bits = 2048; // enough for the exact expanded quadratic
    const Rational c = reference::power2(coordinate_bits);
    const Rational delta = reference::power2(-coordinate_bits - 4);
    const std::array<Rational, 2> roots{c - delta, c + delta};
    const reference::Polynomial polynomial{c * c - delta * delta, -2 * c, 1};
    std::vector<d::Number> coefficients;
    for (const auto &x : polynomial) coefficients.emplace_back(x.str());
    const d::Encloser enclose = [&](acb_poly_t out, slong bits) {
        acb_poly_zero(out);
        d::ComplexBall coefficient;
        for (size_t j = 0; j < coefficients.size(); ++j) {
            coefficients[j].enclose(acb_realref(coefficient.value()), bits);
            acb_poly_set_coeff_acb(out, static_cast<slong>(j), coefficient.value());
        }
    };
    d::ComplexPolynomial p;
    enclose(p.value(), seed_bits);
    d::ComplexBall center;
    d::Number(c.str()).enclose(acb_realref(center.value()), seed_bits);
    d::Float radius;
    arf_set_ui_2exp_si(radius.value(), 1, -coordinate_bits);
    d::Root seed;
    require(d::certify_cluster_box(seed, p.value(), center.value(), radius.value(), 2, 1, seed_bits),
            "coordinate fixture has no genuine initial count");
    d::SearchStatistics statistics{};
    d::Subdivision search({seed}, 1, 64, statistics);
    cps_options options;
    options.goal = CPS_ROOT_ACCURACY;
    options.bits = 53;
    options.max_work_bits = 128;
    require(search.step(options, enclose, 10000, std::numeric_limits<double>::infinity(), accuracy) == CPS_STOP_NONE && search.max_bits() == 128,
            "unrepresentable annulus point did not request precision");
    std::vector<d::Root> cover;
    require(search.pending() && search.complete_cover(cover) && cover.size() == 1 && cover[0].count == 2,
            "annulus precision request lost its certified partial cover");
    // Reach the actual ceiling, so reopening must unblock a stopped component.
    constexpr size_t ceiling_action_limit = 16; // bounded regression guard
    cps_stop_reason reason = CPS_STOP_NONE;
    for (size_t action = 0; action < ceiling_action_limit && reason == CPS_STOP_NONE; ++action)
        reason = search.step(options, enclose, 10000, std::numeric_limits<double>::infinity(), accuracy);
    require(reason == CPS_STOP_PRECISION && search.pending() && search.max_bits() == 128,
            "coordinate fixture did not stop at its precision ceiling");
    require(search.complete_cover(cover) && cover.size() == 1 && cover[0].count == 2,
            "annulus precision ceiling lost its certified partial cover");
    options.max_work_bits = 4096;
    search.reopen(options);
    for (size_t action = 0; action < 1024 && search.pending(); ++action)
        require(search.step(options, enclose, 10000, std::numeric_limits<double>::infinity(), accuracy) == CPS_STOP_NONE, "coordinate fixture failed after reopening");
    require(!search.pending() && search.complete_cover(cover), "coordinate fixture did not finish after reopening");
    check_sample_pair(cover, roots, options.bits);
}

void preserve_useful_uncertain_sample() {
    // Axis seed 13602400673810899049, its 192-bit followup. An uncertain
    // first value can still give a useful count-preserving contraction;
    // replacing it before examining the proposal loses this finite-cap goal.
    const Rational a = 247 * reference::power2(-47);
    const std::vector<reference::Complex> known{
        {a, reference::power2(-147)}, {a, -reference::power2(-147)},
        {reference::power2(-31), reference::power2(-66)},
        {reference::power2(-31), -reference::power2(-66)},
        {Rational("1086317488242689") * reference::power2(-89)},
        {517996543 * reference::power2(-68)}, {-7 * reference::power2(-41)}};
    const auto product = reference::from_roots(known);
    reference::Polynomial polynomial;
    for (const auto &coefficient : product) {
        require(coefficient.imag == 0, "uncertain-sample fixture is not real");
        polynomial.push_back(coefficient.real);
    }
    reference::make_primitive(polynomial);
    certpolysolve::Options options;
    options.method = certpolysolve::Method::secular;
    options.goal = certpolysolve::Goal::approximation();
    options.bits = 53;
    options.max_bits = 192;
    options.work_budget = 1000;
    const auto result = certpolysolve::solve(reference::input(polynomial), options);
    require(result.ok() && result.complete(), "replacing a useful first Newton sample loses a finite-cap goal");
    std::vector<unsigned> visits(known.size());
    for (const auto &root : result.roots()) {
        const auto box = reference::rectangle(root);
        uint64_t count = 0;
        for (size_t j = 0; j < known.size(); ++j)
            if (box.real_lower <= known[j].real && known[j].real <= box.real_upper &&
                box.imag_lower <= known[j].imag && known[j].imag <= box.imag_upper) {
                ++count;
                ++visits[j];
            }
        require(count == root.count, "uncertain-sample output disagrees with exact planted counts");
        const Rational real_lower = box.real_lower > 0 ? box.real_lower
            : (box.real_upper < 0 ? Rational(-box.real_upper) : Rational(0));
        const Rational imag_lower = box.imag_lower > 0 ? box.imag_lower
            : (box.imag_upper < 0 ? Rational(-box.imag_upper) : Rational(0));
        const Rational real_width = box.real_upper - box.real_lower;
        const Rational imag_width = box.imag_upper - box.imag_lower;
        const Rational allowance_squared = (real_lower * real_lower + imag_lower * imag_lower) *
            reference::power2(-2 * static_cast<int>(options.bits));
        require(real_width * real_width + imag_width * imag_width <= 4 * allowance_squared,
                "uncertain-sample output misses its relative approximation goal");
    }
    require(std::all_of(visits.begin(), visits.end(), [](unsigned count) { return count == 1; }),
            "uncertain-sample output loses or repeats a planted root");
    reference::check_disjoint(result.roots());
}

void uncertain_value_sample(d::NewtonAccuracy accuracy = d::NewtonAccuracy::Legacy) {
    // The cluster's center is an exact simple root. Seven small gaps make
    // its derivative tiny. At 7 * 192 + 16 bits that derivative is resolved
    // while cancellation still leaves the value uncertain (checked below).
    // A counted contraction should not require allocating subdivision cells.
    constexpr slong spacing_bits = 192;
    constexpr slong precision = 7 * spacing_bits + 16;
    std::vector<reference::Complex> known;
    for (int j = 0; j < 8; ++j)
        known.emplace_back(1 + j * reference::power2(-spacing_bits), 1);
    for (int j = 0; j < 8; ++j)
        known.emplace_back(4 + j, 4);
    const auto input = reference::input(reference::from_roots(known));
    const d::Encloser enclose = [&](acb_poly_t out, slong bits) {
        acb_poly_zero(out);
        d::ComplexBall coefficient;
        for (const auto &[index, value] : input.get()->terms) {
            value.enclose(coefficient.value(), bits);
            acb_poly_set_coeff_acb(out, static_cast<slong>(index), coefficient.value());
        }
    };
    d::ComplexPolynomial polynomial;
    enclose(polynomial.value(), precision);
    d::ComplexBall center, value, derivative;
    acb_one(center.value());
    arb_one(acb_imagref(center.value()));
    acb_poly_evaluate2(value.value(), derivative.value(), polynomial.value(), center.value(), precision);
    require(acb_contains_zero(value.value()) && !acb_is_zero(value.value()) &&
                !acb_contains_zero(derivative.value()),
            "fixture does not have an uncertain value and a resolved derivative");
    d::Float radius;
    arf_set_ui_2exp_si(radius.value(), 1, -spacing_bits / 4);
    d::Root cluster;
    require(d::certify_cluster_box(cluster, polynomial.value(), center.value(), radius.value(), 8, 1, precision),
            "uncertain-value fixture lacks its initial counted certificate");
    std::vector<d::Root> initial{cluster};
    // The exact product defining the input proves these other point roots.
    for (int j = 0; j < 8; ++j) {
        d::Root root;
        acb_set_si(root.box.value(), 4 + j);
        arb_set_si(acb_imagref(root.box.value()), 4);
        root.count = root.multiplicity = 1;
        root.publish(precision);
        initial.push_back(std::move(root));
    }
    require(d::disjoint(initial), "uncertain-value starting cover overlaps");
    cps_options options;
    options.goal = CPS_CLUSTER_RADIUS;
    options.cluster_radius = d::Number("1@-" + std::to_string(spacing_bits / 2));
    options.max_work_bits = precision;
    d::SearchStatistics statistics{};
    d::Subdivision search(initial, 1, precision, statistics);
    cps_stop_reason reason = CPS_STOP_NONE;
    // A test guard only; zero additional cells and fixed precision are the
    // resource contract. The goal remains far above the internal root gaps.
    for (size_t step = 0; step < 8 * (known.size() + 1) && search.pending() && reason == CPS_STOP_NONE; ++step)
        reason = search.step(options, enclose, 0, std::numeric_limits<double>::infinity(), accuracy);
    std::vector<d::Root> cover;
    require(reason == CPS_STOP_NONE && !search.pending() && search.complete_cover(cover),
            "an uncertain Newton value prevented counted contraction without new cells");
    std::vector<unsigned> visits(known.size());
    const auto epsilon = reference::power2(-spacing_bits / 2);
    for (const auto &root : cover) {
        const auto box = rectangle(root);
        uint64_t count = 0;
        for (size_t j = 0; j < known.size(); ++j)
            if (box.real_lower <= known[j].real && known[j].real <= box.real_upper &&
                box.imag_lower <= known[j].imag && known[j].imag <= box.imag_upper) {
                ++visits[j];
                ++count;
            }
        const auto real_width = box.real_upper - box.real_lower;
        const auto imag_width = box.imag_upper - box.imag_lower;
        require(count == root.count && real_width * real_width + imag_width * imag_width <= 4 * epsilon * epsilon,
                "counted contraction disagrees with planted roots or its absolute-radius goal");
    }
    require(std::all_of(visits.begin(), visits.end(), [](unsigned count) { return count == 1; }),
            "counted contraction loses or repeats a planted root");
}

void resolved_derivative_proposal() {
    // p(z) = z^2 - 2^(1-t) z - 1 has two real roots near +/-1, yet
    // its first Newton proposal at zero is -2^t. All coefficients and
    // that tiny nonzero derivative are exact at the starting precision.
    // Also rotate the variable to exercise the imaginary coordinate.
    constexpr uint64_t goal_bits = 53;
    constexpr uint64_t initial_bits = goal_bits + 32;
    for (int height : {512, 4096})
        for (bool imaginary : {false, true}) {
            const auto delta = reference::power2(-height);
            const reference::Polynomial polynomial{-1, -2 * delta, 1};
            const d::Encloser enclose = [&](acb_poly_t out, slong) {
                acb_poly_zero(out);
                acb_poly_set_coeff_si(out, 0, -1);
                acb_poly_set_coeff_si(out, 2, imaginary ? -1 : 1);
                d::ComplexBall linear;
                auto coordinate = imaginary ? acb_imagref(linear.value()) : acb_realref(linear.value());
                arb_set_si(coordinate, imaginary ? 1 : -1);
                arb_mul_2exp_si(coordinate, coordinate, 1 - height);
                acb_poly_set_coeff_acb(out, 1, linear.value());
            };
            d::ComplexPolynomial enclosed;
            enclose(enclosed.value(), initial_bits);
            d::ComplexBall center, value, derivative;
            acb_poly_evaluate2(value.value(), derivative.value(), enclosed.value(), center.value(), initial_bits);
            require(acb_is_exact(value.value()) && acb_is_exact(derivative.value()) &&
                        !acb_contains_zero(derivative.value()),
                    "resolved-derivative fixture lost its exact nonzero derivative");
            d::Float radius;
            arf_set_ui(radius.value(), 2);
            d::Root initial;
            require(d::certify_cluster_box(initial, enclosed.value(), center.value(), radius.value(), 2, 1, initial_bits),
                    "resolved-derivative fixture lacks its starting certificate");
            cps_options options;
            options.goal = CPS_ROOT_APPROXIMATION;
            options.bits = goal_bits;
            options.max_work_bits = 2 * static_cast<uint64_t>(height);
            d::SearchStatistics statistics{};
            d::Subdivision search({initial}, 1, initial_bits, statistics);
            for (unsigned turn = 0; turn < 2000 && search.pending(); ++turn) {
                require(search.step(options, enclose, 10000) == CPS_STOP_NONE,
                        "resolved-derivative fixture stopped before completing");
                require(search.max_bits() <= initial_bits && statistics.precision_increases == 0,
                        "an outside Newton proposal requested precision for its grid index");
            }
            std::vector<d::Root> cover;
            require(!search.pending() && search.complete_cover(cover) && cover.size() == 2,
                    "resolved-derivative fixture did not complete both roots");
            unsigned negative = 0, positive = 0;
            std::vector<reference::Rectangle> rectangles;
            for (const auto &root : cover) {
                const auto box = rectangle(root);
                const auto &lo = imaginary ? box.imag_lower : box.real_lower;
                const auto &hi = imaginary ? box.imag_upper : box.real_upper;
                require((imaginary ? box.real_lower <= 0 && box.real_upper >= 0
                                   : box.imag_lower <= 0 && box.imag_upper >= 0) &&
                            root.count == 1 && root.multiplicity == 1,
                        "resolved-derivative output has the wrong axis or count");
                const auto lower_value = reference::evaluate(polynomial, lo);
                const auto upper_value = reference::evaluate(polynomial, hi);
                if (hi < 0) {
                    require(lower_value >= 0 && upper_value <= 0, "negative root is not bracketed");
                    ++negative;
                } else {
                    require(lo > delta && lower_value <= 0 && upper_value >= 0, "positive root is not bracketed");
                    ++positive;
                }
                reference::check_accuracy(box, goal_bits);
                rectangles.push_back(box);
            }
            require(negative == 1 && positive == 1, "resolved-derivative fixture lost or repeated a root");
            reference::check_disjoint(rectangles);
        }
}

void stationary_samples(bool integral, bool complex_scalar, cps_goal goal, d::NewtonAccuracy accuracy = d::NewtonAccuracy::Legacy) {
    reference::Polynomial polynomial{1};
    for (const auto &a : squares)
        polynomial = reference::multiply(polynomial, {-a, 0, 1});
    if (integral)
        reference::make_primitive(polynomial);
    const auto derivative = reference::derivative(polynomial);
    for (int point : {0, 1}) {
        require(reference::evaluate(derivative, point) == 0, "fixture sample is not stationary");
        require(reference::evaluate(polynomial, point) != 0, "fixture sample is a root");
    }
    require(reference::gcd(polynomial, derivative).size() == 1, "fixture is not square-free");

    std::vector<d::Number> coefficients;
    for (const auto &coefficient : polynomial)
        coefficients.emplace_back(coefficient.str());
    const d::Encloser enclose = [&](acb_poly_t out, slong bits) {
        acb_poly_zero(out);
        d::ComplexBall coefficient;
        for (size_t j = 0; j < coefficients.size(); ++j) {
            coefficients[j].enclose(acb_realref(coefficient.value()), bits);
            if (complex_scalar)
                arb_set(acb_imagref(coefficient.value()), acb_realref(coefficient.value()));
            acb_poly_set_coeff_acb(out, static_cast<slong>(j), coefficient.value());
        }
    };
    d::ComplexPolynomial p;
    enclose(p.value(), 256);
    d::ComplexBall center;
    d::Float radius;
    arf_set_ui_2exp_si(radius.value(), 1, -2);
    d::Root inner;
    require(d::certify_cluster_box(inner, p.value(), center.value(), radius.value(), 2, 1, 256),
            "the inner stationary fixture has no initial count-two certificate");
    std::vector<d::Root> initial{inner};
    // The factorization independently proves these eight outer singleton
    // enclosures and the completeness of the starting cover.
    for (size_t j = 1; j < squares.size(); ++j) {
        d::RealBall a, root;
        d::Number(squares[j].str()).enclose(a.value(), 256);
        arb_sqrt(root.value(), a.value(), 256);
        for (int sign : {-1, 1}) {
            d::Root region;
            arb_mul_si(acb_realref(region.box.value()), root.value(), sign, 256);
            region.count = region.multiplicity = 1;
            region.publish(128);
            initial.push_back(std::move(region));
        }
    }
    require(d::disjoint(initial), "fixture starting cover overlaps");
    cps_options options;
    options.goal = goal;
    options.bits = 53;
    options.max_work_bits = 1 << 20;
    d::SearchStatistics statistics{};
    d::Subdivision search(initial, 1, 85, statistics);
    uint64_t peak_bits = 85;
    unsigned turns = 0;
    bool checked_symmetric_regrouping = false;
    for (; turns < 2000 && search.pending(); ++turns) {
        peak_bits = std::max(peak_bits, search.max_bits());
        // A test guard, not a production limit: fail before the old retry
        // loop allocates ever larger rational enclosures at unchanged geometry.
        require(peak_bits <= 4096, "stationary Newton samples keep raising precision without progress");
        require(search.step(options, enclose, 10000, std::numeric_limits<double>::infinity(), accuracy) == CPS_STOP_NONE, "stationary fixture stopped early");
        if (!checked_symmetric_regrouping && statistics.boxes_split > 0) {
            std::vector<d::Root> intermediate;
            if (search.complete_cover(intermediate))
                for (const auto &root : intermediate)
                    if (root.count == 2 && !acb_equal(root.box.value(), inner.box.value())) {
                        require(arf_is_zero(arb_midref(acb_realref(root.box.value()))) &&
                                    arf_is_zero(arb_midref(acb_imagref(root.box.value()))),
                                "symmetric regrouping acquired a rounded union center");
                        checked_symmetric_regrouping = true;
                    }
        }
    }
    require(checked_symmetric_regrouping, "stationary fixture did not exercise a new symmetric cluster hull");
    std::vector<d::Root> cover;
    require(!search.pending() && search.complete_cover(cover), "stationary fixture did not finish");
    uint64_t total = 0;
    std::vector<reference::Rectangle> rectangles;
    for (const auto &root : cover) {
        const auto r = rectangle(root);
        unsigned truth = 0;
        if (r.imag_lower <= 0 && r.imag_upper >= 0)
            for (const auto &a : squares) {
                truth += contains_positive_sqrt(r.real_lower, r.real_upper, a);
                truth += contains_positive_sqrt(-r.real_upper, -r.real_lower, a);
            }
        require(truth == root.count, "stationary fixture's output count disagrees with exact factors");
        require(root.count == 1, "stationary fixture did not separate the distinct roots");
        reference::check_accuracy(r, 53);
        rectangles.push_back(r);
        total += root.count;
    }
    reference::check_disjoint(rectangles);
    require(total == 10, "stationary fixture lost roots");
    std::cout << "stationary integral=" << integral << " complex=" << complex_scalar << " goal=" << goal
              << " accuracy_policy=" << static_cast<int>(accuracy) << " turns=" << turns << " peak_bits=" << peak_bits << " splits=" << statistics.boxes_split
              << " contractions=" << statistics.newton_successes << '\n';
}

void outer_certificate_precision() {
    // p = z^2 - 1/64 + 2^-32 z^16. Two small disjoint certified regions
    // merge after the constructor enlarges their support squares. The local
    // search is responsible for these two roots, not the other fourteen.
    constexpr slong degree = 16;
    constexpr slong initial_bits = 85;
    d::ComplexPolynomial exact;
    d::ComplexBall coefficient;
    acb_set_si(coefficient.value(), -1);
    acb_mul_2exp_si(coefficient.value(), coefficient.value(), -6);
    acb_poly_set_coeff_acb(exact.value(), 0, coefficient.value());
    acb_poly_set_coeff_si(exact.value(), 2, 1);
    acb_one(coefficient.value());
    acb_mul_2exp_si(coefficient.value(), coefficient.value(), -2 * degree);
    acb_poly_set_coeff_acb(exact.value(), degree, coefficient.value());

    // On each radius-1/32 circle about +/-1/8, the quadratic has modulus
    // at least (1/32)(1/4-1/32) = 7/1024. The extra term is at most
    // 2^-32 (5/32)^16 < 7/1024, so Rouche independently supplies one root.
    std::vector<d::Root> initial;
    d::Float radius;
    arf_set_ui_2exp_si(radius.value(), 1, -5);
    for (int sign : {-1, 1}) {
        d::ComplexBall center;
        acb_set_si(center.value(), sign);
        acb_mul_2exp_si(center.value(), center.value(), -3);
        d::Root root;
        require(d::certify_cluster_box(root, exact.value(), center.value(), radius.value(), 1, 1, 256),
                "outer-count fixture has no initial singleton certificate");
        initial.push_back(std::move(root));
    }
    require(d::disjoint(initial), "outer-count fixture's initial regions overlap");
    const d::Encloser enclose = [&](acb_poly_t out, slong bits) {
        acb_poly_set_round(out, exact.value(), bits);
        d::ComplexBall error;
        // A deliberately loose but convergent arithmetic enclosure: at the
        // first precision its z^15 error contributes 4 on the radius-2 disk,
        // while contributing only 2^-13 on the inner radius-1 disk.
        arb_add_error_2exp_si(acb_realref(error.value()), initial_bits - (degree - 3) - bits);
        acb_poly_set_coeff_acb(out, degree - 1, error.value());
    };
    d::ComplexPolynomial coarse;
    d::ComplexBall center;
    enclose(coarse.value(), initial_bits);
    arf_one(radius.value());
    const auto inner = d::pellet_count(coarse.value(), center.value(), radius.value(), initial_bits,
                                       d::subdivision_graeffe_steps(degree));
    require(inner.count == 2, "outer-count fixture has no proved inner count at its first precision");
    d::Root refused;
    bool needs_precision = false;
    require(!d::certify_cluster_box(refused, coarse.value(), center.value(), radius.value(), 2, 1,
                                    initial_bits, d::subdivision_graeffe_steps(degree), nullptr,
                                    &needs_precision) && needs_precision,
            "outer-count fixture's outer certificate does not request precision");
    d::SearchStatistics statistics{};
    d::Subdivision search(initial, 1, initial_bits, statistics);
    std::vector<d::Root> cover;
    require(!search.complete_cover(cover), "outer-count fixture did not merge into an uncertified component");
    cps_options options;
    options.goal = CPS_ROOT_APPROXIMATION;
    options.bits = 53;
    options.max_work_bits = 1024;
    require(search.step(options, enclose, 10000) == CPS_STOP_NONE, "outer-count fixture stopped early");
    require(search.min_bits() > initial_bits && statistics.precision_increases == 1 && statistics.boxes_split == 0,
            "a proved inner count discarded the outer certificate's precision request");
    require(search.step(options, enclose, 10000) == CPS_STOP_NONE && search.complete_cover(cover),
            "fresh precision did not recover the two-root certificate");
    require(cover.size() == 1 && cover.front().count == 2 && statistics.boxes_split == 0,
            "outer-certificate precision recovery changed the local count or split unnecessarily");
}
} // namespace

int main() {
    certpolysolve::ThreadScope lifetime;
    try {
        small_cell_budget();
        exact_points_release_support();
        finite_cap_split();
        ring_support_preserves_counts();
        resolved_derivative_proposal();
        preserve_useful_uncertain_sample();
        uncertain_value_sample();
        outer_certificate_precision();
        for (const auto goal : {CPS_ROOT_APPROXIMATION, CPS_ROOT_ACCURACY}) {
            stationary_samples(false, false, goal);
            stationary_samples(false, true, goal);
            stationary_samples(true, false, goal);
        }
        resolved_outside_sample_pair();
        annular_coordinate_precision();
        // Exact arithmetic enclosers below shrink their rounding width with P.
        // Warm-cover coverage is empirical; it does not assert native-start geometry.
        const auto exact = d::NewtonAccuracy::PreparedExact;
        resolved_outside_sample_pair(exact);
        annular_coordinate_precision(exact);
        uncertain_value_sample(exact);
        for (const auto goal : {CPS_ROOT_APPROXIMATION, CPS_ROOT_ACCURACY}) {
            stationary_samples(false, false, goal, exact);
            stationary_samples(false, true, goal, exact);
            stationary_samples(true, false, goal, exact);
        }
        declared_uncertainty_preserves_counted_contraction();
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
