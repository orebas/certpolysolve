#include "numeric/evaluate/structured_source.hpp"
#include "solve/policies.hpp"
#include "solve/local_model_solver.hpp"
// Certificates, continuation ownership and stopping in original coordinates.
#include "numeric/arithmetic/polynomial_arithmetic.hpp"
#include "internal.hpp"
#include "exact_winding_api.hpp"
#include "methods/secular/secular_refinement.hpp"
#include "numeric/exact/exact_complex.hpp"
#include "methods/piecewise/piecewise_roots.hpp"
#include "methods/real/real_refinement.hpp"
#include <cstring>
#include <iostream>
#include <optional>
#include <random>

using namespace reference;
namespace d = certpolysolve::detail;
// A solver snapshot (the internal result) as the header's Result.
static certpolysolve::Result wrap(const cps_result &snapshot) {
    return d::make_result(CPS_OK, d::published_result(snapshot));
}

struct OwnershipObserver {
    bool saw_local = false, reactivated = false;
    static void record(void *context, const d::TraceEvent &event) noexcept {
        auto &self = *static_cast<OwnershipObserver *>(context);
        if (std::strcmp(event.stage, "secular_packet") || std::strcmp(event.phase, "begin"))
            return;
        if (self.saw_local && event.active > 3)
            self.reactivated = true;
        if (event.active <= 3)
            self.saw_local = true;
    }
};

static void inspect(const cps_result &result, cps_goal goal, int bits, bool require_goal) {
    require(result.complete && result.degree == 100, "lost certified degree100 coverage");
    if (require_goal)
        require(result.goal_met, "repaired state did not finish its requested contract");
    const auto exported = wrap(result);
    require(cps_result_goal(exported.get()) == goal, "result changed the session goal");
    const auto &regions = exported.roots();
    check_disjoint(regions);
    uint64_t total = 0;
    for (const auto &region : regions) {
        total += region.count;
        if (result.goal_met && goal == CPS_ROOT_ACCURACY) {
            require(region.count == 1 && region.multiplicity.value_or(0) == 1,
                    "individual isolation retained an unresolved cluster");
            check_accuracy(region, bits);
        } else if (result.goal_met) {
            const Rational lo = endpoint(region.box.real.lo), hi = endpoint(region.box.real.hi);
            const Rational il = endpoint(region.box.imag.lo), ih = endpoint(region.box.imag.hi);
            const Rational re = lo > 0 ? lo : hi < 0 ? -hi : Rational(0);
            const Rational im = il > 0 ? il : ih < 0 ? -ih : Rational(0);
            require(((hi - lo) * (hi - lo) + (ih - il) * (ih - il)) * power2(2 * bits) <=
                        4 * (re * re + im * im), "tightened counted region failed exact accuracy");
        }
    }
    require(total == 100, "session cover count changed");
}

static void cover_ownership() {
        ComplexPolynomial coefficients(101);
        coefficients[0] = Complex(1);
        coefficients[1] = Complex(0, 300);
        coefficients[2] = Complex(-30000);
        coefficients[3] = Complex(0, -1000000);
        coefficients[100] = Complex(1);
        auto polynomial = input(coefficients);
        // Each a fresh solve (one advance a solve state since 18 October 2026; until then one state was advanced a
        // unit at a time, cancelled once local ownership began, and tightened): the whole solve under the
        // observer, where certified ownership is never reactivated; budgeted stops, each a complete cover; a
        // solve cancelled once local ownership begins; and 128 bits.
        for (const auto goal : {CPS_ROOT_ACCURACY, CPS_ROOT_APPROXIMATION}) {
            OwnershipObserver observer;
            auto policy = d::production_policy();
            policy.prepared_source_refinement = true;
            policy.trace = OwnershipObserver::record;
            policy.trace_data = &observer;
            cps_options options;
            options.goal = goal;
            options.bits = 53;
            options.max_work_bits = 8192;
            options.work_units = 100000;
            {
                d::SolveState state(certpolysolve::detail::dense(*polynomial.get()), policy);
                const auto result = state.advance(&options);
                inspect(result, goal, 53, true);
                require(!observer.reactivated, "an inferior cover reactivated certified roots");
                if (goal == CPS_ROOT_ACCURACY) require(observer.saw_local, "control missed local ownership");
            }
            for (const uint64_t units : {1, 2, 4, 8, 16, 32, 64}) {
                observer = {};
                options.work_units = units;
                d::SolveState state(certpolysolve::detail::dense(*polynomial.get()), policy);
                const auto stopped = state.advance(&options);
                inspect(stopped, goal, 53, false);
                require(!observer.reactivated, "an inferior cover reactivated certified roots");
                if (!stopped.goal_met) require(stopped.termination == CPS_WORK_LIMIT, "unexpected refinement stop");
            }
            {
                observer = {};
                options.work_units = 100000;
                options.cancel = [](void *data) -> int { return static_cast<OwnershipObserver *>(data)->saw_local; };
                options.cancel_data = &observer;
                d::SolveState state(certpolysolve::detail::dense(*polynomial.get()), policy);
                const auto stopped = state.advance(&options);
                inspect(stopped, goal, 53, false);
                if (goal == CPS_ROOT_ACCURACY)
                    require(stopped.termination == CPS_CANCELLED, "a cancellation during local ownership was not kept");
                options.cancel = nullptr;
                options.cancel_data = nullptr;
            }
            observer = {};
            options.bits = 128;
            d::SolveState tighter(certpolysolve::detail::dense(*polynomial.get()), policy);
            const auto tightened = tighter.advance(&options);
            inspect(tightened, goal, 128, true);
            require(!observer.reactivated, "a 128-bit solve mixed certificate and ownership states");
            std::cout << "goal=" << goal << " fresh solves: whole, budgeted, cancelled, 128 bits\n";
        }
}

static d::GaussianRational gaussian(const Complex &value) {
    d::GaussianRational result;
    fmpq_set_str(result.real.value(), value.real.str().c_str(), 10);
    fmpq_set_str(result.imag.value(), value.imag.str().c_str(), 10);
    return result;
}

static void enclose(acb_t out, const Complex &value, slong bits) {
    const auto exact = gaussian(value);
    arb_set_fmpq(acb_realref(out), exact.real.value(), bits);
    arb_set_fmpq(acb_imagref(out), exact.imag.value(), bits);
}

static std::unique_ptr<d::PolynomialSource> source(const std::vector<Complex> &roots) {
    d::GaussianPolynomial exact;
    for (const auto &coefficient : from_roots(roots))
        exact.push_back(gaussian(coefficient));
    d::PolynomialSourceOptions options;
    options.input = d::InputSemantics::Exact;
    const bool real = std::all_of(roots.begin(), roots.end(),
                                 [](const Complex &root) { return root.imag == 0; });
    return std::make_unique<d::PolynomialSource>(
        [exact](acb_poly_t out, slong bits) { d::enclose_gaussian(out, exact, bits); },
        static_cast<slong>(roots.size()), real, options);
}

static void check_known_cover(const std::vector<d::Root> &cover,
                             const std::vector<Complex> &known, const cps_options &options) {
    require(d::disjoint(cover), "translated output regions overlap");
    size_t total = 0;
    for (const auto &root : cover) {
        const Rational lo = power2(std::stoi(root.endpoints[0].exponent_text())) * Rational(root.endpoints[0].mantissa_text());
        const Rational hi = power2(std::stoi(root.endpoints[1].exponent_text())) * Rational(root.endpoints[1].mantissa_text());
        const Rational il = power2(std::stoi(root.endpoints[2].exponent_text())) * Rational(root.endpoints[2].mantissa_text());
        const Rational ih = power2(std::stoi(root.endpoints[3].exponent_text())) * Rational(root.endpoints[3].mantissa_text());
        size_t count = 0;
        for (const auto &value : known)
            if (lo <= value.real && value.real <= hi && il <= value.imag && value.imag <= ih)
                ++count;
        require(count == root.count, "translation lost an independently known root");
        total += count;
        const Rational re = lo > 0 ? lo : hi < 0 ? -hi : Rational(0);
        const Rational im = il > 0 ? il : ih < 0 ? -ih : Rational(0);
        const Rational squared_radius = ((hi - lo) * (hi - lo) + (ih - il) * (ih - il)) / 4;
        require(squared_radius * power2(2 * static_cast<int>(options.bits)) <= re * re + im * im,
                "translated enclosure fails exact true-relative accuracy");
        if (options.goal == CPS_ROOT_ACCURACY)
            require(root.count == root.multiplicity, "isolation accepted an unresolved group");
        if (root.realness == CPS_REAL)
            require(std::any_of(known.begin(), known.end(), [&](const Complex &value) {
                return value.imag == 0 && lo <= value.real && value.real <= hi;
            }), "translation preserved invalid local realness");
    }
    require(total == known.size(), "translated cover is incomplete");
}

static void coordinate_stopping() {
    // A local coordinate near one can already have 53 relative bits while
    // its original coordinate near 2^-100 still needs another hundred bits.
    const std::vector<Complex> known = {Complex(power2(-100)), Complex(3 * power2(-100))};
    for (const auto &offset : {Complex(-1), Complex(Rational(-1) / 3, Rational(2) / 5),
                               Complex(-power2(160), power2(120))}) {
        std::vector<Complex> local;
        for (const auto &root : known)
            local.push_back(root - offset);
        const d::RootCoordinates coordinates(gaussian(offset));
        std::vector<d::ComplexBall> nodes(local.size()), residues;
        d::RootVector seeds(static_cast<slong>(local.size()));
        for (size_t j = 0; j < local.size(); ++j) {
            enclose(nodes[j].value(), local[j] + Complex(power2(-130)), 512);
            acb_get_mid(nodes[j].value(), nodes[j].value());
            acb_set(seeds.data() + j, nodes[j].value());
        }
        for (const auto goal : {CPS_ROOT_ACCURACY, CPS_ROOT_APPROXIMATION}) {
            cps_options options;
            options.goal = goal;
            options.bits = 53;
            options.max_work_bits = 4096;
            certpolysolve::detail::CandidateStatistics statistics{};
            certpolysolve::detail::SearchStatistics search{};
            uint64_t maximum = 0;
            auto polynomial = source(local);
            std::vector<slong> bits(local.size(), 512);
            require(d::regenerate_residues(residues, *polynomial, nodes, {}, {}, bits,
                        options, statistics, maximum), "coordinate regression failed regeneration");
            std::vector<d::Root> cover;
            std::vector<uint8_t> inactive;
            require(d::secular_inclusion_cover(cover, inactive, nodes, residues, polynomial->real(),
                        1, static_cast<slong>(maximum), options), "local test cover failed");
            require(std::all_of(inactive.begin(), inactive.end(), [](uint8_t x) { return x != 0; }),
                    "control does not reproduce premature local completion");
            require(d::secular_inclusion_cover(cover, inactive, nodes, residues, polynomial->real(),
                        1, static_cast<slong>(maximum), options, coordinates), "mapped test cover failed");
            require(std::none_of(inactive.begin(), inactive.end(), [](uint8_t x) { return x != 0; }),
                    "local-relative stopping froze an inaccurate original-coordinate root");
            auto absolute = options;
            absolute.absolute_tolerance = d::Number(power2(-110).str());
            require(coordinates.meets(cover.front(), absolute, static_cast<slong>(maximum)),
                    "translation ignored the user's explicit absolute tolerance");

            auto policy = d::production_policy();
            d::SecularRefinement refinement(source(local), seeds.data(), 1, 512, policy, d::solve_local_model, coordinates);
            auto finish = [&] {
                for (unsigned step = 0; step < 64; ++step) {
                    const auto reason = refinement.step(options, search, statistics);
                    if (refinement.cover(cover) && std::all_of(cover.begin(), cover.end(), [&](const d::Root &root) {
                            return d::meets_goal(root, options, static_cast<slong>(refinement.precision()));
                        }))
                        return;
                    require(reason == CPS_STOP_NONE, "original-coordinate refinement stopped prematurely");
                }
                throw std::runtime_error("original-coordinate refinement exhausted its work budget");
            };
            finish();
            check_known_cover(cover, known, options);
            const auto held = cover;
            options.bits = 96;
            options.cancel = [](void *) { return 1; };
            require(refinement.step(options, search, statistics) == CPS_STOP_CANCELLED,
                    "coordinate refinement ignored cancellation");
            options.cancel = nullptr;
            finish();
            check_known_cover(cover, known, options);
            options.bits = 53;
            check_known_cover(held, known, options);
        }
    }
}

static void frozen_duplicates() {
    std::vector<d::ComplexBall> nodes(3);
    std::vector<slong> bits(3, 128);
    for (const auto frozen : {size_t(0), size_t(2)}) {
        std::vector<uint8_t> inactive(3, 0);
        inactive[frozen] = 1;
        for (auto &node : nodes)
            acb_one(node.value());
        const auto held = nodes[frozen];
        d::separate_interpolation_nodes(nodes, inactive, bits);
        require(acb_equal(nodes[frozen].value(), held.value()), "duplicate repair moved a certified point");
        for (size_t i = 0; i < nodes.size(); ++i) {
            require(acb_is_exact(nodes[i].value()), "duplicate repair produced inexact interpolation nodes");
            for (size_t j = 0; j < i; ++j)
                require(!acb_equal(nodes[i].value(), nodes[j].value()), "duplicate repair left equal nodes");
        }
    }
}

// The pairwise scan that separate_interpolation_nodes replaced (O(n^2)), kept
// as the reference for the sorted one.
static void pairwise_separation(std::vector<d::ComplexBall> &proposals, const std::vector<uint8_t> &inactive,
                                const std::vector<slong> &iteration_bits) {
    for (size_t i = 0; i < proposals.size(); ++i)
        for (size_t j = 0; j < i; ++j)
            if (acb_equal(proposals[i].value(), proposals[j].value())) {
                const auto moving = inactive[i] ? j : i;
                if (inactive[moving])
                    continue;
                d::RealBall shift;
                acb_abs(shift.value(), proposals[moving].value(), 64);
                arb_add_ui(shift.value(), shift.value(), 1, 64);
                arb_get_mid_arb(shift.value(), shift.value());
                arb_mul_2exp_si(shift.value(), shift.value(), -iteration_bits[moving] + 8);
                arb_mul_ui(shift.value(), shift.value(), static_cast<ulong>(moving + 1), iteration_bits[moving]);
                arb_add(acb_imagref(proposals[moving].value()), acb_imagref(proposals[moving].value()),
                        shift.value(), iteration_bits[moving]);
                acb_get_mid(proposals[moving].value(), proposals[moving].value());
            }
}

// separate_interpolation_nodes finds equal proposals by sorting: it moves the
// same proposals by the same amounts as the pairwise scan, on random groups of
// equal proposals with frozen members and mixed iteration precisions.
static void duplicate_groups_match_pairwise() {
    std::mt19937_64 random(20261006);
    for (int trial = 0; trial < 300; ++trial) {
        const size_t n = 2 + static_cast<size_t>(random() % 40);
        std::vector<d::ComplexBall> proposals(n);
        std::vector<uint8_t> inactive(n);
        std::vector<slong> bits(n);
        for (size_t i = 0; i < n; ++i) {
            // Values from a set of 18, so that groups of equal proposals form.
            acb_set_si_si(proposals[i].value(), static_cast<slong>(random() % 6) - 3,
                          static_cast<slong>(random() % 3));
            inactive[i] = random() % 4 == 0;
            bits[i] = 64 + static_cast<slong>(random() % 4) * 32;
        }
        auto expected = proposals;
        pairwise_separation(expected, inactive, bits);
        d::separate_interpolation_nodes(proposals, inactive, bits);
        for (size_t i = 0; i < n; ++i)
            require(acb_equal(proposals[i].value(), expected[i].value()),
                    "the sorted duplicate scan moved a proposal differently from the pairwise scan");
    }
}

struct RestartObserver {
    const char *cancel_stage = nullptr;
    bool at_ceiling = false, requested = false, cancelled = false;
    uint64_t ceiling = 0;
    static void trace(void *context, const d::TraceEvent &event) noexcept {
        auto &self = *static_cast<RestartObserver *>(context);
        if (std::strcmp(event.phase, "begin")) return;
        if (!std::strcmp(event.stage, "secular_candidate_restart"))
            self.at_ceiling = event.precision == self.ceiling;
        if (self.cancel_stage && !self.cancelled && self.at_ceiling &&
            !std::strcmp(event.stage, self.cancel_stage))
            self.requested = self.cancelled = true;
    }
};

static void cancelled_restart() {
    for (const auto [exponent, ceiling] : {std::pair{20, 64u}, std::pair{10, 256u}})
        for (const char *stage : {static_cast<const char *>(nullptr), "secular_candidate_restart",
                                  "component_initialization"}) {
            const std::vector<Complex> roots{Complex(1, power2(-exponent)), Complex(1, -power2(-exponent))};
            d::RootVector seeds(2);
            acb_zero(seeds[0]);
            acb_set_si(seeds[1], 2);
            RestartObserver observer{stage, false, false, false, ceiling};
            d::EnginePolicy policy;
            policy.trace = RestartObserver::trace;
            policy.trace_data = &observer;
            d::SecularRefinement state(source(roots), seeds.data(), 1, 32, policy, d::solve_local_model);
            cps_options options;
            options.max_work_bits = ceiling;
            options.cancel = [](void *context) -> int { return static_cast<RestartObserver *>(context)->requested; };
            options.cancel_data = &observer;
            certpolysolve::detail::SearchStatistics search{};
            certpolysolve::detail::CandidateStatistics statistics{};
            bool resuming = false;
            uint64_t raises_at_cancel = 0;
            for (unsigned call = 0; call < 256 && state.active(); ++call) {
                const auto stop = state.step(options, search, statistics);
                if (resuming)
                    require(search.precision_increases == raises_at_cancel && stop != CPS_STOP_PRECISION,
                            "cancelled secular restart consumed a second precision increase");
                resuming = stop == CPS_STOP_CANCELLED;
                if (resuming) {
                    observer.requested = false;
                    raises_at_cancel = search.precision_increases;
                } else require(stop == CPS_STOP_NONE, "secular restart control stopped unexpectedly");
            }
            require(!stage || observer.cancelled, "restart cancellation witness missed its boundary");
            std::vector<d::Root> cover;
            require(!state.active() && state.cover(cover), "cancelled restart failed to resume at its ceiling");
            check_known_cover(cover, roots, options);
        }
}

static void unsupported_query_preserves_state() {
    const std::vector<Complex> known{Complex(1), Complex(2)};
    const auto input = reference::input(from_roots(known));
    d::PiecewiseRootState state(certpolysolve::detail::dense(*input.get()));
    cps_options options;
    options.goal = CPS_ROOT_APPROXIMATION;
    const auto completed = state.advance(options);
    require(completed.complete && completed.goal_met && completed.roots.size() == known.size(),
            "unsupported-query control did not finish");
    for (bool real_only : {false, true}) {
        auto unsupported = options;
        unsupported.bits = 160;
        unsupported.real_only = real_only;
        if (!real_only) {
            unsupported.goal = CPS_CLUSTER_RADIUS;
            unsupported.cluster_radius = d::Number("1@-200");
        }
        const auto rejected = state.advance(unsupported);
        require(rejected.termination == CPS_UNSUPPORTED_QUERY && !rejected.complete && !rejected.goal_met,
                "unsupported piecewise query retained successful flags after mutating its cover");
        auto resume = options;
        resume.work_units = 0;
        const auto unchanged = state.advance(resume);
        require(unchanged.complete && unchanged.goal_met && unchanged.roots.size() == completed.roots.size() &&
                    unchanged.statistics_v2.work_units == completed.statistics_v2.work_units,
                "unsupported query consumed the retained piecewise state");
        for (size_t i = 0; i < completed.roots.size(); ++i) {
            const auto &root = unchanged.roots[i];
            const auto &old = completed.roots[i];
            for (size_t j = 0; j < root.endpoints.size(); ++j)
                require(root.endpoints[j].mantissa_text() == old.endpoints[j].mantissa_text() &&
                            root.endpoints[j].exponent_text() == old.endpoints[j].exponent_text(),
                        "unsupported query changed a retained enclosure");
        }
        check_known_cover(unchanged.roots, known, options);
        check_known_cover(completed.roots, known, options);
    }
}

static void regeneration_precision_memory() {
    std::vector<Complex> roots;
    for(int j=0;j<8;++j)roots.emplace_back(1+j*power2(-30),1+(j*j+1)*power2(-32));
    auto polynomial=source(roots);
    std::vector<d::ComplexBall> nodes(roots.size()),first,second;
    std::vector<slong> bits(roots.size(),64),hints(roots.size(),0);
    for(size_t j=0;j<roots.size();++j)enclose(nodes[j].value(),roots[j]+Complex(power2(-35),power2(-36)),512);
    d::EnginePolicy policy;policy.model_stopping=d::SecularModelPolicy::Continuation;
    policy.regeneration.value_only=true;policy.regeneration.reuse_products=true;
    policy.regeneration.remember_precision=true;policy.study_statistics=std::make_shared<d::StudyStatistics>();
    cps_options options;options.max_work_bits=2048;
    certpolysolve::detail::CandidateStatistics statistics{};uint64_t maximum=0;
    require(d::regenerate_residues(first,*polynomial,nodes,{},{},bits,options,statistics,maximum,&policy,false,hints),
            "initial remembered regeneration failed");
    const auto initial_evaluations=statistics.evaluations;
    require(*std::max_element(hints.begin(),hints.end())>96,"cancellation fixture did not need extra source precision");
    require(d::regenerate_residues(second,*polynomial,nodes,{},{},bits,options,statistics,maximum,&policy,false,hints),
            "remembered regeneration failed");
    require(statistics.evaluations-initial_evaluations<initial_evaluations,
            "successful source precision did not avoid repeated retries");
    for(size_t i=0;i<nodes.size();++i) {
        const auto z=roots[i]+Complex(power2(-35),power2(-36));
        Complex value(1),product(1);
        for(size_t j=0;j<roots.size();++j) {
            value=value*(z-roots[j]);
            if(i!=j)product=product*(z-(roots[j]+Complex(power2(-35),power2(-36))));
        }
        const auto denominator=product.real*product.real+product.imag*product.imag;
        const Complex quotient(-(value.real*product.real+value.imag*product.imag)/denominator,
                               -(value.imag*product.real-value.real*product.imag)/denominator);
        d::ComplexBall expected;enclose(expected.value(),quotient,1024);
        require(acb_contains(first[i].value(),expected.value()) && acb_contains(second[i].value(),expected.value()),
                "precision memory lost independently computed residue identity");
    }
    const auto saved=hints;
    options.cancel=[](void *)->int{return 1;};
    require(!d::regenerate_residues(second,*polynomial,nodes,{},{},bits,options,statistics,maximum,&policy,false,hints) && hints==saved,
            "cancelled regeneration committed precision hints");
}

// One regeneration of the clustered fixture above (differences of 2^-30 at nodes near 1, so that a product
// that lost the low bits of a difference would miss, and values that cancel some 200 bits), with the nodes moved
// by `offset`: its residues enclose the exact ones, and the two integer kernels served what `products` and
// `evaluations` say (every product, every evaluation, or none).
static void clustered_regeneration(const Complex &offset,bool integer_products,bool limb_evaluations,bool products,
                                   bool evaluations,const std::string &what) {
    std::vector<Complex> roots;
    for(int j=0;j<8;++j)roots.emplace_back(1+j*power2(-30),1+(j*j+1)*power2(-32));
    auto polynomial=source(roots);
    std::vector<d::ComplexBall> nodes(roots.size()),residues;
    std::vector<slong> bits(roots.size(),64),hints(roots.size(),0);
    for(size_t j=0;j<roots.size();++j)enclose(nodes[j].value(),roots[j]+offset,512);
    d::EnginePolicy policy;policy.model_stopping=d::SecularModelPolicy::Continuation;
    policy.regeneration.value_only=true;policy.regeneration.reuse_products=true;
    policy.regeneration.remember_precision=true;policy.regeneration.integer_node_products=integer_products;
    policy.regeneration.limb_horner=limb_evaluations;
    policy.study_statistics=std::make_shared<d::StudyStatistics>();
    cps_options options;options.max_work_bits=2048;
    certpolysolve::detail::CandidateStatistics statistics{};uint64_t maximum=0;
    d::RegenerationWork work;
    require(d::regenerate_residues(residues,*polynomial,nodes,{},{},bits,options,statistics,maximum,&policy,false,hints,&work),
            (what+": the regeneration failed").c_str());
    const auto &counts=*policy.study_statistics;
    require(counts.integer_node_products==(products?roots.size():0) && counts.product_builds>=roots.size(),
            (what+": "+std::to_string(counts.integer_node_products)+" products came from integers").c_str());
#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
    require(work.short_product_factors == (products ? roots.size()*(roots.size()-1) : 0),
            "integer node products were charged at a ball precision");
    require((work.limb_horner_work > 0) == evaluations,
            "limb Horner work did not follow the executed evaluator");
    if (!products && !evaluations)
        require(work.service_operations == work.operations, "ball regeneration proxy changed");
#endif
    require(statistics.evaluations>=roots.size() &&
                counts.limb_horner_evaluations==(evaluations?statistics.evaluations:0),
            (what+": "+std::to_string(counts.limb_horner_evaluations)+" of "+std::to_string(statistics.evaluations)+
             " evaluations came from limbs").c_str());
    for(size_t i=0;i<nodes.size();++i) {
        const auto z=roots[i]+offset;
        Complex value(1),product(1);
        for(size_t j=0;j<roots.size();++j) {
            value=value*(z-roots[j]);
            if(i!=j)product=product*(z-(roots[j]+offset));
        }
        const auto denominator=product.real*product.real+product.imag*product.imag;
        const Complex quotient(-(value.real*product.real+value.imag*product.imag)/denominator,
                               -(value.imag*product.real-value.real*product.imag)/denominator);
        d::ComplexBall exact;enclose(exact.value(),quotient,1024);
        require(acb_contains(residues[i].value(),exact.value()),
                (what+": a residue does not enclose the exact one").c_str());
    }
}
// The products of binary64 nodes are formed in integer mantissas (node_products.hpp), one a node, and the
// source's values at nodes of at most two limbs by Horner's rule on limbs (limb_horner_evaluator.hpp), every
// one; the residues enclose the exact ones with either, both or neither. Nodes of up to two limbs keep both
// kernels (the products from 256-bit differences), and nodes of more than 128 bits send every product and every
// evaluation back to ball arithmetic.
static void regeneration_takes_the_integer_kernels() {
    const Complex binary64(power2(-35),power2(-36)),longer(power2(-35)+power2(-70),power2(-36)),
        longest(power2(-35)+power2(-140),power2(-36));
    clustered_regeneration(binary64,true,false,true,false,"binary64 nodes, integer products");
    clustered_regeneration(binary64,false,true,false,true,"binary64 nodes, limb evaluations");
    clustered_regeneration(binary64,true,true,true,true,"binary64 nodes, both kernels");
    clustered_regeneration(binary64,false,false,false,false,"binary64 nodes, ball arithmetic");
    clustered_regeneration(longer,true,true,true,true,"nodes of 71 bits");
    clustered_regeneration(longest,true,true,false,false,"nodes of more than 128 bits");
}

// A partial regeneration of the clustered fixture: some of the eight nodes move, and each unmoved node's residue
// is updated by the moved nodes' factors, prod (b_i - old b_j) / (b_i - b_j). Between binary64 nodes the two
// products come from integer mantissas (ShortNodeProducts::product_over), one a node, however many nodes
// moved; with the kernel off, or from nodes of more than two limbs, the factors are taken one by one in ball
// arithmetic where at most a quarter moved, and the residue is rebuilt where more did. Only the nodes that are
// not updated are evaluated again, and every residue encloses the exact one of the new model.
static void regeneration_updates_unmoved_nodes() {
    std::vector<Complex> roots;
    for(int j=0;j<8;++j)roots.emplace_back(1+j*power2(-30),1+(j*j+1)*power2(-32));
    auto polynomial=source(roots);
    const auto exact_residue=[&](const std::vector<Complex> &at,size_t i) {
        Complex value(1),product(1);
        for(size_t j=0;j<roots.size();++j) {
            value=value*(at[i]-roots[j]);
            if(i!=j)product=product*(at[i]-at[j]);
        }
        const auto denominator=product.real*product.real+product.imag*product.imag;
        return Complex(-(value.real*product.real+value.imag*product.imag)/denominator,
                       -(value.imag*product.real-value.real*product.imag)/denominator);
    };
    // `updates`: the unmoved nodes updated from integers; `rebuilt`: the nodes evaluated again.
    struct Case { const char *what; Complex offset; bool kernels; std::vector<size_t> moved; uint64_t updates,rebuilt; };
    const Complex binary64(power2(-35),power2(-36)),longer(power2(-35)+power2(-70),power2(-36)),
        longest(power2(-35)+power2(-140),power2(-36));
    const Case cases[]={{"binary64 nodes, a quarter moved",binary64,true,{2,5},6,2},
                        {"binary64 nodes, a quarter moved, ball arithmetic",binary64,false,{2,5},0,2},
                        {"nodes of 71 bits, a quarter moved",longer,true,{2,5},6,2},
                        {"nodes of 141 bits, a quarter moved",longest,true,{2,5},0,2},
                        {"binary64 nodes, five moved",binary64,true,{0,2,3,5,6},3,5},
                        {"binary64 nodes, five moved, ball arithmetic",binary64,false,{0,2,3,5,6},0,8},
                        {"binary64 nodes, all but one moved",binary64,true,{0,1,2,3,4,5,6},1,7}};
    for(const auto &c:cases) {
        const std::string what=c.what;
        std::vector<Complex> before,after;
        for(const auto &root:roots)before.push_back(root+c.offset);
        after=before;
        for(const size_t j:c.moved)
            after[j]=after[j]+(j%2?Complex(Rational(0),-power2(-41)):Complex(power2(-40),Rational(0)));
        std::vector<d::ComplexBall> old_nodes(roots.size()),nodes(roots.size()),old_residues,residues;
        for(size_t j=0;j<roots.size();++j) {
            enclose(old_nodes[j].value(),before[j],512);
            enclose(nodes[j].value(),after[j],512);
        }
        std::vector<slong> bits(roots.size(),64),hints(roots.size(),0);
        d::EnginePolicy policy;policy.model_stopping=d::SecularModelPolicy::Continuation;
        policy.regeneration.value_only=true;policy.regeneration.reuse_products=true;
        policy.regeneration.remember_precision=true;
        policy.regeneration.integer_node_products=policy.regeneration.limb_horner=c.kernels;
        policy.study_statistics=std::make_shared<d::StudyStatistics>();
        cps_options options;options.max_work_bits=2048;
        certpolysolve::detail::CandidateStatistics statistics{};uint64_t maximum=0;
        require(d::regenerate_residues(old_residues,*polynomial,old_nodes,{},{},bits,options,statistics,maximum,&policy,false,hints),
                (what+": the first regeneration failed").c_str());
        const auto evaluated=statistics.evaluations;
        require(d::regenerate_residues(residues,*polynomial,nodes,old_nodes,old_residues,bits,options,statistics,maximum,&policy,
                                       true,hints),
                (what+": the partial regeneration failed").c_str());
        require(policy.study_statistics->integer_node_updates==c.updates,
                (what+": "+std::to_string(policy.study_statistics->integer_node_updates)+" updates came from integers, not "+
                 std::to_string(c.updates)).c_str());
        // Each rebuilt node is evaluated once or, where the fixture's cancellation asks, once more at twice the
        // precision; an updated node not at all.
        const auto again=statistics.evaluations-evaluated;
        require(again>=c.rebuilt && again<=2*c.rebuilt,
                (what+": "+std::to_string(again)+" evaluations for "+std::to_string(c.rebuilt)+" rebuilt nodes").c_str());
        for(size_t i=0;i<nodes.size();++i) {
            d::ComplexBall exact;enclose(exact.value(),exact_residue(after,i),1024);
            require(acb_contains(residues[i].value(),exact.value()),
                    (what+": residue "+std::to_string(i)+" of the partial regeneration does not enclose the exact one").c_str());
        }
    }
}

// A local packet that is finite but has not converged still supplies restart
// hints; only a measured deficiency of the local model (non-finite, or a
// leading coefficient containing zero) may raise the working precision
// (principled audit, fix 5). Class F: near_multiple-719 rotated by (3+4i)/5
// doubled from 1,872 to 65,536 bits because a 32-iteration local solve did
// not report convergence, which no precision can change.
static void counted_scale_hints_respect_the_final_request() {
    std::vector<Complex> roots;
    for (int re : {-1, 1})
        for (int im : {-1, 1})
            roots.emplace_back(1 + re * power2(-128), im * power2(-128));
    auto polynomial = source(roots);
    d::Root component;
    acb_one(component.box.value());
    arb_add_error_2exp_si(acb_realref(component.box.value()), -5);
    arb_add_error_2exp_si(acb_imagref(component.box.value()), -5);
    component.count = 4;
    d::EnginePolicy policy;
    cps_options options;
    options.goal = CPS_ROOT_APPROXIMATION;
    options.max_work_bits = 4096;
    const auto seeds = [&](d::RootVector &out, slong final_resolution, bool &unresolved) {
        d::CandidateStatistics statistics{};
        slong working = 2048;
        require(d::component_seeds(out.data(), *polynomial, component, 4, working, options,
                                  statistics, policy, std::nullopt, 64, &unresolved, final_resolution),
                "the counted-hint model did not supply finite hints");
        require(working == 2048, "the determined model unexpectedly changed precision");
        return statistics.candidate_sweeps;
    };
    d::RootVector control(4), counted(4), fine(4), accurate(4);
    bool unresolved = false;
    const auto control_sweeps = seeds(control, 0, unresolved);
    require(control_sweeps > 0, "the control no longer exercises local refinement");
    require(seeds(counted, 64, unresolved) == 0 && unresolved,
            "a counted request refined scale hints already within its resolution");
    for (slong j = 0; j < 4; ++j)
        require(acb_is_finite(counted[j]) && acb_contains(component.box.value(), counted[j]),
                "a counted hint escaped the retained component");
    require(seeds(fine, 256, unresolved) == control_sweeps,
            "the temporary resolution waived refinement needed by the final request");
    options.goal = CPS_ROOT_ACCURACY;
    require(seeds(accurate, 64, unresolved) == control_sweeps,
            "individual accuracy skipped local refinement");
    for (slong j = 0; j < 4; ++j)
        require(acb_equal(control[j], fine[j]) && acb_equal(control[j], accurate[j]),
                "an ineligible request changed its local hints");
}

static void counted_restart_keeps_the_callers_resolution() {
    // KAM4's close pairs fit the temporary 53-bit stage, but not these final
    // requests. Using the stage skipped useful local solves and later raised
    // precision from 730 to 1488, and from 6708 to 11008 bits, respectively.
    ComplexPolynomial polynomial(15);
    const Rational a("1000000000000000000000000");
    polynomial[0] = {4};
    polynomial[4] = {-4 * a};
    polynomial[7] = {4};
    polynomial[8] = {a * a};
    polynomial[11] = {2 * a};
    polynomial[14] = {1};
    const auto exact = input(polynomial);
    for (const uint64_t bits : {uint64_t{333}, uint64_t{3322}}) {
        cps_options options;
        options.goal = CPS_ROOT_APPROXIMATION;
        options.bits = bits;
        options.time_limit = 20;
        d::SolveState state(d::dense(*exact.get()), d::production_policy());
        const auto result = state.advance(&options);
        require(result.complete && result.goal_met, "the final-request regression did not finish");
        // A loose regression budget between the observed good and bad paths,
        // not a promised precision bound for arbitrary inputs.
        require(result.work_bits <= 3 * bits,
                "a staged counted request discarded useful fine-resolution hints");
    }
}

static void local_seeds_do_not_escalate() {
    // Eight roots in one component at 1: four 2^-100 apart inside four 2^-8
    // apart. Separating the inner four needs about 4 x 100 bits, so at 256
    // bits the local packet stays finite but cannot converge: exactly the
    // restart hint secular needs, which only the original source certifies.
    // 72 far roots on the integer grid make the degree exceed
    // small_complex_degree, which exempts small sources from escalation.
    std::vector<Complex> roots;
    for (int x = -4; x <= 4; ++x)
        for (int y = -4; y <= 4; ++y)
            if (std::abs(x - 1) > 1 || std::abs(y) > 1)
                roots.emplace_back(x, y);
    for (int re : {-1, 1})
        for (int im : {-1, 1}) {
            roots.emplace_back(1 + re * power2(-100), im * power2(-100));
            roots.emplace_back(1 + re * power2(-8), im * power2(-8));
        }
    auto polynomial = source(roots);
    d::Root component;
    acb_one(component.box.value());
    arb_add_error_2exp_si(acb_realref(component.box.value()), -5);
    arb_add_error_2exp_si(acb_imagref(component.box.value()), -5);
    component.count = 8;
    cps_options options;
    options.max_work_bits = 65536;
    certpolysolve::detail::CandidateStatistics statistics{};
    d::EnginePolicy policy;
    d::RootVector hints(8);
    slong working = 256;
    const bool seeded = d::component_seeds(hints.data(), *polynomial, component, 8, working, options, statistics,
                                           policy, std::nullopt);
    require(working == 256, "local seeds raised the precision without a measured deficiency");
    require(seeded, "a finite local packet did not supply restart hints");
}

// A restart's centre steps are Newton steps on p^(k-1): they read the
// coefficients k - 1 and k of the source's Taylor expansion about the
// component, which are that derivative's value and slope there
// (PolynomialSource::taylor_pair: one joint evaluation, no expansion); only
// the local model takes a whole shift (component_seeds). Three whole shifts a
// restart for the steps took 27 of near_multiple-769's 35 s
// (measured on 8 October 2026).
static void restart_steps_take_no_expansion() {
    std::vector<Complex> roots;
    for (int x = -4; x <= 4; ++x)
        for (int y = -4; y <= 4; ++y)
            if (std::abs(x - 1) > 1 || std::abs(y) > 1)
                roots.emplace_back(x, y);
    for (int re : {-1, 1})
        for (int im : {-1, 1}) {
            roots.emplace_back(1 + re * power2(-100), im * power2(-100));
            roots.emplace_back(1 + re * power2(-8), im * power2(-8));
        }
    auto polynomial = source(roots);
    const slong working = 256;
    d::ComplexBall center;
    acb_set_d_d(center.value(), 1 + 0x1p-20, 0x1p-30);
    d::ComplexPolynomial whole;
    require(polynomial->taylor(whole.value(), center.value(), working), "no expansion about the component");
    for (const slong k : {slong(1), slong(4), slong(8), static_cast<slong>(roots.size())}) {
        d::ComplexBall lower, upper;
        require(polynomial->taylor_pair(lower.value(), upper.value(), center.value(), working, k) &&
                    acb_overlaps(lower.value(), whole.value()->coeffs + k - 1) &&
                    acb_overlaps(upper.value(), whole.value()->coeffs + k),
                "the pair differs from the whole shift's coefficients k - 1 and k");
        // Not vacuously: the pair is about as tight as the shift's own coefficients.
        require(acb_rel_accuracy_bits(upper.value()) > working / 2, "the pair lost the working precision");
    }
    d::Root component;
    acb_one(component.box.value());
    arb_add_error_2exp_si(acb_realref(component.box.value()), -5);
    arb_add_error_2exp_si(acb_imagref(component.box.value()), -5);
    component.count = 8;
    cps_options options;
    options.max_work_bits = 65536;
    certpolysolve::detail::CandidateStatistics statistics{};
    d::EnginePolicy policy;
    d::RootVector hints(8);
    slong used = working;
    const uint64_t shifts = polynomial->taylor_shifts(), pairs = polynomial->taylor_pairs();
    require(d::component_seeds(hints.data(), *polynomial, component, 8, used, options, statistics, policy,
                               std::nullopt) && used == working,
            "the restart supplied no hints at its precision");
    // Three steps, each a pair (more only where a step's precision had to be raised for it); the model is one
    // whole shift, or the expansion's head.
    require(polynomial->taylor_shifts() - shifts <= 1 && polynomial->taylor_pairs() - pairs >= 3,
            "a restart's centre steps took whole Taylor shifts");
}

// A restart's centre stays in its component's box. Four roots 2^-150 about m with a fifth 2^-130 from it, in a
// box 2^-139 wide: a model at 1,024 bits resolves the four. From the model's noise alone (hints 2^-256 out)
// the three steps aim at 65, 130 and 261 bits, and the box sits at the 139th: the first step, resolved and kept
// to fewer bits than the box has, leaves it; from there Newton's iteration sees five roots, not four, and does
// not return; and no hint lands in the box. Every step aims at the box's bits at least.
static void restart_centres_stay_in_their_boxes() {
    std::vector<Complex> roots;
    for (int x = -4; x <= 4; ++x)
        for (int y = -4; y <= 4; ++y)
            if (std::abs(x - 1) > 1 || std::abs(y) > 1)
                roots.emplace_back(x, y);
    // More bits than a step aimed at 65 keeps.
    const Rational real = Rational(1) + power2(-37) + power2(-71) + power2(-113);
    const Rational imaginary = power2(-41) + power2(-83) + power2(-119);
    constexpr int width = 150, neighbour = 130, box = 139;
    for (int re : {-1, 1})
        for (int im : {-1, 1})
            roots.emplace_back(real + re * power2(-width), imaginary + im * power2(-width));
    roots.emplace_back(real + power2(-neighbour), imaginary);
    auto polynomial = source(roots);
    const slong working = 1024;
    d::Root component;
    // Off the roots' mean by a few of the box's bits.
    enclose(component.box.value(), Complex(real + 3 * power2(-box - 4), imaginary + 5 * power2(-box - 5)), working);
    require(acb_is_exact(component.box.value()), "the box's middle is not exact");
    arb_add_error_2exp_si(acb_realref(component.box.value()), -box);
    arb_add_error_2exp_si(acb_imagref(component.box.value()), -box);
    component.count = 4;
    cps_options options;
    options.max_work_bits = 65536;
    certpolysolve::detail::CandidateStatistics statistics{};
    d::EnginePolicy policy;
    d::RootVector hints(4);
    slong used = working;
    require(d::component_seeds(hints.data(), *polynomial, component, 4, used, options, statistics, policy,
                               std::nullopt) && used == working,
            "a restart's hints left its component's box");
    // On the cluster's own circle, 2^-149.5 from m.
    d::ComplexBall mean, offset;
    enclose(mean.value(), Complex(real, imaginary), working);
    for (slong j = 0; j < 4; ++j) {
        d::Magnitude size;
        acb_sub(offset.value(), hints.data() + j, mean.value(), 4 * working);
        acb_get_mag(size.value(), offset.value());
        const double distance = mag_get_d_log2_approx(size.value());
        require(distance > -width - 1 && distance < -width + 2, "a hint is not on its cluster's circle");
    }
}

// An exact local model has no noise to bound its hints, and its centre takes the working precision. Four roots
// 2^-300 about m in a box 2^-40 wide, expanded exactly at 468 bits: about the roots' mean the model's roots are
// the four. A centre located as a model with noise there needs keeps 154 bits, and m has 200: it stands off
// the mean by far more than the cluster is wide, the terms the model drops then move its roots, and the hints
// came 2^-227 from the mean, resolving nothing of the cluster.
static void exact_models_take_the_whole_centre() {
    std::vector<Complex> roots;
    for (int x : {-1, 1, 3})
        for (int y : {-2, 0, 2})
            if (x != 1 || y != 0)
                roots.emplace_back(x, y);
    // More bits than a centre located for 468 bits of noise keeps.
    const Rational real = Rational(1) + power2(-37) + power2(-71) + power2(-113) + power2(-200);
    const Rational imaginary = power2(-41) + power2(-83) + power2(-119) + power2(-201);
    constexpr int width = 300, box = 40;
    for (int re : {-1, 1})
        for (int im : {-1, 1})
            roots.emplace_back(real + re * power2(-width), imaginary + im * power2(-width));
    d::GaussianPolynomial exact;
    for (const auto &coefficient : from_roots(roots))
        exact.push_back(gaussian(coefficient));
    d::PolynomialSourceOptions source_options;
    source_options.input = d::InputSemantics::Exact;
    // The expansion with no rounding before its last step: the roots and every centre here are dyadic, and
    // `whole` bits hold all their products.
    source_options.exact_taylor = [exact](acb_poly_t out, const acb_t center, slong bits) {
        constexpr slong whole = 32768;
        d::ComplexPolynomial original, shifted;
        d::enclose_gaussian(original.value(), exact, whole);
        acb_poly_taylor_shift(shifted.value(), original.value(), center, whole);
        for (slong j = 0; j < acb_poly_length(shifted.value()); ++j)
            require(acb_is_exact(shifted.value()->coeffs + j), "the test's expansion is not exact");
        acb_poly_set_round(out, shifted.value(), bits);
        return true;
    };
    d::PolynomialSource polynomial([exact](acb_poly_t out, slong bits) { d::enclose_gaussian(out, exact, bits); },
                                   static_cast<slong>(roots.size()), false, source_options);
    require(polynomial.exact_local_models(), "the source has no exact local models");
    const slong working = 468;
    d::Root component;
    enclose(component.box.value(), Complex(real + 3 * power2(-box - 4), imaginary + 5 * power2(-box - 5)), working);
    require(acb_is_exact(component.box.value()), "the box's middle is not exact");
    arb_add_error_2exp_si(acb_realref(component.box.value()), -box);
    arb_add_error_2exp_si(acb_imagref(component.box.value()), -box);
    component.count = 4;
    cps_options options;
    options.max_work_bits = 65536;
    certpolysolve::detail::CandidateStatistics statistics{};
    d::EnginePolicy policy;
    d::RootVector hints(4);
    slong used = working;
    require(d::component_seeds(hints.data(), polynomial, component, 4, used, options, statistics, policy,
                               std::nullopt) && used == working,
            "an exact model supplied no hints at its precision");
    // On the cluster's own circle, 2^-299.5 from m.
    d::ComplexBall mean, offset;
    enclose(mean.value(), Complex(real, imaginary), working);
    for (slong j = 0; j < 4; ++j) {
        d::Magnitude size;
        acb_sub(offset.value(), hints.data() + j, mean.value(), 4 * working);
        acb_get_mag(size.value(), offset.value());
        const double distance = mag_get_d_log2_approx(size.value());
        require(distance > -width - 1 && distance < -width + 2, "an exact model's hint does not resolve its cluster");
    }
}

// A restart's local model is the first k + 1 coefficients of the source's
// expansion about the component. Where FLINT's whole Taylor shift is its
// divide-and-conquer recursion and the head's passes are the fewer
// (heuristics::restart_model_from_head), the head alone is computed
// (PolynomialSource::taylor_head: Horner's rule on k + 1 accumulators with
// disk errors). The whole shift for it was 59% of near_multiple-769's
// instructions (192 restarts of four-root clusters at 936 bits;
// measured on 9 October 2026).
static void restart_model_takes_the_expansion_head() {
    // The rule at its two reference inputs: near_multiple-769's restarts take
    // the head; mignotte_generalized_383 at 44,032 bits (six-root clusters),
    // where FLINT's shift is one convolution and a head was the dearer, keeps
    // its whole shift. No head below FLINT's own switch to Horner's rule, nor
    // where the head's passes exceed an eighth of the length.
    require(d::heuristics::restart_model_from_head(770, 936, 5, 9), "near_multiple-769's restarts take no head");
    require(!d::heuristics::restart_model_from_head(384, 44032, 7, 13),
            "a head where FLINT's one convolution is the cheaper");
    require(!d::heuristics::restart_model_from_head(8, 256, 3, 5) &&
                !d::heuristics::restart_model_from_head(81, 128, 9, 17),
            "a head for a short source, or one of more passes than an eighth of the length");
    // Where the shift is one convolution its passes fall with the precision: a four-root model (five passes,
    // nine when its roots are determined) takes the head at 2,752 bits (the convolution: 12 passes or more)
    // and not at 10,944 (7). A source of few terms, whose head is summed term by term in a fraction of a pass
    // (PolynomialSource::taylor_head_passes), takes the head there too, and at mignotte_generalized_383's
    // 44,032 bits; below FLINT's switch to Horner's rule it has none either.
    require(d::heuristics::restart_model_from_head(720, 2752, 5, 9) &&
                !d::heuristics::restart_model_from_head(720, 10944, 5, 9),
            "the head against FLINT's convolution is not judged by the measured passes");
    require(d::heuristics::restart_model_from_head(720, 10944, 0.25, 0.5) &&
                d::heuristics::restart_model_from_head(384, 44032, 0.25, 0.5) &&
                !d::heuristics::restart_model_from_head(8, 256, 0.25, 0.5),
            "the head of a source of few terms is not judged by its own passes");
    // The cluster of local_seeds_do_not_escalate in a far grid of 160 points:
    // 168 roots, so that at 256 bits FLINT's shift is the recursion (256 is
    // at most twice 169) and the model's 17 passes are below an eighth of the
    // length.
    std::vector<Complex> roots;
    for (int x = -6; x <= 6; ++x)
        for (int y = -6; y <= 6; ++y)
            if (std::abs(x - 1) > 1 || std::abs(y) > 1)
                roots.emplace_back(x, y);
    for (int re : {-1, 1})
        for (int im : {-1, 1}) {
            roots.emplace_back(1 + re * power2(-100), im * power2(-100));
            roots.emplace_back(1 + re * power2(-8), im * power2(-8));
        }
    auto polynomial = source(roots);
    const slong working = 256, k = 8;
    require(d::heuristics::restart_model_from_head(static_cast<slong>(roots.size()) + 1, working, k + 1, 2 * k + 1),
            "the fixture is not one whose model comes from the head");
    // The head is the whole shift's first coefficients, at the cluster and at
    // a centre on the diagonal, where ball Horner would lose half a bit a
    // term (84 of the 256 bits here) and the disk errors lose none.
    for (const auto &[re, im] : {std::pair{1 + 0x1p-20, 0x1p-30}, std::pair{0.75, 0.75}}) {
        d::ComplexBall center;
        acb_set_d_d(center.value(), re, im);
        d::ComplexPolynomial whole, head;
        require(polynomial->taylor(whole.value(), center.value(), working) &&
                    polynomial->taylor_head(head.value(), center.value(), working, k + 1) &&
                    acb_poly_length(head.value()) == k + 1,
                "no expansion head about the centre");
        for (slong j = 0; j <= k; ++j)
            require(acb_overlaps(head.value()->coeffs + j, whole.value()->coeffs + j),
                    "the expansion's head differs from the whole shift's first coefficients");
        require(acb_rel_accuracy_bits(head.value()->coeffs + k) > working - 64,
                "the head lost the working precision to its error rectangles");
        d::ComplexPolynomial value;
        require(polynomial->taylor_head(value.value(), center.value(), working, 1) &&
                    acb_poly_length(value.value()) == 1 &&
                    acb_overlaps(value.value()->coeffs, whole.value()->coeffs),
                "a head of length one is not the value at the centre");
    }
    d::Root component;
    acb_one(component.box.value());
    arb_add_error_2exp_si(acb_realref(component.box.value()), -5);
    arb_add_error_2exp_si(acb_imagref(component.box.value()), -5);
    component.count = static_cast<uint64_t>(k);
    cps_options options;
    options.max_work_bits = 65536;
    certpolysolve::detail::CandidateStatistics statistics{};
    d::EnginePolicy policy;
    d::RootVector hints(k);
    slong used = working;
    const uint64_t shifts = polynomial->taylor_shifts(), heads = polynomial->taylor_heads();
    require(d::component_seeds(hints.data(), *polynomial, component, k, used, options, statistics, policy,
                               std::nullopt) && used == working,
            "the restart supplied no hints from the expansion's head");
    require(polynomial->taylor_shifts() == shifts && polynomial->taylor_heads() > heads,
            "the restart's model took a whole Taylor shift");
    // The outer four roots are the model's to 2^-8 of their distance from the
    // centre at the least (the far roots are a unit away); the inner four it
    // does not resolve at 256 bits, and their hints stay within the outer.
    slong outer = 0, inner = 0;
    for (slong j = 0; j < k; ++j) {
        d::ComplexBall offset;
        acb_sub_ui(offset.value(), hints.data() + j, 1, working);
        d::Magnitude size;
        acb_get_mag(size.value(), offset.value());
        const double log2_size = mag_is_zero(size.value()) ? -HUGE_VAL : mag_get_d_log2_approx(size.value());
        if (log2_size > -8 && log2_size < -7)
            ++outer;
        else if (log2_size < -8)
            ++inner;
    }
    require(outer == 4 && inner == 4, "the head's model does not place the cluster's hints at its roots' scales");
}

// The dropped terms decide for a head's model as they do for a whole shift's
// (place_undetermined_hints): the nested pair of
// undetermined_restart_hints_are_placed_by_scale, whose separation the
// truncated model does not determine (the dropped terms at the pair's
// midpoint exceed the model's value there), in a far grid that makes FLINT's
// shift the recursion at 512 bits. From the head there are no dropped
// coefficients to sum: the bound is what the source's value at the midpoint
// differs from the model's by.
static void head_models_measure_their_dropped_terms() {
    std::vector<Complex> roots{
        Complex(Rational("-1/1024")),
        Complex(Rational("708191041404007847/509014844283922939904")),
        Complex(Rational("177047760351001741/127253711070980734976")),
        Complex(Rational("354095520702004365/254507422141961469952")),
        Complex(Rational("177047760351003507/127253711070980734976")),
        Complex(Rational("708191041404011379/509014844283922939904")),
        Complex(Rational("24917257110736165240970126820211/17909367680024424984345422885552128"))};
    // 280 far roots on the grid of eighths, a quarter away at the least.
    for (int x = -8; x <= 8; ++x)
        for (int y = -8; y <= 8; ++y)
            if (std::abs(x) > 1 || std::abs(y) > 1)
                roots.emplace_back(Rational(x) / 8, Rational(y) / 8);
    auto polynomial = source(roots);
    const slong working = 512, k = 6;
    require(d::heuristics::restart_model_from_head(static_cast<slong>(roots.size()) + 1, working, k + 1, 2 * k + 1),
            "the fixture is not one whose model comes from the head");
    d::Root component;
    acb_set_d(component.box.value(), 708191041404007847.0 / 509014844283922939904.0);
    arb_add_error_2exp_si(acb_realref(component.box.value()), -40);
    arb_add_error_2exp_si(acb_imagref(component.box.value()), -40);
    component.count = static_cast<uint64_t>(k);
    cps_options options;
    options.max_work_bits = 65536;
    certpolysolve::detail::CandidateStatistics statistics{};
    d::EnginePolicy policy;
    uint64_t placed = 0;
    policy.trace_data = &placed;
    policy.trace = [](void *p, const d::TraceEvent &event) {
        if (std::string_view(event.stage) == "local_hints_undetermined")
            *static_cast<uint64_t *>(p) += event.active;
    };
    d::RootVector hints(k);
    slong used = working;
    const uint64_t shifts = polynomial->taylor_shifts(), heads = polynomial->taylor_heads();
    require(d::component_seeds(hints.data(), *polynomial, component, k, used, options, statistics, policy,
                               std::nullopt) && used == working,
            "the restart supplied no hints for the nested pair's cluster");
    require(polynomial->taylor_shifts() == shifts && polynomial->taylor_heads() > heads + 1,
            "the model took a whole Taylor shift, or its dropped terms were not measured at a midpoint");
    require(placed >= 2, "a head's model left the nested pair's hints where its dropped terms do not determine them");
}

// A restart's model is computed at the precision that puts its undetermined
// hints within the iteration's resolution of the centre (component_seeds'
// `resolution`). Four roots 2^-300 apart at 1 among 72 far ones: at 256 bits
// the degree-4 model knows its leading coefficient alone, and its hints go on
// a circle as wide as the fourth root of its noise allows, about 2^-50. Nodes
// iterated at 170 bits settle within 2^-168 of their modulus, so from there
// the iteration would have to carry them in, linearly. With the resolution
// asked for, the model is computed again with four bits for every bit of the
// excess, and the hints lie within 2^-170 of the centre.
// near_multiple-193 at the accuracy goal took 5.2 s restarted at 468 bits
// and 1.8 s at 936.
static void restart_models_resolve_to_the_iteration() {
    std::vector<Complex> roots;
    for (int x = -4; x <= 4; ++x)
        for (int y = -4; y <= 4; ++y)
            if (std::abs(x - 1) > 1 || std::abs(y) > 1)
                roots.emplace_back(x, y);
    for (int re : {-1, 1})
        for (int im : {-1, 1})
            roots.emplace_back(1 + re * power2(-300), im * power2(-300));
    auto polynomial = source(roots);
    d::Root component;
    acb_one(component.box.value());
    arb_add_error_2exp_si(acb_realref(component.box.value()), -5);
    arb_add_error_2exp_si(acb_imagref(component.box.value()), -5);
    component.count = 4;
    cps_options options;
    options.max_work_bits = 65536;
    certpolysolve::detail::CandidateStatistics statistics{};
    d::EnginePolicy policy;
    const slong first = 256, resolution = 170;
    // log2 of the largest distance of a hint from the cluster's centre, 1.
    const auto reach = [&](const d::RootVector &hints) {
        double worst = -HUGE_VAL;
        for (slong j = 0; j < 4; ++j) {
            d::ComplexBall offset;
            acb_sub_ui(offset.value(), hints.data() + j, 1, 4096);
            d::Magnitude size;
            acb_get_mag(size.value(), offset.value());
            if (!mag_is_zero(size.value()))
                worst = std::max(worst, mag_get_d_log2_approx(size.value()));
        }
        return worst;
    };
    d::RootVector plain(4), resolved(4);
    slong used = first;
    require(d::component_seeds(plain.data(), *polynomial, component, 4, used, options, statistics, policy,
                               std::nullopt) && used == first,
            "a restart with no resolution asked for raised its precision");
    require(reach(plain) > -resolution, "the control's hints are within the iteration's resolution at 256 bits");
    used = first;
    require(d::component_seeds(resolved.data(), *polynomial, component, 4, used, options, statistics, policy,
                               std::nullopt, resolution),
            "the restart supplied no hints at the iteration's resolution");
    // Four bits for each bit between the first model's reach and the resolution, and what a second model's
    // remaining excess asks for: below four times the resolution beyond the first precision.
    require(used > first && used <= first + 4 * resolution,
            "the model's precision was not raised by what its hints' excess asks for");
    require(reach(resolved) <= -resolution, "the hints lie beyond the iteration's resolution of the centre");
    // Only the model's noise is raised for. A fifth root at the centre: the four others' circle, 2^-300 out,
    // lies between resolved coefficients at 2,048 bits and is the roots' own at every precision; the one
    // circle the unresolved constant bounds is far inside the resolution asked. The precision stays.
    roots.emplace_back(1);
    auto with_centre = source(roots);
    component.count = 5;
    d::RootVector five(5);
    const slong fine = 2048, finer = 340;
    used = fine;
    require(d::component_seeds(five.data(), *with_centre, component, 5, used, options, statistics, policy,
                               std::nullopt, finer) && used == fine,
            "a model was raised for a circle its roots, not its noise, set");
}

// A local model that the working precision does not resolve gives scales, not
// roots. near_multiple-719 scaled by 3: all 179 restarts built their degree-4
// model from coefficients with no correct bit, and in one of them the noise
// midpoints of the two lowest were exactly zero. The midpoint polynomial then
// has a double root at the centre, two restart hints landed on that one point,
// and the iteration never separated them (a lost completion once fix 5 had
// removed the escalation that hid it).
static void local_seeds_from_an_unresolved_model() {
    const slong working = 256;
    // Coefficient j as the ball [midpoint +- 2^radius]; a missing midpoint is zero.
    const auto noise = [](acb_poly_t p, slong j, std::optional<slong> midpoint, slong radius) {
        d::ComplexBall c;
        if (midpoint)
            arf_set_si_2exp_si(arb_midref(acb_realref(c.value())), 1, *midpoint);
        arb_add_error_2exp_si(acb_realref(c.value()), radius);
        arb_add_error_2exp_si(acb_imagref(c.value()), radius);
        acb_poly_set_coeff_acb(p, j, c.value());
    };
    const auto low_noise = [&](acb_poly_t p) {
        noise(p, 0, std::nullopt, -288);
        noise(p, 1, std::nullopt, -292);
        noise(p, 2, -299, -297);
        noise(p, 3, std::nullopt, -303);
    };
    // 2^lower < |z| < 2^upper, decided on the enclosure of |z|.
    const auto size_between = [](const acb_t z, slong lower, slong upper) {
        d::RealBall size, bound;
        acb_abs(size.value(), z, 64);
        arb_one(bound.value());
        arb_mul_2exp_si(bound.value(), bound.value(), lower);
        if (!arb_gt(size.value(), bound.value()))
            return false;
        arb_mul_2exp_si(bound.value(), bound.value(), upper - lower);
        return arb_lt(size.value(), bound.value()) != 0;
    };
    // The bound polygon of u^4 + (noise): one side from (0, -288 + 1/2) to
    // (4, 0), the half for the two components of the radius. Four roots of
    // size 2^-72 at most.
    const auto on_the_resolution_circle = [&](acb_srcptr points, const char *message) {
        for (slong i = 0; i < 4; ++i) {
            require(size_between(points + i, -73, -71), message);
            for (slong j = 0; j < i; ++j)
                require(!acb_overlaps(points + i, points + j), message);
        }
    };
    {
        d::ComplexPolynomial model;
        acb_poly_set_coeff_si(model.value(), 4, 1);
        low_noise(model.value());
        d::RootVector placed(4), midpoints(4);
        bool determined = true;
        require(d::scale_candidates(placed.data(), model.value(), working, determined) && !determined,
                "a model with unresolved polygon vertices was called determined");
        on_the_resolution_circle(placed.data(), "unresolved coefficients did not place the seeds by their bounds");
        // The midpoint placement reads the zero midpoints as zero roots.
        require(d::initial_candidates(midpoints.data(), model.value(), working) && acb_is_zero(midpoints.data()) &&
                    acb_is_zero(midpoints.data() + 1),
                "the fixture no longer has a double midpoint root");
    }
    {
        // Resolved vertices, and one unresolved coefficient under the polygon
        // (as the centroid step leaves the next-to-leading one): the model is
        // determined and the placement is initial_candidates' own.
        d::ComplexPolynomial model;
        acb_poly_set_coeff_si(model.value(), 4, 1);
        d::ComplexBall constant;
        arf_set_si_2exp_si(arb_midref(acb_realref(constant.value())), -3, -80);
        acb_poly_set_coeff_acb(model.value(), 0, constant.value());
        noise(model.value(), 3, -110, -108);
        d::RootVector placed(4), midpoints(4);
        bool determined = false;
        require(d::scale_candidates(placed.data(), model.value(), working, determined) && determined,
                "an unresolved coefficient under the polygon made the model undetermined");
        require(d::initial_candidates(midpoints.data(), model.value(), working), "midpoint placement failed");
        for (slong i = 0; i < 4; ++i)
            require(acb_equal(placed.data() + i, midpoints.data() + i),
                    "a determined model's placement differs from the midpoint placement");
    }
    // The same through the restart. x^80 + x^4 + (noise): four roots near 0
    // below this precision's resolution, 76 on the unit circle. The centre
    // stays 0 (the cubic coefficient's midpoint is zero), so the local model
    // is the noise above.
    d::PolynomialSource polynomial(
        [&](acb_poly_t out, slong) {
            acb_poly_zero(out);
            acb_poly_set_coeff_si(out, 80, 1);
            acb_poly_set_coeff_si(out, 4, 1);
            low_noise(out);
        },
        80, false);
    d::Root component;
    arb_add_error_2exp_si(acb_realref(component.box.value()), -5);
    arb_add_error_2exp_si(acb_imagref(component.box.value()), -5);
    component.count = 4;
    cps_options options;
    options.max_work_bits = 65536;
    certpolysolve::detail::CandidateStatistics statistics{};
    d::EnginePolicy policy;
    d::RootVector hints(4);
    slong used = working;
    require(d::component_seeds(hints.data(), polynomial, component, 4, used, options, statistics, policy,
                               std::nullopt) && used == working,
            "an unresolved local model did not supply restart hints at its own precision");
    on_the_resolution_circle(hints.data(), "restart hints of an unresolved model are not distinct points at its scale");
}

// A restarted member is iterated on the secular model like every other node,
// so it keeps the round's iteration precision. The maintained secular recipe
// gave it the local model's precision, which is the source precision (14,976
// bits for the restarted members of fresh clustered-277 rotated by (3+4i)/5,
// against 170 for the rest: 54 of that solve's 65 seconds).
// The geometric restart (EnginePolicy::geometric_local_restart, P5 lead 1):
// four clusters of 48 roots about +-1 +- i, each a row of roots 2^-23 apart
// with a jitter (clustered-769's shape at a quarter of its degree). Plain
// secular closes in on each cluster linearly: 26 validations and 318 sweeps.
// A group of nodes that crawls, whose disk's root count is certified exactly,
// is re-seeded from its local model at the roots' mean, and every cluster is
// certified within five validations. Every planted root lies in exactly one
// region, and each region holds exactly its count of them.
// The geometric restart's nested solve (local_model_options) places a model's roots without exact preparation: the
// model is a rounded Taylor polynomial of the input, so its square-free proofs and factors say nothing about the input
// (D1's allowance of 26 October had let the nested solve prepare such models exactly). Its options carry a zero exact
// budget, which refuses every exact step: a solve under them meets its goal with the square-free verdict unknown and
// no exact factor, on real and on complex input, where the same solve with the default budget proves it square-free.
static void local_model_solves_without_exact_preparation() {
    cps_options caller;
    caller.goal = CPS_ROOT_ACCURACY;
    caller.bits = 200;
    caller.max_work_bits = 4096;
    const auto allowed = d::local_model_options(caller, 80, {}, {});
    require(allowed.has_value(), "an ample local-model allowance was refused");
    const auto &nested = *allowed;
    require(nested.exact_budget == 0 && nested.exact_storage_limit() == 0, "the nested solve keeps an exact budget");
    require(nested.goal == CPS_ROOT_APPROXIMATION && nested.method == CPS_METHOD_SECULAR && nested.bits == 80 &&
                nested.max_work_bits == caller.max_work_bits,
            "the nested solve's goal, method or precision cap changed");
    cps_options control = nested;
    control.exact_budget = cps_options{}.exact_budget;
    for (const bool complex : {false, true}) {
        std::vector<Complex> roots;
        for (int j = 0; j < 6; ++j)
            roots.emplace_back(Rational(2 * j + 1) * power2(-3), complex ? Rational(j % 3 + 1) * power2(-2) : Rational(0));
        const auto polynomial = input(from_roots(roots));
        for (const bool zero : {true, false}) {
            const auto result = [&] {
                if (zero) return d::solve_local_model(d::dense(*polynomial.get()), nested);
                d::SolveState state(d::dense(*polynomial.get()), d::maintained_method_policy(CPS_METHOD_SECULAR));
                return state.advance(&control);
            }();
            require(result.complete && result.goal_met, "a model was not solved at the nested goal");
            if (zero)
                require(result.squarefree == CPS_SQUAREFREE_UNKNOWN && result.statistics.exact_factors == 0,
                        "the nested solve prepared its model exactly");
            else
                require(result.squarefree == CPS_SQUAREFREE_PROVED,
                        "the control did not prove the model square-free (the test would not tell)");
        }
    }
}

static void local_model_budget_options() {
    cps_options caller;
    caller.work_units = 7;
    caller.max_iterations = 11;
    d::SearchStatistics search{};
    d::CandidateStatistics candidates{};
    search.work_units = 6;
    candidates.candidate_sweeps = 10;
    auto allowed = d::local_model_options(caller, 80, search, candidates);
    require(allowed && allowed->work_units == 1 && !allowed->work_units_default && allowed->max_iterations == 1,
            "local model did not receive the remaining work and sweep allowances");
    search.work_units = 7;
    require(!d::local_model_options(caller, 80, search, candidates), "spent work became a fresh local allowance");
    search.work_units = 8;
    require(!d::local_model_options(caller, 80, search, candidates), "overspent work wrapped into an allowance");
    search.work_units = 0;
    candidates.candidate_sweeps = 11;
    require(!d::local_model_options(caller, 80, search, candidates), "spent sweeps became an unlimited allowance");
    caller.max_iterations = 0;
    caller.work_units = UINT64_MAX;
    allowed = d::local_model_options(caller, 80, search, candidates);
    require(allowed && !allowed->max_iterations && allowed->work_units == cps_options{}.work_units &&
                !allowed->work_units_default,
            "unlimited outer sweeps or the local work ceiling changed");
}

static thread_local const d::CandidateStatistics *local_parent_statistics = nullptr;
static constexpr uint64_t local_mock_work = 3, local_mock_sweeps = 2;
static thread_local unsigned local_capability_calls = 0;
static thread_local bool misleading_local_capability = false;
static cps_result controlled_local_solver(d::DensePolynomial &&input, const cps_options &options) {
    require(options.work_units == cps_options{}.work_units - local_mock_work * local_capability_calls &&
                !options.work_units_default && local_parent_statistics &&
                options.max_iterations == cps_options{}.work_units - local_parent_statistics->candidate_sweeps,
            "a later local solve received work or sweeps already spent");
    ++local_capability_calls;
    require(options.exact_storage_limit() == 0 && options.method == CPS_METHOD_SECULAR &&
                options.goal == CPS_ROOT_APPROXIMATION && !options.real_only && !options.region.set(),
            "local capability received an unresolved or unbounded query");
    cps_result result;
    result.work_bits = 128;
    result.progress.work_units_this_call = local_mock_work;
    result.statistics_v2.work_units = local_mock_work;
    result.statistics_v3.candidate_sweeps = local_mock_sweeps;
    result.statistics_v3.evaluation_terms = 7;
    if (misleading_local_capability) {
        result.complete = result.goal_met = true;
        for (size_t i = 1; i < input.coefficients.size(); ++i) {
            d::Root hint;
            acb_set_si(hint.box.value(), static_cast<slong>(100 + i));
            hint.count = 1;
            result.roots.push_back(std::move(hint));
        }
    }
    return result;
}

static void local_capability_cannot_publish_source_roots() {
    const std::vector<std::pair<int, int>> centres{{-4, 1}, {0, 2}, {4, 1}};
    std::vector<Complex> roots;
    // Forty-eight members converge slowly enough to trigger geometric restart.
    // Four-member groups finish ordinary refinement before that crawl test fires.
    for (const auto &[x, y] : centres)
        for (int re = -7; re <= 7; re += 2) for (int im = -5; im <= 5; im += 2)
            roots.emplace_back(Rational(x) + re * power2(-30), Rational(y) + im * power2(-30));
    for (bool misleading : {false, true}) {
        local_capability_calls = 0;
        misleading_local_capability = misleading;
        d::RootVector seeds(static_cast<slong>(roots.size()));
        slong index = 0;
        for (const auto &[x, y] : centres)
            for (int re = -7; re <= 7; re += 2) for (int im = -5; im <= 5; im += 2)
                acb_set_d_d(seeds[index++], x + re * 0x1p-8, y + im * 0x1p-8);
        auto policy = d::maintained_method_policy(CPS_METHOD_SECULAR);
        policy.newton_point_refinement = policy.secular_preliminary = false;
        d::SecularRefinement state(source(roots), seeds.data(), 1, 85, policy, controlled_local_solver);
        cps_options options;
        options.goal = CPS_ROOT_ACCURACY;
        options.max_work_bits = 4096;
        d::SearchStatistics search{};
        d::CandidateStatistics statistics{};
        local_parent_statistics = &statistics;
        options.max_iterations = cps_options{}.work_units; // ample, but explicitly shared with every local call
        unsigned after_local = 0;
        for (unsigned call = 0; call < 64 && after_local < 4 && state.active(); ++call) {
            state.step(options, search, statistics);
            if (local_capability_calls) ++after_local;
        }
        if (!local_capability_calls)
            std::cerr << "local fixture: active=" << state.active() << ", sweeps=" << statistics.candidate_sweeps
                      << ", validations=" << statistics.validation_calls << ", bits=" << state.precision() << '\n';
        require(local_capability_calls != 0, "controlled local capability fixture never reached a local solve");
        require(search.work_units == local_mock_work * local_capability_calls &&
                    statistics.candidate_sweeps >= local_mock_sweeps * local_capability_calls,
                "local success or refusal hid its work or sweeps");
        local_parent_statistics = nullptr;
        std::vector<d::Root> cover;
        if (state.cover(cover)) {
            size_t total = 0;
            for (const auto &region : cover) {
                size_t inside = 0;
                for (const auto &root : roots) {
                    d::ComplexBall exact;
                    enclose(exact.value(), root, 256);
                    inside += acb_contains(region.box.value(), exact.value()) != 0;
                }
                require(inside == region.count, "a local hint was published as an original-source root");
                total += inside;
            }
            require(total == roots.size() && d::disjoint(cover), "local refusal lost original-source coverage");
        }
    }
}

static void geometric_restart_ends_a_crawl() {
    constexpr int degree = 192, denominator_bits = 28;
    const Rational unit = power2(-denominator_bits);
    std::vector<Complex> roots;
    for (int j = 0; j < degree; ++j) {
        const int jitter = (7 * j + 3) % 16;
        const Rational re = Rational(j % 2 ? 1 : -1) + Rational(32 * (j + 1) + jitter) * unit;
        const Rational im = Rational(j % 4 < 2 ? 1 : -1) + Rational((5 * j + 1) % 16) * unit;
        roots.emplace_back(re, im);
    }
    const auto polynomial = input(from_roots(roots));
    cps_options options;
    options.goal = CPS_ROOT_APPROXIMATION;
    options.bits = 53;
    options.max_work_bits = 65536;
    // The ample solve reaches eight local models. These bounds leave only one
    // work unit or one sweep at the first model; neither is a fresh allowance
    // for each child. Observe charges at the production local-call boundary.
    struct Charges {
        const d::SearchStatistics *search = nullptr;
        uint64_t units_before = 0, sweeps_before = 0, units = 0, sweeps = 0, calls = 0;
        static void record(void *context, const d::TraceEvent &event) noexcept {
            auto &self = *static_cast<Charges *>(context);
            if (event.search) self.search = event.search;
            if (std::strcmp(event.stage, "local_candidates") || !self.search || !event.candidates) return;
            if (!std::strcmp(event.phase, "begin")) {
                self.units_before = self.search->work_units;
                self.sweeps_before = event.candidates->candidate_sweeps;
                ++self.calls;
            } else if (!std::strcmp(event.phase, "end")) {
                self.units += self.search->work_units - self.units_before;
                self.sweeps += event.candidates->candidate_sweeps - self.sweeps_before;
            }
        }
    };
    for (const auto &[work, sweeps] : std::vector<std::pair<uint64_t, uint64_t>>{{100000, 0}, {8, 0}, {100000, 45}}) {
        Charges charges;
        auto policy = d::maintained_method_policy(CPS_METHOD_SECULAR);
        policy.trace = Charges::record;
        policy.trace_data = &charges;
        options.work_units = work;
        options.work_units_default = false;
        options.max_iterations = sweeps;
        d::SolveState state(certpolysolve::detail::dense(*polynomial.get()), policy);
        const auto result = state.advance(&options);
        require(charges.calls && charges.units, "the production local solve was not reached or hid its work");
        require(result.progress.work_units_this_call == result.statistics_v2.work_units &&
                    result.progress.work_units_this_call <= work &&
                    (!sweeps || result.statistics_v3.candidate_sweeps <= sweeps),
                "nested work was omitted from progress or exceeded the caller's remaining budget");
        require(result.complete, "a bounded local restart lost its certified cover");
        if (work == 100000 && !sweeps) {
            require(result.goal_met && charges.sweeps, "the four clusters were not solved or hid their local sweeps");
            require(result.statistics_v3.validation_calls <= 8, "the restart did not end the crawl onto the clusters");
        } else if (!result.goal_met) {
            require(result.progress.stop_reason == (sweeps ? CPS_STOP_ITERATION_BUDGET : CPS_STOP_WORK_BUDGET),
                    "a bounded local restart lost the reason for stopping");
        }
        const auto exported = wrap(result);
        const auto &regions = exported.roots();
        check_disjoint(regions);
        std::vector<uint64_t> held(regions.size(), 0);
        for (const auto &root : roots) {
            size_t inside = 0, last = 0;
            for (size_t k = 0; k < regions.size(); ++k) {
                const auto &r = regions[k];
                if (endpoint(r.box.real.lo) <= root.real && root.real <= endpoint(r.box.real.hi) &&
                    endpoint(r.box.imag.lo) <= root.imag && root.imag <= endpoint(r.box.imag.hi)) {
                    ++inside;
                    last = k;
                }
            }
            require(inside == 1, "a planted root is not in exactly one region");
            ++held[last];
        }
        for (size_t k = 0; k < regions.size(); ++k)
            require(held[k] == regions[k].count, "a region does not hold exactly its count of planted roots");
    }
}

// Fresh clustered-719's recipe at degree 192, moved off the origin by 3/8:
// four clusters of 48 real roots 2^-94 apart about -19/8, -11/8, -3/8 and 5/8.
// The restart's hints are closer together than binary64 resolves. Iterated
// at the round's 85 bits, each cluster's hints rounded to one point, the next
// packet found the nodes at the arithmetic's floor, and the precision growth
// that followed placed them again (six restarts, 288 nodes placed, 7 s). Now
// the members get the bits their hints need, and each node is placed once.
static void restart_hints_stay_apart() {
    constexpr int degree = 192, separation = 96;
    const Rational step = power2(-(separation + 4));
    std::vector<Complex> roots;
    for (int j = 0; j < degree; ++j)
        roots.emplace_back(Rational(j % 4 - 2) - Rational(3) * power2(-3) +
                               Rational(16 * (j + 1) + (7 * j + 3) % 16) * step,
                           Rational(0));
    const auto polynomial = input(from_roots(roots));
    auto policy = d::maintained_method_policy(CPS_METHOD_SECULAR);
    require(policy.geometric_local_restart, "the maintained secular recipe has no geometric restart");
    policy.study_statistics = std::make_shared<d::StudyStatistics>();
    d::SolveState state(certpolysolve::detail::dense(*polynomial.get()), policy);
    cps_options options;
    options.goal = CPS_ROOT_APPROXIMATION;
    options.bits = 53;
    options.max_work_bits = 65536;
    const auto result = state.advance(&options);
    require(result.complete && result.goal_met, "the four sub-binary64 clusters were not solved");
    require(policy.study_statistics->candidate_restarts >= 4, "the clusters were not restarted from their models");
    require(policy.study_statistics->candidates_replaced <= static_cast<uint64_t>(degree),
            "a restarted node was placed again: the restart's hints did not survive the next packet");
    require(result.statistics_v3.validation_calls <= 10, "the restarted clusters took too many rounds");
}

// Fresh clustered-719's recipe moved off the origin by 3/8: `clusters` rows of real roots 2^-(separation - 2)
// apart, with a jitter, about the integers from -1 up less 3/8 (root j in cluster j mod clusters).
static std::vector<Complex> tight_clusters(int degree, int clusters, int separation) {
    const Rational step = power2(-(separation + 4));
    std::vector<Complex> roots;
    for (int j = 0; j < degree; ++j)
        roots.emplace_back(Rational(j % clusters - 1) - Rational(3) * power2(-3) +
                               Rational(16 * (j + 1) + (7 * j + 3) % 16) * step,
                           Rational(0));
    return roots;
}

// The restart's centre. Three clusters of 72 roots 2^-94 apart (2^-88 wide) under the accuracy goal, which
// isolates every root. The moments give a cluster's mean to one limb, here about 2^-55 of its modulus; a model
// about a point that far from a cluster 2^-88 wide shows the cluster as a point, its hints were a ring of
// rounding noise outside the cluster, and the members closed on it linearly, 2 / 72 of the distance a sweep:
// 594 sweeps and 30 validations here (18.7 s), and the three fresh inputs of degree 719, with clusters of 180,
// never finished at this goal. With the mean located by
// Newton's steps on p^(K-1) the model shows the cluster's shape and its roots are the hints: 76 sweeps and 6
// validations (0.45 s).
// Four times the sweeps the solve takes.
constexpr uint64_t restart_centre_sweeps = 4 * 76;
static void restart_centre_is_the_roots_mean() {
    constexpr int degree = 216, clusters = 3, separation = 96;
    const auto roots = tight_clusters(degree, clusters, separation);
    const auto polynomial = input(from_roots(roots));
    auto policy = d::maintained_method_policy(CPS_METHOD_SECULAR);
    require(policy.geometric_local_restart, "the maintained secular recipe has no geometric restart");
    policy.study_statistics = std::make_shared<d::StudyStatistics>();
    d::SolveState state(certpolysolve::detail::dense(*polynomial.get()), policy);
    cps_options options;
    options.goal = CPS_ROOT_ACCURACY;
    options.bits = 53;
    options.max_work_bits = 65536;
    options.max_iterations = restart_centre_sweeps;
    const auto result = state.advance(&options);
    require(result.complete && result.goal_met,
            "the clusters' roots were not isolated within the sweeps of a restart from their means");
    require(result.roots.size() == static_cast<size_t>(degree), "a cluster was not resolved into its roots");
    require(policy.study_statistics->candidate_restarts >= 1, "no cluster was restarted from its model");
}

// The same recipe with three clusters of 90 roots under the approximation goal, which takes a cluster within
// its tolerance as it stands. A cluster here is 2^-88 of its distance from the others in radius, so a model
// about its mean shows its shape only with 90 x 88 bits and what the evaluation cancels: beyond 7,488, the
// ladder's precision at which the count already dominates on the goal's circle (2^-70 of the modulus), and at
// its next. Taken there, the model made this solve 1.6 times as long and four such clusters 3.6 times (0.69 ->
// 1.11 s and 3.2 -> 11.5 s: every later evaluation at 14,976 bits) for a shape the goal does not ask for. The
// shift stops where the cluster lies within the goal's circle, and the hints go on that circle.
static void restart_within_the_goal_takes_no_finer_model() {
    constexpr int degree = 270, clusters = 3, separation = 96;
    // The bits of a cluster's radius: its roots lie within degree 2^-separation of its first.
    constexpr uint64_t radius_bits = separation - 8;
    static_assert(degree < (1 << 9), "the cluster's extent is below 2^(9 - separation), its radius 2^(8 - separation)");
    constexpr uint64_t shape_bits = (degree / clusters) * radius_bits;
    const auto polynomial = input(from_roots(tight_clusters(degree, clusters, separation)));
    auto policy = d::maintained_method_policy(CPS_METHOD_SECULAR);
    policy.study_statistics = std::make_shared<d::StudyStatistics>();
    d::SolveState state(certpolysolve::detail::dense(*polynomial.get()), policy);
    cps_options options;
    options.goal = CPS_ROOT_APPROXIMATION;
    options.bits = 53;
    options.max_work_bits = 65536;
    const auto result = state.advance(&options);
    require(result.complete && result.goal_met, "the clusters were not solved at the approximation goal");
    require(policy.study_statistics->candidate_restarts >= 1, "no cluster was restarted from its model");
    require(result.work_bits < shape_bits, "a goal that takes a cluster as it stands paid for the cluster's shape");
}

// The active nodes as a whole are a group once everything else is certified:
// 48 roots 2^-40 apart about 1 + i/2 among 60 roots on the integer grid. The
// grid roots are certified first; the cluster's nodes, then the whole active
// set, were never a group and crawled for 46 validations.
static void last_active_cluster_is_a_group() {
    const Rational unit = power2(-40);
    std::vector<Complex> roots;
    for (int j = 0; j < 48; ++j)
        roots.emplace_back(Rational(1) + Rational(3 * j + 1) * unit,
                           power2(-1) + Rational((5 * j + 2) % 7) * unit);
    size_t grid = 0;
    for (int x = -4; x <= 4; ++x)
        for (int y = -3; y <= 3; ++y)
            if (grid < 60 && !(x == 1 && (y == 0 || y == 1))) {
                roots.emplace_back(Rational(x), Rational(y));
                ++grid;
            }
    const auto polynomial = input(from_roots(roots));
    auto policy = d::maintained_method_policy(CPS_METHOD_SECULAR);
    policy.study_statistics = std::make_shared<d::StudyStatistics>();
    d::SolveState state(certpolysolve::detail::dense(*polynomial.get()), policy);
    cps_options options;
    options.goal = CPS_ROOT_APPROXIMATION;
    options.bits = 53;
    options.max_work_bits = 65536;
    const auto result = state.advance(&options);
    require(result.complete && result.goal_met, "the grid and its cluster were not solved");
    require(policy.study_statistics->candidate_restarts >= 1, "the last active cluster was not restarted");
    require(result.statistics_v3.validation_calls <= 12, "the last active cluster crawled");
}

// A node on the wrong cluster does not cost the others their places. Twenty
// clusters of four roots 2^-69 wide on a grid, which binary64 cannot resolve,
// seeded with five points at one and three at another: the five nodes serve
// four roots, their inclusion disks (n times a residue) reach every other
// cluster, and the cover is one component of all 80 nodes, more than a local
// model takes. The precision growth restarted such a component whole, on one
// circle of its own size, the 72 nodes of the eighteen sound clusters with
// it (fresh near_multiple-719, where the binary64 stage ends so when its
// packet sums in lanes: 11 s -> 273 s; near_multiple-193 rotated by i, where
// it ends so with the serial sums: 0.35 s -> 6.8 s). The model's exact
// counts show each cluster's nodes about their own roots: the sound ones are
// restarted from their local models, as certified components of four are,
// and the surplus node moves to the cluster that lacks one
// (SecularRefinement::regroup_restart).
static void a_misplaced_node_keeps_the_other_clusters() {
    constexpr int columns = 5, rows = 4, width = 70;
    std::vector<Complex> roots;
    std::vector<std::pair<int, int>> centres;
    for (int x = 0; x < columns; ++x)
        for (int y = 0; y < rows; ++y) {
            centres.emplace_back(2 * x - 4, 2 * y - 3);
            for (int re : {-1, 1})
                for (int im : {-1, 1})
                    roots.emplace_back(Rational(2 * x - 4) + re * power2(-width), Rational(2 * y - 3) + im * power2(-width));
        }
    const size_t n = roots.size();
    struct Seen {
        uint64_t regrouped = 0, restarted = 0, reseeded = 0;
    };
    // Four seeds 2^-6 about each centre; misplaced, the first cluster takes the second's fourth as well.
    const auto solve = [&](bool misplaced, Seen &seen) {
        d::RootVector seeds(static_cast<slong>(n));
        slong next = 0;
        for (size_t cluster = 0; cluster < centres.size(); ++cluster)
            for (int corner = 0; corner < 4; ++corner) {
                const bool moved = misplaced && cluster == 1 && corner == 3;
                const auto &[x, y] = centres[moved ? 0 : cluster];
                const double reach = moved ? 0x1.8p-6 : 0x1p-6;
                acb_set_d_d(seeds[next++], x + (moved ? reach : (corner & 1 ? reach : -reach)),
                            y + (moved ? 0.0 : (corner & 2 ? reach : -reach)));
            }
        auto policy = d::maintained_method_policy(CPS_METHOD_SECULAR);
        // The certificate attempt and the preliminary phase would place the points first.
        policy.newton_point_refinement = policy.secular_preliminary = false;
        policy.trace_data = &seen;
        policy.trace = [](void *p, const d::TraceEvent &event) {
            auto &s = *static_cast<Seen *>(p);
            if (std::string_view(event.stage) == "secular_component_regroup") {
                ++s.regrouped;
                s.restarted += event.certified;
            } else if (std::string_view(event.stage) == "secular_component_reseed")
                ++s.reseeded;
        };
        d::SecularRefinement state(source(roots), seeds.data(), 1, 85, policy, d::solve_local_model);
        cps_options options;
        options.goal = CPS_ROOT_APPROXIMATION;
        options.bits = 53;
        options.max_work_bits = 4096;
        certpolysolve::detail::SearchStatistics search{};
        certpolysolve::detail::CandidateStatistics statistics{};
        for (unsigned call = 0; call < 4096 && state.active(); ++call)
            state.step(options, search, statistics);
        require(!state.active(), "the clusters were not certified");
        return statistics.evaluations;
    };
    Seen sound, misplaced;
    const uint64_t sound_evaluations = solve(false, sound);
    const uint64_t misplaced_evaluations = solve(true, misplaced);
    require(sound.regrouped == 0 && sound.reseeded == 0, "the control's cover is one component without a misplaced node");
    require(misplaced.regrouped >= 1 && misplaced.reseeded == 0,
            "a component with a misplaced node was re-seeded whole, or its groups were not found");
    require(misplaced.restarted >= centres.size() - 2, "the sound clusters were not restarted from their local models");
    // The misplaced node costs the binary64 rounds that end with it in place, not a second solve: re-seeded
    // whole the run took 22,422 evaluations.
    require(misplaced_evaluations < 2 * sound_evaluations, "a misplaced node cost more than the solve itself");
}

// A restart's precisions follow from what its hints are for. Twenty clusters
// of four roots 2^-300 wide at the accuracy goal: every round's iteration
// doubles its bits r, and the clusters are restarted until r resolves them.
// A model of four roots places its hints 2^(-N / 4) out for noise 2^-N, so
// 4 r bits and the bits the cluster's evaluation cancels put them within the
// iteration's resolution (local_component_restart); and the members' next
// residues need that precision and the accuracy an adequate residue has
// (their hint to the regeneration). The model used to take the
// regeneration's precision, and the regeneration to climb from the
// iteration's bits by doubling: near_multiple-719's last round evaluated at
// 22,016 bits where about 11,000 do.
static void restart_precisions_follow_the_hints() {
    constexpr int columns = 4, rows = 4, width = 300;
    constexpr uint64_t members = 5;
    std::vector<Complex> roots;
    for (int x = 0; x < columns; ++x)
        for (int y = 0; y < rows; ++y) {
            roots.emplace_back(Rational(2 * x - 3), Rational(2 * y - 3));
            for (int re : {-1, 1})
                for (int im : {-1, 1})
                    roots.emplace_back(Rational(2 * x - 3) + re * power2(-width), Rational(2 * y - 3) + im * power2(-width));
        }
    const auto polynomial = input(from_roots(roots));
    auto policy = d::maintained_method_policy(CPS_METHOD_SECULAR);
    struct Round {
        uint64_t resolution = 0, model = 0, residues = 0;
    };
    struct Seen {
        std::vector<Round> rounds;
        bool restarting = false;
    } seen;
    policy.trace_data = &seen;
    policy.trace = [](void *p, const d::TraceEvent &event) {
        auto &s = *static_cast<Seen *>(p);
        const std::string_view stage(event.stage), phase(event.phase);
        if (stage == "secular_candidate_restart" && phase == "begin" && !s.restarting) {
            s.restarting = true;
            s.rounds.push_back({event.precision, 0, 0});
        } else if (stage == "local_candidates" && phase == "begin" && s.restarting)
            s.rounds.back().model = std::max(s.rounds.back().model, event.precision);
        else if (stage == "regeneration" && phase == "end" && s.restarting) {
            s.rounds.back().residues = event.precision;
            s.restarting = false;
        }
    };
    d::SolveState state(certpolysolve::detail::dense(*polynomial.get()), policy);
    cps_options options;
    options.goal = CPS_ROOT_ACCURACY;
    options.bits = 53;
    options.max_work_bits = 65536;
    const auto result = state.advance(&options);
    require(result.goal_met && result.complete, "the clusters were not isolated");
    require(seen.rounds.size() >= 2, "the control no longer restarts its clusters round after round");
    // These clusters cancel few bits (a limb covers them): a model at five bits a bit of the resolution, and
    // residues a limb or two beyond the model (the adequate accuracy, 69 bits here, the guard, and what a
    // measured shortfall adds): 896 and 997 bits in the first round, 1,792 and 1,893 in the second, where
    // the models were at 863 and 2,496 and the residues, by doubling from the iteration's bits, at 1,616
    // and 2,976.
    // A later round's model starts at that sizing with the cancelled bits the first round's raise showed (46
    // here, of which the rounding to whole limbs is most), rounded to whole limbs again: 1,792 bits in the
    // second round, one model, where the model started at the regeneration's precision and was raised to 1,728.
    constexpr uint64_t limb = FLINT_BITS;
    for (const auto &round : seen.rounds) {
        require(round.model >= members * round.resolution && round.model <= members * round.resolution + 2 * limb,
                "a restart's model is not at the precision its hints need");
        require(round.residues >= round.model && round.residues <= round.model + 3 * limb,
                "the restarted members' residues are not at the precision their model showed");
    }
}

// Option 2 (6 October 2026): a node near the original origin of a centred
// source iterates at the bits of the goal at its distance from that origin,
// the origin the goals are relative to. The restart recipe at degree 256
// (four jittered clusters of 64, separation 48) rotated by (3 + 4i) / 5 is
// centred near (3 + 4i) / 5 times -1/2: the cluster at the origin, of moduli
// 2^-48 to 2^-40, needs about 90 bits of the centred coordinates. With the
// rule 11 validations; without it the cluster's nodes converge only as far as
// binary64 holds them and wait for the precision growth (14). On the
// unrotated recipe the rule cost instead (it sent every node of a round to
// the multiprecision packet), which mixed packets remove
// (mixed_rounds_sweep_in_binary64).
static void origin_nodes_get_coordinate_bits() {
    constexpr int degree = 256, separation = 48;
    const Rational step = power2(-(separation + 4));
    const Complex turn(Rational(3, 5), Rational(4, 5));
    std::vector<Complex> roots;
    for (int j = 0; j < degree; ++j)
        roots.push_back(Complex(Rational(j % 4 - 2) + Rational(16 * (j + 1) + (7 * j + 3) % 16) * step) * turn);
    const auto polynomial = input(from_roots(roots));
    auto policy = d::maintained_method_policy(CPS_METHOD_SECULAR);
    require(policy.centred_origin_basis, "the maintained secular recipe has no origin basis");
    d::SolveState state(certpolysolve::detail::dense(*polynomial.get()), policy);
    cps_options options;
    options.goal = CPS_ROOT_APPROXIMATION;
    options.bits = 53;
    options.max_work_bits = 65536;
    const auto result = state.advance(&options);
    require(result.complete && result.goal_met, "the rotated clusters were not solved");
    require(result.statistics_v3.validation_calls <= 12,
            "the cluster at the original origin waited for the precision growth");
}

// Mixed packets (EnginePolicy::mixed_packets): a round in which some active
// nodes iterate above the initial stage still sweeps the others in binary64.
// The restart recipe of origin_nodes_get_coordinate_bits, unrotated, has a
// cluster near the original origin that the coordinate rule raises above the
// stage; without mixed packets every node of those rounds went to the
// multiprecision engine.
static void mixed_rounds_sweep_in_binary64() {
    constexpr int degree = 256, separation = 48;
    const Rational step = power2(-(separation + 4));
    std::vector<Complex> roots;
    for (int j = 0; j < degree; ++j)
        roots.push_back(Complex(Rational(j % 4 - 2) + Rational(16 * (j + 1) + (7 * j + 3) % 16) * step));
    const auto polynomial = input(from_roots(roots));
    cps_options options;
    options.goal = CPS_ROOT_APPROXIMATION;
    options.bits = 53;
    options.max_work_bits = 65536;
    const slong stage = d::SecularMachineWorkspace::stage_bits(options.bits);
    for (bool mixed : {true, false}) {
        auto policy = d::maintained_method_policy(CPS_METHOD_SECULAR);
        require(policy.mixed_packets, "the maintained secular recipe has no mixed packets");
        policy.mixed_packets = mixed;
        policy.study_statistics = std::make_shared<d::StudyStatistics>();
        const auto counts = policy.study_statistics;
        uint64_t machine_before = 0;
        size_t high = 0, high_binary64 = 0;
        policy.snapshot = [&](d::NumericalSnapshot event) {
            if (event.stage == "after_secular_packet") {
                slong bits = 0;
                for (size_t i = 0; i < event.coordinate_precision.size(); ++i)
                    if (i >= event.inactive.size() || !event.inactive[i])
                        bits = std::max(bits, event.coordinate_precision[i]);
                if (bits > stage) {
                    ++high;
                    high_binary64 += counts->machine_packets > machine_before;
                }
            }
            machine_before = counts->machine_packets;
        };
        d::SolveState state(certpolysolve::detail::dense(*polynomial.get()), policy);
        const auto result = state.advance(&options);
        require(result.complete && result.goal_met, "the unrotated clusters were not solved");
        require(high > 0, "no round had a node above the initial stage: the check below proves nothing");
        if (mixed)
            require(high_binary64 > 0, "a round with a node above the initial stage sent every node to the "
                                       "multiprecision engine");
        else
            require(high_binary64 == 0, "without mixed packets such a round ran binary64 sweeps");
    }
}

// A mixed round under the solve's iteration budget (cps_options::max_iterations): each of its sweeps is one
// binary64 sweep of the nodes within the stage and one multiprecision sweep of the others, both counted in
// candidate_sweeps, so a budget of k sweeps allows k counted sweeps, not k pairs (finding S04 of the review
// of 14 October 2026). The fixture of mixed_rounds_sweep_in_binary64: doubling budgets bracket the first mixed
// round, then bisect to an adjacent unmixed/mixed pair and cut at every sweep of the following mixed packet.
// Bracketing avoids hundreds of redundant solves when nested work moves the boundary. Each fresh solve must
// keep its count; a budget of one never starts a mixed round, and
// the unbudgeted solve finishes. (Until 18 October 2026 one state was advanced call by call under each budget;
// a solve state is advanced once since.)
static void mixed_rounds_respect_the_iteration_budget() {
    constexpr int degree = 256, separation = 48;
    const Rational step = power2(-(separation + 4));
    std::vector<Complex> roots;
    for (int j = 0; j < degree; ++j)
        roots.push_back(Complex(Rational(j % 4 - 2) + Rational(16 * (j + 1) + (7 * j + 3) % 16) * step));
    const auto polynomial = input(from_roots(roots));
    cps_options options;
    options.goal = CPS_ROOT_APPROXIMATION;
    options.bits = 53;
    options.max_work_bits = 65536;
    const slong stage = d::SecularMachineWorkspace::stage_bits(options.bits);
    const auto run = [&](uint64_t budget, size_t &high_binary64) {
        options.max_iterations = budget;
        auto policy = d::maintained_method_policy(CPS_METHOD_SECULAR);
        policy.study_statistics = std::make_shared<d::StudyStatistics>();
        const auto counts = policy.study_statistics;
        uint64_t machine_before = 0;
        high_binary64 = 0;
        policy.snapshot = [&](d::NumericalSnapshot event) {
            if (event.stage == "after_secular_packet") {
                slong bits = 0;
                for (size_t i = 0; i < event.coordinate_precision.size(); ++i)
                    if (i >= event.inactive.size() || !event.inactive[i])
                        bits = std::max(bits, event.coordinate_precision[i]);
                if (bits > stage) high_binary64 += counts->machine_packets > machine_before;
            }
            machine_before = counts->machine_packets;
        };
        d::SolveState state(certpolysolve::detail::dense(*polynomial.get()), policy);
        const auto result = state.advance(&options);
        if (budget)
            require(result.statistics_v3.candidate_sweeps <= budget,
                    ("a solve made " + std::to_string(result.statistics_v3.candidate_sweeps) +
                     " sweeps under a budget of " + std::to_string(budget)).c_str());
        return result;
    };
    size_t high = 0;
    uint64_t first_mixed = 0, last_unmixed = 1;
    for (uint64_t budget = 2; budget <= 4096 && !first_mixed; budget *= 2) {
        const auto result = run(budget, high);
        if (high > 0) first_mixed = budget;
        else if (result.goal_met) break;
        else last_unmixed = budget;
    }
    require(first_mixed > 0, "no mixed round under any budget: the checks below prove nothing");
    while (last_unmixed + 1 < first_mixed) {
        const uint64_t middle = last_unmixed + (first_mixed - last_unmixed) / 2;
        (void)run(middle, high);
        if (high) first_mixed = middle;
        else last_unmixed = middle;
    }
    for (uint64_t budget = last_unmixed; budget <= first_mixed + 2 * d::heuristics::secular_packet_sweeps + 2; ++budget)
        (void)run(budget, high);
    // A pair needs two sweeps: under a budget of one every round is an ordinary multiprecision packet.
    (void)run(1, high);
    require(high == 0, "a mixed round started with one sweep left");
    const auto unbudgeted = run(0, high);
    require(unbudgeted.complete && unbudgeted.goal_met, "the unbudgeted solve did not finish the clusters");
}

// An ordinary binary64 packet keeps the point of a node it does not move while
// the goal is the default one, and still writes back the point's rounding at a
// tighter goal (EnginePolicy::machine_keeps_unmoved). The node seeded at
// 1 + 2^-150 sits beyond binary64's bits and rounds to 1; the seed at 11/8
// (its root is 2) keeps the component uncertified, so a packet runs, and the
// binary64 packet is not refused.
static void binary64_packets_keep_points_at_the_default_goal() {
    const std::vector<Complex> roots{Complex(1), Complex(2), Complex(3)};
    for (const uint64_t goal : {cps_options{}.bits, uint64_t(333)}) {
        d::RootVector seeds(3);
        d::ComplexBall offset;
        acb_one(offset.value());
        acb_mul_2exp_si(offset.value(), offset.value(), -150);
        acb_one(seeds[0]);
        acb_add(seeds[0], seeds[0], offset.value(), 256);
        acb_set_d(seeds[1], 1.375);
        acb_set_si(seeds[2], 3);
        d::ComplexBall entry;
        acb_set(entry.value(), seeds[0]);
        auto policy = d::maintained_method_policy(CPS_METHOD_SECULAR);
        require(policy.machine_keeps_unmoved, "the maintained secular recipe does not keep unmoved points");
        // The certificate attempt and the preliminary phase would place the points first.
        policy.newton_point_refinement = policy.secular_preliminary = false;
        policy.study_statistics = std::make_shared<d::StudyStatistics>();
        const auto counts = policy.study_statistics;
        bool seen = false, kept = false, rounded = false;
        policy.snapshot = [&](d::NumericalSnapshot event) {
            if (seen || event.stage != "after_secular_packet") return;
            seen = true;
            kept = acb_equal(event.points[0].value(), entry.value());
            rounded = acb_is_one(event.points[0].value());
        };
        d::SecularRefinement state(source(roots), seeds.data(), 1, 85, policy, d::solve_local_model);
        cps_options options;
        options.bits = goal;
        options.max_work_bits = 4096;
        certpolysolve::detail::SearchStatistics search{};
        certpolysolve::detail::CandidateStatistics statistics{};
        for (unsigned call = 0; call < 64 && !seen && state.active(); ++call)
            state.step(options, search, statistics);
        require(seen, "no secular packet ran");
        require(counts->machine_packets > 0 && counts->continuations == 0, "the binary64 packet did not run");
        if (goal == cps_options{}.bits)
            require(kept, "a binary64 packet at the default goal rounded a node it did not move");
        else
            require(rounded, "a binary64 packet at a tighter goal kept a point it rounds today");
    }
}

// A node exactly on a root has residue zero, but with coefficients that are
// not dyadic its enclosure is rounding noise at every precision: it never
// gains the relative accuracy regeneration asks of it, and the precision
// doubled to the cap. The restart recipe at degree 96 rotated by (3 + 4i) / 5
// is not centred; its roots whose numerator is divisible by 5 are dyadic, and
// nodes reach them exactly (30 s at 65,536 bits). A residue whose enclosure,
// times n, lies below the node's resolution is adequate.
static void exact_roots_do_not_drive_the_precision() {
    constexpr int degree = 96, separation = 48;
    const Rational step = power2(-(separation + 4));
    const Complex turn(Rational(3, 5), Rational(4, 5));
    std::vector<Complex> roots;
    for (int j = 0; j < degree; ++j)
        roots.push_back(Complex(Rational(j % 4 - 2) + Rational(16 * (j + 1) + (7 * j + 3) % 16) * step) * turn);
    const auto polynomial = input(from_roots(roots));
    d::SolveState state(certpolysolve::detail::dense(*polynomial.get()), d::maintained_method_policy(CPS_METHOD_SECULAR));
    cps_options options;
    options.goal = CPS_ROOT_APPROXIMATION;
    options.bits = 53;
    options.max_work_bits = 65536;
    const auto result = state.advance(&options);
    require(result.complete && result.goal_met, "the rotated recipe was not solved");
    require(result.work_bits < options.max_work_bits, "residues at exact roots drove the precision to the cap");
}

static void restarted_members_keep_the_round_precision() {
    // Four roots 2^-60 apart at 1 among 60 on the integer grid; isolating
    // them takes a precision growth and a restart of their component.
    std::vector<Complex> roots;
    for (int x = -4; x <= 4; ++x)
        for (int y = -3; y <= 3; ++y)
            if (std::abs(x - 1) > 1 || std::abs(y) > 1)
                roots.emplace_back(x, y);
    for (int re : {-1, 1})
        for (int im : {-1, 1})
            roots.emplace_back(1 + re * power2(-60), im * power2(-60));
    const auto polynomial = input(from_roots(roots));
    auto policy = d::maintained_method_policy(CPS_METHOD_SECULAR);
    // The certificate attempt's cluster start resolves the four before the growth and restart this checks.
    policy.cluster_start = false;
    slong before_growth = 0, largest = 0, after_restart = 0;
    bool restarted = false;
    policy.snapshot = [&](d::NumericalSnapshot event) {
        slong bits = 0;
        for (size_t i = 0; i < event.coordinate_precision.size(); ++i)
            if (i >= event.inactive.size() || !event.inactive[i])
                bits = std::max(bits, event.coordinate_precision[i]);
        if (event.stage == "before_precision_growth" && !restarted)
            before_growth = bits;
        if (event.stage == "before_candidate_restart")
            restarted = true;
        if (event.stage == "regenerated_model" && restarted && !after_restart)
            after_restart = bits;
        largest = std::max(largest, bits);
    };
    d::SolveState state(certpolysolve::detail::dense(*polynomial.get()), policy);
    cps_options options;
    options.goal = CPS_ROOT_ACCURACY;
    options.bits = 53;
    options.max_work_bits = 4096;
    const auto result = state.advance(&options);
    require(result.goal_met, "the clustered control was not solved");
    require(restarted && before_growth > 0 && after_restart > 0, "the control no longer restarts a component");
    require(after_restart <= 2 * before_growth,
            "restarted members were given more iteration precision than the round's");
}

// Members iterated on the source (the classical recipe) take their restart model's precision for their own
// iteration, so the model is not raised for the round's resolution: no model puts its hints within its own
// precision, and where the raised model still does not resolve the cluster its members crawl at that
// precision (forced classical on mignotte_generalized_767_64_12 at the accuracy goal: 193 sweeps and 56 s
// became 1,856 sweeps and 104 s). x^95 - (3^40 x - 1)^3 has three roots about 2^-2070 apart at 3^-40, which no
// model of the first rounds resolves: every restart's model is at its round's precision, as before the rule.
// The scale is a power of 3 on purpose. About a centre that is a power of 2 (the control had 2^64) a head
// summed term by term is exact, the model then resolves the cluster in its first round and the precision
// follows the hints up, with no raise: 1,488 bits in the round at 680 (13 October 2026).
static void source_iterated_members_ask_no_resolution() {
    constexpr int degree = 95, scale = 40;
    Polynomial coefficients(degree + 1, Rational(0));
    Rational a = 1;
    for (int i = 0; i < scale; ++i)
        a *= 3;
    coefficients[degree] = 1;
    coefficients[3] -= a * a * a;
    coefficients[2] += 3 * a * a;
    coefficients[1] -= 3 * a;
    coefficients[0] += 1;
    const auto polynomial = input(coefficients);
    auto policy = d::maintained_method_policy(CPS_METHOD_CLASSICAL);
    struct Seen {
        uint64_t round = 0, restarts = 0, raised = 0;
    } seen;
    policy.trace_data = &seen;
    policy.trace = [](void *p, const d::TraceEvent &event) {
        auto &s = *static_cast<Seen *>(p);
        const std::string_view stage(event.stage), phase(event.phase);
        if (phase != "begin")
            return;
        if (stage == "secular_candidate_restart")
            s.round = event.precision;
        else if (stage == "local_candidates" && s.round) {
            ++s.restarts;
            s.raised += event.precision > s.round;
        }
    };
    d::SolveState state(certpolysolve::detail::dense(*polynomial.get()), policy);
    cps_options options;
    options.goal = CPS_ROOT_ACCURACY;
    options.bits = 53;
    // Enough for the first restarts, not for the cluster.
    options.max_work_bits = 2048;
    state.advance(&options);
    require(seen.restarts >= 2, "the control no longer restarts its cluster");
    require(seen.raised == 0, "a model of members iterated on the source was raised for the round's resolution");
}

// A restart's model is sized for the goal the caller asked for, of which a default goal staged first is a part.
// Four roots 2^-60 apart among ninety, at a goal of 1,000 bits: the default goal is staged, and the stage's
// restart builds its model at the regeneration's precision, which resolves the four. Sized for the stage's own
// goal it had 320 bits, the four stayed one component within the default tolerance, the certificate attempt at
// the caller's goal could not start from it, and the refinement resolved it through every round on the way
// (unisolve's lsr_200 at 3,322 bits: 0.41 -> 2.2 s).
static void restart_models_are_sized_for_the_callers_goal() {
    std::vector<Complex> roots;
    for (int x = -5; x <= 5; ++x)
        for (int y = -4; y <= 4; ++y)
            if (std::abs(x - 1) > 1 || std::abs(y) > 1)
                roots.emplace_back(x, y);
    for (int re : {-1, 1})
        for (int im : {-1, 1})
            roots.emplace_back(1 + re * power2(-60), im * power2(-60));
    const auto polynomial = input(from_roots(roots));
    auto policy = d::maintained_method_policy(CPS_METHOD_SECULAR);
    constexpr uint64_t goal = 1000;
    struct Seen {
        uint64_t models = 0, least = UINT64_MAX;
    } seen;
    policy.trace_data = &seen;
    policy.trace = [](void *p, const d::TraceEvent &event) {
        auto &s = *static_cast<Seen *>(p);
        if (std::string_view(event.stage) == "local_candidates" && std::string_view(event.phase) == "begin") {
            ++s.models;
            s.least = std::min(s.least, event.precision);
        }
    };
    d::SolveState state(certpolysolve::detail::dense(*polynomial.get()), policy);
    cps_options options;
    options.goal = CPS_ROOT_APPROXIMATION;
    options.bits = goal;
    options.max_work_bits = 16384;
    const auto result = state.advance(&options);
    require(result.goal_met, "the clustered control was not solved at its goal");
    require(seen.models > 0, "the control no longer restarts its cluster");
    require(seen.least >= goal, "a restart's model was sized for the staged default goal, not the caller's");
}

// A model raised once is not started low again. Two clusters of three and of two roots 2^-1200 apart among
// ninety on the integer grid, under the accuracy goal: every precision growth restarts both, and a model that
// starts at the regeneration's precision resolves neither within its round's resolution. Each growth's first
// model was raised to its sizing (three bits of model for each bit of the resolution, and the bits the
// cluster's evaluation cancels, which the growth before had shown and which were forgotten), so every growth
// paid for two models of that cluster (unisolve's lsr1: 14 whole shifts for 9 restarts, 0.79 -> 1.00 s at this
// goal). After the first raise the models start at their sizing.
static void restart_models_are_not_raised_twice() {
    std::vector<Complex> roots;
    for (int x = -5; x <= 5; ++x)
        for (int y = -4; y <= 4; ++y)
            if ((std::abs(x - 1) > 1 || std::abs(y) > 1) && (std::abs(x + 3) > 1 || std::abs(y) > 1))
                roots.emplace_back(x, y);
    const Rational gap = power2(-1200);
    for (int j = -1; j <= 1; ++j)
        roots.emplace_back(Rational(1) + Rational(j) * gap, Rational(0));
    for (int j = 0; j <= 1; ++j)
        roots.emplace_back(Rational(-3) + Rational(j) * gap, Rational(0));
    const auto polynomial = input(from_roots(roots));
    auto policy = d::maintained_method_policy(CPS_METHOD_SECULAR);
    // The certificate attempt's cluster start resolves both clusters before the growths this checks.
    policy.cluster_start = false;
    struct Seen {
        uint64_t models = 0, raises = 0;
    } seen;
    policy.trace_data = &seen;
    policy.trace = [](void *p, const d::TraceEvent &event) {
        auto &s = *static_cast<Seen *>(p);
        const std::string_view stage(event.stage), phase(event.phase);
        if (stage == "local_candidates" && phase == "begin")
            ++s.models;
        else if (stage == "local_model_raised")
            ++s.raises;
    };
    d::SolveState state(certpolysolve::detail::dense(*polynomial.get()), policy);
    cps_options options;
    options.goal = CPS_ROOT_ACCURACY;
    options.bits = 53;
    options.max_work_bits = 65536;
    const auto result = state.advance(&options);
    require(result.goal_met, "the two clusters were not resolved");
    require(seen.models >= 8, "the control no longer restarts its clusters in several growths");
    require(seen.raises <= 1, "a restart's model was started below its sizing again after one had been raised");
}

// What a restart's models record, for the two that follow. The trace's component_initialization scope opens at
// the regeneration's precision; each model in it is a local_candidates scope at the model's.
struct RestartRecord {
    struct Model {
        uint64_t members, entry, precision;
    };
    std::vector<Model> models;
    uint64_t entry = 0, raises = 0;
    static void trace(void *p, const d::TraceEvent &event) {
        auto &s = *static_cast<RestartRecord *>(p);
        const std::string_view stage(event.stage), phase(event.phase);
        if (stage == "component_initialization" && phase == "begin")
            s.entry = event.precision;
        else if (stage == "local_candidates" && phase == "begin")
            s.models.push_back({event.active, s.entry, event.precision});
        else if (stage == "local_model_raised")
            ++s.raises;
    }
};

// Ninety roots on the integer grid with a cluster of three roots `triple` apart and one of two roots 2^-1200
// apart in two of its cells, the three at -3 where `triple_first` (the restart then takes it before the pair),
// all moved by `offset` (3 + 4i): the fixture of restart_models_are_not_raised_twice and its variants.
static std::vector<Complex> grid_with_two_clusters(const Rational &triple, long offset, bool triple_first) {
    std::vector<Complex> roots;
    for (int x = -5; x <= 5; ++x)
        for (int y = -4; y <= 4; ++y)
            if ((std::abs(x - 1) > 1 || std::abs(y) > 1) && (std::abs(x + 3) > 1 || std::abs(y) > 1))
                roots.emplace_back(x, y);
    for (int j = -1; j <= 1; ++j)
        roots.emplace_back(Rational(triple_first ? -3 : 1) + Rational(j) * triple, Rational(0));
    for (int j = 0; j <= 1; ++j)
        roots.emplace_back(Rational(triple_first ? 1 : -3) + Rational(j) * power2(-1200), Rational(0));
    for (auto &root : roots)
        root = Complex(root.real + Rational(3 * offset), root.imag + Rational(4 * offset));
    return roots;
}
static RestartRecord accuracy_restarts(const std::vector<Complex> &roots, bool centred) {
    const auto polynomial = input(from_roots(roots));
    auto policy = d::maintained_method_policy(CPS_METHOD_SECULAR);
    policy.complex_centering = centred;
    // The certificate attempt's cluster start resolves the clusters before the restarts these check.
    policy.cluster_start = false;
    RestartRecord seen;
    policy.trace_data = &seen;
    policy.trace = RestartRecord::trace;
    d::SolveState state(certpolysolve::detail::dense(*polynomial.get()), policy);
    cps_options options;
    options.goal = CPS_ROOT_ACCURACY;
    options.bits = 53;
    options.max_work_bits = 65536;
    const auto result = state.advance(&options);
    require(result.goal_met, "the two clusters were not resolved");
    return seen;
}

// The bits a raise showed are the solve's. The two clusters 768 + 1024 i from the origin and not centred, as an
// input whose centring is declined: an evaluation there cancels some 700 bits, the three-root cluster's first
// model is raised for them, and every later model is sized with them. Kept only for the growth they were shown
// in, they reached the pair as long as every growth raised the three's model anew; once that cluster started at
// its sizing (restart_models_are_not_raised_twice) the pair's next model was bare. Here that model is raised
// in its turn; on unisolve's lsr1 shifted by the square root of 2 it was not (704 bits after 1,152, its hints
// too coarse for the round), the pair was no component in the growth after, and 24 packets carried the two
// nodes in: 438 -> 614 sweeps.
static void restart_models_keep_the_bits_a_raise_showed() {
    constexpr long offset = 256;
    const auto seen = accuracy_restarts(grid_with_two_clusters(power2(-1200), offset, true), false);
    uint64_t pairs = 0, last = 0;
    for (const auto &model : seen.models)
        if (model.members == 2) {
            require(model.precision >= last, "a cluster's model was coarser than its model of the growth before");
            last = model.precision;
            ++pairs;
        }
    require(pairs >= 4, "the control no longer restarts its pair in several growths");
    require(seen.raises <= 1, "a model was sized without the bits a raise of the solve had shown");
}

// A cluster its last model left unresolved starts at its sizing. The three roots 2^-500 apart are resolved a
// growth before the pair, which is then restarted alone: the regeneration's precision is its own, below the
// pair's sizing, where the three's models had kept it above (2,213 bits for 2,816). The pair's models have
// determined none of its hints so far, so a model at 2,213 bits cannot either; started there it was raised, a
// whole shift more on a source with no dense head (unisolve's lsr1 rotated by i: 57 ms of 800).
static void unresolved_restart_components_start_at_their_sizing() {
    const auto seen = accuracy_restarts(grid_with_two_clusters(power2(-500), 0, false), true);
    require(!seen.models.empty() && seen.models.back().members == 2 &&
                seen.models.back().precision > seen.models.back().entry,
            "the control's pair is no longer restarted alone below its sizing");
    require(seen.raises <= 1, "a cluster its last model left unresolved was started below its sizing");
}

// A component with no record starts at the regeneration's precision. Twenty-one roots on a spiral 2^-10 wide
// about -1 and four others, under the accuracy goal: the binary64 stage leaves the cluster one component, and
// the precision of its first restart determines every root of the model. A rule that started such a model at
// its sizing where that was less than twice the regeneration's precision (21 bits of model for each bit of the
// resolution) took 3,584 bits for 1,872 on unisolve's spiral25, 205 ms for 75 of a solve of 0.5 s.
static void unrecorded_restart_models_start_at_the_regenerations_precision() {
    std::vector<Complex> roots{Complex(1, 0), Complex(2, 1), Complex(0, 2), Complex(3, -1)};
    Rational real(1), imaginary(0);
    for (int j = 1; j <= 21; ++j) {
        const Rational turned = (Rational(3) * real - Rational(4) * imaginary) / Rational(5);
        imaginary = (Rational(4) * real + Rational(3) * imaginary) / Rational(5);
        real = turned;
        roots.emplace_back(Rational(-1) + power2(-10 * j) * real, power2(-10 * j) * imaginary);
    }
    const auto polynomial = input(from_roots(roots));
    auto policy = d::maintained_method_policy(CPS_METHOD_SECULAR);
    RestartRecord seen;
    policy.trace_data = &seen;
    policy.trace = RestartRecord::trace;
    d::SolveState state(certpolysolve::detail::dense(*polynomial.get()), policy);
    cps_options options;
    options.goal = CPS_ROOT_ACCURACY;
    options.bits = 53;
    options.max_work_bits = 65536;
    const auto result = state.advance(&options);
    require(result.goal_met, "the spiral was not resolved");
    require(!seen.models.empty() && seen.models.front().members > 10,
            "the control's first restart is no longer its cluster's");
    require(seen.raises == 0, "the control's models are raised: the regeneration's precision no longer resolves it");
    require(seen.models.front().precision <= seen.models.front().entry,
            "a component with no record was started above the regeneration's precision");
}

// Restart hints the local model does not determine are placed by their scale,
// with no conjugate symmetry (place_undetermined_hints). Fuzz seed
// 10295942349935689529: a six-root real cluster 2^-59 wide with a nested pair
// 2^-104 apart, and a seventh root at -1/1024. The local model truncated at
// degree six has a conjugate pair m +- 2^-85 i where the two close roots are
// (the dropped term is 2^-413 on the cluster, the model's value between the
// pair 2^-415), and a Jacobi sweep of a real source keeps a conjugate pair
// conjugate: iterated at the round's precision the pair sat on the bisector of
// the two roots until the work budget (400 units and still running); placed
// on their circle at the placement angles they converge at once.
static void undetermined_restart_hints_are_placed_by_scale() {
    const std::vector<Complex> roots{
        Complex(Rational("-1/1024")),
        Complex(Rational("708191041404007847/509014844283922939904")),
        Complex(Rational("177047760351001741/127253711070980734976")),
        Complex(Rational("354095520702004365/254507422141961469952")),
        Complex(Rational("177047760351003507/127253711070980734976")),
        Complex(Rational("708191041404011379/509014844283922939904")),
        Complex(Rational("24917257110736165240970126820211/17909367680024424984345422885552128"))};
    const auto polynomial = input(from_roots(roots));
    auto policy = d::maintained_method_policy(CPS_METHOD_SECULAR);
    unsigned placements = 0;
    policy.trace_data = &placements;
    policy.trace = [](void *p, const d::TraceEvent &event) {
        if (std::string_view(event.stage) == "local_hints_undetermined")
            ++*static_cast<unsigned *>(p);
    };
    d::SolveState state(certpolysolve::detail::dense(*polynomial.get()), policy);
    cps_options options;
    options.goal = CPS_ROOT_APPROXIMATION;
    options.bits = 100;
    options.max_work_bits = 1024;
    options.work_units = 100; // the solve takes 20; with the conjugate pair it spent 400 and ran on
    const auto result = state.advance(&options);
    require(result.goal_met && result.complete && result.termination == CPS_SUCCESS,
            "a conjugate pair of restart hints for two close real roots was not separated");
    require(placements > 0, "the nested-pair control no longer has restart hints its local model leaves undetermined");
    uint64_t total = 0;
    for (const auto &root : result.roots) total += root.count;
    require(total == roots.size(), "the nested-pair control lost roots");
}

// Two hints are told apart only if no admissible polynomial merges them, and they merge at the model's critical
// point between them, which a third root near moves off their midpoint (place_undetermined_hints). Three real
// roots 2^-60 apart with a fourth 2^-42 away, restarted as one component of four (jacobi-characteristic-3's
// shape): about the four's mean the dropped fifth-degree term is 12.6 gaps cubed on the three, and the model
// has a real root 2.47 gaps out and a conjugate pair at -1.24 +- 1.89 i in their place. Its value is 13.3 at
// the pair's midpoint, beyond the dropped term, and 12.2 at its critical point 0.66 gaps away, within it. Left
// as a conjugate pair, the two are kept conjugate by a real source's sweeps and never reach two real roots.
static void restart_hints_that_would_merge_are_placed_by_scale() {
    std::vector<Complex> roots;
    for (int x = -4; x <= 4; ++x)
        for (int y = -4; y <= 4; ++y)
            if (std::abs(x - 1) > 1 || std::abs(y) > 1)
                roots.emplace_back(x, y);
    constexpr int gap = 60, far = 42;
    for (int k : {-1, 0, 1})
        roots.emplace_back(1 + k * power2(-gap));
    roots.emplace_back(1 - power2(-far));
    auto polynomial = source(roots);
    const slong working = 704;
    d::Root component;
    // About the four's mean, 2^-44 from the three.
    enclose(component.box.value(), Complex(1 - power2(-far - 2)), working);
    require(acb_is_exact(component.box.value()), "the box's middle is not exact");
    arb_add_error_2exp_si(acb_realref(component.box.value()), -far + 1);
    arb_add_error_2exp_si(acb_imagref(component.box.value()), -far + 1);
    component.count = 4;
    cps_options options;
    options.max_work_bits = 65536;
    certpolysolve::detail::CandidateStatistics statistics{};
    d::EnginePolicy policy;
    unsigned placements = 0;
    policy.trace_data = &placements;
    policy.trace = [](void *p, const d::TraceEvent &event) {
        if (std::string_view(event.stage) == "local_hints_undetermined")
            ++*static_cast<unsigned *>(p);
    };
    d::RootVector hints(4);
    slong used = working;
    require(d::component_seeds(hints.data(), *polynomial, component, 4, used, options, statistics, policy,
                               std::nullopt) && used == working,
            "the restart supplied no hints at its precision");
    require(placements == 1, "a pair of hints that an admissible perturbation merges was not placed by scale");
    // No hint is another's conjugate: within a gap's 2^-64, where the model's pair was conjugate to its last bit.
    d::Magnitude tolerance, distance;
    mag_set_ui_2exp_si(tolerance.value(), 1, -gap - 64);
    for (slong i = 0; i < 4; ++i)
        for (slong j = i + 1; j < 4; ++j) {
            d::ComplexBall mirrored;
            acb_conj(mirrored.value(), hints.data() + j);
            acb_sub(mirrored.value(), mirrored.value(), hints.data() + i, 4 * working);
            acb_get_mag(distance.value(), mirrored.value());
            require(arb_is_zero(acb_imagref(hints.data() + i)) || mag_cmp(distance.value(), tolerance.value()) > 0,
                    "two restart hints of a real cluster are a conjugate pair");
        }
}

// A node just outside the counting circle: its series in the count's
// exterior expansion converges at the rate q = rho / |b - c|, so a node at
// rho (1 + 2^-12) asked for about 2^12 times the bits of the expansion's
// resolution in terms (about 120,000 here; a node 2^-40 outside, 10^13: the
// count's time and memory had no bound). Nodes within twice the radius are
// now evaluated directly on each arc, as the members are, and the expansion
// holds only nodes whose series gain a bit a term.
static void exact_count_near_node() {
    const slong precision = 256;
    const auto dyadic = [](d::ComplexBall &z, slong re_num, slong re_exp, slong im_num, slong im_exp) {
        arb_set_si(acb_realref(z.value()), re_num);
        arb_mul_2exp_si(acb_realref(z.value()), acb_realref(z.value()), re_exp);
        arb_set_si(acb_imagref(z.value()), im_num);
        arb_mul_2exp_si(acb_imagref(z.value()), acb_imagref(z.value()), im_exp);
    };
    // Roots: three within 2^-10 of the origin, and 2, 3i, -5/2.
    std::vector<d::ComplexBall> roots(6), nodes(6);
    dyadic(roots[0], 1, -10, 0, 0);
    dyadic(roots[1], -1, -10, 1, -10);
    dyadic(roots[2], 0, 0, -1, -10);
    dyadic(roots[3], 2, 0, 0, 0);
    dyadic(roots[4], 0, 0, 3, 0);
    dyadic(roots[5], -5, -1, 0, 0);
    // Nodes: the cluster's three at their roots (offset by 2^-30), one at
    // 1 + 2^-30 (serving the root 2 from afar), two at 3i and -5/2 (offset).
    for (int k = 0; k < 3; ++k) {
        acb_set(nodes[k].value(), roots[k].value());
        d::ComplexBall nudge;
        dyadic(nudge, 1, -30, 1, -31);
        acb_add(nodes[k].value(), nodes[k].value(), nudge.value(), precision);
    }
    dyadic(nodes[3], (1 << 30) + 1, -30, 0, 0);
    for (int k = 4; k < 6; ++k) {
        acb_set(nodes[k].value(), roots[k].value());
        d::ComplexBall nudge;
        dyadic(nudge, 1, -30, -1, -31);
        acb_add(nodes[k].value(), nodes[k].value(), nudge.value(), precision);
    }
    std::vector<d::ComplexBall> w(nodes.size());
    for (size_t i = 0; i < nodes.size(); ++i) {
        d::ComplexBall value, denominator, factor;
        acb_one(value.value());
        for (const auto &r : roots) {
            acb_sub(factor.value(), nodes[i].value(), r.value(), precision);
            acb_mul(value.value(), value.value(), factor.value(), precision);
        }
        acb_one(denominator.value());
        for (size_t j = 0; j < nodes.size(); ++j) {
            if (j == i) continue;
            acb_sub(factor.value(), nodes[i].value(), nodes[j].value(), precision);
            acb_mul(denominator.value(), denominator.value(), factor.value(), precision);
        }
        acb_div(w[i].value(), value.value(), denominator.value(), precision);
        acb_neg(w[i].value(), w[i].value());
    }
    // The circle |z| = 1 - 2^-12: the node at 1 + 2^-30 lies 2^-12 outside.
    d::ComplexBall centre;
    d::Float radius;
    arf_set_si(radius.value(), (1 << 12) - 1);
    arf_mul_2exp_si(radius.value(), radius.value(), -12);
    std::vector<uint8_t> members{1, 1, 1, 0, 0, 0};
    d::DiskCountScratch scratch;
    require(d::secular_disk_count(nodes, w, members, centre.value(), radius.value(), scratch) ==
                std::optional<slong>(3),
            "the model did not count the three roots inside a circle with a node just outside");
    require(scratch.moments.size() <= 128, "the count's exterior expansion grew with a node near the circle");
}

// Five nodes inside a four-root cluster A and three in a four-root cluster B
// (and six isolated roots, each with its node): the surplus node serves no
// root of A, so its disk reaches B's unserved root and the inclusion cover's
// component spans both clusters, a region no goal tighter than their
// distance accepts. The secular model counts the roots in a disk exactly
// (secular_disk_count: four in A under five nodes; five under five when A
// holds five roots; four under four when balanced), and at a goal the
// clusters meet the cover replaces that component by squares about each
// cluster, certified by those counts (clustered-769 at 10 bits).
static void exact_count_split() {
    const slong precision = 256;
    const auto dyadic = [](d::ComplexBall &z, slong re_num, slong re_exp, slong im_num, slong im_exp) {
        arb_set_si(acb_realref(z.value()), re_num);
        arb_mul_2exp_si(acb_realref(z.value()), acb_realref(z.value()), re_exp);
        arb_set_si(acb_imagref(z.value()), im_num);
        arb_mul_2exp_si(acb_imagref(z.value()), acb_imagref(z.value()), im_exp);
    };
    // Roots: clusters at A = 3 + 2i and B = -5 + 7i, each four roots at
    // distance 2^-20 on the axes; six isolated roots.
    std::vector<d::ComplexBall> roots(14);
    const slong spread = -20;
    const slong axes[4][2] = {{1, 0}, {0, 1}, {-1, 0}, {0, -1}};
    for (int k = 0; k < 4; ++k) {
        dyadic(roots[k], 3 * (1 << 20) + axes[k][0], spread, 2 * (1 << 20) + axes[k][1], spread);
        dyadic(roots[4 + k], -5 * (1 << 20) + axes[k][0], spread, 7 * (1 << 20) + axes[k][1], spread);
    }
    const slong isolated[6][2] = {{9, -4}, {-8, -6}, {1, 11}, {12, 5}, {-11, 2}, {6, -10}};
    for (int k = 0; k < 6; ++k) dyadic(roots[8 + k], isolated[k][0], 0, isolated[k][1], 0);
    // The secular model of the monic polynomial with the given zeros at the
    // nodes: w_i = -p(b_i) / Q'(b_i).
    const auto residues_for = [&](const std::vector<d::ComplexBall> &nodes, const std::vector<d::ComplexBall> &zeros) {
        std::vector<d::ComplexBall> w(nodes.size());
        for (size_t i = 0; i < nodes.size(); ++i) {
            d::ComplexBall value, denominator, factor;
            acb_one(value.value());
            for (const auto &r : zeros) {
                acb_sub(factor.value(), nodes[i].value(), r.value(), precision);
                acb_mul(value.value(), value.value(), factor.value(), precision);
            }
            acb_one(denominator.value());
            for (size_t j = 0; j < nodes.size(); ++j) {
                if (j == i) continue;
                acb_sub(factor.value(), nodes[i].value(), nodes[j].value(), precision);
                acb_mul(denominator.value(), denominator.value(), factor.value(), precision);
            }
            acb_div(w[i].value(), value.value(), denominator.value(), precision);
            acb_neg(w[i].value(), w[i].value());
        }
        return w;
    };
    // Five nodes inside A at the cluster's scale (off its roots), three at
    // three of B's roots (offset by 2^-40), six at the isolated roots.
    std::vector<d::ComplexBall> nodes(14);
    const slong inner[5][2] = {{3, 1}, {-2, 3}, {-3, -2}, {1, -3}, {0, 0}};
    for (int k = 0; k < 5; ++k)
        dyadic(nodes[k], 3 * (1 << 22) + inner[k][0], spread - 2, 2 * (1 << 22) + inner[k][1], spread - 2);
    for (int k = 0; k < 3; ++k) {
        acb_set(nodes[5 + k].value(), roots[4 + k].value());
        d::ComplexBall nudge;
        dyadic(nudge, 1, -40, 1, -41);
        acb_add(nodes[5 + k].value(), nodes[5 + k].value(), nudge.value(), precision);
    }
    for (int k = 0; k < 6; ++k) {
        acb_set(nodes[8 + k].value(), roots[8 + k].value());
        d::ComplexBall nudge;
        dyadic(nudge, 1, -30, -1, -31);
        acb_add(nodes[8 + k].value(), nodes[8 + k].value(), nudge.value(), precision);
    }
    // The same nodes against a polynomial with five roots in A (B's fourth
    // root moved into A's disk, off every node), and a balanced population.
    auto five = roots;
    dyadic(five[7], 3 * (1 << 22) + 2, spread - 2, 2 * (1 << 22) - 1, spread - 2);
    std::vector<d::ComplexBall> balanced(nodes.begin(), nodes.begin() + 4);
    for (int k = 0; k < 4; ++k) {
        balanced.push_back(roots[4 + k]);
        d::ComplexBall nudge;
        dyadic(nudge, 1, -40, 1, -41);
        acb_add(balanced.back().value(), balanced.back().value(), nudge.value(), precision);
    }
    balanced.insert(balanced.end(), nodes.begin() + 8, nodes.end());
    const auto overfull = residues_for(nodes, roots), full = residues_for(nodes, five),
               populated = residues_for(balanced, roots);
    // The count on disks about A's centre from 8 times the cluster's scale up
    // to radius 4 (the nearest other node is about 8.5 away).
    d::ComplexBall centre;
    dyadic(centre, 3, 0, 2, 0);
    std::vector<uint8_t> five_nodes(nodes.size(), 0), four_nodes(balanced.size(), 0);
    std::fill(five_nodes.begin(), five_nodes.begin() + 5, uint8_t(1));
    std::fill(four_nodes.begin(), four_nodes.begin() + 4, uint8_t(1));
    d::DiskCountScratch scratch;
    for (slong exponent = spread + 3; exponent <= 2; ++exponent) {
        d::Float radius;
        arf_one(radius.value());
        arf_mul_2exp_si(radius.value(), radius.value(), exponent);
        require(d::secular_disk_count(nodes, overfull, five_nodes, centre.value(), radius.value(), scratch) ==
                    std::optional<slong>(4),
                "the model did not count four roots under five nodes");
        require(d::secular_disk_count(nodes, full, five_nodes, centre.value(), radius.value(), scratch) ==
                    std::optional<slong>(5),
                "the model did not count five roots under five nodes");
        require(d::secular_disk_count(balanced, populated, four_nodes, centre.value(), radius.value(), scratch) ==
                    std::optional<slong>(4),
                "the model did not count four roots under four nodes");
    }
    // A disk that holds a node outside the group has no count.
    {
        d::Float radius;
        arf_set_si(radius.value(), 9);
        require(!d::secular_disk_count(nodes, overfull, five_nodes, centre.value(), radius.value(), scratch),
                "a disk holding a non-member node was counted");
    }
    // The cover at an 8-bit approximation goal: squares about A and B with
    // four roots each, every region meeting the goal and holding exactly its
    // count of the planted roots.
    cps_options options;
    options.goal = CPS_ROOT_APPROXIMATION;
    options.bits = 8;
    std::vector<d::Root> cover;
    std::vector<uint8_t> inactive;
    require(d::secular_inclusion_cover(cover, inactive, nodes, overfull, false, 1, precision, options),
            "the inclusion cover failed");
    uint64_t total = 0, counted = 0;
    for (const auto &region : cover) {
        require(d::meets_goal(region, options, precision), "a region of the split cover misses the goal");
        uint64_t inside = 0;
        for (const auto &r : roots) inside += acb_contains(region.box.value(), r.value());
        require(inside == region.count, "a region does not hold its count of planted roots");
        total += region.count;
        counted += region.certificate == "secular-model-count" && region.count == 4;
    }
    require(total == roots.size() && counted == 2 &&
                std::all_of(inactive.begin(), inactive.end(), [](uint8_t v) { return v != 0; }),
            "the clusters were not certified by their exact counts");
}

// A component that mixes kept certificates with nodes still uncertified keeps the kept nodes inactive (where they
// were certified) and the others active; the component is one region holding every member's root. Two kept nodes at
// 0 and 1 (regions of radius 2^-60) and two far ones at 5 and 6 whose residues (10) give disks of radius 40.
static void kept_nodes_stay_inactive_in_a_mixed_component() {
    const slong precision = 128;
    std::vector<d::ComplexBall> nodes(4), residues(4);
    const double at[] = {0, 1, 5, 6};
    const double w[] = {0x1p-70, 0x1p-70, 10, 10};
    std::vector<d::Root> kept_roots(4);
    std::vector<uint8_t> kept = {1, 1, 0, 0};
    for (size_t i = 0; i < 4; ++i) {
        acb_set_d(nodes[i].value(), at[i]);
        acb_set_d(residues[i].value(), w[i]);
        if (kept[i]) {
            acb_set_d(kept_roots[i].box.value(), at[i]);
            d::Magnitude r;
            mag_set_d(r.value(), 0x1p-60);
            acb_add_error_mag(kept_roots[i].box.value(), r.value());
            kept_roots[i].count = kept_roots[i].multiplicity = 1;
            kept_roots[i].publish(precision);
        }
    }
    cps_options options;
    std::vector<d::Root> cover;
    std::vector<uint8_t> inactive;
    require(d::secular_inclusion_cover(cover, inactive, nodes, residues, false, 1, precision, options, {}, kept_roots, kept),
            "the inclusion cover with kept regions failed");
    require(cover.size() == 1 && cover[0].count == 4, "the mixed component is not one region of four roots");
    require(inactive == std::vector<uint8_t>({1, 1, 0, 0}),
            "a kept node of a mixed component was made active, or an uncertified one inactive");
}

// The other side of that rule: a kept node whose own disk reaches a converged close pair (their disks small, none
// reaching it) iterates with the component, which restarts whole. Held there, the pair's restart counted 3 roots for
// 2 members and never ran: seven fuzz solves stopped at their caps (28 October 2026).
static void kept_node_not_reached_stays_active() {
    const slong precision = 128;
    std::vector<d::ComplexBall> nodes(3), residues(3);
    const double at[] = {0, 1, 1 + 0x1p-30};
    const double w[] = {1, 0x1p-70, 0x1p-70};
    std::vector<d::Root> kept_roots(3);
    std::vector<uint8_t> kept = {1, 0, 0};
    for (size_t i = 0; i < 3; ++i) {
        acb_set_d(nodes[i].value(), at[i]);
        acb_set_d(residues[i].value(), w[i]);
    }
    acb_set_d(kept_roots[0].box.value(), 0);
    d::Magnitude r;
    mag_set_d(r.value(), 0x1p-60);
    acb_add_error_mag(kept_roots[0].box.value(), r.value());
    kept_roots[0].count = kept_roots[0].multiplicity = 1;
    kept_roots[0].publish(precision);
    cps_options options;
    std::vector<d::Root> cover;
    std::vector<uint8_t> inactive;
    require(d::secular_inclusion_cover(cover, inactive, nodes, residues, false, 1, precision, options, {}, kept_roots, kept),
            "the inclusion cover with a kept region failed");
    require(cover.size() == 1 && cover[0].count == 3, "the kept node and the pair are not one region of three roots");
    require(inactive == std::vector<uint8_t>({0, 0, 0}),
            "a kept node that no uncertified disk reaches was held in its component");
}

// OutwardBounds::apart may only say "apart" of rectangles acb_overlaps finds disjoint, also where an endpoint
// is beyond binary64's range, and must leave touching and unbounded rectangles to the exact test.
static void outward_bounds_never_separate_overlapping_rectangles() {
    const auto ball = [](slong mantissa, slong exponent, slong radius_mantissa, slong radius_exponent) {
        d::ComplexBall z;
        arb_set_si(acb_realref(z.value()), mantissa);
        arb_mul_2exp_si(acb_realref(z.value()), acb_realref(z.value()), exponent);
        d::Magnitude radius;
        mag_set_ui(radius.value(), static_cast<ulong>(radius_mantissa));
        mag_mul_2exp_si(radius.value(), radius.value(), radius_exponent);
        acb_add_error_mag(z.value(), radius.value());
        return z;
    };
    std::vector<d::ComplexBall> boxes;
    for (slong exponent : {-1100, -1074, -60, 0, 60, 1023, 1100})
        for (slong mantissa : {-5, -3, 3, 4, 5})
            for (slong radius : {0, 1, 2})
                boxes.push_back(ball(mantissa, exponent, radius, exponent));
    // Midpoints that agree beyond binary64 with radii that do and do not reach each other.
    for (slong radius_exponent : {-200, -100, -99}) {
        boxes.push_back(ball(1, 0, 1, radius_exponent));
        d::ComplexBall shifted = ball(1, -100, 0, 0);
        arb_add_si(acb_realref(shifted.value()), acb_realref(shifted.value()), 1, 400);
        d::Magnitude radius;
        mag_set_ui_2exp_si(radius.value(), 1, radius_exponent);
        acb_add_error_mag(shifted.value(), radius.value());
        boxes.push_back(std::move(shifted));
    }
    {
        d::ComplexBall unbounded = ball(3, 0, 0, 0);
        mag_inf(arb_radref(acb_realref(unbounded.value())));
        boxes.push_back(std::move(unbounded));
    }
    size_t separated = 0;
    for (const auto &a : boxes)
        for (const auto &b : boxes) {
            // The same test with the axes exchanged covers the imaginary bounds.
            d::ComplexBall ia, ib;
            acb_mul_onei(ia.value(), a.value());
            acb_mul_onei(ib.value(), b.value());
            const bool apart = d::OutwardBounds(a.value()).apart(d::OutwardBounds(b.value()));
            require(!apart || !acb_overlaps(a.value(), b.value()), "binary64 bounds separated overlapping rectangles");
            require(apart == d::OutwardBounds(ia.value()).apart(d::OutwardBounds(ib.value())),
                    "the imaginary bounds disagree with the real ones");
            separated += apart;
        }
    require(separated > boxes.size(), "the bounds separated almost nothing");
}

// A polynomial of few terms: its head is summed term by term (PolynomialSource::taylor_head), agrees with the
// whole shift's first coefficients to the precision, and costs a fraction of a pass.
static void expansion_head_of_few_terms() {
    // Mignotte's x^n - 2 (a x - 1)^2 with n = 257, a = 2^14: four terms and a pair of roots 2^-1800 apart near 1 / a.
    const slong n = 257, working = 4096;
    d::Integer a, coefficient;
    fmpz_one(a.value());
    fmpz_mul_2exp(a.value(), a.value(), 14);
    fmpz_poly_t integer;
    fmpz_poly_init(integer);
    fmpz_poly_set_coeff_si(integer, n, 1);
    fmpz_poly_set_coeff_si(integer, 0, -2);
    fmpz_mul_si(coefficient.value(), a.value(), 4);
    fmpz_poly_set_coeff_fmpz(integer, 1, coefficient.value());
    fmpz_mul(coefficient.value(), a.value(), a.value());
    fmpz_mul_si(coefficient.value(), coefficient.value(), -2);
    fmpz_poly_set_coeff_fmpz(integer, 2, coefficient.value());
    d::PolynomialSource polynomial(
        [integer](acb_poly_t out, slong precision) { acb_poly_set_fmpz_poly(out, integer, precision); }, n, true);
    const auto three = polynomial.taylor_head_passes(working, 3), one = polynomial.taylor_head_passes(working, 1);
    require(three && one && *three < 0.5 && *one < *three, "a four-term polynomial's head is costed in whole passes");
    for (const auto &[re, im] : {std::pair{0x1p-14, 0.0}, std::pair{0x1p-14 + 0x1p-40, 0x1p-45}, std::pair{0.75, 0.75},
                                  // Whole binary64 mantissas: no power of this centre is exact at the precision.
                                  std::pair{1.0 / 3, 0x1p-14 / 3}}) {
        d::ComplexBall center;
        acb_set_d_d(center.value(), re, im);
        d::ComplexPolynomial whole, head;
        const uint64_t heads = polynomial.taylor_heads();
        require(polynomial.taylor(whole.value(), center.value(), working) &&
                    polynomial.taylor_head(head.value(), center.value(), working, 3) &&
                    acb_poly_length(head.value()) == 3 && polynomial.taylor_heads() == heads + 1,
                "no term-by-term head about the centre");
        for (slong j = 0; j < 3; ++j) {
            require(acb_overlaps(head.value()->coeffs + j, whole.value()->coeffs + j),
                    "the term-by-term head differs from the whole shift's first coefficients");
            require(acb_rel_accuracy_bits(head.value()->coeffs + j) >=
                        acb_rel_accuracy_bits(whole.value()->coeffs + j) - 2,
                    "the term-by-term head is less accurate than the whole shift");
        }
    }
    fmpz_poly_clear(integer);
}

// The dropped terms' bound from magnitudes (taylor_tail_bound) against the tail of the whole shift: never below
// it, on a few-term and a dense polynomial, for offsets inside and beyond the geometric case; and within a factor
// of four of it where one term carries the tail (the Mignotte polynomial's x^n beside its pair, the case of
// mignotte_generalized_512_16_16, where the source's value at the midpoint was rounding noise).
static void expansion_tail_bound() {
    const slong n = 257, working = 4096, length = 3;
    d::Integer a, coefficient;
    fmpz_one(a.value());
    fmpz_mul_2exp(a.value(), a.value(), 14);
    fmpz_poly_t mignotte, dense;
    fmpz_poly_init(mignotte);
    fmpz_poly_init(dense);
    fmpz_poly_set_coeff_si(mignotte, n, 1);
    fmpz_poly_set_coeff_si(mignotte, 0, -2);
    fmpz_mul_si(coefficient.value(), a.value(), 4);
    fmpz_poly_set_coeff_fmpz(mignotte, 1, coefficient.value());
    fmpz_mul(coefficient.value(), a.value(), a.value());
    fmpz_mul_si(coefficient.value(), coefficient.value(), -2);
    fmpz_poly_set_coeff_fmpz(mignotte, 2, coefficient.value());
    for (slong i = 0; i <= 40; ++i)
        fmpz_poly_set_coeff_si(dense, i, (i * 7919) % 23 - 11 + (i == 40 ? 12 : 0));
    for (const auto *integer : {mignotte, dense}) {
        const slong degree = fmpz_poly_degree(integer);
        d::PolynomialSource polynomial(
            [integer](acb_poly_t out, slong precision) { acb_poly_set_fmpz_poly(out, integer, precision); }, degree,
            true);
        for (const auto &[re, im] : {std::pair{0x1p-14, 0.0}, std::pair{0x1p-14 + 0x1p-40, 0x1p-45}, std::pair{0.75, 0.75}})
            for (const auto step : {0x1p-40, 0x1p-20, 0.5}) {
                d::ComplexBall center, offset, power;
                acb_set_d_d(center.value(), re, im);
                acb_set_d_d(offset.value(), step, step / 3);
                d::ComplexPolynomial whole;
                d::Magnitude bound, lower, upper, term, size;
                require(polynomial.taylor(whole.value(), center.value(), working) &&
                            polynomial.taylor_tail_bound(bound.value(), center.value(), offset.value(), working, length),
                        "no tail bound about the centre");
                for (slong j = length; j < acb_poly_length(whole.value()); ++j) {
                    acb_pow_ui(power.value(), offset.value(), static_cast<ulong>(j), working);
                    acb_mul(power.value(), power.value(), whole.value()->coeffs + j, working);
                    acb_get_mag_lower(term.value(), power.value());
                    mag_add_lower(lower.value(), lower.value(), term.value());
                    acb_get_mag(term.value(), power.value());
                    mag_add(upper.value(), upper.value(), term.value());
                }
                require(mag_cmp(bound.value(), lower.value()) >= 0, "the tail bound is below the dropped terms");
                if (integer == mignotte && re < 0.5 && step < 0x1p-30) {
                    mag_mul_2exp_si(size.value(), upper.value(), 2);
                    require(mag_cmp(bound.value(), size.value()) <= 0,
                            "the tail bound is loose where one term carries the tail");
                }
            }
    }
    fmpz_poly_clear(mignotte);
    fmpz_poly_clear(dense);
}

// The real engine takes the half with fewer sign variations first, so few halves wait beside a cluster's
// descent. x^31 - 2 (2^20 x - 1)^2 has a pair about 2^-329 apart near 2^-20, in the left halves of its root
// bound's interval, and a third root beyond it; taken left first, every right half waited for the pair, one a
// bisection (mignotte-krs-129-t2048 as a real-only query: 2,192 halves after 2,200 steps, then 4 GB and a kill).
// Taken in Descartes' order, at most a bisection's two halves are pending beside the components.
static void real_isolation_keeps_few_halves_waiting() {
    d::Integer a, coefficient;
    fmpz_one(a.value());
    fmpz_mul_2exp(a.value(), a.value(), 20);
    d::IntegerPolynomial polynomial;
    fmpz_poly_set_coeff_si(polynomial.value(), 31, 1);
    fmpz_poly_set_coeff_si(polynomial.value(), 0, -2);
    fmpz_mul_si(coefficient.value(), a.value(), 4);
    fmpz_poly_set_coeff_fmpz(polynomial.value(), 1, coefficient.value());
    fmpz_mul(coefficient.value(), a.value(), a.value());
    fmpz_mul_si(coefficient.value(), coefficient.value(), -2);
    fmpz_poly_set_coeff_fmpz(polynomial.value(), 2, coefficient.value());
    certpolysolve::detail::SearchStatistics statistics{};
    certpolysolve::detail::CandidateStatistics extended{};
    d::RealRefinement engine(polynomial.value(), 1, 64, statistics, extended);
    cps_options options;
    constexpr slong roots = 3, bisection = 2;
    constexpr unsigned steps = 10000;
    uint64_t most = 0;
    std::vector<d::Root> cover;
    unsigned taken = 0;
    for (; taken < steps && !engine.cover(cover); ++taken) {
        require(engine.step(options) == CPS_STOP_NONE, "the real isolation stopped");
        most = std::max(most, engine.cells());
    }
    require(taken < steps && engine.count() == roots, "the real isolation did not finish");
    require(most <= roots + bisection, "halves waited beside the pair's descent");
}

int main() {
    certpolysolve::ThreadScope lifetime;
    try {
        exact_count_near_node();
        local_model_solves_without_exact_preparation();
        local_model_budget_options();
        local_capability_cannot_publish_source_roots();
        geometric_restart_ends_a_crawl();
        restart_hints_stay_apart();
        restart_centre_is_the_roots_mean();
        restart_within_the_goal_takes_no_finer_model();
        a_misplaced_node_keeps_the_other_clusters();
        restart_precisions_follow_the_hints();
        last_active_cluster_is_a_group();
        origin_nodes_get_coordinate_bits();
        mixed_rounds_sweep_in_binary64();
        mixed_rounds_respect_the_iteration_budget();
        binary64_packets_keep_points_at_the_default_goal();
        exact_roots_do_not_drive_the_precision();
        exact_count_split();
        counted_scale_hints_respect_the_final_request();
        counted_restart_keeps_the_callers_resolution();
        local_seeds_do_not_escalate();
        restart_steps_take_no_expansion();
        restart_centres_stay_in_their_boxes();
        exact_models_take_the_whole_centre();
        restart_models_are_sized_for_the_callers_goal();
        restart_models_are_not_raised_twice();
        restart_models_keep_the_bits_a_raise_showed();
        unresolved_restart_components_start_at_their_sizing();
        unrecorded_restart_models_start_at_the_regenerations_precision();
        restart_hints_that_would_merge_are_placed_by_scale();
        source_iterated_members_ask_no_resolution();
        restart_models_resolve_to_the_iteration();
        restart_model_takes_the_expansion_head();
        head_models_measure_their_dropped_terms();
        restarted_members_keep_the_round_precision();
        undetermined_restart_hints_are_placed_by_scale();
        local_seeds_from_an_unresolved_model();
        regeneration_precision_memory();
        regeneration_takes_the_integer_kernels();
        regeneration_updates_unmoved_nodes();
        kept_nodes_stay_inactive_in_a_mixed_component();
    kept_node_not_reached_stays_active();
        cancelled_restart();
        unsupported_query_preserves_state();
        cover_ownership();
        coordinate_stopping();
        frozen_duplicates();
        duplicate_groups_match_pairwise();
        outward_bounds_never_separate_overlapping_rectangles();
        expansion_head_of_few_terms();
        expansion_tail_bound();
        real_isolation_keeps_few_halves_waiting();
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
