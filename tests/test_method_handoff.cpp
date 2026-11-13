#include "solve/policies.hpp"
#include "numeric/arithmetic/polynomial_arithmetic.hpp"
#include "exact_winding_api.hpp"
#include "api_test_support.hpp"
#include "internal.hpp"
#include "core/dovetail.hpp"
#include "methods/subdivision/subdivision.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <random>

namespace d = certpolysolve::detail;
using namespace reference;

int main() {
    certpolysolve::ThreadScope lifetime;
    try {
        // The dovetail share schedule: a fixed share, each Subdivision step
        // paid for before it runs (the last step's cost predicts the next),
        // every step charged at least one operation.
        {
            using S = d::DovetailSchedule;
            // A 1/16 share, engaged after 100 operations of iteration work,
            // with a first step expected to cost 10.
            S schedule(4, 100, 10);
            require(schedule.share() == 1.0 / 16 && !schedule.subdivision_due(),
                    "Subdivision was due before its share covered the first step");
            schedule.charge_iteration(149); // 10 of 159 is above 1/16
            require(!schedule.subdivision_due(), "Subdivision was due before its share covered the first step");
            schedule.charge_iteration(1); // 10 of 160
            require(schedule.subdivision_due(), "the share covered the first step, but it was not due");
            schedule.charge_subdivision(40); // four times the prediction
            // The next step is predicted at four times the last, 160:
            // 40 + 160 of 3,199 is above 1/16.
            schedule.charge_iteration(2849);
            require(!schedule.subdivision_due(), "a step that cost more than predicted was not paid off first");
            schedule.charge_iteration(1); // 200 of 3,200
            require(schedule.subdivision_due(), "a paid-off share did not make the next step due");
            schedule.charge_subdivision(0);
            require(schedule.subdivision() == 41, "a step that computes nothing was not charged");
            // Four iteration steps so far, 3,000 operations: a turn of
            // Subdivision is sized to their mean.
            require(schedule.iteration_step() == 750, "a turn was not sized to the mean iteration step");
            require(S(4, 100).iteration_step() == 100,
                    "before the first iteration step, a turn was not sized to the earlier work");
            // Records. Engagement counts as the first, with 100 operations
            // before it; a gap is overdue once it exceeds all earlier work.
            S records(4, 100);
            records.charge_iteration(100);
            require(!records.overdue(), "a gap equal to all earlier work was overdue");
            records.charge_iteration(1);
            require(records.overdue(), "a gap above all earlier work was not overdue");
            records.observe(S::Structure{4, 2, -3}, 0); // a record, 201 operations in
            require(!records.overdue(), "a record did not clear the overdue state");
            records.charge_iteration(201);
            require(!records.overdue(), "a gap equal to all earlier work was overdue");
            records.charge_iteration(1);
            require(records.overdue(), "a gap above all earlier work was not overdue");
            records.observe(S::Structure{4, 2, -3}, -0.5);
            require(records.overdue(), "a partial radius decrease counted as a record");
            records.observe(S::Structure{4, 2, -3}, -1);
            require(!records.overdue(), "a halved radius did not count as a record");
            records.charge_iteration(1000);
            require(records.overdue(), "a long gap was not overdue");
            // The iteration's own convergence: a halved largest correction
            // is a record; a smaller contraction is not.
            records.charge_iteration(5000);
            require(records.overdue(), "a long gap was not overdue");
            records.observe(S::Structure{4, 2, -3}, -1, -10);
            require(!records.overdue(), "a first correction did not count as a record");
            records.charge_iteration(20000);
            records.observe(S::Structure{4, 2, -3}, -1, -10.5);
            require(records.overdue(), "a partial contraction counted as a record");
            records.observe(S::Structure{4, 2, -3}, -1, -11);
            require(!records.overdue(), "a halved correction did not count as a record");
            // At the precision budget the share stops being the junior one:
            // Subdivision runs while its work, the next step included, is
            // within the iteration's (a share of one half), from whatever
            // share it had.
            S budget(4, 100, 10);
            budget.charge_iteration(100); // 10 of 110: not due at 1/16
            require(!budget.subdivision_due(), "Subdivision was due at a sixteenth before the budget");
            budget.at_budget(true);
            require(budget.at_budget_share() && budget.share() == 0.5 && budget.subdivision_due(),
                    "the budget did not make the two engines equal partners");
            budget.charge_subdivision(60); // 60 + 240 predicted of 400: not within the iteration's 100
            require(!budget.subdivision_due(), "an equal share let Subdivision run ahead of the iteration");
            budget.charge_iteration(200); // 300 >= 60 + 240
            require(budget.subdivision_due(), "an equal share did not let Subdivision catch up");
            budget.at_budget(false); // the budget raised: 300 of 600 is above 1/16
            require(!budget.at_budget_share() && budget.share() == 1.0 / 16 && !budget.subdivision_due(),
                    "a raised budget did not restore the share");
            // A step stopped before an operation of known cost (a disk test
            // under the share's allowance): available() is the largest step
            // the share covers, and after wait_for(c) the next step is due
            // exactly when the share covers c.
            S waiting(4, 100, 10);
            waiting.charge_iteration(1500); // a sixteenth of 1,600 is 100 above Subdivision's 0
            require(waiting.available() == 100, "the share's allowance was not the largest covered step");
            waiting.wait_for(150);
            require(!waiting.subdivision_due(), "a step waiting for 150 was due with 100 covered");
            waiting.charge_iteration(750); // 2,250 + 150 = 2,400, of which 150 is a sixteenth
            require(waiting.subdivision_due() && waiting.available() == 150,
                    "a step waiting for 150 was not due once the share covered it");
            S half(1, 100, 10);
            half.at_budget(true);
            require(half.share() == 0.5 && !half.at_budget_share(),
                    "the budget changed a share that was already one half");
        }
        // The meter's charge for Subdivision's steps follows FLINT's block
        // multiplication (dovetail.hpp): a flat exponent profile is one block
        // and one quasi-linear product; a profile whose exponents vary more
        // than 3 prec + 512 bits between neighbours is cut into blocks, and
        // every block pair is charged; a Graeffe step's squared polynomial
        // has twice the profile's curvature, hence more blocks.
        {
            const slong n = 769, prec = 2720;
            d::ComplexPolynomial flat, steep;
            acb_poly_fit_length(flat.value(), n + 1);
            acb_poly_fit_length(steep.value(), n + 1);
            for (slong j = 0; j <= n; ++j) {
                acb_set_ui(flat.value()->coeffs + j, 3); // one block
                // Exponent 2 j^2 bits: the curvature exceeds the height
                // allowance every few coefficients.
                acb_set_ui(steep.value()->coeffs + j, 1);
                acb_mul_2exp_si(steep.value()->coeffs + j, steep.value()->coeffs + j, 2 * j * j);
            }
            _acb_poly_set_length(flat.value(), n + 1);
            _acb_poly_set_length(steep.value(), n + 1);
            const double length = n + 1;
            const double one_block = 3 * (2 * length) * std::log2(4 * length);
            require(d::flint_block_lengths(flat.value(), prec).size() == 1 &&
                        d::block_product_operations(flat.value(), prec) == one_block,
                    "a flat profile is not one block charged its FFT");
            require(d::flint_block_lengths(steep.value(), prec).size() > 10 &&
                        d::block_product_operations(steep.value(), prec) > 4 * one_block &&
                        d::block_product_operations(steep.value(), prec) <= length * length,
                    "a curved profile is not charged for its blocks, within the schoolbook count");
            // The Taylor shift's branches: a convolution's product and n
            // scalings above twice the length, the divide-and-conquer
            // recursion below it, Horner at short lengths.
            require(d::block_taylor_shift_operations(flat.value(), prec) == one_block + length &&
                        d::block_taylor_shift_operations(flat.value(), 170) == std::log2(length) / 2 * one_block,
                    "a Taylor shift is not charged by FLINT's algorithm selection");
            d::ComplexPolynomial shortp;
            acb_poly_set_coeff_si(shortp.value(), 20, 1);
            require(d::block_taylor_shift_operations(shortp.value(), prec) == 21 * 21,
                    "a short Taylor shift is not charged as Horner's");
            // A degree-769 flat shift at 2,720 bits: far below n^2 (measured 18x).
            require(length * length / d::block_taylor_shift_operations(flat.value(), prec) > 8,
                    "the block charge is not far below n^2 for a flat profile");
        }
        // At secular's precision stop Subdivision continues from whichever
        // cover is ahead: by structure first, then by the largest region still
        // unresolved. (By structure alone, a dovetailed search whose
        // unresolved cluster was still wider than secular's was kept, and the
        // precision stop took 106 work units in place of 42: fuzz seed
        // 7664058671117322999 at a 160-bit cap.)
        {
            // Regions: `resolved` single roots of radius 2^-40, and one
            // cluster of `count` roots of radius 2^radius_log2.
            const auto cover = [](int resolved, uint64_t count, slong radius_log2) {
                std::vector<d::Root> roots(static_cast<size_t>(resolved) + 1);
                for (size_t i = 0; i < roots.size(); ++i) {
                    const bool cluster = i + 1 == roots.size();
                    acb_set_si(roots[i].box.value(), static_cast<slong>(i));
                    mag_set_ui_2exp_si(arb_radref(acb_realref(roots[i].box.value())), 1, cluster ? radius_log2 : -40);
                    roots[i].count = cluster ? count : 1;
                    roots[i].multiplicity = cluster ? 0 : 1;
                }
                return roots;
            };
            require(d::cover_ahead(cover(9, 7, -60), cover(9, 7, -20)) && !d::cover_ahead(cover(9, 7, -20), cover(9, 7, -60)),
                    "with equal structure, the cover with the smaller unresolved region was not ahead");
            require(!d::cover_ahead(cover(9, 7, -20), cover(9, 7, -20)), "a cover was ahead of itself");
            require(d::cover_ahead(cover(12, 4, -10), cover(9, 7, -60)) && !d::cover_ahead(cover(9, 7, -60), cover(12, 4, -10)),
                    "a better structure was not ahead of a smaller radius");
        }
        // The regeneration meter. A Horner step at a node with a short
        // mantissa costs far less than the same precision's full
        // multiply-add; the table's corners are reproduced; a mantissa longer
        // than the precision costs what the full one does; binary64 is 1.
        {
            require(std::abs(d::horner_weight(65536, 53) - 4701) < 1 &&
                        std::abs(d::horner_weight(65536, 65536) - 203445) < 1 &&
                        std::abs(d::horner_weight(64, 53) - 107) < 1,
                    "the Horner weight does not reproduce its calibration");
            require(d::horner_weight(30000, 53) < d::horner_weight(30000, 30000) / 16,
                    "a binary64 node at 30,000 bits was charged like a full-length one");
            require(d::horner_weight(4096, 128) < d::horner_weight(4096, 1024) &&
                        d::horner_weight(4096, 1024) < d::horner_weight(16384, 1024),
                    "the Horner weight is not increasing in the precision and the mantissa");
            require(d::horner_weight(256, 4096) == d::horner_weight(256, 256) && d::horner_weight(53, 53) == 1,
                    "the Horner weight's limits are wrong");
        }
        // The meter's other interpolations, in base 2 by portable::log2 and portable::exp2 since 21 October 2026
        // (the same bits on every platform): the calibration's points are reproduced to rounding, and the operation
        // weight rises with the precision.
        {
            constexpr std::array<std::array<double, 2>, 9> limbs_weights{{{1, 148.0},
                                                                         {2, 188.2},
                                                                         {4, 230.6},
                                                                         {8, 312.4},
                                                                         {16, 684.3},
                                                                         {32, 2059.3},
                                                                         {64, 4976.1},
                                                                         {256, 46535.6},
                                                                         {1024, 227839.8}}};
            constexpr double rounding = 0x1p-40;
            for (const auto &[limbs, weight] : limbs_weights)
                require(std::abs(d::operation_weight(static_cast<int64_t>(64 * limbs)) / weight - 1) < rounding,
                        "the operation weight does not reproduce its calibration");
            require(d::operation_weight(53) == 1 && d::operation_weight(64) == 148,
                    "the operation weight's first points are not exact");
            double previous = d::operation_weight(53);
            for (int64_t bits = 54; bits <= (int64_t{1} << 20); bits += 61) {
                const double weight = d::operation_weight(bits);
                require(weight >= previous, "the operation weight is not increasing in the precision");
                previous = weight;
            }
            constexpr std::array<std::array<double, 2>, 8> shift_passes{
                {{256, 28.4}, {512, 25.8}, {1024, 20.4}, {2048, 13.6}, {4096, 11.3}, {8192, 8.6}, {16384, 5.3}, {32768, 3.8}}};
            for (const auto &[bits, passes] : shift_passes)
                require(std::abs(d::heuristics::convolution_shift_passes(static_cast<slong>(bits)) / passes - 1) < rounding,
                        "the convolution shift's passes do not reproduce their measurement");
        }
        // Explicit secular first tries classical's machine/Newton certificate
        // (docs/method-composition-design-20260928.md). On a sparse polynomial
        // that certifies everything at 85 working bits; secular's own route
        // needs 170 here, so the work-bit cap makes the regression deterministic.
        {
            constexpr int degree = 256, terms = 10;
            std::mt19937_64 random(20260928);
            ComplexPolynomial sparse(degree + 1, Complex(0));
            const auto coefficient = [&] {
                return Rational(static_cast<long long>(random() % (2u << 20)) - (1 << 20)) / power2(20);
            };
            sparse[0] = {coefficient(), coefficient()};
            for (int k = 0; k < terms; ++k)
                sparse[1 + random() % (degree - 1)] = {coefficient(), coefficient()};
            sparse[degree] = {Rational(1)};
            certpolysolve::Options options;
            options.method = certpolysolve::Method::secular;
            options.goal = certpolysolve::Goal::approximation();
            options.bits = 53;
            options.max_bits = 85;
            const auto result = certpolysolve::solve(reference::input(sparse), options);
            require(result.complete() && result.goal_met(), "secular did not use the machine/Newton certificate");
            const auto regions = result.roots();
            check_disjoint(regions);
            uint64_t total = 0;
            for (const auto &region : regions) total += region.count;
            require(total == degree, "machine/Newton certificate miscounted the sparse roots");
        }
        // MPSolve's lsr1 rotated by i under the maintained secular recipe
        // (research/lsr1-rotate-i-20260928): x^498 (x - 20000i)^2 dominates
        // the top terms, a double root at any practical precision, and a
        // near-triple cluster sits at 1e-4 i. After the first growth the pair
        // converges linearly into its double root while nothing else can be
        // certified. (Until 7 October 2026 this also compared the private
        // restart_motion_growth rule, since removed, against the recipe.)
        {
            constexpr int degree = 500;
            ComplexPolynomial lsr1(degree + 1, Complex(0));
            const std::pair<int, const char *> terms[] = {
                {0, "1"}, {1, "30000"}, {2, "300000000"}, {3, "1000000000000"}, {200, "1"},
                {298, "400000000"}, {299, "12000000040000"}, {300, "120000001200000001"},
                {301, "400000012000000030000"}, {302, "40000000300000000"}, {303, "1000000000000"},
                {498, "400000000"}, {499, "40000"}, {500, "1"}};
            // p(i x): the coefficient of x^k gains the factor i^k.
            for (const auto &[k, value] : terms) {
                const Rational a(value);
                const Rational unit[4][2] = {{1, 0}, {0, 1}, {-1, 0}, {0, -1}};
                lsr1[static_cast<size_t>(k)] = {a * unit[k % 4][0], a * unit[k % 4][1]};
            }
            const auto input = reference::input(lsr1);
            const auto policy = d::maintained_method_policy(CPS_METHOD_SECULAR);
            cps_options options;
            options.goal = CPS_ROOT_APPROXIMATION;
            options.bits = 53;
            d::SolveState solver(certpolysolve::detail::dense(*input.get()), policy);
            const auto result = solver.advance(&options);
            require(result.goal_met && result.complete && d::disjoint(result.roots),
                    "lsr1 rotated by i did not complete");
            uint64_t total = 0;
            bool pair = false, triple = false;
            for (const auto &root : result.roots) {
                auto endpoint = [&](size_t i) { return reference::endpoint(certpolysolve::Dyadic{
                    root.endpoints[i].mantissa_text(), root.endpoints[i].exponent_text()}); };
                auto contains = [&](const Complex &z) {
                    return endpoint(0) <= z.real && z.real <= endpoint(1) &&
                           endpoint(2) <= z.imag && z.imag <= endpoint(3);
                };
                total += root.count;
                pair = pair || (root.count == 2 && contains({0, 20000}));
                triple = triple || (root.count == 3 && contains({0, Rational(1, 10000)}));
            }
            require(total == degree && pair && triple, "lsr1 rotated by i lost its count-2 or count-3 region");
        }
        // Secular has no global convergence proof, so a factor whose precision
        // budget it has spent continues with the complete method from its
        // certified cover (principled audit, fix 2). x^8 - 2 (64 x - 1)^2 has
        // two roots sqrt(2) 2^-30 apart near 1/64. At an 82-bit budget the
        // plain evaluation of p at the pair's nodes cancels past it, so
        // secular stops on the budget; Subdivision finishes from its cover
        // within the same budget. (The research build's control
        // exact_residues_at_budget, off by default, would let the recipe
        // finish this input itself at 80 bits.) The
        // window is narrow: from 84 bits secular finishes alone, and below 82
        // Subdivision's last box about one of the pair can land at radius
        // exactly 2^-59, which misses the goal's relative radius at 1/64 by
        // 2^-112, with the budget allowing no further bisection. Where in the
        // window that happens depends on the cover handed over: 80 bits
        // sufficed while the certificate attempt's certificates for the six
        // other roots were thrown away (the pair refused it), not since they
        // are kept.
        {
            constexpr int degree = 8;
            ComplexPolynomial mignotte(degree + 1, Complex(0));
            const Rational a = power2(6);
            mignotte[0] = {Rational(-2)};
            mignotte[1] = {Rational(4) * a};
            mignotte[2] = {Rational(-2) * a * a};
            mignotte[degree] = {Rational(1)};
            const auto input = reference::input(mignotte);
            auto policy = d::maintained_method_policy(CPS_METHOD_SECULAR);
            // The limb evaluation is held off so that the plain evaluation cancels: at 82 bits it works in a
            // window of three limbs with these short exact coefficients, which holds the cancellation, and
            // secular finishes alone.
            policy.regeneration.limb_horner = false;
            // The compensated tier proves the pair at the 53-bit goal from one binary64 evaluation, where the ball
            // certificate's 82-bit evaluation cancels, and secular then finishes alone: the window is the ball
            // certificate's, so the tier is held off. The handoff still serves wherever the tier declines.
            policy.compensated_certificates = false;
            // The handoff alone: dovetailed, Subdivision runs as secular's
            // equal partner at the budget (dovetail.hpp, at_budget) and
            // finishes this small factor before secular stops.
            policy.dovetail_subdivision = false;
            unsigned handoffs = 0;
            policy.trace_data = &handoffs;
            policy.trace = [](void *p, const d::TraceEvent &event) {
                if (std::string_view(event.stage) == "secular_subdivision_handoff")
                    ++*static_cast<unsigned *>(p);
            };
            cps_options options;
            options.goal = CPS_ROOT_APPROXIMATION;
            options.bits = 53;
            options.max_work_bits = 82;
            d::SolveState solver(certpolysolve::detail::dense(*input.get()), policy);
            auto result = solver.advance(&options);
            require(handoffs == 1, "secular did not hand its spent precision budget to subdivision");
            require(result.complete && result.goal_met && result.termination == CPS_SUCCESS &&
                        d::disjoint(result.roots),
                    "subdivision did not finish secular's factor within the same budget");
            const auto published = certpolysolve_test::published(std::move(result));
            uint64_t total = 0;
            for (const auto &region : published.roots()) {
                const auto exact = rectangle_count(mignotte, region);
                require(exact && static_cast<uint64_t>(*exact) == region.count,
                        "subdivision's count disagrees with the exact winding count");
                total += region.count;
            }
            require(total == degree, "subdivision lost roots after the secular handoff");
            // With the tier on, the same budget suffices: secular proves the pair alone, and every region's count
            // is still the exact winding count.
            auto tiered = d::maintained_method_policy(CPS_METHOD_SECULAR);
            tiered.regeneration.limb_horner = false;
            tiered.dovetail_subdivision = false;
            unsigned tiered_handoffs = 0;
            tiered.trace_data = &tiered_handoffs;
            tiered.trace = policy.trace;
            d::SolveState alone(certpolysolve::detail::dense(*input.get()), tiered);
            auto finished = alone.advance(&options);
            require(tiered_handoffs == 0 && finished.complete && finished.goal_met &&
                        finished.termination == CPS_SUCCESS && d::disjoint(finished.roots),
                    "with the compensated tier secular did not finish the pair alone at the same budget");
            const auto proved = certpolysolve_test::published(std::move(finished));
            uint64_t counted = 0;
            for (const auto &region : proved.roots()) {
                const auto exact = rectangle_count(mignotte, region);
                require(exact && static_cast<uint64_t>(*exact) == region.count,
                        "with the compensated tier a count disagrees with the exact winding count");
                counted += region.count;
            }
            require(counted == degree, "with the compensated tier secular lost roots");
            // At a 72-bit budget no method separates the pair. A split whose
            // direct counts that budget cannot decide is abandoned, so the
            // solve stops promptly for precision with its complete certified
            // cover, instead of multiplying undecided cells until the cell
            // budget.
            handoffs = 0;
            options.max_work_bits = 72;
            d::SolveState limited(certpolysolve::detail::dense(*input.get()), policy);
            auto stopped = limited.advance(&options);
            require(handoffs == 1 && stopped.termination == CPS_PRECISION_LIMIT &&
                        stopped.progress.stop_reason == CPS_STOP_PRECISION && stopped.complete &&
                        !stopped.goal_met && d::disjoint(stopped.roots),
                    "subdivision did not stop for the precision budget");
            const auto partial = certpolysolve_test::published(std::move(stopped));
            total = 0;
            for (const auto &region : partial.roots()) {
                const auto exact = rectangle_count(mignotte, region);
                require(exact && static_cast<uint64_t>(*exact) == region.count,
                        "a precision-limited count disagrees with the exact winding count");
                total += region.count;
            }
            require(total == degree, "the precision-limited cover lost roots");
        }
        // At its precision budget regeneration accepts residue enclosures with
        // no correct bits. Such a node is a root at that precision, so the
        // model settles it at once (residue_unresolved); before, the model kept
        // proposing uninformed steps until the work budget ran out, and the
        // spent precision budget never reached subdivision. With the plain
        // evaluation alone (the research build's exact_residues_at_budget is
        // off by default, so the plain evaluation's cancellation leaves the
        // model uninformed), Wilkinson's
        // polynomial at 64 bits stalls secular in binary64 packets and
        // subdivision then separates its roots within the same budget; a
        // 2^-89 cluster at 256 bits stalls it in multiprecision packets, and
        // neither engine separates the cluster from that evaluation, so the
        // solve stops for precision with a complete cover. (With the research
        // build's exact residues, secular finishes both itself: the block after this one.)
        {
            struct Case {
                const char *name;
                std::vector<Complex> roots;
                cps_goal goal;
                uint64_t bits, budget;
            };
            std::vector<Complex> wilkinson, cluster{Complex(Rational(-1, 63))};
            for (int j = 1; j <= 20; ++j) wilkinson.emplace_back(Rational(j));
            for (int j : {0, 1, 2, 5}) cluster.emplace_back(Rational(-2, 3) + j * power2(-89));
            for (const auto &c : {Case{"Wilkinson-20 at 64 bits", wilkinson, CPS_ROOT_APPROXIMATION, 53, 64},
                                  Case{"a 2^-89 cluster at 256 bits", cluster, CPS_ROOT_ACCURACY, 16, 256}}) {
                const auto fail = [&](const char *what) { return std::string(c.name) + ": " + what; };
                const auto polynomial = from_roots(c.roots);
                const auto input = reference::input(polynomial);
                auto policy = d::maintained_method_policy(CPS_METHOD_SECULAR);
                // The handoff with no dovetailed search to keep: a kept search
                // has precisions of its own (dovetailing is tested below).
                policy.dovetail_subdivision = false;
                struct Progress { unsigned handoffs = 0; bool cancel = false; } progress;
                policy.trace_data = &progress;
                policy.trace = [](void *p, const d::TraceEvent &event) {
                    if (std::string_view(event.stage) == "secular_subdivision_handoff") {
                        auto &state = *static_cast<Progress *>(p);
                        ++state.handoffs;
                        state.cancel = true;
                    }
                };
                cps_options options;
                options.goal = c.goal;
                options.bits = c.bits;
                options.max_work_bits = c.budget;
                // Far more than either solve needs; the uninformed steps
                // exhausted it.
                options.work_units = 2000;
                options.cancel_data = &progress;
                options.cancel = [](void *p) -> int { return static_cast<Progress *>(p)->cancel; };
                d::SolveState solver(certpolysolve::detail::dense(*input.get()), policy);
                // Subdivision starts at the goal's working precision, not at
                // the precision secular stopped at.
                const auto paused = solver.advance(&options);
                const uint64_t start = std::min<uint64_t>(c.budget, std::max<uint64_t>(64, c.bits + 32));
                require(progress.handoffs == 1 && paused.termination == CPS_CANCELLED &&
                            paused.progress.min_pending_bits == start && paused.progress.max_pending_bits == start,
                        fail("subdivision did not start at the goal's working precision after the handoff").c_str());
                // The whole solve, without the pause: a fresh solve (one advance a solve state since 18 October
                // 2026; until then the paused state was resumed).
                options.cancel = nullptr;
                progress = Progress{};
                d::SolveState whole(certpolysolve::detail::dense(*input.get()), policy);
                auto result = whole.advance(&options);
                const bool separable = c.budget == 64;
                require(progress.handoffs == 1 && result.complete && d::disjoint(result.roots),
                        fail("secular did not hand its precision budget to subdivision at an uninformed model").c_str());
                require(separable ? result.goal_met && result.termination == CPS_SUCCESS
                                  : !result.goal_met && result.termination == CPS_PRECISION_LIMIT &&
                                        result.progress.stop_reason == CPS_STOP_PRECISION,
                        fail("an uninformed model did not end in the budget's result").c_str());
                const auto published = certpolysolve_test::published(std::move(result));
                uint64_t total = 0;
                for (const auto &region : published.roots()) {
                    const auto exact = rectangle_count(polynomial, region);
                    require(exact && static_cast<uint64_t>(*exact) == region.count,
                        fail("a count after an uninformed model disagrees with the exact winding count").c_str());
                    total += region.count;
                }
                require(total == c.roots.size(),
                        fail("the cover after an uninformed model lost roots").c_str());
            }
        }
        // The research build's control exact_residues_at_budget
        // (engine_policy.hpp): at the budget a residue the plain evaluation
        // leaves without its demanded accuracy is taken from the factor's
        // exact value at the node, rounded once. With the binary64 machine
        // retired at the budget and the restart circle's placement (both
        // always on), each case stops for precision or hands off with the
        // control off (every build: the maintained stop and hand-off) and,
        // in the research build, is finished by secular alone with it on, no
        // handoff.
        {
            struct Case {
                const char *name;
                std::vector<Complex> roots;
                cps_goal goal;
                uint64_t bits, budget;
                const char *cluster_radius; // for the cluster goal
                bool without_model; // stops for precision with the exact residues off
                uint64_t units; // at most, with them: a crawl guard, far above the measured 8-23
            };
            std::vector<Complex> wilkinson, cluster{Complex(Rational(-1, 63))};
            for (int j = 1; j <= 20; ++j) wilkinson.emplace_back(Rational(j));
            for (int j : {0, 1, 2, 5}) cluster.emplace_back(Rational(-2, 3) + j * power2(-89));
            // x^8 - 2 (64 x - 1)^2: the pair sqrt(2) 2^-30 apart near 1/64.
            ComplexPolynomial mignotte(9, Complex(0));
            mignotte[0] = {Rational(-2)};
            mignotte[1] = {Rational(256)};
            mignotte[2] = {Rational(-8192)};
            mignotte[8] = {Rational(1)};
            // Fuzz seed 12915955158621147063 under the cluster goal: five roots
            // within 2^-33 of -0.068 (a nested pair 2^-49 apart), a double
            // root and x^5 beside them, radius 2^-29 at a 96-bit budget. The
            // plain residues of the cluster have no correct bit at 96 bits.
            const std::vector<Complex> nested{
                Complex(Rational("-142955/2097152")), Complex(Rational("-142955/2097152")),
                Complex(Rational("-2342174723/34359738368")), Complex(Rational("-585543681/8589934592")),
                Complex(Rational("-2342174719/34359738368")), Complex(Rational("-1171087361/17179869184")),
                Complex(Rational("-38374190645247/562949953421312")),
                Complex(0), Complex(0), Complex(0), Complex(0), Complex(0)};
            // Fuzz seed 538619078387263787 (cap2x): four roots 2^-117 apart
            // near -2^-25, a 16-bit accuracy goal at a 128-bit budget. With
            // exact residues the model is resolved; the restart that follows
            // the binary64 floor placed the four nodes at the unit roots of
            // order four about the real centre, a conjugate-symmetric start
            // from which two nodes stayed on the bisector of two real roots
            // for 110 rounds (candidates.hpp, placement_angle).
            const std::vector<Complex> symmetric{
                Complex(Rational("-3089/137438953472")),
                Complex(Rational("-3734371856789589520667377665/166153499473114484112975882535043072")),
                Complex(Rational("-3734371856789589520667377663/166153499473114484112975882535043072")),
                Complex(Rational("-1867185928394794760333688831/83076749736557242056487941267521536"))};
            const auto cases = {
                Case{"Wilkinson-20 at 64 bits", wilkinson, CPS_ROOT_APPROXIMATION, 53, 64, nullptr, false, 100},
                Case{"a 2^-89 cluster at 256 bits", cluster, CPS_ROOT_ACCURACY, 16, 256, nullptr, true, 100},
                Case{"the pair near 1/64 at 80 bits", {}, CPS_ROOT_APPROXIMATION, 53, 80, nullptr, false, 100},
                Case{"the nested cluster under the cluster goal at 96 bits", nested, CPS_CLUSTER_RADIUS, 53, 96, "1@-29", true, 100},
                Case{"four roots 2^-117 apart at 128 bits", symmetric, CPS_ROOT_ACCURACY, 16, 128, nullptr, true, 100}};
            for (const auto &c : cases) {
                const auto fail = [&](const char *what) { return std::string(c.name) + ": " + what; };
                const auto polynomial = c.roots.empty() ? mignotte : from_roots(c.roots);
                const auto input = reference::input(polynomial);
#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
                const std::initializer_list<bool> arms{false, true}; // the exact residues off, then on
#else
                const std::initializer_list<bool> arms{false}; // the default build has no exact residues
#endif
                for (const bool model : arms) {
                    auto policy = d::maintained_method_policy(CPS_METHOD_SECULAR);
#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
                    policy.exact_residues_at_budget = model;
#endif
                    // The cases are about what the exact residues add to the plain evaluation at the
                    // budget. The limb evaluation is held off: its window of whole limbs holds more than
                    // the budget's bits of short exact coefficients, and resolves the pair near 1/64 at 80
                    // bits by itself.
                    policy.regeneration.limb_horner = false;
                    // The compensated tier, likewise, proves the pair near 1/64 from one binary64 evaluation.
                    policy.compensated_certificates = false;
                    policy.dovetail_subdivision = false;
                    unsigned handoffs = 0;
                    policy.trace_data = &handoffs;
                    policy.trace = [](void *p, const d::TraceEvent &event) {
                        if (std::string_view(event.stage) == "secular_subdivision_handoff")
                            ++*static_cast<unsigned *>(p);
                    };
                    cps_options options;
                    options.goal = c.goal;
                    options.bits = c.bits;
                    options.max_work_bits = c.budget;
                    options.work_units = 2000;
                    if (c.cluster_radius)
                        require(cps_options_set_cluster_radius(&options, c.cluster_radius) == CPS_OK,
                                fail("the cluster radius was refused").c_str());
                    d::SolveState solver(certpolysolve::detail::dense(*input.get()), policy);
                    const auto result = solver.advance(&options);
                    require(result.complete && d::disjoint(result.roots),
                            fail("the solve did not end with a disjoint complete cover").c_str());
                    uint64_t total = 0;
                    for (const auto &root : result.roots) total += root.count;
                    require(total == static_cast<uint64_t>(polynomial.size() - 1), fail("the cover lost roots").c_str());
                    if (!model) {
                        if (c.without_model)
                            require(!result.goal_met && result.termination == CPS_PRECISION_LIMIT &&
                                        result.progress.stop_reason == CPS_STOP_PRECISION,
                                    fail("stopped otherwise than for precision without the exact residues; "
                                         "the case no longer shows what they add").c_str());
                        else
                            require(handoffs == 1,
                                    fail("finished without a handoff without the exact residues; "
                                         "the case no longer shows what they add").c_str());
                        continue;
                    }
                    require(result.goal_met && result.termination == CPS_SUCCESS,
                            fail("not finished with the exact residues at the budget").c_str());
                    require(handoffs == 0, fail("secular handed off a factor it finishes with the exact residues").c_str());
                    require(result.progress.work_units_this_call <= c.units,
                            fail("finished with the exact residues, but slowly").c_str());
                }
            }
        }
#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
        // The exact residues under a centred study: the source is the centred
        // factor and the nodes live in its coordinates, so the exact value
        // must be taken from that polynomial (a mismatch would be a false
        // certificate, not a failure). Wilkinson-64 is centred at 65/2
        // (centering-selected) and at a 128-bit budget its plain residues
        // have no correct bit, so every node takes the exact value; every
        // region is checked against the exact winding count.
        {
            std::vector<Complex> wilkinson;
            for (int j = 1; j <= 64; ++j) wilkinson.emplace_back(Rational(j));
            const auto polynomial = from_roots(wilkinson);
            const auto input = reference::input(polynomial);
            auto policy = d::maintained_method_policy(CPS_METHOD_SECULAR);
            policy.exact_residues_at_budget = true;
            struct Seen { unsigned centred = 0, exact = 0; } seen;
            policy.trace_data = &seen;
            policy.trace = [](void *p, const d::TraceEvent &event) {
                auto &s = *static_cast<Seen *>(p);
                if (std::string_view(event.stage) == "exact-centering" && std::string_view(event.reason) == "centering-selected")
                    ++s.centred;
                if (std::string_view(event.stage) == "exact_regeneration" && std::string_view(event.phase) == "end")
                    ++s.exact;
            };
            cps_options options;
            options.goal = CPS_ROOT_APPROXIMATION;
            options.bits = 53;
            options.max_work_bits = 128;
            d::SolveState solver(certpolysolve::detail::dense(*input.get()), policy);
            const auto result = solver.advance(&options);
            require(seen.centred == 1 && seen.exact > 0,
                    "Wilkinson-64 at 128 bits no longer exercises the exact residues under a centred study");
            require(result.goal_met && result.complete && result.termination == CPS_SUCCESS && d::disjoint(result.roots),
                    "Wilkinson-64 at 128 bits did not finish with exact residues under a centred study");
            const auto published = certpolysolve_test::published(std::move(result));
            uint64_t total = 0;
            for (const auto &region : published.roots()) {
                const auto exact = rectangle_count(polynomial, region);
                require(exact && static_cast<uint64_t>(*exact) == region.count && region.count == 1,
                        "a region of the centred exact-residue solve disagrees with the exact winding count");
                total += region.count;
            }
            require(total == 64, "the centred exact-residue solve lost roots");
        }
#endif
        // Subdivision at the precision budget retries a failed contraction at
        // half its exponent before it bisects (subdivision.cpp, the step's
        // Newton branch). Bisected at once, the children's direct tests are
        // undecidable at the budget, the split is abandoned and the search
        // stops for precision short of a resolution it can certify.
        {
            // The step that fails for geometry: two roots 2^-58 apart at
            // modulus 2^59 (one of them the integer -15 2^55), a 100-bit
            // approximation goal at 128 bits, started as the handoff starts
            // it, from one count-2 region (radius 2^-2, centred off both
            // roots). The contractions by 2, 4, 8 and 16 bits succeed, the one
            // by 32 bits (to 2^-62) cannot hold both roots, and the retry by
            // 16 bits (to 2^-46) meets the goal. (The factor of fuzz seed
            // 2687697434157657036 that holds the pair.)
            const std::vector<Complex> pair{
                Complex(Rational("-540431955284459520")),
                Complex(Rational("-155768905756044828855914889876602879/288230376151711744"))};
            const auto input = reference::input(from_roots(pair));
            const d::Encloser enclose = [&](acb_poly_t out, slong bits) {
                d::ComplexBall coefficient;
                acb_poly_zero(out);
                for (const auto &[i, c] : input.get()->terms) {
                    c.enclose(coefficient.value(), bits);
                    acb_poly_set_coeff_acb(out, static_cast<slong>(i), coefficient.value());
                }
            };
            d::Root region;
            // Centre -15 2^55 + 1/8: both roots are within 1/8 + 2^-58 of it.
            arb_set_si(acb_realref(region.box.value()), -540431955284459520L);
            arb_mul_2exp_si(acb_realref(region.box.value()), acb_realref(region.box.value()), 3);
            arb_add_ui(acb_realref(region.box.value()), acb_realref(region.box.value()), 1, 256);
            arb_mul_2exp_si(acb_realref(region.box.value()), acb_realref(region.box.value()), -3);
            mag_set_ui_2exp_si(arb_radref(acb_realref(region.box.value())), 1, -2);
            mag_set_ui_2exp_si(arb_radref(acb_imagref(region.box.value())), 1, -2);
            region.count = 2;
            cps_options options;
            options.goal = CPS_ROOT_APPROXIMATION;
            options.bits = 100;
            options.max_work_bits = 128;
            certpolysolve::detail::SearchStatistics statistics{};
            d::Subdivision search({region}, 1, options.max_work_bits, statistics);
            cps_stop_reason reason = CPS_STOP_NONE;
            // Far more steps than the solve needs (about ten); a guard only.
            for (unsigned step = 0; step < 1000 && reason == CPS_STOP_NONE && search.pending(); ++step)
                reason = search.step(options, enclose, 65536);
            std::vector<d::Root> cover;
            require(reason == CPS_STOP_NONE && !search.pending() && search.complete_cover(cover),
                    "Subdivision at the budget stopped where a gentler contraction reaches the goal");
            uint64_t total = 0;
            for (const auto &root : cover) {
                require(d::meets_goal(root, options, static_cast<slong>(options.max_work_bits)),
                        "a region of the retried contraction misses the goal");
                total += root.count;
            }
            require(total == 2, "the retried contraction lost a root");
            // The same search with every step under a small work allowance
            // (dovetailing's share, Subdivision::step): a disk test stops
            // before each Taylor shift or Graeffe step the allowance does not
            // cover, and the next step, given just the pending operation's
            // cost, resumes it. The cover is the unbounded search's.
            {
                certpolysolve::detail::SearchStatistics bounded_statistics{};
                d::Subdivision bounded({region}, 1, options.max_work_bits, bounded_statistics);
                cps_stop_reason bounded_reason = CPS_STOP_NONE;
                double allowance = 1;
                unsigned suspensions = 0;
                for (unsigned step = 0; step < 100000 && bounded_reason == CPS_STOP_NONE && bounded.pending(); ++step) {
                    bounded_reason = bounded.step(options, enclose, 65536, allowance);
                    if (bounded.suspended()) {
                        ++suspensions;
                        allowance = bounded.pending_work();
                    } else {
                        allowance = 1;
                    }
                }
                std::vector<d::Root> bounded_cover;
                require(bounded_reason == CPS_STOP_NONE && !bounded.pending() && bounded.complete_cover(bounded_cover),
                        "a search under small allowances did not finish");
                require(suspensions > 0, "small allowances never stopped a disk test");
                require(bounded_cover.size() == cover.size() &&
                            std::equal(cover.begin(), cover.end(), bounded_cover.begin(),
                                       [](const d::Root &a, const d::Root &b) {
                                           return a.count == b.count && acb_equal(a.box.value(), b.box.value());
                                       }),
                        "a search under small allowances found another cover");
            }
            // The step that fails for precision: 39 roots on the unit circle
            // (one pair repeated), 100-bit accuracy at 128 bits, through the
            // secular recipe. Secular's final cover leaves six regions at
            // radius 2^-99; from there the contraction to 2^-107 is below what
            // 128 bits resolve and the one to 2^-103 certifies. (Fuzz seed
            // 12447402984292181221.)
            const std::vector<std::pair<const char *, const char *>> circle{
                {"12/37", "35/37"}, {"12/37", "-35/37"}, {"60/61", "11/61"}, {"60/61", "-11/61"}, {"-3/5", "4/5"},
                {"-3/5", "-4/5"}, {"-15/17", "-8/17"}, {"-15/17", "8/17"}, {"-33/65", "-56/65"}, {"-33/65", "56/65"},
                {"8/17", "15/17"}, {"8/17", "-15/17"}, {"15/17", "-8/17"}, {"15/17", "8/17"}, {"-55/73", "48/73"},
                {"-55/73", "-48/73"}, {"-45/53", "-28/53"}, {"-45/53", "28/53"}, {"-9/41", "40/41"},
                {"-9/41", "-40/41"}, {"45/53", "28/53"}, {"45/53", "-28/53"}, {"-60/61", "11/61"},
                {"-60/61", "-11/61"}, {"28/53", "-45/53"}, {"28/53", "45/53"}, {"-21/29", "-20/29"},
                {"-21/29", "20/29"}, {"-56/65", "33/65"}, {"-56/65", "-33/65"}, {"-20/29", "-21/29"},
                {"-20/29", "21/29"}, {"40/41", "9/41"}, {"40/41", "-9/41"}, {"-45/53", "-28/53"},
                {"-45/53", "28/53"}, {"-8/17", "15/17"}, {"-8/17", "-15/17"}, {"-1", "0"}};
            std::vector<Complex> roots;
            for (const auto &[re, im] : circle) roots.emplace_back(Rational(re), Rational(im));
            const auto circle_input = reference::input(from_roots(roots));
            options.goal = CPS_ROOT_ACCURACY;
            d::SolveState solver(certpolysolve::detail::dense(*circle_input.get()), d::maintained_method_policy(CPS_METHOD_SECULAR));
            const auto result = solver.advance(&options);
            require(result.goal_met && result.complete && result.termination == CPS_SUCCESS && d::disjoint(result.roots),
                    "roots on the unit circle at 128 bits stopped for precision one contraction short of the goal");
            total = 0;
            for (const auto &root : result.roots) total += root.count;
            require(total == roots.size(), "the solve on the unit circle lost roots");
        }
        // A finished factor and a factor blocked by the precision budget: the
        // driver's stall counter counts both, and the exit taken on a finished
        // factor's turn left the stop reason unset (fuzz F1, x^2 (x - r) at 72
        // bits). The blocked factor's budget is why the solve stops.
        {
            const std::vector<Complex> roots{Complex(Rational("179/39212225875655018250662308278270164992")),
                                             Complex(0), Complex(0)};
            const auto input = reference::input(from_roots(roots));
            cps_options options;
            options.goal = CPS_ROOT_ACCURACY;
            options.bits = 100;
            options.max_work_bits = 72;
            d::SolveState solver(certpolysolve::detail::dense(*input.get()), d::maintained_method_policy(CPS_METHOD_CLASSICAL));
            const auto result = solver.advance(&options);
            require(result.complete && !result.goal_met && result.termination == CPS_PRECISION_LIMIT &&
                        result.progress.stop_reason == CPS_STOP_PRECISION,
                    "a precision-limited solve with a finished factor did not report the precision budget");
        }
        // At the precision limit the model steers no node whose residue lacks
        // correct bits, but the cover's radius still depends on where the
        // members sit. Local seeds land at the numerical roots, bunched
        // inside the pseudo-zero set, where the certificate is worst; the
        // certificate re-seed places them to minimise it. Fuzz case
        // 17285531581314345458: i and an 8-root cluster near 1.052i (a triple
        // root among them) within 2^-103, at a 256-bit budget. Settling the
        // cluster without the re-seed stopped this solve on its 100-unit work
        // budget; the promoted build finished it by uninformed jitter.
        {
            std::vector<Complex> roots{Complex(0, Rational(1))};
            for (const char *imaginary :
                 {"821437588947892652169863677083571/780872769740589311321969174511616",
                  "821437588947892652169863677083417/780872769740589311321969174511616",
                  "205359397236973163042465919270989/195218192435147327830492293627904",
                  "821437588947892652169863677083879/780872769740589311321969174511616",
                  "821437588947892652169863677083725/780872769740589311321969174511616", "81/77",
                  "410718794473946326084931838541901/390436384870294655660984587255808",
                  "410718794473946326084931838541747/390436384870294655660984587255808",
                  "205359397236973163042465919270989/195218192435147327830492293627904",
                  "205359397236973163042465919270989/195218192435147327830492293627904"})
                roots.emplace_back(Rational(0), Rational(imaginary));
            const auto polynomial = from_roots(roots);
            const auto input = reference::input(polynomial);
            cps_options options;
            options.goal = CPS_ROOT_APPROXIMATION;
            options.bits = 53;
            options.max_work_bits = 256;
            options.work_units = 100;
            options.live_cells = 256;
            d::SolveState solver(certpolysolve::detail::dense(*input.get()), d::maintained_method_policy(CPS_METHOD_SECULAR));
            auto result = solver.advance(&options);
            require(result.goal_met && result.complete && result.termination == CPS_SUCCESS &&
                        d::disjoint(result.roots),
                    "the certificate re-seed did not recover a cluster at the precision limit");
            const auto published = certpolysolve_test::published(std::move(result));
            uint64_t total = 0;
            for (const auto &region : published.roots()) {
                const auto exact = rectangle_count(polynomial, region);
                require(exact && static_cast<uint64_t>(*exact) == region.count,
                        "a count after the certificate re-seed disagrees with the exact winding count");
                total += region.count;
            }
            require(total == roots.size(), "the cover after the certificate re-seed lost roots");
            // A solve under a 96-bit budget stops for precision (the cluster's 2^-103 separation is beyond it under
            // either setting of exact_residues_at_budget; 128 bits is not, with the control on), and one under 256
            // bits recovers the cluster. (Until 18 October 2026 a session stopped at 96 bits was resumed at 256;
            // a solve state is advanced once since.)
            options.max_work_bits = 96;
            options.work_units = 1000;
            d::SolveState limited_state(certpolysolve::detail::dense(*input.get()), d::maintained_method_policy(CPS_METHOD_SECULAR));
            const auto limited = limited_state.advance(&options);
            require(limited.complete && !limited.goal_met && limited.termination == CPS_PRECISION_LIMIT,
                    "the cluster did not stop for precision at 96 bits");
            options.max_work_bits = 256;
            d::SolveState raised_state(certpolysolve::detail::dense(*input.get()), d::maintained_method_policy(CPS_METHOD_SECULAR));
            const auto raised = raised_state.advance(&options);
            require(raised.goal_met && raised.complete && d::disjoint(raised.roots),
                    "a solve under a 256-bit budget did not recover the cluster");
        }
        // The candidate route's handoff starts Subdivision the same way: at the
        // goal's working precision, not at the budget its candidates reached.
        // Fuzz case 16333413833024016003: thirteen roots within 2^-60 of 1,
        // repeated with multiplicities 1, 2 and 6, 16-bit accuracy under a
        // 72-bit budget. The candidates spend the budget on the factors whose
        // roots are 2^-60 apart and hand them over.
        {
            const Rational e = power2(-62);
            std::vector<Complex> roots{Complex(1), Complex(1), Complex(1 + e), Complex(1 + 4 * e), Complex(1 - 4 * e)};
            for (int j = 0; j < 2; ++j) roots.emplace_back(1 - 3 * e);
            for (int j = 0; j < 6; ++j) roots.emplace_back(1 + 3 * e);
            const auto input = reference::input(from_roots(roots));
            auto policy = d::maintained_method_policy(CPS_METHOD_CLASSICAL);
            struct Progress { unsigned handoffs = 0; bool cancel = false; } progress;
            policy.trace_data = &progress;
            policy.trace = [](void *p, const d::TraceEvent &event) {
                if (std::string_view(event.stage) == "candidate_subdivision_handoff") {
                    auto &state = *static_cast<Progress *>(p);
                    ++state.handoffs;
                    state.cancel = true;
                }
            };
            cps_options options;
            options.goal = CPS_ROOT_ACCURACY;
            options.bits = 16;
            options.max_work_bits = 72;
            options.work_units = 1000;
            options.cancel_data = &progress;
            options.cancel = [](void *p) -> int { return static_cast<Progress *>(p)->cancel; };
            d::SolveState solver(certpolysolve::detail::dense(*input.get()), policy);
            const auto paused = solver.advance(&options);
            const uint64_t start = std::max<uint64_t>(64, options.bits + 32); // below the 72-bit budget
            require(progress.handoffs == 1 && paused.termination == CPS_CANCELLED,
                    "the candidate route did not hand its spent precision budget to subdivision");
            bool started = false;
            for (const auto &component : paused.components)
                started = started || (component.state == CPS_COMPONENT_ACTIVE && component.work_bits == start);
            require(started, "subdivision did not start at the goal's working precision after the candidate handoff");
        }
        // Dovetailing (dovetail.hpp): secular has no global convergence proof,
        // so a factor must finish through Subdivision's share whatever the
        // iteration does. Here secular keeps working but never moves a node
        // after its first certified cover; only the share can finish it.
        {
            // x^8 - 3x^5 + 2x^2 + 5: no rational roots (+-1, +-5 are not), so
            // exact preparation cannot certify any root on its own.
            ComplexPolynomial polynomial(9, Complex(0));
            polynomial[0] = {Rational(5)};
            polynomial[2] = {Rational(2)};
            polynomial[5] = {Rational(-3)};
            polynomial[8] = {Rational(1)};
            const uint64_t degree = 8;
            const auto input = reference::input(polynomial);
            for (const bool dovetail : {false, true}) {
                auto policy = d::maintained_method_policy(CPS_METHOD_SECULAR);
                policy.test_frozen_secular = true;
                policy.dovetail_subdivision = dovetail;
                // The guarantee holds for any share; one half keeps this test
                // short (at the production share of 1/16 the same solve takes
                // 11,432 work units).
                policy.dovetail_share_log2 = 1;
                // The Newton certificate and the preliminary ordinary phase,
                // both on the original polynomial, would finish these simple
                // roots before the secular model iterates at all.
                policy.newton_point_refinement = policy.secular_preliminary = false;
                cps_options options;
                options.goal = CPS_ROOT_ACCURACY;
                options.bits = 120; // beyond the initial binary64 certificate
                options.max_work_bits = 4096;
                // A work unit is one step, and the frozen iteration's steps are cheap: a request for the model
                // it already has computes nothing (SecularRefinement::step), so the share is earned over many
                // steps. The dovetailed solve took 23,708 (1,072 while such a request rebuilt the model, in the
                // same time); the undovetailed one is stopped at 3,000.
                options.work_units = dovetail ? 30000 : 3000;
                d::SolveState solver(certpolysolve::detail::dense(*input.get()), policy);
                auto result = solver.advance(&options);
                if (!dovetail) {
                    require(!result.goal_met && result.progress.stop_reason == CPS_STOP_WORK_BUDGET,
                            "the frozen secular iteration finished on its own; the test no longer exercises dovetailing");
                    continue;
                }
                require(result.goal_met && result.complete && result.termination == CPS_SUCCESS &&
                            d::disjoint(result.roots),
                        "dovetailed Subdivision did not finish a factor whose secular iteration never converges");
                const auto published = certpolysolve_test::published(std::move(result));
                uint64_t total = 0;
                for (const auto &region : published.roots()) {
                    const auto exact = rectangle_count(polynomial, region);
                    require(exact && static_cast<uint64_t>(*exact) == region.count,
                            "a dovetailed count disagrees with the exact winding count");
                    total += region.count;
                }
                require(total == degree, "the dovetailed cover lost roots");
            }
        }
        // A guard from the meter's qualification. Fuzz seed 510668236044846158
        // under Auto: 26 roots, a pair 2^-45 apart, at a 72-bit cap, solved
        // in 66 work units. Under a meter that charged Subdivision's steps
        // too little (a one-block cost table), the dovetailed search, engaged
        // early on a coarse cover and blocked at 72 bits by the cells that
        // cover forced, was kept at secular's stop and stopped for precision
        // at 86 units where a fresh start from secular's final cover
        // finishes. The block-structure meter removed the case; this keeps
        // it from returning.
        {
            std::vector<Complex> roots{
                {Rational("-37/65536"), Rational("-2079/4194304")},
                {Rational("5368709120"), Rational("4634017792")},
                {Rational("-3217143/32768")},
                {Rational("10240")},
                {Rational("-524288"), Rational("101711872/163")},
                {Rational("239685/262144"), Rational("761/1024")},
                {Rational("-5087/131072"), Rational("22273/262144")},
                {Rational("59/348127232"), Rational("-287849/4398046511104")},
                {Rational("2218745745375232"), Rational("1776810790486016")},
                {Rational("-14423/33554432"), Rational("772895/268435456")},
                {Rational("3356703/17592186044416")},
                {Rational("26212357206179840/879"), Rational("-31182452424704")},
                {Rational("372776480")},
                {Rational("14677101/72057594037927936")},
                {Rational("9292065/35184372088832")},
                {Rational("3925/1099511627776"), Rational("-151691/70368744177664")},
                {Rational("-1030147/302231454903657293676544"), Rational("4239017/2417851639229258349412352")},
                {Rational("460277/274877906944")},
                {Rational("1305/140737488355328"), Rational("605/70368744177664")},
                {Rational("26327/1099511627776"), Rational("-119/536870912")},
                {Rational("-122159/2048"), Rational("4")},
                {Rational("-10335409301094400"), Rational("-37735239065272320")},
                {Rational("55/77594624"), Rational("1/262144")},
                {Rational("9292063/35184372088832")},
                {Rational("9292061/35184372088832"), Rational("-1/35184372088832")},
                {Rational("0")}
            };
            const auto polynomial = from_roots(roots);
            const auto input = reference::input(polynomial);
            cps_options options;
            options.goal = CPS_ROOT_APPROXIMATION;
            options.bits = 53;
            options.max_work_bits = 72;
            options.work_units = 100000;
            options.live_cells = 65536;
            // With Subdivision as an equal partner at the budget (the cap is
            // the starting precision here; the recipe's rule since 1 October
            // 2026), its stop comes first.
            auto policy = d::maintained_method_policy(CPS_METHOD_SECULAR);
            d::SolveState solver(certpolysolve::detail::dense(*input.get()), policy);
            const auto result = solver.advance(&options);
            require(result.goal_met && result.complete && result.termination == CPS_SUCCESS &&
                        d::disjoint(result.roots),
                    "a dovetailed search kept at secular's stop blocked a solve a fresh start finishes");
            require(result.progress.work_units_this_call < 200, "the fresh start from secular's cover took long");
            uint64_t total = 0;
            for (const auto &root : result.roots) total += root.count;
            require(total == roots.size(), "the solve at a 72-bit cap lost roots");
        }
        // When the complete method itself has spent the precision budget, its
        // share guarantees nothing more, and a secular iteration that has
        // stopped making records must not run on to the work budget. Fuzz case
        // 12828861162913203822: one root and five within 2^-100 of each other,
        // under the accuracy goal at 100 bits. A 96-bit budget cannot
        // separate them (nor can a 192-bit one, unless the private control
        // exact_residues_at_budget is on); 384 bits can. The solve must stop
        // for precision long before its 100,000 work units (it used to spend
        // them all), and must continue to success when the budget is raised.
        {
            // (Secular now stops this solve itself: its model cannot steer the
            // five nodes. The dovetail rule on its own is checked below.)
            std::vector<Complex> roots{Complex(Rational("215539/262144"), Rational("-1/11"))};
            const char *near = "389022657345430751387049787391/633825300114114700748351602688";
            roots.emplace_back(Rational("-31057439705591620336669228531557/49121460758843889307997249208320"),
                               Rational(near));
            roots.emplace_back(Rational("-98/155"), Rational(near));
            roots.emplace_back(Rational("-124229758822366481346676914127003/196485843035375557231988996833280"),
                               Rational(near));
            roots.emplace_back(Rational("-124229758822366481346676914127313/196485843035375557231988996833280"),
                               Rational("1257/2048"));
            roots.emplace_back(Rational("-31057439705591620336669228531557/49121460758843889307997249208320"),
                               Rational("778045314690861502774099574783/1267650600228229401496703205376"));
            const auto polynomial = from_roots(roots);
            const auto input = reference::input(polynomial);
            cps_options options;
            options.goal = CPS_ROOT_ACCURACY;
            options.bits = 100;
            options.max_work_bits = 96;
            options.work_units = 100000;
            d::SolveState solver(certpolysolve::detail::dense(*input.get()), d::maintained_method_policy(CPS_METHOD_SECULAR));
            auto stopped = solver.advance(&options);
            require(stopped.complete && !stopped.goal_met && stopped.termination == CPS_PRECISION_LIMIT &&
                        stopped.progress.stop_reason == CPS_STOP_PRECISION,
                    "a factor whose complete method spent the precision budget did not stop for precision");
            // One tenth of the budget: far above what the stop needs, far below
            // the spin that spent all of it.
            require(stopped.progress.work_units_this_call < options.work_units / 10,
                    "the precision stop came only after a long spin");
            // Under a 384-bit budget a fresh solve succeeds (until 18 October 2026 the stopped state was resumed
            // with the budget raised; a solve state is advanced once since).
            options.max_work_bits = 384;
            d::SolveState raised(certpolysolve::detail::dense(*input.get()), d::maintained_method_policy(CPS_METHOD_SECULAR));
            auto result = raised.advance(&options);
            require(result.goal_met && result.complete && result.termination == CPS_SUCCESS &&
                        d::disjoint(result.roots),
                    "a solve under a 384-bit budget did not succeed");
            const auto published = certpolysolve_test::published(std::move(result));
            uint64_t total = 0;
            for (const auto &region : published.roots()) {
                const auto exact = rectangle_count(polynomial, region);
                require(exact && static_cast<uint64_t>(*exact) == region.count,
                        "a count after the raised budget disagrees with the exact winding count");
                total += region.count;
            }
            require(total == roots.size(), "the cover after the raised budget lost roots");
            // The dovetail rule alone: a secular iteration that keeps working,
            // never moves a node and never stops by itself (frozen after its
            // first certified cover, its residues accurate), beside a
            // Subdivision that spends the 192-bit budget on the five roots.
            // Once secular has gone longer without a record than all its work
            // before, the factor stops for precision.
            auto policy = d::maintained_method_policy(CPS_METHOD_SECULAR);
            policy.test_frozen_secular = true;
            policy.newton_point_refinement = policy.secular_preliminary = false;
            policy.dovetail_share_log2 = 1; // as above: any share; one half keeps the test short
            options.max_work_bits = 192;
            // In steps, as above: this stop took 19,076 (the frozen iteration's steps on the model it
            // already has are cheap, and the stop is earned by metered work).
            options.work_units = 40000;
            d::SolveState frozen(certpolysolve::detail::dense(*input.get()), policy);
            const auto blocked = frozen.advance(&options);
            require(blocked.complete && blocked.termination == CPS_PRECISION_LIMIT &&
                        blocked.progress.stop_reason == CPS_STOP_PRECISION,
                    "a frozen iteration beside a precision-blocked Subdivision did not stop for precision");
            // The same stop with the budget at the frozen iteration's own
            // precision, so the factor is at the budget throughout: the two
            // engines are then equal partners whatever the configured share
            // (dovetail.hpp, at_budget), and the stop at the production share
            // of a sixteenth costs what it costs at one half: 26,091 work
            // units either way (1,258 while a request for the same model
            // rebuilt it, where a sixteenth held to its share took 16,358).
            options.max_work_bits = blocked.work_bits;
            uint64_t units[2];
            for (const unsigned share_log2 : {1u, 4u}) {
                auto shared = policy;
                shared.dovetail_share_log2 = share_log2;
                d::SolveState at_budget(certpolysolve::detail::dense(*input.get()), shared);
                const auto stopped_at_budget = at_budget.advance(&options);
                require(stopped_at_budget.complete && stopped_at_budget.termination == CPS_PRECISION_LIMIT &&
                            stopped_at_budget.progress.stop_reason == CPS_STOP_PRECISION &&
                            stopped_at_budget.work_bits == options.max_work_bits,
                        "the frozen iteration at the budget did not stop for precision");
                units[share_log2 == 1 ? 0 : 1] = stopped_at_budget.progress.work_units_this_call;
            }
            require(units[1] <= 2 * units[0],
                    "at the budget the production share starved the precision stop of Subdivision's turns");
        }
        // Zero work returns without mathematical progress, and no precision beyond the caller's cap is reported: a
        // forced-piecewise solve with zero work used to enclose the input at max(64, bits), certify a stop cover there
        // and report that precision (review of 14 October 2026, S01: 64 bits under a 32-bit cap, 1,000 under a 64-bit
        // cap). Zero work is the options' internal setting since stage 3b-3 (a work budget of 0 selects the automatic
        // one); an internal solve whose budget is spent runs with it. (Moved from test_methods.)
        {
            certpolysolve::Polynomial p;
            p.set(0, certpolysolve::Number("1")).set(2, certpolysolve::Number("1")).set(5, certpolysolve::Number("3"));
            for (auto [bits, cap] : {std::pair<uint64_t, uint64_t>{53, 32}, {1000, 64}}) {
                const auto o = certpolysolve_test::raw_options();
                require(cps_options_set_method(o.get(), CPS_METHOD_PIECEWISE) == CPS_OK &&
                            cps_options_set_goal(o.get(), CPS_ROOT_APPROXIMATION) == CPS_OK &&
                            cps_options_set_precision(o.get(), bits) == CPS_OK &&
                            cps_options_set_max_work_bits(o.get(), cap) == CPS_OK,
                        "zero-work options refused");
                o.get()->work_units = 0;
                o.get()->work_units_default = false;
                const auto zero = certpolysolve_test::solve(p, o.get());
                const auto search = zero.statistics();
                require(cps_result_termination(zero.get()) == CPS_WORK_LIMIT && zero.statistics().work_units == 0 &&
                            search.work_units == 0 && zero.roots().size() == 0 && !zero.complete(),
                        "zero work made progress");
                require(zero.work_bits() <= cap, "zero work reported precision beyond the cap");
            }
        }
        std::cout << "owned classical/secular continuation preserves source, counts, coordinates and cancellation;"
                     " secular uses the machine/Newton certificate; restart motion does not hold the precision;"
                     " secular's spent precision budget continues in subdivision, also from an uninformed model\n";
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
