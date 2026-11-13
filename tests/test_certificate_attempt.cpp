// The certificate attempt through the public API. Both maintained recipes
// begin with it (a binary64 Aberth packet, then a certificate at each point,
// with further packets on the points it refuses), and a solve it finishes
// makes no FLINT validation at all. The fixture needs more than the pass's
// precision on at least one recipe: the attempt must read the exact input at
// the precision it raises to (the classical recipe once gave it a fixed
// enclosure), and the result must report that precision.
#include "exact_winding_api.hpp"
#include <certpolysolve/certpolysolve.hpp>
#include <iostream>

using namespace reference;

int main() {
    certpolysolve::ThreadScope lifetime;
    try {
        // 100 rational points on each of the circles of radius 1 and 5/4:
        // simple roots, the closest about 0.02 apart, whose evaluation loses
        // up to 32 bits. Binary64 leaves the worst of them a few per cent
        // from their roots, and 85-bit certificates are refused for rounding.
        std::vector<Complex> planted;
        for (const auto &radius : {Rational(1), Rational(5) / 4})
            for (int k = 1; k <= 25; ++k) {
                const Rational t = Rational(k) / 26, scale = radius / (1 + t * t);
                const Rational x = (1 - t * t) * scale, y = 2 * t * scale;
                for (const auto &z : {Complex(x, y), Complex(-x, y), Complex(x, -y), Complex(-x, -y)})
                    planted.push_back(z);
            }
        const auto input = reference::input(from_roots(planted));
        uint64_t raised = 0;
        // Both accuracy goals: the attempt's cover isolates every root, which
        // the accuracy goal (the library's default) asks for as well.
        for (const auto &goal : {certpolysolve::Goal::approximation(), certpolysolve::Goal::accuracy()})
        for (const auto method : {certpolysolve::Method::classical, certpolysolve::Method::secular, certpolysolve::Method::automatic}) {
            certpolysolve::Options options;
            options.method = method;
            options.goal = goal;
            options.bits = 53;
            options.max_bits = 4096;
            const auto result = certpolysolve::solve(input, options);
            require(result.goal_met() && cps_result_termination(result.get()) == CPS_SUCCESS && result.complete(),
                    "a polynomial with 200 simple roots was not solved");
            uint64_t total = 0;
            for (const auto &region : result.roots()) {
                const Rational real_lower = endpoint(region.box.real.lo), real_upper = endpoint(region.box.real.hi);
                const Rational imag_lower = endpoint(region.box.imag.lo), imag_upper = endpoint(region.box.imag.hi);
                uint64_t inside = 0;
                for (const auto &root : planted)
                    inside += real_lower <= root.real && root.real <= real_upper && imag_lower <= root.imag &&
                              root.imag <= imag_upper;
                require(inside == region.count && region.count == 1, "a region does not hold exactly its planted root");
                total += region.count;
            }
            require(total == planted.size(), "the cover lost roots");
            require(result.engine_counters().validation_calls == 0,
                    "the solve did not end in the certificate attempt (a FLINT validation ran)");
            require(result.work_bits() >= 85 && result.work_bits() <= 4096,
                    "the result does not report the precision the certificates were computed at");
            raised += result.work_bits() > 85;
        }
        // Since the compensated tier's Newton retry (28 October 2026) the secular recipe's points are proved at the
        // pass's precision; the classical recipe's still need a raise, which its result must report.
        require(raised > 0, "no solve of the fixture reported a raised certificate precision");
        // The first point the attempt refuses may need two of its packets.
        // This polynomial of degree 33 (the squared distance from a point to
        // a Bezier curve, from the engineering controls) has eight points
        // that binary64 cannot settle; the first of them is not certified
        // after one packet of four sweeps and is after two. A round on that
        // one point used to end the attempt ("at least half" of one point),
        // and the solve went on to the ladder and its validations; it now
        // keeps the point for the next round, and the attempt finishes.
        {
            reference::Polynomial bezier;
            for (const char *coefficient :
                 {"-5797/256", "34459/32", "-11522243/448", "6074219/16", "-3048789095/768", "1114970075/32",
                  "-115108423903/384", "1704896110381/672", "-9732423576993/512", "78922019804485/672",
                  "-226124382322267/384", "154959667768217/64", "-3197461380268961/384", "588156183355117/24",
                  "-7078600181002695/112", "49038358056302095/336", "-549174402904911681/1792",
                  "9413600381794683/16", "-65778771880772409/64", "103356459534876645/64",
                  "-144077329184080845/64", "351774970180304957/128", "-371876273209725427/128",
                  "84264431190548739/32", "-1037813787952705725/512", "335893051512472831/256",
                  "-90419066273084445/128", "9993704305925445/32", "-14293867850618767/128",
                  "507066367676325/16", "-27823718194197/4", "72376923288413/64", "-127066474939569/1024",
                  "7195223171093/1024"})
                bezier.push_back(Rational(coefficient));
            const auto input = reference::input(bezier);
            for (const auto &goal : {certpolysolve::Goal::approximation(), certpolysolve::Goal::accuracy()})
            for (const auto method : {certpolysolve::Method::classical, certpolysolve::Method::secular, certpolysolve::Method::automatic}) {
                certpolysolve::Options options;
                options.method = method;
            options.goal = goal;
            options.bits = 53;
            options.max_bits = 4096;
                const auto result = certpolysolve::solve(input, options);
                require(result.goal_met() && cps_result_termination(result.get()) == CPS_SUCCESS && result.complete(),
                        "the Bezier distance polynomial was not solved");
                uint64_t total = 0;
                for (const auto &region : result.roots()) total += region.count;
                require(total == bezier.size() - 1, "the Bezier distance polynomial's cover lost roots");
                require(result.engine_counters().validation_calls == 0,
                        "the attempt ended at a first refused point that needed a second packet");
            }
        }
        // The rounds come after the pass. MPSolve's chrma86 (degree 85):
        // binary64 leaves a point among the first few that the certificate
        // refuses, and a round then, set against the few certificates held,
        // was refused by the balance and ended the attempt; the solve went
        // on to the ladder and its validations. After the pass the round is
        // set against every certificate, and the attempt finishes.
        {
            reference::Polynomial chrma86;
            for (const char *coefficient :
                 {"-37323563944774684366413", "1421488044269587510553502", "-26736753027456629518691835",
                  "331131917982836213928263988", "-3037787541450031336549542855", "22018752959335124409905815062",
                  "-131348126008031306354473934673", "663247089279227279850889124112",
                  "-2893969841637200007409915485946", "11084431184857667532985950195340",
                  "-37733325581815351761156010417102", "115316863319364359670373342858520",
                  "-319023440352616323904295490736326", "804536900111299256164858897853900",
                  "-1860584098829215600236203470430274", "3966135801140013878425374322412992",
                  "-7828063640729597190937835338161389", "14362670950002180848565058293977694",
                  "-24583835110302873412208118637823811", "39380784031550800205531673731033908",
                  "-59210764237552010680390560543626527", "83783223904243710783570912272460674",
                  "-111848058903521517223938772237552989", "141196005694514322375884524945567648",
                  "-168924678100849441586888690527074296", "191934355684789350658715635188069488",
                  "-207529802799946780673414419809639384", "213958040975255756704422401895883072",
                  "-210729340913306303176850451120179464", "198642266669751717285042389936645296",
                  "-179530671395226377420085992417733080", "155830840583051867732289999886407448",
                  "-130102534106192830626941611459245026", "104623685251901350570308210774003732",
                  "-81130911771324285222331733833967214", "60721722751669841823234237554944080",
                  "-43890422571558290599065897455493430", "30648028199839735564175746760112260",
                  "-20675584960953223583060935969704890", "13472210443671907508240815261549264",
                  "-8475294779706316125174305742084868", "5144487693046667206325008305832416",
                  "-3010834527207528735366020874011474", "1697638629644615452705976169070808",
                  "-921426113335525764896515559785512", "481034063400901647743131872281000",
                  "-241347328851242444803448993231742", "116285169746791043529510051319304",
                  "-53764891560915180763283384684810", "23837215437419433715128912120988",
                  "-10127294077652926185677801840720", "4120180398235404845479768511272",
                  "-1604086200506194392837821143706", "597208944236533091195912455324",
                  "-212471068653559962033678740264", "72180488770290051910410589216",
                  "-23395718702945191146889883064", "7228905698367478566499054336",
                  "-2127251879937399703739813620", "595563471353731861856952208", "-158456039418808447335444048",
                  "40014750144061982090660208", "-9577771674280587490372388", "2169604691922977357956540",
                  "-464338998088665475693365", "93715883286149103976346", "-17799469352226551798303",
                  "3173970799687484152320", "-529986455049921656391", "82625566344162047514",
                  "-11986905090030193013", "1612129192377644928", "-200130039240030010", "22817831845317932",
                  "-2375516438647068", "224276849409560", "-19045648342274", "1440364708500", "-95821759548",
                  "5520463004", "-269833373", "10880130", "-347461", "8240", "-129", "1"})
                chrma86.push_back(Rational(coefficient));
            const auto input = reference::input(chrma86);
            for (const auto &goal : {certpolysolve::Goal::approximation(), certpolysolve::Goal::accuracy()})
            for (const auto method : {certpolysolve::Method::classical, certpolysolve::Method::secular, certpolysolve::Method::automatic}) {
                certpolysolve::Options options;
                options.method = method;
            options.goal = goal;
            options.bits = 53;
            options.max_bits = 4096;
                const auto result = certpolysolve::solve(input, options);
                require(result.goal_met() && cps_result_termination(result.get()) == CPS_SUCCESS && result.complete(),
                        "chrma86 was not solved");
                uint64_t total = 0;
                for (const auto &region : result.roots()) total += region.count;
                require(total == chrma86.size() - 1, "chrma86's cover lost roots");
                require(result.engine_counters().validation_calls == 0,
                        "an early refusal's round ended the attempt (a FLINT validation ran)");
            }
        }
        // A goal above the certificate ladder's first rung that the attempt
        // refuses from initial candidates is met in two steps of one solve:
        // secular iterates to the library's default goal, then retries the
        // attempt from those points at the caller's goal, starting at the
        // caller's working precision. Chebyshev's T_40 (leading coefficient
        // 2^39): binary64 leaves its points too far from its roots for a
        // 1,024-bit certificate, and the refusals' allowance runs out after a
        // few; secular then iterated at four times the goal's working
        // precision (4,352 bits, 8 work units, 7 FLINT validations). With the
        // two steps the caller's goal costs two work units beyond the default
        // goal's: the refused attempt and the retry, which certifies every
        // point (2,112 bits). A work budget bounds both steps.
        {
            reference::Polynomial previous{Rational(1)}, chebyshev{Rational(0), Rational(1)};
            for (int n = 1; n < 40; ++n) {
                reference::Polynomial next(chebyshev.size() + 1, Rational(0));
                for (size_t k = 0; k < chebyshev.size(); ++k) next[k + 1] += 2 * chebyshev[k];
                for (size_t k = 0; k < previous.size(); ++k) next[k] -= previous[k];
                previous = std::move(chebyshev);
                chebyshev = std::move(next);
            }
            const auto input = reference::input(chebyshev);
            constexpr uint64_t goal = 1024;
            const uint64_t working = goal + 32; // the driver's guard bits
            certpolysolve::Options options;
            options.goal = certpolysolve::Goal::accuracy();
            options.bits = 53;
            options.max_bits = 1U << 16;
            const auto initial = certpolysolve::solve(input, options);
            require(initial.goal_met(), "T_40 was not solved at the default goal");
            options.bits = goal;
            options.max_bits = 1U << 16;
            const auto result = certpolysolve::solve(input, options);
            require(result.goal_met() && cps_result_termination(result.get()) == CPS_SUCCESS && result.complete(),
                    "T_40 was not solved at 1,024 bits");
            uint64_t total = 0;
            for (const auto &region : result.roots()) total += region.count;
            require(total == chebyshev.size() - 1, "T_40's cover lost roots");
            require(result.work_bits() < 4 * working, "a 1,024-bit goal iterated at four times its working precision");
            require(result.statistics().work_units <= initial.statistics().work_units + 2,
                    "the caller's goal cost more than two attempts beyond the default goal");
            options.work_budget = 3;
            options.cell_budget = CPS_DEFAULT_LIVE_CELLS;
            const auto first = certpolysolve::solve(input, options);
            require(cps_result_termination(first.get()) == CPS_WORK_LIMIT && first.complete() && !first.goal_met(),
                    "a spent budget did not return a complete cover");
        }
        // A refused attempt keeps its certificates. Two roots 2^-100 apart at
        // 1/64 (a cluster at a 53-bit approximation goal) and 28 simple roots,
        // rational points on the unit circle: the attempt certifies the 28 and
        // refuses the pair, which no region of one root splits at its
        // precision. Thrown away, the 28 were iterated again by secular at the
        // precision the pair needs, and their regions were secular's; kept,
        // they are the attempt's in the cover, and secular iterates the pair
        // alone. Under the accuracy goal secular splits the pair. Every
        // region holds exactly its count of the planted roots.
        {
            std::vector<Complex> roots{Complex(Rational(1, 64)), Complex(Rational(1, 64) + power2(-100))};
            for (int k = 1; k <= 7; ++k) {
                const Rational t = Rational(k) / 8, scale = Rational(1) / (1 + t * t);
                const Rational x = (1 - t * t) * scale, y = 2 * t * scale;
                for (const auto &z : {Complex(x, y), Complex(-x, y), Complex(x, -y), Complex(-x, -y)})
                    roots.push_back(z);
            }
            const auto input = reference::input(from_roots(roots));
            for (const auto &goal : {certpolysolve::Goal::approximation(), certpolysolve::Goal::accuracy()}) {
                certpolysolve::Options options;
                options.method = certpolysolve::Method::secular;
                options.goal = goal;
                options.bits = 53;
                options.max_bits = 4096;
                const auto result = certpolysolve::solve(input, options);
                require(result.goal_met() && cps_result_termination(result.get()) == CPS_SUCCESS && result.complete(),
                        "a pair 2^-100 apart beside 28 simple roots was not solved");
                uint64_t total = 0, kept = 0;
                for (const auto &region : result.roots()) {
                    const Rational real_lower = endpoint(region.box.real.lo), real_upper = endpoint(region.box.real.hi);
                    const Rational imag_lower = endpoint(region.box.imag.lo), imag_upper = endpoint(region.box.imag.hi);
                    uint64_t inside = 0;
                    for (const auto &root : roots)
                        inside += real_lower <= root.real && root.real <= real_upper && imag_lower <= root.imag &&
                                  root.imag <= imag_upper;
                    require(inside == region.count, "a region does not hold exactly its count of planted roots");
                    total += region.count;
                    kept += region.certificate.rfind("original-coefficient", 0) == 0;
                }
                require(total == roots.size(), "the cover lost roots");
                require(kept == roots.size() - 2, "the refused attempt's 28 certificates were not kept");
            }
        }
        std::cout << "certificate attempt tests passed\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "FAILED: " << error.what() << '\n';
        return 1;
    }
}
