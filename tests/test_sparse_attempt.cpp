// A sparse polynomial with many terms through the public API. The certificate
// attempt both maintained recipes begin with must finish it, at its first
// precision: degree 1,600 with 320 nonzero coefficients, whose value ball
// Horner over the terms encloses with no correct bit at 85 bits (about half a
// bit lost a term). The attempt used to refuse its first point there and
// leave the solve to the ladders (2.5-4.5 s at 170 bits where the attempt
// needs 0.35 s).
//
// Every region is then checked without the library: value and slope at the
// region's centre, term by term at 256 bits with FLINT's own powering, and
// Rouche's theorem for the linear Taylor model on the disk inside the region
// and on the disk around it, with the coefficient majorant of p'' as the
// remainder bound. (A real root's region is an interval of the axis; there
// the disk on the interval serves for both.)
#include <flint/acb.h>
#include <flint/arb.h>
#include <flint/fmpz.h>
#include <certpolysolve/certpolysolve.hpp>
#include <iostream>
#include <map>
#include <set>

namespace {
void require(bool value, const char *why) {
    if (!value) throw std::runtime_error(why);
}
struct Acb {
    acb_t v;
    Acb() { acb_init(v); }
    ~Acb() { acb_clear(v); }
    Acb(const Acb &) = delete;
    Acb &operator=(const Acb &) = delete;
};
struct Arb {
    arb_t v;
    Arb() { arb_init(v); }
    ~Arb() { arb_clear(v); }
    Arb(const Arb &) = delete;
    Arb &operator=(const Arb &) = delete;
};
struct Generator {
    uint64_t state;
    uint64_t next() {
        state = state * 6364136223846793005ULL + 1442695040888963407ULL;
        return state >> 33;
    }
};
constexpr slong precision = 256;
void endpoint(arb_t out, const certpolysolve::Dyadic &x) {
    fmpz_t mantissa, exponent;
    fmpz_init(mantissa);
    fmpz_init(exponent);
    require(!fmpz_set_str(mantissa, x.mantissa().c_str(), 10) && !fmpz_set_str(exponent, x.exponent().c_str(), 10),
            "an endpoint is not a dyadic number");
    arb_set_fmpz(out, mantissa);
    arb_mul_2exp_fmpz(out, out, exponent);
    fmpz_clear(mantissa);
    fmpz_clear(exponent);
}
// |p(c)| + (r^2 / 2) max |p''| < |p'(c)| r on the disk |z - c| <= r, with
// max |p''| <= sum k (k - 1) |a_k| (|c| + r)^(k - 2): p and its linear model
// p(c) + p'(c)(z - c) then have the same number of zeros in the disk, and the
// model has one.
bool one_root(const std::map<slong, slong> &terms, const acb_t center, const arb_t radius) {
    Acb value, slope, power, term;
    Arb majorant, size, rho, left, right;
    acb_abs(rho.v, center, precision);
    arb_add(rho.v, rho.v, radius, precision);
    for (const auto &[exponent, coefficient] : terms) {
        acb_pow_ui(power.v, center, static_cast<ulong>(exponent), precision);
        acb_mul_si(term.v, power.v, coefficient, precision);
        acb_add(value.v, value.v, term.v, precision);
        acb_mul_si(term.v, term.v, exponent, precision);
        acb_add(slope.v, slope.v, term.v, precision); // z p'(z)
        if (exponent >= 2) {
            arb_pow_ui(size.v, rho.v, static_cast<ulong>(exponent - 2), precision);
            arb_mul_si(size.v, size.v, coefficient < 0 ? -coefficient : coefficient, precision);
            arb_mul_si(size.v, size.v, exponent, precision);
            arb_mul_si(size.v, size.v, exponent - 1, precision);
            arb_add(majorant.v, majorant.v, size.v, precision);
        }
    }
    acb_div(slope.v, slope.v, center, precision);
    acb_abs(left.v, value.v, precision);
    arb_mul(majorant.v, majorant.v, radius, precision);
    arb_mul(majorant.v, majorant.v, radius, precision);
    arb_mul_2exp_si(majorant.v, majorant.v, -1);
    arb_add(left.v, left.v, majorant.v, precision);
    acb_abs(right.v, slope.v, precision);
    arb_mul(right.v, right.v, radius, precision);
    return arb_lt(left.v, right.v);
}
} // namespace

int main() {
    certpolysolve::ThreadScope lifetime;
    try {
        constexpr slong degree = 1600, count = 320;
        Generator random{3};
        std::set<slong> support{0, degree};
        while (static_cast<slong>(support.size()) < count)
            support.insert(1 + static_cast<slong>(random.next() % static_cast<uint64_t>(degree - 1)));
        std::map<slong, slong> terms;
        certpolysolve::Polynomial input;
        for (const slong exponent : support) {
            const auto drawn = static_cast<slong>(random.next() % 19) - 9;
            const slong coefficient = drawn ? drawn : 1;
            terms[exponent] = coefficient;
            input.set(static_cast<uint64_t>(exponent), certpolysolve::Number(std::to_string(coefficient)));
        }
        bool checked = false;
        for (const auto method : {certpolysolve::Method::classical, certpolysolve::Method::secular, certpolysolve::Method::automatic}) {
            certpolysolve::Options options;
            options.method = method;
            options.goal = certpolysolve::Goal::approximation();
            options.bits = 53;
            options.max_bits = 4096;
            const auto result = certpolysolve::solve(input, options);
            require(result.goal_met() && cps_result_termination(result.get()) == CPS_SUCCESS && result.complete(),
                    "a sparse polynomial of 320 terms was not solved");
            const auto &regions = result.roots();
            require(regions.size() == static_cast<size_t>(degree), "the cover does not have one region for each root");
            for (const auto &region : regions)
                require(region.count == 1, "a region counts more than one root");
            require(result.engine_counters().validation_calls == 0,
                    "the solve did not end in the certificate attempt (a FLINT validation ran)");
            require(result.work_bits() == 85, "the certificates were not computed at the attempt's first precision");
            if (checked)
                continue;
            // The regions of one route, each against the polynomial itself.
            checked = true;
            for (const auto &region : regions) {
                Arb real_lower, real_upper, imag_lower, imag_upper, inner, outer, half;
                Acb center;
                endpoint(real_lower.v, region.box.real.lo);
                endpoint(real_upper.v, region.box.real.hi);
                endpoint(imag_lower.v, region.box.imag.lo);
                endpoint(imag_upper.v, region.box.imag.hi);
                arb_add(acb_realref(center.v), real_lower.v, real_upper.v, ARF_PREC_EXACT);
                arb_add(acb_imagref(center.v), imag_lower.v, imag_upper.v, ARF_PREC_EXACT);
                acb_mul_2exp_si(center.v, center.v, -1);
                // The half-widths: the smaller is the radius of the disk inside
                // the region, their sum bounds the radius of the disk around it.
                arb_sub(inner.v, real_upper.v, real_lower.v, ARF_PREC_EXACT);
                arb_sub(half.v, imag_upper.v, imag_lower.v, ARF_PREC_EXACT);
                arb_mul_2exp_si(inner.v, inner.v, -1);
                arb_mul_2exp_si(half.v, half.v, -1);
                arb_add(outer.v, inner.v, half.v, ARF_PREC_EXACT);
                if (arb_is_zero(half.v)) {
                    // A root proved real is exported as an interval of the
                    // axis. The polynomial is real: one root in the disk on
                    // that interval is its own conjugate, so it is real and
                    // lies in the interval.
                    require(region.realness == certpolysolve::Realness::real && arb_is_zero(acb_imagref(center.v)),
                            "a region without height is not a real root's interval");
                } else {
                    require(region.realness != certpolysolve::Realness::real, "a real root's region has height");
                    arb_min(inner.v, inner.v, half.v, precision);
                }
                require(arb_is_positive(inner.v), "a region has no interior");
                require(one_root(terms, center.v, inner.v), "a region is not shown to hold a root");
                require(one_root(terms, center.v, outer.v), "a region is not shown to hold only one root");
            }
        }
        std::cout << "sparse certificate attempt test passed\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "FAILED: " << error.what() << '\n';
        return 1;
    }
}
