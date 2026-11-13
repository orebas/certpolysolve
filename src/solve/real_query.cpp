#include "solve/real_query.hpp"
#include <algorithm>

namespace certpolysolve::detail::real_query {
bool exact_input(const DensePolynomial &input) {
    return std::all_of(input.coefficients.begin(), input.coefficients.end(),
                       [](const Coefficient &c) { return c.exact() && c.is_real(); });
}
bool exact_input(const cps_polynomial &input) {
    return std::all_of(input.terms.begin(), input.terms.end(),
                       [](const auto &t) { return t.second.exact() && t.second.is_real(); });
}
bool applies(const cps_options &o, const DensePolynomial &input) {
    return o.real_only && o.method == CPS_METHOD_AUTO && exact_input(input);
}
bool applies(const cps_options &o, const cps_polynomial &input) {
    return o.real_only && o.method == CPS_METHOD_AUTO && exact_input(input);
}

cps_options complex_options(cps_options o) {
    o.bits = o.goal == CPS_ROOT_ISOLATION ? cps_options{}.bits : o.bits;
    o.real_only = false;
    // A stage of the real-only question, which the route answers from this cover and may refine further: the
    // certificate attempt's cluster start stays out (cps_options::staged_goal; the fuzz's d1 case 10185732176165118585,
    // a cluster of three 2^-89 apart, took 17 driver steps with it where main takes 10).
    o.staged_goal = true;
    // An exact real factor whose roots are all real starts from the all-real seed pass's seeds instead of crawling
    // from the circles (docs/real-line-seed-pass-design-20261020.md).
    o.real_line_seeds = true;
    o.goal = CPS_ROOT_ACCURACY;
    o.region = {};
    return o;
}

std::unique_ptr<RealRegionFilter> region(const cps_options &o, const DensePolynomial &input) {
    if (!o.region.set()) return nullptr;
    std::vector<SparseRealTerm> terms;
    for (size_t j = 0; j < input.coefficients.size(); ++j)
        if (!input.coefficients[j].zero()) terms.push_back({static_cast<slong>(j), input.coefficients[j].real});
    return std::make_unique<RealRegionFilter>(o, std::move(terms));
}
std::unique_ptr<RealRegionFilter> region(const cps_options &o, const cps_polynomial &input) {
    if (!o.region.set()) return nullptr;
    std::vector<SparseRealTerm> terms;
    for (const auto &[exponent, coefficient] : input.terms) terms.push_back({static_cast<slong>(exponent), coefficient.real});
    return std::make_unique<RealRegionFilter>(o, std::move(terms));
}

} // namespace certpolysolve::detail::real_query
