#pragma once
// For tests on the 0.5 C++ header that need what the header leaves out on purpose: certpolysolve_testing.h's poll
// counter, a region with imaginary bounds, an exact allowance alone (set_exact_allowance), a secular equation's
// resource limits. The test sets those on C options and solves through the C API; the answer comes back as the
// header's own Result, so the rest of the test reads it as any other.
#include <certpolysolve/certpolysolve.hpp>
#include <certpolysolve/certpolysolve_testing.h>
#include <memory>
#include <new>
#include <string>
#include <type_traits>
#include <utility>

namespace certpolysolve_test {
struct OptionsDeleter {
    void operator()(cps_options *o) const noexcept { cps_options_destroy(o); }
};
using RawOptions = std::unique_ptr<cps_options, OptionsDeleter>;
// Fresh C options with every default.
inline RawOptions raw_options() {
    cps_options *o = nullptr;
    if (cps_options_create(&o) != CPS_OK) throw std::bad_alloc();
    return RawOptions(o);
}
// What the exact-budget setter stage 3b-3 removed did, through the resource limits (C15): `bits` of exact storage,
// a compact exponent's expansion allowed, every other limit kept.
inline cps_error set_exact_allowance(cps_options *o, uint64_t bits) {
    cps_resource_limits limits;
    if (const cps_error e = cps_options_get_resource_limits(o, &limits); e != CPS_OK) return e;
    limits.exact_materialization_bits = bits;
    limits.flags |= CPS_RESOURCE_ALLOW_EXPONENT_EXPANSION;
    return cps_options_set_resource_limits(o, &limits);
}
// A solve with C options (NULL: the defaults), answered as certpolysolve::solve answers: an input that recorded a bad
// value gives its error result without a solve.
inline certpolysolve::Result solve(const certpolysolve::Polynomial &p, const cps_options *o) {
    if (!p.ok()) return certpolysolve::detail::error_result(p.status());
    cps_result *r = nullptr;
    const cps_error e = cps_solve(p.get(), o, &r);
    return certpolysolve::detail::make_result(e, r);
}
inline certpolysolve::Result solve(const certpolysolve::SecularEquation &s, const cps_options *o) {
    if (!s.ok()) return certpolysolve::detail::error_result(s.status());
    cps_result *r = nullptr;
    const cps_error e = cps_solve_secular(s.get(), o, &r);
    return certpolysolve::detail::make_result(e, r);
}
inline certpolysolve::Result solve(const certpolysolve::Program &p, const cps_options *o) {
    if (!p.ok()) return certpolysolve::detail::error_result(p.status());
    cps_result *r = nullptr;
    const cps_error e = cps_solve_program(p.get(), o, &r);
    return certpolysolve::detail::make_result(e, r);
}
// The C options the header's own solve makes from `o` (and from real_only and a region), to add what Options leaves
// out before solving through the C API: keep the object for the solve, since it owns the options and the link of
// o.stop to their token. status() holds the first refused setting (the header's solve would return it as an error).
class COptions {
    std::unique_ptr<certpolysolve::detail::Settings> settings_;

  public:
    explicit COptions(const certpolysolve::Options &o, bool real_only = false, const certpolysolve::Region *region = nullptr)
        : settings_(certpolysolve::detail::settings(o, real_only, region)) {}
    cps_options *get() const noexcept { return settings_->options.get(); }
    bool ok() const noexcept { return settings_->error.empty(); }
    const std::string &status() const noexcept { return settings_->error; }
};
// A C call's answer (its code and result, which the Result then owns) as the header's Result.
inline certpolysolve::Result from_c(cps_error e, cps_result *r) { return certpolysolve::detail::make_result(e, r); }
} // namespace certpolysolve_test
namespace certpolysolve::detail {
cps_result *published_result(cps_result result); // the library's internals (query.hpp): the answer and its text
} // namespace certpolysolve::detail
namespace certpolysolve_test {
// An answer a test made with the library's internal classes (a cps_result from a solve state; the type is complete
// only where internal.hpp is included), as the header's Result, the way the C entry points publish one.
template <class InternalResult> certpolysolve::Result published(InternalResult &&internal) {
    return certpolysolve::detail::make_result(CPS_OK,
                                          certpolysolve::detail::published_result(std::forward<InternalResult>(internal)));
}
} // namespace certpolysolve_test
