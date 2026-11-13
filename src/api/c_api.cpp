// The public C API: every function that include/certpolysolve/certpolysolve.h and certpolysolve_testing.h declare. Each checks its
// arguments, converts them to the solver's types and hands back the contract's objects, results and error codes.
//
// Every function that returns cps_error, except cps_result_get_approximation, runs in guarded(), which clears the
// thread's message and maps std::bad_alloc, std::invalid_argument (thrown by require()), std::out_of_range and anything
// else to CPS_OUT_OF_MEMORY, CPS_INVALID_ARGUMENT, CPS_OUT_OF_RANGE and CPS_INTERNAL_ERROR, with a message for
// cps_last_error. Handle outputs are NULL and node outputs CPS_NO_NODE unless the call succeeds, a failed setter or
// builder leaves its object unchanged, and NULL options mean the defaults (or_defaults). Literals become exact Numbers,
// kept as a polynomial's nonzero terms (store) or made exact rationals under the resource limits (cps_secular_create).
// A solve (cps_solve, cps_solve_secular, cps_solve_program) runs in the default floating-point environment, then restores
// the caller's (DefaultFloatingPoint), refuses the options no input of its kind can answer (validate_query), and fixes
// its deadline on a copy of the options (timed_options). It advances one SessionState, which borrows the polynomial's
// terms for the call (cps_solve, via certpolysolve::detail::solve) or copies its input (the other two). finish_call puts the
// regions in exact midpoint order (MidpointOrder), gives a cancelled stop the stop reason CPS_STOP_TIME_LIMIT when the
// deadline has passed and no token was requested (name_time_limit), and records the call's wall time.
//
// State: the thread's last error (thread_local error_message, 512 bytes) and the const defaults (default_options). A
// solve only reads its input and options; its deadline and testing poll counter live in its copy and frame. Each object
// is the caller's until its destroy deletes it (layouts: internal.hpp, secular.hpp, program.hpp), a result's views live
// as long as it, and a token shares its std::atomic<bool> flag, asserted lock-free, with the options naming it.
// Entry points: cps_abi_version, cps_version, cps_backend_version, cps_last_error, cps_cleanup_thread (flint_cleanup);
// cps_polynomial_*, cps_secular_*, cps_program_*; cps_options_*, cps_resource_limits_default, cps_cancel_token_*; the solves;
// cps_result_*, cps_termination_name, cps_stop_reason_name; cps_dyadic_*; cps_testing_*. They call certpolysolve::detail::solve
// and compare_numbers (solver.cpp), SessionState (session.cpp), RealRegionFilter::exact_ends (real_region.hpp),
// polynomial_add and kin (polynomial_arithmetic.cpp), Program::append and parse_program (program.cpp), parse_dyadic and
// kin (exact_text.cpp), and statistics::carried_stage_fields (engine_statistics.hpp).
// Not here: routing, the goal per input kind, input-dependent refusals, a program's output check (session.cpp); literal
// parsing (number.cpp); the cancellation polls (cps_options::cancelled, internal.hpp); the C++ API (certpolysolve.hpp).
// Design notes: docs/api-spec-0.5.md (errors 2.1, NULL rules 2.3, usage errors 4.3 and 5.2, region order 6.3,
// statistics 7, threads 8; the comments' Cn labels are its change list, section 10), docs/api.md (the guide to the
// built library) and docs/thread-safety-audit-20261018.md (C3 and E4, the floating-point environment).
#include "platform/floating_environment.hpp"
#include "platform/deadline.hpp"
#include "core/cover.hpp"
#include "platform/parallel.hpp"
#include "input/coefficients/polynomial_input.hpp"
#include "core/query.hpp"
#include "solve/solve_state.hpp"
#include "numeric/arithmetic/polynomial_arithmetic.hpp"
#include "numeric/arithmetic/exact_text.hpp"
#include "input/program/program.hpp"
#include "input/secular/secular.hpp"
#include "certify/real_region.hpp"
#include "core/session_statistics.hpp"
#include <certpolysolve/certpolysolve_testing.h>
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <thread>
#include <type_traits>

namespace {
thread_local char error_message[512] = "";
void error_text(const char *text) noexcept {
    std::snprintf(error_message, sizeof(error_message), "%s", text);
}
template <class Function> cps_error guarded(Function &&function) noexcept {
    // The thread's message describes its last failed call only: a call that succeeds leaves it empty.
    error_message[0] = '\0';
    try {
        function();
        return CPS_OK;
    } catch (const std::bad_alloc &) {
        error_text("allocation failed");
        return CPS_OUT_OF_MEMORY;
    } catch (const std::invalid_argument &e) {
        error_text(e.what());
        return CPS_INVALID_ARGUMENT;
    } catch (const std::out_of_range &e) {
        error_text(e.what());
        return CPS_OUT_OF_RANGE;
    } catch (const std::exception &e) {
        error_text(e.what());
        return CPS_INTERNAL_ERROR;
    } catch (...) {
        error_text("unknown C++ exception");
        return CPS_INTERNAL_ERROR;
    }
}
void require(bool predicate, const char *message) {
    if (!predicate)
        throw std::invalid_argument(message);
}
// An enum argument of the C API as the integer the caller passed. A C caller may pass any int, and the contract
// refuses one outside the enumeration; in C++ reading the parameter as its enum type is undefined for such a value
// (UBSan's enum check reports it), so the bytes are copied into the enumeration's underlying integer type and the
// range is checked there. Used before any read of the parameter.
// History: docs/code-history-20261021.md, c_api.cpp 1.
template <class Enum> long long enum_argument(const Enum &argument) {
    using Underlying = std::underlying_type_t<Enum>;
    static_assert(std::is_enum_v<Enum> && sizeof(Underlying) < sizeof(long long)); // widened without loss
    Underlying value;
    std::memcpy(&value, &argument, sizeof value);
    return value;
}
// An exponent a polynomial can hold (C40: any, up to the largest a degree can name).
void valid_exponent(uint64_t index) {
    if (index > certpolysolve::detail::max_exponent)
        throw std::out_of_range("an exponent above 2^63 - 1 cannot be a polynomial's degree");
}
// c as the coefficient of x^index: an exact zero removes the term, since the object keeps its nonzero terms only.
void store(cps_polynomial &p, uint64_t index, certpolysolve::detail::Coefficient &&c) {
    if (c.zero()) p.terms.erase(index);
    else p.terms.insert_or_assign(index, std::move(c));
}

void validate_resources(const cps_resource_limits &r) {
    require(!(r.flags & ~CPS_RESOURCE_ALLOW_EXPONENT_EXPANSION), "unknown resource flag");
}
// Appends a node to the program, transactionally (Program::append checks before it changes anything). These
// helpers return and take C++ types, so they live outside the extern "C" block.
cps_error program_node(cps_program *p, uint64_t *out, const certpolysolve::detail::ProgramNode &node) {
    if (out)
        *out = CPS_NO_NODE;
    return guarded([&] {
        require(p && out, "null program or node output");
        *out = p->program.append(node);
    });
}
static_assert(std::atomic<bool>::is_always_lock_free, "a cancellation request must be a lock-free store");
static_assert(CPS_DEFAULT_WORK_UNITS_PER_ROOT == certpolysolve::detail::heuristics::default_work_units_per_root,
              "the header's per-root default is the solver's");
// certpolysolve_testing.h's poll counter for one solve call: from its stop_at-th poll on the call is cancelled, as by a
// callback that counts its calls; before that it asks the options' own callback, if any (the internal one: no public
// call sets it). It is installed as the call's callback (timed_options), so the solver polls it exactly where it
// polls a callback, through every internal wrapper and clear of the callback, and it lives in the call's frame for
// as long as the call's options. History: docs/code-history-20261021.md, c_api.cpp 2.
struct TestingPolls {
    uint64_t polls = 0, stop_at = 0;
    certpolysolve::detail::CancelCallback chained = nullptr;
    void *chained_data = nullptr;
    static int poll(void *data) {
        auto &t = *static_cast<TestingPolls *>(data);
        if (++t.polls >= t.stop_at) return 1;
        return t.chained ? t.chained(t.chained_data) : 0;
    }
};
// The options a solve call works with: the caller's (or the defaults), with the call's deadline fixed from the time
// limit at its entry. Every solve the call runs copies these options and so inherits that deadline. With the testing
// counter on, `polls` becomes the call's callback; off, the options are the caller's as they are.
cps_options timed_options(const cps_options *o, TestingPolls &polls) {
    cps_options timed = o ? *o : cps_options{};
    if (timed.time_limit > 0 && std::isfinite(timed.time_limit))
        timed.deadline = certpolysolve::detail::deadline_after(timed.time_limit, std::chrono::steady_clock::now());
    if (timed.testing_stop_at_poll) {
        polls.stop_at = timed.testing_stop_at_poll;
        polls.chained = timed.cancel;
        polls.chained_data = timed.cancel_data;
        timed.cancel = TestingPolls::poll;
        timed.cancel_data = &polls;
    }
    return timed;
}
// The threads a solve call starts for its split loops (W9, parallel.hpp): the options' count, at most the machine's
// hardware threads (more cannot run at once, and the answer is the same at every count; an unknown count, 0, sets no
// bound), or one where a callback is installed (the testing poll counter), which no two threads may call at once.
unsigned call_workers(const cps_options &timed) {
    if (timed.cancel || timed.workers <= 1) return 1; // the default asks the system nothing
    const unsigned hardware = std::thread::hardware_concurrency();
    return hardware ? std::min<unsigned>(timed.workers, hardware) : timed.workers;
}
// The usage errors (docs/api-spec-0.5.md sections 4.3 and 5.2; C7, C8): settings that no input of the entry point's
// kind can answer, and a region whose ends are out of order or not exact within the exact storage policy. The
// options alone decide them, before any work; each is CPS_INVALID_ARGUMENT from the solve, with no result. (A refusal
// that depends on the input, a forced method on ball input for one, is an answer: CPS_UNSUPPORTED_QUERY.)
enum class EntryKind { Coefficients, Secular, Program };
void validate_query(const cps_options &o, EntryKind kind) {
    const bool forced = o.method != CPS_METHOD_AUTO, cluster = o.goal == CPS_CLUSTER_RADIUS,
               isolation = o.goal == CPS_ROOT_ISOLATION, region = o.region.set();
    switch (kind) {
    case EntryKind::Coefficients:
        require(o.real_only || !region, "a region is for real-only queries");
        require(!(forced && (o.real_only || isolation)),
                "real-only and isolation queries are answered by CPS_METHOD_AUTO only");
        require(!(forced && cluster), "the cluster-radius goal is answered by CPS_METHOD_AUTO only");
        require(!(o.real_only && cluster), "a real-only query does not take the cluster-radius goal");
        break;
    case EntryKind::Secular:
        require(o.method == CPS_METHOD_AUTO || o.method == CPS_METHOD_SECULAR,
                "native secular input is solved by CPS_METHOD_AUTO or CPS_METHOD_SECULAR");
        require(!isolation && !region, "native secular input takes no isolation goal and no region");
        require(!(forced && (o.real_only || cluster)),
                "forced secular answers all complex roots, at the accuracy or the approximation goal");
        require(!(o.real_only && cluster), "a real-only query does not take the cluster-radius goal");
        break;
    case EntryKind::Program:
        require(o.method != CPS_METHOD_PIECEWISE,
                "program input is solved by CPS_METHOD_AUTO, CPS_METHOD_CLASSICAL or CPS_METHOD_SECULAR");
        require(!o.real_only && !isolation && !region,
                "program input answers all complex roots: no real-only query, isolation goal or region");
        break;
    }
    if (!region) return;
    // The ends' order, decided exactly within the precision cap, and the real ends' exact form.
    const auto ordered = [&](const std::optional<certpolysolve::detail::Number> &lower,
                             const std::optional<certpolysolve::detail::Number> &upper) {
        if (!lower || !upper) return true;
        const auto order = certpolysolve::detail::compare_numbers(*lower, *upper, o);
        require(order.has_value(), "the order of a region's ends could not be decided within the precision cap");
        return *order <= 0;
    };
    require(ordered(o.region.real_lower, o.region.real_upper) && ordered(o.region.imag_lower, o.region.imag_upper),
            "a region's lower end exceeds its upper end");
    std::array<std::optional<certpolysolve::detail::Rational>, 2> ends;
    if (const char *refused = certpolysolve::detail::RealRegionFilter::exact_ends(o, ends))
        throw std::invalid_argument(refused);
}
// The size protocol of cps_statistics and cps_engine_counters: the caller's size must hold the 0.5 layout; the
// library fills that layout, leaves any bytes beyond it as they are and reports the bytes it filled.
template <class Block> void fill_sized(Block *out, const Block &value) {
    require(out, "null statistics output");
    require(out->size >= sizeof(Block), "the size field is smaller than this library's layout (set it to sizeof)");
    Block filled = value;
    filled.size = sizeof(Block);
    std::memcpy(out, &filled, sizeof(Block));
}
// A literal as the getters give it back: surrounding white space removed.
std::string trimmed(const char *text) {
    std::string s(text);
    const auto first = s.find_first_not_of(" \t\n\r\f\v");
    if (first == std::string::npos) return {};
    return s.substr(first, s.find_last_not_of(" \t\n\r\f\v") - first + 1);
}
// A solve stopped by its checkpoints says which question stopped it: the time limit when the deadline has passed
// and the token was not requested (a callback and a limit together are told apart by the clock alone).
void name_time_limit(cps_result &result, const cps_options &o) {
    if (result.progress.stop_reason != CPS_STOP_CANCELLED || !o.deadline) return;
    if (o.cancel_flag && o.cancel_flag->load(std::memory_order_relaxed)) return;
    if (std::chrono::steady_clock::now() < *o.deadline) return;
    result.progress.stop_reason = CPS_STOP_TIME_LIMIT;
    result.diagnostics = "the time limit was reached; " + result.diagnostics;
}
// The published order of the regions (C26, O13): increasing exact midpoint, by its real part and then its imaginary
// part, on every route. The regions are disjoint closed rectangles, so no two share a midpoint, and the sum of a
// region's two ends on an axis orders as its midpoint does. Reading an end back from its decimal text costs as much as
// a small solve's last few percent, so the order is built cheaply and proved:
// 1. sort (indices) by the boxes' centres; a route's list is often in that order already;
// 2. prove each neighbouring pair's order from the boxes, which hold the exported rectangles (Root::publish and the
//    later edits keep them so) and so their exact midpoints: a part's boxes lie apart when the centres' difference,
//    rounded toward zero, exceeds the sum of the radii, rounded up; the pair is in order when its real parts' boxes
//    lie apart, or when the two real intervals are the same (the same canonical ends) and their imaginary parts'
//    boxes lie apart;
// 3. a pair the boxes cannot order (a conjugate pair's real intervals differ in their last bits) is compared exactly,
//    from the ends read back, the real ends first and the imaginary ones only on a tie;
// 4. if a pair is out of order, the list is sorted by the exact keys.
// A chain of neighbours each in exact order is the exact order; a region moves at most once, at the end.
class MidpointOrder {
    using Root = certpolysolve::detail::Root;
    std::vector<Root> &roots_;
    // The exact sums of a region's ends on each axis, read on demand.
    struct Sums {
        arf_t real, imag;
        bool have_real = false, have_imag = false;
        Sums() { arf_init(real), arf_init(imag); }
        ~Sums() { arf_clear(real), arf_clear(imag); }
        Sums(const Sums &) = delete;
        Sums &operator=(const Sums &) = delete;
    };
    std::unique_ptr<Sums[]> sums_;
    static void read(const certpolysolve::detail::Endpoint &e, arf_t out) { arf_set(out, e.value()); }
    const arf_struct *sum(size_t i, bool real) {
        if (!sums_) sums_.reset(new Sums[roots_.size()]);
        Sums &s = sums_[i];
        bool &have = real ? s.have_real : s.have_imag;
        arf_ptr out = real ? s.real : s.imag;
        if (!have) {
            const auto &e = roots_[i].endpoints;
            arf_t lo, hi;
            arf_init(lo);
            arf_init(hi);
            read(e[real ? 0 : 2], lo), read(e[real ? 1 : 3], hi);
            arf_add(out, lo, hi, ARF_PREC_EXACT, ARF_RND_DOWN);
            arf_clear(lo);
            arf_clear(hi);
            have = true;
        }
        return out;
    }
    bool exactly_before(size_t a, size_t b) {
        const int real = arf_cmp(sum(a, true), sum(b, true));
        return real ? real < 0 : arf_cmp(sum(a, false), sum(b, false)) < 0;
    }
    // Every point of x lies below every point of y: y's centre less x's, rounded toward zero, is positive and above
    // the radii's sum, rounded up.
    static bool apart(const arb_t x, const arb_t y) {
        arf_t difference;
        arf_init(difference);
        arf_sub(difference, arb_midref(y), arb_midref(x), MAG_BITS, ARF_RND_DOWN);
        mag_t radii;
        mag_init(radii);
        mag_add(radii, arb_radref(x), arb_radref(y));
        const bool below = arf_sgn(difference) > 0 && arf_cmpabs_mag(difference, radii) > 0;
        mag_clear(radii);
        arf_clear(difference);
        return below;
    }
    static bool same(const certpolysolve::detail::Endpoint &x, const certpolysolve::detail::Endpoint &y) {
        return arf_equal(x.value(), y.value());
    }
    // The pair's order, from the boxes where they decide it, else exactly.
    bool before(size_t a, size_t b) {
        const Root &x = roots_[a], &y = roots_[b];
        if (apart(acb_realref(x.box.value()), acb_realref(y.box.value()))) return true;
        if (same(x.endpoints[0], y.endpoints[0]) && same(x.endpoints[1], y.endpoints[1]) &&
            apart(acb_imagref(x.box.value()), acb_imagref(y.box.value())))
            return true;
        return exactly_before(a, b);
    }

  public:
    explicit MidpointOrder(std::vector<Root> &roots) : roots_(roots) {}
    void apply() {
        const size_t n = roots_.size();
        if (n < 2) return;
        std::vector<size_t> order(n);
        std::iota(order.begin(), order.end(), size_t(0));
        const auto by_centre = [&](size_t a, size_t b) {
            const auto &x = roots_[a].box, &y = roots_[b].box;
            const int real = arf_cmp(arb_midref(acb_realref(x.value())), arb_midref(acb_realref(y.value())));
            return real ? real < 0 : arf_cmp(arb_midref(acb_imagref(x.value())), arb_midref(acb_imagref(y.value()))) < 0;
        };
        if (!std::is_sorted(order.begin(), order.end(), by_centre)) std::stable_sort(order.begin(), order.end(), by_centre);
        bool ordered = true;
        for (size_t i = 1; i < n && ordered; ++i) ordered = before(order[i - 1], order[i]);
        if (!ordered) std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) { return exactly_before(a, b); });
        bool identity = true;
        for (size_t i = 0; i < n && identity; ++i) identity = order[i] == i;
        if (identity) return;
        std::vector<Root> sorted;
        sorted.reserve(n);
        for (const size_t i : order) sorted.push_back(std::move(roots_[i]));
        roots_ = std::move(sorted);
    }
};
// What every solve entry does once its answer is formed: the regions put in their published order (MidpointOrder),
// the time-limit stop named, the call's wall time and the testing counter's polls recorded.
void finish_call(cps_result &result, const cps_options &timed, const TestingPolls &polls,
                 std::chrono::steady_clock::time_point started) {
    MidpointOrder(result.roots).apply();
    certpolysolve::detail::publish_text(result);
    name_time_limit(result, timed);
    result.total_ns = static_cast<uint64_t>(
        std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - started).count());
    result.testing_polls = polls.polls;
    result.testing_split_jobs = certpolysolve::detail::split_jobs();
}
certpolysolve::detail::ProgramNode program_op(certpolysolve::detail::ProgramOp op, uint64_t a = 0, uint64_t b = 0,
                                          uint64_t value = 0) {
    certpolysolve::detail::ProgramNode node;
    node.op = op;
    node.a = a;
    node.b = b;
    node.value = value;
    return node;
}
// A new polynomial from an arithmetic call (C39): *out stays NULL until the result is complete.
template <class Make> cps_error polynomial_result(cps_polynomial **out, Make &&make) {
    if (out)
        *out = nullptr;
    return guarded([&] {
        require(out, "null polynomial output");
        auto result = std::make_unique<cps_polynomial>(make());
        *out = result.release();
    });
}
// Text into a caller's buffer: at most size bytes with the NUL; *length (if not NULL) gets the full length without
// the NUL, so size 0 asks for it.
void write_text(const std::string &text, char *buffer, size_t size, size_t *length) {
    require(buffer || !size, "null buffer with a nonzero size");
    if (length)
        *length = text.size();
    if (size) {
        const size_t n = std::min(size - 1, text.size());
        std::memcpy(buffer, text.data(), n);
        buffer[n] = '\0';
    }
}
// Every getter's value for NULL options is the setting's default: a default-made cps_options says what that is.
const cps_options &default_options() {
    static const cps_options defaults;
    return defaults;
}
const cps_options &or_defaults(const cps_options *o) { return o ? *o : default_options(); }
} // namespace

extern "C" {
uint32_t cps_abi_version(void) { return CPS_ABI_VERSION; }
const char *cps_version(void) { return CPS_VERSION_STRING; }
const char *cps_backend_version(void) { return flint_version; }
const char *cps_last_error(void) { return error_message; }
void cps_cleanup_thread(void) { flint_cleanup(); }

size_t cps_result_component_count(const cps_result *result) { return result ? result->components.size() : 0; }

cps_error cps_polynomial_create(cps_polynomial **out) {
    if (out)
        *out = nullptr;
    return guarded([&] {
        require(out, "null polynomial output");
        *out = new cps_polynomial;
    });
}
void cps_polynomial_destroy(cps_polynomial *p) { delete p; }
cps_error cps_secular_create(const cps_secular_term_view *terms, size_t count, const cps_resource_limits *limits,
                           cps_secular **out) {
    if (out)
        *out = nullptr;
    return guarded([&] {
        require(out && (count == 0 || terms), "null secular input or output");
        cps_resource_limits resources;
        cps_resource_limits_default(&resources);
        if (limits) resources = *limits;
        validate_resources(resources);
        if (count > static_cast<size_t>(std::numeric_limits<slong>::max()) ||
            count > std::vector<certpolysolve::detail::SecularTerm>().max_size())
            throw std::out_of_range("secular input exceeds addressable storage");
        auto source = std::make_unique<cps_secular>();
        using namespace certpolysolve::detail;
        const auto metadata = resources::multiply(count,
            sizeof(SecularTerm) + sizeof(size_t) + 4 * sizeof(Number));
        if (!resources::storage_fits(metadata, resources.temporary_bytes))
            throw std::out_of_range("secular input exceeds estimated temporary storage allowance");
        uint64_t remaining = std::min(resources.exact_materialization_bits,
            certpolysolve::detail::resources::multiply(resources.temporary_bytes - metadata, 8));
        if (count > resources.exact_operations / 4)
            throw std::out_of_range("secular normalization exceeds exact operation allowance");
        std::vector<std::array<Number, 4>> parsed;
        parsed.reserve(count);
        uint64_t owned = 0, expanded = 0;
        for (size_t i = 0; i < count; ++i) {
            const auto &v = terms[i];
            require(v.residue_real && v.pole_real, "null secular literal");
            // A NULL imaginary part is 0 (C19).
            const auto imaginary = [](const char *text) { return text ? Number(text) : Number(); };
            parsed.push_back({Number(v.residue_real), imaginary(v.residue_imag), Number(v.pole_real),
                              imaginary(v.pole_imag)});
            for (const auto &n : parsed.back()) {
                owned = certpolysolve::detail::resources::add(owned, n.owned_storage_bits());
                expanded = certpolysolve::detail::resources::add(expanded, n.expanded_storage_bits());
            }
        }
        if (!(resources.flags & CPS_RESOURCE_ALLOW_EXPONENT_EXPANSION) &&
            (expanded == UINT64_MAX || expanded > owned))
            throw std::out_of_range("secular normalization requires explicit exponent expansion");
        for (size_t i = 0; i < count; ++i) {
            certpolysolve::detail::SecularTerm term;
            fmpq *targets[] = {term.residue.real.value(), term.residue.imag.value(), term.pole.real.value(),
                               term.pole.imag.value()};
            for (size_t j = 0; j < 4; ++j) {
                if (!parsed[i][j].rational_with_budget(targets[j], remaining, true))
                    throw std::out_of_range("secular normalization needs a larger exact allowance or explicit exponent expansion");
            }
            require(!term.residue.zero(), "secular residues must be nonzero");
            source->terms.push_back(std::move(term));
        }
        // Distinct poles: an index sort (the terms keep their supplied order)
        // and a check of neighbours, n log n comparisons, not a pairwise scan.
        std::vector<size_t> order(count);
        std::iota(order.begin(), order.end(), size_t(0));
        std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
            const auto &x = source->terms[a].pole;
            const auto &y = source->terms[b].pole;
            const int real = fmpq_cmp(x.real.value(), y.real.value());
            return real ? real < 0 : fmpq_cmp(x.imag.value(), y.imag.value()) < 0;
        });
        for (size_t i = 1; i < count; ++i) {
            const auto &x = source->terms[order[i - 1]].pole;
            const auto &y = source->terms[order[i]].pole;
            require(!(fmpq_equal(x.real.value(), y.real.value()) && fmpq_equal(x.imag.value(), y.imag.value())),
                    "secular poles must be distinct");
        }
        *out = source.release();
    });
}
void cps_secular_destroy(cps_secular *p) { delete p; }
cps_error cps_solve_secular(const cps_secular *source, const cps_options *options, cps_result **out) {
    if (out)
        *out = nullptr;
    return guarded([&] {
        const auto started = std::chrono::steady_clock::now();
        const certpolysolve::detail::DefaultFloatingPoint environment;
        require(source && out, "null secular input or output");
        validate_query(or_defaults(options), EntryKind::Secular);
        TestingPolls polls;
        const cps_options timed = timed_options(options, polls);
        const certpolysolve::detail::TeamScope team(call_workers(timed));
        certpolysolve::detail::SessionState state(*source);
        auto result = std::make_unique<cps_result>(state.advance(&timed));
        finish_call(*result, timed, polls, started);
        *out = result.release();
    });
}
cps_error cps_polynomial_set_ball_coefficient(cps_polynomial *p, uint64_t index, const char *real,
                                            const char *imag, const char *real_radius,
                                            const char *imag_radius) {
    return guarded([&] {
        require(p && real, "null polynomial or coefficient");
        valid_exponent(index);
        using certpolysolve::detail::Number;
        // A NULL imaginary part or radius is 0 (C19).
        const auto part = [](const char *text) { return text ? Number(text) : Number(); };
        certpolysolve::detail::Coefficient c{Number(real), part(imag), part(real_radius), part(imag_radius)};
        require(c.real_radius.sign() >= 0 && c.imag_radius.sign() >= 0,
                "coefficient radii must be nonnegative");
        store(*p, index, std::move(c));
    });
}
cps_error cps_polynomial_set_coefficient(cps_polynomial *p, uint64_t index, const char *real, const char *imag) {
    return cps_polynomial_set_ball_coefficient(p, index, real, imag, "0", "0");
}
cps_error cps_polynomial_set_coefficient_d(cps_polynomial *p, uint64_t index, double real, double imag) {
    return guarded([&] {
        require(p, "null polynomial");
        valid_exponent(index);
        certpolysolve::detail::Coefficient c;
        c.real = certpolysolve::detail::Number::from_double(real);
        c.imag = certpolysolve::detail::Number::from_double(imag);
        store(*p, index, std::move(c));
    });
}
cps_error cps_polynomial_copy(const cps_polynomial *p, cps_polynomial **out) {
    return polynomial_result(out, [&] {
        require(p, "null polynomial");
        return *p;
    });
}
cps_error cps_polynomial_add(const cps_polynomial *a, const cps_polynomial *b, cps_polynomial **out) {
    return polynomial_result(out, [&] {
        require(a && b, "null polynomial");
        return certpolysolve::detail::polynomial_add(*a, *b);
    });
}
cps_error cps_polynomial_subtract(const cps_polynomial *a, const cps_polynomial *b, cps_polynomial **out) {
    return polynomial_result(out, [&] {
        require(a && b, "null polynomial");
        return certpolysolve::detail::polynomial_subtract(*a, *b);
    });
}
cps_error cps_polynomial_multiply(const cps_polynomial *a, const cps_polynomial *b, cps_polynomial **out) {
    return polynomial_result(out, [&] {
        require(a && b, "null polynomial");
        return certpolysolve::detail::polynomial_multiply(*a, *b);
    });
}
cps_error cps_polynomial_scale(const cps_polynomial *a, const char *real, const char *imag, cps_polynomial **out) {
    return polynomial_result(out, [&] {
        require(a && real, "null polynomial or number");
        using certpolysolve::detail::Number;
        const certpolysolve::detail::Coefficient factor{Number(real), imag ? Number(imag) : Number(), Number(), Number()};
        return certpolysolve::detail::polynomial_scale(*a, factor);
    });
}
cps_error cps_polynomial_power(const cps_polynomial *a, uint64_t exponent, cps_polynomial **out) {
    return polynomial_result(out, [&] {
        require(a, "null polynomial");
        return certpolysolve::detail::polynomial_power(*a, exponent);
    });
}
cps_error cps_polynomial_derivative(const cps_polynomial *a, cps_polynomial **out) {
    return polynomial_result(out, [&] {
        require(a, "null polynomial");
        return certpolysolve::detail::polynomial_derivative(*a);
    });
}
cps_error cps_polynomial_from_roots(const char *const *real, const char *const *imag, const uint64_t *multiplicity,
                                  size_t count, cps_polynomial **out) {
    return polynomial_result(out, [&] {
        require(!count || real, "null root list");
        using certpolysolve::detail::Number;
        std::vector<certpolysolve::detail::Coefficient> roots(count);
        std::vector<uint64_t> multiplicities(count, 1);
        for (size_t k = 0; k < count; ++k) {
            require(real[k], "null root");
            roots[k].real = Number(real[k]);
            if (imag && imag[k]) roots[k].imag = Number(imag[k]);
            if (multiplicity) multiplicities[k] = multiplicity[k];
        }
        return certpolysolve::detail::polynomial_from_roots(roots, multiplicities);
    });
}
int64_t cps_polynomial_degree(const cps_polynomial *p) {
    if (!p || p->terms.empty())
        return -1;
    return static_cast<int64_t>(p->terms.rbegin()->first); // at most max_exponent
}
cps_error cps_polynomial_get_coefficient(const cps_polynomial *p, uint64_t index, cps_coefficient_part part, char *buffer,
                                       size_t size, size_t *length) {
    return guarded([&] {
        require(p, "null polynomial");
        const auto which = enum_argument(part); // CPS_PART_REAL = 0 to CPS_PART_IMAG_RADIUS
        require(which >= CPS_PART_REAL && which <= CPS_PART_IMAG_RADIUS, "unknown coefficient part");
        const certpolysolve::detail::Number zero;
        const certpolysolve::detail::Number *n = &zero;
        if (const auto term = p->terms.find(index); term != p->terms.end()) {
            const auto &c = term->second;
            n = which == CPS_PART_REAL ? &c.real : which == CPS_PART_IMAG ? &c.imag
                : which == CPS_PART_REAL_RADIUS ? &c.real_radius : &c.imag_radius;
        }
        write_text(certpolysolve::detail::number_text(*n), buffer, size, length);
    });
}
cps_error cps_polynomial_terms(const cps_polynomial *p, uint64_t *exponents, size_t capacity, size_t *count) {
    return guarded([&] {
        require(p, "null polynomial");
        require(exponents || !capacity, "null exponent array with a nonzero capacity");
        if (count)
            *count = p->terms.size();
        size_t k = 0;
        for (auto term = p->terms.begin(); term != p->terms.end() && k < capacity; ++term)
            exponents[k++] = term->first;
    });
}
cps_error cps_dyadic_compare(cps_dyadic_view a, cps_dyadic_view b, int *out) {
    return guarded([&] {
        require(out, "null comparison output");
        certpolysolve::detail::Float x, y;
        certpolysolve::detail::parse_dyadic(x.value(), a);
        certpolysolve::detail::parse_dyadic(y.value(), b);
        const int c = arf_cmp(x.value(), y.value());
        *out = (c > 0) - (c < 0);
    });
}
cps_error cps_dyadic_to_double(cps_dyadic_view x, cps_rounding rounding, double *out) {
    return guarded([&] {
        require(out, "null double output");
        const auto mode = enum_argument(rounding); // CPS_ROUND_DOWN = 0 to CPS_ROUND_UP
        require(mode >= CPS_ROUND_DOWN && mode <= CPS_ROUND_UP, "unknown rounding");
        certpolysolve::detail::Float v;
        certpolysolve::detail::parse_dyadic(v.value(), x);
        *out = arf_get_d(v.value(), mode == CPS_ROUND_DOWN ? ARF_RND_FLOOR
                                    : mode == CPS_ROUND_UP ? ARF_RND_CEIL
                                                          : ARF_RND_NEAR);
    });
}
cps_error cps_dyadic_format(cps_dyadic_view x, unsigned digits, char *buffer, size_t size, size_t *length) {
    return guarded([&] {
        require(digits >= 1, "at least one significant digit");
        certpolysolve::detail::Float v;
        certpolysolve::detail::parse_dyadic(v.value(), x);
        write_text(certpolysolve::detail::dyadic_format(v.value(), digits), buffer, size, length);
    });
}
cps_error cps_dyadic_rational(cps_dyadic_view x, char *buffer, size_t size, size_t *length) {
    return guarded([&] {
        certpolysolve::detail::Float v;
        certpolysolve::detail::parse_dyadic(v.value(), x);
        write_text(certpolysolve::detail::dyadic_rational(v.value()), buffer, size, length);
    });
}
cps_error cps_dyadic_midpoint(cps_dyadic_view a, cps_dyadic_view b, char *buffer, size_t size, size_t *length) {
    return guarded([&] {
        certpolysolve::detail::Float x, y, m;
        certpolysolve::detail::parse_dyadic(x.value(), a);
        certpolysolve::detail::parse_dyadic(y.value(), b);
        arf_add(m.value(), x.value(), y.value(), ARF_PREC_EXACT, ARF_RND_DOWN);
        arf_mul_2exp_si(m.value(), m.value(), -1);
        write_text(certpolysolve::detail::dyadic_text(m.value()), buffer, size, length);
    });
}
cps_error cps_dyadic_width(cps_dyadic_view a, cps_dyadic_view b, char *buffer, size_t size, size_t *length) {
    return guarded([&] {
        certpolysolve::detail::Float x, y, w;
        certpolysolve::detail::parse_dyadic(x.value(), a);
        certpolysolve::detail::parse_dyadic(y.value(), b);
        arf_sub(w.value(), y.value(), x.value(), ARF_PREC_EXACT, ARF_RND_DOWN);
        write_text(certpolysolve::detail::dyadic_text(w.value()), buffer, size, length);
    });
}
const char *cps_termination_name(cps_termination t) {
    switch (enum_argument(t)) { // any value a C caller passes; one outside the enumeration is "unknown"
    case CPS_SUCCESS: return "success";
    case CPS_PRECISION_LIMIT: return "precision limit";
    case CPS_CANCELLED: return "cancelled";
    case CPS_ZERO_POLYNOMIAL: return "zero polynomial";
    case CPS_INDETERMINATE_DEGREE: return "indeterminate degree";
    case CPS_UNSUPPORTED_QUERY: return "unsupported query";
    case CPS_WORK_LIMIT: return "work limit";
    }
    return "unknown";
}
const char *cps_stop_reason_name(cps_stop_reason s) {
    switch (enum_argument(s)) {
    case CPS_STOP_NONE: return "none";
    case CPS_STOP_WORK_BUDGET: return "work budget";
    case CPS_STOP_CELL_BUDGET: return "cell budget";
    case CPS_STOP_PRECISION: return "precision cap";
    case CPS_STOP_CANCELLED: return "cancelled";
    case CPS_STOP_INPUT: return "input";
    case CPS_STOP_ITERATION_BUDGET: return "iteration budget";
    case CPS_STOP_TIME_LIMIT: return "time limit";
    }
    return "unknown";
}
cps_error cps_program_create(cps_program **out) {
    if (out)
        *out = nullptr;
    return guarded([&] {
        require(out, "null program output");
        *out = new cps_program;
    });
}
void cps_program_destroy(cps_program *p) { delete p; }
cps_error cps_program_variable(cps_program *p, uint64_t *out) {
    return program_node(p, out, program_op(certpolysolve::detail::ProgramOp::Variable));
}
cps_error cps_program_constant(cps_program *p, const char *real, const char *imag, uint64_t *out) {
    if (out)
        *out = CPS_NO_NODE;
    return guarded([&] {
        require(p && out && real, "null program, node output or constant");
        certpolysolve::detail::Coefficient c;
        c.real = certpolysolve::detail::Number(real);
        c.imag = certpolysolve::detail::Number(imag ? imag : "0");
        // append_constant adds the constant and its node together; a refusal adds neither.
        const auto constants = p->program.constants.size();
        try {
            *out = p->program.append_constant(std::move(c));
        } catch (...) {
            p->program.constants.resize(constants);
            throw;
        }
    });
}
cps_error cps_program_add(cps_program *p, uint64_t a, uint64_t b, uint64_t *out) {
    return program_node(p, out, program_op(certpolysolve::detail::ProgramOp::Add, a, b));
}
cps_error cps_program_subtract(cps_program *p, uint64_t a, uint64_t b, uint64_t *out) {
    return program_node(p, out, program_op(certpolysolve::detail::ProgramOp::Subtract, a, b));
}
cps_error cps_program_multiply(cps_program *p, uint64_t a, uint64_t b, uint64_t *out) {
    return program_node(p, out, program_op(certpolysolve::detail::ProgramOp::Multiply, a, b));
}
cps_error cps_program_negate(cps_program *p, uint64_t a, uint64_t *out) {
    return program_node(p, out, program_op(certpolysolve::detail::ProgramOp::Negate, a));
}
cps_error cps_program_power(cps_program *p, uint64_t a, uint64_t exponent, uint64_t *out) {
    return program_node(p, out, program_op(certpolysolve::detail::ProgramOp::Power, a, 0, exponent));
}
cps_error cps_program_set_output(cps_program *p, uint64_t node) {
    return guarded([&] {
        require(p, "null program");
        require(node < p->program.nodes.size(), "program output is not a node");
        p->program.output = node;
    });
}
cps_error cps_program_parse(const char *text, cps_program **out) {
    if (out)
        *out = nullptr;
    return guarded([&] {
        require(text && out, "null program text or output");
        auto parsed = std::make_unique<cps_program>();
        parsed->program = certpolysolve::detail::parse_program(text);
        *out = parsed.release();
    });
}
uint64_t cps_program_size(const cps_program *p) { return p ? static_cast<uint64_t>(p->program.nodes.size()) : 0; }
uint64_t cps_program_output(const cps_program *p) {
    return p && p->program.output ? *p->program.output : CPS_NO_NODE;
}
cps_error cps_solve_program(const cps_program *program, const cps_options *options, cps_result **out) {
    if (out)
        *out = nullptr;
    return guarded([&] {
        const auto started = std::chrono::steady_clock::now();
        const certpolysolve::detail::DefaultFloatingPoint environment;
        require(program && out, "null program or result output");
        validate_query(or_defaults(options), EntryKind::Program);
        TestingPolls polls;
        const cps_options timed = timed_options(options, polls);
        const certpolysolve::detail::TeamScope team(call_workers(timed));
        certpolysolve::detail::SessionState state(*program);
        auto result = std::make_unique<cps_result>(state.advance(&timed));
        finish_call(*result, timed, polls, started);
        *out = result.release();
    });
}
cps_error cps_options_create(cps_options **out) {
    if (out)
        *out = nullptr;
    return guarded([&] {
        require(out, "null options output");
        *out = new cps_options;
    });
}
void cps_options_destroy(cps_options *o) { delete o; }
void cps_resource_limits_default(cps_resource_limits *out) {
    if (out) *out = {UINT64_MAX, UINT64_MAX, UINT64_MAX, 0};
}
cps_error cps_options_set_resource_limits(cps_options *o, const cps_resource_limits *r) {
    return guarded([&] {
        require(o && r, "null options or resource settings");
        validate_resources(*r);
        o->exact_budget = r->exact_materialization_bits;
        o->exact_operations = r->exact_operations;
        o->temporary_bytes = r->temporary_bytes;
        o->allow_exponent_expansion = r->flags & CPS_RESOURCE_ALLOW_EXPONENT_EXPANSION;
    });
}
cps_error cps_options_get_resource_limits(const cps_options *o, cps_resource_limits *r) {
    return guarded([&] {
        require(o && r, "null options or resource output");
        *r = {o->exact_budget, o->exact_operations, o->temporary_bytes,
              o->allow_exponent_expansion ? CPS_RESOURCE_ALLOW_EXPONENT_EXPANSION : UINT64_C(0)};
    });
}
cps_error cps_options_set_precision(cps_options *o, uint64_t bits) {
    return guarded([&] {
        require(o && bits >= 1, "requested bits must be positive");
        if (bits > certpolysolve::detail::resources::precision_max)
            throw std::out_of_range("precision exceeds backend representability");
        o->bits = bits;
    });
}
cps_error cps_options_set_method(cps_options *o, int method) {
    return guarded([&] {
        require(o && method >= CPS_METHOD_AUTO && method <= CPS_METHOD_PIECEWISE, "invalid solver method");
        o->method = static_cast<cps_method>(method);
    });
}
cps_error cps_options_set_absolute_tolerance(cps_options *o, const char *tolerance) {
    return guarded([&] {
        require(o && tolerance, "null options or tolerance");
        certpolysolve::detail::Number value(tolerance);
        require(value.sign() >= 0, "absolute tolerance must be nonnegative");
        o->absolute_tolerance_text = trimmed(tolerance);
        o->absolute_tolerance = std::move(value);
    });
}
cps_error cps_options_set_goal(cps_options *o, int goal) {
    return guarded([&] {
        require(o && (goal == CPS_ROOT_ACCURACY || goal == CPS_ROOT_APPROXIMATION || goal == CPS_CLUSTER_RADIUS ||
                      goal == CPS_ROOT_ISOLATION || goal == CPS_GOAL_DEFAULT),
                "invalid options or goal");
        require(goal != CPS_CLUSTER_RADIUS || o->cluster_radius.sign() > 0,
                "cluster-radius goal requires a positive cluster radius");
        if (goal != CPS_CLUSTER_RADIUS) {
            o->cluster_radius = certpolysolve::detail::Number();
            o->cluster_radius_text = "0";
        }
        o->goal = static_cast<cps_goal>(goal);
    });
}
cps_error cps_options_set_cluster_radius(cps_options *o, const char *radius) {
    return guarded([&] {
        require(o && radius, "null options or cluster radius");
        certpolysolve::detail::Number value(radius);
        require(value.sign() >= 0, "cluster radius must be nonnegative");
        o->cluster_radius_text = value.zero() ? "0" : trimmed(radius);
        o->cluster_radius = std::move(value);
        o->goal = o->cluster_radius.zero() ? CPS_GOAL_DEFAULT : CPS_CLUSTER_RADIUS;
    });
}
cps_error cps_options_set_max_iterations(cps_options *o, uint64_t sweeps) {
    return guarded([&] {
        require(o, "null options");
        o->max_iterations = sweeps;
    });
}
cps_error cps_options_set_workers(cps_options *o, uint32_t workers) {
    return guarded([&] {
        require(o, "null options");
        require(workers >= 1, "a solve needs at least one worker");
        o->workers = workers;
    });
}
cps_error cps_options_set_region(cps_options *o, const char *real_lower, const char *real_upper,
                               const char *imag_lower, const char *imag_upper) {
    return guarded([&] {
        require(o, "null options");
        // Each end's literal only (C8): the ends' order and their exact form are the solve's checks.
        const auto end = [](const char *text) -> std::optional<certpolysolve::detail::Number> {
            if (!text || !*text) return std::nullopt;
            return certpolysolve::detail::Number(text);
        };
        cps_options::Region region{end(real_lower), end(real_upper), end(imag_lower), end(imag_upper)};
        const auto text = [](const char *end) { return end ? trimmed(end) : std::string(); };
        std::array<std::string, 4> texts{text(real_lower), text(real_upper), text(imag_lower), text(imag_upper)};
        o->region = std::move(region);
        o->region_text = o->region.set() ? std::move(texts) : std::array<std::string, 4>{};
    });
}
cps_error cps_options_set_real_only(cps_options *o, int enabled) {
    return guarded([&] {
        require(o, "null options");
        o->real_only = enabled != 0;
    });
}
cps_error cps_cancel_token_create(cps_cancel_token **out) {
    if (out)
        *out = nullptr;
    return guarded([&] {
        require(out, "null token output");
        *out = new cps_cancel_token();
    });
}
void cps_cancel_token_request(cps_cancel_token *token) {
    if (token) token->flag->store(true, std::memory_order_relaxed);
}
int cps_cancel_token_requested(const cps_cancel_token *token) {
    return token && token->flag->load(std::memory_order_relaxed) ? 1 : 0;
}
void cps_cancel_token_destroy(cps_cancel_token *token) { delete token; }
cps_error cps_options_set_cancel_token(cps_options *o, cps_cancel_token *token) {
    return guarded([&] {
        require(o, "null options");
        o->cancel_flag = token ? token->flag : nullptr;
    });
}
cps_error cps_options_set_time_limit(cps_options *o, double seconds) {
    return guarded([&] {
        require(o, "null options");
        require(!std::isnan(seconds) && seconds >= 0, "the time limit must be 0 (none), +infinity (none) or positive");
        o->time_limit = std::isinf(seconds) ? 0 : seconds;
    });
}
cps_error cps_solve(const cps_polynomial *p, const cps_options *o, cps_result **out) {
    if (out)
        *out = nullptr;
    return guarded([&] {
        const auto started = std::chrono::steady_clock::now();
        const certpolysolve::detail::DefaultFloatingPoint environment;
        require(p && out, "null polynomial or result output");
        validate_query(or_defaults(o), EntryKind::Coefficients);
        TestingPolls polls;
        const cps_options timed = timed_options(o, polls);
        const certpolysolve::detail::TeamScope team(call_workers(timed));
        auto result = std::make_unique<cps_result>(certpolysolve::detail::solve(*p, timed));
        finish_call(*result, timed, polls, started);
        *out = result.release();
    });
}
void cps_result_destroy(cps_result *r) { delete r; }
int64_t cps_result_degree(const cps_result *r) { return r ? r->degree : -1; }
cps_termination cps_result_termination(const cps_result *r) { return r ? r->termination : CPS_UNSUPPORTED_QUERY; }
int cps_result_complete(const cps_result *r) { return r && r->complete; }
int cps_result_goal_met(const cps_result *r) { return r && r->goal_met; }
cps_goal cps_result_goal(const cps_result *r) { return r ? r->goal : CPS_ROOT_ACCURACY; }
uint64_t cps_result_work_bits(const cps_result *r) { return r ? r->work_bits : 0; }
size_t cps_result_size(const cps_result *r) { return r ? r->roots.size() : 0; }
const char *cps_result_diagnostics(const cps_result *r) { return r ? r->diagnostics.c_str() : "null result"; }
cps_error cps_result_get_root(const cps_result *r, size_t index, cps_root_view *out) {
    return guarded([&] {
        require(r && out, "null result or root output");
        if (index >= r->roots.size())
            throw std::out_of_range("root index out of range");
        require(r->exported.size() == r->roots.size(), "a result's regions were not published as text");
        const auto &root = r->roots[index];
        const auto &text = r->exported[index];
        *out = {text[0].view(),
                text[1].view(),
                text[2].view(),
                text[3].view(),
                root.count,
                root.multiplicity,
                root.realness,
                root.certificate.c_str()};
    });
}
cps_error cps_result_get_approximation(const cps_result *r, size_t index, double *real, double *imag) {
    error_message[0] = '\0'; // as guarded() does: a call that succeeds leaves the thread's message empty (C2)
    if (!r || !real || !imag) {
        error_text("null result or approximation output");
        return CPS_INVALID_ARGUMENT;
    }
    if (index >= r->roots.size()) {
        error_text("root index out of range");
        return CPS_OUT_OF_RANGE;
    }
    const auto &root = r->roots[index];
    if (!r->complete || !r->goal_met || root.multiplicity != root.count) {
        error_text("accuracy unmet or cluster unresolved; inspect certified enclosures");
        return CPS_UNRESOLVED;
    }
    // An isolating interval may be as wide as the root's distance to its neighbours: its midpoint is no
    // approximation of the root.
    if (r->goal == CPS_ROOT_ISOLATION) {
        error_text("an isolation answer promises no width; inspect the exact bounds or tighten to an accuracy goal");
        return CPS_UNRESOLVED;
    }
    const auto re = arb_midref(acb_realref(root.box.value()));
    const auto im = arb_midref(acb_imagref(root.box.value()));
    const double re_double = arf_get_d(re, ARF_RND_NEAR), im_double = arf_get_d(im, ARF_RND_NEAR);
    if (!std::isfinite(re_double) || !std::isfinite(im_double) || (re_double == 0 && !arf_is_zero(re)) ||
        (im_double == 0 && !arf_is_zero(im))) {
        error_text("binary64 conversion would overflow or underflow; use exact dyadic bounds");
        return CPS_OUT_OF_RANGE;
    }
    *real = re_double;
    *imag = im_double;
    return CPS_OK;
}

// The 0.5 calls on regions, options and results (docs/api-spec-0.5.md).
cps_error cps_options_get_region(const cps_options *o, const char **real_lower, const char **real_upper,
                               const char **imag_lower, const char **imag_upper) {
    return guarded([&] {
        require(o && real_lower && real_upper && imag_lower && imag_upper, "null options or region output");
        const auto end = [&](size_t k) { return o->region_text[k].empty() ? nullptr : o->region_text[k].c_str(); };
        *real_lower = end(0), *real_upper = end(1), *imag_lower = end(2), *imag_upper = end(3);
    });
}
cps_error cps_options_copy(const cps_options *source, cps_options **out) {
    if (out)
        *out = nullptr;
    return guarded([&] {
        require(source && out, "null options or copy output");
        auto copy = std::make_unique<cps_options>(*source);
        // A call's state, never a caller's setting: the call that owns it fixes it on its own copy.
        copy->deadline.reset();
        *out = copy.release();
    });
}
cps_method cps_options_method(const cps_options *o) { return or_defaults(o).method; }
uint64_t cps_options_precision(const cps_options *o) { return or_defaults(o).bits; }
uint64_t cps_options_max_work_bits(const cps_options *o) { return or_defaults(o).max_work_bits; }
cps_goal cps_options_goal(const cps_options *o) { return or_defaults(o).goal; }
const char *cps_options_cluster_radius(const cps_options *o) { return or_defaults(o).cluster_radius_text.c_str(); }
const char *cps_options_absolute_tolerance(const cps_options *o) {
    return or_defaults(o).absolute_tolerance_text.c_str();
}
int cps_options_real_only(const cps_options *o) { return or_defaults(o).real_only ? 1 : 0; }
uint64_t cps_options_work_budget(const cps_options *o) {
    const auto &v = or_defaults(o);
    return v.work_units_default ? 0 : v.work_units;
}
uint64_t cps_options_live_cells(const cps_options *o) { return or_defaults(o).live_cells; }
uint64_t cps_options_max_iterations(const cps_options *o) { return or_defaults(o).max_iterations; }
double cps_options_time_limit(const cps_options *o) { return or_defaults(o).time_limit; }
uint32_t cps_options_workers(const cps_options *o) { return or_defaults(o).workers; }
cps_error cps_options_set_max_work_bits(cps_options *o, uint64_t bits) {
    return guarded([&] {
        require(o && bits >= CPS_MIN_MAX_WORK_BITS, "the maximum working precision is below CPS_MIN_MAX_WORK_BITS");
        if (bits > certpolysolve::detail::resources::precision_max)
            throw std::out_of_range("precision exceeds backend representability");
        o->max_work_bits = bits;
    });
}
cps_error cps_options_set_work_budget(cps_options *o, uint64_t units) {
    return guarded([&] {
        require(o, "null options");
        o->work_units_default = units == 0;
        o->work_units = units ? units : default_options().work_units;
    });
}
cps_error cps_options_set_live_cells(cps_options *o, uint64_t cells) {
    return guarded([&] {
        require(o, "null options");
        o->live_cells = cells ? cells : default_options().live_cells;
    });
}
cps_error cps_result_get_statistics(const cps_result *r, cps_statistics *out) {
    using certpolysolve::detail::resources::add;
    return guarded([&] {
        require(r, "null result");
        cps_statistics s{};
        s.work_units = r->statistics_v2.work_units;
        s.sweeps = r->statistics_v3.candidate_sweeps;
        s.precision_increases = r->statistics_v2.precision_increases;
        s.final_work_bits = r->work_bits;
        s.highest_work_bits = std::max(r->progress.max_outer_bits, r->work_bits);
        s.regions = r->roots.size();
        for (const auto &root : r->roots) {
            s.roots += root.count;
            s.resolved_regions += root.multiplicity != 0;
        }
        s.exact_factors = r->statistics.exact_factors;
        s.total_ns = r->total_ns;
        s.preparation_ns = add(r->statistics.preprocessing_ns, r->carried_v1.preprocessing_ns);
        fill_sized(out, s);
    });
}
cps_error cps_result_get_engine_counters(const cps_result *r, cps_engine_counters *out) {
    using certpolysolve::detail::resources::add;
    return guarded([&] {
        require(r, "null result");
        // The v1 counters with those of the engines that handed over (the old v1 getter has the last engine's).
        certpolysolve::detail::StageStatistics a = r->statistics;
        certpolysolve::detail::statistics::add_saturating(a, r->carried_v1, certpolysolve::detail::statistics::carried_stage_fields);
        const auto &b = r->statistics_v2;
        const auto &c = r->statistics_v3;
        const auto &p = r->progress;
        cps_engine_counters e{};
        e.input_enclosure_ns = a.input_enclosure_ns, e.isolation_ns = a.isolation_ns;
        e.coverage_validation_ns = a.coverage_validation_ns, e.backend_ns = b.backend_ns;
        e.taylor_shift_ns = b.taylor_shift_ns, e.counting_ns = b.counting_ns, e.refinement_ns = b.refinement_ns;
        e.search_ns = b.search_ns, e.candidate_ns = c.candidate_ns, e.evaluation_ns = c.evaluation_ns;
        e.interaction_ns = c.interaction_ns, e.validation_ns = c.validation_ns, e.regeneration_ns = c.regeneration_ns;
        e.precision_attempts = a.precision_attempts, e.isolation_passes = a.isolation_passes;
        e.disk_tests = b.disk_tests, e.graeffe_steps = b.graeffe_steps, e.taylor_shifts = b.taylor_shifts;
        e.taylor_cache_hits = b.taylor_cache_hits, e.boxes_split = b.boxes_split, e.boxes_excluded = b.boxes_excluded;
        e.newton_attempts = b.newton_attempts, e.newton_successes = b.newton_successes;
        e.backend_calls = b.backend_calls, e.evaluations = c.evaluations, e.evaluation_terms = c.evaluation_terms;
        e.root_interactions = c.root_interactions, e.validation_calls = c.validation_calls;
        e.secular_regenerations = c.secular_regenerations, e.double_seed_calls = c.double_seed_calls;
        e.candidate_stagnations = c.candidate_stagnations, e.real_refinement_steps = c.real_refinement_steps;
        e.real_bisections = c.real_bisections;
        e.pending_components = p.pending_components, e.live_cells = p.live_cells;
        e.min_pending_bits = p.min_pending_bits, e.max_pending_bits = p.max_pending_bits;
        e.input_uncertainty = p.input_uncertainty;
        fill_sized(out, e);
    });
}
cps_error cps_result_get_component(const cps_result *r, size_t index, cps_component_progress *out) {
    return guarded([&] {
        require(r && out, "null result or component output");
        if (index >= r->components.size())
            throw std::out_of_range("component index out of range");
        const auto &c = r->components[index];
        *out = {c.factor_index, c.work_bits, c.live_cells, c.certified_count, c.state};
    });
}
cps_stop_reason cps_result_stop_reason(const cps_result *r) { return r ? r->progress.stop_reason : CPS_STOP_INPUT; }
int cps_result_preparation_limited(const cps_result *r) { return r && r->progress.preprocessing_limited ? 1 : 0; }
uint64_t cps_result_sufficient_work_bits(const cps_result *r) { return r ? r->sufficient_work_bits : 0; }
cps_squarefree cps_result_squarefree(const cps_result *r) { return r ? r->squarefree : CPS_SQUAREFREE_UNKNOWN; }
cps_handoff cps_result_handoff(const cps_result *r) { return r ? r->handoff : CPS_HANDOFF_NONE; }
// certpolysolve_testing.h.
cps_error cps_testing_options_stop_at_poll(cps_options *o, uint64_t k) {
    return guarded([&] {
        require(o, "null options");
        o->testing_stop_at_poll = k;
    });
}
uint64_t cps_testing_result_polls(const cps_result *r) { return r ? r->testing_polls : 0; }
static_assert(certpolysolve::detail::split_sites == CPS_TESTING_SPLIT_SITES, "the testing header numbers every split site");
uint64_t cps_testing_result_split_jobs(const cps_result *r, unsigned site) {
    return r && site < certpolysolve::detail::split_sites ? r->testing_split_jobs[site] : 0;
}
} // extern "C"
