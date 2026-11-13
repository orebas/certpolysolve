#pragma once
// CertPolySolve's C++ interface (the guide: docs/api.md; a first tour: docs/quickstart-cpp.md).
//
// The types are values, and each copies. Exact ends, intervals, boxes and polynomials compare with ==; numbers,
// coefficients, polynomials, exact ends, intervals, roots and stops print with std::format. Exact numbers go in
// (integers, doubles at their exact binary value, text such as "1/3", "0.1" or "3@-2"), and exact root boxes come out. Nothing throws but
// std::bad_alloc: an input remembers its first bad value (ok(), status()), and a solve reports it in its Result,
// as it reports every stop (ok(), stop(), message()). Braces and parentheses mean the same thing for every type:
// Polynomial{3} and Polynomial(3) are both the constant 3.
//
// C++ traps this header cannot remove (docs/quickstart-cpp.md, "C++ pitfalls"): 1/3 in code is integer division
// (0), so write "1/3"; 0.1 in code is a binary double, not one tenth, so write "0.1" for the exact decimal; and
// `Polynomial p(Number(x));` declares a function, so write `auto p = Polynomial(x);` or use braces.
#include <certpolysolve/certpolysolve.h>
#include <charconv>
#include <chrono>
#include <cmath>
#include <compare>
#include <complex>
#include <concepts>
#include <cstdint>
#include <format>
#include <initializer_list>
#include <memory>
#include <new>
#include <optional>
#include <ranges>
#include <stop_token>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

namespace certpolysolve {
inline const char *version() noexcept { return cps_version(); }
inline uint32_t abi_version() noexcept { return cps_abi_version(); }
inline const char *backend_version() noexcept { return cps_backend_version(); }

// Construct at a worker thread's entry: its destructor frees the thread's FLINT caches (cps_cleanup_thread; a
// short-lived thread that skips it leaks them). Objects made on the thread may outlive it and move between threads.
class ThreadScope {
  public:
    ThreadScope() = default;
    ThreadScope(const ThreadScope &) = delete;
    ThreadScope &operator=(const ThreadScope &) = delete;
    ~ThreadScope() { cps_cleanup_thread(); }
};

// Specialize to let a type of your own (or a library's) be a Number: `static std::string text(const T &)` returns
// its exact value in the coefficient grammar. <certpolysolve/gmp.hpp> does this for mpz_class and mpq_class.
template <class T> struct number_traits;

namespace detail {
template <class T> concept character = std::same_as<std::remove_cv_t<T>, char> || std::same_as<std::remove_cv_t<T>, wchar_t> ||
    std::same_as<std::remove_cv_t<T>, char8_t> || std::same_as<std::remove_cv_t<T>, char16_t> ||
    std::same_as<std::remove_cv_t<T>, char32_t> || std::same_as<std::remove_cv_t<T>, bool>;
template <class T> concept traited = requires(const T &v) {
    { ::certpolysolve::number_traits<std::remove_cvref_t<T>>::text(v) } -> std::convertible_to<std::string>;
};
inline void memory(cps_error e) {
    if (e == CPS_OUT_OF_MEMORY) throw std::bad_alloc();
}
// The exact binary value of a double in the coefficient grammar: "m@e" (or text the library rejects for NaN and
// the infinities, so the value is remembered as bad).
inline std::string double_text(double v) {
    if (std::isnan(v)) return "nan";
    if (std::isinf(v)) return v > 0 ? "inf" : "-inf";
    if (v == 0) return "0";
    int e = 0;
    const double f = std::frexp(v, &e); // v = f 2^e, 1/2 <= |f| < 1: f 2^53 is an exact integer
    auto m = static_cast<long long>(std::ldexp(f, 53));
    e -= 53;
    while (m % 2 == 0) m /= 2, ++e; // an odd mantissa: one text for one value (0.5 is "1@-1")
    return std::to_string(m) + "@" + std::to_string(e);
}
// A C call that writes text (buffer, size, length): the whole text.
template <class Call> std::string text(Call &&call) {
    size_t length = 0;
    memory(call(nullptr, 0, &length));
    std::string s(length + 1, '\0');
    memory(call(s.data(), s.size(), &length));
    s.resize(length);
    return s;
}
} // namespace detail

// An exact number you give: an integer, a double (its exact binary value), text in the coefficient grammar
// (integer, "p/q", decimal or scientific, "m@e" for m 2^e), a std::complex<double>, or Number::complex(re, im).
// Characters and bools are refused at compile time ('3' would be 51). The library checks text when it is used.
class Number {
    std::string real_ = "0", imag_ = "0";

  public:
    Number() = default;
    template <std::integral T> requires(!detail::character<T>) Number(T v) : real_(std::to_string(v)) {}
    template <class T> requires detail::character<T> Number(T) = delete;
    Number(double v) : real_(detail::double_text(v)) {}
    Number(const char *text) : real_(text ? text : "") {}
    Number(std::string_view text) : real_(text) {}
    Number(const std::string &text) : real_(text) {}
    Number(std::complex<double> z) : real_(detail::double_text(z.real())), imag_(detail::double_text(z.imag())) {}
    template <detail::traited T> Number(const T &v) : real_(number_traits<std::remove_cvref_t<T>>::text(v)) {}
    // re + i im, from two real numbers.
    static Number complex(const Number &re, const Number &im) {
        Number n;
        n.real_ = re.real_;
        n.imag_ = im.real_;
        return n;
    }
    const std::string &real() const noexcept { return real_; }
    const std::string &imag() const noexcept { return imag_; }
    bool is_real() const noexcept { return imag_ == "0"; }
};

// A polynomial's coefficient: an exact number, or a ball (a number known within a radius; ball() makes one).
struct Coefficient {
    Number value;
    Number real_radius, imag_radius;
    Coefficient() = default;
    template <class T> requires std::constructible_from<Number, T> Coefficient(T &&v) : value(std::forward<T>(v)) {}
    bool exact() const noexcept { return real_radius.real() == "0" && imag_radius.real() == "0"; }
};
// value +/- radius (for a complex value, radius on each part), or with a radius for each part.
inline Coefficient ball(Number value, Number radius) {
    Coefficient c(std::move(value));
    c.real_radius = radius;
    if (!c.value.is_real()) c.imag_radius = std::move(radius);
    return c;
}
inline Coefficient ball(Number value, Number real_radius, Number imag_radius) {
    Coefficient c(std::move(value));
    c.real_radius = std::move(real_radius);
    c.imag_radius = std::move(imag_radius);
    return c;
}

enum class Rounding { down, nearest, up };

// An exact root end m 2^e: a value you keep, compare exactly (with another Dyadic, an integer or a double), round
// to a double, or write exactly ("p/q") or in decimal.
class Dyadic {
    std::string mantissa_ = "0", exponent_ = "0";

  public:
    Dyadic() = default;
    Dyadic(std::string mantissa, std::string exponent) : mantissa_(std::move(mantissa)), exponent_(std::move(exponent)) {}
    explicit Dyadic(cps_dyadic_view v) : mantissa_(v.mantissa ? v.mantissa : "0"), exponent_(v.exponent ? v.exponent : "0") {}
    template <std::integral T> requires(!detail::character<T>) Dyadic(T v) : mantissa_(std::to_string(v)) {}
    Dyadic(double v) {
        const auto t = detail::double_text(v);
        const auto at = t.find('@');
        if (at == std::string::npos) mantissa_ = t;
        else mantissa_ = t.substr(0, at), exponent_ = t.substr(at + 1);
    }
    // "m@e" or "m" (the text the library's midpoint and width give).
    static Dyadic parse(std::string_view t) {
        const auto at = t.find('@');
        return at == std::string_view::npos ? Dyadic(std::string(t), "0") : Dyadic(std::string(t.substr(0, at)), std::string(t.substr(at + 1)));
    }
    cps_dyadic_view view() const noexcept { return {mantissa_.c_str(), exponent_.c_str()}; }
    const std::string &mantissa() const noexcept { return mantissa_; }
    const std::string &exponent() const noexcept { return exponent_; }
    double to_double(Rounding r = Rounding::nearest) const {
        double d = 0;
        detail::memory(cps_dyadic_to_double(view(), static_cast<cps_rounding>(r), &d));
        return d;
    }
    // The exact value, "p/q" or "p" in lowest terms.
    std::string to_string() const {
        return detail::text([&](char *b, size_t n, size_t *l) { return cps_dyadic_rational(view(), b, n, l); });
    }
    // Decimal, at most `digits` significant digits, as std::format("{:.N}") writes a double (cps_dyadic_format).
    std::string format(unsigned digits) const {
        return detail::text([&](char *b, size_t n, size_t *l) { return cps_dyadic_format(view(), digits ? digits : 1, b, n, l); });
    }
    friend std::strong_ordering operator<=>(const Dyadic &a, const Dyadic &b) {
        int c = 0;
        detail::memory(cps_dyadic_compare(a.view(), b.view(), &c));
        return c < 0 ? std::strong_ordering::less : c > 0 ? std::strong_ordering::greater : std::strong_ordering::equal;
    }
    friend bool operator==(const Dyadic &a, const Dyadic &b) { return (a <=> b) == 0; }
};

// A closed interval [lo, hi] with exact ends; Box: its real and imaginary intervals.
struct Interval {
    Dyadic lo, hi;
    Dyadic mid() const {
        return Dyadic::parse(detail::text([&](char *b, size_t n, size_t *l) { return cps_dyadic_midpoint(lo.view(), hi.view(), b, n, l); }));
    }
    Dyadic width() const {
        return Dyadic::parse(detail::text([&](char *b, size_t n, size_t *l) { return cps_dyadic_width(lo.view(), hi.view(), b, n, l); }));
    }
    bool contains(const Dyadic &x) const { return lo <= x && x <= hi; }
    bool is_point() const { return lo == hi; }
    friend bool operator==(const Interval &, const Interval &) = default;
};
struct Box {
    Interval real, imag;
    friend bool operator==(const Box &, const Box &) = default;
};

enum class Realness { unknown, real, not_real };

// One region of an answer: its exact box, the number of roots in it (counted with multiplicity), the multiplicity
// when the region is proved to hold one root of it, whether the roots are proved real or not real, and its value as
// a double where the answer gives one. Printing gives the centre: {} to double precision, {:.N} with N digits from
// the exact box.
struct Root {
    Box box;
    // The root as a double, exactly the number cps_result_get_approximation gives for this region, and only where it
    // gives one: the result complete with its goal met, the region one root of proved multiplicity, a goal that asks
    // for a width (not isolation), and a centre inside binary64's range. Empty otherwise: after a stop, for an
    // unresolved cluster, and for an isolation answer. C and C++ never give different values.
    std::optional<std::complex<double>> value;
    uint64_t count = 1;
    std::optional<uint64_t> multiplicity;
    Realness realness = Realness::unknown;
    std::string certificate;
    // The box's centre rounded to doubles, for any region: what printing shows. Not a root's value unless value is
    // set (then the two are equal); a cluster's centre, or overflowing to infinity, is possible.
    std::complex<double> centre() const { return {box.real.mid().to_double(), box.imag.mid().to_double()}; }
};

// Why a solve stopped; Stop::none with ok() when its goal is met.
enum class Stop {
    none, precision_cap, input_uncertainty, work_budget, cell_budget, iteration_budget, time_limit, cancelled,
    unsupported, zero_polynomial, indeterminate_degree, error
};
inline const char *name(Stop s) noexcept {
    switch (s) {
    case Stop::none: return "none";
    case Stop::precision_cap: return "precision cap";
    case Stop::input_uncertainty: return "input uncertainty";
    case Stop::work_budget: return "work budget";
    case Stop::cell_budget: return "cell budget";
    case Stop::iteration_budget: return "iteration budget";
    case Stop::time_limit: return "time limit";
    case Stop::cancelled: return "cancelled";
    case Stop::unsupported: return "unsupported";
    case Stop::zero_polynomial: return "zero polynomial";
    case Stop::indeterminate_degree: return "indeterminate degree";
    case Stop::error: return "error";
    }
    return "unknown";
}

// Whether every root is simple: proved, disproved (a root of multiplicity above one is proved), or not known.
enum class Squarefree { unknown = CPS_SQUAREFREE_UNKNOWN, proved = CPS_SQUAREFREE_PROVED, disproved = CPS_SQUAREFREE_DISPROVED };
// The piecewise method handed the solve to the secular method within the call, and why: it spent more than a global
// solve's predicted cost, its storage was refused, the coefficients' exponents exceed its range, or it reached its
// precision budget.
enum class Handoff {
    none = CPS_HANDOFF_NONE, regime_exit = CPS_HANDOFF_REGIME_EXIT, storage = CPS_HANDOFF_STORAGE,
    range = CPS_HANDOFF_RANGE, precision = CPS_HANDOFF_PRECISION
};

class Result;
namespace detail {
::certpolysolve::Result make_result(cps_error e, cps_result *r);
::certpolysolve::Result error_result(std::string message);
} // namespace detail

// What a solve returns: never throws, always tells why it stopped.
class [[nodiscard]] Result {
    std::shared_ptr<const cps_result> handle_;
    std::vector<Root> roots_;
    Stop stop_ = Stop::error;
    std::string message_;
    bool complete_ = false, goal_met_ = false;
    int64_t degree_ = -1;
    friend Result detail::make_result(cps_error, cps_result *);
    friend Result detail::error_result(std::string);

  public:
    // The goal is met: the roots answer the question asked.
    bool ok() const noexcept { return stop_ == Stop::none && goal_met_; }
    Stop stop() const noexcept { return stop_; }
    // The error, or the solver's diagnostics.
    const std::string &message() const noexcept { return message_; }
    // The regions hold every root (also after a stop, when the cover is complete).
    bool complete() const noexcept { return complete_; }
    bool goal_met() const noexcept { return goal_met_; }
    const std::vector<Root> &roots() const noexcept { return roots_; }
    // The roots as numbers: each region's value (Root::value, the numbers cps_result_get_approximation gives),
    // repeated as many times as the region holds roots (a double root gives two copies), in the regions' order.
    // Doubles, not enclosures: the proved boxes are in roots(). A region without a value is left out, so the list has
    // degree() numbers exactly when every region has one; after a stop or an error it is empty.
    std::vector<std::complex<double>> values() const {
        std::vector<std::complex<double>> v;
        for (const Root &r : roots_)
            if (r.value) v.insert(v.end(), static_cast<size_t>(r.count), *r.value);
        return v;
    }
    int64_t degree() const noexcept { return degree_; }
    uint64_t work_bits() const noexcept { return handle_ ? cps_result_work_bits(handle_.get()) : 0; }
    // Exact preparation was refused or limited by the resource limits, so multiplicities may be unknown.
    bool preparation_limited() const noexcept { return handle_ && cps_result_preparation_limited(handle_.get()); }
    // After a stop at the precision cap on exact real input: a cap at which the roots can be isolated at the
    // requested precision (Mahler's bound; a bound, not a promise that the solve then meets its goal). 0 otherwise.
    uint64_t sufficient_work_bits() const noexcept { return handle_ ? cps_result_sufficient_work_bits(handle_.get()) : 0; }
    Squarefree squarefree() const noexcept {
        return handle_ ? static_cast<Squarefree>(cps_result_squarefree(handle_.get())) : Squarefree::unknown;
    }
    Handoff handoff() const noexcept {
        return handle_ ? static_cast<Handoff>(cps_result_handoff(handle_.get())) : Handoff::none;
    }
    // For experts: the C result, what the call spent (a promise: work units, sweeps, precision, regions, times), the
    // engines' counters (diagnostics, which may change in any release) and the search's components (zeros and
    // empty for an error).
    const cps_result *get() const noexcept { return handle_.get(); }
    cps_statistics statistics() const noexcept {
        cps_statistics s{};
        s.size = sizeof s;
        if (handle_) (void)cps_result_get_statistics(handle_.get(), &s);
        return s;
    }
    cps_engine_counters engine_counters() const noexcept {
        cps_engine_counters c{};
        c.size = sizeof c;
        if (handle_) (void)cps_result_get_engine_counters(handle_.get(), &c);
        return c;
    }
    std::vector<cps_component_progress> components() const {
        std::vector<cps_component_progress> out;
        if (!handle_) return out;
        out.resize(cps_result_component_count(handle_.get()));
        for (size_t i = 0; i < out.size(); ++i) (void)cps_result_get_component(handle_.get(), i, &out[i]);
        return out;
    }
};

// The goal: every root alone in its region to the accuracy asked (accuracy, the default for coefficient and secular
// input), clusters allowed when each region's centre approximates its roots (approximation, the default for
// programs), distinct roots each alone in a region of any width (isolation, exact coefficient input), or clusters of radius r.
class Goal {
    cps_goal goal_ = CPS_ROOT_ACCURACY;
    std::optional<Number> radius_;
    explicit Goal(cps_goal g, std::optional<Number> r = {}) : goal_(g), radius_(std::move(r)) {}

  public:
    static Goal accuracy() { return Goal(CPS_ROOT_ACCURACY); }
    static Goal approximation() { return Goal(CPS_ROOT_APPROXIMATION); }
    static Goal isolation() { return Goal(CPS_ROOT_ISOLATION); }
    static Goal clusters(Number radius) { return Goal(CPS_CLUSTER_RADIUS, std::move(radius)); }
    cps_goal value() const noexcept { return goal_; }
    const std::optional<Number> &radius() const noexcept { return radius_; }
};

enum class Method { automatic = CPS_METHOD_AUTO, classical = CPS_METHOD_CLASSICAL, secular = CPS_METHOD_SECULAR, piecewise = CPS_METHOD_PIECEWISE };

// Advanced: the exact-arithmetic allowances (docs/resource-policy.md), each with no fixed ceiling by default; the
// constructor fills in the defaults. flags |= CPS_RESOURCE_ALLOW_EXPONENT_EXPANSION lets a compact exponent
// ("1@1000000") expand within them.
struct ResourceLimits : cps_resource_limits {
    ResourceLimits() noexcept : cps_resource_limits{} { cps_resource_limits_default(this); }
};

// A solve's settings, each with its default; set them by name (`Options o; o.bits = 64;`) or with designated
// initializers in this order ({.bits = 64, .time_limit = 60s}).
struct Options {
    uint64_t bits = CPS_DEFAULT_BITS;            // relative accuracy of every root, in bits
    uint64_t max_bits = CPS_DEFAULT_MAX_WORK_BITS; // the most working precision the solver may use
    std::chrono::duration<double> time_limit{}; // zero: none
    std::optional<Goal> goal{};                 // empty: the input kind's default
    Method method = Method::automatic;          // advanced
    std::stop_token stop{};                     // advanced: cancellation from another thread
    std::optional<uint64_t> work_budget{};      // advanced: work units (empty or 0: the automatic budget)
    std::optional<uint64_t> cell_budget{};      // advanced: live cells (empty or 0: CPS_DEFAULT_LIVE_CELLS)
    std::optional<uint64_t> iteration_budget{}; // advanced: sweeps (empty: none)
    std::optional<Number> absolute_tolerance{}; // advanced: also accept a region within this distance of its roots
    std::optional<ResourceLimits> resources{};  // advanced: exact-arithmetic allowances (empty: the defaults)
    uint32_t workers = 1;                       // advanced: the threads a solve may use (1: none started; the same
                                                // answer at every count; cps_options_set_workers)
};

// A closed region of the real line for real_roots: between(a, b), at_least(a), at_most(b).
struct Region {
    std::optional<Number> lower, upper;
};
inline Region between(Number a, Number b) { return {std::move(a), std::move(b)}; }
inline Region at_least(Number a) { return {std::move(a), std::nullopt}; }
inline Region at_most(Number b) { return {std::nullopt, std::move(b)}; }

// A polynomial with exact (or ball) coefficients: Polynomial{a0, a1, ...} with the constant term first (also
// Polynomial::lowest_first), Polynomial::highest_first with the highest power first, a range of coefficients, a
// constant, its terms (Polynomial::sparse), or built with arithmetic from Polynomial::x(), from_roots and monomial.
// It keeps its nonzero terms only, so x^1000000000 + 1 is two terms. Reading: degree(), p[i], coefficients(), terms();
// printing with std::format and std::println. A bad value is remembered (ok(), status()).
class Polynomial {
    struct Deleter {
        void operator()(cps_polynomial *p) const noexcept { cps_polynomial_destroy(p); }
    };
    using Handle = std::unique_ptr<cps_polynomial, Deleter>;
    Handle handle_;
    std::string status_;
    explicit Polynomial(Handle owned) noexcept : handle_(std::move(owned)) {}
    static cps_polynomial *fresh() {
        cps_polynomial *p = nullptr;
        detail::memory(cps_polynomial_create(&p));
        return p;
    }
    void fail(cps_error e) {
        detail::memory(e);
        if (e != CPS_OK && status_.empty()) status_ = cps_last_error();
    }
    // The result of an arithmetic call: its polynomial, or a zero polynomial that remembers the error. The message is
    // read before anything else calls the library (a call that succeeds clears it).
    static Polynomial from_call(cps_error e, cps_polynomial *out) {
        detail::memory(e);
        // Adopt a successful result before doing anything that could allocate.
        // Constructing an empty polynomial first could throw and leak `out`.
        if (e == CPS_OK) return Polynomial(Handle(out));
        std::string message = cps_last_error();
        Polynomial p;
        p.status_ = message.empty() ? "polynomial arithmetic failed" : std::move(message);
        return p;
    }
    template <class Call> static Polynomial binary(const Polynomial &a, const Polynomial &b, Call call) {
        if (!a.ok()) return a;
        if (!b.ok()) return b;
        cps_polynomial *out = nullptr;
        const cps_error e = call(a.get(), b.get(), &out);
        return from_call(e, out);
    }
    std::string part(uint64_t exponent, cps_coefficient_part which) const {
        return detail::text([&](char *b, size_t n, size_t *l) { return cps_polynomial_get_coefficient(get(), exponent, which, b, n, l); });
    }

  public:
    Polynomial() : handle_(fresh()) {}
    template <class T> requires std::constructible_from<Coefficient, T> Polynomial(T &&constant) : Polynomial() {
        set(0, Coefficient(std::forward<T>(constant)));
    }
    Polynomial(std::initializer_list<Coefficient> ascending) : Polynomial() {
        uint64_t i = 0;
        for (const auto &c : ascending) set(i++, c);
    }
    template <std::ranges::input_range R>
        requires std::constructible_from<Coefficient, std::ranges::range_reference_t<R>> && (!std::constructible_from<Coefficient, R>)
    explicit Polynomial(R &&ascending) : Polynomial() {
        uint64_t i = 0;
        for (auto &&c : ascending) set(i++, Coefficient(c));
    }
    Polynomial(const Polynomial &o) : status_(o.status_) {
        cps_polynomial *p = nullptr;
        detail::memory(cps_polynomial_copy(o.get(), &p));
        handle_.reset(p);
    }
    Polynomial &operator=(const Polynomial &o) {
        if (this != &o) *this = Polynomial(o);
        return *this;
    }
    Polynomial(Polynomial &&) noexcept = default;
    Polynomial &operator=(Polynomial &&) noexcept = default;
    ~Polynomial() = default;

    // A list either way round, named. lowest_first({-2, 0, 1}) is x^2 - 2 with the constant term first, as the
    // brace form reads a list; highest_first({1, 0, -2}) is the same polynomial with the highest power first, as
    // MATLAB's roots and numpy.roots read one. Zeros may stand anywhere; the degree is that of the highest nonzero
    // coefficient. Print the polynomial to see which one you built.
    static Polynomial lowest_first(std::initializer_list<Coefficient> coefficients) { return Polynomial(coefficients); }
    template <std::ranges::input_range R>
        requires std::constructible_from<Coefficient, std::ranges::range_reference_t<R>> && (!std::constructible_from<Coefficient, R>)
    static Polynomial lowest_first(R &&coefficients) {
        return Polynomial(std::forward<R>(coefficients));
    }
    static Polynomial highest_first(std::initializer_list<Coefficient> coefficients) {
        Polynomial p;
        uint64_t exponent = coefficients.size();
        for (const auto &c : coefficients) p.set(--exponent, c);
        return p;
    }
    template <std::ranges::forward_range R>
        requires std::constructible_from<Coefficient, std::ranges::range_reference_t<R>> && (!std::constructible_from<Coefficient, R>)
    static Polynomial highest_first(R &&coefficients) {
        Polynomial p;
        uint64_t exponent = static_cast<uint64_t>(std::ranges::distance(coefficients));
        for (auto &&c : coefficients) p.set(--exponent, Coefficient(c));
        return p;
    }

    // x, and c x^exponent.
    static Polynomial x() { return Polynomial{0, 1}; }
    static Polynomial monomial(uint64_t exponent, const Coefficient &c = 1) {
        Polynomial p;
        p.set(exponent, c);
        return p;
    }
    // The product of (x - r): from a range of numbers, or of pairs (root, multiplicity).
    template <std::ranges::input_range R> static Polynomial from_roots(R &&roots) {
        std::vector<Number> values;
        std::vector<uint64_t> multiplicities;
        for (auto &&r : roots) {
            if constexpr (requires { r.first; r.second; }) {
                values.emplace_back(r.first);
                multiplicities.push_back(static_cast<uint64_t>(r.second));
            } else {
                values.emplace_back(r);
                multiplicities.push_back(1);
            }
        }
        std::vector<const char *> re, im;
        for (const auto &v : values) re.push_back(v.real().c_str()), im.push_back(v.imag().c_str());
        cps_polynomial *out = nullptr;
        const cps_error e = cps_polynomial_from_roots(re.data(), im.data(), multiplicities.data(), values.size(), &out);
        return from_call(e, out);
    }
    static Polynomial from_roots(std::initializer_list<Number> roots) { return from_roots(std::vector<Number>(roots)); }

    Polynomial &set(uint64_t exponent, const Coefficient &c) {
        fail(cps_polynomial_set_ball_coefficient(handle_.get(), exponent, c.value.real().c_str(), c.value.imag().c_str(),
                                                c.real_radius.real().c_str(), c.imag_radius.real().c_str()));
        return *this;
    }
    // The coefficient of x^exponent, exact (zero beyond the degree).
    Coefficient coefficient(uint64_t exponent) const {
        Coefficient c(Number::complex(part(exponent, CPS_PART_REAL), part(exponent, CPS_PART_IMAG)));
        c.real_radius = part(exponent, CPS_PART_REAL_RADIUS);
        c.imag_radius = part(exponent, CPS_PART_IMAG_RADIUS);
        return c;
    }
    Coefficient operator[](uint64_t exponent) const { return coefficient(exponent); }
    // -1 for the zero polynomial.
    int64_t degree() const noexcept { return cps_polynomial_degree(get()); }
    bool is_zero() const noexcept { return degree() < 0; }
    // degree() + 1 coefficients, ascending (for a polynomial of few terms and a high degree, read terms()).
    std::vector<Coefficient> coefficients() const {
        std::vector<Coefficient> c;
        for (int64_t k = 0; k <= degree(); ++k) c.push_back(coefficient(static_cast<uint64_t>(k)));
        return c;
    }
    struct Term {
        uint64_t exponent;
        Coefficient coefficient;
    };
    // A polynomial from its terms, in any order: Polynomial::sparse({{1000000, 1}, {7, -3}, {0, "1/2"}}).
    static Polynomial sparse(std::initializer_list<Term> terms) {
        Polynomial p;
        for (const auto &t : terms) p.set(t.exponent, t.coefficient);
        return p;
    }
    template <std::ranges::input_range R> requires std::convertible_to<std::ranges::range_reference_t<R>, Term>
    static Polynomial sparse(R &&terms) {
        Polynomial p;
        for (auto &&t : terms) {
            const Term term = t;
            p.set(term.exponent, term.coefficient);
        }
        return p;
    }
    // The nonzero terms, ascending.
    std::vector<Term> terms() const {
        size_t count = 0;
        detail::memory(cps_polynomial_terms(get(), nullptr, 0, &count));
        std::vector<uint64_t> exponents(count);
        detail::memory(cps_polynomial_terms(get(), exponents.data(), exponents.size(), &count));
        std::vector<Term> t;
        t.reserve(exponents.size());
        for (const uint64_t e : exponents) t.push_back({e, coefficient(e)});
        return t;
    }
    bool ok() const noexcept { return status_.empty(); }
    const std::string &status() const noexcept { return status_; }
    const cps_polynomial *get() const noexcept { return handle_.get(); }

    Polynomial derivative() const {
        if (!ok()) return *this;
        cps_polynomial *out = nullptr;
        const cps_error e = cps_polynomial_derivative(get(), &out);
        return from_call(e, out);
    }
    friend Polynomial operator+(const Polynomial &a, const Polynomial &b) { return binary(a, b, cps_polynomial_add); }
    friend Polynomial operator-(const Polynomial &a, const Polynomial &b) { return binary(a, b, cps_polynomial_subtract); }
    friend Polynomial operator*(const Polynomial &a, const Polynomial &b) { return binary(a, b, cps_polynomial_multiply); }
    friend Polynomial operator-(const Polynomial &a) { return Polynomial() - a; }
    Polynomial &operator+=(const Polynomial &b) { return *this = *this + b; }
    Polynomial &operator-=(const Polynomial &b) { return *this = *this - b; }
    Polynomial &operator*=(const Polynomial &b) { return *this = *this * b; }
    friend Polynomial pow(const Polynomial &a, uint64_t exponent) {
        if (!a.ok()) return a;
        cps_polynomial *out = nullptr;
        const cps_error e = cps_polynomial_power(a.get(), exponent, &out);
        return from_call(e, out);
    }
    // Equal coefficients (exact values, and the same radii for balls): the same terms, since the library writes one
    // text for one value.
    friend bool operator==(const Polynomial &a, const Polynomial &b) {
        if (a.degree() != b.degree()) return false;
        const auto ta = a.terms(), tb = b.terms();
        if (ta.size() != tb.size()) return false;
        for (size_t k = 0; k < ta.size(); ++k) {
            const Coefficient &x = ta[k].coefficient, &y = tb[k].coefficient;
            if (ta[k].exponent != tb[k].exponent || x.value.real() != y.value.real() || x.value.imag() != y.value.imag() ||
                x.real_radius.real() != y.real_radius.real() || x.imag_radius.real() != y.imag_radius.real())
                return false;
        }
        return true;
    }
};

// A native secular equation, a sum of residue / (x - pole) terms, minus one: SecularEquation{{.residue = 1, .pole =
// 2}, ...}. Immutable once made; copies share it.
class SecularEquation {
    std::shared_ptr<const cps_secular> handle_;
    std::string status_;

  public:
    struct Term {
        Number residue;
        Number pole;
    };
    SecularEquation(std::initializer_list<Term> terms) : SecularEquation(std::vector<Term>(terms)) {}
    template <std::ranges::input_range R> requires std::convertible_to<std::ranges::range_reference_t<R>, Term>
    explicit SecularEquation(R &&terms) {
        std::vector<Term> held(std::ranges::begin(terms), std::ranges::end(terms));
        std::vector<cps_secular_term_view> views;
        for (const auto &t : held)
            views.push_back({t.residue.real().c_str(), t.residue.imag().c_str(), t.pole.real().c_str(), t.pole.imag().c_str()});
        cps_secular *s = nullptr;
        const cps_error e = cps_secular_create(views.data(), views.size(), nullptr, &s);
        detail::memory(e);
        if (e == CPS_OK) handle_.reset(s, cps_secular_destroy);
        else status_ = cps_last_error();
    }
    bool ok() const noexcept { return status_.empty(); }
    const std::string &status() const noexcept { return status_; }
    const cps_secular *get() const noexcept { return handle_.get(); }
};

// Program input: a polynomial given by the program that computes it.
// Program::from(f) calls f once with a ProgramValue, whose operators record the program; constants are plain
// numbers. A bad value is remembered (ok(), status()).
class Program;
class ProgramValue {
    struct Recording {
        cps_program *program = nullptr;
        std::string status;
    };
    std::shared_ptr<Recording> recording_;
    uint64_t node_ = 0;
    friend class Program;
    ProgramValue(std::shared_ptr<Recording> r, uint64_t node) : recording_(std::move(r)), node_(node) {}
    // One recorded operation: call(&node) adds it to the program and writes its node. The call is made here, before the
    // node is read. (As record(call(..., &n), n) the two were arguments of one call, whose order C++ leaves open: GCC
    // read n first, and every program recorded through Program::from became the variable alone; 19 October 2026.)
    template <class Call> ProgramValue record(Call &&call) const {
        uint64_t node = 0;
        const cps_error e = call(&node);
        detail::memory(e);
        if (e != CPS_OK && recording_->status.empty()) recording_->status = cps_last_error();
        return ProgramValue(recording_, e == CPS_OK ? node : node_);
    }
    ProgramValue constant(const Number &c) const {
        return record([&](uint64_t *node) {
            return cps_program_constant(recording_->program, c.real().c_str(), c.imag().c_str(), node);
        });
    }
    // The recording's first mistake, kept; the value stays this one.
    ProgramValue refuse(const char *why) const {
        if (recording_->status.empty()) recording_->status = why;
        return *this;
    }
    template <class Call> ProgramValue apply(Call call, const ProgramValue &b) const {
        // A node number means something only in its own program: an operand recorded by another one is refused.
        if (b.recording_ != recording_) return refuse("an operand belongs to another program's recording");
        return record([&](uint64_t *node) { return call(recording_->program, node_, b.node_, node); });
    }
    template <class T> ProgramValue lift(const T &c) const {
        if constexpr (std::is_same_v<T, ProgramValue>) return c;
        else return constant(Number(c));
    }

  public:
    friend ProgramValue operator+(const ProgramValue &a, const ProgramValue &b) { return a.apply(cps_program_add, b); }
    friend ProgramValue operator-(const ProgramValue &a, const ProgramValue &b) { return a.apply(cps_program_subtract, b); }
    friend ProgramValue operator*(const ProgramValue &a, const ProgramValue &b) { return a.apply(cps_program_multiply, b); }
    friend ProgramValue operator-(const ProgramValue &a) {
        return a.record([&](uint64_t *node) { return cps_program_negate(a.recording_->program, a.node_, node); });
    }
    template <class T> requires std::constructible_from<Number, T> friend ProgramValue operator+(const ProgramValue &a, const T &b) { return a + a.lift(b); }
    template <class T> requires std::constructible_from<Number, T> friend ProgramValue operator+(const T &a, const ProgramValue &b) { return b.lift(a) + b; }
    template <class T> requires std::constructible_from<Number, T> friend ProgramValue operator-(const ProgramValue &a, const T &b) { return a - a.lift(b); }
    template <class T> requires std::constructible_from<Number, T> friend ProgramValue operator-(const T &a, const ProgramValue &b) { return b.lift(a) - b; }
    template <class T> requires std::constructible_from<Number, T> friend ProgramValue operator*(const ProgramValue &a, const T &b) { return a * a.lift(b); }
    template <class T> requires std::constructible_from<Number, T> friend ProgramValue operator*(const T &a, const ProgramValue &b) { return b.lift(a) * b; }
    template <class T> ProgramValue &operator+=(const T &b) { return *this = *this + b; }
    template <class T> ProgramValue &operator-=(const T &b) { return *this = *this - b; }
    template <class T> ProgramValue &operator*=(const T &b) { return *this = *this * b; }
    // a^0 is the constant 1 and a^1 is a (the program format's power takes an exponent of 2 or more).
    friend ProgramValue pow(const ProgramValue &a, uint64_t exponent) {
        if (exponent == 0) return a.constant(Number(1));
        if (exponent == 1) return a;
        return a.record([&](uint64_t *node) { return cps_program_power(a.recording_->program, a.node_, exponent, node); });
    }
};

class Program {
    std::shared_ptr<ProgramValue::Recording> recording_;
    std::shared_ptr<cps_program> handle_;
    void own(cps_program *p) {
        handle_.reset(p, cps_program_destroy);
        recording_ = std::make_shared<ProgramValue::Recording>();
        recording_->program = p;
    }

  public:
    Program() {
        cps_program *p = nullptr;
        detail::memory(cps_program_create(&p));
        own(p);
    }
    // A program from its text form ("PSPROGRAM 1", then one operation a line; cps_program_parse).
    explicit Program(std::string_view text) {
        const std::string held(text);
        cps_program *p = nullptr;
        const cps_error e = cps_program_parse(held.c_str(), &p);
        detail::memory(e);
        if (e == CPS_OK) own(p);
        else {
            std::string message = cps_last_error(); // before the next library call clears it
            *this = Program();
            recording_->status = message.empty() ? "the program text was refused" : std::move(message);
        }
    }
    template <class Callable> static Program from(Callable &&f) {
        Program program;
        uint64_t x = 0;
        const cps_error e = cps_program_variable(program.handle_.get(), &x);
        detail::memory(e);
        if (e != CPS_OK) {
            program.recording_->status = cps_last_error();
            return program;
        }
        const ProgramValue variable(program.recording_, x);
        const ProgramValue result = variable.lift(std::forward<Callable>(f)(variable));
        if (result.recording_ != program.recording_) {
            if (program.recording_->status.empty()) program.recording_->status = "the program's result was not recorded by it";
            return program;
        }
        const cps_error out = cps_program_set_output(program.handle_.get(), result.node_);
        detail::memory(out);
        if (out != CPS_OK && program.recording_->status.empty()) program.recording_->status = cps_last_error();
        return program;
    }
    bool ok() const noexcept { return recording_->status.empty(); }
    const std::string &status() const noexcept { return recording_->status; }
    const cps_program *get() const noexcept { return handle_.get(); }
};

namespace detail {
inline ::certpolysolve::Result error_result(std::string message) {
    ::certpolysolve::Result r;
    r.stop_ = ::certpolysolve::Stop::error;
    r.message_ = std::move(message);
    return r;
}
inline ::certpolysolve::Stop stop_of(cps_termination t, cps_stop_reason s) {
    switch (t) {
    case CPS_SUCCESS: return ::certpolysolve::Stop::none;
    case CPS_PRECISION_LIMIT: return s == CPS_STOP_INPUT ? ::certpolysolve::Stop::input_uncertainty : ::certpolysolve::Stop::precision_cap;
    case CPS_WORK_LIMIT:
        return s == CPS_STOP_CELL_BUDGET ? ::certpolysolve::Stop::cell_budget : s == CPS_STOP_ITERATION_BUDGET ? ::certpolysolve::Stop::iteration_budget : ::certpolysolve::Stop::work_budget;
    case CPS_CANCELLED: return s == CPS_STOP_TIME_LIMIT ? ::certpolysolve::Stop::time_limit : ::certpolysolve::Stop::cancelled;
    case CPS_UNSUPPORTED_QUERY: return ::certpolysolve::Stop::unsupported;
    case CPS_ZERO_POLYNOMIAL: return ::certpolysolve::Stop::zero_polynomial;
    case CPS_INDETERMINATE_DEGREE: return ::certpolysolve::Stop::indeterminate_degree;
    }
    return ::certpolysolve::Stop::error;
}
inline ::certpolysolve::Result make_result(cps_error e, cps_result *raw) {
    memory(e);
    if (e != CPS_OK) return error_result(cps_last_error());
    ::certpolysolve::Result r;
    r.handle_.reset(raw, cps_result_destroy);
    r.stop_ = stop_of(cps_result_termination(raw), cps_result_stop_reason(raw));
    r.message_ = cps_result_diagnostics(raw);
    r.complete_ = cps_result_complete(raw) != 0;
    r.goal_met_ = cps_result_goal_met(raw) != 0;
    r.degree_ = cps_result_degree(raw);
    for (size_t i = 0; i < cps_result_size(raw); ++i) {
        cps_root_view v{};
        memory(cps_result_get_root(raw, i, &v));
        ::certpolysolve::Root root;
        root.box = {{::certpolysolve::Dyadic(v.real_lower), ::certpolysolve::Dyadic(v.real_upper)}, {::certpolysolve::Dyadic(v.imag_lower), ::certpolysolve::Dyadic(v.imag_upper)}};
        root.count = v.count;
        // The value as C gives it, or none where C refuses (the refusal's message is the C call's own).
        double re = 0, im = 0;
        if (cps_result_get_approximation(raw, i, &re, &im) == CPS_OK) root.value = std::complex<double>(re, im);
        if (v.multiplicity) root.multiplicity = v.multiplicity;
        root.realness = v.realness == CPS_REAL ? ::certpolysolve::Realness::real : v.realness == CPS_NONREAL ? ::certpolysolve::Realness::not_real : ::certpolysolve::Realness::unknown;
        root.certificate = v.certificate ? v.certificate : "";
        r.roots_.push_back(std::move(root));
    }
    return r;
}
// The C options of a solve, its token linked to the caller's std::stop_token.
struct Settings {
    struct OptionsDeleter {
        void operator()(cps_options *o) const noexcept { cps_options_destroy(o); }
    };
    struct TokenDeleter {
        void operator()(cps_cancel_token *t) const noexcept { cps_cancel_token_destroy(t); }
    };
    struct Request {
        cps_cancel_token *token;
        void operator()() const noexcept { cps_cancel_token_request(token); }
    };
    std::unique_ptr<cps_options, OptionsDeleter> options;
    std::unique_ptr<cps_cancel_token, TokenDeleter> token;
    std::optional<std::stop_callback<Request>> link;
    std::string error;
    void check(cps_error e) {
        memory(e);
        if (e != CPS_OK && error.empty()) error = cps_last_error();
    }
};
inline std::unique_ptr<Settings> settings(const ::certpolysolve::Options &o, bool real_only, const ::certpolysolve::Region *region) {
    auto s = std::make_unique<Settings>();
    cps_options *raw = nullptr;
    memory(cps_options_create(&raw));
    s->options.reset(raw);
    s->check(cps_options_set_precision(raw, o.bits));
    s->check(cps_options_set_max_work_bits(raw, o.max_bits));
    if (o.time_limit.count() != 0) s->check(cps_options_set_time_limit(raw, o.time_limit.count()));
    if (o.goal) {
        if (o.goal->radius()) s->check(cps_options_set_cluster_radius(raw, o.goal->radius()->real().c_str()));
        else s->check(cps_options_set_goal(raw, o.goal->value()));
    }
    s->check(cps_options_set_method(raw, static_cast<int>(o.method)));
    if (o.work_budget) s->check(cps_options_set_work_budget(raw, *o.work_budget));
    if (o.cell_budget) s->check(cps_options_set_live_cells(raw, *o.cell_budget));
    if (o.iteration_budget) s->check(cps_options_set_max_iterations(raw, *o.iteration_budget));
    if (o.absolute_tolerance) s->check(cps_options_set_absolute_tolerance(raw, o.absolute_tolerance->real().c_str()));
    if (o.resources) s->check(cps_options_set_resource_limits(raw, &*o.resources));
    if (o.workers != 1) s->check(cps_options_set_workers(raw, o.workers));
    if (real_only) s->check(cps_options_set_real_only(raw, 1));
    if (region)
        s->check(cps_options_set_region(raw, region->lower ? region->lower->real().c_str() : "",
                                          region->upper ? region->upper->real().c_str() : "", "", ""));
    if (o.stop.stop_possible()) {
        cps_cancel_token *t = nullptr;
        memory(cps_cancel_token_create(&t));
        s->token.reset(t);
        s->check(cps_options_set_cancel_token(raw, t));
        s->link.emplace(o.stop, Settings::Request{t});
    }
    return s;
}
template <class Input, class Solve> ::certpolysolve::Result run(const Input &input, const ::certpolysolve::Options &o, bool real_only, const ::certpolysolve::Region *region, Solve solve) {
    if (!input.ok()) return error_result(input.status());
    const auto s = settings(o, real_only, region);
    if (!s->error.empty()) return error_result(s->error);
    cps_result *raw = nullptr;
    const cps_error e = solve(input.get(), s->options.get(), &raw);
    return make_result(e, raw);
}
} // namespace detail

// Every root (complex roots), or the real roots, in a region if one is given.
[[nodiscard]] inline Result solve(const Polynomial &p, const Options &o = {}) { return detail::run(p, o, false, nullptr, cps_solve); }
[[nodiscard]] inline Result solve(const SecularEquation &s, const Options &o = {}) { return detail::run(s, o, false, nullptr, cps_solve_secular); }
[[nodiscard]] inline Result solve(const Program &p, const Options &o = {}) { return detail::run(p, o, false, nullptr, cps_solve_program); }
[[nodiscard]] inline Result real_roots(const Polynomial &p, const Options &o = {}) { return detail::run(p, o, true, nullptr, cps_solve); }
[[nodiscard]] inline Result real_roots(const Polynomial &p, const Region &r, const Options &o = {}) {
    return detail::run(p, o, true, &r, cps_solve);
}
[[nodiscard]] inline Result real_roots(const SecularEquation &s, const Options &o = {}) { return detail::run(s, o, true, nullptr, cps_solve_secular); }

namespace detail {
// A coefficient part for printing: a dyadic with a small exponent as its fraction, any other text as given (a large
// exponent, or text the library has not checked yet: nothing here throws).
inline std::string printable(const std::string &t) {
    constexpr long long fraction_exponents = 64; // |e| up to this prints as p/q; beyond, m@e
    const auto at = t.find('@');
    if (at == std::string::npos) return t;
    long long e = 0;
    const char *last = t.data() + t.size();
    const auto read = std::from_chars(t.data() + at + 1, last, e);
    if (read.ec != std::errc{} || read.ptr != last || e < -fraction_exponents || e > fraction_exponents) return t;
    return ::certpolysolve::Dyadic::parse(t).to_string();
}
inline bool negative(const std::string &t) { return !t.empty() && t[0] == '-'; }
// A number or a coefficient on its own: "3", "-1/3", "1/3 - 2i", "3/4 +/- 1/1000".
inline std::string coefficient_text(const ::certpolysolve::Coefficient &c) {
    const auto re = printable(c.value.real()), im = printable(c.value.imag());
    std::string s;
    if (im == "0") s = re;
    else if (re == "0") s = im + "i";
    else s = re + (negative(im) ? " - " + im.substr(1) : " + " + im) + "i";
    if (!c.exact()) {
        const auto rr = printable(c.real_radius.real()), ri = printable(c.imag_radius.real());
        s += " +/- " + (ri == "0" ? rr : rr + " + " + ri + "i");
    }
    return s;
}
// One term of a printed polynomial: its sign, and its coefficient as it stands before x^k ("3", "1/3" alone,
// "(1/3)" before x, "(1 - 2i)", "(3/4 +/- 1/1000)"; empty for a unit before x).
struct TermText {
    bool minus = false;
    std::string coefficient;
};
inline TermText term_text(const ::certpolysolve::Coefficient &c, uint64_t exponent) {
    TermText t;
    if (!c.exact() || printable(c.value.imag()) != "0") {
        t.coefficient = "(" + coefficient_text(c) + ")";
        return t;
    }
    t.coefficient = printable(c.value.real());
    if (negative(t.coefficient)) t.minus = true, t.coefficient.erase(0, 1);
    if (exponent == 0) return t;
    if (t.coefficient == "1") t.coefficient.clear();
    else if (t.coefficient.find_first_of("/@") != std::string::npos) t.coefficient = "(" + t.coefficient + ")";
    return t;
}
} // namespace detail
} // namespace certpolysolve

// Printing: a Dyadic in decimal ({} 17 significant digits, {:.N} N digits); an Interval "[lo, hi]"; a Root its centre
// ({} to double precision, {:.N} with N digits from the exact box); a Stop its name; a Number or a Coefficient
// exactly ("-1/3", "1 - 2i", "3/4 +/- 1/1000"); a Polynomial descending, exact ("x^3 - (1/2)x - 7/3").
namespace certpolysolve::detail {
struct digits_spec {
    std::optional<unsigned> digits;
    constexpr auto parse(std::format_parse_context &ctx) {
        auto it = ctx.begin();
        if (it != ctx.end() && *it == '.') {
            unsigned d = 0;
            for (++it; it != ctx.end() && *it >= '0' && *it <= '9'; ++it) d = 10 * d + static_cast<unsigned>(*it - '0');
            digits = d;
        }
        if (it != ctx.end() && *it != '}') throw std::format_error("certpolysolve: the format is {} or {:.N}");
        return it;
    }
};
} // namespace certpolysolve::detail
template <> struct std::formatter<certpolysolve::Dyadic> : certpolysolve::detail::digits_spec {
    auto format(const certpolysolve::Dyadic &d, std::format_context &ctx) const {
        return std::format_to(ctx.out(), "{}", d.format(digits.value_or(17)));
    }
};
template <> struct std::formatter<certpolysolve::Interval> : certpolysolve::detail::digits_spec {
    auto format(const certpolysolve::Interval &i, std::format_context &ctx) const {
        const unsigned d = digits.value_or(17);
        return std::format_to(ctx.out(), "[{}, {}]", i.lo.format(d), i.hi.format(d));
    }
};
template <> struct std::formatter<certpolysolve::Root> : certpolysolve::detail::digits_spec {
    auto format(const certpolysolve::Root &r, std::format_context &ctx) const {
        const bool real = r.box.imag.lo == 0 && r.box.imag.hi == 0;
        if (!digits) {
            const std::complex<double> c = r.centre();
            if (real) return std::format_to(ctx.out(), "{}", c.real());
            return std::format_to(ctx.out(), "{} {} {}i", c.real(), c.imag() < 0 ? '-' : '+', std::abs(c.imag()));
        }
        const auto re = r.box.real.mid().format(*digits);
        if (real) return std::format_to(ctx.out(), "{}", re);
        const auto im = r.box.imag.mid().format(*digits);
        if (!im.empty() && im[0] == '-') return std::format_to(ctx.out(), "{} - {}i", re, im.substr(1));
        return std::format_to(ctx.out(), "{} + {}i", re, im);
    }
};
template <> struct std::formatter<certpolysolve::Stop> : std::formatter<std::string_view> {
    auto format(certpolysolve::Stop s, std::format_context &ctx) const {
        return std::formatter<std::string_view>::format(certpolysolve::name(s), ctx);
    }
};
template <> struct std::formatter<certpolysolve::Number> {
    constexpr auto parse(std::format_parse_context &ctx) { return ctx.begin(); }
    auto format(const certpolysolve::Number &n, std::format_context &ctx) const {
        return std::format_to(ctx.out(), "{}", certpolysolve::detail::coefficient_text(certpolysolve::Coefficient(n)));
    }
};
template <> struct std::formatter<certpolysolve::Coefficient> {
    constexpr auto parse(std::format_parse_context &ctx) { return ctx.begin(); }
    auto format(const certpolysolve::Coefficient &c, std::format_context &ctx) const {
        return std::format_to(ctx.out(), "{}", certpolysolve::detail::coefficient_text(c));
    }
};
template <> struct std::formatter<certpolysolve::Polynomial> {
    constexpr auto parse(std::format_parse_context &ctx) { return ctx.begin(); }
    auto format(const certpolysolve::Polynomial &p, std::format_context &ctx) const {
        if (!p.ok()) return std::format_to(ctx.out(), "<bad polynomial: {}>", p.status());
        const auto terms = p.terms();
        if (terms.empty()) return std::format_to(ctx.out(), "0");
        std::string s;
        for (auto t = terms.rbegin(); t != terms.rend(); ++t) {
            const auto term = certpolysolve::detail::term_text(t->coefficient, t->exponent);
            if (s.empty()) s += term.minus ? "-" : "";
            else s += term.minus ? " - " : " + ";
            s += term.coefficient;
            if (t->exponent != 0) s += t->exponent == 1 ? "x" : "x^" + std::to_string(t->exponent);
        }
        return std::format_to(ctx.out(), "{}", s);
    }
};
