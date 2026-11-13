// Exact square-free decomposition over Q(i) of an exact polynomial with a nonreal coefficient (exact_complex_factors;
// no factors where an allowance or cancellation stops it), and Gaussian-rational ball enclosures (enclose_gaussian).
//
// exact_complex_factors admits exact, not all real coefficients within the allowances, splits off the exact zero roots
// as x^m (with_zero_factor), expands the rest, f, within exact_storage_limit, and refuses when exact_operations is
// below f's coefficient count. Where the compact-exponent rule refuses the expansion, the split prime's test still
// runs on the coefficients as written, and a pass is returned as squarefree_proved, with no factors. Two one-sided proofs follow, the split prime's (squarefree_modular), else the extension
// field's (modular_squarefree, ModularSquarefreeTest: F_p[t]/(t^2 + 1), a field as p = 2^31 - 1 is 3 mod 4): an image
// keeping the degree, with gcd(f, f') = 1 and no denominator divisible by p, proves f square-free (it is returned
// whole); a failure proves nothing. Otherwise Arithmetic decomposes the monic f: u = gcd(f, f'), w = f/u, then at
// multiplicity k = 1, 2, ... y = gcd(w, u) leaves the factor w/y, w = y and u = u/y. Divisions are checked exact
// (exact_divide) and the factors' powers multiplied back to f (multiply, equal); a mismatch throws std::logic_error.
//
// Arithmetic::gcd is multimodular: at the primes p = 1 (mod 4) above 2^61 (next_split_prime), the two maps Z[i] -> F_p,
// i -> +-r with r^2 = -1, give images (images) whose gcds (nmod_poly_gcd) of the least degree seen are lifted by
// Chinese remaindering and rational reconstruction; a candidate unchanged by one more prime and dividing both inputs is
// the gcd (proof at the function). Work is charged in words to exact_operations, sizes held to exact_storage_limit
// (charge, guard, storage); exhaustion or cancellation (BudgetExceeded, Cancelled) leaves no factors.
//
// State: none; Arithmetic (the remaining exact_operations) and ModularSquarefreeTest live for one call.
// Entry points: exact_complex_factors, called by SolveState::Impl::preprocess (solver.cpp) and prepare_squarefree_input
// (exact_preprocessing.cpp); enclose_gaussian, by solver.cpp, factor_state.hpp, factor_state.cpp and centering.cpp;
// extension_proof_words, read by tests/test_exact_complex.cpp.
// Not here: the split-prime proof at 998244353 and modular_proof_admitted (modular_squarefree.cpp), the real input's
// square-free factorization (exact_preprocessing.cpp), solving the factors (solver.cpp, FactorState::Gaussian).
// Design notes: docs/mathematics.md (Individual root certificates), docs/resource-policy.md (the allowances),
// docs/preparation-selection-findings-20260924.md (the order of the two proofs).
#include "numeric/exact/exact_complex.hpp"
#include "core/telemetry.hpp"
#include "numeric/exact/modular_squarefree.hpp"
#include <flint/fq_nmod.h>
#include <flint/fq_nmod_poly.h>
#include <flint/nmod.h>
#include <flint/nmod_poly.h>
#include <flint/ulong_extras.h>
#include <algorithm>
#include <optional>
#include <stdexcept>

namespace certpolysolve::detail {
namespace {
struct BudgetExceeded {};
struct Cancelled {};
uint64_t size(const Rational &q) {
    return resources::add(fmpz_bits(fmpq_numref(q.value())), fmpz_bits(fmpq_denref(q.value())));
}
void trim(GaussianPolynomial &p) {
    while (!p.empty() && p.back().zero())
        p.pop_back();
}
GaussianRational one() {
    GaussianRational z;
    fmpq_one(z.real.value());
    return z;
}
bool equal(const GaussianPolynomial &a, const GaussianPolynomial &b) {
    if (a.size() != b.size())
        return false;
    for (size_t i = 0; i < a.size(); ++i)
        if (!fmpq_equal(a[i].real.value(), b[i].real.value()) || !fmpq_equal(a[i].imag.value(), b[i].imag.value()))
            return false;
    return true;
}

// A word-sized prime p = 1 (mod 4) with a square root of -1: the two ring maps
// Z[i] -> F_p send i to that root and to its negative.
struct SplitPrime {
    nmod_t field;
    ulong unit; // unit^2 = -1
};
// Primes just above 2^61, in increasing order: each image carries 61 bits,
// and nmod arithmetic keeps two spare bits in the word.
constexpr ulong split_prime_floor = UINT64_C(1) << 61;
SplitPrime next_split_prime(ulong after) {
    ulong p = after;
    do
        p = n_nextprime(p, 1);
    while (p % 4 != 1);
    SplitPrime result;
    nmod_init(&result.field, p);
    result.unit = n_sqrtmod(p - 1, p);
    if (nmod_mul(result.unit, result.unit, result.field) != p - 1)
        throw std::logic_error("no square root of -1 modulo a prime congruent to 1 mod 4");
    return result;
}
// One nmod polynomial; the modulus is set per prime.
class ResiduePolynomial {
    nmod_poly_t value_;
  public:
    explicit ResiduePolynomial(ulong prime) { nmod_poly_init(value_, prime); }
    ~ResiduePolynomial() { nmod_poly_clear(value_); }
    ResiduePolynomial(const ResiduePolynomial &) = delete;
    ResiduePolynomial &operator=(const ResiduePolynomial &) = delete;
    nmod_poly_struct *value() & { return value_; }
    const nmod_poly_struct *value() const & { return value_; }
};

// Conservative guards for this exact field implementation. They bound
// operand/result sizes and the work charged, NOT process peak RAM or time
// inside one FLINT call. No global generic-ring tables or arithmetic/thread
// settings are initialized or changed.
//
// The unit of exact_operations is one word operation, as in the real path's
// estimate (exact_preprocessing.cpp, squarefree_effort: length^2 words). A
// rational operation is charged one unit plus its operands' size in words, so
// a budget bounds time as well as count: a count alone does not bound the
// time of operations on coefficients that grow.
// History: docs/code-history-20261021.md, exact_complex.cpp 1.
struct Arithmetic {
    const cps_options &options;
    const EnginePolicy *trace;
    uint64_t operations;
    uint64_t calls = 0;
    explicit Arithmetic(const cps_options &o, const EnginePolicy *t)
        : options(o), trace(t), operations(o.exact_operations) {}
    void checkpoint() const {
        if (options.cancelled())
            throw Cancelled{};
    }
    void charge(uint64_t words) {
        if (operations < words)
            throw BudgetExceeded{};
        // An unlimited budget stays unlimited.
        if (operations != UINT64_MAX)
            operations -= words;
    }
    void guard(const Rational &a, const Rational &b) {
        const auto bits = resources::add(size(a), size(b));
        charge(resources::add(1, bits / FLINT_BITS));
        // Cancellation is polled every 64 operations.
        if ((++calls & 63) == 0)
            checkpoint();
        // Uncancelled numerator/denominator products and a carry bit fit this
        // estimate for addition, subtraction, multiplication and division.
        if (resources::add(bits, 4) > options.exact_storage_limit())
            throw BudgetExceeded{};
    }
    Rational add(const Rational &a, const Rational &b, bool subtract = false) {
        guard(a, b);
        Rational c;
        if (subtract)
            fmpq_sub(c.value(), a.value(), b.value());
        else
            fmpq_add(c.value(), a.value(), b.value());
        return c;
    }
    Rational mul(const Rational &a, const Rational &b, bool divide = false) {
        guard(a, b);
        Rational c;
        if (divide) {
            if (fmpq_is_zero(b.value()))
                throw std::logic_error("zero divisor in Q(i) preprocessing");
            fmpq_div(c.value(), a.value(), b.value());
        } else
            fmpq_mul(c.value(), a.value(), b.value());
        return c;
    }
    GaussianRational sub(const GaussianRational &a, const GaussianRational &b) {
        return {add(a.real, b.real, true), add(a.imag, b.imag, true)};
    }
    GaussianRational product(const GaussianRational &a, const GaussianRational &b) {
        return {add(mul(a.real, b.real), mul(a.imag, b.imag), true),
                add(mul(a.real, b.imag), mul(a.imag, b.real))};
    }
    GaussianRational inverse(const GaussianRational &a) {
        auto norm = add(mul(a.real, a.real), mul(a.imag, a.imag));
        GaussianRational result{mul(a.real, norm, true), mul(a.imag, norm, true)};
        fmpq_neg(result.imag.value(), result.imag.value());
        return result;
    }
    void storage(const GaussianPolynomial &p) const {
        uint64_t remaining = options.exact_storage_limit();
        for (const auto &z : p) {
            const auto bits = resources::add(size(z.real), size(z.imag));
            if (bits > remaining)
                throw BudgetExceeded{};
            remaining -= bits;
        }
    }
    GaussianPolynomial monic(GaussianPolynomial p) {
        trim(p);
        if (!p.empty()) {
            const auto inv = inverse(p.back());
            for (auto &z : p)
                z = product(z, inv);
            storage(p);
        }
        return p;
    }
    GaussianPolynomial derivative(const GaussianPolynomial &p) {
        GaussianPolynomial d(p.size() > 1 ? p.size() - 1 : 0);
        for (size_t i = 1; i < p.size(); ++i) {
            GaussianRational factor;
            fmpq_set_ui(factor.real.value(), static_cast<ulong>(i), 1);
            d[i - 1] = product(p[i], factor);
        }
        storage(d);
        trim(d);
        return d;
    }
    std::pair<GaussianPolynomial, GaussianPolynomial> divide(GaussianPolynomial p,
                                                             const GaussianPolynomial &q) {
        checkpoint();
        if (q.empty())
            throw std::logic_error("zero polynomial divisor in Q(i) preprocessing");
        GaussianPolynomial quotient(p.size() >= q.size() ? p.size() - q.size() + 1 : 0);
        const auto inv = inverse(q.back());
        while (!p.empty() && p.size() >= q.size()) {
            const auto shift = p.size() - q.size();
            const auto multiplier = product(p.back(), inv);
            quotient[shift] = multiplier;
            for (size_t j = 0; j < q.size(); ++j)
                p[shift + j] = sub(p[shift + j], product(multiplier, q[j]));
            trim(p);
            storage(p);
            storage(quotient);
        }
        trim(quotient);
        return {std::move(quotient), std::move(p)};
    }
    GaussianPolynomial exact_divide(const GaussianPolynomial &p, const GaussianPolynomial &q) {
        auto [quotient, remainder] = divide(p, q);
        if (!remainder.empty())
            throw std::logic_error("nonexact division in Q(i) squarefree decomposition");
        return quotient;
    }
    // The images of p under the two maps at a split prime, or false when the
    // prime is unusable for p: a denominator is divisible by it, or a leading
    // coefficient's image vanishes (the image would lose degree).
    bool images(nmod_poly_t plus, nmod_poly_t minus, const GaussianPolynomial &p, const SplitPrime &prime) {
        const ulong modulus = prime.field.n;
        nmod_poly_zero(plus);
        nmod_poly_zero(minus);
        uint64_t words = 0;
        for (size_t k = 0; k < p.size(); ++k) {
            ulong part[2];
            const Rational *source[2] = {&p[k].real, &p[k].imag};
            for (int t = 0; t < 2; ++t) {
                const auto denominator = fmpz_fdiv_ui(fmpq_denref(source[t]->value()), modulus);
                if (!denominator)
                    return false;
                part[t] = nmod_mul(fmpz_fdiv_ui(fmpq_numref(source[t]->value()), modulus),
                                   n_invmod(denominator, modulus), prime.field);
                words = resources::add(words, 1 + size(*source[t]) / FLINT_BITS);
            }
            const ulong twist = nmod_mul(prime.unit, part[1], prime.field);
            nmod_poly_set_coeff_ui(plus, static_cast<slong>(k), nmod_add(part[0], twist, prime.field));
            nmod_poly_set_coeff_ui(minus, static_cast<slong>(k), nmod_sub(part[0], twist, prime.field));
        }
        charge(words);
        const auto degree = static_cast<slong>(p.size()) - 1;
        return nmod_poly_degree(plus) == degree && nmod_poly_degree(minus) == degree;
    }
    // The monic gcd in Q(i)[x], from its images at split primes.
    //
    // A Euclidean remainder sequence over Q(i) lets coefficients grow along
    // the sequence, about n^5 in all. Here each prime costs two word-sized
    // gcds, and the number of primes follows the size of the answer.
    // History: docs/code-history-20261021.md, exact_complex.cpp 2.
    //
    // For a prime p = 1 (mod 4) and either map i -> +-r, r^2 = -1, call the
    // prime good when no denominator of a or b is divisible by p and both
    // leading coefficients have nonzero images. The monic gcd G then has
    // p-integral coefficients (Gauss's lemma in the local ring), and its
    // image divides gcd(image a, image b): every good image's gcd has degree
    // at least deg G. Let d be the smallest degree seen, and g the monic
    // polynomial of degree d whose coefficients' real and imaginary parts are
    // recovered by the Chinese remainder theorem and rational reconstruction
    // from the images of that degree (an image pair u, v of a coefficient
    // alpha + i beta is alpha + r beta, alpha - r beta). If g divides a and b
    // exactly, then g divides G, and deg g = d >= deg G: g = G. Exact division
    // is the whole proof; wrong images, unlucky primes and a premature
    // reconstruction can only make it fail and cost more primes.
    GaussianPolynomial gcd(GaussianPolynomial a, GaussianPolynomial b) {
        trim(a);
        trim(b);
        if (b.empty())
            return monic(std::move(a));
        if (a.empty())
            return monic(std::move(b));
        if (a.size() == 1 || b.size() == 1)
            return {one()};
        TraceScope phase(trace, "gaussian_modular_gcd_images", 0, std::min(a.size(), b.size()));
        slong degree = static_cast<slong>(std::min(a.size(), b.size())) - 1;
        std::vector<Integer> residues; // real and imaginary parts of the d low coefficients
        Integer modulus;
        GaussianPolynomial previous;
        ulong after = split_prime_floor;
        for (;;) {
            checkpoint();
            const auto prime = next_split_prime(after);
            after = prime.field.n;
            ResiduePolynomial a_plus(after), a_minus(after), b_plus(after), b_minus(after), g_plus(after),
                g_minus(after);
            if (!images(a_plus.value(), a_minus.value(), a, prime) ||
                !images(b_plus.value(), b_minus.value(), b, prime))
                continue;
            // Two word-sized gcds of these lengths.
            charge(resources::multiply(2, resources::multiply(a.size(), b.size())));
            nmod_poly_gcd(g_plus.value(), a_plus.value(), b_plus.value());
            nmod_poly_gcd(g_minus.value(), a_minus.value(), b_minus.value());
            const auto found = nmod_poly_degree(g_plus.value());
            if (found != nmod_poly_degree(g_minus.value()) || found > degree)
                continue; // one of the two maps is unlucky, or both are
            if (found == 0)
                return {one()}; // a good image of degree 0 bounds deg G by 0
            if (found < degree || residues.empty()) {
                degree = found;
                residues.assign(2 * static_cast<size_t>(degree), Integer{});
                fmpz_one(modulus.value());
                previous.clear();
            }
            const ulong half = n_invmod(2, after);
            const ulong half_over_unit = nmod_mul(half, n_invmod(prime.unit, after), prime.field);
            for (slong k = 0; k < degree; ++k) {
                const ulong u = nmod_poly_get_coeff_ui(g_plus.value(), k),
                            v = nmod_poly_get_coeff_ui(g_minus.value(), k);
                const ulong real = nmod_mul(nmod_add(u, v, prime.field), half, prime.field);
                const ulong imag = nmod_mul(nmod_sub(u, v, prime.field), half_over_unit, prime.field);
                auto &re = residues[2 * static_cast<size_t>(k)], &im = residues[2 * static_cast<size_t>(k) + 1];
                fmpz_CRT_ui(re.value(), re.value(), modulus.value(), real, after, 0);
                fmpz_CRT_ui(im.value(), im.value(), modulus.value(), imag, after, 0);
            }
            fmpz_mul_ui(modulus.value(), modulus.value(), after);
            const auto modulus_words = 1 + fmpz_bits(modulus.value()) / FLINT_BITS;
            charge(resources::multiply(2 * static_cast<uint64_t>(degree), modulus_words));
            // Rational reconstruction of every part; stop at the first miss.
            GaussianPolynomial candidate(static_cast<size_t>(degree) + 1);
            candidate.back() = one();
            bool reconstructed = true;
            for (slong k = 0; k < degree && reconstructed; ++k) {
                auto &z = candidate[static_cast<size_t>(k)];
                reconstructed =
                    fmpq_reconstruct_fmpz(z.real.value(), residues[2 * static_cast<size_t>(k)].value(),
                                          modulus.value()) &&
                    fmpq_reconstruct_fmpz(z.imag.value(), residues[2 * static_cast<size_t>(k) + 1].value(),
                                          modulus.value());
            }
            charge(resources::multiply(2 * static_cast<uint64_t>(degree), modulus_words));
            if (!reconstructed) {
                previous.clear();
                continue;
            }
            // Divide only once one more prime has left the candidate as it
            // was: a reconstruction from too few primes is usually wrong, and
            // its large coefficients make the trial division expensive.
            if (!equal(candidate, previous)) {
                previous = std::move(candidate);
                continue;
            }
            storage(previous);
            if (divide(a, previous).second.empty() && divide(b, previous).second.empty())
                return previous;
            previous.clear();
        }
    }
    GaussianPolynomial multiply(const GaussianPolynomial &a, const GaussianPolynomial &b) {
        checkpoint();
        if (a.empty() || b.empty())
            return {};
        GaussianPolynomial p(a.size() + b.size() - 1);
        for (size_t i = 0; i < a.size(); ++i) {
            for (size_t j = 0; j < b.size(); ++j) {
                const auto v = product(a[i], b[j]);
                p[i + j] = {add(p[i + j].real, v.real), add(p[i + j].imag, v.imag)};
            }
            storage(p);
        }
        trim(p);
        return p;
    }
};

// F_p(i), p = 2^31-1 (prime, p == 3 mod 4), so x^2+1 is irreducible.
// Admissible denominator and leading-coefficient reductions plus gcd(f,f')=1
// prove squarefreeness over Q(i). A failed test proves NOTHING about repeats;
// it merely selects exact arithmetic. FLINT supplies its fast polynomial GCD;
// the prime, field, admissibility checks and one-sided proof are unchanged.
constexpr uint64_t prime = UINT64_C(2147483647);
uint64_t power(uint64_t x, uint64_t exponent) {
    uint64_t y = 1;
    for (; exponent; exponent >>= 1, x = x * x % prime)
        if (exponent & 1)
            y = y * x % prime;
    return y;
}
// All context-dependent values are scoped together. No borrowed FLINT handle
// escapes, and destruction clears polynomials before their field context.
class ModularSquarefreeTest {
    fq_nmod_ctx_t field_;
    fq_nmod_poly_t polynomial_, derivative_, gcd_;
    nmod_poly_t coefficient_;
  public:
    ModularSquarefreeTest() {
        nmod_poly_init(coefficient_, static_cast<ulong>(prime));
        nmod_poly_set_coeff_ui(coefficient_, 0, 1);
        nmod_poly_set_coeff_ui(coefficient_, 2, 1);
        fq_nmod_ctx_init_modulus(field_, coefficient_, "i"); // i² + 1, p ≡ 3 (mod 4)
        nmod_poly_zero(coefficient_);
        fq_nmod_poly_init(polynomial_, field_);
        fq_nmod_poly_init(derivative_, field_);
        fq_nmod_poly_init(gcd_, field_);
    }
    ~ModularSquarefreeTest() {
        fq_nmod_poly_clear(gcd_, field_);
        fq_nmod_poly_clear(derivative_, field_);
        fq_nmod_poly_clear(polynomial_, field_);
        fq_nmod_ctx_clear(field_);
        nmod_poly_clear(coefficient_);
    }
    ModularSquarefreeTest(const ModularSquarefreeTest &) = delete;
    ModularSquarefreeTest &operator=(const ModularSquarefreeTest &) = delete;
    bool proves(const GaussianPolynomial &p, const cps_options &options) {
        for (size_t i = 0; i < p.size(); ++i) {
            if (options.cancelled()) throw Cancelled{};
            const auto reduce = [](const Rational &q, ulong &out) {
                const auto den = fmpz_fdiv_ui(fmpq_denref(q.value()), static_cast<ulong>(prime));
                if (!den) return false;
                out = static_cast<ulong>(fmpz_fdiv_ui(fmpq_numref(q.value()), static_cast<ulong>(prime)) *
                                         power(den, prime - 2) % prime);
                return true;
            };
            ulong re, im;
            if (!reduce(p[i].real, re) || !reduce(p[i].imag, im)) return false;
            // fq_nmod uses the nmod polynomial basis re + im*i. Degree < 2
            // is already reduced by the explicitly chosen field modulus.
            nmod_poly_set_coeff_ui(coefficient_, 0, re);
            nmod_poly_set_coeff_ui(coefficient_, 1, im);
            fq_nmod_poly_set_coeff(polynomial_, static_cast<slong>(i), coefficient_, field_);
        }
        if (p.empty() || fq_nmod_poly_degree(polynomial_, field_) != static_cast<slong>(p.size()-1))
            return false;
        fq_nmod_poly_derivative(derivative_, polynomial_, field_);
        if (options.cancelled()) throw Cancelled{};
        // One backend call is atomic with respect to cooperative cancellation.
        // A cancellation observed on return discards the proof transactionally.
        fq_nmod_poly_gcd(gcd_, polynomial_, derivative_, field_);
        if (options.cancelled()) throw Cancelled{};
        return fq_nmod_poly_degree(gcd_, field_) == 0;
    }
};
ModularSquarefreeResult modular_squarefree(const GaussianPolynomial &p, const cps_options &options) {
    if (!modular_proof_admitted(p.size(), extension_proof_words, options)) return ModularSquarefreeResult::Refused;
    return ModularSquarefreeTest().proves(p, options) ? ModularSquarefreeResult::Proved
                                                      : ModularSquarefreeResult::Inconclusive;
}
} // namespace

const uint64_t extension_proof_words = 3 * (sizeof(fq_nmod_struct) + 2 * sizeof(ulong)) / sizeof(ulong);

GaussianDecomposition exact_complex_factors(std::span<const Coefficient> coefficients,
                                            const cps_options &options, const EnginePolicy *trace) {
    GaussianDecomposition result;
    result.status = "Q(i) preprocessing not applicable";
    if (coefficients.size() < 2 || std::all_of(coefficients.begin(), coefficients.end(),
                                               [](const Coefficient &c) { return c.is_real(); }))
        return result;
    if (std::any_of(coefficients.begin(), coefficients.end(),
                    [](const Coefficient &c) { return !c.exact(); })) {
        result.status = "Q(i) preprocessing requires exact coefficients";
        return result;
    }
    if (!options.exact_storage_limit() || !options.exact_operations ||
        coefficients.size() > GaussianPolynomial().max_size() ||
        !resources::storage_fits(resources::multiply(coefficients.size(), sizeof(GaussianRational)), options.temporary_bytes)) {
        result.status = "Q(i) preprocessing storage/budget guard";
        return result;
    }
    if (!options.allow_exponent_expansion &&
        !compact_expansion_fits(coefficients, options.exact_expansion_allowance())) {
        // The split prime's test needs no write-out (squarefree_modular reduces each m@e number directly), so it
        // still runs: a pass keeps the square-free proof for the ball path, and anything else says nothing.
        result.status = "compact exponent expansion needs explicit allowance";
        if (modular_proof_admitted(coefficients.size(), split_proof_words, options)) {
            const auto proof = squarefree_modular(coefficients, options);
            result.cancelled = proof == ModularSquarefreeResult::Cancelled;
            result.squarefree_proved = proof == ModularSquarefreeResult::Proved;
            if (result.squarefree_proved)
                result.status += "; squarefreeness proved modulo 998244353 on the coefficients as written";
        }
        return result;
    }
    // Exact zero roots first: p = x^m q with q(0) != 0. Both modular
    // certificates and the exact decomposition then concern q alone, and x^m
    // is an exact factor of known multiplicity, as FLINT's integer squarefree
    // factorization returns it on the real path. Without this split, any zero
    // root makes gcd(p, p') nontrivial, so both certificates fail by
    // construction and the exact Euclidean fallback runs at full degree.
    size_t zero_order = 0;
    while (zero_order + 1 < coefficients.size() && coefficients[zero_order].zero()) ++zero_order;
    coefficients = coefficients.subspan(zero_order);
    const auto with_zero_factor = [&](GaussianDecomposition &decomposition) {
        if (zero_order) decomposition.factors.push_back({GaussianPolynomial{GaussianRational{}, one()}, zero_order});
    };
    bool modular_refused = false;
    try {
        Arithmetic arithmetic(options, trace);
        arithmetic.checkpoint();
        GaussianPolynomial polynomial(coefficients.size());
        uint64_t remaining = options.exact_storage_limit();
        {
        TraceScope phase(trace, "gaussian_expand", 0, coefficients.size());
        for (size_t i = 0; i < coefficients.size(); ++i) {
            if (!coefficients[i].real.rational_with_budget(polynomial[i].real.value(), remaining) ||
                !coefficients[i].imag.rational_with_budget(polynomial[i].imag.value(), remaining))
                throw BudgetExceeded{};
        }
        }
        trim(polynomial);
        arithmetic.storage(polynomial);
        if (polynomial.size() < 2) {
            // A monomial c x^m: only the exact zero root remains.
            with_zero_factor(result);
            if (zero_order) result.status = "Q(i) monomial: exact zero root";
            return result;
        }
        const size_t size = polynomial.size();
        // Below the coefficient count the exact_operations allowance admits neither modular proof
        // (modular_proof_admitted: one unit an imported element) nor the exact decomposition, which charges at
        // least six word operations a coefficient to make the polynomial monic: refused before any stage, so
        // that no proof runs uncharged. A refusal says nothing about the input.
        // History: docs/code-history-20261021.md, exact_complex.cpp 3.
        if (!options.exact_effort_fits(size)) {
            result.status = "Q(i) preprocessing refused: the exact_operations allowance is below the coefficient "
                            "count, which admits neither a modular proof nor the exact decomposition (budget guard)";
            return result;
        }
        // Each modular proof runs only where its residue polynomials fit the temporary storage allowance. A
        // proof refused for that is not inconclusive at its prime: the refusal is kept for the diagnostics.
        auto split = ModularSquarefreeResult::Refused;
        if (modular_proof_admitted(size, split_proof_words, options)) {
            TraceScope phase(trace, "gaussian_split_modular_gcd", 30, size);
            split = squarefree_modular(coefficients.first(size), options);
            if (split == ModularSquarefreeResult::Cancelled) throw Cancelled{};
        }
        arithmetic.checkpoint();
        auto extension = ModularSquarefreeResult::Refused;
        if (split != ModularSquarefreeResult::Proved && modular_proof_admitted(size, extension_proof_words, options)) {
            TraceScope phase(trace, "gaussian_modular_gcd", 31, size);
            extension = modular_squarefree(polynomial, options);
        }
        modular_refused = split == ModularSquarefreeResult::Refused ||
                          (split != ModularSquarefreeResult::Proved && extension == ModularSquarefreeResult::Refused);
        if (split == ModularSquarefreeResult::Proved || extension == ModularSquarefreeResult::Proved) {
            result.factors.push_back({std::move(polynomial), 1});
            with_zero_factor(result);
            result.status = split == ModularSquarefreeResult::Proved
                                ? "Q(i) squarefreeness proved modulo 998244353 with i mapped to F_p"
                                : "Q(i) squarefreeness proved modulo 2147483647";
            return result;
        }
        // Characteristic-zero squarefree decomposition: w=f/gcd(f,f'), then
        // successively remove gcd(w,g). Each exact division is checked.
        TraceScope exact_phase(trace, "gaussian_exact_decomposition", 0, polynomial.size());
        polynomial = arithmetic.monic(std::move(polynomial));
        auto repeated = arithmetic.gcd(polynomial, arithmetic.derivative(polynomial));
        auto remaining_factor = arithmetic.exact_divide(polynomial, repeated);
        std::vector<GaussianFactor> factors;
        for (uint64_t multiplicity = 1; remaining_factor.size() > 1; ++multiplicity) {
            if (multiplicity >= polynomial.size())
                throw std::logic_error("Q(i) multiplicity exceeds degree");
            auto shared = arithmetic.gcd(remaining_factor, repeated);
            auto factor = arithmetic.exact_divide(remaining_factor, shared);
            if (factor.size() > 1)
                factors.push_back({std::move(factor), multiplicity});
            remaining_factor = std::move(shared);
            repeated = arithmetic.exact_divide(repeated, remaining_factor);
        }
        GaussianPolynomial reconstructed{one()};
        for (const auto &factor : factors)
            for (uint64_t i = 0; i < factor.multiplicity; ++i)
                reconstructed = arithmetic.multiply(reconstructed, factor.polynomial);
        if (!equal(reconstructed, polynomial))
            throw std::logic_error("Q(i) factor reconstruction failed");
        result.factors = std::move(factors);
        with_zero_factor(result);
        result.status = "Q(i) squarefree decomposition; exact reconstruction verified";
    } catch (const BudgetExceeded &) {
        result.status = modular_refused ? "Q(i) preprocessing arithmetic/storage budget exhausted; the modular "
                                          "proofs were refused by the temporary_bytes allowance"
                                        : "Q(i) preprocessing arithmetic/storage budget exhausted";
    } catch (const Cancelled &) {
        result.cancelled = true;
        result.status = "Q(i) preprocessing cancelled";
    }
    return result;
}

void enclose_gaussian(acb_poly_t out, const GaussianPolynomial &polynomial, slong precision) {
    acb_poly_zero(out);
    ComplexBall coefficient;
    for (size_t i = 0; i < polynomial.size(); ++i) {
        arb_set_fmpq(acb_realref(coefficient.value()), polynomial[i].real.value(), precision);
        arb_set_fmpq(acb_imagref(coefficient.value()), polynomial[i].imag.value(), precision);
        acb_poly_set_coeff_acb(out, static_cast<slong>(i), coefficient.value());
    }
}
} // namespace certpolysolve::detail
