#pragma once
// Optional interop with GMP's C++ classes (link gmpxx and gmp): mpz_class and mpq_class are Numbers, so they go into
// coefficient lists and from_roots as they are, and an exact root end converts to mpq_class.
//   std::vector<mpz_class> c = ...;  certpolysolve::Polynomial p(c);
//   mpq_class lo = certpolysolve::to_mpq(result.roots()[0].box.real.lo);
#include <certpolysolve/certpolysolve.hpp>
#include <gmpxx.h>
#include <charconv>
#include <new>
#include <string>

namespace certpolysolve {
template <> struct number_traits<mpz_class> {
    static std::string text(const mpz_class &v) { return v.get_str(); }
};
template <> struct number_traits<mpq_class> {
    static std::string text(const mpq_class &v) {
        mpq_class canonical(v);
        canonical.canonicalize();
        return canonical.get_str();
    }
};
// The exact value of a root end, m 2^e (an end as a result gives it: decimal integers). An exponent beyond long long
// is a number of more than 2^63 bits, which no memory holds: std::bad_alloc, as for any other allocation that fails.
inline mpq_class to_mpq(const Dyadic &d) {
    mpq_class q(mpz_class(d.mantissa()), 1);
    long long e = 0;
    const std::string &text = d.exponent();
    const auto read = std::from_chars(text.data(), text.data() + text.size(), e);
    if (read.ec == std::errc::result_out_of_range) throw std::bad_alloc();
    if (read.ec != std::errc{} || read.ptr != text.data() + text.size())
        throw std::invalid_argument("a root end's exponent is not a decimal integer");
    const unsigned long long magnitude = e >= 0 ? static_cast<unsigned long long>(e) : 0ULL - static_cast<unsigned long long>(e);
    if (magnitude > static_cast<unsigned long long>(static_cast<mp_bitcnt_t>(-1))) throw std::bad_alloc();
    if (e >= 0) mpz_mul_2exp(q.get_num_mpz_t(), q.get_num_mpz_t(), static_cast<mp_bitcnt_t>(magnitude));
    else mpz_mul_2exp(q.get_den_mpz_t(), q.get_den_mpz_t(), static_cast<mp_bitcnt_t>(magnitude));
    q.canonicalize();
    return q;
}
} // namespace certpolysolve
