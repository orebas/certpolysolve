// <certpolysolve/gmp.hpp>: wilk30's exact coefficients built with mpz_class and given as they are; an exact root end
// back as mpq_class; mpq_class coefficients.
#include <certpolysolve/gmp.hpp>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

int main() {
    certpolysolve::ThreadScope lifetime;
    try {
        std::vector<mpz_class> c{1};
        for (int k = 1; k <= 30; ++k) { // multiply by (x - k)
            c.push_back(0);
            for (size_t i = c.size() - 1; i > 0; --i) c[i] = c[i - 1] - k * c[i];
            c[0] *= -k;
        }
        const certpolysolve::Polynomial w30(c);
        if (!(w30 == certpolysolve::Polynomial::from_roots(std::vector<int>{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,
                                                                         17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30})))
            throw std::runtime_error("mpz_class coefficients differ from from_roots");
        const auto result = certpolysolve::real_roots(w30, {.bits = 64});
        if (!result.ok() || result.roots().size() != 30) throw std::runtime_error("wilk30 not solved: " + result.message());
        for (const auto &root : result.roots()) {
            const mpq_class lo = certpolysolve::to_mpq(root.box.real.lo), hi = certpolysolve::to_mpq(root.box.real.hi);
            int inside = 0;
            for (int k = 1; k <= 30; ++k) inside += lo <= k && k <= hi;
            if (inside != 1) throw std::runtime_error("a region does not hold exactly one integer");
        }
        const certpolysolve::Polynomial q{mpq_class(1, 3), mpq_class(-4, 6)};
        if (q[1].value.real() != "-2/3") throw std::runtime_error("mpq_class: " + q[1].value.real());
        // to_mpq: an end's exact value; an exponent beyond long long is a number no memory holds (std::bad_alloc).
        if (certpolysolve::to_mpq(certpolysolve::Dyadic("3", "-2")) != mpq_class(3, 4) ||
            certpolysolve::to_mpq(certpolysolve::Dyadic("-5", "3")) != mpq_class(-40))
            throw std::runtime_error("to_mpq of small ends");
        bool refused = false;
        try {
            static_cast<void>(certpolysolve::to_mpq(certpolysolve::Dyadic("1", "100000000000000000000")));
        } catch (const std::bad_alloc &) {
            refused = true;
        }
        if (!refused) throw std::runtime_error("to_mpq of a 2^(10^20) end did not report bad_alloc");
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
    std::cout << "gmp interop tests passed\n";
    return 0;
}
