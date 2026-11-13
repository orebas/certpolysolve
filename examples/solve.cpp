#include <certpolysolve/certpolysolve.hpp>
#include <algorithm>
#include <charconv>
#include <format>
#include <iostream>
#include <optional>
#include <string>
#include <vector>
#include "secular_text.hpp"

int main(int argc, char **argv) {
    certpolysolve::ThreadScope thread_lifetime;
    try {
        certpolysolve::Polynomial polynomial;
        certpolysolve::Options options;
        std::optional<certpolysolve::SecularEquation> secular;
        bool real_only = false, highest_first = false;
        std::vector<std::string> coefficients;
        for (int i = 1; i < argc; ++i) {
            const std::string argument(argv[i]);
            if (argument == "--real") real_only = true;
            else if (argument == "--method" && i + 1 < argc) {
                const std::string method(argv[++i]);
                if (method == "auto") options.method = certpolysolve::Method::automatic;
                else if (method == "classical") options.method = certpolysolve::Method::classical;
                else if (method == "secular") options.method = certpolysolve::Method::secular;
                else if (method == "piecewise") options.method = certpolysolve::Method::piecewise;
                else throw std::invalid_argument("--method expects auto, classical, secular or piecewise");
            }
            else if (argument == "--approximation") options.goal = certpolysolve::Goal::approximation();
            else if (argument == "--clusters" && i + 1 < argc) options.goal = certpolysolve::Goal::clusters(certpolysolve::Number(argv[++i]));
            else if (argument == "--highest-first" || argument == "--descending") highest_first = true;
            else if (argument == "--secular" && i + 1 < argc) {
                if (secular) throw std::invalid_argument("only one secular input is allowed");
                std::vector<certpolysolve::SecularEquation::Term> terms;
                for (const auto &t : read_secular_terms(std::string(argv[++i])))
                    terms.push_back({.residue = certpolysolve::Number::complex(certpolysolve::Number(t.residue_real),
                                                                           certpolysolve::Number(t.residue_imag)),
                                     .pole = certpolysolve::Number::complex(certpolysolve::Number(t.pole_real),
                                                                        certpolysolve::Number(t.pole_imag))});
                secular.emplace(terms);
                if (!secular->ok()) throw std::invalid_argument(secular->status());
            }
            else if (argument == "--bits" && i + 1 < argc) {
                const std::string text(argv[++i]);
                uint64_t bits = 0;
                const auto parsed = std::from_chars(text.data(), text.data() + text.size(), bits);
                if (parsed.ec != std::errc{} || parsed.ptr != text.data() + text.size() || bits > 1000000)
                    throw std::invalid_argument("--bits expects an integer in 1..1,000,000");
                options.bits = bits;
                options.max_bits = std::min<uint64_t>(1000000, std::max<uint64_t>(4096, bits + 64));
            } else if (argument == "--help") {
                std::cout << "Usage: certpolysolve_cli [--real] [--approximation | --clusters R] [--bits N] [--highest-first] a0 a1 ... an\n"
                             "       certpolysolve_cli [--real] [--approximation | --clusters R] [--bits N] --secular FILE\n"
                             "Coefficients with the constant term first (x^2 - 2 is: -2 0 1), or with the highest power\n"
                             "first after --highest-first (1 0 -2); exact decimal, rational, or mantissa@binary_exponent.\n"
                             "The polynomial as read is echoed on stderr.\n"
                             "--approximation permits sufficiently accurate counted clusters; --clusters R asks for boxes of\n"
                             "radius at most R (a number, e.g. 1e-6), each with its count.\n"
                             "--method auto|classical|secular|piecewise selects a solver.\n"
                             "Secular files require modern Secular; Degree=N; headers.\n";
                return 0;
            } else coefficients.push_back(argument);
        }
        if (coefficients.empty() == !secular) throw std::invalid_argument("supply coefficients or one secular file; use --help");
        if (!coefficients.empty()) {
            const std::vector<certpolysolve::Number> numbers(coefficients.begin(), coefficients.end());
            polynomial = highest_first ? certpolysolve::Polynomial::highest_first(numbers) : certpolysolve::Polynomial::lowest_first(numbers);
            if (!polynomial.ok()) throw std::invalid_argument(polynomial.status());
            std::cerr << "p(x) = " << std::format("{}", polynomial) << '\n';
        }
        auto result = secular ? (real_only ? certpolysolve::real_roots(*secular, options) : certpolysolve::solve(*secular, options))
                              : (real_only ? certpolysolve::real_roots(polynomial, options) : certpolysolve::solve(polynomial, options));
        if (result.stop() == certpolysolve::Stop::error) throw std::runtime_error(result.message());
        std::cout << "degree=" << result.degree() << " complete=" << result.complete()
                  << " goal_met=" << result.goal_met() << " termination=" << cps_result_termination(result.get())
                  << " outer_work_bits=" << result.work_bits() << '\n';
        auto endpoint = [](const certpolysolve::Dyadic& d) { return d.mantissa() + " * 2^(" + d.exponent() + ")"; };
        for (const auto& root : result.roots()) {
            std::cout << "real=[" << endpoint(root.box.real.lo) << ", " << endpoint(root.box.real.hi)
                      << "] imag=[" << endpoint(root.box.imag.lo) << ", " << endpoint(root.box.imag.hi)
                      << "] count=" << root.count << " multiplicity=" << root.multiplicity.value_or(0)
                      << " certificate=" << root.certificate << '\n';
        }
        std::cout << result.message() << '\n';
        return result.goal_met() ? 0 : 2;
    } catch (const std::exception& error) {
        std::cerr << "certpolysolve: " << error.what() << '\n';
        return 1;
    }
}
