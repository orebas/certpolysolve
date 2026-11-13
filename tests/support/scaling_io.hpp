#pragma once
#include "case_io.hpp"
#include <chrono>
#include <iostream>
#include <certpolysolve/certpolysolve.h>
#include <sstream>
#include <stdexcept>

namespace scaling_io {
// Explicit private extension. PSCASE 1 remains strict and unchanged elsewhere.
// Version 2 admits exact binary mantissa@exponent fields and up to 2^20 terms.
inline corpus_io::Coefficients read(const std::string &path) {
    std::ifstream stream(path);
    std::string magic, version, count;
    if (!(stream >> magic >> version >> count) || magic != "PSCASE" || (version != "1" && version != "2"))
        throw std::invalid_argument("expected PSCASE 1 or 2");
    corpus_io::Coefficients result(static_cast<size_t>(corpus_io::integer(count, 1, 1048576)));
    auto number = [&](const std::string &s) {
        if (corpus_io::rational(s))
            return true;
        if (version != "2" || s.size() > 100000)
            return false;
        const auto at = s.find('@');
        return at != s.npos && s.find('/') == s.npos && corpus_io::rational(s.substr(0, at)) &&
               corpus_io::rational(s.substr(at + 1));
    };
    for (auto &row : result)
        for (auto &field : row)
            if (!(stream >> field) || !number(field))
                throw std::invalid_argument("invalid exact component");
    std::string extra;
    if (stream >> extra)
        throw std::invalid_argument("trailing coefficient data");
    return result;
}
// A solve's answer as one JSON line, read through the C API, so a worker on either C++ header (or on the C API
// alone) prints it alike: the fields and their order as they were when this took the pre-0.5 header's RootResult.
inline void result(const cps_result *result, const std::string &query, uint64_t elapsed,
                   const std::string &extra = "", std::ostream &stream = std::cout) {
    if (!result) throw std::invalid_argument("no result to export");
    const auto start = std::chrono::steady_clock::now();
    std::ostringstream output;
    output << "{\"schema\":1,\"query\":" << corpus_io::quote(query) << ",\"degree\":" << cps_result_degree(result)
           << ",\"goal\":" << static_cast<int>(cps_result_goal(result))
           << ",\"complete\":" << (cps_result_complete(result) ? "true" : "false")
           << ",\"goal_met\":" << (cps_result_goal_met(result) ? "true" : "false")
           << ",\"work_bits\":" << cps_result_work_bits(result)
           << ",\"termination\":" << static_cast<int>(cps_result_termination(result)) << ",\"elapsed_ns\":" << elapsed
           << ",\"diagnostics\":" << corpus_io::quote(cps_result_diagnostics(result)) << extra << ",\"regions\":[";
    auto endpoint = [](const cps_dyadic_view &d) {
        return corpus_io::quote(std::string(d.mantissa) + "@" + d.exponent);
    };
    for (size_t i = 0; i < cps_result_size(result); ++i) {
        cps_root_view root{};
        if (cps_result_get_root(result, i, &root) != CPS_OK) throw std::runtime_error(cps_last_error());
        if (i)
            output << ',';
        output << "{\"count\":" << root.count << ",\"multiplicity\":" << root.multiplicity << ",\"real\":["
               << endpoint(root.real_lower) << ',' << endpoint(root.real_upper) << "],\"imag\":["
               << endpoint(root.imag_lower) << ',' << endpoint(root.imag_upper)
               << "],\"certificate\":" << corpus_io::quote(root.certificate) << '}';
    }
    const auto export_ns =
        std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - start)
            .count();
    output << "],\"export_ns\":" << export_ns << "}\n";
    stream << output.str();
}
} // namespace scaling_io
