// Owned program construction, dependency/liveness metadata and the PSPROGRAM text reader.
// append checks operands before publishing a node; parse_program reports the first invalid line.
// Called by the C API and program transformations. Evaluation lives in numeric/evaluate.
#include "input/program/program.hpp"
#include <algorithm>
#include <charconv>
#include <limits>
#include <stdexcept>
#include <string>

namespace certpolysolve::detail {

uint64_t Program::append(const ProgramNode &node) {
    const auto size = static_cast<uint64_t>(nodes.size());
    const auto operand = [&](uint64_t index) {
        if (index >= size)
            throw std::invalid_argument("program operand " + std::to_string(index) + " is not an earlier node");
    };
    switch (node.op) {
    case ProgramOp::Variable: break;
    case ProgramOp::Constant:
        if (node.value >= constants.size()) throw std::invalid_argument("program constant out of range");
        break;
    case ProgramOp::Add:
    case ProgramOp::Subtract:
    case ProgramOp::Multiply: operand(node.a); operand(node.b); break;
    case ProgramOp::Negate: operand(node.a); break;
    case ProgramOp::Power:
        operand(node.a);
        if (node.value < 2) throw std::invalid_argument("program power exponent must be at least 2");
        break;
    }
    nodes.push_back(node);
    return size;
}

uint64_t Program::append_constant(Coefficient constant) {
    if (!constant.exact()) throw std::invalid_argument("program constants must be exact");
    constants.push_back(std::move(constant));
    ProgramNode node;
    node.op = ProgramOp::Constant;
    node.value = constants.size() - 1;
    return append(node);
}

std::vector<uint64_t> Program::last_uses() const {
    const auto none = static_cast<uint64_t>(nodes.size());
    std::vector<uint64_t> last(nodes.size(), none);
    for (uint64_t i = 0; i < nodes.size(); ++i) {
        const auto &node = nodes[i];
        if (node.op == ProgramOp::Variable || node.op == ProgramOp::Constant) continue;
        last[node.a] = i;
        if (node.op != ProgramOp::Negate && node.op != ProgramOp::Power) last[node.b] = i;
    }
    return last;
}

std::vector<uint64_t> Program::degree_bounds() const {
    constexpr uint64_t top = std::numeric_limits<uint64_t>::max();
    const auto sum = [](uint64_t a, uint64_t b) { return a > top - b ? top : a + b; };
    const auto product = [](uint64_t a, uint64_t b) { return a && b > top / a ? top : a * b; };
    std::vector<uint64_t> degree(nodes.size(), 0);
    for (size_t i = 0; i < nodes.size(); ++i) {
        const auto &node = nodes[i];
        switch (node.op) {
        case ProgramOp::Variable: degree[i] = 1; break;
        case ProgramOp::Constant: degree[i] = 0; break;
        case ProgramOp::Add:
        case ProgramOp::Subtract: degree[i] = std::max(degree[node.a], degree[node.b]); break;
        case ProgramOp::Multiply: degree[i] = sum(degree[node.a], degree[node.b]); break;
        case ProgramOp::Negate: degree[i] = degree[node.a]; break;
        case ProgramOp::Power: degree[i] = product(degree[node.a], node.value); break;
        }
    }
    return degree;
}

namespace {
std::string_view trim(std::string_view s) {
    while (!s.empty() && (s.front() == ' ' || s.front() == '\t' || s.front() == '\r')) s.remove_prefix(1);
    while (!s.empty() && (s.back() == ' ' || s.back() == '\t' || s.back() == '\r')) s.remove_suffix(1);
    return s;
}
std::vector<std::string_view> fields(std::string_view line) {
    std::vector<std::string_view> out;
    while (!(line = trim(line)).empty()) {
        const auto end = line.find_first_of(" \t");
        out.push_back(line.substr(0, end));
        if (end == std::string_view::npos) break;
        line.remove_prefix(end);
    }
    return out;
}
uint64_t index_of(std::string_view text) {
    uint64_t value = 0;
    const auto [end, error] = std::from_chars(text.data(), text.data() + text.size(), value);
    if (error != std::errc{} || end != text.data() + text.size())
        throw std::invalid_argument("expected a node index or an exponent, found '" + std::string(text) + "'");
    return value;
}
} // namespace

Program parse_program(std::string_view text) {
    Program program;
    bool header = false;
    size_t number = 0;
    while (!text.empty()) {
        const auto end = text.find('\n');
        auto line = text.substr(0, end);
        text.remove_prefix(end == std::string_view::npos ? text.size() : end + 1);
        ++number;
        if (const auto comment = line.find('#'); comment != std::string_view::npos) line = line.substr(0, comment);
        const auto f = fields(line);
        if (f.empty()) continue;
        try {
            if (!header) {
                if (f.size() != 2 || f[0] != "PSPROGRAM" || f[1] != "1")
                    throw std::invalid_argument("expected the header 'PSPROGRAM 1'");
                header = true;
                continue;
            }
            const auto arity = [&](size_t n) {
                if (f.size() != n + 1)
                    throw std::invalid_argument("'" + std::string(f[0]) + "' takes " + std::to_string(n) + " fields");
            };
            ProgramNode node;
            if (f[0] == "x") {
                arity(0);
                node.op = ProgramOp::Variable;
                program.append(node);
            } else if (f[0] == "c") {
                if (f.size() != 2 && f.size() != 3) throw std::invalid_argument("'c' takes a real and an optional imaginary part");
                Coefficient c;
                c.real = Number(f[1]);
                c.imag = f.size() == 3 ? Number(f[2]) : Number("0");
                program.append_constant(std::move(c));
            } else if (f[0] == "add" || f[0] == "sub" || f[0] == "mul") {
                arity(2);
                node.op = f[0] == "add" ? ProgramOp::Add : f[0] == "sub" ? ProgramOp::Subtract : ProgramOp::Multiply;
                node.a = index_of(f[1]);
                node.b = index_of(f[2]);
                program.append(node);
            } else if (f[0] == "neg") {
                arity(1);
                node.op = ProgramOp::Negate;
                node.a = index_of(f[1]);
                program.append(node);
            } else if (f[0] == "pow") {
                arity(2);
                node.op = ProgramOp::Power;
                node.a = index_of(f[1]);
                node.value = index_of(f[2]);
                program.append(node);
            } else if (f[0] == "out") {
                arity(1);
                const auto output = index_of(f[1]);
                if (output >= program.nodes.size()) throw std::invalid_argument("output is not a node");
                program.output = output;
            } else
                throw std::invalid_argument("unknown operation '" + std::string(f[0]) + "'");
        } catch (const std::exception &error) {
            throw std::invalid_argument("program line " + std::to_string(number) + ": " + error.what());
        }
    }
    if (!header) throw std::invalid_argument("program text has no 'PSPROGRAM 1' header");
    if (!program.output) throw std::invalid_argument("program text has no 'out' line");
    return program;
}

} // namespace certpolysolve::detail
