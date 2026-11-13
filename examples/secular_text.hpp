#pragma once
// A secular-equation file (MPSolve's modern secular subset) read into exact terms, without any CertPolySolve API.
#include <algorithm>
#include <cctype>
#include <charconv>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <set>
#include <sstream>
#include <string_view>
#include <limits>

inline bool secular_literal_matches(std::string_view value, std::string_view kind) {
    // Scan iteratively: a recursive regular-expression matcher can exhaust
    // the stack on a valid, long coefficient well within the library's limit.
    size_t cursor = 0;
    auto sign = [&] {
        if (cursor < value.size() && (value[cursor] == '+' || value[cursor] == '-')) ++cursor;
    };
    auto digits = [&] {
        const auto first = cursor;
        while (cursor < value.size() && value[cursor] >= '0' && value[cursor] <= '9') ++cursor;
        return cursor != first;
    };
    sign();
    const bool integral = digits();
    if (kind == "integer") return integral && cursor == value.size();
    if (kind == "rational") {
        if (!integral) return false;
        if (cursor < value.size() && value[cursor] == '/') {
            ++cursor;
            sign();
            if (!digits()) return false;
        }
        return cursor == value.size();
    }
    bool fractional = false;
    if (cursor < value.size() && value[cursor] == '.') {
        ++cursor;
        fractional = digits();
    }
    if (!integral && !fractional) return false;
    if (cursor < value.size() && (value[cursor] == 'e' || value[cursor] == 'E')) {
        ++cursor;
        sign();
        if (!digits()) return false;
    }
    return cursor == value.size();
}

// One term of a secular file, as exact text in the coefficient grammar: residue/(x - pole).
struct SecularTermText {
    std::string residue_real, residue_imag = "0", pole_real = "0", pole_imag = "0";
};
// Explicit modern MPSolve secular subset; no representation guessing or input
// program execution. The library itself accepts owned numeric terms, not files:
// callers build their equation from these terms.
inline std::vector<SecularTermText> read_secular_terms(std::istream &stream) {
    std::string text, line;
    while (std::getline(stream, line))
        text += line.substr(0, line.find('!')) + '\n';
    std::set<std::string> seen;
    bool declared = false, complex = true;
    int degree = -1;
    std::string kind = "floatingpoint";
    auto trim = [](std::string value) {
        const auto first = value.find_first_not_of(" \r\n\t");
        if (first == std::string::npos)
            return std::string{};
        return value.substr(first, value.find_last_not_of(" \r\n\t") - first + 1);
    };
    while (text.find(';') != std::string::npos) {
        const auto end = text.find(';');
        auto option = trim(text.substr(0, end));
        text.erase(0, end + 1);
        std::transform(option.begin(), option.end(), option.begin(),
                       [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
        auto key = option;
        const auto equal = option.find('=');
        if (equal != std::string::npos) {
            key = trim(option.substr(0, equal));
            const auto value = trim(option.substr(equal + 1));
            unsigned count = 0;
            const auto parsed = std::from_chars(value.data(), value.data() + value.size(), count);
            if (parsed.ec != std::errc{} || parsed.ptr != value.data() + value.size())
                throw std::invalid_argument("invalid secular header integer");
            if (key == "degree" && count <= static_cast<unsigned>(std::numeric_limits<int>::max()))
                degree = static_cast<int>(count);
            else if (!(key == "precision" && count == 0))
                throw std::invalid_argument("unsupported secular degree/precision/header");
        } else if (option == "secular") {
            declared = true;
            key = "representation";
        } else if (option == "real" || option == "complex") {
            complex = option == "complex";
            key = "domain";
        } else if (option == "integer" || option == "rational" || option == "floatingpoint") {
            kind = option;
            key = "kind";
        } else if (option == "dense")
            key = "density";
        else
            throw std::invalid_argument("unsupported secular header option");
        if (!seen.insert(key).second)
            throw std::invalid_argument("duplicate secular header option");
    }
    if (!declared || degree < 0)
        throw std::invalid_argument("secular file requires Secular; and Degree=N;");
    std::istringstream values(text);
    auto component = [&] {
        std::string value;
        if (!(values >> value))
            throw std::invalid_argument("truncated secular data");
        if (!secular_literal_matches(value, kind))
            throw std::invalid_argument("literal does not match secular numeric declaration");
        return value;
    };
    std::vector<SecularTermText> terms;
    for (int i = 0; i < degree; ++i) {
        SecularTermText term;
        term.residue_real = component();
        if (complex)
            term.residue_imag = component();
        term.pole_real = component();
        if (complex)
            term.pole_imag = component();
        terms.push_back(std::move(term));
    }
    if (values >> line)
        throw std::invalid_argument("trailing secular data");
    return terms;
}

inline std::vector<SecularTermText> read_secular_terms(const std::string &path) {
    std::ifstream stream(path);
    if (!stream) throw std::invalid_argument("cannot open secular input");
    return read_secular_terms(stream);
}
