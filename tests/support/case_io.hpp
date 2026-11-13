#pragma once
#include <array>
#include <charconv>
#include <fstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

// Deliberately small process protocol, not a public library serialization API.
namespace corpus_io {
inline int integer(std::string_view text, int minimum, int maximum) {
    int value = 0;
    const auto parsed = std::from_chars(text.data(), text.data() + text.size(), value);
    if (parsed.ec != std::errc{} || parsed.ptr != text.data() + text.size() || value < minimum || value > maximum)
        throw std::invalid_argument("integer option out of range");
    return value;
}
inline bool digits(std::string_view text) {
    if (text.empty()) return false;
    for (char c : text) if (c < '0' || c > '9') return false;
    return true;
}
inline bool rational(std::string_view text) {
    if (text.empty() || text.size() > 100000) return false;
    if (text.front() == '-') text.remove_prefix(1);
    // The protocol uses normalized decimal integers. In particular, forbid
    // leading zeros: Boost's string constructor otherwise recognizes octal,
    // while GMP with explicit base 10 does not. Never let checkers disagree
    // silently about the meaning of a replayed coefficient or endpoint.
    auto decimal = [](std::string_view value) {
        return digits(value) && (value.size() == 1 || value.front() != '0');
    };
    const auto slash = text.find('/');
    if (slash == std::string_view::npos) return decimal(text);
    const auto denominator = text.substr(slash + 1);
    return decimal(text.substr(0, slash)) && decimal(denominator) && denominator != "0";
}
using Coefficients = std::vector<std::array<std::string, 4>>;
inline Coefficients read_fields(std::istream& stream) {
    std::string magic, version, count;
    if (!(stream >> magic >> version >> count) || magic != "PSCASE" || version != "1")
        throw std::invalid_argument("expected PSCASE 1 input");
    Coefficients result(static_cast<size_t>(integer(count, 1, 10001)));
    for (auto& row : result) for (auto& value : row)
        if (!(stream >> value) || !rational(value)) throw std::invalid_argument("invalid rational component");
    return result;
}
inline Coefficients read(const std::string& path) {
    std::ifstream stream(path);
    auto result=read_fields(stream);
    std::string extra;
    if (stream >> extra) throw std::invalid_argument("trailing PSCASE data");
    return result;
}
inline std::string quote(std::string_view text) {
    std::string result = "\"";
    for (unsigned char c : text) {
        if (c == '\\' || c == '"') { result += '\\'; result += static_cast<char>(c); }
        else if (c == '\n') result += "\\n";
        else if (c == '\r') result += "\\r";
        else if (c == '\t') result += "\\t";
        else if (c < 32) result += '?';
        else result += static_cast<char>(c);
    }
    return result + '"';
}
} // namespace corpus_io
