// Reading the answer: a proof that a closed loop is stable, or that it is not.
// The loop (s + 1)^3 closed with gain K has the characteristic polynomial s^3 + 3s^2 + 3s + 1 + K.
#include <certpolysolve/certpolysolve.hpp>
#include <print>
#include <string_view>

std::string_view verdict(double gain) {
    certpolysolve::Polynomial p{1 + gain, 3, 3, 1}; // s^3 + 3s^2 + 3s + 1 + K
    certpolysolve::Result result = certpolysolve::solve(p);
    if (!result.ok()) return "no answer";

    // Each box holds roots, and together the boxes hold every root, so this is a proof. Its ends are exact numbers,
    // and comparing them with 0 is exact too.
    bool all_left = true;
    for (const certpolysolve::Root &root : result.roots()) {
        if (root.box.real.lo > 0) return "unstable (a root is proved right of the axis)";
        if (root.box.real.hi >= 0) all_left = false; // this box touches the axis
    }
    return all_left ? "stable (every root is proved left of the axis)" : "undecided (a root box touches the axis)";
}

int main() {
    for (double gain : {7.9, 7.999999, 8.0, 8.1}) std::println("K = {}: {}", gain, verdict(gain));
}
