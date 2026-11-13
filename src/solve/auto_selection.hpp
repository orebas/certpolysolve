#pragma once
#include <optional>
#include "input/coefficients/polynomial_input.hpp"
#include "core/query.hpp"

namespace certpolysolve::detail {
// Auto route from the selector campaign (docs/selector-campaign-20260928.md),
// the public Auto.
// History: docs/code-history-20261021.md, auto_selection.hpp 1.
//
// Queries that the explicit methods support (exact coefficients, a complex
// approximation or accuracy goal) go to the maintained secular recipe, or to
// the maintained piecewise recipe when piecewise's own work estimate is below
// the global candidate stage's minimum sweep budget. Every other query, and
// native secular input, keeps the D5 route. No names, hashes or family labels.
enum class AutoRoute { D5, Secular, Piecewise };

struct AutoSelection {
    AutoRoute route = AutoRoute::D5;
    const char *reason = "";
    uint64_t elapsed_ns = 0;
    // With route Piecewise: the piecewise recipe begins with the certificate
    // attempt the global recipes begin with (PiecewiseRootConfig::
    // certificate_attempt). See select_auto_route for the rule.
    bool attempt_first = false;
    // With route Piecewise: piecewise's predicted cost, its admission's work in complex coefficient operations
    // (exact), which bounds the attempt before it (PiecewiseRootConfig::certificate_attempt_fallback).
    uint64_t piecewise_operations = 0;
};

// prove_squarefree false leaves out piecewise admission's square-free proof (PiecewiseAdmissionPolicy): the route
// then follows its estimate alone.
AutoSelection select_auto_route(const DensePolynomial &, const cps_options &, bool prove_squarefree = true);
// The same selection from the input object's terms where it needs no coefficient list (the query outside the
// explicit methods, inexact coefficients, a cap below piecewise's start, or an admission its structure refuses: a
// power composite's), or nothing where it does (C40's second step).
std::optional<AutoSelection> select_auto_route(const cps_polynomial &, const cps_options &, bool prove_squarefree = true);

inline const char *name(AutoRoute route) {
    return route == AutoRoute::Secular ? "secular" : route == AutoRoute::Piecewise ? "piecewise" : "d5";
}
} // namespace certpolysolve::detail
