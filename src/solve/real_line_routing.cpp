// The all-real seed pass (real_line_seeds.cpp) under a call's options. real_line_seed_limits gives its limits (the
// precision cap, temporary storage, cell share and cancellation, borrowing the options) to a factor's start
// (solver.cpp, seed_from_real_line) and to real_line_check, Auto's check before routing (session.cpp): the pass on
// the primitive part of the exact real coefficients' integer import (streamed_real_polynomial), as a
// RealLineSeedHint whatever it found (none below heuristics::real_line_seed_min_degree); every_root_seeded keeps it
// only when it seeds every root. State: none. Design notes: docs/real-line-seed-pass-design-20261020.md (the pass),
// docs/allroots-seed-routing-note-20261020.md (the check), dev/docs/slowdowns-20261029.md (the kept stopped pass).
#include "solve/real_line_routing.hpp"
#include "numeric/exact/exact_preprocessing.hpp"
#include "core/heuristics.hpp"
#include <algorithm>

namespace certpolysolve::detail {
RealLineSeedLimits real_line_seed_limits(const cps_options &o) {
    RealLineSeedLimits limits;
    limits.max_precision = static_cast<slong>(std::min<uint64_t>(o.max_work_bits, static_cast<uint64_t>(WORD_MAX)));
    limits.max_bytes = o.temporary_bytes;
    limits.min_cell_fraction = heuristics::real_line_seed_min_cell_fraction;
    limits.cancelled = [&o] { return o.cancelled(); };
    return limits;
}

std::optional<RealLineSeedHint> real_line_check(std::span<const Coefficient> coefficients, const cps_options &o) {
    RealLineSeedHint hint;
    // The nonzero part, which the solver prepares once it has set the exact zero roots apart (SolveState::Impl::
    // preprocess), so that a square-free nonzero part's one factor equals the hint.
    size_t zeros = 0;
    while (zeros + 1 < coefficients.size() && coefficients[zeros].zero()) ++zeros;
    // The solver's own import (exact_preprocessing.hpp), under the same budgets its preparation charges.
    if (!streamed_real_polynomial(hint.polynomial, coefficients.subspan(zeros), o)) return std::nullopt;
    // Primitive with a positive leading coefficient: the form a square-free factor takes (fmpz_poly_factor_squarefree),
    // so that a square-free input's one factor equals it.
    fmpz_poly_primitive_part(hint.polynomial.value(), hint.polynomial.value());
    const slong degree = fmpz_poly_degree(hint.polynomial.value());
    if (degree < heuristics::real_line_seed_min_degree) return std::nullopt;
    hint.pass = real_line_seeds(hint.polynomial.value(), real_line_seed_limits(o));
    return hint;
}

bool seeds_every_root(const RealLineSeedHint &hint) {
    return hint.pass.statistics.stop == RealLineSeedStop::Seeds &&
           static_cast<slong>(hint.pass.seeds.size()) == fmpz_poly_degree(hint.polynomial.value());
}

std::optional<RealLineSeedHint> every_root_seeded(std::span<const Coefficient> coefficients, const cps_options &o) {
    auto hint = real_line_check(coefficients, o);
    if (!hint || !seeds_every_root(*hint)) return std::nullopt;
    return hint;
}
} // namespace certpolysolve::detail
