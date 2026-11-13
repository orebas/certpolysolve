#pragma once
// Rouche's test on a chart's model (certify_root_model, from goal_half_width's rectangle), the chart's root count
// (count_root_model), what a refused certificate reports (ChartFailure, ChartCertificateDiagnostic) and the
// jet a caller keeps across the certificates at one point (CertificateJet, the jet orders). Defined in
// chart_certificate.cpp; used by the piecewise state and the certificate attempt. No method state is retained.
#include <certpolysolve/certpolysolve.h>
#include "core/cover.hpp"
#include "core/number.hpp"
#include "numeric/models/root_model.hpp"
#include "core/wide_int.hpp"
#include <array>
#include <cstdint>
#include <optional>
#include <vector>

namespace certpolysolve::detail {
struct ChartCertificateStatistics;
struct ChartCertificateWork;
struct ChartCertificateDiagnostic {
    double error_log2 = 0, margin_log2 = 0;
    bool component_terms_available = false;
    std::array<double, 5> component_log2{};
};
enum class ChartFailure { Domain, ModelRemainder, Rounding, Candidate, Width };
// A candidate's Taylor jet and the scratch that lets it grow, kept by a caller across the certificates it proves
// at one point of one model at one precision (certify_root_model's shared_jet). `jet` is the longest computed so
// far; `top`, while `extensible`, the last accumulator's column (certify_root_model's ladder).
struct CertificateJet {
    std::vector<ComplexBall> jet, top, below;
    bool extensible = false;
    void clear() {
        jet.clear();
        extensible = false;
    }
};
// Jet orders of a certificate: a final rectangle (a relative goal, a tiny box) needs a cubic jet; a guard, whose
// rectangle is a fraction of the chart, up to sixteen orders to keep the cancellation near the root. Both bound
// every further order by a majorant.
inline constexpr slong certificate_jet_order = 3, guard_jet_order = 16;
// The remainder bound's binomials C(n, m) without an fmpz where they fit in a word: every C(n, m) with n up to
// word_binomial_degree is below 2^63 (C(66, 33) is the largest), and each step's product C * n in 128 bits. An fmpz
// above 62 bits takes an mpz, and the first in a process initializes FLINT's pool (about 1.1 million instructions, a
// third of a one-shot degree-8 solve of a rounded double root; dev/docs/easy-panel-rerun-20261027.md). arb_mul_ui and
// arb_mul_fmpz of the same value give the same ball: both multiply by the exact arf of it (FLINT 3.6, arb/mul.c).
inline constexpr slong word_binomial_degree = 66;
inline uint128 word_binomial(ulong n, ulong m) {
    if (m > n) return 0;  // as fmpz_bin_uiui
    uint128 c = 1;
    for (ulong i = 1; i <= m; ++i) c = c * (n - m + i) / i;  // C(n - m + i, i), exact at each step
    return c;
}
// The half-width of the rectangle certify_root_model tries first around a centre: the goal's relative width with
// five guard bits, or an eighth of the absolute tolerance where that is larger; the certificate attempt's too
// (certificate_attempt.cpp).
void goal_half_width(arb_t radius, const acb_t center, const cps_options &options, slong precision);
bool certify_root_model(Root &, const RootModel &, const acb_t local_point, const cps_options &,
                        slong precision, ChartFailure *failure = nullptr,
                        ChartCertificateStatistics *statistics = nullptr,
                        double *correction_log2 = nullptr,
                        ChartCertificateDiagnostic *diagnostic = nullptr,
                        slong experimental_jet_order = 0, bool defer_guard_jet = false,
                        ChartCertificateWork *work = nullptr, CertificateJet *shared_jet = nullptr);
// A disk count only guides candidate work. Exported rectangles need their own
// certificate. An empty optional is inconclusive, not a zero count.
std::optional<uint64_t> count_root_model(const RootModel &, slong precision);
} // namespace certpolysolve::detail
