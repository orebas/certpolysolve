#pragma once
#include <cstdint>

#ifndef CERTPOLYSOLVE_EXPERIMENT_DENSE_WORK_LIMIT
#define CERTPOLYSOLVE_EXPERIMENT_DENSE_WORK_LIMIT (1ULL << 29)
#endif

namespace certpolysolve::detail {
// Experimental controls: deliberately absent from the installed ABI.
enum class IterationMethod {
    Automatic,
    Aberth,
    Secular,
    DurandKerner,
    StructuredQR,
    StructuredQZ,
    JenkinsTraub,
    Laguerre,
    DenseEigen
};
enum class ArithmeticTier {
    Automatic,
    Binary64,
    ScaledBinary64,
    DoubleDouble,
    Adaptive,
    Binary128,
    NFloat,
    MPFR
};
#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
enum class RealMethod { Automatic, Descartes, Sturm }; // the Portfolio engine's real-only choice
// The research build's evaluator of a coefficient source (PolynomialSourceOptions::evaluation): the coefficients, or
// the local Taylor cache that CandidateEngine::Series selects (LocalPolynomialEvaluator). The default build
// evaluates the coefficients. History: docs/code-history-20261021.md, internal.hpp 4.
enum class EvaluationMethod { Coefficients, TaylorCache };
#endif
#ifdef CERTPOLYSOLVE_BUILD_RESEARCH
// The research build's Portfolio engine's controls and the names of the methods and tiers it reports.
namespace controller_config {
// Frozen with benchmarks/campaign-c3/protocol.json before numerical tuning.
inline constexpr unsigned sweeps_per_packet = 2;
inline constexpr unsigned stagnant_packets = 3;
inline constexpr double improvement_bits = 0.5;
inline constexpr unsigned subdivision_period = 4;
// Retention policy only: eviction never makes an evaluator ineligible.
inline constexpr uint64_t jit_retained_cache_bytes = 64ULL * 1024 * 1024;
// Cubic-work experiment admission; independent of caller storage restrictions.
inline constexpr uint64_t dense_work_limit = CERTPOLYSOLVE_EXPERIMENT_DENSE_WORK_LIMIT;
} // namespace controller_config
const char *name(IterationMethod);
const char *name(ArithmeticTier);
#endif
} // namespace certpolysolve::detail
