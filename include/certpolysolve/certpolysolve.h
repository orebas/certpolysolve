#ifndef CERTPOLYSOLVE_CERTPOLYSOLVE_H
#define CERTPOLYSOLVE_CERTPOLYSOLVE_H

#include <stddef.h>
#include <stdint.h>

#if defined(_WIN32)
# if defined(CERTPOLYSOLVE_BUILDING_LIBRARY)
#  define CPS_API __declspec(dllexport)
# else
#  define CPS_API __declspec(dllimport)
# endif
#else
# define CPS_API __attribute__((visibility("default")))
#endif

#ifdef __cplusplus
extern "C" {
#endif

#define CPS_ABI_VERSION 2u
/* The version of these headers; cps_version() reports the library's. No compatibility of any kind is promised
   before 1.0.0. */
#define CPS_VERSION_MAJOR 0
#define CPS_VERSION_MINOR 5
#define CPS_VERSION_PATCH 1
#define CPS_VERSION_STRING "0.5.1"
typedef struct cps_polynomial cps_polynomial;
typedef struct cps_options cps_options;
typedef struct cps_result cps_result;
typedef struct cps_secular cps_secular;
typedef struct cps_program cps_program;

/* Advanced: most callers leave these unspecified. Optional per-stage
   estimates, NOT peak-RSS or wall-time guarantees.
   UINT64_MAX means unspecified (no fixed policy ceiling); zero is an explicit
   zero allowance. Initialize with cps_resource_limits_default before use.
   Default exponent policy retains compact input when expanded coefficient
   storage would exceed the owned input representation. The flag explicitly
   permits that expansion, subject to exact_materialization_bits. */
#define CPS_RESOURCE_ALLOW_EXPONENT_EXPANSION UINT64_C(1)
typedef struct cps_resource_limits {
    uint64_t exact_materialization_bits;
    uint64_t exact_operations;
    uint64_t temporary_bytes;
    uint64_t flags;
} cps_resource_limits;
/* Every allowance unspecified, no flag. NULL does nothing. */
CPS_API void cps_resource_limits_default(cps_resource_limits *);

/* S(x) = sum residue_i/(x-pole_i) - 1. Strings are borrowed during create
   only and denote exact values in the coefficient literal grammar; a NULL
   imaginary part is 0 (the real parts must not be NULL). */
typedef struct cps_secular_term_view {
    const char *residue_real, *residue_imag, *pole_real, *pole_imag;
} cps_secular_term_view;

typedef enum cps_error {
    CPS_OK = 0, CPS_INVALID_ARGUMENT = 1, CPS_OUT_OF_MEMORY = 2,
    CPS_INTERNAL_ERROR = 3, CPS_OUT_OF_RANGE = 4, CPS_UNRESOLVED = 5
} cps_error;
typedef enum cps_termination {
    CPS_SUCCESS = 0, CPS_PRECISION_LIMIT = 1, CPS_CANCELLED = 2,
    CPS_ZERO_POLYNOMIAL = 3, CPS_INDETERMINATE_DEGREE = 4,
    CPS_UNSUPPORTED_QUERY = 5, CPS_WORK_LIMIT = 6
} cps_termination;
typedef enum cps_realness { CPS_REAL_UNKNOWN = 0, CPS_REAL = 1, CPS_NONREAL = 2 } cps_realness;
/* CPS_GOAL_DEFAULT: the options' initial goal, resolved per input kind at the solve: CPS_ROOT_ACCURACY for coefficient
   and secular input, CPS_ROOT_APPROXIMATION for a program. cps_result_goal reports the resolved goal. */
typedef enum cps_goal {
    CPS_ROOT_ACCURACY = 0, CPS_CLUSTER_RADIUS = 1, CPS_ROOT_APPROXIMATION = 2,
    CPS_ROOT_ISOLATION = 3, CPS_GOAL_DEFAULT = 4
} cps_goal;
typedef enum cps_method {
    CPS_METHOD_AUTO = 0, CPS_METHOD_CLASSICAL = 1,
    CPS_METHOD_SECULAR = 2, CPS_METHOD_PIECEWISE = 3
} cps_method;

/* Exact endpoint = mantissa * 2^exponent; both are decimal integer strings.
   All returned string pointers are borrowed until the result is destroyed. */
typedef struct cps_dyadic_view { const char *mantissa; const char *exponent; } cps_dyadic_view;
/* Exact operations on a root end m 2^e (decimal mantissa and exponent, as results give them), done inside the
   library, so a caller needs no big-number library (the C++ layer is built on them). A malformed or NULL end is
   CPS_INVALID_ARGUMENT. */
typedef enum cps_rounding { CPS_ROUND_DOWN = 0, CPS_ROUND_NEAREST = 1, CPS_ROUND_UP = 2 } cps_rounding;
/* -1, 0 or 1 in *out as a < b, a = b, a > b. */
CPS_API cps_error cps_dyadic_compare(cps_dyadic_view a, cps_dyadic_view b, int *out);
/* The binary64 number in *out, rounded toward -infinity, to nearest or toward +infinity, correctly also where the
   exponent range overflows or underflows. */
CPS_API cps_error cps_dyadic_to_double(cps_dyadic_view x, cps_rounding rounding, double *out);
/* Text into buffer, at most size bytes with the terminating NUL; *length (if not NULL) receives the full length
   without the NUL, so size 0 asks for it. cps_dyadic_format: decimal as C's "%.*g" writes a double with `digits`
   (1 or more) as its precision, and as C++'s std::format("{:.N}") does: at most `digits` significant digits,
   correctly rounded to nearest with ties to even, trailing zeros dropped, the exponent form ("7e-30") below 10^-4
   and from 10^digits, '.' as the decimal point whatever the host's locale; a value beyond MPFR's exponent range
   (about 2^(+/-2^30)) is CPS_OUT_OF_RANGE.
   cps_dyadic_rational: the exact value, "p/q" or "p" in lowest terms (as long as the exponent makes it).
   cps_dyadic_midpoint and cps_dyadic_width: (a + b)/2 and b - a, exact, as "m@e" (the coefficient grammar's exact
   dyadic). */
CPS_API cps_error cps_dyadic_format(cps_dyadic_view x, unsigned digits, char *buffer, size_t size, size_t *length);
CPS_API cps_error cps_dyadic_rational(cps_dyadic_view x, char *buffer, size_t size, size_t *length);
CPS_API cps_error cps_dyadic_midpoint(cps_dyadic_view a, cps_dyadic_view b, char *buffer, size_t size, size_t *length);
CPS_API cps_error cps_dyadic_width(cps_dyadic_view a, cps_dyadic_view b, char *buffer, size_t size, size_t *length);
typedef struct cps_root_view {
    cps_dyadic_view real_lower, real_upper, imag_lower, imag_upper;
    uint64_t count;
    uint64_t multiplicity; /* 0 = not established; count is NOT multiplicity. */
    cps_realness realness;
    const char *certificate;
} cps_root_view;

/* CPS_STOP_ITERATION_BUDGET: the solve made the sweeps cps_options_set_max_iterations
   allows it; the termination is CPS_WORK_LIMIT. */
typedef enum cps_stop_reason {
    CPS_STOP_NONE = 0, CPS_STOP_WORK_BUDGET = 1, CPS_STOP_CELL_BUDGET = 2,
    CPS_STOP_PRECISION = 3, CPS_STOP_CANCELLED = 4, CPS_STOP_INPUT = 5,
    CPS_STOP_ITERATION_BUDGET = 6, CPS_STOP_TIME_LIMIT = 7
} cps_stop_reason;
/* The search's components when the call ended (cps_result_get_component), diagnostics: an index identifies an
   entry of this result only, not a region of the answer and not a persistent handle. */
typedef enum cps_component_state {
    CPS_COMPONENT_ACTIVE = 0, CPS_COMPONENT_FINISHED = 1,
    CPS_COMPONENT_PRECISION_LIMIT = 2, CPS_COMPONENT_CELL_LIMIT = 4
} cps_component_state;
typedef struct cps_component_progress {
    uint64_t factor_index, work_bits, live_cells, certified_count;
    cps_component_state state;
} cps_component_progress;

/* What a solve call spent, a promise. The caller sets size to
   sizeof(cps_statistics); a smaller size is CPS_INVALID_ARGUMENT and leaves *out unchanged; otherwise the library
   fills the fields that fit, leaves the bytes beyond them untouched and writes the bytes it filled to size. */
typedef struct cps_statistics {
    size_t size;
    uint64_t work_units;          /* controller work units the call spent, every engine and internal solve included */
    uint64_t sweeps;              /* sweeps of the simultaneous iterations: what the iteration budget counts */
    uint64_t precision_increases; /* raises of the outer working precision */
    uint64_t final_work_bits;     /* the last outer working precision (cps_result_work_bits) */
    uint64_t highest_work_bits;   /* the highest outer working precision reached */
    uint64_t regions;             /* exported regions (cps_result_size) */
    uint64_t roots;               /* the sum of the regions' counts */
    uint64_t resolved_regions;    /* regions whose multiplicity is established */
    uint64_t exact_factors;       /* factors produced by exact preparation */
    uint64_t total_ns;            /* wall time of the call, on the solve's monotonic clock */
    uint64_t preparation_ns;      /* wall time of exact preparation, every engine included */
} cps_statistics;
/* The engines' counters: diagnostics, whose fields may be added, removed or renamed in any release. The size
   protocol is cps_statistics'. Stage times overlap. */
typedef struct cps_engine_counters {
    size_t size;
    uint64_t input_enclosure_ns, isolation_ns, coverage_validation_ns, backend_ns, taylor_shift_ns, counting_ns,
        refinement_ns, search_ns, candidate_ns, evaluation_ns, interaction_ns, validation_ns, regeneration_ns;
    uint64_t precision_attempts, isolation_passes, disk_tests, graeffe_steps, taylor_shifts, taylor_cache_hits,
        boxes_split, boxes_excluded, newton_attempts, newton_successes, backend_calls, evaluations, evaluation_terms,
        root_interactions, validation_calls, secular_regenerations, double_seed_calls, candidate_stagnations,
        real_refinement_steps, real_bisections;
    /* the search's state at the end */
    uint64_t pending_components, live_cells, min_pending_bits, max_pending_bits, input_uncertainty;
} cps_engine_counters;

/* Advanced: a cancellation flag the library owns (no caller code runs during a solve). cps_cancel_token_request may be called from any thread at any time, also while solves read the token,
   and is lock-free; a request is sticky. Every solve whose options name the token stops at its next cooperative
   checkpoint with CPS_CANCELLED, stop reason CPS_STOP_CANCELLED and its certified partial answer; a token requested
   before the solve starts stops it at once. A solve whose goal is met reports CPS_SUCCESS, whatever is requested
   afterwards. Options that name a token share it, so the caller may destroy its
   handle at any time. NULL does nothing for request and destroy, and reads as not requested. */
typedef struct cps_cancel_token cps_cancel_token;

CPS_API uint32_t cps_abi_version(void);
CPS_API const char *cps_version(void);
CPS_API const char *cps_backend_version(void);
/* Thread-local diagnostic, valid until the next error-producing call on this thread. */
CPS_API const char *cps_last_error(void);
/* Release the calling thread's FLINT caches, before the thread exits (a short-lived
   thread that skips it leaks them). Objects made on the thread stay valid and may
   be used and destroyed on any thread. Coordinate with any other FLINT users on the
   thread. */
CPS_API void cps_cleanup_thread(void);

CPS_API cps_error cps_polynomial_create(cps_polynomial **out);
CPS_API void cps_polynomial_destroy(cps_polynomial *polynomial);
/* Ascending coefficient order. Grammar: integer, p/q, decimal/scientific,
   or integer@integer for an exact dyadic. No NaN/Inf, hex, or expressions.
   The real part must not be NULL; a NULL imaginary part or radius is 0.
   Setters are transactional: invalid inputs do not change the polynomial.
   The polynomial keeps its nonzero terms only, so any exponent up to 2^63 - 1 can be set and a polynomial of few
   terms stays small whatever its degree (an exact zero removes a term; a ball whose value is zero is a term). A
   solve reads every coefficient up to the degree: a degree beyond what memory holds fails there, with
   CPS_OUT_OF_MEMORY or CPS_OUT_OF_RANGE, not when the coefficient is set. */
CPS_API cps_error cps_polynomial_set_coefficient(cps_polynomial *, uint64_t index,
    const char *real, const char *imag);
CPS_API cps_error cps_polynomial_set_coefficient_d(cps_polynomial *, uint64_t index,
    double real, double imag);
CPS_API cps_error cps_polynomial_set_ball_coefficient(cps_polynomial *, uint64_t index,
    const char *real, const char *imag, const char *real_radius, const char *imag_radius);
/* A new polynomial with the same coefficients (balls included) in *out. */
CPS_API cps_error cps_polynomial_copy(const cps_polynomial *, cps_polynomial **out);
/* Reading a polynomial. cps_polynomial_degree: the highest power whose coefficient is not exactly zero (a ball
   counts as nonzero), -1 for the zero polynomial or NULL. cps_polynomial_get_coefficient: one part of the
   coefficient of x^index (zero beyond the degree) as exact text in the coefficient grammar, so that the text sets
   the same value again ("p/q", "p", or "m@e" for a binary exponent), into buffer as cps_dyadic_format writes.
   cps_polynomial_terms: the exponents of the terms (the powers whose coefficients are not exactly zero), ascending,
   at most capacity of them into exponents; *count (if not NULL) receives how many terms there are, so capacity 0
   asks for it. */
typedef enum cps_coefficient_part {
    CPS_PART_REAL = 0, CPS_PART_IMAG = 1, CPS_PART_REAL_RADIUS = 2, CPS_PART_IMAG_RADIUS = 3
} cps_coefficient_part;
CPS_API int64_t cps_polynomial_degree(const cps_polynomial *);
CPS_API cps_error cps_polynomial_get_coefficient(const cps_polynomial *, uint64_t index, cps_coefficient_part part,
    char *buffer, size_t size, size_t *length);
CPS_API cps_error cps_polynomial_terms(const cps_polynomial *, uint64_t *exponents, size_t capacity, size_t *count);
/* Exact polynomial arithmetic, for building inputs: products, powers, recurrences, a polynomial from its roots.
   Each call makes a new polynomial in *out (NULL on error; destroy it with cps_polynomial_destroy); the operands
   stay unchanged and may be the same object. The arithmetic is exact, over the Gaussian rationals: every
   coefficient of an operand must be exact (a ball coefficient is CPS_INVALID_ARGUMENT), and a coefficient's
   exponent is expanded as written ("1@1000000" is a million-bit integer). The arithmetic works on the terms: a
   product (or a step of a power) pairs terms when that pairs no more of them than the dense result has
   coefficients, and is dense otherwise, so (x^1000000 + 1)^2 stays three terms. A result whose degree exceeds
   2^63 - 1 is CPS_OUT_OF_RANGE. Expanding a product or a recurrence gives up the structure the program route uses
   (a product's factors are solved one by one there); a program keeps it. */
CPS_API cps_error cps_polynomial_add(const cps_polynomial *a, const cps_polynomial *b, cps_polynomial **out);
CPS_API cps_error cps_polynomial_subtract(const cps_polynomial *a, const cps_polynomial *b, cps_polynomial **out);
CPS_API cps_error cps_polynomial_multiply(const cps_polynomial *a, const cps_polynomial *b, cps_polynomial **out);
/* a times the number real + i imag (the coefficient grammar above; a NULL imag means 0). */
CPS_API cps_error cps_polynomial_scale(const cps_polynomial *a, const char *real, const char *imag, cps_polynomial **out);
/* a to the power exponent; any polynomial to the power 0 is 1. */
CPS_API cps_error cps_polynomial_power(const cps_polynomial *a, uint64_t exponent, cps_polynomial **out);
/* The derivative of a, exactly. */
CPS_API cps_error cps_polynomial_derivative(const cps_polynomial *a, cps_polynomial **out);
/* The product of (x - r_k)^m_k over count roots r_k = real[k] + i imag[k] (the coefficient grammar; imag NULL means
   every root real, a NULL imag[k] that root), with m_k = multiplicity[k] (NULL means 1 each; 0 leaves the root out).
   count 0 gives the polynomial 1. */
CPS_API cps_error cps_polynomial_from_roots(const char *const *real, const char *const *imag,
    const uint64_t *multiplicity, size_t count, cps_polynomial **out);

/* Owned immutable input, distinct poles, nonzero residues; term storage and
   arithmetic are checked against representable sizes and the resource limits
   (NULL: the defaults of cps_resource_limits_default). Empty sum is the nonzero
   constant -1. Secular input currently requires exact normalization;
   refused expansion returns CPS_OUT_OF_RANGE with an explanatory diagnostic. */
CPS_API cps_error cps_secular_create(const cps_secular_term_view *, size_t count, const cps_resource_limits *limits,
    cps_secular **out);
CPS_API void cps_secular_destroy(cps_secular *);
CPS_API cps_error cps_solve_secular(const cps_secular *, const cps_options *, cps_result **out);

/* Program input: a polynomial given by a
   straight-line program of exact constants, the variable, +, -, x, negation
   and integer powers. Each builder call appends one node and writes its index
   to *out (CPS_NO_NODE when the call fails); nodes are numbered from 0 in order, and an operand must be an
   earlier node. Constants use the coefficient literal grammar (imag may be
   null for 0). The text form ("PSPROGRAM 1", one operation a line) is parsed
   by cps_program_parse. A failed call leaves the program unchanged.

   Solving a program (cps_solve_program): the program is never expanded; its
   degree and leading coefficient are proved from its top coefficients, its roots are bounded by a
   majorant, and the structured route iterates on evaluations of the program
   (Aberth mode for Auto and classical, secular mode for secular; piecewise
   is a usage error, CPS_INVALID_ARGUMENT). All complex roots: a real-only query,
   the isolation goal and a region are usage errors. A region off the real axis is reported non-real.
   When every constant of the program is real, a region holding one distinct root whose reflection in
   the real axis meets no other region is reported real (by conjugation; its imaginary interval is 0);
   any other region meeting the axis is left unknown.
   Goals: CPS_GOAL_DEFAULT is CPS_ROOT_APPROXIMATION (counted clusters) for a
   program; CPS_CLUSTER_RADIUS with any of the three methods (the one route
   serves them all, so, unlike coefficient input, explicit classical and
   secular do not refuse it); CPS_ROOT_ACCURACY is a best
   effort: a cluster that does not separate within the precision cap ends
   the call with CPS_PRECISION_LIMIT and the counted cluster, and the
   diagnostics say whether the roots are known distinct ("distinct roots,
   unresolved at the cap", when the image of the program modulo a word prime
   proved the polynomial square-free) or "possibly a multiple root". A
   multiplicity above one is certified only through a power in the program
   (the output c q^k, or a factor q^k of a product): q's regions, each
   count and each multiplicity q's region proves, times k; or when every
   root is zero (the program computes c x^n, proved from its majorant and
   its constants' denominators): the point 0, count and multiplicity n.
   Refusals:
   CPS_ZERO_POLYNOMIAL for a program that computes zero,
   CPS_INDETERMINATE_DEGREE when its leading terms cancel beyond what the
   library resolves, CPS_UNSUPPORTED_QUERY with CPS_STOP_INPUT when its length
   or storage exceeds the resource limits; a nonzero constant is answered as
   degree 0 with no roots. Budgets, cancellation and the iteration budget are
   as for coefficients. */
CPS_API cps_error cps_program_create(cps_program **out);
CPS_API void cps_program_destroy(cps_program *);
CPS_API cps_error cps_program_variable(cps_program *, uint64_t *out);
CPS_API cps_error cps_program_constant(cps_program *, const char *real, const char *imag, uint64_t *out);
CPS_API cps_error cps_program_add(cps_program *, uint64_t a, uint64_t b, uint64_t *out);
CPS_API cps_error cps_program_subtract(cps_program *, uint64_t a, uint64_t b, uint64_t *out);
CPS_API cps_error cps_program_multiply(cps_program *, uint64_t a, uint64_t b, uint64_t *out);
CPS_API cps_error cps_program_negate(cps_program *, uint64_t a, uint64_t *out);
CPS_API cps_error cps_program_power(cps_program *, uint64_t a, uint64_t exponent, uint64_t *out);
CPS_API cps_error cps_program_set_output(cps_program *, uint64_t node);
CPS_API cps_error cps_program_parse(const char *text, cps_program **out);
/* Not a node: what a failed builder call writes to *out, and cps_program_output without an output. */
#define CPS_NO_NODE UINT64_MAX
/* The number of nodes, and the output node (CPS_NO_NODE if none is set). */
CPS_API uint64_t cps_program_size(const cps_program *);
CPS_API uint64_t cps_program_output(const cps_program *);
/* The program must have an output (CPS_INVALID_ARGUMENT otherwise); it is
   copied, so it may be destroyed after the call. Null options selects the
   defaults. */
CPS_API cps_error cps_solve_program(const cps_program *, const cps_options *, cps_result **out);

CPS_API cps_error cps_options_create(cps_options **out);
CPS_API void cps_options_destroy(cps_options *options);
/* Auto preserves the public selection policy. Explicit methods support exact
   coefficients, all-complex ROOT_ACCURACY/ROOT_APPROXIMATION queries. Native
   secular input supports Auto/Secular. Program input supports Auto, classical
   and secular (never real-only; see the program note above). Other explicit
   combinations are usage errors: the solve returns CPS_INVALID_ARGUMENT and no result;
   a forced method on ball input is answered with
   CPS_UNSUPPORTED_QUERY (a refusal that depends on the input). Shared exact preparation/direct formulas and each
   method's internal refinement stages remain permitted. Explicit piecewise
   hands the stops that are not the caller's (its precision budget, a storage
   or coefficient-range refusal, its regime exit) to the maintained secular
   recipe on the original input, within the same call and work budget; the
   diagnostics name the handoff. Unknown integer values leave options
   unchanged. */
CPS_API cps_error cps_options_set_method(cps_options *, int method);
/* The default maximum outer working precision. A solve that needs more stops
   for precision with a complete cover, and for exact real input its diagnostics
   and cps_result_sufficient_work_bits name a cap that suffices (from Mahler's
   root-separation bound); a caller who needs more raises the cap. */
#define CPS_DEFAULT_MAX_WORK_BITS 65536
/* The smallest maximum working precision a setter accepts. */
#define CPS_MIN_MAX_WORK_BITS 32
/* The default relative accuracy of every root, in bits: binary64's 53-bit significand, so a root's region is about
   as narrow as a double can tell (its midpoint is not promised to be the root correctly rounded). */
#define CPS_DEFAULT_BITS 53
/* Defaults: CPS_DEFAULT_BITS requested bits, CPS_DEFAULT_MAX_WORK_BITS maximum outer working
   bits (cps_options_set_max_work_bits), 0 absolute tolerance (disabled), automatic iterations, unspecified
   resource allowances. The requested bits: 1 or more (0 is CPS_INVALID_ARGUMENT); a setter's range is limited by
   backend representability (CPS_OUT_OF_RANGE beyond it), not a fixed bit cap. */
CPS_API cps_error cps_options_set_precision(cps_options *, uint64_t bits);
CPS_API cps_error cps_options_set_absolute_tolerance(cps_options *, const char *tolerance);
/* CPS_ROOT_ACCURACY (CPS_GOAL_DEFAULT's goal for coefficient and secular input) requires disjoint individual-root regions with
   established multiplicities. CPS_ROOT_APPROXIMATION also accepts disjoint
   counted clusters when every root in the exported rectangle is approximated
   by its midpoint to the requested true-relative accuracy (or explicit absolute
   tolerance). Count does not assert multiplicity. Real-only queries currently
   refine individual roots for either accuracy goal.
   CPS_ROOT_ISOLATION, under CPS_METHOD_AUTO on exact coefficient input, asks
   for every distinct complex root alone in a disjoint closed rectangle with
   its multiplicity established. A real-only query asks for every real root
   alone in a closed real interval and requires exact real coefficients.
   This goal promises no width: requested bits are not used. Ball input is
   refused with CPS_UNSUPPORTED_QUERY. A forced method, native secular input
   or program input is a usage error (CPS_INVALID_ARGUMENT from the solve).
   Refused exact factor preparation reports CPS_UNSUPPORTED_QUERY with cps_result_preparation_limited true.
   Approximations are refused for isolation results; the exact bounds are the answer.
   Every goal but CPS_CLUSTER_RADIUS clears cluster_radius (CPS_GOAL_DEFAULT
   included: it restores the initial goal). Selecting
   CPS_CLUSTER_RADIUS requires a positive radius previously set by
   cps_options_set_cluster_radius. Pass a cps_goal constant; other integer
   values are rejected without changing the options. The integer parameter
   allows validation before constructing a C++ enum value. */
CPS_API cps_error cps_options_set_goal(cps_options *, int goal);
/* A positive radius selects certified cluster coverage with each exported
   rectangle at Euclidean radius <= this absolute bound from its midpoint.
   Separation/multiplicity are NOT required. Zero clears the radius and
   restores CPS_GOAL_DEFAULT, including after CPS_ROOT_APPROXIMATION was selected.
   A real-only query with this goal, or a forced method with it on coefficient
   input, is a usage error (CPS_INVALID_ARGUMENT from the solve). */
CPS_API cps_error cps_options_set_cluster_radius(cps_options *, const char *radius);
/* The iteration budget of a solve: at most this many sweeps of the
   simultaneous iterations (the Aberth and secular packets of the classical and
   secular methods, at every precision, the binary64 ones included), counted by
   cps_statistics.sweeps. Zero, the default, sets no limit. When the
   budget is spent the solve returns its last certified cover with termination
   CPS_WORK_LIMIT and stop reason CPS_STOP_ITERATION_BUDGET. The value changes no
   engine's own packet lengths.
   Not under the budget: Subdivision (the complete method a solve hands a
   factor to) and the real engine's exact isolation, which make no sweeps (a
   real-only query on exact real input is answered by the complex methods,
   whose sweeps are under the budget, or by the piecewise method's axis
   charts, which are not), and the piecewise
   method's own work (its solves of small local models do appear in
   cps_statistics.sweeps, but they are not sweeps of the input). Their limits are
   the work and cell budgets, the time limit and the cancellation token. A solve the
   piecewise method hands to the secular method is under the budget from there
   on. This is not a time bound: use the time limit for a deadline. */
CPS_API cps_error cps_options_set_max_iterations(cps_options *, uint64_t sweeps);
/* The search budget's defaults: the work units and live cells a solve may use when the caller sets none
   (cps_options_set_work_budget and cps_options_set_live_cells below). */
#define CPS_DEFAULT_WORK_UNITS UINT64_C(100000)
/* The automatic work budget grows with the input: at least this many units a root (the degree). */
#define CPS_DEFAULT_WORK_UNITS_PER_ROOT UINT64_C(64)
#define CPS_DEFAULT_LIVE_CELLS UINT64_C(65536)
CPS_API cps_error cps_options_set_real_only(cps_options *, int enabled);
/* A closed rectangle for the roots asked for: ends in the coefficient grammar,
   NULL or "" for an unbounded side. A real-only query answers exactly its
   real roots x with real_lower <= x <= real_upper when imag_lower <= 0 <=
   imag_upper, and none when the rectangle misses the real axis (an empty,
   complete answer); complete and goal_met refer to those roots. Ends are
   closed: a root exactly at an end is decided exactly and is in the answer
   (an exact point when the end is dyadic); a root whose interval holds an
   end is placed on its side by the input's signs and its interval cut to
   that side, except a root of even multiplicity, which is refined until it
   lies on one side and may refine every root.
   Only real-only queries under CPS_METHOD_AUTO on exact real polynomial input
   take a region (CPS_UNSUPPORTED_QUERY on inexact or nonreal input; any other
   query with a region is a usage error, CPS_INVALID_ARGUMENT from the solve).
   The setter checks each end's literal only; four unbounded sides clear the
   region. The solve checks the rest before any work, each failure a usage error:
   on each axis the lower end must not exceed the upper (decided exactly within
   the precision cap), and each real end must expand to an exact rational within
   the exact storage policy (the resource limits: exact_materialization_bits, and
   CPS_RESOURCE_ALLOW_EXPONENT_EXPANSION for a compact exponent).
   cps_options_get_region writes NULL for an unbounded side. */
CPS_API cps_error cps_options_set_region(cps_options *, const char *real_lower, const char *real_upper,
                                      const char *imag_lower, const char *imag_upper);
CPS_API cps_error cps_options_get_region(const cps_options *, const char **real_lower, const char **real_upper,
                                      const char **imag_lower, const char **imag_upper);
CPS_API cps_error cps_cancel_token_create(cps_cancel_token **out);
CPS_API void cps_cancel_token_request(cps_cancel_token *token);
CPS_API int cps_cancel_token_requested(const cps_cancel_token *token);
CPS_API void cps_cancel_token_destroy(cps_cancel_token *token);
/* Names the token that solves with these options read; NULL clears it. */
CPS_API cps_error cps_options_set_cancel_token(cps_options *, cps_cancel_token *token);
/* A wall-clock limit for one solve, in seconds, on a monotonic clock from the start of the solve call, covering
   everything the call does (the real-only route's solves and a method's hand-off included). 0 (the default) or
   +infinity means none; NaN or a negative value is CPS_INVALID_ARGUMENT. A solve that reaches it stops at its next
   cooperative checkpoint with CPS_CANCELLED, stop reason CPS_STOP_TIME_LIMIT and its certified partial answer; the
   overrun is bounded by the longest step that cannot be interrupted. A solve whose goal is met before the limit is
   seen reports CPS_SUCCESS. A time-limited answer depends on the machine; the work and iteration budgets are
   deterministic. A finite limit beyond the clock's range uses its latest representable deadline. */
CPS_API cps_error cps_options_set_time_limit(cps_options *, double seconds);
/* The threads a solve may use. 1, the default, starts none. A larger count lets each solve call start up to workers
   - 1 threads, at most the machine's hardware threads, and join them before it returns; where they cannot be started
   the solve runs on its own thread. The answer, its counters and its diagnostics are the same at every count (a time
   limit or a token stops where the machine's speed puts it, as with one). The workers never call caller code, never
   use FLINT's thread pool and start with every signal blocked. 0 is CPS_INVALID_ARGUMENT. */
CPS_API cps_error cps_options_set_workers(cps_options *, uint32_t workers);
/* cps_options_copy: an independent copy of every setting, naming the same token as the source.
   cps_options_set_max_work_bits: the maximum outer working precision, CPS_MIN_MAX_WORK_BITS or more (it may be below
   the requested bits).
   Advanced, the search budget: most callers keep the defaults, CPS_DEFAULT_WORK_UNITS controller work units a solve or
   CPS_DEFAULT_WORK_UNITS_PER_ROOT a root of the input if that is more (a degree-8,192 input gets 524,288), and
   CPS_DEFAULT_LIVE_CELLS live search squares,
   and use the time limit for a deadline. cps_options_set_work_budget: 0 selects that automatic budget, any other
   value is used as given; cps_options_set_live_cells: 0 selects CPS_DEFAULT_LIVE_CELLS. A unit is a
   degree-proof/preprocessing pass, backend call, bounded local certification/refinement attempt, disk exclusion
   attempt, square split, or frontier transition. The work budget is a soft, deterministic limit: it is checked
   between units, the unit in progress is finished (at high precision one unit can take long), and nothing is
   promised about how far a solve gets within it. A solve that stops for it returns the certified cover it has,
   with CPS_STOP_WORK_BUDGET. These are not wall-time, byte-allocation, or backend-internal-operation guarantees.
   cps_options_set_resource_limits / _get_: the resource limits above (NULL is CPS_INVALID_ARGUMENT, so is an unknown
   flag); they set the exact storage policy for exact preparation and for a region's ends.
   A getter returns the setting (NULL options: its default); a string getter returns the literal as last set,
   trimmed ("0" when none was set; NULL for an unbounded region side), valid until the options change or go. */
CPS_API cps_error cps_options_copy(const cps_options *source, cps_options **out);
CPS_API cps_error cps_options_set_max_work_bits(cps_options *, uint64_t bits);
CPS_API cps_error cps_options_set_work_budget(cps_options *, uint64_t units);
CPS_API cps_error cps_options_set_live_cells(cps_options *, uint64_t cells);
CPS_API cps_error cps_options_set_resource_limits(cps_options *, const cps_resource_limits *);
CPS_API cps_error cps_options_get_resource_limits(const cps_options *, cps_resource_limits *out);
CPS_API cps_method cps_options_method(const cps_options *);
CPS_API uint64_t cps_options_precision(const cps_options *);
CPS_API uint64_t cps_options_max_work_bits(const cps_options *);
CPS_API cps_goal cps_options_goal(const cps_options *);
CPS_API const char *cps_options_cluster_radius(const cps_options *);
CPS_API const char *cps_options_absolute_tolerance(const cps_options *);
CPS_API int cps_options_real_only(const cps_options *);
CPS_API uint64_t cps_options_work_budget(const cps_options *);
CPS_API uint64_t cps_options_live_cells(const cps_options *);
CPS_API uint64_t cps_options_max_iterations(const cps_options *);
CPS_API double cps_options_time_limit(const cps_options *);
CPS_API uint32_t cps_options_workers(const cps_options *);

/* Null options selects defaults. Mathematical non-success returns CPS_OK with
   a result describing the termination; errors denote invalid API use/runtime errors.
   On error *out is null. Polynomial and options are read-only during solving. */
CPS_API cps_error cps_solve(const cps_polynomial *, const cps_options *, cps_result **out);
CPS_API void cps_result_destroy(cps_result *result);
/* -1 for zero; otherwise the degree upper bound if leading nonzero is unresolved. */
CPS_API int64_t cps_result_degree(const cps_result *);
CPS_API cps_termination cps_result_termination(const cps_result *);
CPS_API int cps_result_complete(const cps_result *);
CPS_API int cps_result_goal_met(const cps_result *);
CPS_API cps_goal cps_result_goal(const cps_result *);
CPS_API uint64_t cps_result_work_bits(const cps_result *); /* Last outer precision, NOT backend peak. */
CPS_API size_t cps_result_size(const cps_result *);
CPS_API const char *cps_result_diagnostics(const cps_result *);
/* Static English names: "success", "precision limit", "cancelled", ...; "none", "work budget", "time limit", ...
   An unknown value gives "unknown". */
CPS_API const char *cps_termination_name(cps_termination);
CPS_API const char *cps_stop_reason_name(cps_stop_reason);
CPS_API size_t cps_result_component_count(const cps_result *);
/* cps_result_get_statistics and cps_result_get_engine_counters follow their structs' size protocol;
   NULL result or output is CPS_INVALID_ARGUMENT. cps_result_get_component: CPS_OUT_OF_RANGE for an index at or past
   cps_result_component_count. cps_result_stop_reason: why the solve stopped (NULL: CPS_STOP_INPUT).
   cps_result_preparation_limited: 1 when exact preparation was refused or limited by the resource limits, so that
   multiplicities may be unknown (NULL: 0). */
CPS_API cps_error cps_result_get_statistics(const cps_result *, cps_statistics *out);
CPS_API cps_error cps_result_get_engine_counters(const cps_result *, cps_engine_counters *out);
CPS_API cps_error cps_result_get_component(const cps_result *, size_t index, cps_component_progress *out);
CPS_API cps_stop_reason cps_result_stop_reason(const cps_result *);
CPS_API int cps_result_preparation_limited(const cps_result *);
/* Facts about the solve (the diagnostics may repeat them as text; these getters are the reference).
   cps_result_sufficient_work_bits: 0 unless the solve stopped at its precision cap (CPS_PRECISION_LIMIT with
   CPS_STOP_PRECISION) on exact real input, whose distinct factors are then integer polynomials; then
   2 (separation bits + size bits + requested bits) from Mahler's root-separation bound, a cap at which the roots
   can be isolated at the requested precision (a bound, not a promise that a solve under it meets its goal within
   its other budgets). NULL: 0.
   cps_result_squarefree: PROVED, every root is simple; DISPROVED, a root of multiplicity above one is proved;
   UNKNOWN, neither (also where exact preparation was limited, or on ball input). A program's verdict comes from
   its image modulo a word-size prime, so it is PROVED or UNKNOWN. NULL: UNKNOWN.
   cps_result_handoff: the piecewise method, forced or chosen by Auto, handed the solve to the secular method within
   the call, and why: it spent more than a global solve's predicted cost (REGIME_EXIT), its storage was refused
   (STORAGE), the coefficients' exponents exceed its working range (RANGE), or it reached its precision budget
   (PRECISION); the answer then comes from the secular method. NULL: NONE. */
typedef enum cps_squarefree { CPS_SQUAREFREE_UNKNOWN = 0, CPS_SQUAREFREE_PROVED = 1, CPS_SQUAREFREE_DISPROVED = 2 } cps_squarefree;
typedef enum cps_handoff {
    CPS_HANDOFF_NONE = 0, CPS_HANDOFF_REGIME_EXIT = 1, CPS_HANDOFF_STORAGE = 2, CPS_HANDOFF_RANGE = 3,
    CPS_HANDOFF_PRECISION = 4
} cps_handoff;
CPS_API uint64_t cps_result_sufficient_work_bits(const cps_result *);
CPS_API cps_squarefree cps_result_squarefree(const cps_result *);
CPS_API cps_handoff cps_result_handoff(const cps_result *);
/* The regions come in increasing order of their exact midpoint's real part, then its imaginary part, on every
   route. */
CPS_API cps_error cps_result_get_root(const cps_result *, size_t index, cps_root_view *out);
/* Explicit lossy midpoint conversion; refuses unmet goals, unresolved clusters,
   isolation results (no width was asked), non-finite results and rounding a
   nonzero midpoint component to zero.
   Nonzero subnormals are allowed and may lose substantial relative accuracy;
   the requested precision applies to the exact enclosure, not this conversion. */
CPS_API cps_error cps_result_get_approximation(const cps_result *, size_t index, double *real, double *imag);

#ifdef __cplusplus
}
#endif
#endif
