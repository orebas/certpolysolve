/* The C API's contract where no other test reaches it (docs/api-spec-0.5.md, change C35): the version macros and
   calls, cps_last_error's rule (a failed call describes itself; a call that returns cps_error and succeeds leaves the
   message empty; a call of another return type leaves it alone), every NULL fallback of section 2.3, the outputs a
   failed call leaves (a handle NULL, a node CPS_NO_NODE), and the component accessors. C11, the public header only. */
#include <certpolysolve/certpolysolve.h>
#include <stdio.h>
#include <string.h>

#define CHECK(test)                                                                                          \
    do {                                                                                                     \
        if (!(test)) {                                                                                       \
            fprintf(stderr, "FAIL line %d: %s (last error: \"%s\")\n", __LINE__, #test, cps_last_error());   \
            return 1;                                                                                        \
        }                                                                                                    \
    } while (0)
/* A failed call's message is not empty; the next successful cps_error call empties it. */
#define FAILS(call) CHECK((call) != CPS_OK && cps_last_error()[0] != '\0')
#define CLEARS(call) CHECK((call) == CPS_OK && cps_last_error()[0] == '\0')

static int versions(void) {
    char expected[32];
    snprintf(expected, sizeof expected, "%d.%d.%d", CPS_VERSION_MAJOR, CPS_VERSION_MINOR, CPS_VERSION_PATCH);
    CHECK(strcmp(cps_version(), CPS_VERSION_STRING) == 0);
    CHECK(strncmp(CPS_VERSION_STRING, expected, strlen(expected)) == 0); /* "0.5.0" or "0.5.0-dev" */
    CHECK(cps_abi_version() == CPS_ABI_VERSION && CPS_ABI_VERSION == 2u);
    CHECK(cps_backend_version() != NULL && cps_backend_version()[0] != '\0');
    return 0;
}

static int last_error(void) {
    cps_options *o = NULL;
    cps_polynomial *p = NULL;
    cps_result *r = NULL;
    CLEARS(cps_options_create(&o));
    CLEARS(cps_polynomial_create(&p));
    /* Each failure sets a message, and the next success of any kind of cps_error call clears it. */
    FAILS(cps_options_set_goal(o, 99));
    CLEARS(cps_options_set_goal(o, CPS_ROOT_APPROXIMATION));
    FAILS(cps_polynomial_set_coefficient(p, 0, "not a number", NULL));
    CLEARS(cps_polynomial_set_coefficient(p, 0, "-2", NULL));
    CLEARS(cps_polynomial_set_coefficient(p, 2, "1", NULL));
    FAILS(cps_options_set_precision(o, 0));
    CLEARS(cps_options_set_goal(o, CPS_ROOT_ACCURACY)); /* every root resolved, so its approximation is given */
    CLEARS(cps_solve(p, o, &r));
    cps_root_view view;
    FAILS(cps_result_get_root(r, 99, &view));
    CLEARS(cps_result_get_root(r, 0, &view));
    double re = 0, im = 0;
    FAILS(cps_result_get_approximation(r, 99, &re, &im));
    CLEARS(cps_result_get_approximation(r, 0, &re, &im));
    int order = 0;
    const cps_dyadic_view a = {"1", "0"}, b = {"3", "-1"};
    FAILS(cps_dyadic_compare(a, (cps_dyadic_view){"x", "0"}, &order));
    CLEARS(cps_dyadic_compare(a, b, &order));
    CHECK(order == -1);
    /* A call of another return type leaves the message as it was. */
    FAILS(cps_options_set_max_work_bits(o, 1));
    const char *kept = cps_last_error();
    char copy[512];
    snprintf(copy, sizeof copy, "%s", kept);
    (void)cps_result_size(r);
    (void)cps_options_precision(o);
    (void)cps_result_termination(NULL);
    CHECK(strcmp(cps_last_error(), copy) == 0);
    cps_result_destroy(r);
    cps_polynomial_destroy(p);
    cps_options_destroy(o);
    return 0;
}

static int null_fallbacks(void) {
    /* Getters of a NULL result (section 2.3). */
    CHECK(cps_result_degree(NULL) == -1);
    CHECK(cps_result_termination(NULL) == CPS_UNSUPPORTED_QUERY);
    CHECK(cps_result_stop_reason(NULL) == CPS_STOP_INPUT);
    CHECK(cps_result_complete(NULL) == 0 && cps_result_goal_met(NULL) == 0 && cps_result_preparation_limited(NULL) == 0);
    CHECK(cps_result_goal(NULL) == CPS_ROOT_ACCURACY);
    CHECK(cps_result_work_bits(NULL) == 0 && cps_result_size(NULL) == 0 && cps_result_component_count(NULL) == 0);
    CHECK(cps_result_sufficient_work_bits(NULL) == 0);
    CHECK(strcmp(cps_result_diagnostics(NULL), "null result") == 0);
    CHECK(cps_result_squarefree(NULL) == CPS_SQUAREFREE_UNKNOWN);
    CHECK(cps_result_handoff(NULL) == CPS_HANDOFF_NONE);
    /* Programs and tokens. */
    CHECK(cps_program_size(NULL) == 0 && cps_program_output(NULL) == CPS_NO_NODE);
    CHECK(cps_cancel_token_requested(NULL) == 0);
    cps_cancel_token_request(NULL);
    cps_resource_limits_default(NULL);
    /* Options getters: the setting's default. */
    CHECK(cps_options_method(NULL) == CPS_METHOD_AUTO && cps_options_precision(NULL) == CPS_DEFAULT_BITS);
    CHECK(cps_options_max_work_bits(NULL) == CPS_DEFAULT_MAX_WORK_BITS && cps_options_goal(NULL) == CPS_GOAL_DEFAULT);
    CHECK(strcmp(cps_options_cluster_radius(NULL), "0") == 0 && strcmp(cps_options_absolute_tolerance(NULL), "0") == 0);
    CHECK(cps_options_real_only(NULL) == 0 && cps_options_work_budget(NULL) == 0);
    CHECK(cps_options_live_cells(NULL) == CPS_DEFAULT_LIVE_CELLS && cps_options_max_iterations(NULL) == 0);
    CHECK(cps_options_time_limit(NULL) == 0.0);
    /* Destroying NULL does nothing. */
    cps_polynomial_destroy(NULL);
    cps_options_destroy(NULL);
    cps_result_destroy(NULL);
    cps_secular_destroy(NULL);
    cps_program_destroy(NULL);
    cps_cancel_token_destroy(NULL);
    /* A NULL handle or output of a cps_error function is CPS_INVALID_ARGUMENT. */
    CHECK(cps_options_set_method(NULL, CPS_METHOD_AUTO) == CPS_INVALID_ARGUMENT);
    CHECK(cps_polynomial_create(NULL) == CPS_INVALID_ARGUMENT && cps_options_create(NULL) == CPS_INVALID_ARGUMENT);
    CHECK(cps_program_create(NULL) == CPS_INVALID_ARGUMENT && cps_cancel_token_create(NULL) == CPS_INVALID_ARGUMENT);
    cps_statistics statistics;
    statistics.size = sizeof statistics;
    CHECK(cps_result_get_statistics(NULL, &statistics) == CPS_INVALID_ARGUMENT);
    /* NULL options of a solve and NULL limits of a secular equation are the defaults. */
    cps_polynomial *p = NULL;
    cps_result *r = NULL;
    CHECK(cps_polynomial_create(&p) == CPS_OK && cps_polynomial_set_coefficient(p, 1, "1", NULL) == CPS_OK &&
          cps_polynomial_set_coefficient(p, 0, "-1/3", NULL) == CPS_OK);
    CHECK(cps_solve(p, NULL, &r) == CPS_OK && cps_result_goal_met(r) && cps_result_goal(r) == CPS_ROOT_ACCURACY);
    cps_result_destroy(r);
    cps_polynomial_destroy(p);
    const cps_secular_term_view term = {"1", NULL, "2", NULL};
    cps_secular *s = NULL;
    CHECK(cps_secular_create(&term, 1, NULL, &s) == CPS_OK && s);
    cps_secular_destroy(s);
    return 0;
}

static int failed_outputs(void) {
    /* A failed call sets its handle output to NULL and its node output to CPS_NO_NODE. */
    cps_polynomial *p = (cps_polynomial *)&p;
    CHECK(cps_polynomial_copy(NULL, &p) == CPS_INVALID_ARGUMENT && p == NULL);
    cps_polynomial *sum = (cps_polynomial *)&sum;
    CHECK(cps_polynomial_add(NULL, NULL, &sum) == CPS_INVALID_ARGUMENT && sum == NULL);
    cps_options *copy = (cps_options *)&copy;
    CHECK(cps_options_copy(NULL, &copy) == CPS_INVALID_ARGUMENT && copy == NULL);
    cps_result *r = (cps_result *)&r;
    CHECK(cps_solve(NULL, NULL, &r) == CPS_INVALID_ARGUMENT && r == NULL);
    cps_program *program = (cps_program *)&program;
    CHECK(cps_program_parse("not a program", &program) == CPS_INVALID_ARGUMENT && program == NULL);
    cps_secular *s = (cps_secular *)&s;
    const cps_secular_term_view zero = {"0", NULL, "1", NULL}; /* a zero residue */
    CHECK(cps_secular_create(&zero, 1, NULL, &s) == CPS_INVALID_ARGUMENT && s == NULL);
    CHECK(cps_program_create(&program) == CPS_OK);
    uint64_t node = 7;
    CHECK(cps_program_negate(program, 5, &node) == CPS_INVALID_ARGUMENT && node == CPS_NO_NODE);
    CHECK(cps_program_size(program) == 0);
    cps_program_destroy(program);
    return 0;
}

static int components(void) {
    /* A solve whose search keeps components: forced classical on a polynomial with a double root. */
    cps_polynomial *p = NULL;
    cps_options *o = NULL;
    cps_result *r = NULL;
    CHECK(cps_polynomial_create(&p) == CPS_OK && cps_options_create(&o) == CPS_OK);
    CHECK(cps_polynomial_set_coefficient(p, 0, "1", NULL) == CPS_OK && cps_polynomial_set_coefficient(p, 1, "-2", NULL) == CPS_OK &&
          cps_polynomial_set_coefficient(p, 2, "1", NULL) == CPS_OK && cps_polynomial_set_coefficient(p, 3, "1", NULL) == CPS_OK);
    CHECK(cps_options_set_method(o, CPS_METHOD_CLASSICAL) == CPS_OK);
    CHECK(cps_solve(p, o, &r) == CPS_OK);
    const size_t n = cps_result_component_count(r);
    cps_component_progress c;
    for (size_t i = 0; i < n; ++i) {
        CHECK(cps_result_get_component(r, i, &c) == CPS_OK);
        CHECK(c.state == CPS_COMPONENT_ACTIVE || c.state == CPS_COMPONENT_FINISHED || c.state == CPS_COMPONENT_PRECISION_LIMIT ||
              c.state == CPS_COMPONENT_CELL_LIMIT);
    }
    CHECK(cps_result_get_component(r, n, &c) == CPS_OUT_OF_RANGE);
    CHECK(cps_result_get_component(r, 0, NULL) == CPS_INVALID_ARGUMENT);
    CHECK(cps_result_get_component(NULL, 0, &c) == CPS_INVALID_ARGUMENT);
    cps_result_destroy(r);
    cps_options_destroy(o);
    cps_polynomial_destroy(p);
    return 0;
}

/* An enum argument outside its enumeration is refused (CPS_INVALID_ARGUMENT), or named "unknown". Tested from C, where passing such a
   value is legal (a C enum is an integer type); in C++ the cast to the enum type is itself undefined for a value
   outside the enumeration's range. The library reads these arguments as integers before any use (c_api.cpp,
   enum_argument), so the refusal is well-defined for every C caller. */
static int enum_arguments(void) {
    cps_polynomial *p = NULL;
    char buffer[32];
    size_t length = 0;
    double d = 0;
    const cps_dyadic_view one = {"1", "0"};
    CHECK(cps_polynomial_create(&p) == CPS_OK && cps_polynomial_set_coefficient(p, 0, "1", NULL) == CPS_OK);
    CHECK(cps_polynomial_get_coefficient(p, 0, (cps_coefficient_part)7, buffer, sizeof buffer, &length) == CPS_INVALID_ARGUMENT);
    CHECK(cps_polynomial_get_coefficient(p, 0, (cps_coefficient_part)-1, buffer, sizeof buffer, &length) == CPS_INVALID_ARGUMENT);
    CHECK(cps_polynomial_get_coefficient(p, 0, CPS_PART_IMAG_RADIUS, buffer, sizeof buffer, &length) == CPS_OK);
    CHECK(cps_dyadic_to_double(one, (cps_rounding)9, &d) == CPS_INVALID_ARGUMENT);
    CHECK(cps_dyadic_to_double(one, (cps_rounding)-1, &d) == CPS_INVALID_ARGUMENT);
    CHECK(cps_dyadic_to_double(one, CPS_ROUND_UP, &d) == CPS_OK && d == 1.0);
    /* The name functions answer "unknown" for a value outside the enumeration, and the names for the enumerators. */
    CHECK(strcmp(cps_termination_name((cps_termination)99), "unknown") == 0 &&
          strcmp(cps_termination_name((cps_termination)-1), "unknown") == 0 &&
          strcmp(cps_termination_name(CPS_WORK_LIMIT), "work limit") == 0);
    CHECK(strcmp(cps_stop_reason_name((cps_stop_reason)99), "unknown") == 0 &&
          strcmp(cps_stop_reason_name(CPS_STOP_TIME_LIMIT), "time limit") == 0);
    cps_polynomial_destroy(p);
    return 0;
}

int main(void) {
    if (versions() || last_error() || null_fallbacks() || failed_outputs() || components() || enum_arguments()) return 1;
    cps_cleanup_thread();
    puts("C contract tests passed");
    return 0;
}
