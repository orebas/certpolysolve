#include <certpolysolve/certpolysolve.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#define CHECK(test)                                                                                          \
    do {                                                                                                     \
        if (!(test)) {                                                                                       \
            fprintf(stderr, "FAIL line %d: %s\n", __LINE__, #test);                                          \
            return 1;                                                                                        \
        }                                                                                                    \
    } while (0)

int main(void) {
    CHECK(cps_abi_version() == CPS_ABI_VERSION);
    CHECK(strcmp(cps_version(), CPS_VERSION_STRING) == 0);
    {
        /* The cancellation token and the time limit from C (the 0.5 API's stage 1). */
        cps_cancel_token *token = NULL;
        cps_options *timed = NULL;
        CHECK(cps_cancel_token_create(NULL) == CPS_INVALID_ARGUMENT);
        CHECK(cps_cancel_token_create(&token) == CPS_OK && !cps_cancel_token_requested(token));
        CHECK(cps_options_create(&timed) == CPS_OK);
        CHECK(cps_options_set_cancel_token(timed, token) == CPS_OK);
        cps_cancel_token_request(token);
        CHECK(cps_cancel_token_requested(token));
        cps_cancel_token_destroy(token); /* the options keep the request */
        CHECK(cps_options_set_cancel_token(NULL, NULL) == CPS_INVALID_ARGUMENT);
        CHECK(cps_options_set_time_limit(timed, -1) == CPS_INVALID_ARGUMENT);
        CHECK(cps_options_set_time_limit(timed, 2.5) == CPS_OK && strcmp(cps_last_error(), "") == 0);
        cps_polynomial *quadratic = NULL;
        cps_result *stopped = NULL;
        CHECK(cps_polynomial_create(&quadratic) == CPS_OK);
        CHECK(cps_polynomial_set_coefficient(quadratic, 0, "-2", "0") == CPS_OK &&
              cps_polynomial_set_coefficient(quadratic, 2, "1", "0") == CPS_OK);
        CHECK(cps_solve(quadratic, timed, &stopped) == CPS_OK);
        CHECK(cps_result_termination(stopped) == CPS_CANCELLED && !cps_result_goal_met(stopped));
        cps_cancel_token_request(NULL);
        cps_cancel_token_destroy(NULL);
        CHECK(!cps_cancel_token_requested(NULL));
        cps_result_destroy(stopped);
        cps_polynomial_destroy(quadratic);
        cps_options_destroy(timed);
    }
    CHECK(cps_polynomial_create(NULL) == CPS_INVALID_ARGUMENT);
    cps_polynomial *p = NULL;
    CHECK(cps_polynomial_create(&p) == CPS_OK);
    CHECK(cps_polynomial_set_coefficient(p, 0, "-2", "0") == CPS_OK);
    CHECK(cps_polynomial_set_coefficient(p, 2, "1", "0") == CPS_OK);
    CHECK(cps_polynomial_set_coefficient(p, 2, "1/0", "0") == CPS_INVALID_ARGUMENT);
    CHECK(cps_polynomial_set_coefficient_d(p, 2, NAN, 0) == CPS_INVALID_ARGUMENT);
    const char *invalid[] = {"", ".", "++1", "1e", "1e 2", "1/2.0", "0x3", "2**4", "Inf", "nan", "1@2@3"};
    for (size_t i = 0; i < sizeof(invalid) / sizeof(invalid[0]); ++i)
        CHECK(cps_polynomial_set_coefficient(p, 2, invalid[i], "0") == CPS_INVALID_ARGUMENT);
    CHECK(cps_polynomial_set_ball_coefficient(p, 0, "0", "0", "-1", "0") == CPS_INVALID_ARGUMENT);
    CHECK(cps_polynomial_set_coefficient(p, UINT64_MAX, "0", "0") == CPS_OUT_OF_RANGE);
    cps_result *r = NULL;
    CHECK(cps_solve(p, NULL, &r) == CPS_OK);
    CHECK(cps_result_complete(r) && cps_result_goal_met(r));
    CHECK(cps_result_goal(r) == CPS_ROOT_ACCURACY);
    CHECK(cps_result_size(r) == 2);
    cps_statistics statistics;
    statistics.size = sizeof statistics;
    cps_engine_counters counters;
    counters.size = sizeof counters;
    CHECK(cps_result_get_statistics(NULL, &statistics) == CPS_INVALID_ARGUMENT);
    CHECK(cps_result_get_statistics(r, NULL) == CPS_INVALID_ARGUMENT);
    CHECK(cps_result_get_statistics(r, &statistics) == CPS_OK);
    CHECK(cps_result_get_engine_counters(r, &counters) == CPS_OK);
    /* A quadratic takes Auto's closed forms: no exact preparation, one precision tried. */
    CHECK(counters.precision_attempts >= 1 && counters.isolation_passes == 0);
    CHECK(statistics.exact_factors == 0);
    {
        /* A quintic takes the general route, whose engines fill the counters (the closed forms answer degree 3). */
        cps_polynomial *quintic = NULL;
        cps_result *general = NULL;
        CHECK(cps_polynomial_create(&quintic) == CPS_OK);
        CHECK(cps_polynomial_set_coefficient(quintic, 0, "-2", "0") == CPS_OK &&
              cps_polynomial_set_coefficient(quintic, 5, "1", "0") == CPS_OK);
        CHECK(cps_solve(quintic, NULL, &general) == CPS_OK && cps_result_goal_met(general));
        CHECK(cps_result_get_statistics(general, &statistics) == CPS_OK);
        CHECK(cps_result_get_engine_counters(general, &counters) == CPS_OK);
        CHECK(counters.precision_attempts >= 1 && counters.isolation_passes >= 1);
        CHECK(statistics.exact_factors == 1);
        cps_result_destroy(general);
        cps_polynomial_destroy(quintic);
    }
    cps_root_view v;
    CHECK(cps_result_get_root(r, 0, &v) == CPS_OK);
    CHECK(v.count == 1 && v.multiplicity == 1 && v.realness == CPS_REAL);
    CHECK(strlen(v.real_lower.mantissa) > 0);
    CHECK(cps_result_get_root(r, 2, &v) == CPS_OUT_OF_RANGE);
    double re, im;
    CHECK(cps_result_get_approximation(r, 0, &re, &im) == CPS_OK);
    CHECK(re < -1.4 && im == 0);
    cps_result_destroy(r);
    cps_options *o = NULL;
    CHECK(cps_options_create(&o) == CPS_OK);
    CHECK(cps_options_set_method(NULL, CPS_METHOD_CLASSICAL) == CPS_INVALID_ARGUMENT);
    CHECK(cps_options_set_method(o, -1) == CPS_INVALID_ARGUMENT);
    CHECK(cps_options_set_method(o, 99) == CPS_INVALID_ARGUMENT);
    for (int method = CPS_METHOD_CLASSICAL; method <= CPS_METHOD_PIECEWISE; ++method) {
        CHECK(cps_options_set_method(o, method) == CPS_OK);
        CHECK(cps_solve(p, o, &r) == CPS_OK && cps_result_goal_met(r));
        CHECK(cps_result_size(r) == 2);
        cps_result_destroy(r);
    }
    CHECK(cps_options_set_method(o, CPS_METHOD_AUTO) == CPS_OK);
    CHECK(cps_options_set_goal(NULL, CPS_ROOT_APPROXIMATION) == CPS_INVALID_ARGUMENT);
    CHECK(cps_options_set_goal(o, 99) == CPS_INVALID_ARGUMENT);
    CHECK(cps_options_set_goal(o, -1) == CPS_INVALID_ARGUMENT);
    CHECK(cps_options_set_goal(o, CPS_CLUSTER_RADIUS) == CPS_INVALID_ARGUMENT);
    CHECK(cps_options_set_goal(o, CPS_ROOT_APPROXIMATION) == CPS_OK);
    CHECK(cps_options_set_goal(o, 99) == CPS_INVALID_ARGUMENT);
    CHECK(cps_solve(p, o, &r) == CPS_OK);
    CHECK(cps_result_goal(r) == CPS_ROOT_APPROXIMATION && cps_result_goal_met(r));
    cps_result_destroy(r);
    /* Exact coefficient isolation under Auto: all complex roots or real-only;
       the regions promise no width and are not approximations. */
    CHECK(cps_options_set_goal(o, CPS_ROOT_ISOLATION) == CPS_OK);
    CHECK(cps_solve(p, o, &r) == CPS_OK);
    CHECK(cps_result_goal(r) == CPS_ROOT_ISOLATION && cps_result_complete(r) && cps_result_goal_met(r));
    CHECK(cps_result_size(r) == 2);
    CHECK(cps_result_get_root(r, 0, &v) == CPS_OK && v.count == 1 && v.multiplicity == 1);
    CHECK(cps_result_get_approximation(r, 0, &re, &im) == CPS_UNRESOLVED);
    cps_result_destroy(r);
    CHECK(cps_options_set_real_only(o, 1) == CPS_OK);
    CHECK(cps_solve(p, o, &r) == CPS_OK);
    CHECK(cps_result_goal(r) == CPS_ROOT_ISOLATION && cps_result_complete(r) && cps_result_goal_met(r));
    CHECK(cps_result_size(r) == 2);
    CHECK(cps_result_get_root(r, 0, &v) == CPS_OK && v.count == 1 && v.multiplicity == 1 && v.realness == CPS_REAL);
    CHECK(cps_result_get_approximation(r, 0, &re, &im) == CPS_UNRESOLVED);
    cps_result_destroy(r);
    CHECK(cps_options_set_method(o, CPS_METHOD_SECULAR) == CPS_OK);
    CHECK(cps_solve(p, o, &r) == CPS_INVALID_ARGUMENT && r == NULL);
    CHECK(cps_options_set_method(o, CPS_METHOD_AUTO) == CPS_OK);
    /* A region: the setter's argument checks (each end's literal); ends out of order are the solve's usage error;
       [0, +inf) holds sqrt(2) alone; four unbounded sides clear it. */
    CHECK(cps_options_set_region(NULL, "0", "2", NULL, NULL) == CPS_INVALID_ARGUMENT);
    CHECK(cps_options_set_region(o, "2", "0", NULL, NULL) == CPS_OK);
    CHECK(cps_solve(p, o, &r) == CPS_INVALID_ARGUMENT && r == NULL && strstr(cps_last_error(), "lower end exceeds"));
    CHECK(cps_options_set_region(o, "x", "", NULL, NULL) == CPS_INVALID_ARGUMENT);
    CHECK(cps_options_set_region(o, "0", "", "", NULL) == CPS_OK);
    CHECK(cps_solve(p, o, &r) == CPS_OK && cps_result_goal_met(r) && cps_result_size(r) == 1);
    CHECK(cps_result_get_root(r, 0, &v) == CPS_OK && v.realness == CPS_REAL && v.real_lower.mantissa[0] != '-');
    cps_result_destroy(r);
    CHECK(cps_options_set_region(o, NULL, NULL, NULL, NULL) == CPS_OK);
    CHECK(cps_solve(p, o, &r) == CPS_OK && cps_result_goal_met(r) && cps_result_size(r) == 2);
    cps_result_destroy(r);
    CHECK(cps_options_set_real_only(o, 0) == CPS_OK);
    CHECK(cps_options_set_goal(o, CPS_ROOT_APPROXIMATION) == CPS_OK);
    CHECK(cps_options_set_precision(o, 0) == CPS_INVALID_ARGUMENT);
    {
        /* C19: a NULL imaginary part or radius is 0; a NULL real part is refused. */
        cps_polynomial *q = NULL;
        char text[16];
        CHECK(cps_polynomial_create(&q) == CPS_OK);
        CHECK(cps_polynomial_set_coefficient(q, 0, "-2", NULL) == CPS_OK);
        CHECK(cps_polynomial_set_ball_coefficient(q, 2, "1", NULL, NULL, NULL) == CPS_OK);
        CHECK(cps_polynomial_set_coefficient(q, 1, NULL, "0") == CPS_INVALID_ARGUMENT);
        CHECK(cps_polynomial_set_ball_coefficient(q, 1, NULL, NULL, NULL, NULL) == CPS_INVALID_ARGUMENT);
        CHECK(cps_polynomial_degree(q) == 2);
        CHECK(cps_polynomial_get_coefficient(q, 0, CPS_PART_IMAG, text, sizeof text, NULL) == CPS_OK && strcmp(text, "0") == 0);
        CHECK(cps_polynomial_get_coefficient(q, 2, CPS_PART_REAL_RADIUS, text, sizeof text, NULL) == CPS_OK &&
              strcmp(text, "0") == 0);
        cps_result *exact = NULL;
        CHECK(cps_solve(q, NULL, &exact) == CPS_OK && cps_result_goal_met(exact) && cps_result_size(exact) == 2);
        cps_result_destroy(exact);
        cps_polynomial_destroy(q);
        cps_secular_term_view real_term = {"1/3", NULL, "2", NULL};
        cps_secular *s = NULL;
        CHECK(cps_secular_create(&real_term, 1, NULL, &s) == CPS_OK && s);
        cps_secular_destroy(s);
        s = NULL;
        real_term.pole_real = NULL;
        CHECK(cps_secular_create(&real_term, 1, NULL, &s) == CPS_INVALID_ARGUMENT && !s);
    }
    CHECK(cps_options_set_cluster_radius(o, "-1") == CPS_INVALID_ARGUMENT);
    CHECK(cps_options_set_cluster_radius(o, NULL) == CPS_INVALID_ARGUMENT);
    CHECK(cps_options_set_cluster_radius(NULL, "1") == CPS_INVALID_ARGUMENT);
    CHECK(cps_options_set_cluster_radius(o, "1e3") == CPS_OK);
    CHECK(cps_options_set_goal(o, CPS_CLUSTER_RADIUS) == CPS_OK);
    CHECK(cps_options_set_cluster_radius(o, "invalid") == CPS_INVALID_ARGUMENT);
    CHECK(cps_solve(p, o, &r) == CPS_OK);
    CHECK(cps_result_goal(r) == CPS_CLUSTER_RADIUS && cps_result_goal_met(r));
    CHECK(cps_result_size(r) == 1);
    CHECK(cps_result_get_root(r, 0, &v) == CPS_OK && v.count == 2 && v.multiplicity == 0);
    CHECK(cps_result_get_approximation(r, 0, &re, &im) == CPS_UNRESOLVED);
    cps_result_destroy(r);
    CHECK(cps_options_set_goal(o, CPS_ROOT_APPROXIMATION) == CPS_OK);
    CHECK(cps_options_set_goal(o, CPS_CLUSTER_RADIUS) == CPS_INVALID_ARGUMENT);
    CHECK(cps_solve(p, o, &r) == CPS_OK);
    CHECK(cps_result_goal(r) == CPS_ROOT_APPROXIMATION && cps_result_goal_met(r));
    cps_result_destroy(r);
    CHECK(cps_options_set_cluster_radius(o, "0") == CPS_OK);
    /* A token requested before the solve: no work, a cancelled answer. */
    cps_cancel_token *stop = NULL;
    CHECK(cps_cancel_token_create(&stop) == CPS_OK);
    cps_cancel_token_request(stop);
    CHECK(cps_options_set_cancel_token(o, stop) == CPS_OK);
    CHECK(cps_solve(p, o, &r) == CPS_OK);
    CHECK(cps_result_termination(r) == CPS_CANCELLED);
    CHECK(!cps_result_complete(r));
    CHECK(cps_result_get_engine_counters(r, &counters) == CPS_OK);
    CHECK(counters.precision_attempts == 0 && counters.isolation_passes == 0);
    cps_result_destroy(r);
    r = (cps_result *)p;
    CHECK(cps_solve(NULL, o, &r) == CPS_INVALID_ARGUMENT && r == NULL);
    CHECK(cps_options_set_cancel_token(o, NULL) == CPS_OK);
    /* 0 selects a budget's default; a one-unit work budget stops after that unit (a quintic: the closed forms answer
       degree 3 and below within one unit). */
    CHECK(cps_options_set_live_cells(o, 0) == CPS_OK && cps_options_live_cells(o) == CPS_DEFAULT_LIVE_CELLS);
    CHECK(cps_options_set_work_budget(o, 0) == CPS_OK && cps_options_work_budget(o) == 0);
    CHECK(cps_options_set_work_budget(o, 1) == CPS_OK);
    {
        cps_polynomial *quintic = NULL;
        CHECK(cps_polynomial_create(&quintic) == CPS_OK);
        CHECK(cps_polynomial_set_coefficient(quintic, 0, "-2", "0") == CPS_OK &&
              cps_polynomial_set_coefficient(quintic, 5, "1", "0") == CPS_OK);
        CHECK(cps_solve(quintic, o, &r) == CPS_OK);
        cps_polynomial_destroy(quintic);
    }
    CHECK(cps_result_termination(r) == CPS_WORK_LIMIT && cps_result_stop_reason(r) == CPS_STOP_WORK_BUDGET);
    CHECK(cps_result_get_statistics(r, &statistics) == CPS_OK && statistics.work_units == 1);
    CHECK(cps_result_get_engine_counters(r, NULL) == CPS_INVALID_ARGUMENT);
    CHECK(cps_result_get_engine_counters(NULL, &counters) == CPS_INVALID_ARGUMENT);
    cps_result_destroy(r);
    CHECK(cps_solve(p, o, &r) == CPS_OK && cps_result_goal_met(r));
    CHECK(cps_result_get_statistics(r, &statistics) == CPS_OK && statistics.work_units == 1);
    cps_result_destroy(r);
    CHECK(cps_options_set_work_budget(o, 100000) == CPS_OK);
    CHECK(cps_solve(p, o, &r) == CPS_OK && cps_result_goal_met(r));
    CHECK(cps_result_get_root(r, 0, &v) == CPS_OK && v.count == 1);
    cps_result_destroy(r);
    cps_secular *secular = NULL;
    cps_secular_term_view term = {"1/3", "2/7", "2", "-1"};
    CHECK(cps_secular_create(&term, 1, NULL, &secular) == CPS_OK);
    CHECK(cps_solve_secular(secular, NULL, &r) == CPS_OK && cps_result_goal_met(r));
    cps_secular_destroy(secular);
    CHECK(cps_result_get_engine_counters(r, &counters) == CPS_OK);
    cps_result_destroy(r);
    CHECK(cps_secular_create(NULL, 0, NULL, &secular) == CPS_OK);
    CHECK(cps_solve_secular(secular, NULL, &r) == CPS_OK && cps_result_goal_met(r));
    CHECK(cps_result_degree(r) == 0 && cps_result_size(r) == 0);
    cps_result_destroy(r);
    cps_secular_destroy(secular);
    CHECK(cps_secular_create(NULL, 1, NULL, &secular) == CPS_INVALID_ARGUMENT && !secular);
    CHECK(cps_secular_create(&term, SIZE_MAX, NULL, &secular) == CPS_OUT_OF_RANGE && !secular);
    cps_secular_term_view large_terms[2049];
    char large_poles[2049][16];
    for (size_t i = 0; i < 2049; ++i) {
        snprintf(large_poles[i], sizeof(large_poles[i]), "%zu", i);
        large_terms[i] = (cps_secular_term_view){"1", "0", large_poles[i], "0"};
    }
    CHECK(cps_secular_create(large_terms, 2049, NULL, &secular) == CPS_OK && secular);
    CHECK(cps_options_set_cancel_token(o, stop) == CPS_OK);
    CHECK(cps_solve_secular(secular, o, &r) == CPS_OK && cps_result_termination(r) == CPS_CANCELLED);
    cps_cancel_token_destroy(stop);
    cps_secular_destroy(secular);
    cps_result_destroy(r);
    term.residue_real = term.residue_imag = "0";
    CHECK(cps_secular_create(&term, 1, NULL, &secular) == CPS_INVALID_ARGUMENT && !secular);
    cps_secular_term_view repeated_poles[2] = {{"1", "0", "1/2", "0"}, {"2", "0", "0.5", "0"}};
    CHECK(cps_secular_create(repeated_poles, 2, NULL, &secular) == CPS_INVALID_ARGUMENT && !secular);
    CHECK(cps_solve_secular(NULL, NULL, &r) == CPS_INVALID_ARGUMENT && !r);
    cps_secular_destroy(NULL);
    /* Program input: x^2 - 2 built node by node; solves and the refusals. */
    cps_program *program = NULL;
    uint64_t x = 0, two = 0, square = 0, difference = 0;
    CHECK(cps_program_create(&program) == CPS_OK && program);
    CHECK(cps_program_variable(program, &x) == CPS_OK && x == 0);
    CHECK(cps_program_constant(program, "2", NULL, &two) == CPS_OK && two == 1);
    CHECK(cps_program_multiply(program, x, x, &square) == CPS_OK && square == 2);
    CHECK(cps_program_subtract(program, square, two, &difference) == CPS_OK && difference == 3);
    /* C20: no output is CPS_NO_NODE, and so is the node a failed builder call writes; the program is unchanged. */
    CHECK(cps_program_output(program) == CPS_NO_NODE);
    uint64_t failed = 0;
    CHECK(cps_program_add(program, x, 99, &failed) == CPS_INVALID_ARGUMENT && failed == CPS_NO_NODE);
    failed = 0;
    CHECK(cps_program_constant(program, "x", NULL, &failed) == CPS_INVALID_ARGUMENT && failed == CPS_NO_NODE);
    CHECK(cps_program_size(program) == 4);
    CHECK(cps_solve_program(program, NULL, &r) == CPS_INVALID_ARGUMENT && r == NULL); /* no output yet */
    CHECK(cps_program_set_output(program, difference) == CPS_OK && cps_program_output(program) == difference);
    CHECK(cps_solve_program(program, NULL, &r) == CPS_OK && cps_result_complete(r) && cps_result_goal_met(r));
    CHECK(cps_result_goal(r) == CPS_ROOT_APPROXIMATION && cps_result_degree(r) == 2 && cps_result_size(r) == 2);
    /* Every constant real: the roots +-sqrt(2) are proved real by conjugation (the program note, certpolysolve.h). */
    CHECK(cps_result_get_root(r, 0, &v) == CPS_OK && v.count == 1 && v.realness == CPS_REAL);
    CHECK(cps_result_get_approximation(r, 0, &re, &im) == CPS_OK && re < -1.414 && re > -1.4143 && fabs(im) < 1e-15);
    CHECK(cps_result_get_approximation(r, 1, &re, &im) == CPS_OK && re > 1.414 && re < 1.4143);
    cps_result_destroy(r);
    cps_options *po = NULL;
    CHECK(cps_options_create(&po) == CPS_OK);
    CHECK(cps_options_set_goal(po, CPS_ROOT_ACCURACY) == CPS_OK);
    CHECK(cps_solve_program(program, po, &r) == CPS_OK && cps_result_goal_met(r) && cps_result_goal(r) == CPS_ROOT_ACCURACY);
    CHECK(cps_result_get_root(r, 1, &v) == CPS_OK && v.count == 1 && v.multiplicity == 1);
    cps_result_destroy(r);
    /* Real-only and piecewise programs: usage errors, no result. */
    CHECK(cps_options_set_real_only(po, 1) == CPS_OK);
    CHECK(cps_solve_program(program, po, &r) == CPS_INVALID_ARGUMENT && r == NULL);
    CHECK(cps_options_set_real_only(po, 0) == CPS_OK);
    CHECK(cps_options_set_method(po, CPS_METHOD_PIECEWISE) == CPS_OK);
    CHECK(cps_solve_program(program, po, &r) == CPS_INVALID_ARGUMENT && r == NULL);
    CHECK(cps_options_set_method(po, CPS_METHOD_SECULAR) == CPS_OK);
    CHECK(cps_solve_program(program, po, &r) == CPS_OK && cps_result_goal_met(r) && cps_result_size(r) == 2);
    cps_result_destroy(r);
    cps_program_destroy(program);
    CHECK(cps_solve_program(NULL, NULL, &r) == CPS_INVALID_ARGUMENT && !r);
    cps_options_destroy(po);
    cps_options_destroy(o);
    cps_polynomial_destroy(p);
    cps_result_destroy(NULL);
    cps_options_destroy(NULL);
    cps_polynomial_destroy(NULL);
    cps_cleanup_thread();
    puts("C ABI tests passed");
    return 0;
}
