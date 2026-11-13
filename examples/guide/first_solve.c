/* The first solve of docs/api.md, section 1, in C. The suite builds and runs it (guide_first_solve_c). */
#include <certpolysolve/certpolysolve.h>
#include <stdio.h>

int main(void) {
    cps_polynomial *p = NULL;
    cps_options *o = NULL;
    cps_result *r = NULL;
    int status = 1;
    if (cps_polynomial_create(&p) != CPS_OK ||
        cps_polynomial_set_coefficient(p, 0, "-2", NULL) != CPS_OK || /* x^2 - 2 */
        cps_polynomial_set_coefficient(p, 2, "1", NULL) != CPS_OK)
        goto done;
    if (cps_options_create(&o) != CPS_OK || cps_options_set_precision(o, 128) != CPS_OK ||
        cps_options_set_time_limit(o, 10.0) != CPS_OK)
        goto done;
    if (cps_solve(p, o, &r) != CPS_OK) { /* a usage error: no result */
        fprintf(stderr, "%s\n", cps_last_error());
        goto done;
    }
    if (!cps_result_goal_met(r)) { /* a stop: the regions are proved, the goal is not met */
        fprintf(stderr, "%s: %s\n", cps_termination_name(cps_result_termination(r)), cps_result_diagnostics(r));
        goto done;
    }
    for (size_t i = 0; i < cps_result_size(r); ++i) {
        cps_root_view v;
        char lower[64], upper[64];
        if (cps_result_get_root(r, i, &v) != CPS_OK ||
            cps_dyadic_format(v.real_lower, 40, lower, sizeof lower, NULL) != CPS_OK ||
            cps_dyadic_format(v.real_upper, 40, upper, sizeof upper, NULL) != CPS_OK)
            goto done;
        printf("[%s, %s] holds %llu root(s)\n", lower, upper, (unsigned long long)v.count);
    }
    status = 0;
done:
    cps_result_destroy(r);
    cps_options_destroy(o);
    cps_polynomial_destroy(p);
    return status;
}
