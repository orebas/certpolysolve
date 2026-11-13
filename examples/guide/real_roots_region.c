/* The API guide's real roots in a region (docs/api.md, section 7): the real roots of 16x^5 - 20x^3 + 5x - 1/2
   between 0 and 1, each in an interval proved to hold it. */
#include <certpolysolve/certpolysolve.h>
#include <stdio.h>

int main(void) {
    cps_polynomial *p = NULL;
    cps_options *o = NULL;
    cps_result *r = NULL;
    int status = 1;
    if (cps_polynomial_create(&p) != CPS_OK ||
        cps_polynomial_set_coefficient(p, 0, "-1/2", NULL) != CPS_OK ||
        cps_polynomial_set_coefficient(p, 1, "5", NULL) != CPS_OK ||
        cps_polynomial_set_coefficient(p, 3, "-20", NULL) != CPS_OK ||
        cps_polynomial_set_coefficient(p, 5, "16", NULL) != CPS_OK)
        goto done;
    if (cps_options_create(&o) != CPS_OK || cps_options_set_real_only(o, 1) != CPS_OK ||
        cps_options_set_region(o, "0", "1", NULL, NULL) != CPS_OK) /* 0 <= x <= 1; NULL: no bound */
        goto done;
    if (cps_solve(p, o, &r) != CPS_OK) {
        fprintf(stderr, "%s\n", cps_last_error());
        goto done;
    }
    if (!cps_result_goal_met(r)) {
        fprintf(stderr, "%s: %s\n", cps_termination_name(cps_result_termination(r)), cps_result_diagnostics(r));
        goto done;
    }
    for (size_t i = 0; i < cps_result_size(r); ++i) {
        cps_root_view v;
        char lower[32], upper[32];
        if (cps_result_get_root(r, i, &v) != CPS_OK ||
            cps_dyadic_format(v.real_lower, 20, lower, sizeof lower, NULL) != CPS_OK ||
            cps_dyadic_format(v.real_upper, 20, upper, sizeof upper, NULL) != CPS_OK)
            goto done;
        printf("[%s, %s]\n", lower, upper);
    }
    status = 0;
done:
    cps_result_destroy(r);
    cps_options_destroy(o);
    cps_polynomial_destroy(p);
    return status;
}
