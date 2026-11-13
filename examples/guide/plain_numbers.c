/* The API guide's plain numbers (docs/api.md, section 8): the roots of x^2 - 2 as doubles, as the README's first
   program prints them in C++. */
#include <certpolysolve/certpolysolve.h>
#include <stdio.h>

int main(void) {
    cps_polynomial *p = NULL;
    cps_result *r = NULL;
    int status = 1;
    if (cps_polynomial_create(&p) != CPS_OK ||
        cps_polynomial_set_coefficient(p, 0, "-2", NULL) != CPS_OK || /* constant term first: x^2 - 2 */
        cps_polynomial_set_coefficient(p, 2, "1", NULL) != CPS_OK)
        goto done;
    if (cps_solve(p, NULL, &r) != CPS_OK) { /* NULL options: the defaults */
        fprintf(stderr, "%s\n", cps_last_error());
        goto done;
    }
    if (!cps_result_goal_met(r)) {
        fprintf(stderr, "%s: %s\n", cps_termination_name(cps_result_termination(r)), cps_result_diagnostics(r));
        goto done;
    }
    for (size_t i = 0; i < cps_result_size(r); ++i) {
        cps_root_view v;
        double re, im;
        if (cps_result_get_root(r, i, &v) != CPS_OK || cps_result_get_approximation(r, i, &re, &im) != CPS_OK)
            goto done;
        for (uint64_t k = 0; k < v.count; ++k) { /* a double root prints twice */
            if (im == 0)
                printf("%.17g\n", re);
            else
                printf("%.17g %+.17gi\n", re, im);
        }
    }
    status = 0;
done:
    cps_result_destroy(r);
    cps_polynomial_destroy(p);
    return status;
}
