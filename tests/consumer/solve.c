#include <certpolysolve/certpolysolve.h>

// The installed C interface: a piecewise solve and a native secular solve, each one call.
int main(void) {
    cps_polynomial *polynomial = NULL;
    cps_options *options = NULL;
    cps_result *result = NULL;
    cps_secular *secular = NULL;
    int success = 0;
    if (cps_polynomial_create(&polynomial) != CPS_OK || cps_options_create(&options) != CPS_OK)
        goto cleanup;
    if (cps_polynomial_set_coefficient(polynomial, 0, "-2", "0") != CPS_OK ||
        cps_polynomial_set_coefficient(polynomial, 2, "1", "0") != CPS_OK)
        goto cleanup;
    if (cps_options_set_method(options, CPS_METHOD_PIECEWISE) != CPS_OK ||
        cps_options_set_work_budget(options, 100000) != CPS_OK || cps_options_set_live_cells(options, 65536) != CPS_OK ||
        cps_solve(polynomial, options, &result) != CPS_OK)
        goto cleanup;
    cps_statistics statistics;
    statistics.size = sizeof statistics;
    cps_engine_counters counters;
    counters.size = sizeof counters;
    success = cps_result_goal_met(result) && cps_result_complete(result) &&
              cps_result_get_statistics(result, &statistics) == CPS_OK && statistics.work_units > 0 &&
              cps_result_stop_reason(result) == CPS_STOP_NONE;
    if (!success) goto cleanup;
    cps_result_destroy(result);
    result = NULL;
    const cps_secular_term_view term = {"1/3", "0", "0", "0"};
    if (cps_secular_create(&term, 1, NULL, &secular) != CPS_OK || cps_solve_secular(secular, NULL, &result) != CPS_OK) {
        success = 0;
        goto cleanup;
    }
    success = cps_result_goal_met(result) && cps_result_size(result) == 1 &&
              cps_result_get_engine_counters(result, &counters) == CPS_OK;
cleanup:
    cps_secular_destroy(secular);
    cps_result_destroy(result);
    cps_polynomial_destroy(polynomial);
    cps_options_destroy(options);
    cps_cleanup_thread();
    return success ? 0 : 1;
}
