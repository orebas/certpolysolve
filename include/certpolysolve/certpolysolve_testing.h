#ifndef CERTPOLYSOLVE_CERTPOLYSOLVE_TESTING_H
#define CERTPOLYSOLVE_CERTPOLYSOLVE_TESTING_H

/* For the library's own tests and fuzz harness, not for applications (not installed): the stops at an exact
   cancellation poll that a callback counting its calls gave until stage 3b-3 removed the callback
   (docs/api-spec-0.5.md). */
#include <certpolysolve/certpolysolve.h>

#ifdef __cplusplus
extern "C" {
#endif

/* A solve call with these options counts its cancellation polls from zero, and from its k-th poll on (k >= 1)
   answers "cancelled" (CPS_CANCELLED with CPS_STOP_CANCELLED), as a callback that returned nonzero from its k-th call
   on did. A requested token or a reached time limit answers a poll first, and that poll is not counted. k = 0, the
   default, turns counting off; k = UINT64_MAX counts without stopping. */
CPS_API cps_error cps_testing_options_stop_at_poll(cps_options *options, uint64_t k);
/* The cancellation polls the call made while counting was on; 0 when it was off (NULL: 0). */
CPS_API uint64_t cps_testing_result_polls(const cps_result *result);
/* With poll counting on (above), a solve call runs on one thread whatever cps_options_set_workers asks: the counter
   is not shared between threads. */
/* The jobs a split loop ran on the call's workers (cps_options_set_workers), by site: 0 the multiprecision secular
   packet's sweep, 1 the binary64 engine's evaluations, 2 a Taylor row, 3 the fixed-point transform, 4 the
   finisher's blocked evaluation. 0 with one worker, for NULL, and for a site past the last. The tests read it to
   show that an input reaches the split it means to. */
#define CPS_TESTING_SPLIT_SITES 5u
CPS_API uint64_t cps_testing_result_split_jobs(const cps_result *result, unsigned site);

#ifdef __cplusplus
}
#endif
#endif
