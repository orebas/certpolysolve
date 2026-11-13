/* ELF interposition test only; never linked into the installed library. */
#include <flint/acb_poly.h>
#include <flint/arb_fmpz_poly.h>
#include <flint/fmpz_poly.h>
#include <stdio.h>
#include <stdlib.h>
static _Noreturn void forbidden(const char *name) {
    fprintf(stderr, "forbidden backend root solver: %s\n", name);
    _Exit(90);
}
slong acb_poly_find_roots(acb_ptr a, const acb_poly_t b, acb_srcptr c, slong d, slong e) {
    (void)a;
    (void)b;
    (void)c;
    (void)d;
    (void)e;
    forbidden("acb_poly_find_roots");
}
slong _acb_poly_find_roots(acb_ptr a, acb_srcptr b, acb_srcptr c, slong d, slong e, slong f) {
    (void)a;
    (void)b;
    (void)c;
    (void)d;
    (void)e;
    (void)f;
    forbidden("_acb_poly_find_roots");
}
double _acb_poly_find_roots_double(acb_ptr a, acb_srcptr b, acb_srcptr c, slong d, slong e, slong f) {
    (void)a;
    (void)b;
    (void)c;
    (void)d;
    (void)e;
    (void)f;
    forbidden("_acb_poly_find_roots_double");
}
double cd_poly_find_roots(double *a, const double *b, const double *c, slong d, slong e, double f) {
    (void)a;
    (void)b;
    (void)c;
    (void)d;
    (void)e;
    (void)f;
    forbidden("cd_poly_find_roots");
}
slong arb_fmpz_poly_real_roots(arb_ptr a, const fmpz_poly_t b, int c, slong d) {
    (void)a;
    (void)b;
    (void)c;
    (void)d;
    forbidden("arb_fmpz_poly_real_roots");
}
void arb_fmpz_poly_complex_roots(acb_ptr a, const fmpz_poly_t b, int c, slong d) {
    (void)a;
    (void)b;
    (void)c;
    (void)d;
    forbidden("arb_fmpz_poly_complex_roots");
}
void fmpz_poly_isolate_real_roots(fmpq *a, slong *b, fmpz *c, slong *d, slong *e, const fmpz_poly_t f) {
    (void)a;
    (void)b;
    (void)c;
    (void)d;
    (void)e;
    (void)f;
    forbidden("fmpz_poly_isolate_real_roots");
}
