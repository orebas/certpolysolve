#include <flint/acb_poly.h>
int main(void) {
    acb_poly_t polynomial;
    acb_t root;
    acb_poly_init(polynomial);
    acb_init(root);
    acb_poly_set_coeff_si(polynomial, 0, -1);
    acb_poly_set_coeff_si(polynomial, 1, 1);
    acb_poly_find_roots(root, polynomial, 0, 1, 64);
    acb_clear(root);
    acb_poly_clear(polynomial);
    flint_cleanup();
    return 0; /* With the shim this program MUST instead exit 90. */
}
