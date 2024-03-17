#include <assert.h>
#include <stdint.h>
#include <string.h>

#include "GF.h"
#include "poly.h"
#include "utils.h"

int GF_elem_prod(GF_elem_t *res, GF_elem_t a, GF_elem_t b) {
  if (!res || a.poly->len != b.poly->len) {
    return 1;
  }

  // Different fields.
  if (!GF_eq(res->GF, a.GF) && !GF_eq(res->GF, b.GF)) {
    return 1;
  }

  poly_t *res_poly = poly_create_zero(a.poly->len + b.poly->len - 2);
  
  if (!res_poly) {
    return 1;
  }
  
  poly_carryless_mul(res_poly, *a.poly, *b.poly, *res->GF->I, res->GF->p);

  memmove(res->poly->coeff, res_poly->coeff, sizeof(*res_poly->coeff) * res->poly->len);

  poly_destroy(res_poly);

  GF_elem_normalize(res);
  
  return 0;
}
