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

  poly_t *tmp = poly_create_zero(a.poly->deg + b.poly->deg + 1);

  if (!tmp) {
    return 1;
  }

  poly_carryless_mul(tmp, *a.poly, *b.poly, res->GF->p);
  poly_normalize_deg(tmp);
  poly_carryless_div(tmp, *res->GF->I, res->GF->p);

  memcpy(res->poly->coeff, tmp->coeff, sizeof(*tmp->coeff) * (tmp->deg + 1));
  res->poly->deg = tmp->deg;

  poly_destroy(tmp);

  return 0;
}
