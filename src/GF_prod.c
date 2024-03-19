#include <assert.h>
#include <stdint.h>
#include <string.h>

#include "GF.h"
#include "poly.h"
#include "utils.h"

void GF_elem_prod(GF_elem_t *res, GF_elem_t a, GF_elem_t b) {
  if (!res) {
    return;
  }

  // Different fields.
  if (!GF_eq(res->GF, a.GF) && !GF_eq(res->GF, b.GF)) {
    return;
  }

  poly_t *tmp = poly_create_zero(a.poly->deg + b.poly->deg + 1);
  if (!tmp) {
    return;
  }

  poly_carryless_mul(tmp, *a.poly, *b.poly, res->GF->p);
  poly_normalize_deg(tmp);
  poly_carryless_div(tmp, *res->GF->I, res->GF->p);

  memcpy(res->poly->coeff, tmp->coeff, sizeof(*tmp->coeff) * (tmp->deg + 1));
  res->poly->deg = tmp->deg;

  poly_destroy(tmp);
}
