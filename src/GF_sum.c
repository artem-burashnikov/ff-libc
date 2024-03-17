#include <assert.h>
#include <stddef.h>

#include "GF.h"
#include "poly.h"
#include "utils.h"

int GF_elem_sum(GF_elem_t *res, GF_elem_t a, GF_elem_t b) {
  if (!res) {
    return 1;
  }

  // Different fields.
  if (!GF_eq(res->GF, a.GF) && !GF_eq(res->GF, b.GF)) {
    return 1;
  }

  for (size_t i = 0; i < res->GF->n; ++i) {
    res->poly->coeff[i] = (a.poly->coeff[i] + b.poly->coeff[i]) % res->GF->p;
  }

  poly_normalize_coeff(res->poly, res->GF->p);
  poly_normalize_deg(res->poly);
  return 0;
}
