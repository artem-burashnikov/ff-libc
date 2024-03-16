#include <assert.h>
#include <stdint.h>

#include "GF.h"
#include "poly.h"
#include "utils.h"

int GF_elem_sum(GF_elem_t *res, GF_elem_t a, GF_elem_t b) {
  GF_elem_t h = *res;

  // Different fields.
  if (!GF_eq(h.GF, a.GF) && !GF_eq(h.GF, b.GF)) {
    return 1;
  }

  // Dimension of the field extension.
  int8_t n = a.GF->n;

  // Characteristic of the field.
  int8_t p = h.GF->p;

  int8_t *w = h.poly->coeff;
  int8_t *u = a.poly->coeff;
  int8_t *v = b.poly->coeff;

  for (int8_t i = 0; i < n; ++i) {
    assert(u[i] >= 0);
    assert(v[i] >= 0);
    w[i] = (u[i] + v[i]) % p;
  }

  poly_normalize_deg(res->poly);
  return 0;
}
