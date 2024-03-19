#include <assert.h>
#include <stddef.h>

#include "GF.h"
#include "poly.h"
#include "utils.h"

void GF_elem_sum(GF_elem_t *res, GF_elem_t a, GF_elem_t b) {
  if (!res) {
    return;
  }

  // Different fields.
  if (!GF_eq(res->GF, a.GF) && !GF_eq(res->GF, b.GF)) {
    return;
  }

  poly_carryless_sum(res->poly, *a.poly, *b.poly, res->GF->p);

  GF_elem_normalize(res);
}
