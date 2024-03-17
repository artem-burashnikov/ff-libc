#include <assert.h>
#include <stddef.h>

#include "GF.h"
#include "poly.h"
#include "utils.h"

int GF_elem_sum(GF_elem_t *res, GF_elem_t a, GF_elem_t b) {
  if (!res || res->poly->len!= a.poly->len || a.poly->len != b.poly->len) {
    return 1;
  }

  // Different fields.
  if (!GF_eq(res->GF, a.GF) && !GF_eq(res->GF, b.GF)) {
    return 1;
  }

  poly_carryless_sum(res->poly, *a.poly, *b.poly, res->GF->p);
  
  return 0;
}
