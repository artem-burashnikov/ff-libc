#include <assert.h>
#include <stdint.h>
#include <string.h>

#include "GF.h"
#include "poly.h"
#include "utils.h"

int GF_elem_div(GF_elem_t *res, GF_elem_t a, GF_elem_t b) {
  if (!res || ((b.poly->deg == 0) && (*b.poly->coeff == 0))) {
    return 1;
  }

  // Different fields.
  if (!GF_eq(res->GF, a.GF) && !GF_eq(res->GF, b.GF)) {
    return 1;
  }

  GF_elem_t *inv_b = GF_elem_get_inverse(b);
  if (!inv_b) {
    return 1;
  }

  GF_elem_prod(res, a, *inv_b);

  GF_elem_destroy(inv_b);

  return 0;
}
