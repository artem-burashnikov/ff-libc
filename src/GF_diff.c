#include <assert.h>
#include <stdint.h>
#include <string.h>

#include "GF.h"
#include "poly.h"
#include "utils.h"

int GF_elem_diff(GF_elem_t *res, GF_elem_t a, GF_elem_t b) {
  // Copy b coefficients to res.
  memcpy(res->poly->coeff, b.poly->coeff, b.poly->len);

  // Set res coefficients to their complements mod p.
  if (!res->poly || GF_elem_get_complement(res, b)) {
    return 1;
  }

  // Now calculate res = a + complement(b) which is res = a + res.
  return GF_elem_sum(res, a, *res);
}
