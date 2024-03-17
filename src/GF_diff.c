#include <assert.h>
#include <stdint.h>
#include <string.h>

#include "GF.h"
#include "poly.h"
#include "utils.h"

int GF_elem_diff(GF_elem_t *res, GF_elem_t a, GF_elem_t b) {
  if (!res) {
    return 1;
  }
  // Set b coefficients to their complements mod p.
  GF_elem_get_complement(&b, b);

  // Now calculate res = a + ~b.
  GF_elem_sum(res, a, b);

  // Complement b again.
  GF_elem_get_complement(&b, b);

  poly_normalize_deg(res);
  return 0;
}
