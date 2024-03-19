#include <assert.h>
#include <stdint.h>
#include <string.h>

#include "GF.h"
#include "poly.h"
#include "utils.h"

void GF_elem_diff(GF_elem_t *res, GF_elem_t a, GF_elem_t b) {
  if (!res) {
    return;
  }
  // Get b's complement.
  GF_elem_t *negb = GF_elem_get_complement(b);

  // Now calculate res = a + ~b. res is normalized in sum.
  GF_elem_sum(res, a, *negb);

  GF_elem_destroy(negb);
}
