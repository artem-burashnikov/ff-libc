#include "poly.h"

#include <stddef.h>
#include <stdlib.h>

int poly_init(poly_t *poly, int8_t degree, int8_t *coeff) {
  if (!coeff || (degree <= 0)) {
    return 1;
  }

  poly->deg = degree;
  poly->coeff = coeff;

  return 0;
}

void poly_destroy(poly_t *poly) {
  free(poly->coeff);
  free(poly);
  return;
}
