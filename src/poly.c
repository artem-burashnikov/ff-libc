#include "poly.h"

#include <stddef.h>
#include <stdlib.h>

poly_t *poly_init(int8_t degree, int8_t *coeff) {
  poly_t *poly = malloc(sizeof(*poly));

  if (!poly || !coeff || (degree < 0)) {
    return NULL;
  }

  poly->deg = degree;
  poly->coeff = coeff;

  return poly;
}

void poly_destroy(poly_t *poly) {
  free(poly->coeff);
  free(poly);
  return;
}

int poly_eq(poly_t *a, poly_t *b) {
  if (!a || !b || (a->deg != b->deg)) {
    return 0;
  }

  for (size_t i = 0; i < a->deg; ++i) {
    if (a->coeff[i] == b->coeff[i]) {
      continue;
    }
    return 0;
  }

  return 1;
}
