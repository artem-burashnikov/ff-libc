#include "poly.h"

#include <stddef.h>
#include <stdlib.h>

poly_t *poly_init(int8_t degree, int8_t *coeff) {
  poly_t *poly = malloc(sizeof(*poly));

  if (!poly || !coeff) {
    return NULL;
  }
  
  poly->deg = (degree > 0) ? degree : MINUS_INFINITY;
  poly->coeff = coeff;

  return poly;
}

void poly_destroy(poly_t *poly) {
  free(poly->coeff);
  free(poly);
  return;
}
