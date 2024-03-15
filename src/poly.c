#include "poly.h"

#include <stddef.h>
#include <stdlib.h>

#define MIN(A, B) ((A) < (B) ? (A) : (B))

poly_t *poly_init(int8_t degree, int8_t *coeff, int8_t len) {
  poly_t *poly = malloc(sizeof(*poly));

  if (!poly || !coeff || (len < 1) || (degree < 0) || (len < (degree + 1))) {
    poly_destroy(poly);
    return NULL;
  }

  poly->deg = degree;
  poly->coeff = coeff;
  poly->len = len;

  return poly;
}

void poly_destroy(poly_t *poly) {
  free(poly->coeff);
  free(poly);
  return;
}

int poly_eq(const poly_t *a, const poly_t *b) {
  if (!a || !b || (a->deg != b->deg)) {
    return 0;
  }

  int8_t n = MIN(a->len, b->len);

  for (size_t i = 0; i < n; ++i) {
    if (a->coeff[i] == b->coeff[i]) {
      continue;
    }
    return 0;
  }

  return 1;
}
