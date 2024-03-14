#include "GF.h"

#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include "poly.h"

GF_t *GF_init(int8_t p, int8_t n, poly_t *I) {
  GF_t *GF = malloc(sizeof(*GF));

  if (!GF || !I || (n < 2) || (n > 100) || (p > 11) || (p < 2)) {
    return NULL;
  }

  GF->p = p;
  GF->n = n;
  GF->I = I;

  return GF;
}

void GF_destroy(GF_t *GF) {
  poly_destroy(GF->I);
  free(GF);
  return;
}

void GFelement_destroy(GFelement_t *a) {
  poly_destroy(a->poly);
  return;
}

int GF_eq(GF_t *f, GF_t *k) {
  // Irreducible polynomials must match.
  int8_t ret = poly_eq(f->I, k->I);

  // Characteristics of fields and dimensions of extensions must match.
  if (f->p != k->p || f->n != k->n) {
    ret = 0;
  }

  return ret;
}

int GFelement_add(GFelement_t *res, GFelement_t *a, GFelement_t *b) {
  int8_t p, n;
  int8_t *u, *v, *w;

  // Invalid input.
  if (!res || !a || !b) {
    return 1;
  }

  // Different fields.
  if (!GF_eq(res->GF, a->GF) && !GF_eq(res->GF, b->GF)) {
    return 1;
  }

  // Dimension of the field extension.
  n = a->GF->n;

  // Characteristic of the field.
  p = res->GF->p;

  w = res->poly->coeff;
  u = a->poly->coeff;
  v = b->poly->coeff;

  for (size_t i = 0; i < n; ++i) {
    w[i] = (u[i] + v[i]) % p;
  }

  return 0;
}

GFelement_t *GFelement_get_neutral(GF_t *GF) {
  GFelement_t *zero;
  poly_t *poly;
  int8_t *coeff;
  int8_t dim;

  if (!GF) {
    return NULL;
  }

  dim = GF->n - 1;
  coeff = calloc(dim, sizeof(*coeff));
  poly = poly_init(0, coeff);

  if (!coeff || !poly) {
    return NULL;
  }

  zero->GF = GF;
  zero->poly = poly;

  return zero;
}

GFelement_t *GFelement_get_unity(GF_t *GF) {
  GFelement_t *unity;

  /* Get neutral and set the least significant digit to one. */
  unity = GFelement_get_neutral(GF);

  if (!unity || (!GF->n < 1)) {
    return NULL;
  }

  *unity->poly->coeff = 1;

  return unity;
}
