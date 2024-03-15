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

void GFelement_destroy(GFelement_t *a) {
  poly_destroy(a->poly);
  free(a);
  return;
}

int GF_eq(const GF_t *f, const GF_t *k) {
  // Irreducible polynomials must match.
  int8_t ret = poly_eq(f->I, k->I);

  // Characteristics of fields and dimensions of extensions must match.
  if (f->p != k->p || f->n != k->n) {
    ret = 0;
  }

  return ret;
}

int GFelement_add(GFelement_t *res, const GFelement_t *a,
                  const GFelement_t *b) {
  // Invalid input.
  if (!res || !a || !b) {
    return 1;
  }

  // Different fields.
  if (!GF_eq(res->GF, a->GF) && !GF_eq(res->GF, b->GF)) {
    return 1;
  }

  // Dimension of the field extension.
  int8_t n = a->GF->n;

  // Characteristic of the field.
  int8_t p = res->GF->p;

  int8_t *w = res->poly->coeff;
  int8_t *u = a->poly->coeff;
  int8_t *v = b->poly->coeff;

  for (size_t i = 0; i < n; ++i) {
    w[i] = (u[i] + v[i]) % p;
  }

  return 0;
}

GFelement_t *GFelement_get_neutral(GF_t *GF) {
  if (!GF) {
    return NULL;
  }

  int8_t dim = GF->n;

  GFelement_t *ret = malloc(sizeof(*ret));
  int8_t *coeff = calloc(dim, sizeof(*coeff));
  poly_t *poly = poly_init(0, coeff, dim);

  if (!ret || !coeff || !poly) {
    free(ret);
    free(coeff);
    poly_destroy(poly);
    return NULL;
  }

  ret->GF = GF;
  ret->poly = poly;

  return ret;
}

GFelement_t *GFelement_get_unity(GF_t *GF) {
  /* Get neutral and set the least significant digit to one. */
  GFelement_t *unity = GFelement_get_neutral(GF);

  if (!unity || !GF || (GF->n < 1)) {
    GFelement_destroy(unity);
    return NULL;
  }

  unity->poly->coeff[0] = 1;

  return unity;
}
