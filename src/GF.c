#include "GF.h"

#include <assert.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include "poly.h"
#include "utils.h"

GF_t *GF_init(int8_t p, int8_t n, poly_t *I) {
  GF_t *GF = malloc(sizeof(*GF));
  if (!GF) {
    return NULL;
  }

  // Artificial constraints.
  if ((n < 2) || (n > 100) || (p > 11) || (p < 2)) {
    return NULL;
  }

  if (!I || !I->coeff || (I->deg >= I->len)) {
    free(GF);
    return NULL;
  }

  GF->p = p;
  GF->n = n;
  GF->I = I;

  return GF;
}

void GF_elem_destroy(GF_elem_t *a) {
  poly_destroy(a->poly);
  free(a);
}

void GF_elem_normalize(GF_elem_t *a) {
  if (!a) {
    return;
  }
  // First set coefficients mod p.
  poly_normalize_coeff(a->GF->p, a->poly);
  // Then correct the polynomial degree if necessary.
  poly_normalize_deg(a->poly);
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

GF_elem_t *GF_elem_get_neutral(GF_t *GF) {
  if (!GF) {
    return NULL;
  }

  int8_t dim = GF->n;

  GF_elem_t *neutral = malloc(sizeof(*neutral));
  int8_t *coeff = calloc(dim, sizeof(*coeff));
  poly_t *poly = poly_from_array(0, coeff, dim);

  if (!neutral || !coeff || !poly) {
    free(neutral);
    free(coeff);
    poly_destroy(poly);
    return NULL;
  }

  neutral->GF = GF;
  neutral->poly = poly;

  return neutral;
}

GF_elem_t *GF_elem_get_unity(GF_t *GF) {
  /* Get neutral and set the least significant digit to one. */
  GF_elem_t *unity = GF_elem_get_neutral(GF);

  if (!unity || !GF || (GF->n < 1)) {
    GF_elem_destroy(unity);
    return NULL;
  }

  *unity->poly->coeff = 1;

  return unity;
}

int GF_elem_get_complement(GF_elem_t *res, GF_elem_t a) {
  if (!res) {
    return 1;
  }
  assert(res->poly->len == a.poly->len);
  for (int8_t i = 0; i < a.poly->len; ++i) {
    res->poly->coeff[i] = get_complement_mod_p(a.poly->coeff[i], a.GF->p);
  }
  return 0;
}
