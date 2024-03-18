#include "GF.h"

#include <assert.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "poly.h"
#include "utils.h"

GF_t *GF_init(int8_t p, size_t n, poly_t *I) {
  GF_t *GF = malloc(sizeof(*GF));
  if (!GF) {
    return NULL;
  }

  // Artificial constraints.
  if ((n < 2) || (n > 100) || (p > 11) || (p < 2)) {
    free(GF);
    return NULL;
  }

  // Sanity check.
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

int GF_elem_normalize(GF_elem_t *a) {
  if (!a || !a->GF || !a->poly || (a->poly->len != a->GF->n)) {
    return 1;
  }

  // Normalize coefficients mod p and deg.
  poly_normalize_coeff(a->poly, a->GF->p);
  poly_normalize_deg(a->poly);

  /* Then normalize a over GF(p)[X]/(I) setting a = a mod I. */
  // If deg a < deg I: q = 0, r = a.
  if (a->poly->deg < a->GF->I->deg) {
    return 0;
  }

  // At this point deg a >= deg I
  poly_carryless_div(a->poly, *a->GF->I, a->GF->p);

  return 0;
}

int GF_eq(const GF_t *F, const GF_t *K) {
  // Irreducible polynomials must match.
  int8_t ret = poly_eq(F->I, K->I);

  // Characteristics of fields and dimensions of extensions must match.
  if (F->p != K->p || F->n != K->n) {
    ret = 0;
  }

  return ret;
}

GF_elem_t *GF_elem_from_array(int8_t *coeff, size_t len, GF_t *GF) {
  if (!coeff || (!len) || !GF) {
    return NULL;
  }

  if (!GF->I || !GF->I->coeff || (GF->I->deg < 2)) {
    return NULL;
  }

  GF_elem_t *a = malloc(sizeof(*a));
  poly_t *poly = poly_from_array(len - 1, coeff, len);
  if (!a || !poly) {
    free(a);
    poly_destroy(poly);
    return NULL;
  }

  poly_normalize_coeff(poly, GF->p);
  poly_normalize_deg(poly);

  if (poly->deg >= GF->I->deg) {
    poly_carryless_div(poly, *GF->I, GF->p);
  }

  // Normalize array length.
  if (len != GF->n) {
    int8_t *tmp = realloc(poly->coeff, sizeof(*tmp) * GF->n);
    if (!tmp) {
      free(a);
      poly_destroy(poly);
      return NULL;
    }

    // If we added any blocks, then need to initialize.
    if (len < GF->n) {
      memset(tmp + len, 0, sizeof(*tmp) * (GF->n - len));
    }

    poly->coeff = tmp;
    poly->len = GF->n;
  }

  a->GF = GF;
  a->poly = poly;

  return a;
}

GF_elem_t *GF_elem_get_neutral(GF_t *GF) {
  if (!GF) {
    return NULL;
  }

  GF_elem_t *neutral = malloc(sizeof(*neutral));
  int8_t coeff[GF->n];
  memset(coeff, 0, sizeof(*coeff) * GF->n);
  poly_t *poly = poly_from_array(0, coeff, GF->n);

  if (!neutral || !poly) {
    free(neutral);
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

GF_elem_t *GF_elem_get_complement(GF_elem_t a) {
  GF_elem_t *res = GF_elem_get_neutral(a.GF);
  if (!res) {
    return NULL;
  }

  // Assume coefficients are non-negative since a is a normalized element of the
  // field.
  for (size_t i = 0; i < a.GF->n; ++i) {
    res->poly->coeff[i] = get_complement_mod_p(a.poly->coeff[i], a.GF->p);
  }

  poly_normalize_deg(res->poly);

  return res;
}

GF_elem_t *GF_elem_get_inverse(GF_elem_t a) {
  if (a.poly->deg == 0 && (*a.poly->coeff == 0)) {
    return NULL;
  }
  /* GF(p)[x]/(I) \ {0} forms a multiplicative group G* of order (p^n - 1)
     So for any a in GF(p)[x]/(I): a^(|G*|) = 1. That implies a^(|G*| - 1) =
     inverse(a). */
  uint64_t mul_group_ord = fpow(a.GF->p, a.GF->n) - 2;

  GF_elem_t *res = GF_elem_get_neutral(a.GF);
  if (!res) {
    return NULL;
  }

  poly_fpowm(res->poly, *a.poly, mul_group_ord, *a.GF->I, a.GF->p);

  return res;
}
