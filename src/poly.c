#include "poly.h"

#include <assert.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include "utils.h"

poly_t *poly_from_array(uint8_t deg, int8_t *coeff, uint8_t len) {
  if (!coeff || !len || (deg >= len)) {
    return NULL;
  }

  poly_t *poly = malloc(sizeof(*poly));
  int8_t *poly_coeff = malloc(sizeof(*poly->coeff) * len);

  if (!poly || !poly_coeff) {
    free(poly);
    free(poly_coeff);
    return NULL;
  }

  poly->deg = deg;
  poly->coeff = memcpy(poly_coeff, coeff, len);
  poly->len = len;

  return poly;
}

void poly_destroy(poly_t *poly) {
  if (poly) {
    free(poly->coeff);
    free(poly);
  }
}

int poly_eq(const poly_t *a, const poly_t *b) {
  if (!a || !b || (a->deg != b->deg)) {
    return 0;
  }

  /* Polynomials of the same degree may have different array lengths.
     Leading zero coefficients don't matter in that case.
     So we only check coefficients up to the minimal length. */
  int8_t n = MIN(a->len, b->len);

  for (size_t i = 0; i < n; ++i) {
    if (a->coeff[i] == b->coeff[i]) {
      continue;
    }
    return 0;
  }

  return 1;
}

poly_t *poly_cpy(const poly_t *a) {
  if (!a) {
    return NULL;
  }

  /* Initialize a new polynomial then copy the array of coefficients from the
     old one. */
  int8_t *res_coeff = malloc(sizeof(*res_coeff) * a->len);
  if (!res_coeff) {
    return NULL;
  }

  poly_t *res = poly_from_array(a->deg, res_coeff, a->len);
  if (!res) {
    free(res_coeff);
    return NULL;
  }

  memcpy(res->coeff, a->coeff, a->len);

  return res;
}

poly_t *poly_create_zero(uint8_t len) {
  if (!len) {
    return NULL;
  }
  // A zero polynomial of degree 0 is a 0-filled array of the given length.
  int8_t res_coeff[len];
  memset(res_coeff, 0, len);
  poly_t *res = poly_from_array(0, res_coeff, len);

  return res;
}

void poly_normalize_deg(poly_t *a) {
  if (!a || (a->deg >= a->len)) {
    return;
  }
  while ((a->deg > 0) && a->coeff[a->deg] == 0) {
    a->deg--;
  }
}

void poly_normalize_coeff(poly_t *a, int8_t p) {
  if (!a || (a->deg >= a->len) || (p < 2)) {
    return;
  }

  for (size_t i = 0; i < a->len; ++i) {
    a->coeff[i] = eu_mod(a->coeff[i], p);
  }
}

int poly_long_div(poly_t *a, poly_t b, int8_t p) {
  if (!a) {
    return 1;
  }

  // Assume deg a >= deg b.
  int8_t n = a->deg;
  int8_t m = b.deg;

  int8_t *u = a->coeff;
  int8_t *v = b.coeff;

  for (int8_t k = n - m; k >= 0; --k) {
    int8_t q = find_q_mod_p(u[k + m], v[m], p);
    for (int8_t i = m + k; i >= k; --i) {
      u[i] = eu_mod(u[i] - ((q * v[i - k]) % p), p);
    }
  }

  poly_normalize_deg(a);
  assert(a->deg == (b.deg - 1));
  return 0;
}
