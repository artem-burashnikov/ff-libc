#include "poly.h"

#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include "utils.h"

poly_t *poly_from_array(int8_t degree, int8_t *coeff, int8_t len) {
  poly_t *poly = malloc(sizeof(*poly));

  if (!poly || !coeff || !len || (len < 1) || (degree < 0) || (degree >= len)) {
    poly_destroy(poly);
    return NULL;
  }

  poly->deg = degree;
  poly->coeff = coeff;
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

poly_t *poly_create_zero(int8_t len) {
  if (!len) {
    return NULL;
  }

  int8_t *res_coeff = calloc(len, sizeof(*res_coeff));
  if (!res_coeff) {
    return NULL;
  }

  poly_t *res = poly_from_array(0, res_coeff, len);

  if (!res) {
    free(res_coeff);
    return NULL;
  }

  return res;
}

void poly_normalize_deg(poly_t *a) {
  if (!a || a->deg >= a->len) {
    return;
  }

  int8_t deg = a->deg;
  while ((deg > 0) && a->coeff[deg] == 0) {
    --deg;
  }
  a->deg = deg;
  return;
}

void poly_normalize_coeff(int8_t p, poly_t *a) {
  if (!a) {
    return;
  }

  for (int8_t i = 0; i < a->len; ++i) {
    a->coeff[i] = eu_mod(a->coeff[i], p);
  }

  return;
}
