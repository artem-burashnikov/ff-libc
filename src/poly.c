#include "poly.h"

#include <assert.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "utils.h"

poly_t *poly_from_array(uint8_t deg, int8_t *coeff, size_t len) {
  if (!coeff || (!len) || (deg >= len)) {
    return NULL;
  }

  poly_t *poly = malloc(sizeof(*poly));
  int8_t *tmp = malloc(sizeof(*poly->coeff) * len);

  if (!poly || !tmp) {
    free(poly);
    free(tmp);
    return NULL;
  }

  poly->deg = deg;
  poly->coeff = memcpy(tmp, coeff, len * sizeof(*tmp));
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

  /* At this point polynomials have the same degree.
     Polynomials of the same degree may have different array lengths.
     Leading zero coefficients don't matter in that case.
     So we only check coefficients up to the degree.
     Not using memcmp, since coefficeints may be negative. */
  assert(a->deg == b->deg);
  for (size_t i = 0; i <= a->deg; ++i) {
    if (a->coeff[i] != b->coeff[i]) {
      return 0;
    }
  }

  return 1;
}

poly_t *poly_create_zero(size_t len) {
  if (!len) {
    return NULL;
  }
  // A zero polynomial of degree 0 is a 0-filled array of the given length.
  int8_t *tmp = calloc(len, sizeof(*tmp));

  poly_t *res = poly_from_array(0, tmp, len);

  free(tmp);

  return res;
}

void poly_normalize_deg(poly_t *a) {
  if (!a || (a->deg >= a->len)) {
    return;
  }
  size_t k = a->len - 1;
  while ((k > 0) && a->coeff[k] == 0) {
    k--;
  }
  a->deg = k;
}

void poly_normalize_coeff(poly_t *a, int8_t p) {
  if (!a || (a->deg >= a->len) || (p < 2)) {
    return;
  }

  for (size_t i = 0; i < a->len; ++i) {
    a->coeff[i] = eu_mod(a->coeff[i], p);
  }
}

/* Set res = a + b, where a and b are polynomials over Fp.
   Assume max(a.deg, b.deg) < min(a.len, b.len) */
void poly_carryless_sum(poly_t *res, poly_t a, poly_t b, int8_t p) {
  if (!res) {
    return;
  }

  int8_t w;
  size_t max_deg = MAX(a.deg, b.deg);
  for (size_t i = 0; i <= max_deg; ++i) {
    w = 0;
    if (i <= a.deg) {
      w += a.coeff[i];
    }
    if (i <= b.deg) {
      w += b.coeff[i];
    }
    res->coeff[i] = eu_mod(w, p);
  }
}

/* Set a = a mod b, where a and b are polynomials over Fp. */
void poly_carryless_div(poly_t *a, poly_t b, int8_t p) {
  if (!a) {
    return;
  }

  // Assume deg a >= deg b.
  int8_t n = a->deg;
  int8_t m = b.deg;

  int8_t *u = a->coeff;
  int8_t *v = b.coeff;

  int8_t q;
  for (int8_t k = n - m; k >= 0; --k) {
    q = find_q_mod_p(u[k + m], v[m], p);
    for (int8_t i = m + k; i >= k; --i) {
      u[i] = eu_mod(u[i] - (q * v[i - k]), p);
    }
  }
  poly_normalize_deg(a);
  assert(a->deg < b.deg);
}

// Set res = a * b mod p. Must be guaranteed res is niether a nor b.
void poly_carryless_mul(poly_t *res, poly_t a, poly_t b, int8_t p) {
  if (!res || (res->len < (a.deg + b.deg + 1))) {
    return;
  }

  memset(res->coeff, 0, res->len * sizeof(*res->coeff));

  for (size_t i = 0; i <= a.deg; ++i) {
    for (size_t j = 0; j <= b.deg; ++j) {
      if (a.coeff[i] == 0) {
        break;
      }
      res->coeff[i + j] =
          eu_mod(res->coeff[i + j] + a.coeff[i] * b.coeff[j], p);
    }
  }
}

void poly_fpowm(poly_t *res, poly_t a, uint64_t exp, poly_t I, int8_t p) {
  if (!res) {
    return;
  }

  int8_t *tmp = NULL;
  poly_t *base = poly_create_zero(I.deg + I.deg);
  base->deg = a.deg;
  memcpy(base->coeff, a.coeff, (a.deg + 1) * sizeof(*base->coeff));

  // Temporary buffer
  poly_t *buff = poly_create_zero(I.deg + I.deg);

  // Set prod equal to 1. Prod holds the result.
  poly_t *prod = poly_create_zero(I.deg + I.deg);
  *prod->coeff = 1;

  while (exp > 0) {
    if ((exp % 2) != 0) {
      // Set buff = prod * base
      poly_carryless_mul(buff, *prod, *base, p);
      poly_normalize_deg(buff);
      poly_carryless_div(buff, I, p);
      exp = exp - 1;
      // Swap buff and prod.
      tmp = prod->coeff;
      prod->coeff = buff->coeff;
      prod->deg = buff->deg;
      prod->len = buff->len;
      buff->coeff = tmp;
    }
    // Set buff = base * base;
    poly_carryless_mul(buff, *base, *base, p);
    poly_normalize_deg(buff);
    poly_carryless_div(buff, I, p);
    exp = exp / 2;
    // Swap buff and base.
    tmp = base->coeff;
    base->coeff = buff->coeff;
    base->deg = buff->deg;
    base->len = buff->len;
    buff->coeff = tmp;
  }

  memcpy(res->coeff, prod->coeff, sizeof(*prod->coeff) * (prod->deg + 1));
  res->deg = prod->deg;

  // Clean up.
  poly_destroy(prod);
  poly_destroy(buff);
  poly_destroy(base);
}
