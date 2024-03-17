#pragma once

#include <stdint.h>

// Polynomial.
typedef struct {
  uint8_t deg;    // Degree of polynomial.
  int8_t *coeff;  // Array of coefficients.
  uint8_t len;    // Length of the array of coefficients.
} poly_t;

// Initialize a polynomial.
poly_t *poly_from_array(uint8_t deg, int8_t *coeff, uint8_t len);

// Destroy a given polynomial.
void poly_destroy(poly_t *a);

/* Check if two polynomials are equal.
   Meaning they have the same degree and matching coefficients.
   Return 1 if the degree and corresponding coefficients match. */
int poly_eq(const poly_t *a, const poly_t *b);

/* Return a copy of the given polynomial. */
poly_t *poly_cpy(const poly_t *a);

/* Return a zero polynomial of the given length.*/
poly_t *poly_create_zero(uint8_t len);

/* Set a = a mod b, where a and b are polynomials over Fp. */
int poly_long_div(poly_t *a, poly_t b, int8_t p);

/* Normalize the degree (find the greatest non zero coefficient index)
   of the given polynomial. */
void poly_normalize_deg(poly_t *a);

/* Normalize coefficients modulo p of the given polynomial. */
void poly_normalize_coeff(poly_t *a, int8_t p);
