#pragma once

#include <stdint.h>

// Polynomial.
typedef struct {
  int8_t deg;     // Degree of polynomial.
  int8_t *coeff;  // Array of coefficients.
  int8_t len;     // Length of the array of coefficients.
} poly_t;

// Initialize a polynomial.
poly_t *poly_init(int8_t degree, int8_t *coeff, int8_t len);

// Destroy a given polynomial.
void poly_destroy(poly_t *poly);

/* Check if two polynomials are equal.
   Return 1 if the degree and corresponding coefficients match. */
int poly_eq(const poly_t *a, const poly_t *b);
