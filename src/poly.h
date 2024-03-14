#pragma once

#include <stdint.h>

#define MINUS_INFINITY -1

// Polynomial.
typedef struct {
  int8_t deg;     // Degree of polynomial.
  int8_t *coeff;  // Array of coefficients.
} poly_t;

/* Initialize a polynomial.
   No checking is done to ensure that degree matches the array of coefficients.
 */
poly_t *poly_init(int8_t degree, int8_t *coeff);

// Destroy a given polynomial.
void poly_destroy(poly_t *poly);

/* Check if two polynomials are equal.
   Return 1 if the degree and corresponding coefficients match. */
int poly_eq(poly_t *a, poly_t *b);
