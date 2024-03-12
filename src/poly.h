#pragma once

#include <stdint.h>

// Polynomial.
typedef struct {
  int8_t deg;    // Degree of polynomial.
  int8_t *coeff; // Array of coefficients.
} poly_t;

// Initialize a polynomial.
int poly_init(poly_t *poly, int8_t degree, int8_t *coeff);

// Destroy a given polynomial.
void poly_destroy(poly_t *poly);
