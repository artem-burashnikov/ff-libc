#pragma once

#include <stdint.h>

#include "poly.h"

// Galois field.
typedef struct {
  int8_t p;  // Characteristic of the field GF(p).
  int8_t n;  // Dimension of the field extension.
  poly_t *I; // Irreducible polynomial over GF(p)[X] of degree n.
} GF_t;

// Element of the Galois field.
typedef struct {
  GF_t *GF;     // Galois field.
  poly_t *poly; // Element of the GF(p)/(I) quotient field.
} GFelement_t;
