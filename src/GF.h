#pragma once

#include <stdint.h>

#include "poly.h"

// Galois field.
typedef struct {
  int8_t p;   // Characteristic of the field GF(p).
  int8_t n;   // Dimension of the field extension.
  poly_t *I;  // Irreducible polynomial over GF(p)[X] of degree n.
} GF_t;

// Element of the Galois field.
typedef struct {
  GF_t *GF;      // Galois field.
  poly_t *poly;  // Element of the GF(p)/(I) quotient field.
} GFelement_t;

// Initialize a quotient field.
GF_t *GF_init(int8_t p, int8_t n, poly_t *I);

// Destroy a given Galois field.
void GF_destroy(GF_t *GF);

/* Destroy a given element of the Galois Field.
   The field itself is left untouched. */
void GFelement_destroy(GFelement_t *a);

/* Return 1 if the given fields are equal.
   Meaning they have the same characteristic and irreducible polynomials. */
int GF_eq(GF_t *f, GF_t *k);

/* Calculate the sum in res of a and b.
   Feilds must match. Return 0 on success. */
int GFelement_add(GFelement_t *res, GFelement_t *a, GFelement_t *b);

/* Calculate the product in res of a and b.
   Feilds must match. Return 0 on success. */
int GFelement_mul(GFelement_t *res, GFelement_t *a, GFelement_t *b);

/* Calculate q and r such that a = bq + r, where deg r < deg b.
   Arguments q and/or r can be NULL, meaning the result for that variable won't
   be stored. Feilds must match. Return 0 on success. */
int GFelement_divmod(GFelement_t *q, GFelement_t *r, GFelement_t *a,
                     GFelement_t *b);

/* Return neutral element in the given finite field. */
GFelement_t *GFelement_get_neutral(GF_t *GF);

/* Return unity element in the given finite field. */
GFelement_t *GFelement_get_unity(GF_t *GF);
