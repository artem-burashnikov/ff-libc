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
} GF_elem_t;

// Initialize a quotient field.
GF_t *GF_init(int8_t p, int8_t n, poly_t *I);

/* Destroy a given element of the Galois Field.
   The field itself is left untouched. */
void GF_elem_destroy(GF_elem_t *a);

/* Set coefficients mod p. */
void GF_elem_normalize(GF_elem_t *a);

/* Return 1 if the given fields are equal.
   Meaning they have the same characteristic and irreducible polynomials. */
int GF_eq(const GF_t *f, const GF_t *k);

/* res = a + b. Return 0 on success. */
int GF_elem_sum(GF_elem_t *res, GF_elem_t a, GF_elem_t b);

/* res = a - b. Return 0 on success. */
int GF_elem_diff(GF_elem_t *res, GF_elem_t a, GF_elem_t b);

/* res = a * b in GF(p)/(I). Return 0 on success. */
int GF_elem_prod(GF_elem_t *res, GF_elem_t a, GF_elem_t b);

/* Find q and r over GF(p)[x]/(I) such that a = bq + r, where deg r < deg b.
   One of q or r can be NULL, meaning the result for that variable won't
   be stored. Return 0 on success. */
int GF_elem_div(GF_elem_t *q, GF_elem_t *r, GF_elem_t a, GF_elem_t b);

/* Return neutral element of the given finite field. */
GF_elem_t *GF_elem_get_neutral(GF_t *GF);

/* Return unity element of the given finite field. */
GF_elem_t *GF_elem_get_unity(GF_t *GF);

/* Calculate res = -a mod p. Return 0 on success. */
int GF_elem_get_complement(GF_elem_t *res, GF_elem_t a);

/* Calculate res s.t. res * a = 1 mod (I) */
int GF_elem_get_inverse(GF_elem_t *res, GF_elem_t a);
