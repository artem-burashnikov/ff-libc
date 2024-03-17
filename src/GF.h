#pragma once

#include <stddef.h>
#include <stdint.h>

#include "poly.h"

// Galois field.
typedef struct {
  uint8_t p;  // Characteristic of the field GF(p).
  size_t n;   // Dimension of the field extension.
  poly_t *I;  // Irreducible polynomial over GF(p)[X] of degree n.
} GF_t;

// Element of the Galois field.
typedef struct {
  GF_t *GF;      // Galois field.
  poly_t *poly;  // Element of the GF(p)/(I) quotient field.
} GF_elem_t;

/* Initialize a quotient field. */
GF_t *GF_init(uint8_t p, size_t n, poly_t *I);

/* Given an array of coefficients of the specified length and GF(p)[x]/(I),
   Return an element over that field. */
GF_elem_t *GF_elem_from_array(int8_t *coeff, size_t len, GF_t *GF);

/* Destroy a given element of the Galois Field.
   The field itself is left untouched. */
void GF_elem_destroy(GF_elem_t *a);

/* Return 1 if the given fields are equal.
   Meaning they have the same characteristic and irreducible polynomials. */
int GF_eq(const GF_t *F, const GF_t *K);

/* res = a + b mod (I). Return 0 on success. */
int GF_elem_sum(GF_elem_t *res, GF_elem_t a, GF_elem_t b);

/* res = a - b mod (I). Return 0 on success. */
int GF_elem_diff(GF_elem_t *res, GF_elem_t a, GF_elem_t b);

/* res = a * b mod (I). Return 0 on success. */
int GF_elem_prod(GF_elem_t *res, GF_elem_t a, GF_elem_t b);

/* Calculate res: a = b * res mod (I). Return 0 on success. */
int GF_elem_div(GF_elem_t *res, GF_elem_t a, GF_elem_t b);

/* Calculate res = -a mod p. Return 0 on success. */
int GF_elem_get_complement(GF_elem_t *res, GF_elem_t a);

/* Calculate res: res * a = 1 mod (I). Return 0 on success. */
int GF_elem_get_inverse(GF_elem_t *res, GF_elem_t a);

/* Return neutral element of the given finite field. */
GF_elem_t *GF_elem_get_neutral(GF_t *GF);

/* Return unity element of the given finite field. */
GF_elem_t *GF_elem_get_unity(GF_t *GF);
