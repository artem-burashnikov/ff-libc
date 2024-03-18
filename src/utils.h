#pragma once

#include <stdint.h>

#define MAX(A, B) ((A) > (B) ? (A) : (B))

// Return integer abs.
int8_t iabs(int8_t x);

// Return euclidean x mod y.
int8_t eu_mod(int8_t x, int8_t y);

// Return complement(a) : a + complement(a) = 0 mod p.
int8_t get_complement_mod_p(uint8_t a, uint8_t p);

// Find q such that y*q = x mod p.
int8_t find_q_mod_p(int8_t x, int8_t y, int8_t p);

#if 0
// Find inverse of d mod p: d * inverse(d) = 1 mod p.
int8_t get_inv_mod_p(int8_t d, int8_t p);
#endif

// Fast base to the exponent.
uint8_t fpow(uint8_t base, uint8_t exp);
