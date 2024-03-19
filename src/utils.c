#include "utils.h"

#include <assert.h>
#include <stdint.h>
#include <stdlib.h>

#include "GF.h"
#include "poly.h"

uint8_t complement(uint8_t a, uint8_t p) {
  assert(a < p);
  return (p - a) % p;
}

// Asuume p is prime.
int8_t inverse(int8_t a, int8_t p) {
  assert(a != 0);
  int8_t t = 0;
  int8_t new_t = 1;
  int8_t r = p;
  int8_t new_r = a;

  int8_t q;
  int8_t tmp;
  while (new_r != 0) {
    q = r / new_r;

    tmp = new_t;
    new_t = t - q * new_t;
    t = tmp;

    tmp = new_r;
    new_r = r - q * new_r;
    r = tmp;
  }
  t = (t < 0) ? (t + p) : t;
  return t;
}

// Assume p is prime.
uint8_t x_div_y_mod_p(uint8_t x, uint8_t y, uint8_t p) {
  uint8_t q = 0;
  while (((y * q) % p) != x) {
    q += 1;
  }
  return q;
}

uint64_t fpow(uint8_t base, uint8_t exp) {
  uint64_t res = 1;
  while (exp > 0) {
    if ((exp % 2) != 0) {
      res *= base;
    }
    base *= base;
    exp = exp / 2;
  }
  return res;
}
