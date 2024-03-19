#include "utils.h"

#include <assert.h>
#include <stdint.h>

#include "GF.h"
#include "poly.h"

int8_t iabs(int8_t x) { return (x < 0) ? -x : x; }

int8_t eu_mod(int8_t x, int8_t y) {
  assert(y != 0);
  int8_t r;
  r = x % y;
  if (r < 0) {
    r += iabs(y);
  }
  return r;
}

int8_t get_complement_mod_p(uint8_t a, uint8_t p) {
  assert(p > 1);
  assert(a < p);
  return (p - a) % p;
}

int8_t find_q_mod_p(int8_t x, int8_t y, int8_t p) {
  assert(p > 1);
  assert((y >= 0) && (x >= 0));
  int8_t q = 0;
  while (((y * q) % p) != x) {
    q += 1;
  }
  assert(q < p);
  return q;
}

uint64_t fpow(uint8_t base, uint8_t exp) {
  uint64_t res = 1;
  while (exp > 0) {
    if ((exp % 2) != 0) {
      res *= base;
      exp = exp - 1;
    }
    base *= base;
    exp = exp / 2;
  }
  return res;
}
