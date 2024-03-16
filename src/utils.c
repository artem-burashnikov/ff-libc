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

int8_t get_complement_mod_p(int8_t a, int8_t p) {
  assert(a > 0);
  assert(p > 1);
  assert(a < p);
  return (p - a) % p;
}
