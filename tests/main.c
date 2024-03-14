#include <stdlib.h>

#include "gf.h"
#include "minunit.h"
#include "poly.h"

MU_TEST(poly_eq_test) {
  int8_t coeff_a[] = {1, 0, 0, 40, 21, 105};
  poly_t a = {.deg = 5, .coeff = coeff_a};

  int8_t coeff_b[] = {1, 0, 0, -70, 11, 32};
  poly_t b = {.deg = 5, .coeff = coeff_b};

  poly_t c = a;

  mu_check(poly_eq(&a, &c));
  mu_check(!poly_eq(&a, &b));
}

MU_TEST(gf_eq_test) {
  int8_t p = 2;
  int8_t n = 2;
  int8_t coeff[] = {1, 1, 1};
  poly_t *I = poly_init(2, coeff);

  GF_t *GF2_a = GF_init(p, n, I);
  GF_t *GF2_b = GF_init(p, n, I);
  GF_t *GF2_c = GF_init(p + 1, n, I);
  GF_t *GF2_d = GF_init(p, n + 1, I);

  mu_check(GF_eq(GF2_a, GF2_b));
  mu_check(!GF_eq(GF2_a, GF2_c));
  mu_check(!GF_eq(GF2_c, GF2_d));

  free(I);
  free(GF2_a);
  free(GF2_b);
  free(GF2_c);
  free(GF2_d);
}

int main() {
  MU_RUN_TEST(poly_eq_test);
  MU_RUN_TEST(gf_eq_test);
  MU_REPORT();
  return MU_EXIT_CODE;
}
