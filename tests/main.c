#include <stdlib.h>

#include "../src/GF.h"
#include "../src/poly.h"
#include "../src/utils.h"
#include "minunit.h"

const int8_t IGF2_2_coeff[3] = {1, 1, 1};
const poly_t IGF2_2 = {.coeff = IGF2_2_coeff, .deg = 2, .len = 3};
const GF_t GF2_2 = {.p = 2, .n = 2, .I = &IGF2_2};

const int8_t IGF2_5_coeff[6] = {1, 0, 1, 0, 0, 1};
const poly_t IGF2_5 = {.coeff = IGF2_5_coeff, .deg = 5, .len = 6};
const GF_t GF2_5 = {.p = 2, .n = 5, .I = &IGF2_5};

MU_TEST(poly_init_finit_test) {
  int8_t coeff_a[6] = {1, 0, 0, 40, 21, 105};
  uint8_t len_a = 6;
  uint8_t deg_a = 5;
  
  int8_t coeff_b[3] = {1, 0, 0};
  uint8_t len_b = 3;
  uint8_t deg_b = 2;

  poly_t *a = poly_from_array(deg_a, coeff_a, len_a);
  mu_check(a->deg == deg_a);
  mu_check(a->len == len_a);
  for (size_t i = 0; i < a->len; i++) {
    mu_check(a->coeff[i] == coeff_a[i]);
  }

  poly_t *b = poly_from_array(deg_b, coeff_b, len_b);
  mu_check(b->deg == deg_b);
  mu_check(b->len == len_b);
  for (size_t i = 0; i < b->len; i++) {
    mu_check(b->coeff[i] == coeff_b[i]);
  }

  poly_t *c = poly_create_zero(5);
  for (size_t i = 0; i < 5; i++) {
    mu_check(c->coeff[i] == 0);
  }

  poly_normalize_deg(b);
  mu_check(b->deg == 0);
  mu_check(b->len == 3);

  poly_normalize_deg(a);
  mu_check(a->deg == 5);
  mu_check(a->len == 6);

  poly_normalize_coeff(20, a);
  for (size_t i = 0; i < a->len; ++i) {
    mu_check(a->coeff[i] == coeff_a[i] % 20);
  }

  free(a);
  free(b);
  poly_destroy(c);
}

MU_TEST(poly_eq_test) {
  int8_t coeff_a[] = {1, 0, 0, 40, 21, 105};
  poly_t a = {.deg = 5, .coeff = coeff_a, .len = 6};

  int8_t coeff_b[] = {1, 0, 0, -70, 11, 32};
  poly_t b = {.deg = 5, .coeff = coeff_b, .len = 6};

  poly_t c = a;

  mu_check(poly_eq(&a, &c));
  mu_check(!poly_eq(&a, &b));
}

MU_TEST(gf_eq_test) {
  int8_t p = 2;
  int8_t n = 2;
  int8_t coeff[] = {1, 1, 1};
  poly_t *I = poly_from_array(2, coeff, 3);

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

MU_TEST(gf_get_neutral_test) {
  int8_t zero_coeff_2_2[2] = {0};
  poly_t zero_2_2 = {.coeff = zero_coeff_2_2, .deg = 0, .len = 2};

  int8_t zero_coeff_2_5[5] = {0};
  poly_t zero_2_5 = {.coeff = zero_coeff_2_5, .deg = 0, .len = 5};

  GF_elem_t *neutral_2_2 = GF_elem_get_neutral(&GF2_2);
  mu_check(neutral_2_2->GF->p == 2);
  mu_check(neutral_2_2->GF->n == 2);
  mu_check(neutral_2_2->GF->I->deg == 2);
  mu_check(neutral_2_2->GF->I->len == 3);
  mu_check(poly_eq(neutral_2_2->GF->I, &IGF2_2));
  mu_check(neutral_2_2->poly->deg == 0);
  mu_check(neutral_2_2->poly->len == 2);
  mu_check(poly_eq(neutral_2_2->poly, &zero_2_2));
  GF_elem_destroy(neutral_2_2);

  GF_elem_t *neutral_2_5 = GF_elem_get_neutral(&GF2_5);
  mu_check(neutral_2_5->GF->p == 2);
  mu_check(neutral_2_5->GF->n == 5);
  mu_check(neutral_2_5->GF->I->deg == 5);
  mu_check(neutral_2_5->GF->I->len == 6);
  mu_check(poly_eq(neutral_2_5->GF->I, &IGF2_5));
  mu_check(neutral_2_5->poly->deg == 0);
  mu_check(neutral_2_5->poly->len == 5);
  mu_check(poly_eq(neutral_2_5->poly, &zero_2_5));
  GF_elem_destroy(neutral_2_5);
}

MU_TEST(gf_get_unity_test) {
  int8_t one_coeff_2_5[5] = {1};
  poly_t one_2_5 = {.coeff = one_coeff_2_5, .deg = 0, .len = 5};

  GF_elem_t *unity_2_5 = GF_elem_get_unity(&GF2_5);
  mu_check(poly_eq(unity_2_5->GF->I, &IGF2_5));
  mu_check(unity_2_5->poly->deg == 0);
  mu_check(unity_2_5->poly->len == 5);
  mu_check(poly_eq(unity_2_5->poly, &one_2_5));
  GF_elem_destroy(unity_2_5);
}

int main() {
  MU_RUN_TEST(poly_eq_test);
  MU_RUN_TEST(gf_eq_test);
  MU_RUN_TEST(gf_get_neutral_test);
  MU_RUN_TEST(gf_get_unity_test);
  MU_RUN_TEST(poly_init_finit_test);
  MU_REPORT();
  return MU_EXIT_CODE;
}
