#include <stdlib.h>

#include "../src/GF.h"
#include "../src/poly.h"
#include "../src/utils.h"
#include "minunit.h"

// x^2 + x + 1
int8_t IGF2_2_coeff[3] = {1, 1, 1};
poly_t IGF2_2 = {.coeff = IGF2_2_coeff, .deg = 2, .len = 3};
GF_t GF2_2 = {.p = 2, .n = 2, .I = &IGF2_2};

// x^5 + x^2 + 1
int8_t IGF2_5_coeff[6] = {1, 0, 1, 0, 0, 1};
poly_t IGF2_5 = {.coeff = IGF2_5_coeff, .deg = 5, .len = 6};
GF_t GF2_5 = {.p = 2, .n = 5, .I = &IGF2_5};

// x^4 + 5x^2 + 4x + 3
int8_t IGF7_4_coeff[5] = {3, 4, 5, 0, 1};
poly_t IGF7_4 = {.coeff = IGF7_4_coeff, .deg = 4, .len = 5};
GF_t GF7_4 = {.p = 7, .n = 4, .I = &IGF7_4};

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

  poly_normalize_coeff(a, 20);
  for (size_t i = 0; i < a->len; ++i) {
    mu_check(a->coeff[i] == coeff_a[i] % 20);
  }

  poly_destroy(a);
  poly_destroy(b);
  poly_destroy(c);
}

MU_TEST(poly_eq_test) {
  int8_t coeff_a[] = {1, 0, 0, 40, 21, 105};
  poly_t *a = poly_from_array(5, coeff_a, 6);

  int8_t coeff_b[] = {1, 0, 0, -70, 11, 32};
  poly_t *b = poly_from_array(5, coeff_b, 6);

  mu_check(poly_eq(a, a));
  mu_check(!poly_eq(a, b));

  poly_destroy(a);
  poly_destroy(b);
}

MU_TEST(poly_div_test_case1) {
  /* Dimension. */
  uint8_t n = 10;

  /* Characteristic. */
  uint8_t p = 10;

  // 1 + x
  int8_t arr1[10] = {1, 1, 0, 0, 0, 0, 0, 0, 0, 0};
  poly_t *a = poly_from_array(1, arr1, n);

  // x^2 + x^4 + x^8
  int8_t arr2[10] = {0, 0, 1, 0, 1, 0, 0, 0, 1, 0};
  poly_t *b = poly_from_array(8, arr2, n);

  poly_carryless_div(b, *a, p);
  mu_check(b->deg == 0);
  mu_check(*b->coeff == 3);

  poly_destroy(a);
  poly_destroy(b);
}

MU_TEST(poly_div_test_case2) {
  /* Dimension. */
  uint8_t n = 10;

  /* Characteristic. */
  uint8_t p = 5;

  // 2 + 3x^2
  int8_t arr1[10] = {2, 0, 3, 0, 0, 0, 0, 0, 0, 0};
  poly_t *a = poly_from_array(2, arr1, n);

  // 1 + x + 4x^2
  int8_t arr2[10] = {1, 1, 4, 0, 0, 0, 0, 0, 0, 0};
  poly_t *b = poly_from_array(2, arr2, n);

  poly_carryless_div(a, *b, p);
  mu_check(a->deg == 1);
  mu_check(a->coeff[0] == 0);
  mu_check(a->coeff[1] == 3);

  poly_destroy(a);
  poly_destroy(b);
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

  poly_destroy(I);
  free(GF2_a);
  free(GF2_b);
  free(GF2_c);
  free(GF2_d);
}

MU_TEST(gf_elem_from_array_test1) {
  int8_t coeff[] = {5, 11, 20, 0, 20};
  GF_elem_t *a = GF_elem_from_array(coeff, 5, &GF2_5);

  mu_check(GF_eq(a->GF, &GF2_5));
  mu_check(a->poly->len == 5);
  mu_check(a->poly->deg == 1);
  for (size_t i = 0; i < 5; ++i) {
    mu_check(a->poly->coeff[i] == coeff[i] % 2);
  }

  GF_elem_destroy(a);
}

MU_TEST(gf_elem_from_array_test2) {
  int8_t coeff[] = {1, 1, 0, 0, 0, 0, 0, 1};
  GF_elem_t *a = GF_elem_from_array(coeff, 8, &GF2_5);

  mu_check(GF_eq(a->GF, &GF2_5));
  mu_check(a->poly->len == 5);
  mu_check(a->poly->deg == 4);
  mu_check(a->poly->coeff[0] == 1);
  mu_check(a->poly->coeff[1] == 1);
  mu_check(a->poly->coeff[2] == 1);
  mu_check(a->poly->coeff[3] == 0);
  mu_check(a->poly->coeff[4] == 1);

  GF_elem_destroy(a);
}

MU_TEST(gf_elem_from_array_test3) {
  int8_t coeff[] = {1, 1};
  GF_elem_t *a = GF_elem_from_array(coeff, 2, &GF2_5);

  mu_check(GF_eq(a->GF, &GF2_5));
  mu_check(a->poly->len == 5);
  mu_check(a->poly->deg == 1);
  mu_check(a->poly->coeff[0] == 1);
  mu_check(a->poly->coeff[1] == 1);
  mu_check(a->poly->coeff[2] == 0);
  mu_check(a->poly->coeff[3] == 0);
  mu_check(a->poly->coeff[4] == 0);

  GF_elem_destroy(a);
}

MU_TEST(gf_elem_from_array_test4) {
  int8_t coeff[] = {1, 1, 0, 1, 0};
  GF_elem_t *a = GF_elem_from_array(coeff, 5, &GF2_5);

  mu_check(GF_eq(a->GF, &GF2_5));
  mu_check(a->poly->len == 5);
  mu_check(a->poly->deg == 3);
  mu_check(a->poly->coeff[0] == 1);
  mu_check(a->poly->coeff[1] == 1);
  mu_check(a->poly->coeff[2] == 0);
  mu_check(a->poly->coeff[3] == 1);
  mu_check(a->poly->coeff[4] == 0);

  GF_elem_destroy(a);
}

MU_TEST(gf_get_neutral_test) {
  poly_t *zero_2_2 = poly_create_zero(2);
  poly_t *zero_2_5 = poly_create_zero(5);

  GF_elem_t *neutral_2_2 = GF_elem_get_neutral(&GF2_2);
  mu_check(neutral_2_2->GF->p == 2);
  mu_check(neutral_2_2->GF->n == 2);
  mu_check(neutral_2_2->GF->I->deg == 2);
  mu_check(neutral_2_2->GF->I->len == 3);
  mu_check(poly_eq(neutral_2_2->GF->I, &IGF2_2));
  mu_check(neutral_2_2->poly->deg == 0);
  mu_check(neutral_2_2->poly->len == 2);
  mu_check(poly_eq(neutral_2_2->poly, zero_2_2));

  GF_elem_destroy(neutral_2_2);
  poly_destroy(zero_2_2);

  GF_elem_t *neutral_2_5 = GF_elem_get_neutral(&GF2_5);
  mu_check(neutral_2_5->GF->p == 2);
  mu_check(neutral_2_5->GF->n == 5);
  mu_check(neutral_2_5->GF->I->deg == 5);
  mu_check(neutral_2_5->GF->I->len == 6);
  mu_check(poly_eq(neutral_2_5->GF->I, &IGF2_5));
  mu_check(neutral_2_5->poly->deg == 0);
  mu_check(neutral_2_5->poly->len == 5);
  mu_check(poly_eq(neutral_2_5->poly, zero_2_5));

  GF_elem_destroy(neutral_2_5);
  poly_destroy(zero_2_5);
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

MU_TEST(gf_sum_test1) {
  int8_t coeff_a[] = {1, 1, 1, 0, 0};
  GF_elem_t *a = GF_elem_from_array(coeff_a, 5, &GF2_5);
  mu_check(a->poly->deg == 2);

  int8_t coeff_b[] = {0, 1, 1, 0, 1};
  GF_elem_t *b = GF_elem_from_array(coeff_b, 5, &GF2_5);
  mu_check(b->poly->deg == 4);

  GF_elem_sum(a, *a, *b);
  mu_check(a->poly->len == 5);
  mu_check(a->poly->deg == 4);
  mu_check(a->poly->coeff[0] == 1);
  mu_check(a->poly->coeff[1] == 0);
  mu_check(a->poly->coeff[2] == 0);
  mu_check(a->poly->coeff[3] == 0);
  mu_check(a->poly->coeff[4] == 1);

  GF_elem_destroy(a);
  GF_elem_destroy(b);
}

MU_TEST(gf_sum_test2) {
  int8_t coeff_a[] = {1, 1, 1, 1, 1};
  GF_elem_t *a = GF_elem_from_array(coeff_a, 5, &GF2_5);

  int8_t coeff_b[] = {1, 1, 1, 1, 1};
  GF_elem_t *b = GF_elem_from_array(coeff_b, 5, &GF2_5);

  GF_elem_sum(a, *a, *b);
  mu_check(a->poly->len == 5);
  mu_check(a->poly->deg == 0);
  mu_check(a->poly->coeff[0] == 0);
  mu_check(a->poly->coeff[1] == 0);
  mu_check(a->poly->coeff[2] == 0);
  mu_check(a->poly->coeff[3] == 0);
  mu_check(a->poly->coeff[4] == 0);

  GF_elem_destroy(a);
  GF_elem_destroy(b);
}

MU_TEST(gf_diff_test1) {
  int8_t coeff_a[5] = {1, 1, 1, 1, 1};
  GF_elem_t *a = GF_elem_from_array(coeff_a, 5, &GF2_5);

  int8_t coeff_b[5] = {1, 1, 1, 1, 1};
  GF_elem_t *b = GF_elem_from_array(coeff_b, 5, &GF2_5);

  GF_elem_diff(a, *a, *b);
  mu_check(a->poly->len == 5);
  mu_check(a->poly->deg == 0);
  mu_check(a->poly->coeff[0] == 0);
  mu_check(a->poly->coeff[1] == 0);
  mu_check(a->poly->coeff[2] == 0);
  mu_check(a->poly->coeff[3] == 0);
  mu_check(a->poly->coeff[4] == 0);

  GF_elem_destroy(a);
  GF_elem_destroy(b);
}

MU_TEST(gf_diff_test2) {
  int8_t coeff_a[] = {0, 5, 3, 1};
  GF_elem_t *a = GF_elem_from_array(coeff_a, 4, &GF7_4);

  int8_t coeff_b[] = {3, 0, 6, 0};
  GF_elem_t *b = GF_elem_from_array(coeff_b, 4, &GF7_4);

  GF_elem_diff(a, *a, *b);
  mu_check(a->poly->len == 4);
  mu_check(a->poly->deg == 3);
  mu_check(a->poly->coeff[0] == 4);
  mu_check(a->poly->coeff[1] == 5);
  mu_check(a->poly->coeff[2] == 4);
  mu_check(a->poly->coeff[3] == 1);

  GF_elem_destroy(a);
  GF_elem_destroy(b);
}

MU_TEST(poly_mul_no_carry_test1) {
  int8_t coeff[] = {1, 1};
  poly_t *a = poly_from_array(1, coeff, 2);
  poly_t *res = poly_create_zero(3);
  poly_carryless_mul(res, *a, *a, 3);

  mu_check(res->deg = 2);
  mu_check(res->len = 3);
  mu_check(res->coeff[0] == 1);
  mu_check(res->coeff[1] == 2);
  mu_check(res->coeff[2] == 1);

  poly_destroy(res);
  poly_destroy(a);
}

MU_TEST(poly_mul_no_carry_test2) {
  int8_t coeff_a[] = {1, 0, 1, 5, 15};
  poly_t *a = poly_from_array(4, coeff_a, 5);

  int8_t coeff_b[] = {1, 0, 1, 0, 0};
  poly_t *b = poly_from_array(2, coeff_b, 5);

  poly_t *res = poly_create_zero(7);
  poly_carryless_mul(res, *a, *b, 17);

  mu_check(res->deg = 6);
  mu_check(res->len = 7);
  mu_check(res->coeff[0] == 1);
  mu_check(res->coeff[1] == 0);
  mu_check(res->coeff[2] == 2);
  mu_check(res->coeff[3] == 5);
  mu_check(res->coeff[4] == 16);
  mu_check(res->coeff[5] == 5);
  mu_check(res->coeff[6] == 15);

  poly_destroy(res);
  poly_destroy(a);
  poly_destroy(b);
}

MU_TEST(gf_prod_test1) {
  int8_t coeff_a[] = {1, 0, 1, 5, 15};
  GF_elem_t *a = GF_elem_from_array(coeff_a, 5, &GF2_5);

  int8_t coeff_b[] = {1, 0, 1, 0};
  GF_elem_t *b = GF_elem_from_array(coeff_b, 4, &GF2_5);

  GF_elem_t *res = GF_elem_get_neutral(&GF2_5);

  // 15x^6 + 5x^5 + 16x^4 + 5x^3 + 2x^2 + 1
  GF_elem_prod(res, *a, *b);

  mu_check(res->poly->deg = 2);
  mu_check(res->poly->len = 5);
  mu_check(res->poly->coeff[0] == 0);
  mu_check(res->poly->coeff[1] == 1);
  mu_check(res->poly->coeff[2] == 1);
  mu_check(res->poly->coeff[3] == 0);
  mu_check(res->poly->coeff[4] == 0);

  GF_elem_destroy(a);
  GF_elem_destroy(b);
  GF_elem_destroy(res);
}

MU_TEST(gf_prod_test2) {
  int8_t coeff[] = {1, 1};
  GF_elem_t *a = GF_elem_from_array(coeff, 2, &GF2_5);
  GF_elem_t *res = GF_elem_get_neutral(&GF2_5);

  GF_elem_prod(res, *a, *a);

  mu_check(res->poly->deg = 2);
  mu_check(res->poly->len = 5);
  mu_check(res->poly->coeff[0] == 1);
  mu_check(res->poly->coeff[1] == 0);
  mu_check(res->poly->coeff[2] == 1);
  mu_check(res->poly->coeff[3] == 0);
  mu_check(res->poly->coeff[4] == 0);

  GF_elem_destroy(a);
  GF_elem_destroy(res);
}

MU_TEST(poly_fpowm_test1) {
  int8_t coeff[2] = {1, 1, 0};
  poly_t *a = poly_from_array(1, coeff, 3);
  poly_t *res = poly_create_zero(3);
  poly_fpowm(res, *a, 2, IGF2_2, 2);

  mu_check(res->deg = 1);
  mu_check(res->len = 3);
  mu_check(res->coeff[0] == 0);
  mu_check(res->coeff[1] == 1);

  poly_destroy(res);
  poly_destroy(a);
}

MU_TEST_SUITE(gf_tests) {
  MU_RUN_TEST(gf_eq_test);
  MU_RUN_TEST(gf_get_neutral_test);
  MU_RUN_TEST(gf_get_unity_test);
  MU_RUN_TEST(gf_elem_from_array_test1);
  MU_RUN_TEST(gf_elem_from_array_test2);
  MU_RUN_TEST(gf_elem_from_array_test3);
  MU_RUN_TEST(gf_elem_from_array_test4);
}

MU_TEST_SUITE(gf_arithmetic_tests) {
  MU_RUN_TEST(poly_div_test_case1);
  MU_RUN_TEST(poly_div_test_case2);
  MU_RUN_TEST(poly_mul_no_carry_test1);
  MU_RUN_TEST(poly_mul_no_carry_test2);
  MU_RUN_TEST(gf_sum_test1);
  MU_RUN_TEST(gf_sum_test2);
  MU_RUN_TEST(gf_diff_test1);
  MU_RUN_TEST(gf_diff_test2);
  MU_RUN_TEST(gf_prod_test1);
  MU_RUN_TEST(gf_prod_test2);
  MU_RUN_TEST(poly_fpowm_test1);
}

int main() {
  MU_RUN_TEST(poly_eq_test);
  MU_RUN_TEST(poly_init_finit_test);
  MU_RUN_SUITE(gf_tests);
  MU_RUN_SUITE(gf_arithmetic_tests);
  MU_REPORT();
  return MU_EXIT_CODE;
}
