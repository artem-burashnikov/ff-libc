#include <stdlib.h>

#include "GF.h"
#include "minunit.h"
#include "poly.h"
#include "utils.h"

// x^2 + x + 1
uint8_t IGF2_2_coeff[3] = {1, 1, 1};
poly_t IGF2_2 = {.coeff = IGF2_2_coeff, .deg = 2};
GF_t GF2_2 = {.p = 2, .I = &IGF2_2};

// x^3 + x + 1
uint8_t IGF2_3_coeff[4] = {1, 1, 0, 1};
poly_t IGF2_3 = {.coeff = IGF2_3_coeff, .deg = 3};
GF_t GF2_3 = {.p = 2, .I = &IGF2_3};

// x^5 + x^2 + 1
uint8_t IGF2_5_coeff[6] = {1, 0, 1, 0, 0, 1};
poly_t IGF2_5 = {.coeff = IGF2_5_coeff, .deg = 5};
GF_t GF2_5 = {.p = 2, .I = &IGF2_5};

// x^4 + 5x^2 + 4x + 3
uint8_t IGF7_4_coeff[5] = {3, 4, 5, 0, 1};
poly_t IGF7_4 = {.coeff = IGF7_4_coeff, .deg = 4};
GF_t GF7_4 = {.p = 7, .I = &IGF7_4};

MU_TEST(poly_init_finit_test) {
  uint8_t coeff_a[6] = {1, 0, 0, 40, 21, 105};
  uint8_t deg_a = 5;

  uint8_t coeff_b[3] = {1, 0, 0};
  uint8_t deg_b = 0;

  poly_t *a = poly_from_array(deg_a, coeff_a);
  mu_check(a->deg == deg_a);
  for (size_t i = 0; i <= a->deg; i++) {
    mu_check(a->coeff[i] == coeff_a[i]);
  }

  poly_t *b = poly_from_array(deg_b, coeff_b);
  mu_check(b->deg == deg_b);
  for (size_t i = 0; i <= b->deg; i++) {
    mu_check(b->coeff[i] == coeff_b[i]);
  }

  poly_t *c = poly_create_zero(5);
  mu_check(c->deg == 0);
  for (size_t i = 0; i < 5; i++) {
    mu_check(c->coeff[i] == 0);
  }

  poly_destroy(a);
  poly_destroy(b);
  poly_destroy(c);
}

MU_TEST(poly_eq_test) {
  uint8_t coeff_a[6] = {1, 0, 0, 40, 21, 105};
  poly_t *a = poly_from_array(5, coeff_a);

  uint8_t coeff_b[6] = {1, 0, 0, 70, 11, 32};
  poly_t *b = poly_from_array(5, coeff_b);

  mu_check(poly_eq(a, a));
  mu_check(!poly_eq(a, b));

  poly_destroy(a);
  poly_destroy(b);
}

MU_TEST(poly_div_test_case1) {
  // 1 + x
  uint8_t arr1[2] = {1, 1};
  poly_t *a = poly_from_array(1, arr1);

  // x^2 + x^4 + x^8
  uint8_t arr2[9] = {0, 0, 1, 0, 1, 0, 0, 0, 1};
  poly_t *b = poly_from_array(8, arr2);

  poly_t *res = poly_create_zero(9);

  poly_div(res, *b, *a, 10);
  mu_check(res->deg == 0);
  mu_check(*res->coeff == 3);

  poly_destroy(a);
  poly_destroy(b);
  poly_destroy(res);
}

MU_TEST(poly_div_test_case2) {
  // 2 + 3x^2
  uint8_t arr1[10] = {2, 0, 3, 0, 0, 0, 0, 0, 0, 0};
  poly_t *a = poly_from_array(2, arr1);

  // 1 + x + 4x^2
  uint8_t arr2[10] = {1, 1, 4, 0, 0, 0, 0, 0, 0, 0};
  poly_t *b = poly_from_array(2, arr2);

  poly_div(a, *a, *b, 5);
  mu_check(a->deg == 1);
  mu_check(a->coeff[0] == 0);
  mu_check(a->coeff[1] == 3);

  poly_destroy(a);
  poly_destroy(b);
}

MU_TEST(gf_eq_test) {
  uint8_t coeff[] = {1, 1, 1};
  poly_t *I = poly_from_array(2, coeff);

  GF_t *GF2_a = GF_init_field(2, *I);
  GF_t *GF2_b = GF_init_field(2, *I);
  GF_t *GF2_c = GF_init_field(2, IGF2_5);

  mu_check(GF_eq(GF2_a, GF2_b));
  mu_check(!GF_eq(GF2_a, GF2_c));

  poly_destroy(I);
  GF_destroy_field(GF2_a);
  GF_destroy_field(GF2_b);
  GF_destroy_field(GF2_c);
}

MU_TEST(gf_elem_from_array_test1) {
  uint8_t coeff[] = {5, 11, 20, 0, 20};
  GF_elem_t *a = GF_elem_from_array(4, coeff, &GF2_5);

  mu_check(GF_eq(a->GF, &GF2_5));
  mu_check(a->poly->deg == 1);
  for (size_t i = 0; i < 5; ++i) {
    mu_check(a->poly->coeff[i] == coeff[i] % 2);
  }

  GF_elem_destroy(a);
}

MU_TEST(gf_elem_from_array_test2) {
  uint8_t coeff[] = {1, 1, 0, 0, 0, 0, 0, 1};
  GF_elem_t *a = GF_elem_from_array(7, coeff, &GF2_5);

  mu_check(GF_eq(a->GF, &GF2_5));
  mu_check(a->poly->deg == 4);
  mu_check(a->poly->coeff[0] == 1);
  mu_check(a->poly->coeff[1] == 1);
  mu_check(a->poly->coeff[2] == 1);
  mu_check(a->poly->coeff[3] == 0);
  mu_check(a->poly->coeff[4] == 1);

  GF_elem_destroy(a);
}

MU_TEST(gf_elem_from_array_test3) {
  uint8_t coeff[] = {1, 1};
  GF_elem_t *a = GF_elem_from_array(1, coeff, &GF2_5);

  mu_check(GF_eq(a->GF, &GF2_5));
  mu_check(a->poly->deg == 1);
  mu_check(a->poly->coeff[0] == 1);
  mu_check(a->poly->coeff[1] == 1);
  mu_check(a->poly->coeff[2] == 0);
  mu_check(a->poly->coeff[3] == 0);
  mu_check(a->poly->coeff[4] == 0);

  GF_elem_destroy(a);
}

MU_TEST(gf_elem_from_array_test4) {
  uint8_t coeff[] = {1, 1, 0, 1, 0};
  GF_elem_t *a = GF_elem_from_array(3, coeff, &GF2_5);

  mu_check(GF_eq(a->GF, &GF2_5));
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
  mu_check(neutral_2_2->GF->I->deg == 2);
  mu_check(poly_eq(neutral_2_2->GF->I, &IGF2_2));
  mu_check(neutral_2_2->poly->deg == 0);
  mu_check(poly_eq(neutral_2_2->poly, zero_2_2));

  GF_elem_destroy(neutral_2_2);
  poly_destroy(zero_2_2);

  GF_elem_t *neutral_2_5 = GF_elem_get_neutral(&GF2_5);
  mu_check(neutral_2_5->GF->p == 2);
  mu_check(neutral_2_5->GF->I->deg == 5);
  mu_check(poly_eq(neutral_2_5->GF->I, &IGF2_5));
  mu_check(neutral_2_5->poly->deg == 0);
  mu_check(poly_eq(neutral_2_5->poly, zero_2_5));

  GF_elem_destroy(neutral_2_5);
  poly_destroy(zero_2_5);
}

MU_TEST(gf_get_unity_test) {
  uint8_t one_coeff_2_5[5] = {1};
  poly_t one_2_5 = {.deg = 0, .coeff = one_coeff_2_5};

  GF_elem_t *unity_2_5 = GF_elem_get_unity(&GF2_5);
  mu_check(poly_eq(unity_2_5->GF->I, &IGF2_5));
  mu_check(unity_2_5->poly->deg == 0);
  mu_check(poly_eq(unity_2_5->poly, &one_2_5));
  GF_elem_destroy(unity_2_5);
}

MU_TEST(gf_sum_test1) {
  uint8_t coeff_a[] = {1, 1, 1, 0, 0};
  GF_elem_t *a = GF_elem_from_array(2, coeff_a, &GF2_5);
  mu_check(a->poly->deg == 2);

  uint8_t coeff_b[] = {0, 1, 1, 0, 1};
  GF_elem_t *b = GF_elem_from_array(4, coeff_b, &GF2_5);
  mu_check(b->poly->deg == 4);

  GF_elem_sum(a, *a, *b);
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
  uint8_t coeff_a[] = {1, 1, 1, 1, 1};
  GF_elem_t *a = GF_elem_from_array(4, coeff_a, &GF2_5);

  GF_elem_sum(a, *a, *a);
  mu_check(a->poly->deg == 0);
  mu_check(a->poly->coeff[0] == 0);
  mu_check(a->poly->coeff[1] == 0);
  mu_check(a->poly->coeff[2] == 0);
  mu_check(a->poly->coeff[3] == 0);
  mu_check(a->poly->coeff[4] == 0);

  GF_elem_destroy(a);
}

MU_TEST(gf_diff_test1) {
  uint8_t coeff_a[5] = {1, 1, 1, 1, 1};
  GF_elem_t *a = GF_elem_from_array(4, coeff_a, &GF2_5);

  GF_elem_diff(a, *a, *a);
  mu_check(a->poly->deg == 0);
  mu_check(a->poly->coeff[0] == 0);
  mu_check(a->poly->coeff[1] == 0);
  mu_check(a->poly->coeff[2] == 0);
  mu_check(a->poly->coeff[3] == 0);
  mu_check(a->poly->coeff[4] == 0);

  GF_elem_destroy(a);
}

MU_TEST(gf_diff_test2) {
  uint8_t coeff_a[] = {0, 5, 3, 1};
  GF_elem_t *a = GF_elem_from_array(3, coeff_a, &GF7_4);

  uint8_t coeff_b[] = {3, 0, 6, 0};
  GF_elem_t *b = GF_elem_from_array(2, coeff_b, &GF7_4);

  GF_elem_diff(a, *a, *b);
  mu_check(a->poly->deg == 3);
  mu_check(a->poly->coeff[0] == 4);
  mu_check(a->poly->coeff[1] == 5);
  mu_check(a->poly->coeff[2] == 4);
  mu_check(a->poly->coeff[3] == 1);

  GF_elem_destroy(a);
  GF_elem_destroy(b);
}

MU_TEST(poly_mul_test1) {
  uint8_t coeff[] = {1, 1};
  poly_t *a = poly_from_array(1, coeff);
  poly_t *res = poly_create_zero(3);
  poly_mul(res, *a, *a, 3);

  mu_check(res->deg = 2);
  mu_check(res->coeff[0] == 1);
  mu_check(res->coeff[1] == 2);
  mu_check(res->coeff[2] == 1);

  poly_destroy(res);
  poly_destroy(a);
}

MU_TEST(poly_mul_test2) {
  uint8_t coeff_a[] = {1, 0, 1, 5, 15};
  poly_t *a = poly_from_array(4, coeff_a);

  uint8_t coeff_b[] = {1, 0, 1, 0, 0};
  poly_t *b = poly_from_array(2, coeff_b);

  poly_t *res = poly_create_zero(7);
  poly_mul(res, *a, *b, 17);

  mu_check(res->deg = 6);
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
  uint8_t coeff_a[] = {1, 0, 1, 5, 15};
  GF_elem_t *a = GF_elem_from_array(4, coeff_a, &GF2_5);

  uint8_t coeff_b[] = {1, 0, 1, 0};
  GF_elem_t *b = GF_elem_from_array(2, coeff_b, &GF2_5);

  GF_elem_t *res = GF_elem_get_neutral(&GF2_5);

  // 15x^6 + 5x^5 + 16x^4 + 5x^3 + 2x^2 + 1
  GF_elem_prod(res, *a, *b);

  mu_check(res->poly->deg = 2);
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
  uint8_t coeff[] = {1, 1};
  GF_elem_t *a = GF_elem_from_array(2, coeff, &GF2_5);
  GF_elem_t *res = GF_elem_get_neutral(&GF2_5);

  GF_elem_prod(res, *a, *a);

  mu_check(res->poly->deg = 2);
  mu_check(res->poly->coeff[0] == 1);
  mu_check(res->poly->coeff[1] == 0);
  mu_check(res->poly->coeff[2] == 1);
  mu_check(res->poly->coeff[3] == 0);
  mu_check(res->poly->coeff[4] == 0);

  GF_elem_destroy(a);
  GF_elem_destroy(res);
}

MU_TEST(poly_fpowm_test1) {
  uint8_t coeff[2] = {1, 1};
  poly_t *a = poly_from_array(1, coeff);
  poly_t *res = poly_create_zero(2);
  poly_fpowm(res, *a, 2, IGF2_2, 2);

  mu_check(res->deg = 1);
  mu_check(res->coeff[0] == 0);
  mu_check(res->coeff[1] == 1);

  poly_destroy(res);
  poly_destroy(a);
}

MU_TEST(poly_fpowm_test2) {
  uint8_t coeff[2] = {1, 1};
  poly_t *a = poly_from_array(1, coeff);
  poly_t *res = poly_create_zero(2);
  poly_fpowm(res, *a, 2, IGF2_2, 2);

  mu_check(res->deg = 1);
  mu_check(res->coeff[0] == 0);
  mu_check(res->coeff[1] == 1);

  poly_destroy(res);
  poly_destroy(a);
}

MU_TEST(gf_get_inverse) {
  uint8_t coeff_a[3] = {1, 0, 0};
  GF_elem_t *a = GF_elem_from_array(1, coeff_a, &GF2_3);
  GF_elem_t *inv_a = GF_elem_get_inverse(*a);
  mu_check(inv_a->poly->deg == 0);
  mu_check(inv_a->poly->coeff[0] == 1);
  mu_check(inv_a->poly->coeff[1] == 0);
  mu_check(inv_a->poly->coeff[2] == 0);
  GF_elem_destroy(a);
  GF_elem_destroy(inv_a);

  uint8_t coeff_b[3] = {0, 1, 0};
  GF_elem_t *b = GF_elem_from_array(1, coeff_b, &GF2_3);
  GF_elem_t *inv_b = GF_elem_get_inverse(*b);
  mu_check(inv_b->poly->deg == 2);
  mu_check(inv_b->poly->coeff[0] == 1);
  mu_check(inv_b->poly->coeff[1] == 0);
  mu_check(inv_b->poly->coeff[2] == 1);
  GF_elem_destroy(b);
  GF_elem_destroy(inv_b);

  uint8_t coeff_c[3] = {0, 0, 1};
  GF_elem_t *c = GF_elem_from_array(2, coeff_c, &GF2_3);
  GF_elem_t *inv_c = GF_elem_get_inverse(*c);
  mu_check(inv_c->poly->deg == 2);
  mu_check(inv_c->poly->coeff[0] == 1);
  mu_check(inv_c->poly->coeff[1] == 1);
  mu_check(inv_c->poly->coeff[2] == 1);
  GF_elem_destroy(c);
  GF_elem_destroy(inv_c);

  uint8_t coeff_d[3] = {1, 1, 0};
  GF_elem_t *d = GF_elem_from_array(1, coeff_d, &GF2_3);
  GF_elem_t *inv_d = GF_elem_get_inverse(*d);
  mu_check(inv_d->poly->deg == 2);
  mu_check(inv_d->poly->coeff[0] == 0);
  mu_check(inv_d->poly->coeff[1] == 1);
  mu_check(inv_d->poly->coeff[2] == 1);
  GF_elem_destroy(d);
  GF_elem_destroy(inv_d);

  uint8_t coeff_e[3] = {1, 0, 1};
  GF_elem_t *e = GF_elem_from_array(2, coeff_e, &GF2_3);
  GF_elem_t *inv_e = GF_elem_get_inverse(*e);
  mu_check(inv_e->poly->deg == 1);
  mu_check(inv_e->poly->coeff[0] == 0);
  mu_check(inv_e->poly->coeff[1] == 1);
  mu_check(inv_e->poly->coeff[2] == 0);
  GF_elem_destroy(e);
  GF_elem_destroy(inv_e);

  uint8_t coeff_f[3] = {0, 1, 1};
  GF_elem_t *f = GF_elem_from_array(2, coeff_f, &GF2_3);
  GF_elem_t *inv_f = GF_elem_get_inverse(*f);
  mu_check(inv_f->poly->deg == 1);
  mu_check(inv_f->poly->coeff[0] == 1);
  mu_check(inv_f->poly->coeff[1] == 1);
  mu_check(inv_f->poly->coeff[2] == 0);
  GF_elem_destroy(f);
  GF_elem_destroy(inv_f);

  uint8_t coeff_g[3] = {1, 1, 1};
  GF_elem_t *g = GF_elem_from_array(2, coeff_g, &GF2_3);
  GF_elem_t *inv_g = GF_elem_get_inverse(*g);
  mu_check(inv_g->poly->deg == 2);
  mu_check(inv_g->poly->coeff[0] == 0);
  mu_check(inv_g->poly->coeff[1] == 0);
  mu_check(inv_g->poly->coeff[2] == 1);
  GF_elem_destroy(g);
  GF_elem_destroy(inv_g);
}

MU_TEST_SUITE(gf_tests) {
  MU_RUN_TEST(gf_eq_test);
  MU_RUN_TEST(gf_get_neutral_test);
  MU_RUN_TEST(gf_get_unity_test);
  MU_RUN_TEST(gf_elem_from_array_test1);
  MU_RUN_TEST(gf_elem_from_array_test2);
  MU_RUN_TEST(gf_elem_from_array_test3);
  MU_RUN_TEST(gf_elem_from_array_test4);
  MU_RUN_TEST(gf_get_inverse);
}

MU_TEST_SUITE(gf_arithmetic_tests) {
  MU_RUN_TEST(poly_div_test_case1);
  MU_RUN_TEST(poly_div_test_case2);
  MU_RUN_TEST(poly_mul_test1);
  MU_RUN_TEST(poly_mul_test2);
  MU_RUN_TEST(gf_sum_test1);
  MU_RUN_TEST(gf_sum_test2);
  MU_RUN_TEST(gf_diff_test1);
  MU_RUN_TEST(gf_diff_test2);
  MU_RUN_TEST(gf_prod_test1);
  MU_RUN_TEST(gf_prod_test2);
  MU_RUN_TEST(poly_fpowm_test1);
  MU_RUN_TEST(poly_fpowm_test2);
}

int main() {
  MU_RUN_TEST(poly_eq_test);
  MU_RUN_TEST(poly_init_finit_test);
  MU_RUN_SUITE(gf_tests);
  MU_RUN_SUITE(gf_arithmetic_tests);
  MU_REPORT();
  return MU_EXIT_CODE;
}
