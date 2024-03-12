#include "minunit.h"

MU_TEST(check) { mu_check(1 == 1); }

int main() {
  MU_RUN_TEST(check);
  MU_REPORT();
  return MU_EXIT_CODE;
}
