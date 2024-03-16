#pragma once

#include <stdint.h>

#define MIN(A, B) ((A) < (B) ? (A) : (B))

int8_t iabs(int8_t x);

int8_t eu_mod(int8_t x, int8_t y);

int8_t get_complement_mod_p(int8_t a, int8_t p);
