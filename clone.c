#include <stdint.h>
#include <stdio.h>

#define MIX_CONST1 0xff51afd7ed558ccdULL
#define MIX_CONST2 0xc4ceb9fe1a85ec53ULL

/* PlantSDS encoding: A = 001, C = 110, G = 011, T = 100 */
static const int8_t BASE_LOOKUP[256] = {
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, 1,  -1, 6,  -1, -1, -1, 3,  -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, 4,  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, 1,  -1, 6,  -1, -1, -1, 3,  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, 4,  -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1};

static uint64_t mix_hash(__uint128_t hash, uint64_t seed) {
  __uint128_t mixed_hash =
      (hash ^ MIX_CONST1) * ((__uint128_t)seed ^ MIX_CONST2);
  return (uint64_t)(mixed_hash ^ (mixed_hash >> 64));
}

static uint64_t reverse_bits64(uint64_t hash) {
#if defined(__aarch64__)
  uint64_t r;
  __asm__("rbit %0, %1" : "=r"(r) : "r"(hash));
  return r;
#else
  uint64_t r = __builtin_bswap64(n);
  r = ((r & 0x5555555555555555ULL) << 1) | ((r & 0xAAAAAAAAAAAAAAAAULL) >> 1);
  r = ((r & 0x3333333333333333ULL) << 2) | ((r & 0xCCCCCCCCCCCCCCCCULL) >> 2);
  r = ((r & 0x0F0F0F0F0F0F0F0FULL) << 4) | ((r & 0xF0F0F0F0F0F0F0F0ULL) >> 4);
  return r;
#endif
}

static __uint128_t reverse_bits128(__uint128_t hash) {
  return ((__uint128_t)reverse_bits64((uint64_t)hash) << 64) |
         reverse_bits64((uint64_t)(hash >> 64));
}
