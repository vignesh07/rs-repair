// GF(256) with the ISA-L / Hadoop / Backblaze polynomial x^8+x^4+x^3+x^2+1 (0x11D), generator 2.
#ifndef GF_H
#define GF_H
#include <stdint.h>
static uint8_t gf_exp[512], gf_log[256];
static uint8_t gf_mul_tab[256][256];
static void gf_init(void) {
    int x = 1;
    for (int i = 0; i < 255; i++) { gf_exp[i] = (uint8_t)x; gf_log[x] = (uint8_t)i; x <<= 1; if (x & 0x100) x ^= 0x11D; }
    for (int i = 255; i < 512; i++) gf_exp[i] = gf_exp[i - 255];
    for (int a = 0; a < 256; a++) for (int b = 0; b < 256; b++)
        gf_mul_tab[a][b] = (a && b) ? gf_exp[gf_log[a] + gf_log[b]] : 0;
}
static inline uint8_t gf_mul(uint8_t a, uint8_t b) { return gf_mul_tab[a][b]; }
static inline uint8_t gf_inv(uint8_t a) { return gf_exp[255 - gf_log[a]]; }
// rank over GF(2) of m bytes viewed as vectors in GF(2)^8
static inline int rank_bytes(const uint8_t *v, int m) {
    uint8_t basis[8] = {0}; int rk = 0;
    for (int i = 0; i < m; i++) {
        uint8_t x = v[i];
        while (x) { int b = 31 - __builtin_clz((unsigned)x); if (basis[b]) x ^= basis[b]; else { basis[b] = x; rk++; break; } }
    }
    return rk;
}
#endif
