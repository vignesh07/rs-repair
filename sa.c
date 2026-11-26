// Simulated annealing over ALL linear repair schemes for one node of an RS code over GF(256), base field GF(2).
// Scheme: g_u(x) = u + sum_{m=1}^{r-1} L_m(u) (x-a*)^m, u in F, L_m GF(2)-linear. Stored as c[i][m] = L_m(e_i), e_i = 1<<i.
// Bandwidth = sum over helpers j of rank_GF(2){ g_{e_i}(a_j) : i=0..7 }.
// usage: sa n r node iters seed T0 T1 restarts  p0 p1 ... p(n-1)
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "gf.h"
static uint64_t s[2];
static inline uint64_t rng(void) { uint64_t a = s[0], b = s[1]; s[0] = b; a ^= a << 23; s[1] = a ^ b ^ (a >> 17) ^ (b >> 26); return s[1] + b; }
static inline double rnd(void) { return (rng() >> 11) * (1.0 / 9007199254740992.0); }
int main(int argc, char **argv) {
    gf_init();
    int n = atoi(argv[1]), r = atoi(argv[2]), node = atoi(argv[3]);
    long iters = atol(argv[4]); uint64_t seed = strtoull(argv[5], 0, 10);
    double T0 = atof(argv[6]), T1 = atof(argv[7]); int restarts = atoi(argv[8]);
    uint8_t A[64]; for (int i = 0; i < n; i++) A[i] = (uint8_t)atoi(argv[9 + i]);
    int H = n - 1; uint8_t dpow[64][8];
    for (int j = 0, h = 0; j < n; j++) if (j != node) { uint8_t d = A[j] ^ A[node], p = 1; for (int m = 1; m < r; m++) { p = gf_mul(p, d); dpow[h][m] = p; } h++; }
    s[0] = seed * 0x9E3779B97F4A7C15ULL + 1; s[1] = seed ^ 0xD1B54A32D192ED03ULL; for (int i = 0; i < 20; i++) rng();
    int best = 1 << 30; uint8_t bestc[8][8];
    for (int rs = 0; rs < restarts; rs++) {
        uint8_t c[8][8], v[8][64]; int rk[64], cur = 0;
        for (int i = 0; i < 8; i++) for (int m = 1; m < r; m++) c[i][m] = (uint8_t)(rng() & 255);
        for (int i = 0; i < 8; i++) for (int h = 0; h < H; h++) { uint8_t x = (uint8_t)(1 << i); for (int m = 1; m < r; m++) x ^= gf_mul(c[i][m], dpow[h][m]); v[i][h] = x; }
        for (int h = 0; h < H; h++) { uint8_t col[8]; for (int i = 0; i < 8; i++) col[i] = v[i][h]; rk[h] = rank_bytes(col, 8); cur += rk[h]; }
        int lbest = cur;
        for (long it = 0; it < iters; it++) {
            double T = T0 * pow(T1 / T0, (double)it / (double)iters);
            int i = (int)(rng() & 7), m = 1 + (int)(rng() % (uint64_t)(r - 1)); uint8_t nv = (uint8_t)(rng() & 255), old = c[i][m];
            if (nv == old) continue;
            uint8_t diff = nv ^ old, nvv[64]; int nrk[64], nw = 0;
            for (int h = 0; h < H; h++) {
                nvv[h] = v[i][h] ^ gf_mul(diff, dpow[h][m]);
                uint8_t col[8]; for (int q = 0; q < 8; q++) col[q] = (q == i) ? nvv[h] : v[q][h];
                nrk[h] = rank_bytes(col, 8); nw += nrk[h];
            }
            int d = nw - cur;
            if (d <= 0 || rnd() < exp(-d / T)) {
                c[i][m] = nv; cur = nw; for (int h = 0; h < H; h++) { v[i][h] = nvv[h]; rk[h] = nrk[h]; }
                if (cur < lbest) lbest = cur;
                if (cur < best) { best = cur; memcpy(bestc, c, sizeof(c)); }
            }
        }
        fprintf(stderr, "restart %d: local best %d, global best %d\n", rs, lbest, best);
    }
    printf("BEST n=%d r=%d node=%d bw=%d  c=", n, r, node, best);
    for (int i = 0; i < 8; i++) for (int m = 1; m < r; m++) printf("%d%s", bestc[i][m], (i == 7 && m == r - 1) ? "\n" : ",");
    return 0;
}
