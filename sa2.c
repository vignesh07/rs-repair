// General-space simulated annealing with structured proposals.
// Same scheme representation as sa.c (c[i][m] = L_m(e_i)); proposals:
//   (1) change one byte;  (2) pick u != 0 and r-1 helpers R, force g_u := u * P_R(x)/P_R(a*) by a rank-one update
//   lambda (x) Delta with lambda(u)=1;  (3) pick u and r-2 helpers plus a random extra factor.
// usage: sa2 n r node iters seed T0 T1 restarts pstruct  p0 ... p(n-1)   [init c bytes via env SA_INIT="a,b,c,..."]
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "gf.h"
static uint64_t s[2];
static inline uint64_t rng(void) { uint64_t a = s[0], b = s[1]; s[0] = b; a ^= a << 23; s[1] = a ^ b ^ (a >> 17) ^ (b >> 26); return s[1] + b; }
static inline double rnd(void) { return (rng() >> 11) * (1.0 / 9007199254740992.0); }
static int n, r, node, H; static uint8_t A[64], hd[64], dpow[64][8];
static int eval(uint8_t c[8][8], uint8_t v[8][64], int *rk) {
    int tot = 0;
    for (int i = 0; i < 8; i++) for (int h = 0; h < H; h++) { uint8_t x = (uint8_t)(1 << i); for (int m = 1; m < r; m++) x ^= gf_mul(c[i][m], dpow[h][m]); v[i][h] = x; }
    for (int h = 0; h < H; h++) { uint8_t col[8]; for (int i = 0; i < 8; i++) col[i] = v[i][h]; rk[h] = rank_bytes(col, 8); tot += rk[h]; }
    return tot;
}
int main(int argc, char **argv) {
    gf_init();
    n = atoi(argv[1]); r = atoi(argv[2]); node = atoi(argv[3]);
    long iters = atol(argv[4]); uint64_t seed = strtoull(argv[5], 0, 10);
    double T0 = atof(argv[6]), T1 = atof(argv[7]); int restarts = atoi(argv[8]); double pstruct = atof(argv[9]);
    for (int i = 0; i < n; i++) A[i] = (uint8_t)atoi(argv[10 + i]);
    H = 0; for (int j = 0; j < n; j++) if (j != node) { uint8_t d = A[j] ^ A[node], p = 1; hd[H] = d; for (int m = 1; m < r; m++) { p = gf_mul(p, d); dpow[H][m] = p; } H++; }
    s[0] = seed * 0x9E3779B97F4A7C15ULL + 1; s[1] = seed ^ 0xD1B54A32D192ED03ULL; for (int i = 0; i < 20; i++) rng();
    int best = 1 << 30; uint8_t bestc[8][8]; memset(bestc, 0, sizeof(bestc));
    const char *init = getenv("SA_INIT");
    for (int rs = 0; rs < restarts; rs++) {
        uint8_t c[8][8], v[8][64], nc[8][8], nv[8][64]; int rk[64], nrk[64];
        memset(c, 0, sizeof(c));
        if (init) { const char *p = init; for (int i = 0; i < 8; i++) for (int m = 1; m < r; m++) { c[i][m] = (uint8_t)strtol(p, (char **)&p, 10); if (*p == ',') p++; } }
        else for (int i = 0; i < 8; i++) for (int m = 1; m < r; m++) c[i][m] = (uint8_t)(rng() & 255);
        int cur = eval(c, v, rk), lbest = cur;
        if (cur < best) { best = cur; memcpy(bestc, c, sizeof(c)); }
        for (long it = 0; it < iters; it++) {
            double T = T0 * pow(T1 / T0, (double)it / (double)iters);
            memcpy(nc, c, sizeof(c));
            if (rnd() >= pstruct) {
                int i = (int)(rng() & 7), m = 1 + (int)(rng() % (uint64_t)(r - 1)); nc[i][m] = (uint8_t)(rng() & 255);
            } else {
                uint8_t u = (uint8_t)(1 + rng() % 255);
                // target polynomial t(x) = u * prod_{h in R}(x - a_h) / prod_{h in R}(a* - a_h), expressed in powers of y = x - a*:
                // prod (y + d_h) / prod d_h ; coefficients of y^m, m = 0..r-1 (constant term = 1 -> times u)
                int R[8], cnt = 0; int want = r - 1; if (r > 2 && (rng() & 3) == 0) want = r - 2;
                while (cnt < want) { int h = (int)(rng() % (uint64_t)H), dup = 0; for (int q = 0; q < cnt; q++) if (R[q] == h) dup = 1; if (!dup) R[cnt++] = h; }
                uint8_t poly[8] = {1, 0, 0, 0, 0, 0, 0, 0}; int deg = 0;
                for (int q = 0; q < cnt; q++) { uint8_t d = hd[R[q]], di = gf_inv(d);  // multiply by (y + d)/d = 1 + y/d
                    for (int m = deg + 1; m >= 1; m--) poly[m] ^= gf_mul(poly[m - 1], di); deg++; }
                if (want == r - 2) { uint8_t e = (uint8_t)(rng() & 255);                     // extra factor (1 + e*y), e arbitrary
                    for (int m = deg + 1; m >= 1; m--) poly[m] ^= gf_mul(poly[m - 1], e); deg++; }
                uint8_t target[8], curg[8]; for (int m = 1; m < r; m++) { target[m] = gf_mul(u, poly[m]); curg[m] = 0; }
                for (int i = 0; i < 8; i++) if (u >> i & 1) for (int m = 1; m < r; m++) curg[m] ^= c[i][m];
                uint8_t lam; do { lam = (uint8_t)(1 + rng() % 255); } while ((__builtin_popcount(lam & u) & 1) == 0);
                for (int i = 0; i < 8; i++) if (lam >> i & 1) for (int m = 1; m < r; m++) nc[i][m] ^= (uint8_t)(target[m] ^ curg[m]);
            }
            int nw = eval(nc, nv, nrk), d = nw - cur;
            if (d <= 0 || rnd() < exp(-d / T)) {
                memcpy(c, nc, sizeof(c)); memcpy(v, nv, sizeof(v)); memcpy(rk, nrk, sizeof(int) * H); cur = nw;
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
