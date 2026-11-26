// r = 4 search restricted to check polynomials g(y) = a*y^3 + b*y + c (no y^2 term), y = x - a*.
// For points in a GF(2)-subspace these make g(d)*d GF(2)-linear in d. Scheme stored as a[i], b[i] for basis e_i (c = e_i).
// Proposals: (1) change one byte; (2) force g_u to vanish on {d1, d2, d1^d2} (a 2-dim subspace minus 0);
//            (3) force g_u to vanish at one helper d1 with a random cubic coefficient.
// usage: sa3 n node iters seed T0 T1 restarts pstruct p0..p(n-1)   [env SA_INIT = "a0,b0,a1,b1,...,a7,b7"]
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "gf.h"
static uint64_t s[2];
static inline uint64_t rng(void) { uint64_t a = s[0], b = s[1]; s[0] = b; a ^= a << 23; s[1] = a ^ b ^ (a >> 17) ^ (b >> 26); return s[1] + b; }
static inline double rnd(void) { return (rng() >> 11) * (1.0 / 9007199254740992.0); }
static int n, node, H; static uint8_t A[64], hd[64], d1p[64], d3p[64];
static int eval(const uint8_t *a, const uint8_t *b, int *rk) {
    int tot = 0;
    for (int h = 0; h < H; h++) { uint8_t col[8]; for (int i = 0; i < 8; i++) col[i] = (uint8_t)((1 << i) ^ gf_mul(a[i], d3p[h]) ^ gf_mul(b[i], d1p[h])); rk[h] = rank_bytes(col, 8); tot += rk[h]; }
    return tot;
}
int main(int argc, char **argv) {
    gf_init();
    n = atoi(argv[1]); node = atoi(argv[2]);
    long iters = atol(argv[3]); uint64_t seed = strtoull(argv[4], 0, 10);
    double T0 = atof(argv[5]), T1 = atof(argv[6]); int restarts = atoi(argv[7]); double pstruct = atof(argv[8]);
    for (int i = 0; i < n; i++) A[i] = (uint8_t)atoi(argv[9 + i]);
    H = 0; for (int j = 0; j < n; j++) if (j != node) { uint8_t d = A[j] ^ A[node]; hd[H] = d; d1p[H] = d; d3p[H] = gf_mul(d, gf_mul(d, d)); H++; }
    s[0] = seed * 0x9E3779B97F4A7C15ULL + 1; s[1] = seed ^ 0xD1B54A32D192ED03ULL; for (int i = 0; i < 20; i++) rng();
    int best = 1 << 30; uint8_t ba[8], bb[8]; memset(ba, 0, 8); memset(bb, 0, 8);
    const char *init = getenv("SA_INIT");
    for (int rs = 0; rs < restarts; rs++) {
        uint8_t a[8], b[8], na[8], nb[8]; int rk[64], nrk[64];
        if (init) { const char *p = init; for (int i = 0; i < 8; i++) { a[i] = (uint8_t)strtol(p, (char **)&p, 10); if (*p == ',') p++; b[i] = (uint8_t)strtol(p, (char **)&p, 10); if (*p == ',') p++; } }
        else for (int i = 0; i < 8; i++) { a[i] = (uint8_t)(rng() & 255); b[i] = (uint8_t)(rng() & 255); }
        int cur = eval(a, b, rk), lbest = cur;
        if (cur < best) { best = cur; memcpy(ba, a, 8); memcpy(bb, b, 8); }
        for (long it = 0; it < iters; it++) {
            double T = T0 * pow(T1 / T0, (double)it / (double)iters);
            memcpy(na, a, 8); memcpy(nb, b, 8);
            double x = rnd();
            if (x >= pstruct) { int i = (int)(rng() & 7); if (rng() & 1) na[i] = (uint8_t)(rng() & 255); else nb[i] = (uint8_t)(rng() & 255); }
            else {
                uint8_t u = (uint8_t)(1 + rng() % 255), ta, tb;
                if (rng() % 4) {                       // roots d1, d2, d1^d2 : g = u*(y^3 + A y + B)/B with y^4+Ay^2+By = prod (y - t), t in T
                    uint8_t d1 = hd[rng() % (uint64_t)H], d2; do { d2 = hd[rng() % (uint64_t)H]; } while (d2 == d1);
                    uint8_t d3 = d1 ^ d2;
                    uint8_t B = gf_mul(d1, gf_mul(d2, d3));                                   // product of the three nonzero roots
                    uint8_t Acoef = (uint8_t)(gf_mul(d1, d2) ^ gf_mul(d1, d3) ^ gf_mul(d2, d3)); // e2 of the roots (e1 = 0)
                    uint8_t Bi = gf_inv(B);
                    ta = gf_mul(u, Bi); tb = gf_mul(u, gf_mul(Acoef, Bi));
                } else {                               // one prescribed root d1, random cubic coefficient
                    uint8_t d1 = hd[rng() % (uint64_t)H]; ta = (uint8_t)(rng() & 255);
                    tb = gf_mul((uint8_t)(u ^ gf_mul(ta, gf_mul(d1, gf_mul(d1, d1)))), gf_inv(d1));
                }
                uint8_t ca = 0, cb = 0; for (int i = 0; i < 8; i++) if (u >> i & 1) { ca ^= a[i]; cb ^= b[i]; }
                uint8_t lam; do { lam = (uint8_t)(1 + rng() % 255); } while ((__builtin_popcount(lam & u) & 1) == 0);
                for (int i = 0; i < 8; i++) if (lam >> i & 1) { na[i] ^= (uint8_t)(ta ^ ca); nb[i] ^= (uint8_t)(tb ^ cb); }
            }
            int nw = eval(na, nb, nrk), d = nw - cur;
            if (d <= 0 || rnd() < exp(-d / T)) {
                memcpy(a, na, 8); memcpy(b, nb, 8); cur = nw;
                if (cur < lbest) lbest = cur;
                if (cur < best) { best = cur; memcpy(ba, a, 8); memcpy(bb, b, 8); }
            }
        }
        fprintf(stderr, "restart %d: local best %d, global best %d\n", rs, lbest, best);
    }
    printf("BEST n=%d r=4 node=%d bw=%d  ab=", n, node, best);
    for (int i = 0; i < 8; i++) printf("%d,%d%s", ba[i], bb[i], i == 7 ? "\n" : ",");
    return 0;
}
