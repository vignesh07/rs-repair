// adj.c — subspace family for r = 4 over GF(256): g(y) = a*y^3 + b*y + c, y = x - a*, so g(d)*d = a d^4 + b d^2 + c d is GF(2)-linear in d.
//   S_Y = { (a,b,c) : a d^4 + b d^2 + c d in Y for all d in V },  V = 4-dim point subspace, Y = 4-dim subspace.
//
// Claim tested here ("adjoint-space criterion"):
//   dim S_Y            = 8 + delta_0,   delta_0 = dim { lam : lam * Uhat_V <= Yperp }
//   rank at helper d   = 4 + delta_0 - delta_d,   delta_d = dim { lam : lam * U_d <= Yperp }
// where W_d = { v^2 + d v : v in V } (3-dim), U_d = span of the 2x2 Moore minors of a basis of W_d
// (= products of the nonzero elements of the 2-dim subspaces of W_d), Uhat_V = same with 3x3 minors of V,
// and Yperp is the trace-dual of Y.
//
//   adj A          : compare the criterion with direct computation for every 4-dim Y, V = {0..15}
//   adj B          : best scheme in the family over every 4-dim V containing 1 (all V up to scaling), n = 16
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "gf.h"

static uint8_t TR[256], FV[256];
static uint8_t SUB[210000][4];
static int NSUB;
static int FILT;                 // adj B v1 v2 v3 v4 : only the point subspace spanned by these four elements
static uint8_t FILTV[4];

static void init_all(void) {
    gf_init();
    for (int x = 0; x < 256; x++) {
        uint8_t t = (uint8_t)x, s = (uint8_t)x;
        for (int i = 1; i < 8; i++) { t = gf_mul(t, t); s ^= t; }
        if (s > 1) { fprintf(stderr, "trace not in GF(2)\n"); exit(1); }
        TR[x] = s;
    }
    for (int w = 0; w < 256; w++) {          // Tr(w*lam) = parity(FV[w] & lam)
        uint8_t f = 0;
        for (int k = 0; k < 8; k++) if (TR[gf_mul((uint8_t)w, (uint8_t)(1 << k))]) f |= (uint8_t)(1 << k);
        FV[w] = f;
    }
}

static int nullspace(const uint32_t *rows_in, int nrows, int ncols, uint32_t *out) {
    uint32_t r[32]; int pivcol[32], rank = 0, ispiv[32] = {0};
    memcpy(r, rows_in, (size_t)nrows * sizeof(uint32_t));
    for (int c = 0; c < ncols && rank < nrows; c++) {
        int p = -1;
        for (int i = rank; i < nrows; i++) if (r[i] >> c & 1) { p = i; break; }
        if (p < 0) continue;
        uint32_t t = r[p]; r[p] = r[rank]; r[rank] = t;
        for (int i = 0; i < nrows; i++) if (i != rank && (r[i] >> c & 1)) r[i] ^= r[rank];
        pivcol[rank++] = c;
    }
    for (int i = 0; i < rank; i++) ispiv[pivcol[i]] = 1;
    int nd = 0;
    for (int c = 0; c < ncols; c++) if (!ispiv[c]) {
        uint32_t v = 1u << c;
        for (int i = 0; i < rank; i++) if (r[i] >> c & 1) v |= 1u << pivcol[i];
        out[nd++] = v;
    }
    return nd;
}

static int perp(const uint8_t *B, int m, uint8_t *out) {
    uint32_t rows[8], o[8];
    for (int i = 0; i < m; i++) rows[i] = FV[B[i]];
    int d = nullspace(rows, m, 8, o);
    for (int i = 0; i < d; i++) out[i] = (uint8_t)o[i];
    return d;
}

static void elements(const uint8_t *B, int k, uint8_t *el) {     // the 2^k - 1 nonzero elements
    int n = 0;
    for (int m = 1; m < (1 << k); m++) { uint8_t x = 0; for (int i = 0; i < k; i++) if (m >> i & 1) x ^= B[i]; el[n++] = x; }
}

static void enum_sub4(void) {
    NSUB = 0;
    for (int mask = 0; mask < 256; mask++) {
        if (__builtin_popcount((unsigned)mask) != 4) continue;
        int p[4], np = 0, fpos[4][8], nf[4], tot = 0;
        for (int b = 7; b >= 0; b--) if (mask >> b & 1) p[np++] = b;
        for (int i = 0; i < 4; i++) { nf[i] = 0; for (int b = 0; b < p[i]; b++) if (!(mask >> b & 1)) fpos[i][nf[i]++] = b; tot += nf[i]; }
        for (uint32_t x = 0; x < (1u << tot); x++) {
            uint32_t y = x;
            for (int i = 0; i < 4; i++) {
                uint8_t r = (uint8_t)(1 << p[i]);
                for (int t = 0; t < nf[i]; t++) { if (y & 1) r |= (uint8_t)(1 << fpos[i][t]); y >>= 1; }
                SUB[NSUB][i] = r;
            }
            NSUB++;
        }
    }
}

typedef struct { int dimS, valid, rk[15], tot; uint32_t S[24]; } Res;

// direct computation from the definition
static void evalY(const uint8_t *Vb, const uint8_t *Vel, const uint8_t *Yp, Res *R) {
    uint32_t rows[16]; int n = 0; uint8_t tmp[24];
    for (int i = 0; i < 4; i++) {
        uint8_t v = Vb[i], v2 = gf_mul(v, v), v4 = gf_mul(v2, v2);
        for (int j = 0; j < 4; j++) {
            uint8_t z = Yp[j];
            rows[n++] = (uint32_t)FV[gf_mul(z, v4)] | ((uint32_t)FV[gf_mul(z, v2)] << 8) | ((uint32_t)FV[gf_mul(z, v)] << 16);
        }
    }
    R->dimS = nullspace(rows, 16, 24, R->S);
    for (int s = 0; s < R->dimS; s++) tmp[s] = (uint8_t)(R->S[s] >> 16);
    R->valid = (rank_bytes(tmp, R->dimS) == 8);
    R->tot = 0;
    for (int e = 0; e < 15; e++) {
        uint8_t d = Vel[e], d3 = gf_mul(d, gf_mul(d, d));
        for (int s = 0; s < R->dimS; s++) {
            uint8_t a = (uint8_t)(R->S[s] & 255), b = (uint8_t)((R->S[s] >> 8) & 255), c = (uint8_t)(R->S[s] >> 16);
            tmp[s] = gf_mul(a, d3) ^ gf_mul(b, d) ^ c;
        }
        R->rk[e] = rank_bytes(tmp, R->dimS); R->tot += R->rk[e];
    }
}

static void Ud(const uint8_t *Vb, uint8_t d, uint8_t *u) {
    uint8_t w[4]; int nw = 0;
    for (int i = 0; i < 4 && nw < 3; i++) {
        w[nw] = gf_mul(Vb[i], Vb[i]) ^ gf_mul(d, Vb[i]);
        if (rank_bytes(w, nw + 1) == nw + 1) nw++;
    }
    if (nw != 3) { fprintf(stderr, "W_d is not 3-dimensional\n"); exit(1); }
    u[0] = gf_mul(gf_mul(w[1], w[2]), w[1] ^ w[2]);
    u[1] = gf_mul(gf_mul(w[0], w[2]), w[0] ^ w[2]);
    u[2] = gf_mul(gf_mul(w[0], w[1]), w[0] ^ w[1]);
}

static void UV(const uint8_t *Vb, uint8_t *u) {
    for (int i = 0; i < 4; i++) {
        uint8_t h[3], p = 1; int n = 0;
        for (int j = 0; j < 4; j++) if (j != i) h[n++] = Vb[j];
        for (int m = 1; m < 8; m++) { uint8_t x = 0; for (int t = 0; t < 3; t++) if (m >> t & 1) x ^= h[t]; p = gf_mul(p, x); }
        u[i] = p;
    }
}

// dim { lam : lam*u_i in Yperp for all i } , Yb = basis of Y
static int delta(const uint8_t *u, int nu, const uint8_t *Yb) {
    uint8_t rows[16]; int n = 0;
    for (int i = 0; i < nu; i++) for (int j = 0; j < 4; j++) rows[n++] = FV[gf_mul(u[i], Yb[j])];
    return 8 - rank_bytes(rows, n);
}

static void modeA(void) {
    uint8_t Vb[4] = {1, 2, 4, 8}, Vel[15], ud[15][3], uv[4];
    if (FILT) memcpy(Vb, FILTV, 4);          // adj A v1 v2 v3 v4 : same test on another point subspace
    elements(Vb, 4, Vel);
    for (int e = 0; e < 15; e++) Ud(Vb, Vel[e], ud[e]);
    UV(Vb, uv);
    long cnt[16] = {0}, mism = 0, n9 = 0, ninj = 0, nlow = 0, checks = 0; int best = 999;
    for (int s = 0; s < NSUB; s++) {
        const uint8_t *Yb = SUB[s]; uint8_t Yp[8]; Res R;
        if (perp(Yb, 4, Yp) != 4) { fprintf(stderr, "perp dim\n"); exit(1); }
        evalY(Vb, Vel, Yp, &R);
        int d0 = delta(uv, 4, Yb);
        checks++; if (R.dimS != 8 + d0) mism++;
        for (int e = 0; e < 15; e++) { int de = delta(ud[e], 3, Yb); checks++; if (R.rk[e] != 4 + d0 - de) mism++; }
        if (R.dimS > 8) { n9++; continue; }
        if (!R.valid) { ninj++; continue; }
        int n3 = 0;
        for (int e = 0; e < 15; e++) { if (R.rk[e] == 3) n3++; if (R.rk[e] < 3) nlow++; }
        cnt[n3]++;
        if (R.tot < best) best = R.tot;
    }
    printf("V = span(%d,%d,%d,%d): %d subspaces Y, %ld predictions checked, %ld mismatches\n", Vb[0], Vb[1], Vb[2], Vb[3], NSUB, checks, mism);
    printf("dim S_Y > 8: %ld   c not injective: %ld   helpers below 3 bits: %ld\n", n9, ninj, nlow);
    for (int i = 0; i < 16; i++) if (cnt[i]) printf("valid Y with %d three-bit helpers: %ld\n", i, cnt[i]);
    printf("best total: %d bits\n", best);
}

// dim S = 9 or 10: pick the best 8-dim D inside S on which c is injective (every complement of ker c is tried)
static int best_hyper(const Res *R, const uint8_t *Vel, uint32_t *Dout, int *rkout) {
    int ds = R->dimS, k = ds - 8, nk = 0, n = 0;
    uint32_t ker[4], bs[8]; uint8_t cv[9];
    for (int m = 1; m < (1 << ds); m++) {                  // elements of S with c = 0
        uint32_t x = 0; for (int i = 0; i < ds; i++) if (m >> i & 1) x ^= R->S[i];
        if ((x >> 16) == 0) {
            int indep = 1;                                   // keep a basis of the kernel
            for (int q = 1; q < (1 << nk); q++) { uint32_t y = 0; for (int i = 0; i < nk; i++) if (q >> i & 1) y ^= ker[i]; if (y == x) indep = 0; }
            if (indep && nk < k) ker[nk++] = x;
        }
    }
    if (nk != k) return -1;
    for (int i = 0; i < ds && n < 8; i++) { cv[n] = (uint8_t)(R->S[i] >> 16); if (rank_bytes(cv, n + 1) == n + 1) bs[n++] = R->S[i]; }
    if (n != 8) return -1;
    uint8_t val[8][15], kv[4][15];                          // kv[q] = value of kernel element number q (q = 0 is zero)
    uint32_t kel[4];
    for (int q = 0; q < (1 << k); q++) { kel[q] = 0; for (int i = 0; i < k; i++) if (q >> i & 1) kel[q] ^= ker[i]; }
    for (int e = 0; e < 15; e++) {
        uint8_t d = Vel[e], d3 = gf_mul(d, gf_mul(d, d));
        for (int i = 0; i < 8; i++) val[i][e] = gf_mul((uint8_t)(bs[i] & 255), d3) ^ gf_mul((uint8_t)((bs[i] >> 8) & 255), d) ^ (uint8_t)(bs[i] >> 16);
        for (int q = 0; q < (1 << k); q++) kv[q][e] = gf_mul((uint8_t)(kel[q] & 255), d3) ^ gf_mul((uint8_t)((kel[q] >> 8) & 255), d) ^ (uint8_t)(kel[q] >> 16);
    }
    int best = 999;
    for (long phi = 0; phi < (1L << (8 * k)); phi++) {
        int tot = 0, rk[15];
        for (int e = 0; e < 15; e++) {
            uint8_t t[8]; for (int i = 0; i < 8; i++) t[i] = val[i][e] ^ kv[(phi >> (k * i)) & ((1 << k) - 1)][e];
            rk[e] = rank_bytes(t, 8); tot += rk[e];
        }
        if (tot < best) { best = tot; for (int i = 0; i < 8; i++) Dout[i] = bs[i] ^ kel[(phi >> (k * i)) & ((1 << k) - 1)]; memcpy(rkout, rk, sizeof(rk)); }
    }
    return best;
}

static void modeB(void) {
    long hist[130] = {0}, nV = 0, n10 = 0, n10ok = 0, nUbad = 0, nev = 0;
    int best = 999; uint8_t bV[4], bYp[4]; uint32_t bS[8]; int brk[15];
    for (int s = 0; s < NSUB; s++) {
        const uint8_t *Vb = SUB[s]; uint8_t t[5] = {0}, Vel[15];
        if (FILT) { uint8_t t8[8]; memcpy(t8, Vb, 4); memcpy(t8 + 4, FILTV, 4); if (rank_bytes(t8, 8) != 4) continue; }
        else { memcpy(t, Vb, 4); t[4] = 1; if (rank_bytes(t, 5) != 4) continue; }
        nV++;
        elements(Vb, 4, Vel);
        int vbest = 999;
        for (int e = 0; e <= 15; e++) {
            uint8_t u[4], inU[256] = {0}, seen[256] = {0}, Uel[7];
            int nz;
            if (e < 15) {
                Ud(Vb, Vel[e], u);
                if (rank_bytes(u, 3) != 3) { nUbad++; continue; }
                elements(u, 3, Uel); inU[0] = 1; for (int i = 0; i < 7; i++) inU[Uel[i]] = 1;
                nz = 256;
            } else {
                UV(Vb, u);
                if (rank_bytes(u, 4) != 4) { nUbad++; continue; }
                nz = 2;                       // single candidate Yperp = Uhat_V
            }
            for (int z = 1; z < nz; z++) {
                uint8_t Yp[4]; Res R; uint32_t D[8]; int rk[15], tot;
                if (e < 15) {
                    if (inU[z] || seen[z]) continue;
                    seen[z] = 1; for (int i = 0; i < 7; i++) seen[z ^ Uel[i]] = 1;
                    Yp[0] = u[0]; Yp[1] = u[1]; Yp[2] = u[2]; Yp[3] = (uint8_t)z;
                } else memcpy(Yp, u, 4);
                evalY(Vb, Vel, Yp, &R); nev++;
                if (!R.valid) continue;
                if (R.dimS == 8) { tot = R.tot; memcpy(D, R.S, sizeof(D)); memcpy(rk, R.rk, sizeof(rk)); }
                else if (R.dimS <= 10) { tot = best_hyper(&R, Vel, D, rk); if (tot < 0) continue; if (R.dimS == 10) n10ok++; }
                else { n10++; continue; }
                if (tot < vbest) vbest = tot;
                if (tot < best) { best = tot; memcpy(bV, Vb, 4); memcpy(bYp, Yp, 4); memcpy(bS, D, sizeof(bS)); memcpy(brk, rk, sizeof(brk)); }
            }
        }
        if (vbest < 130) hist[vbest]++;
    }
    printf("point subspaces V containing 1: %ld   exact evaluations: %ld   dim S = 10 handled: %ld   dim S > 10 skipped: %ld   degenerate U: %ld\n", nV, nev, n10ok, n10, nUbad);
    for (int i = 0; i < 130; i++) if (hist[i]) printf("best total %d bits: %ld subspaces V\n", i, hist[i]);
    printf("BEST %d  V=%d,%d,%d,%d  Yperp=%d,%d,%d,%d\n", best, bV[0], bV[1], bV[2], bV[3], bYp[0], bYp[1], bYp[2], bYp[3]);
    uint8_t Vel[15]; elements(bV, 4, Vel);
    printf("offsets:"); for (int e = 0; e < 15; e++) printf(" %d", Vel[e]); printf("\nranks:  ");
    for (int e = 0; e < 15; e++) printf(" %d", brk[e]);
    printf("\nabc:");
    for (int i = 0; i < 8; i++) printf(" %u,%u,%u", bS[i] & 255, (bS[i] >> 8) & 255, bS[i] >> 16);
    printf("\n");
}

int main(int argc, char **argv) {
    init_all(); enum_sub4();
    if (NSUB != 200787) { fprintf(stderr, "subspace count %d\n", NSUB); return 1; }
    if (argc == 6) { FILT = 1; for (int i = 0; i < 4; i++) FILTV[i] = (uint8_t)atoi(argv[2 + i]); }
    if (argc > 1 && argv[1][0] == 'B') modeB(); else modeA();
    return 0;
}
