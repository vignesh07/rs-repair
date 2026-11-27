/* Chooses, for one failed node, the change of coordinates on the lost byte that makes the rebuild XOR network smallest.
 *
 * The rebuild is  lost = sum_i D_i * y_i  over GF(2), y_i = the r_i bits helper i sends per byte.  We may send any basis
 * of the helper's bits and compute any invertible image A*lost first, so the cost of helper i is the weight of the
 * lightest basis of A*V_i (V_i = column space of D_i), and we minimise the total over A by hill climbing with restarts.
 *
 * stdin:  N, then per node: H, then H lines "r d1 .. dr" (d = columns of D_i as bytes)
 * stdout: per node: "cost a0 .. a7"   (a_b = row b of A as a mask over the bits of the lost byte)                       */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXH 40
static int H, R[MAXH], nspan[MAXH];
static uint8_t span[MAXH][255];
static uint64_t rng_state = 0x9E3779B97F4A7C15ULL;
static uint64_t rnd(void) { rng_state ^= rng_state << 13; rng_state ^= rng_state >> 7; rng_state ^= rng_state << 17; return rng_state; }

static int cost_of(const uint8_t *T)
{
    int total = 0;
    for (int h = 0; h < H; h++) {
        uint8_t byw[9][255]; int cnt[9] = {0};
        for (int i = 0; i < nspan[h]; i++) { uint8_t v = T[span[h][i]]; int w = __builtin_popcount(v); byw[w][cnt[w]++] = v; }
        uint8_t piv[8] = {0}; int got = 0;
        for (int w = 1; w <= 8 && got < R[h]; w++)
            for (int i = 0; i < cnt[w] && got < R[h]; i++) {
                uint8_t v = byw[w][i];
                for (int b = 7; b >= 0 && v; b--) if (v >> b & 1) { if (piv[b]) v ^= piv[b]; else { piv[b] = v; got++; total += w; v = 0; } }
            }
    }
    return total;
}

static void table_of(const uint8_t *A, uint8_t *T)
{
    for (int x = 0; x < 256; x++) { uint8_t z = 0; for (int b = 0; b < 8; b++) z |= (uint8_t)((__builtin_popcount(A[b] & x) & 1) << b); T[x] = z; }
}

int main(void)
{
    int N;
    if (scanf("%d", &N) != 1) return 1;
    for (int node = 0; node < N; node++) {
        if (scanf("%d", &H) != 1) return 1;
        for (int h = 0; h < H; h++) {
            int d[8];
            if (scanf("%d", &R[h]) != 1) return 1;
            for (int j = 0; j < R[h]; j++) if (scanf("%d", &d[j]) != 1) return 1;
            nspan[h] = 0;
            for (int m = 1; m < (1 << R[h]); m++) { uint8_t v = 0; for (int j = 0; j < R[h]; j++) if (m >> j & 1) v ^= (uint8_t)d[j]; span[h][nspan[h]++] = v; }
        }
        uint8_t bestA[8]; int best = 1 << 30;
        for (int restart = 0; restart < 60; restart++) {
            uint8_t A[8], T[256];
            for (int b = 0; b < 8; b++) A[b] = (uint8_t)(1 << b);
            if (restart) for (int s = 0; s < 64; s++) { int a = rnd() % 8, b = rnd() % 8; if (a != b) A[b] ^= A[a]; }
            table_of(A, T);
            int cur = cost_of(T);
            for (;;) {
                int ba = -1, bb = -1, bc = cur;
                for (int a = 0; a < 8; a++) for (int b = 0; b < 8; b++) {
                    if (a == b) continue;
                    uint8_t T2[256];
                    for (int x = 0; x < 256; x++) T2[x] = T[x] ^ (uint8_t)(((T[x] >> a) & 1) << b);
                    int c = cost_of(T2);
                    if (c < bc) { bc = c; ba = a; bb = b; }
                }
                if (ba < 0) break;
                A[bb] ^= A[ba]; table_of(A, T); cur = bc;
            }
            if (cur < best) { best = cur; memcpy(bestA, A, 8); }
        }
        printf("%d", best);
        for (int b = 0; b < 8; b++) printf(" %d", bestA[b]);
        printf("\n");
    }
    return 0;
}
