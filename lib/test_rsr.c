/* Self-contained correctness test: every supported code, both flavours, every failed node, many lengths, every code
 * path.  The stripes are produced here by a plain reference encoder (interop with the real libraries is tested
 * separately).  Buffers are allocated at their exact size so a sanitizer build catches any out-of-bounds access. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "rsrepair.h"

static uint8_t gexp[512], glog[256];
static void gf_init(void) { int x = 1; for (int i = 0; i < 255; i++) { gexp[i] = (uint8_t)x; glog[x] = (uint8_t)i; x <<= 1; if (x & 0x100) x ^= 0x11D; } for (int i = 255; i < 512; i++) gexp[i] = gexp[i - 255]; }
static uint8_t mul(uint8_t a, uint8_t b) { return (a && b) ? gexp[glog[a] + glog[b]] : 0; }
static uint8_t inv(uint8_t a) { return gexp[255 - glog[a]]; }

static uint64_t s = 88172645463325252ULL;
static uint8_t rnd(void) { s ^= s << 13; s ^= s >> 7; s ^= s << 17; return (uint8_t)(s >> 24); }

#define LMAX 4351
static const size_t lens[] = {0, 1, 7, 8, 9, 63, 64, 65, 511, 512, 513, 1024, 1553, LMAX};

int main(void)
{
    gf_init();
    static uint8_t chunk[RSR_MAX_N][LMAX], coef[RSR_MAX_N][LMAX];
    long repairs = 0, plans = 0;
    int nb = rsr_set_backend(RSR_BACKEND_AUTO);
    printf("best backend on this CPU: %s\n", rsr_backend_name());
    for (int n = 4; n <= 20; n++) for (int k = 1; k < n; k++) {
        if (!rsr_supported(n, k)) continue;
        for (int flavour = 0; flavour < 2; flavour++) {
            for (int j = 0; j < k; j++) for (int t = 0; t < LMAX; t++) coef[j][t] = rnd();
            for (int i = 0; i < n; i++) for (int t = 0; t < LMAX; t++) {
                uint8_t acc = 0;
                if (flavour == RSR_PLAIN) { uint8_t pw = 1; for (int j = 0; j < k; j++) { acc ^= mul(coef[j][t], pw); pw = mul(pw, (uint8_t)i); } }
                else if (i < k) acc = coef[i][t];
                else for (int j = 0; j < k; j++) acc ^= mul(coef[j][t], inv((uint8_t)(i ^ j)));
                chunk[i][t] = acc;
            }
            for (int failed = 0; failed < n; failed++) {
                rsr_plan plan;
                if (rsr_plan_init(&plan, n, k, failed, flavour)) { printf("plan init failed\n"); return 1; }
                plans++;
                for (size_t li = 0; li < sizeof lens / sizeof *lens; li++) {
                    size_t len = lens[li];
                    uint8_t *ref[RSR_MAX_N] = {0};
                    for (int be = 1; be <= nb; be++) {
                        rsr_set_backend(be);
                        uint8_t *pay[RSR_MAX_N] = {0};
                        for (int i = 0; i < n; i++) {
                            if (i == failed) continue;
                            size_t pl = rsr_payload_len(&plan, i, len);
                            uint8_t *src = malloc(len ? len : 1);           /* exact-size copies: no slack to hide overruns */
                            memcpy(src, chunk[i], len);
                            pay[i] = malloc(pl ? pl : 1);
                            if (rsr_helper(&plan, i, src, len, pay[i]) != pl) { printf("payload length mismatch\n"); return 1; }
                            free(src);
                            if (be == 1) { ref[i] = malloc(pl ? pl : 1); memcpy(ref[i], pay[i], pl); }
                            else if (memcmp(ref[i], pay[i], pl)) { printf("payload differs between code paths: RS(%d,%d) helper %d len %zu backend %d\n", n, k, i, len, be); return 1; }
                        }
                        uint8_t *out = malloc(len ? len : 1);
                        rsr_rebuild(&plan, (const uint8_t *const *)pay, len, out);
                        if (memcmp(out, chunk[failed], len)) { printf("WRONG REBUILD: RS(%d,%d) flavour %d failed %d len %zu backend %d\n", n, k, flavour, failed, len, be); return 1; }
                        repairs++;
                        free(out);
                        for (int i = 0; i < n; i++) free(pay[i]);
                    }
                    for (int i = 0; i < n; i++) free(ref[i]);
                    if (len == LMAX && rsr_total_payload_len(&plan, 512 * 8) * 8 != (size_t)plan.total_bits * 512 * 8) { printf("payload accounting wrong\n"); return 1; }
                }
            }
        }
    }
    printf("OK: %ld plans (code x flavour x failed node), %ld rebuilds byte-identical, payloads identical across %d code paths\n", plans, repairs, nb);
    return 0;
}
