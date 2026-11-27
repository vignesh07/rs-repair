/* Speed and traffic of rsrepair against Intel ISA-L's own single-chunk reconstruction, same machine, one pinned thread.
 * Stripes are encoded by ISA-L with its Cauchy matrix.  Every figure is the best time seen over three interleaved
 * rounds (ISA-L, then rsrepair, then ISA-L again ...), so a slow spell on a shared machine cannot favour one side.
 * usage: ./bench n k [n k ...]   -> JSON on stdout */
#define _GNU_SOURCE
#include <sched.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "erasure_code.h"
#include "rsrepair.h"

#define ROUNDS 3
static double now(void) { struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts); return ts.tv_sec + ts.tv_nsec * 1e-9; }
/* keep the best time of stmt in *slot: at least 3 runs and at least `window` seconds */
#define BEST(slot, window, stmt) do { double t0_ = now(); int it_ = 0; \
        while (it_ < 3 || now() - t0_ < (window)) { double a_ = now(); stmt; double d_ = now() - a_; if (d_ < (slot)) (slot) = d_; it_++; } } while (0)

static void *xalloc(size_t n) { void *p = NULL; if (posix_memalign(&p, 64, n ? n : 64)) exit(1); memset(p, 0, n ? n : 64); return p; }

int main(int argc, char **argv)
{
    cpu_set_t set; CPU_ZERO(&set); CPU_SET(0, &set); sched_setaffinity(0, sizeof set, &set);
    static const size_t sizes[] = {4 << 10, 64 << 10, 1 << 20, 16 << 20};
    static const char *bname[] = {"", "portable", "avx2", "avx512+gfni"};
    const int best = rsr_set_backend(0);
    printf("[\n");
    for (int ai = 1; ai + 1 < argc; ai += 2) {
        int n = atoi(argv[ai]), k = atoi(argv[ai + 1]);
        if (!rsr_supported(n, k)) { fprintf(stderr, "RS(%d,%d) unsupported\n", n, k); continue; }
        for (size_t si = 0; si < sizeof sizes / sizeof *sizes; si++) {
            size_t len = sizes[si];
            uint8_t *chunk[32], *out = xalloc(len), *pay[32] = {0};
            static uint8_t a[32 * 32], tbls[32 * 32 * 32], tb_dec[32 * 32], tb_base[32 * 32], b[32 * 32], binv[32 * 32];
            for (int i = 0; i < n; i++) chunk[i] = xalloc(len);
            uint64_t s = 0x9E3779B97F4A7C15ULL + len;
            for (int i = 0; i < k; i++) for (size_t t = 0; t < len; t += 8) { s ^= s << 13; s ^= s >> 7; s ^= s << 17; memcpy(chunk[i] + t, &s, 8); }
            gf_gen_cauchy1_matrix(a, n, k);
            ec_init_tables(k, n - k, &a[k * k], tbls);
            ec_encode_data((int)len, k, n - k, tbls, chunk, &chunk[k]);

            /* normal repair with ISA-L: lost chunk 0 from the next k chunks (row 0 of the inverted matrix) */
            uint8_t *srcs[32];
            for (int i = 0; i < k; i++) { srcs[i] = chunk[i + 1]; memcpy(&b[k * i], &a[k * (i + 1)], k); }
            if (gf_invert_matrix(b, binv, k)) return 1;
            ec_init_tables(k, 1, binv, tb_dec);          /* table format of the kernel ISA-L picks on this CPU */
            ec_init_tables_base(k, 1, binv, tb_base);    /* classic format, for its AVX2 kernel */

            double t_isal = 1e30, t_isal_avx2 = 1e30, reb[4][32], hlp[4][32];
            for (int be = 0; be < 4; be++) for (int i = 0; i < 32; i++) reb[be][i] = hlp[be][i] = 1e30;
            size_t traffic[32] = {0};
            for (int round = 0; round < ROUNDS; round++) {
                BEST(t_isal, 0.05, ec_encode_data((int)len, k, 1, tb_dec, srcs, &out));
                if (memcmp(out, chunk[0], len)) { fprintf(stderr, "ISA-L decode wrong?\n"); return 1; }
                BEST(t_isal_avx2, 0.05, ec_encode_data_avx2((int)len, k, 1, tb_base, srcs, &out));
                if (memcmp(out, chunk[0], len)) { fprintf(stderr, "ISA-L avx2 decode wrong?\n"); return 1; }
                for (int be = best; be >= 1; be--) {
                    rsr_set_backend(be);
                    for (int failed = 0; failed < (be == 1 ? 1 : n); failed++) {      /* portable path: one lost chunk is enough */
                        rsr_plan plan;
                        if (rsr_plan_init(&plan, n, k, failed, RSR_CAUCHY)) return 1;
                        traffic[failed] = rsr_total_payload_len(&plan, len);
                        for (int i = 0; i < n; i++) if (i != failed) {
                            pay[i] = xalloc(rsr_payload_len(&plan, i, len));
                            if (failed == 0) BEST(hlp[be][i], 0.01, rsr_helper(&plan, i, chunk[i], len, pay[i]));
                            else rsr_helper(&plan, i, chunk[i], len, pay[i]);
                        }
                        BEST(reb[be][failed], 0.05, rsr_rebuild(&plan, (const uint8_t *const *)pay, len, out));
                        if (memcmp(out, chunk[failed], len)) { fprintf(stderr, "WRONG rebuild\n"); return 1; }
                        for (int i = 0; i < n; i++) { free(pay[i]); pay[i] = NULL; }
                    }
                }
            }
            double mb = len / 1048576.0, tmin = 1e30, tmax = 0, tsum = 0;
            for (int f = 0; f < n; f++) { if (traffic[f] < tmin) tmin = (double)traffic[f]; if (traffic[f] > tmax) tmax = (double)traffic[f]; tsum += (double)traffic[f]; }
            printf("  {\"n\": %d, \"k\": %d, \"chunk_bytes\": %zu, \"isal_MBps\": %.0f, \"isal_avx2_MBps\": %.0f,\n", n, k, len, mb / t_isal, mb / t_isal_avx2);
            for (int be = best; be >= 1; be--) {
                int cnt = be == 1 ? 1 : n;
                double lo = 1e30, hi = 0, sum = 0, hsum = 0;
                for (int f = 0; f < cnt; f++) { if (reb[be][f] < lo) lo = reb[be][f]; if (reb[be][f] > hi) hi = reb[be][f]; sum += reb[be][f]; }
                for (int i = 1; i < n; i++) hsum += hlp[be][i];
                printf("   \"%s\": {\"rebuild_MBps_slowest\": %.0f, \"rebuild_MBps_avg\": %.0f, \"rebuild_MBps_fastest\": %.0f, \"helper_MBps\": %.0f},\n",
                       bname[be], mb / hi, mb / (sum / cnt), mb / lo, mb / (hsum / (n - 1)));
            }
            printf("   \"download_vs_normal\": {\"best\": %.4f, \"avg\": %.4f, \"worst\": %.4f}}%s\n",
                   tmin / ((double)k * len), tsum / n / ((double)k * len), tmax / ((double)k * len),
                   (ai + 3 < argc || si + 1 < sizeof sizes / sizeof *sizes) ? "," : "");
            fflush(stdout);
            for (int i = 0; i < n; i++) free(chunk[i]);
            free(out);
        }
    }
    printf("]\n");
    return 0;
}
