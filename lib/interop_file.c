/* Interop test against stripes written by another library.  Each file stripe_<n>_<k>.bin holds the n shards of one
 * stripe, back to back, exactly as that library encoded them.  For every shard: pretend it is lost, build the
 * helper payloads from the other shards, rebuild, and compare with the library's own bytes.
 * usage: ./interop_file <dir> <plain|cauchy> <label> */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "rsrepair.h"

int main(int argc, char **argv)
{
    if (argc < 4) return 2;
    int flavour = strcmp(argv[2], "cauchy") ? RSR_PLAIN : RSR_CAUCHY;
    long rebuilt = 0, bytes = 0; int codes = 0;
    for (int n = 4; n <= 20; n++) for (int k = 1; k < n; k++) {
        if (!rsr_supported(n, k)) continue;
        char path[512];
        snprintf(path, sizeof path, "%s/stripe_%d_%d.bin", argv[1], n, k);
        FILE *f = fopen(path, "rb");
        if (!f) { printf("missing %s\n", path); return 1; }
        fseek(f, 0, SEEK_END); long sz = ftell(f); fseek(f, 0, SEEK_SET);
        if (sz <= 0 || sz % n) { printf("bad size %s\n", path); return 1; }
        size_t len = (size_t)sz / (size_t)n;
        uint8_t *all = malloc((size_t)sz);
        if (fread(all, 1, (size_t)sz, f) != (size_t)sz) return 1;
        fclose(f);
        codes++;
        for (int failed = 0; failed < n; failed++) {
            rsr_plan plan;
            if (rsr_plan_init(&plan, n, k, failed, flavour)) return 1;
            uint8_t *pay[RSR_MAX_N] = {0}, *out = malloc(len);
            for (int i = 0; i < n; i++) if (i != failed) { pay[i] = malloc(rsr_payload_len(&plan, i, len)); rsr_helper(&plan, i, all + (size_t)i * len, len, pay[i]); }
            rsr_rebuild(&plan, (const uint8_t *const *)pay, len, out);
            if (memcmp(out, all + (size_t)failed * len, len)) { printf("MISMATCH %s shard %d\n", path, failed); return 1; }
            rebuilt++; bytes += (long)len;
            for (int i = 0; i < n; i++) free(pay[i]);
            free(out);
        }
        free(all);
    }
    printf("OK: %s: %d codes, %ld lost shards rebuilt byte-identical (%ld bytes)\n", argv[3], codes, rebuilt, bytes);
    return 0;
}
