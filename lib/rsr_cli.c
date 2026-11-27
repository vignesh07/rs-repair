/* Command-line front end for rsrepair, working on shard files.
 *
 *   rsr info
 *   rsr helper  <plain|cauchy> <n> <k> <lost> <i> <shard file> <payload file>
 *   rsr rebuild <plain|cauchy> <n> <k> <lost> <shard bytes> <out file> <payload 0> ... <payload n-1>   ("-" for the lost one)
 *
 * Files are streamed in 1 MiB pieces, so shards of any size work in a few MiB of memory. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "rsrepair.h"

#define PIECE (1u << 20) /* a multiple of the 512-byte block */

static int flavour_of(const char *s) { return strcmp(s, "cauchy") ? RSR_PLAIN : RSR_CAUCHY; }

static int info(void)
{
    printf("bits downloaded per lost byte\n");
    printf("%-10s %7s   %-28s %s\n", "code", "normal", "rsrepair: best / avg / worst", "saving (avg)");
    for (int r = 2; r <= 4; r++)
        for (int n = 4; n <= 20; n++) {
            int k = n - r, lo = 1 << 30, hi = 0, sum = 0;
            char name[16], range[40];
            if (!rsr_supported(n, k)) continue;
            for (int f = 0; f < n; f++) {
                rsr_plan p;
                if (rsr_plan_init(&p, n, k, f, RSR_PLAIN)) return 1;
                if (p.total_bits < lo) lo = p.total_bits;
                if (p.total_bits > hi) hi = p.total_bits;
                sum += p.total_bits;
            }
            snprintf(name, sizeof name, "RS(%d,%d)", n, k);
            snprintf(range, sizeof range, "%d / %.2f / %d", lo, (double)sum / n, hi);
            printf("%-10s %7d   %-28s %.1f%%\n", name, 8 * k, range, 100.0 * (1.0 - (double)sum / n / (8.0 * k)));
        }
    return 0;
}

int main(int argc, char **argv)
{
    if (argc >= 2 && !strcmp(argv[1], "info")) return info();
    if (argc == 9 && !strcmp(argv[1], "helper")) {
        rsr_plan p;
        int n = atoi(argv[3]), k = atoi(argv[4]), lost = atoi(argv[5]), i = atoi(argv[6]);
        if (rsr_plan_init(&p, n, k, lost, flavour_of(argv[2])) || i < 0 || i >= n || i == lost) { fprintf(stderr, "unsupported code or chunk number\n"); return 1; }
        FILE *in = fopen(argv[7], "rb"), *out = fopen(argv[8], "wb");
        if (!in || !out) { perror("open"); return 1; }
        uint8_t *buf = malloc(PIECE), *pay = malloc(PIECE);
        size_t got;
        while ((got = fread(buf, 1, PIECE, in)) > 0) {
            size_t pl = rsr_helper(&p, i, buf, got, pay);
            if (fwrite(pay, 1, pl, out) != pl) { perror("write"); return 1; }
        }
        fclose(in);
        if (fclose(out)) { perror("close"); return 1; }
        return 0;
    }
    if (argc >= 8 && !strcmp(argv[1], "rebuild")) {
        rsr_plan p;
        int n = atoi(argv[3]), k = atoi(argv[4]), lost = atoi(argv[5]);
        unsigned long long total = strtoull(argv[6], NULL, 10);
        if (rsr_plan_init(&p, n, k, lost, flavour_of(argv[2])) || argc != 8 + n) { fprintf(stderr, "unsupported code, or expected %d payload files\n", n); return 1; }
        FILE *out = fopen(argv[7], "wb"), *in[RSR_MAX_N] = {0};
        uint8_t *pay[RSR_MAX_N] = {0}, *buf = malloc(PIECE);
        if (!out) { perror("open"); return 1; }
        for (int i = 0; i < n; i++) {
            if (i == lost) continue;
            in[i] = fopen(argv[8 + i], "rb");
            if (!in[i]) { perror(argv[8 + i]); return 1; }
            pay[i] = malloc(PIECE);
        }
        while (total) {
            size_t len = total < PIECE ? (size_t)total : PIECE;
            for (int i = 0; i < n; i++) {
                if (i == lost) continue;
                size_t pl = rsr_payload_len(&p, i, len);
                if (fread(pay[i], 1, pl, in[i]) != pl) { fprintf(stderr, "payload %d is too short\n", i); return 1; }
            }
            rsr_rebuild(&p, (const uint8_t *const *)pay, len, buf);
            if (fwrite(buf, 1, len, out) != len) { perror("write"); return 1; }
            total -= len;
        }
        if (fclose(out)) { perror("close"); return 1; }
        return 0;
    }
    fprintf(stderr, "usage:\n  rsr info\n  rsr helper  <plain|cauchy> n k lost i shard payload\n  rsr rebuild <plain|cauchy> n k lost shard_bytes out payload_0 ... payload_{n-1}   (\"-\" for the lost one)\n");
    return 2;
}
