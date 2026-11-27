/* Interop test against Intel ISA-L itself: ISA-L builds the Cauchy matrix and encodes the stripe; one chunk is
 * dropped; rsrepair rebuilds it from the helper payloads; the result must equal ISA-L's chunk byte for byte. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "erasure_code.h"
#include "rsrepair.h"

int main(void)
{
    const size_t len = 8192 + 513;          /* not a multiple of the block size on purpose */
    static uint8_t a[32 * 32], tbls[32 * 32 * 32];
    uint8_t *chunk[32];
    for (int i = 0; i < 32; i++) chunk[i] = malloc(len);
    long rebuilt = 0; int codes = 0;
    srand(12345);
    for (int n = 4; n <= 20; n++) for (int k = 1; k < n; k++) {
        if (!rsr_supported(n, k)) continue;
        codes++;
        gf_gen_cauchy1_matrix(a, n, k);
        ec_init_tables(k, n - k, &a[k * k], tbls);
        for (int i = 0; i < k; i++) for (size_t t = 0; t < len; t++) chunk[i][t] = (uint8_t)rand();
        ec_encode_data((int)len, k, n - k, tbls, chunk, &chunk[k]);
        for (int failed = 0; failed < n; failed++) {
            rsr_plan plan;
            if (rsr_plan_init(&plan, n, k, failed, RSR_CAUCHY)) return 1;
            uint8_t *pay[32] = {0}, *out = malloc(len);
            for (int i = 0; i < n; i++) if (i != failed) { pay[i] = malloc(rsr_payload_len(&plan, i, len)); rsr_helper(&plan, i, chunk[i], len, pay[i]); }
            rsr_rebuild(&plan, (const uint8_t *const *)pay, len, out);
            if (memcmp(out, chunk[failed], len)) { printf("MISMATCH RS(%d,%d) chunk %d\n", n, k, failed); return 1; }
            rebuilt++;
            for (int i = 0; i < n; i++) free(pay[i]);
            free(out);
        }
    }
    printf("OK: %d codes encoded by ISA-L (gf_gen_cauchy1_matrix + ec_encode_data), %ld lost chunks rebuilt byte-identical (%s)\n", codes, rebuilt, rsr_backend_name());
    return 0;
}
