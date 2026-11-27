/* rsrepair - low-traffic repair of one lost chunk of a Reed-Solomon stripe.
 *
 * Normal repair downloads k whole chunks.  Here every surviving node ("helper") sends only a few bits per byte of its
 * chunk, and the lost chunk is rebuilt from those bits.  The stripe itself is ordinary Reed-Solomon over GF(256) as
 * written by existing libraries - nothing is re-encoded.
 *
 * Supported codes: RS(n,k) with n-k in {2,3,4} and n <= 16, plus RS(20,17), in two flavours
 *   RSR_PLAIN   chunk i holds f(i) for a polynomial f of degree < k
 *               (Backblaze JavaReedSolomon, klauspost/reedsolomon default matrix)
 *   RSR_CAUCHY  systematic code with parity rows 1/(i xor j)
 *               (Intel ISA-L gf_gen_cauchy1_matrix, Hadoop HDFS RS coder)
 * both over GF(256) with polynomial 0x11D.  Chunks are numbered 0..n-1 in the library's own order (data, then parity).
 *
 * Typical use
 *   rsr_plan plan;
 *   rsr_plan_init(&plan, 14, 10, lost, RSR_CAUCHY);
 *   on every surviving node i:  rsr_helper(&plan, i, chunk_i, len, payload_i);     // send payload_i
 *   on the node rebuilding:     rsr_rebuild(&plan, payloads, len, out);            // payloads[lost] is ignored
 *
 * Wire format of a payload (identical on every CPU):  the chunk is cut into blocks of 512 bytes (the last may be
 * shorter).  For each block of m bytes the helper sends r "planes" of ceil(m/8) bytes; bit t of plane j (byte t/8,
 * bit t%8) is the j-th parity bit computed from chunk byte t.  To stream, cut the chunk at multiples of 512 bytes.
 */
#ifndef RSREPAIR_H
#define RSREPAIR_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define RSR_MAX_N 32

enum { RSR_PLAIN = 0, RSR_CAUCHY = 1 };
enum { RSR_BACKEND_AUTO = 0, RSR_BACKEND_PORTABLE = 1, RSR_BACKEND_AVX2 = 2, RSR_BACKEND_AVX512 = 3 };

typedef struct rsr_plan {
    int n, k, failed;
    int total_bits;                 /* bits downloaded per lost byte (normal repair: 8*k) */
    uint8_t bits[RSR_MAX_N];        /* bits helper i sends per chunk byte (0 for the failed node) */
    uint8_t hmask[RSR_MAX_N][8];    /* bit j sent by helper i = parity(hmask[i][j] & chunk byte) */
    uint8_t dmask[RSR_MAX_N][8];    /* that bit feeds the intermediate bits selected by dmask[i][j] */
    uint8_t fin[8];                 /* lost byte = XOR of fin[b] over the set intermediate bits b */
    uint64_t h_sel[RSR_MAX_N];      /* hmask packed for the vector code */
    uint64_t fin_mat;               /* fin as a GF(2) matrix for the vector code */
    uint8_t fin_lut[256];
    uint8_t h_lut[RSR_MAX_N][256];  /* chunk byte -> the r bits to send */
    /* rebuild network: intermediate plane b = XOR of payload planes term[term_start[b] .. term_start[b+1]) */
    uint16_t term_start[9];
    struct { uint8_t helper, plane; } term[64 * (RSR_MAX_N - 1)];
} rsr_plan;

/* 0 on success, -1 if the code or node is not supported. */
int rsr_plan_init(rsr_plan *plan, int n, int k, int failed, int flavour);
/* Same for a code whose chunk i holds mult[i] * f(i) (generalised Reed-Solomon); mult = n nonzero bytes. */
int rsr_plan_init_grs(rsr_plan *plan, int n, int k, int failed, const uint8_t *mult);
/* 1 if RS(n,k) is in the built-in table. */
int rsr_supported(int n, int k);

/* Size of the payload helper i sends for a chunk of chunk_len bytes, and the sum over all helpers. */
size_t rsr_payload_len(const rsr_plan *plan, int helper, size_t chunk_len);
size_t rsr_total_payload_len(const rsr_plan *plan, size_t chunk_len);

/* Helper side: squeeze a chunk into its payload. Returns the payload length. */
size_t rsr_helper(const rsr_plan *plan, int helper, const uint8_t *chunk, size_t chunk_len, uint8_t *payload);
/* Rebuild side: payloads[i] for every i != failed, each produced from chunk_len bytes. Writes chunk_len bytes. */
void rsr_rebuild(const rsr_plan *plan, const uint8_t *const *payloads, size_t chunk_len, uint8_t *out);

/* Pick the code path (default: best the CPU supports). Returns the backend now in use. */
int rsr_set_backend(int backend);
const char *rsr_backend_name(void);

#ifdef __cplusplus
}
#endif
#endif
