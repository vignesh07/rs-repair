/* rsrepair - low-traffic repair of one lost chunk of a Reed-Solomon stripe. See rsrepair.h. */
#include "rsrepair.h"

#include <string.h>

#include "rsr_tables.h"

#define SB 512                      /* chunk bytes per block */
#define LSB8 0x0101010101010101ULL
#define GATHER 0x0102040810204080ULL /* (x & LSB8) * GATHER >> 56 collects bit 0 of byte t into bit t */
#define MAX_TERMS (64 * (RSR_MAX_N - 1))

/* ------------------------------------------------------------------ GF(256), polynomial 0x11D */

static uint8_t gf_mul(uint8_t a, uint8_t b)
{
    uint8_t r = 0;
    while (b) {
        if (b & 1) r ^= a;
        a = (uint8_t)((a << 1) ^ ((a & 0x80) ? 0x1D : 0));
        b >>= 1;
    }
    return r;
}

static uint8_t gf_inv(uint8_t a)
{
    uint8_t r = 1;
    for (int i = 0; i < 254; i++) r = gf_mul(r, a);
    return r;
}

static unsigned parity8(unsigned x)
{
    x ^= x >> 4;
    x ^= x >> 2;
    x ^= x >> 1;
    return x & 1;
}

/* ------------------------------------------------------------------ plans */

int rsr_supported(int n, int k)
{
    for (int i = 0; i < RSR_NCODES; i++)
        if (rsr_codes[i].n == n && rsr_codes[i].k == k) return 1;
    return 0;
}

int rsr_plan_init_grs(rsr_plan *p, int n, int k, int failed, const uint8_t *mult)
{
    int ci = -1;
    for (int i = 0; i < RSR_NCODES; i++)
        if (rsr_codes[i].n == n && rsr_codes[i].k == k) ci = i;
    if (ci < 0 || failed < 0 || failed >= n || n > RSR_MAX_N) return -1;
    if (mult)
        for (int i = 0; i < n; i++)
            if (!mult[i]) return -1;

    memset(p, 0, sizeof *p);
    p->n = n;
    p->k = k;
    p->failed = failed;
    const uint8_t *rec = rsr_blob + rsr_node_off[rsr_codes[ci].first + failed];
    uint8_t us = mult ? mult[failed] : 1;
    for (int b = 0; b < 8; b++) p->fin[b] = gf_mul(us, rec[b]);
    rec += 8;
    for (int i = 0; i < n; i++) {
        if (i == failed) {
            rec++; /* 0xFF marker */
            continue;
        }
        int r = *rec++;
        uint8_t uinv = mult ? gf_inv(mult[i]) : 1;
        p->bits[i] = (uint8_t)r;
        p->total_bits += r;
        for (int j = 0; j < r; j++) {
            /* the table's mask acts on f(i) = chunk byte / mult[i]; fold the division into the mask */
            uint8_t hm = 0;
            for (int e = 0; e < 8; e++) hm |= (uint8_t)(parity8(rec[j] & gf_mul(uinv, (uint8_t)(1 << e))) << e);
            p->hmask[i][j] = hm;
            p->dmask[i][j] = rec[r + j];
            p->h_sel[i] |= (uint64_t)hm << (8 * j);
        }
        rec += 2 * r;
        for (int c = 0; c < 256; c++) {
            uint8_t v = 0;
            for (int j = 0; j < r; j++) v |= (uint8_t)(parity8(p->hmask[i][j] & (unsigned)c) << j);
            p->h_lut[i][c] = v;
        }
    }
    for (int z = 0; z < 256; z++) {
        uint8_t v = 0;
        for (int b = 0; b < 8; b++)
            if (z >> b & 1) v ^= p->fin[b];
        p->fin_lut[z] = v;
    }
    for (int e = 0; e < 8; e++) {
        uint64_t row = 0;
        for (int b = 0; b < 8; b++) row |= (uint64_t)((p->fin[b] >> e) & 1) << b;
        p->fin_mat |= row << (8 * (7 - e));
    }
    int T = 0;
    for (int b = 0; b < 8; b++) {
        p->term_start[b] = (uint16_t)T;
        for (int i = 0; i < n; i++)
            for (int j = 0; j < p->bits[i]; j++)
                if (p->dmask[i][j] >> b & 1) {
                    p->term[T].helper = (uint8_t)i;
                    p->term[T].plane = (uint8_t)j;
                    T++;
                }
    }
    p->term_start[8] = (uint16_t)T;
    return 0;
}

int rsr_plan_init(rsr_plan *p, int n, int k, int failed, int flavour)
{
    if (flavour == RSR_PLAIN) return rsr_plan_init_grs(p, n, k, failed, NULL);
    if (flavour != RSR_CAUCHY || n > RSR_MAX_N || k <= 0 || k >= n) return -1;
    /* [I ; 1/(i^j)] is the code { u_i f(i) } with u_i = 1 / prod_{t<k, t!=i} (i ^ t) */
    uint8_t u[RSR_MAX_N];
    for (int i = 0; i < n; i++) {
        uint8_t prod = 1;
        for (int t = 0; t < k; t++)
            if (t != i) prod = gf_mul(prod, (uint8_t)(i ^ t));
        u[i] = gf_inv(prod);
    }
    return rsr_plan_init_grs(p, n, k, failed, u);
}

size_t rsr_payload_len(const rsr_plan *p, int helper, size_t len)
{
    if (helper < 0 || helper >= p->n) return 0;
    return (size_t)p->bits[helper] * (64 * (len / SB) + ((len % SB) + 7) / 8);
}

size_t rsr_total_payload_len(const rsr_plan *p, size_t len)
{
    size_t s = 0;
    for (int i = 0; i < p->n; i++) s += rsr_payload_len(p, i, len);
    return s;
}

/* ------------------------------------------------------------------ portable kernels (one block of m <= 512 bytes) */

static void helper_block(const uint8_t *lut, int r, const uint8_t *src, size_t m, uint8_t *dst)
{
    size_t pb = (m + 7) >> 3;
    for (size_t u = 0; u < pb; u++) {
        size_t cnt = m - 8 * u < 8 ? m - 8 * u : 8;
        uint64_t v = 0;
        for (size_t t = 0; t < cnt; t++) v |= (uint64_t)lut[src[8 * u + t]] << (8 * t);
        for (int j = 0; j < r; j++) dst[(size_t)j * pb + u] = (uint8_t)((((v >> j) & LSB8) * GATHER) >> 56);
    }
}

static void rebuild_block(const rsr_plan *p, const uint8_t *const *blk, size_t m, uint8_t *out)
{
    size_t pb = (m + 7) >> 3;
    uint8_t z[8][64];
    memset(z, 0, sizeof z);
    for (int i = 0; i < p->n; i++)
        for (int j = 0; j < p->bits[i]; j++) {
            const uint8_t *pl = blk[i] + (size_t)j * pb;
            unsigned mask = p->dmask[i][j];
            for (int b = 0; b < 8; b++)
                if (mask >> b & 1)
                    for (size_t u = 0; u < pb; u++) z[b][u] ^= pl[u];
        }
    for (size_t u = 0; u < pb; u++) {
        size_t cnt = m - 8 * u < 8 ? m - 8 * u : 8;
        uint64_t w = 0;
        for (int b = 0; b < 8; b++) w |= (uint64_t)z[b][u] << (8 * b);
        for (size_t v = 0; v < cnt; v++) out[8 * u + v] = p->fin_lut[(uint8_t)((((w >> v) & LSB8) * GATHER) >> 56)];
    }
}

/* ------------------------------------------------------------------ x86 vector kernels (whole 512-byte blocks) */

#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
#include <immintrin.h>
#define HAVE_X86 1

/* addresses of the payload planes of block 0, in the order the rebuild network reads them, and their step per block */
static inline int load_terms(const rsr_plan *p, const uint8_t *const *pay, const uint8_t **tp, size_t *ts)
{
    const int T = p->term_start[8];
    for (int t = 0; t < T; t++) {
        const int i = p->term[t].helper;
        tp[t] = pay[i] + 64 * (size_t)p->term[t].plane;
        ts[t] = 64 * (size_t)p->bits[i];
    }
    return T;
}

#define T512 __attribute__((target("avx512f,avx512bw,avx512vl,avx512vbmi,gfni")))

/* r[q] <- { r[0].qword[q], r[1].qword[q], ..., r[7].qword[q] } */
T512 static inline void transpose8(__m512i *r)
{
    __m512i t0 = _mm512_unpacklo_epi64(r[0], r[1]), t1 = _mm512_unpackhi_epi64(r[0], r[1]);
    __m512i t2 = _mm512_unpacklo_epi64(r[2], r[3]), t3 = _mm512_unpackhi_epi64(r[2], r[3]);
    __m512i t4 = _mm512_unpacklo_epi64(r[4], r[5]), t5 = _mm512_unpackhi_epi64(r[4], r[5]);
    __m512i t6 = _mm512_unpacklo_epi64(r[6], r[7]), t7 = _mm512_unpackhi_epi64(r[6], r[7]);
    __m512i u0 = _mm512_shuffle_i64x2(t0, t2, 0x88), u1 = _mm512_shuffle_i64x2(t1, t3, 0x88);
    __m512i u2 = _mm512_shuffle_i64x2(t0, t2, 0xdd), u3 = _mm512_shuffle_i64x2(t1, t3, 0xdd);
    __m512i u4 = _mm512_shuffle_i64x2(t4, t6, 0x88), u5 = _mm512_shuffle_i64x2(t5, t7, 0x88);
    __m512i u6 = _mm512_shuffle_i64x2(t4, t6, 0xdd), u7 = _mm512_shuffle_i64x2(t5, t7, 0xdd);
    r[0] = _mm512_shuffle_i64x2(u0, u4, 0x88);
    r[1] = _mm512_shuffle_i64x2(u1, u5, 0x88);
    r[2] = _mm512_shuffle_i64x2(u2, u6, 0x88);
    r[3] = _mm512_shuffle_i64x2(u3, u7, 0x88);
    r[4] = _mm512_shuffle_i64x2(u0, u4, 0xdd);
    r[5] = _mm512_shuffle_i64x2(u1, u5, 0xdd);
    r[6] = _mm512_shuffle_i64x2(u2, u6, 0xdd);
    r[7] = _mm512_shuffle_i64x2(u3, u7, 0xdd);
}

T512 static size_t helper_avx512(const rsr_plan *p, int i, const uint8_t *src, size_t len, uint8_t *dst)
{
    const int r = p->bits[i];
    const size_t nsb = len / SB;
    uint8_t idx[64];
    for (int a = 0; a < 8; a++)
        for (int b = 0; b < 8; b++) idx[8 * a + b] = (uint8_t)(8 * b + a);
    const __m512i TR = _mm512_loadu_si512(idx);
    const __m512i X = _mm512_set1_epi64((long long)p->h_sel[i]);
    const __m512i BITREV = _mm512_set1_epi64((long long)0x8040201008040201ULL);
    for (size_t s = 0; s < nsb; s++) {
        __m512i u[8];
        /* with the chunk bytes as the matrix operand, one instruction yields, per 8 chunk bytes, byte j = their 8 j-th parity bits */
        for (int m = 0; m < 8; m++)
            u[m] = _mm512_permutexvar_epi8(TR, _mm512_gf2p8affine_epi64_epi8(X, _mm512_loadu_si512(src + SB * s + 64 * m), 0));
        transpose8(u);
        for (int j = 0; j < r; j++)
            _mm512_storeu_si512(dst + (s * (size_t)r + (size_t)j) * 64, _mm512_gf2p8affine_epi64_epi8(u[j], BITREV, 0));
    }
    return nsb * SB;
}

/* XOR of cnt payload planes.  The jump into the unrolled chain has one fixed target per call site, so the branch
 * predictor never misses, whatever the lengths of the eight lists. */
#define XT(i) a = _mm512_xor_si512(a, _mm512_loadu_si512(q[i]));
T512 static inline __attribute__((always_inline)) __m512i xor_planes(const uint8_t *const *q, int cnt)
{
    __m512i a = _mm512_setzero_si512();
    while (cnt > 16) {
        XT(0) XT(1) XT(2) XT(3) XT(4) XT(5) XT(6) XT(7) XT(8) XT(9) XT(10) XT(11) XT(12) XT(13) XT(14) XT(15)
        q += 16;
        cnt -= 16;
    }
    switch (cnt) {
    case 16: XT(15) /* fall through */
    case 15: XT(14) /* fall through */
    case 14: XT(13) /* fall through */
    case 13: XT(12) /* fall through */
    case 12: XT(11) /* fall through */
    case 11: XT(10) /* fall through */
    case 10: XT(9) /* fall through */
    case 9: XT(8) /* fall through */
    case 8: XT(7) /* fall through */
    case 7: XT(6) /* fall through */
    case 6: XT(5) /* fall through */
    case 5: XT(4) /* fall through */
    case 4: XT(3) /* fall through */
    case 3: XT(2) /* fall through */
    case 2: XT(1) /* fall through */
    case 1: XT(0) /* fall through */
    default: break;
    }
    return a;
}
#undef XT

T512 static size_t rebuild_avx512(const rsr_plan *p, const uint8_t *const *pay, size_t len, uint8_t *out)
{
    const size_t nsb = len / SB;
    if (!nsb) return 0;
    const uint8_t *tp[MAX_TERMS];
    size_t ts[MAX_TERMS];
    const int T = load_terms(p, pay, tp, ts);
    const uint16_t *st = p->term_start;
    uint8_t idx[64];
    for (int l = 0; l < 8; l++)
        for (int c = 0; c < 8; c++) idx[8 * l + c] = (uint8_t)(8 * (7 - c) + l);
    const __m512i IDX = _mm512_loadu_si512(idx);
    const __m512i SEL = _mm512_set1_epi64((long long)0x8040201008040201ULL);
    const __m512i FIN = _mm512_set1_epi64((long long)p->fin_mat);
    for (size_t s = 0; s < nsb; s++) {
        __m512i z[8];
        z[0] = xor_planes(tp + st[0], st[1] - st[0]);
        z[1] = xor_planes(tp + st[1], st[2] - st[1]);
        z[2] = xor_planes(tp + st[2], st[3] - st[2]);
        z[3] = xor_planes(tp + st[3], st[4] - st[3]);
        z[4] = xor_planes(tp + st[4], st[5] - st[4]);
        z[5] = xor_planes(tp + st[5], st[6] - st[5]);
        z[6] = xor_planes(tp + st[6], st[7] - st[6]);
        z[7] = xor_planes(tp + st[7], st[8] - st[7]);
        for (int t = 0; t < T; t++) tp[t] += ts[t];
        transpose8(z); /* z[m] = the 8 intermediate plane words of output bytes 64m..64m+63 */
        for (int m = 0; m < 8; m++) {
            __m512i w = _mm512_permutexvar_epi8(IDX, z[m]);
            w = _mm512_gf2p8affine_epi64_epi8(SEL, w, 0); /* planes -> one intermediate byte per output byte */
            _mm512_storeu_si512(out + SB * s + 64 * m, _mm512_gf2p8affine_epi64_epi8(w, FIN, 0));
        }
    }
    return nsb * SB;
}

#define T256 __attribute__((target("avx2")))

T256 static size_t helper_avx2(const rsr_plan *p, int i, const uint8_t *src, size_t len, uint8_t *dst)
{
    const int r = p->bits[i];
    const size_t nsb = len / SB;
    uint8_t hi16[16];
    for (int x = 0; x < 16; x++) hi16[x] = p->h_lut[i][x << 4];
    const __m256i LO = _mm256_broadcastsi128_si256(_mm_loadu_si128((const __m128i *)p->h_lut[i]));
    const __m256i HI = _mm256_broadcastsi128_si256(_mm_loadu_si128((const __m128i *)hi16));
    const __m256i M4 = _mm256_set1_epi8(0x0f);
    for (size_t s = 0; s < nsb; s++)
        for (int q = 0; q < 16; q++) {
            __m256i x = _mm256_loadu_si256((const __m256i *)(src + SB * s + 32 * q));
            __m256i y = _mm256_xor_si256(_mm256_shuffle_epi8(LO, _mm256_and_si256(x, M4)),
                                         _mm256_shuffle_epi8(HI, _mm256_and_si256(_mm256_srli_epi16(x, 4), M4)));
            for (int j = 0; j < r; j++) {
                uint32_t mm = (uint32_t)_mm256_movemask_epi8(_mm256_sll_epi16(y, _mm_cvtsi32_si128(7 - j)));
                memcpy(dst + (s * (size_t)r + (size_t)j) * 64 + 4 * q, &mm, 4);
            }
        }
    return nsb * SB;
}

#define XT(i) a0 = _mm256_xor_si256(a0, _mm256_loadu_si256((const __m256i *)q[i])); a1 = _mm256_xor_si256(a1, _mm256_loadu_si256((const __m256i *)(q[i] + 32)));
T256 static inline __attribute__((always_inline)) void xor_planes256(const uint8_t *const *q, int cnt, uint8_t *z)
{
    __m256i a0 = _mm256_setzero_si256(), a1 = a0;
    while (cnt > 16) {
        XT(0) XT(1) XT(2) XT(3) XT(4) XT(5) XT(6) XT(7) XT(8) XT(9) XT(10) XT(11) XT(12) XT(13) XT(14) XT(15)
        q += 16;
        cnt -= 16;
    }
    switch (cnt) {
    case 16: XT(15) /* fall through */
    case 15: XT(14) /* fall through */
    case 14: XT(13) /* fall through */
    case 13: XT(12) /* fall through */
    case 12: XT(11) /* fall through */
    case 11: XT(10) /* fall through */
    case 10: XT(9) /* fall through */
    case 9: XT(8) /* fall through */
    case 8: XT(7) /* fall through */
    case 7: XT(6) /* fall through */
    case 6: XT(5) /* fall through */
    case 5: XT(4) /* fall through */
    case 4: XT(3) /* fall through */
    case 3: XT(2) /* fall through */
    case 2: XT(1) /* fall through */
    case 1: XT(0) /* fall through */
    default: break;
    }
    _mm256_storeu_si256((__m256i *)z, a0);
    _mm256_storeu_si256((__m256i *)(z + 32), a1);
}
#undef XT

T256 static size_t rebuild_avx2(const rsr_plan *p, const uint8_t *const *pay, size_t len, uint8_t *out)
{
    const size_t nsb = len / SB;
    if (!nsb) return 0;
    const uint8_t *tp[MAX_TERMS];
    size_t ts[MAX_TERMS];
    const int T = load_terms(p, pay, tp, ts);
    const uint16_t *st = p->term_start;
    uint8_t e[32], bs[32];
    for (int t = 0; t < 32; t++) {
        e[t] = (uint8_t)(t >> 3);
        bs[t] = (uint8_t)(1 << (t & 7));
    }
    const __m256i EXP = _mm256_loadu_si256((const __m256i *)e), BITSEL = _mm256_loadu_si256((const __m256i *)bs);
    __m256i F[8];
    for (int b = 0; b < 8; b++) F[b] = _mm256_set1_epi8((char)p->fin[b]);
    for (size_t s = 0; s < nsb; s++) {
        uint8_t z[8][64];
        xor_planes256(tp + st[0], st[1] - st[0], z[0]);
        xor_planes256(tp + st[1], st[2] - st[1], z[1]);
        xor_planes256(tp + st[2], st[3] - st[2], z[2]);
        xor_planes256(tp + st[3], st[4] - st[3], z[3]);
        xor_planes256(tp + st[4], st[5] - st[4], z[4]);
        xor_planes256(tp + st[5], st[6] - st[5], z[5]);
        xor_planes256(tp + st[6], st[7] - st[6], z[6]);
        xor_planes256(tp + st[7], st[8] - st[7], z[7]);
        for (int t = 0; t < T; t++) tp[t] += ts[t];
        for (int q = 0; q < 16; q++) {
            __m256i acc = _mm256_setzero_si256();
            for (int b = 0; b < 8; b++) {
                uint32_t w;
                memcpy(&w, z[b] + 4 * q, 4);
                __m256i x = _mm256_shuffle_epi8(_mm256_set1_epi32((int)w), EXP); /* byte t = plane byte t/8 */
                x = _mm256_cmpeq_epi8(_mm256_and_si256(x, BITSEL), BITSEL);        /* all-ones where bit t is set */
                acc = _mm256_xor_si256(acc, _mm256_and_si256(x, F[b]));
            }
            _mm256_storeu_si256((__m256i *)(out + SB * s + 32 * q), acc);
        }
    }
    return nsb * SB;
}

static int detect_backend(void)
{
    __builtin_cpu_init();
    if (__builtin_cpu_supports("avx512f") && __builtin_cpu_supports("avx512bw") && __builtin_cpu_supports("avx512vl") &&
        __builtin_cpu_supports("avx512vbmi") && __builtin_cpu_supports("gfni"))
        return RSR_BACKEND_AVX512;
    if (__builtin_cpu_supports("avx2")) return RSR_BACKEND_AVX2;
    return RSR_BACKEND_PORTABLE;
}
#else
static int detect_backend(void) { return RSR_BACKEND_PORTABLE; }
#endif

/* ------------------------------------------------------------------ public entry points */

static int g_backend = 0;

static int backend(void)
{
    if (!g_backend) g_backend = detect_backend();
    return g_backend;
}

int rsr_set_backend(int b)
{
    int best = detect_backend();
    g_backend = (b <= 0 || b > best) ? best : b;
    return g_backend;
}

const char *rsr_backend_name(void)
{
    switch (backend()) {
    case RSR_BACKEND_AVX512: return "avx512+gfni";
    case RSR_BACKEND_AVX2: return "avx2";
    default: return "portable";
    }
}

size_t rsr_helper(const rsr_plan *p, int i, const uint8_t *chunk, size_t len, uint8_t *payload)
{
    const int r = (i >= 0 && i < p->n) ? p->bits[i] : 0;
    if (!r) return 0;
    size_t done = 0;
#ifdef HAVE_X86
    if (backend() == RSR_BACKEND_AVX512) done = helper_avx512(p, i, chunk, len, payload);
    else if (backend() == RSR_BACKEND_AVX2) done = helper_avx2(p, i, chunk, len, payload);
#endif
    for (size_t off = done; off < len; off += SB) {
        size_t m = len - off < SB ? len - off : SB;
        helper_block(p->h_lut[i], r, chunk + off, m, payload + (off / SB) * 64 * (size_t)r);
    }
    return rsr_payload_len(p, i, len);
}

void rsr_rebuild(const rsr_plan *p, const uint8_t *const *payloads, size_t len, uint8_t *out)
{
    size_t done = 0;
#ifdef HAVE_X86
    if (backend() == RSR_BACKEND_AVX512) done = rebuild_avx512(p, payloads, len, out);
    else if (backend() == RSR_BACKEND_AVX2) done = rebuild_avx2(p, payloads, len, out);
#endif
    for (size_t off = done; off < len; off += SB) {
        size_t m = len - off < SB ? len - off : SB;
        const uint8_t *blk[RSR_MAX_N];
        for (int i = 0; i < p->n; i++) blk[i] = p->bits[i] ? payloads[i] + (off / SB) * 64 * (size_t)p->bits[i] : NULL;
        rebuild_block(p, blk, m, out + off);
    }
}
