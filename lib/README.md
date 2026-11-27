# rsrepair: rebuild a lost Reed–Solomon chunk with 25–44% less download, at ISA-L speed

When one chunk of a Reed–Solomon stripe is lost, normal repair downloads **k whole chunks**. With this library every
surviving node instead sends **a few bits per byte** of its chunk (2 to 7, usually 4), and the lost chunk is rebuilt
from those bits. The stripe is ordinary Reed–Solomon as existing libraries write it. Nothing is re-encoded.

It is the repair schemes of this repository turned into working code: a C library (`rsrepair.c`, under 500 lines),
a command-line tool, tests against four real erasure-coding libraries, and measurements against Intel ISA-L.

Status: 26 Nov 2026. A library and its measurements, not yet wired into any storage system. Unreviewed.

## Results

| Code | Where these parameters are used | Download to rebuild one chunk | Rebuild speed, 1 MiB chunks (GB/s, one core) |
|---|---|---|---|
| RS(14,10) | HDFS policy RS-10-4; MinIO default on 14 drives; Facebook f4's parameters | 10 chunks → 6.2 (**−38%**) | ISA-L 2.5 → rsrepair **3.8** |
| RS(9,6) | HDFS policy RS-6-3 | 6 chunks → 4.0 (**−33%**) | ISA-L 4.0 → rsrepair **5.4** |
| RS(5,3) | HDFS policy RS-3-2; MinIO default on 5 drives | 3 chunks → 2.2 (**−28%**) | ISA-L 6.4 → rsrepair **8.4** |
| RS(16,12) | MinIO default on 16 drives | 12 chunks → 7.1 (**−41%**) | ISA-L 2.2 → rsrepair **3.6** |
| RS(12,8) | MinIO default on 12 drives | 8 chunks → 5.2 (**−35%**) | ISA-L 2.8 → rsrepair **4.0** |
| RS(8,4) | MinIO default on 8 drives | 4 chunks → 2.2 (**−44%**) | ISA-L 6.3 → rsrepair **10.5** |
| RS(6,3) | MinIO default on 6 drives | 3 chunks → 2.0 (**−33%**) | ISA-L 6.8 → rsrepair **9.2** |
| RS(4,2) | MinIO default on 4 drives | 2 chunks → 1.5 (**−25%**) | ISA-L 10.1 → rsrepair **13.6** |
| RS(20,17) | Backblaze Vaults (17 data + 3 parity) | 17 chunks → 12.6 (**−26%**) | ISA-L 1.6 → rsrepair **2.1** |
| RS(6,4) | common 4+2 layout | 4 chunks → 3.0 (**−25%**) | ISA-L 5.4 → rsrepair **6.7** |
| RS(10,8) | common 8+2 layout | 8 chunks → 6.0 (**−25%**) | ISA-L 3.4 → rsrepair **4.4** |

"rsrepair" is its AVX-512 + GFNI path, averaged over every possible lost chunk. "ISA-L" is ISA-L's own reconstruction
of one chunk from k chunks (`ec_encode_data`, which picks its AVX-512 + GFNI kernel on this machine). Same machine, one
core. All sizes and the other code paths are in [the full table](#all-measurements).

In short, on the eleven codes above:

- **Download**: 25% to 44% less than normal repair. These are real payload bytes, counted by the library.
- **Rebuild speed**: faster than ISA-L at 1 MiB chunks (1.2–1.7×) and 16 MiB chunks (1.1–2.2×); between 0.8× and 1.8×
  at 64 KiB; between 0.7× and 1.3× at 4 KiB, where ISA-L wins for the narrow codes. The library reads fewer bytes than
  normal repair, which is why it pulls ahead once buffers no longer fit in the cache.
- **Helper cost** (new work on every surviving node): 15–80 GB/s per core, far above any disk or network.
- **End to end on files**: `demo.sh` splits a 256 MiB file with klauspost's own example encoder into 10+4 shards,
  deletes one, and rebuilds it byte-identical from 156.8 MiB of payloads instead of 256 MiB.

## What it costs and where it does not apply

- **More disk reads.** All n−1 surviving nodes read their whole chunk, where normal repair reads k chunks. For
  RS(14,10) that is 13 chunks read instead of 10, to send 6.2 chunks' worth instead of 10. It pays off where the
  network, or a cap on repair traffic, is the limit, and not where disk reads are.
- **More CPU in total, though little.** Normal repair does all its arithmetic on the rebuilding node. Here each of the
  n−1 helpers also makes one pass over its chunk. Summed over all nodes (rebuild time plus n−1 helper times from the
  table below), the coding arithmetic costs 1.1–1.4× the CPU of normal repair at 1 MiB chunks, 1.5–2.5× at 16 MiB,
  1.8–3.1× at 64 KiB and 3–5× at 4 KiB. In absolute terms that is between 0.07 and 2.1 ms of CPU per MiB rebuilt.
  The CPU spent on disk and network I/O is not measured.
- **Every surviving chunk must be reachable.** If a second node is down or slow, fall back to normal repair.
- **One lost chunk at a time.** Helpers must be told which chunk is lost.
- **Codes covered**: n−k of 2, 3 or 4 with n ≤ 16, plus RS(20,17); 37 codes (`./rsr info`). Two families: chunk i
  holds f(i) (Backblaze JavaReedSolomon, klauspost/reedsolomon default) or ISA-L's Cauchy matrix
  (`gf_gen_cauchy1_matrix`, also Hadoop's RS coder). **Not** ISA-L's `gf_gen_rs_matrix`, and Jerasure is untested.
- **CPUs**: the fast path needs AVX-512 with GFNI and VBMI (Intel Ice Lake, AMD Zen 4, and newer). There is an AVX2
  path (1–8 GB/s depending on code and buffer size; slower than ISA-L's AVX2 kernel for narrow codes on small
  buffers, down to 0.4×) and portable C (0.5–1.6 GB/s rebuild, 1.4–2.3 GB/s helper). No ARM NEON path yet.
- **Tiny buffers.** At 4 KiB per chunk ISA-L is up to 1.5× faster for the narrow codes.
- Measured on one cloud VM with two virtual CPUs. Expect different absolute numbers elsewhere.

## Tested

- `make test`: every code, both families, every lost chunk, 14 lengths from 0 to 4351 bytes, all three code paths;
  33,348 rebuilds byte-identical and helper payloads identical across code paths, under the address and
  undefined-behaviour sanitizers with exact-size buffers.
- Against the real libraries, 37 codes and 397 lost shards each, all byte-identical ([interop/README.md](interop/README.md)):
  **Intel ISA-L** v2.32.1; **klauspost/reedsolomon** v1.14.2 and v1.12.4 (the version MinIO pins; MinIO calls it with
  the default matrix); **Backblaze JavaReedSolomon**; **Apache Hadoop**'s RS coder classes (`RSUtil`, `GF256`), which
  serve the HDFS policies RS-3-2, RS-6-3 and RS-10-4.
- The table itself: `gen_tables.py` re-derives it from the scheme files one directory up and repairs random codewords
  of both families with it before writing it out.
- An independent decoder: `check_wire.py` rebuilds shards in Python from the payload files the C tool writes, using
  only the table and the wire format as documented, and confirms the payload sizes.

## Use

```c
#include "rsrepair.h"

rsr_plan plan;
rsr_plan_init(&plan, 14, 10, lost, RSR_CAUCHY);            /* or RSR_PLAIN */

/* on each surviving node i */
size_t plen = rsr_helper(&plan, i, chunk_i, len, payload_i);   /* plen == rsr_payload_len(&plan, i, len) */

/* on the node rebuilding the chunk; payloads[lost] is ignored */
rsr_rebuild(&plan, payloads, len, out);
```

```
make            # librsrepair.a and the rsr tool
make test
./rsr info      # codes and download per lost byte
./rsr helper  plain 14 10 3 0 shard.0 payload.0
./rsr rebuild plain 14 10 3 <shard bytes> shard.3 payload.0 payload.1 payload.2 - payload.4 ... payload.13
./demo.sh <klauspost simple-encoder> <klauspost simple-decoder>
```

Chunks are numbered 0..n−1 in the library's own order, data first. To stream, cut a chunk at multiples of 512 bytes.

## How it works

1. **The scheme.** For each code and each lost chunk, the scheme files of this repository give every helper a small
   set of parity functions: helper i computes r_i fixed XOR-combinations of the bits of each of its bytes. The lost
   byte is a fixed XOR-combination of all the bits received. (This is Guruswami–Wootters trace repair with the
   best-known schemes for these codes, e.g. 50 bits per lost byte for RS(14,10) against 62 in the 2022 tables for the
   same points, and 80 for normal repair.)
2. **Wire format: bit planes.** A chunk is cut into 512-byte blocks. For each block a helper sends r planes of 64
   bytes (fewer for a final partial block); bit t of plane j is its j-th parity bit of byte t. The format is the same
   on every CPU.
3. **Helper.** One GFNI instruction with the *chunk bytes as the matrix operand* yields, for 8 chunk bytes at a time,
   their parity bits already grouped by plane; a byte shuffle and a transpose lay the planes out.
4. **Rebuild.** Because the payloads are planes, the rebuild needs no unpacking: each of 8 intermediate planes is the
   XOR of a fixed list of payload planes, 64 bytes per instruction. One shuffle and two GFNI instructions per 64 bytes
   then turn the 8 planes back into bytes.
5. **Two choices that matter for speed.** `opt.c` picks the basis each helper sends in, and the coordinates of the
   intermediate byte, to make those XOR lists as short as possible (RS(14,10): from 129 plane reads per block to 92).
   And the XOR lists are run as unrolled chains entered by a jump whose target never changes, so the branch predictor
   never misses; this took RS(14,10) from 8.4 to 12 GB/s at 64 KiB chunks.

## All measurements

GB/s of rebuilt data on one core. Machine and versions: [results/machine.txt](results/machine.txt). Raw numbers:
[results/bench.json](results/bench.json). Each figure is the best time over three interleaved rounds (`bench.c`).
"Helper" is the speed at which one surviving node turns its chunk into a payload (AVX-512 path).

| Code | Chunk | ISA-L | klauspost | rsrepair (slowest–fastest lost chunk) | rsrepair / ISA-L | Helper | AVX2 only: ISA-L | AVX2 only: rsrepair | Portable C |
|---|---|---|---|---|---|---|---|---|---|
| RS(14,10) | 4 KiB | 10.9 | 9.6 | **10.7** (10.4–11.3) | 0.99× | 49.7 | 3.9 | 4.8 | 0.9 |
| RS(14,10) | 64 KiB | 8.2 | 12.1 | **12.1** (11.7–12.8) | 1.47× | 70.9 | 3.9 | 5.2 | 0.9 |
| RS(14,10) | 1 MiB | 2.5 | 2.3 | **3.8** (3.7–3.9) | 1.49× | 59.0 | 2.3 | 3.2 | 0.9 |
| RS(14,10) | 16 MiB | 2.0 | 1.6 | **3.7** (3.6–3.7) | 1.83× | 19.1 | 1.9 | 3.3 | 0.8 |
| RS(9,6) | 4 KiB | 22.3 | 14.7 | **16.8** (16.5–17.2) | 0.75× | 49.5 | 6.9 | 5.9 | 1.2 |
| RS(9,6) | 64 KiB | 15.2 | 19.6 | **18.1** (17.8–18.5) | 1.19× | 70.5 | 7.8 | 6.3 | 1.2 |
| RS(9,6) | 1 MiB | 4.0 | 3.6 | **5.4** (5.2–5.7) | 1.35× | 55.8 | 3.7 | 4.5 | 1.2 |
| RS(9,6) | 16 MiB | 3.3 | 3.1 | **5.0** (4.9–5.2) | 1.50× | 17.6 | 3.0 | 4.1 | 1.1 |
| RS(5,3) | 4 KiB | 36.2 | 19.4 | **25.7** (25.2–26.2) | 0.71× | 48.7 | 11.6 | 6.9 | 1.4 |
| RS(5,3) | 64 KiB | 28.8 | 31.6 | **25.2** (24.7–25.7) | 0.88× | 67.0 | 13.2 | 7.2 | 1.4 |
| RS(5,3) | 1 MiB | 6.4 | 6.5 | **8.4** (8.0–8.6) | 1.31× | 58.7 | 5.7 | 6.7 | 1.4 |
| RS(5,3) | 16 MiB | 6.0 | 5.9 | **7.4** (7.2–7.6) | 1.23× | 16.0 | 5.6 | 6.1 | 1.3 |
| RS(16,12) | 4 KiB | 7.7 | 2.9 | **9.1** (9.0–9.3) | 1.19× | 48.1 | 3.0 | 4.4 | 0.8 |
| RS(16,12) | 64 KiB | 6.5 | 3.8 | **10.5** (10.3–10.7) | 1.61× | 70.6 | 3.2 | 4.7 | 0.9 |
| RS(16,12) | 1 MiB | 2.2 | 1.0 | **3.6** (3.6–3.7) | 1.64× | 66.3 | 2.2 | 3.4 | 0.8 |
| RS(16,12) | 16 MiB | 1.1 | 0.9 | **2.5** (2.3–3.1) | 2.21× | 16.6 | 1.1 | 2.5 | 0.7 |
| RS(12,8) | 4 KiB | 16.7 | 11.1 | **11.8** (11.5–12.2) | 0.71× | 48.2 | 4.7 | 5.1 | 1.0 |
| RS(12,8) | 64 KiB | 10.9 | 14.5 | **13.6** (13.0–14.4) | 1.24× | 69.6 | 4.8 | 5.5 | 1.0 |
| RS(12,8) | 1 MiB | 2.8 | 2.8 | **4.0** (4.0–4.1) | 1.43× | 60.9 | 2.7 | 3.3 | 1.0 |
| RS(12,8) | 16 MiB | 1.9 | 2.1 | **3.4** (3.3–3.5) | 1.85× | 16.4 | 1.7 | 2.6 | 0.8 |
| RS(8,4) | 4 KiB | 30.5 | 16.8 | **22.4** (21.9–22.7) | 0.74× | 54.8 | 9.3 | 6.6 | 1.4 |
| RS(8,4) | 64 KiB | 23.8 | 25.5 | **23.8** (23.4–24.0) | 1.00× | 81.3 | 11.1 | 7.2 | 1.4 |
| RS(8,4) | 1 MiB | 6.3 | 5.0 | **10.5** (10.1–11.1) | 1.68× | 77.8 | 6.2 | 7.0 | 1.4 |
| RS(8,4) | 16 MiB | 5.6 | 4.6 | **8.7** (8.7–8.8) | 1.56× | 21.9 | 5.5 | 6.4 | 1.3 |
| RS(6,3) | 4 KiB | 35.8 | 18.9 | **25.5** (24.9–26.2) | 0.71× | 52.5 | 11.4 | 6.8 | 1.4 |
| RS(6,3) | 64 KiB | 29.8 | 31.4 | **26.8** (25.7–28.1) | 0.90× | 76.7 | 13.8 | 7.4 | 1.4 |
| RS(6,3) | 1 MiB | 6.8 | 6.5 | **9.2** (9.1–9.6) | 1.35× | 67.9 | 6.1 | 7.1 | 1.4 |
| RS(6,3) | 16 MiB | 6.4 | 5.9 | **8.1** (8.0–8.3) | 1.27× | 18.5 | 5.8 | 6.4 | 1.3 |
| RS(4,2) | 4 KiB | 45.4 | 21.9 | **29.9** (29.8–30.0) | 0.66× | 48.8 | 15.3 | 7.3 | 1.5 |
| RS(4,2) | 64 KiB | 39.2 | 39.9 | **32.3** (32.2–32.4) | 0.82× | 69.2 | 20.1 | 7.7 | 1.6 |
| RS(4,2) | 1 MiB | 10.1 | 8.9 | **13.6** (13.3–13.9) | 1.35× | 57.6 | 9.0 | 7.6 | 1.5 |
| RS(4,2) | 16 MiB | 8.5 | 7.7 | **9.3** (9.3–9.5) | 1.10× | 17.1 | 7.7 | 6.8 | 1.4 |
| RS(20,17) | 4 KiB | 5.4 | 2.3 | **6.8** (6.6–7.0) | 1.26× | 47.4 | 2.3 | 3.3 | 0.6 |
| RS(20,17) | 64 KiB | 4.3 | 2.6 | **7.7** (7.5–7.8) | 1.77× | 65.7 | 2.2 | 3.6 | 0.6 |
| RS(20,17) | 1 MiB | 1.6 | 1.2 | **2.1** (2.1–2.4) | 1.33× | 52.2 | 1.6 | 2.0 | 0.6 |
| RS(20,17) | 16 MiB | 0.7 | 0.6 | **1.2** (1.1–1.3) | 1.61× | 14.6 | 0.7 | 1.1 | 0.5 |
| RS(6,4) | 4 KiB | 30.5 | 17.0 | **21.8** (21.6–22.3) | 0.71× | 48.5 | 9.3 | 6.6 | 1.3 |
| RS(6,4) | 64 KiB | 23.8 | 25.5 | **22.3** (21.9–23.0) | 0.94× | 67.7 | 11.1 | 7.0 | 1.3 |
| RS(6,4) | 1 MiB | 5.4 | 5.0 | **6.7** (6.6–6.9) | 1.24× | 50.8 | 4.9 | 5.7 | 1.3 |
| RS(6,4) | 16 MiB | 4.9 | 4.5 | **6.1** (6.0–6.2) | 1.24× | 15.9 | 4.3 | 5.0 | 1.1 |
| RS(10,8) | 4 KiB | 16.3 | 11.4 | **12.3** (11.9–12.8) | 0.76× | 46.6 | 4.7 | 5.2 | 1.0 |
| RS(10,8) | 64 KiB | 10.9 | 15.0 | **13.5** (12.9–13.9) | 1.25× | 62.1 | 4.7 | 5.5 | 1.0 |
| RS(10,8) | 1 MiB | 3.4 | 2.8 | **4.4** (4.2–4.6) | 1.29× | 51.3 | 3.4 | 4.0 | 1.0 |
| RS(10,8) | 16 MiB | 2.4 | 2.0 | **3.4** (2.9–3.7) | 1.43× | 16.3 | 2.3 | 2.5 | 0.8 |

## Earlier work

- V. Guruswami, M. Wootters, *Repairing Reed–Solomon codes*, 2017: the trace-repair framework.
- T. X. Dinh, L. Y. N. Nguyen, L. J. Mohan, S. Boztas, T. T. Luong, S. H. Dau, *Practical Considerations in Repairing
  Reed–Solomon Codes*, ISIT 2022, arXiv:2205.11015: the first tables of schemes for these codes, and a C prototype on
  top of ISA-L (github.com/dausonhoang/tracerepair). They report the prototype as 2.8× slower than ISA-L for RS(9,6)
  (paper, Section V) and 1.8–2.4× slower for five other codes (Dinh's thesis, Section 5.1.4), and leave further
  optimisation open. That prototype is the only public implementation we found.
- The schemes used here are the best-known ones from the parent directory of this repository.

## Files

| file | what it is |
|---|---|
| `rsrepair.h`, `rsrepair.c` | the library |
| `rsr_tables.h` | scheme table, generated |
| `gen_tables.py`, `opt.c` | generate and self-check the table from `../schemes_*.json` |
| `rsr_cli.c`, `demo.sh` | command-line tool and the end-to-end file demo |
| `test_rsr.c` | self-contained correctness test |
| `check_wire.py` | independent Python decoder of the wire format |
| `interop_isal.c`, `interop_file.c`, `interop/` | tests against ISA-L, klauspost, Backblaze, Hadoop |
| `bench.c`, `results/` | measurements |
