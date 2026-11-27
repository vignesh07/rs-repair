# Interop tests against the real libraries

Each test lets the other library encode one random stripe per supported code (shards of 8705 bytes), then
`interop_file` drops every shard in turn, rebuilds it from the helper payloads and compares it byte for byte with what
that library wrote. 37 codes and 397 lost shards per library.

| library | version tested | how the stripes are made | flavour |
|---|---|---|---|
| Intel ISA-L | v2.32.1 (commit af6a939) | `interop_isal.c` calls `gf_gen_cauchy1_matrix`, `ec_init_tables`, `ec_encode_data` directly | cauchy |
| klauspost/reedsolomon | v1.14.2, and v1.12.4 (the version MinIO pins) | `kp/main.go gen`, `reedsolomon.New(k, n-k)` with default options | plain |
| Backblaze JavaReedSolomon | commit d3c481d | `java/BBStripes.java`, `ReedSolomon.create(k, n-k).encodeParity` | plain |
| Apache Hadoop RS coder | commit 3b1dab9 | `java/HdfsStripes.java` drives the unmodified `RSUtil` and `GF256` classes the way `RSRawEncoder` does | cauchy |

```
# ISA-L
git clone https://github.com/intel/isa-l && make -C isa-l -f Makefile.unx
make interop_isal ISAL=/path/to/isa-l && ./interop_isal

# klauspost/reedsolomon  (CODES = "4 2 5 3 ... 20 17", every pair printed by `./rsr info`)
cd interop/kp && go mod init kp && go get github.com/klauspost/reedsolomon@v1.14.2 && go build -o kp .
mkdir stripes && ./kp gen stripes 8705 $CODES && ../../interop_file stripes plain klauspost

# Backblaze JavaReedSolomon
git clone https://github.com/Backblaze/JavaReedSolomon
javac -d classes JavaReedSolomon/src/main/java/com/backblaze/erasure/*.java interop/java/BBStripes.java
mkdir stripes_bb && java -cp classes BBStripes stripes_bb 8705 $CODES && ./interop_file stripes_bb plain backblaze

# Hadoop (only the three coder source files are needed; the stub replaces one annotation class)
U=hadoop/hadoop-common-project/hadoop-common/src/main/java/org/apache/hadoop/io/erasurecode/rawcoder/util
javac -d classes interop/java/stub/org/apache/hadoop/classification/InterfaceAudience.java \
      $U/GF256.java $U/GaloisField.java $U/RSUtil.java interop/java/HdfsStripes.java
mkdir stripes_hdfs && java -cp classes HdfsStripes stripes_hdfs 8705 $CODES && ./interop_file stripes_hdfs cauchy hadoop
```

`kp/main.go bench n k ...` also times klauspost's own one-shard `ReconstructSome` (single goroutine); those are the
`klauspost_MBps` figures in `../results/bench.json`.

Not covered: Intel ISA-L's other matrix, `gf_gen_rs_matrix` (it is not the same family of codes on the points 0..n-1;
checked numerically for seven codes), and Jerasure's matrices (not tested).
