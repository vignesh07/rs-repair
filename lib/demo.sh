#!/bin/sh
# End-to-end demo on real files.
#   1. klauspost/reedsolomon's own example encoder splits a random file into 10 data + 4 parity shards.
#   2. One shard is deleted.
#   3. Every other shard is squeezed by `rsr helper`; the lost shard is rebuilt by `rsr rebuild` from those payloads.
#   4. The rebuilt shard is compared with the original, and klauspost's own decoder reassembles the file.
# usage: demo.sh <simple-encoder> <simple-decoder> [file size in MiB] [lost shard]
set -e
ENC=$1; DEC=$2; MB=${3:-256}; LOST=${4:-3}; N=14; K=10
HERE=$(cd "$(dirname "$0")" && pwd)
W=$(mktemp -d); trap 'rm -rf "$W"' EXIT
mkdir "$W/shards"
head -c $((MB * 1024 * 1024)) /dev/urandom > "$W/file"
"$ENC" -data $K -par $((N - K)) -out "$W/shards" "$W/file" > /dev/null
SH=$(wc -c < "$W/shards/file.0" | tr -d " ")
mv "$W/shards/file.$LOST" "$W/lost.orig"
echo "file: $MB MiB, RS($N,$K) encoded by klauspost/reedsolomon, shard size $SH bytes; shard $LOST deleted"

ARGS=""; TOTAL=0; i=0
while [ $i -lt $N ]; do
    if [ $i -eq $LOST ]; then ARGS="$ARGS -"; else
        "$HERE/rsr" helper plain $N $K $LOST $i "$W/shards/file.$i" "$W/pay.$i"
        TOTAL=$((TOTAL + $(wc -c < "$W/pay.$i" | tr -d " ")))
        ARGS="$ARGS $W/pay.$i"
    fi
    i=$((i + 1))
done
"$HERE/rsr" rebuild plain $N $K $LOST $SH "$W/shards/file.$LOST" $ARGS
cmp "$W/shards/file.$LOST" "$W/lost.orig" && echo "rebuilt shard is byte-identical to the deleted one"
"$DEC" -data $K -par $((N - K)) -out "$W/file.back" "$W/shards/file" > /dev/null
# klauspost's example decoder does not know the file length and pads the output to a multiple of the shard size
[ "$(sha256sum < "$W/file")" = "$(head -c $((MB * 1024 * 1024)) "$W/file.back" | sha256sum)" ] && echo "klauspost's own decoder accepts the 14 shards and returns the original file (sha256 matches)"
awk -v t=$TOTAL -v n=$((K * SH)) 'BEGIN {
    printf "downloaded to rebuild: %.1f MiB instead of %.1f MiB with normal repair (%.1f%% less)\n", t / 1048576, n / 1048576, 100 * (1 - t / n) }'
