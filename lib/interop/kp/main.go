// Stripes and baseline timings from klauspost/reedsolomon (default matrix, as shipped).
//
//	kp gen   <outdir> <shardlen> n k [n k ...]   write every shard of one random stripe per code
//	kp bench n k [n k ...]                       time the library's own one-shard reconstruction, single goroutine
package main

import (
	"fmt"
	"math/rand"
	"os"
	"path/filepath"
	"strconv"
	"time"

	"github.com/klauspost/reedsolomon"
)

func stripe(n, k, l int, rng *rand.Rand, opts ...reedsolomon.Option) (reedsolomon.Encoder, [][]byte) {
	enc, err := reedsolomon.New(k, n-k, opts...)
	if err != nil {
		panic(err)
	}
	shards := make([][]byte, n)
	for i := range shards {
		shards[i] = make([]byte, l)
		if i < k {
			rng.Read(shards[i])
		}
	}
	if err := enc.Encode(shards); err != nil {
		panic(err)
	}
	if ok, _ := enc.Verify(shards); !ok {
		panic("verify failed")
	}
	return enc, shards
}

func main() {
	rng := rand.New(rand.NewSource(20261126))
	switch os.Args[1] {
	case "gen":
		dir := os.Args[2]
		l, _ := strconv.Atoi(os.Args[3])
		for a := 4; a+1 < len(os.Args); a += 2 {
			n, _ := strconv.Atoi(os.Args[a])
			k, _ := strconv.Atoi(os.Args[a+1])
			_, shards := stripe(n, k, l, rng)
			f, err := os.Create(filepath.Join(dir, fmt.Sprintf("stripe_%d_%d.bin", n, k)))
			if err != nil {
				panic(err)
			}
			for _, s := range shards {
				f.Write(s)
			}
			f.Close()
		}
	case "bench":
		fmt.Println("[")
		first := true
		for a := 2; a+1 < len(os.Args); a += 2 {
			n, _ := strconv.Atoi(os.Args[a])
			k, _ := strconv.Atoi(os.Args[a+1])
			for _, l := range []int{4 << 10, 64 << 10, 1 << 20, 16 << 20} {
				enc, shards := stripe(n, k, l, rng, reedsolomon.WithMaxGoroutines(1))
				want := append([]byte(nil), shards[0]...)
				required := make([]bool, n)
				required[0] = true
				best := time.Duration(1 << 62)
				start := time.Now()
				for it := 0; it < 5 || time.Since(start) < 150*time.Millisecond; it++ {
					shards[0] = shards[0][:0] // lost; the library reuses the capacity
					t0 := time.Now()
					if err := enc.ReconstructSome(shards, required); err != nil {
						panic(err)
					}
					if d := time.Since(t0); d < best {
						best = d
					}
				}
				if string(shards[0]) != string(want) {
					panic("wrong reconstruction")
				}
				if !first {
					fmt.Println(",")
				}
				first = false
				fmt.Printf("  {\"n\": %d, \"k\": %d, \"chunk_bytes\": %d, \"klauspost_MBps\": %.0f}", n, k, l, float64(l)/1048576/best.Seconds())
			}
		}
		fmt.Println("\n]")
	}
}
