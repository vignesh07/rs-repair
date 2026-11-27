# rsrepair - see README.md
CC      ?= cc
CFLAGS  ?= -O3 -Wall -Wextra
ISAL    ?= ../../isa-l          # an Intel ISA-L source tree built with `make -f Makefile.unx` (only for interop_isal / bench)

all: librsrepair.a rsr

librsrepair.a: rsrepair.o
	ar rcs $@ $^

rsrepair.o: rsrepair.c rsrepair.h rsr_tables.h
	$(CC) $(CFLAGS) -c -o $@ rsrepair.c

rsr: rsr_cli.c librsrepair.a
	$(CC) $(CFLAGS) -o $@ rsr_cli.c librsrepair.a

# every code, flavour, lost chunk, length and code path, under the address and undefined-behaviour sanitizers
test: test_rsr.c rsrepair.c rsrepair.h rsr_tables.h
	$(CC) -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer -o test_rsr test_rsr.c rsrepair.c
	./test_rsr

# rebuild the scheme table from the json files in the parent directory (deterministic; about half a minute)
tables: opt
	python3 gen_tables.py > rsr_tables.h.new && mv rsr_tables.h.new rsr_tables.h

opt: opt.c
	$(CC) -O2 -o $@ opt.c

# rebuild shards of stripes written by another library (see interop/README.md for how the stripes are made)
interop_file: interop_file.c librsrepair.a
	$(CC) $(CFLAGS) -o $@ interop_file.c librsrepair.a

interop_isal: interop_isal.c librsrepair.a
	$(CC) $(CFLAGS) -I$(ISAL)/include -o $@ interop_isal.c librsrepair.a $(ISAL)/bin/isa-l.a

bench: bench.c librsrepair.a
	$(CC) $(CFLAGS) -Wno-deprecated-declarations -I$(ISAL)/include -o $@ bench.c librsrepair.a $(ISAL)/bin/isa-l.a

clean:
	rm -f *.o librsrepair.a rsr test_rsr opt interop_file interop_isal bench

.PHONY: all test tables clean
