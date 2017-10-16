CFLAGS = -Wall -shared -fPIC

all: local make_kalign

local:
	gcc $(CFLAGS) align.c -o libalign.so
	gcc $(CFLAGS) swalign.c -o libswalign.so -lm
	gcc $(CFLAGS) seqtools.c -o libseqtools.so
	gcc $(CFLAGS) consensus.c -o libconsensus.so

make_kalign:
	cd kalign && make

clean: clean_local clean_kalign

clean_kalign:
	cd kalign && make clean

clean_local:
	rm -f libalign.so libswalign.so libseqtools.so libconsensus.so
