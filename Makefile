CFLAGS = -Wall -shared -fPIC

.PHONY: all clean local kalign clean_kalign clean_local

all: local kalign

local:
	gcc $(CFLAGS) align.c -o libalign.so
	gcc $(CFLAGS) swalign.c -o libswalign.so -lm
	gcc $(CFLAGS) seqtools.c -o libseqtools.so
	gcc $(CFLAGS) consensus.c -o libconsensus.so

kalign:
	if [ -f kalign/Makefile ]; then make -C kalign; fi

clean: clean_local clean_kalign

clean_kalign:
	if [ -f kalign/Makefile ]; then make -C kalign clean; fi

clean_local:
	rm -f libalign.so libswalign.so libseqtools.so libconsensus.so
