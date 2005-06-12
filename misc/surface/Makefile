

all: geodesic

geodesic: geodesic.o

geodesic.o: geo-lib.h

geo-lib.o: geo-lib.h

.c.o:
	cc -O2 -c $<

.o:
	cc -o $* $^ geo-lib.o -lm

