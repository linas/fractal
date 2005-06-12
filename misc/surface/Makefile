

all: geodesic period

geodesic: geodesic.o geo-lib.o
period: period.o geo-lib.o

geodesic.o: geo-lib.h

geo-lib.o: geo-lib.h

.c.o:
	cc -O2 -c $<

.o:
	cc -o $* $^ -lm

