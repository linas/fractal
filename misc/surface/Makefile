

all: geodesic

geodesic: geodesic.o


.c.o:
	cc -O2 -c $<

.o:
	cc -o $* $^ -lm

