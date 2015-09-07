
INCLUDES =  -I../../tools/inc -I../../misc/bignum \
   -I ../../generate

FUNC=../../tools/lib/libfunc.a
HIST=../../generate/brat.o ../../generate/util.o


all: lytic-1d lytic-parts lytic-2d

lytic-1d: lytic-1d.o lytic.o ../sum/sum.o
lytic-parts: lytic-parts.o lytic.o ../sum/sum.o
lytic-2d: lytic-2d.o lytic.o ../sum/sum.o $(HIST)

.o:
	cc -o $@ $^ -lm -lstdc++


.C.o:
	cc -c $^  -Wall -O2 $(INCLUDES)
