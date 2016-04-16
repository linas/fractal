
INCLUDES =  -I../../tools/inc -I../../misc/bignum \
   -I ../../generate

FUNC=../../tools/lib/libfunc.a
HIST=../../generate/brat.o ../../generate/util.o


all: gpf-gen gpf-2d

gpf-gen: gpf-gen.o $(FUNC)
gpf-2d: gpf-2d.o $(HIST)

.o:
	cc -o $@ $^ -lm -lstdc++

.c.o:
	cc -c $^  -std=c99 -Wall -O2 $(INCLUDES)

.C.o:
	cc -c $^  -Wall -O2 $(INCLUDES)
