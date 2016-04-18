
ANANT=../anant-git/src

INCLUDES =  -I../../tools/inc -I$(ANANT) \
   -I ../../generate

FUNC=../../tools/lib/libfunc.a
HIST=../../generate/brat.o ../../generate/util.o


all: gpf-gen gpf-2d gpf-zero

gpf-2d: gpf-2d.o gpf-gen-bignum.o $(HIST) $(FUNC)
gpf-gen: gpf-gen.o gpf-gen-bignum.o $(FUNC)
gpf-zero: gpf-zero.o gpf-gen-bignum.o $(FUNC)

gpf-gen-bignum.o: gpf-gen-bignum.h

.o:
	cc -o $@ $^ -L$(ANANT) -lanant -lgmp -ldb -lm -lstdc++

.c.o:
	cc -c $^  -std=gnu11 -Wall -O2 $(INCLUDES)

.C.o:
	cc -c $^  -std=gnu++11  -Wall -O2 $(INCLUDES)
