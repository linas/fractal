
ANANT=../anant-git/src

INCLUDES =  -I../../tools/inc -I$(ANANT) \
   -I ../../generate

FUNC=../../tools/lib/libfunc.a
HIST=../../generate/brat.o ../../generate/util.o


all: distrib gpf-gen gpf-2d gpf-zero scribe

distrib: distrib.o $(FUNC)
gpf-2d: gpf-2d.o gpf-gen-bignum.o $(HIST) $(FUNC)
gpf-gen: gpf-gen.o gpf-gen-bignum.o $(FUNC)
gpf-zero: gpf-zero.o gpf-gen-bignum.o $(FUNC)

gpf-gen-bignum.o: gpf-gen-bignum.h
gpf-2d.o: ../../generate/brat.h
scribe.o: ../../generate/brat.h

scribe: scribe.o $(HIST)

.o:
	cc -o $@ $^ -L$(ANANT) -lanant -lgmp -ldb -lpthread -lm -lstdc++

.c.o:
	cc -c -g $^  -std=gnu11 -Wall -O2 $(INCLUDES)

.C.o:
	c++ -c $^  -std=gnu++11 -Wall -g -O2 $(INCLUDES)
