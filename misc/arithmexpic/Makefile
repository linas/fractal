
ANANT=../anant-git/src

INCLUDES =  -I../../tools/inc -I$(ANANT) \
   -I ../../generate

FUNC=../../tools/lib/libfunc.a
HIST=../../generate/brat.o ../../generate/util.o


all: genfunc-2d totient_ord_phase

genfunc-2d: genfunc-2d.o genfunc.o $(HIST) $(FUNC)

totient_ord_phase: genfunc-2d.o genfunc.o
	ln -f genfunc-2d totient_ord_phase
	ln -f genfunc-2d totient_exp_phase
	ln -f genfunc-2d totient_exp_mag

genfunc.o: genfunc.h

.o:
	cc -o $@ $^ -L$(ANANT) -lanant -lgmp -ldb -lpthread -lm -lstdc++

.c.o:
	cc -c -g $^  -std=gnu11 -Wall -O2 $(INCLUDES)

.C.o:
	cc -c $^  -std=gnu++11 -Wall -O2 $(INCLUDES)
