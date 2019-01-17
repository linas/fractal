
ANANT=../anant-git/src

INCLUDES =  -I../../tools/inc -I$(ANANT) \
   -I ../../generate

FUNC=../../tools/lib/libfunc.a
HIST=../../generate/brat.o ../../generate/util.o


all: genfunc-2d unity

unity.o: unity.c

genfunc-2d: genfunc-2d.o $(HIST) $(FUNC)
unity: $(FUNC)

totient_ord_phase: genfunc-2d
	ln -f genfunc-2d mobius_exp_mag
	ln -f genfunc-2d mobius_big
	ln -f genfunc-2d divisor_exp_mag
	ln -f genfunc-2d divisor_big

.o:
	cc -o $@ $^ -L$(ANANT) -lanant -lgmp -ldb -lpthread -lm -lstdc++

.c.o:
	cc -c -g $^  -std=gnu11 -Wall -O2 $(INCLUDES)

.C.o:
	cc -c $^  -std=gnu++11 -Wall -O2 $(INCLUDES)
