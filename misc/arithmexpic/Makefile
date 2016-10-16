
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
	ln -f genfunc-2d mobius_exp_mag
	ln -f genfunc-2d divisor_exp_mag
	ln -f genfunc-2d divisor_big
	ln -f genfunc-2d liouv_omega_exp_mag
	ln -f genfunc-2d mertens_m
	ln -f genfunc-2d thue_morse
	ln -f genfunc-2d thue_morse_big
	ln -f genfunc-2d xperiment

genfunc.o: genfunc.h
genfunc-2d.o: genfunc.h

.o:
	cc -o $@ $^ -L$(ANANT) -lanant -lgmp -ldb -lpthread -lm -lstdc++

.c.o:
	cc -c -g $^  -std=gnu11 -Wall -O2 $(INCLUDES)

.C.o:
	cc -c $^  -std=gnu++11 -Wall -O2 $(INCLUDES)
