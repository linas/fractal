
ANANT=../anant-git/src

INCLUDES =  -I../../tools/inc -I$(ANANT) \
   -I ../../generate

FUNC=../../tools/lib/libfunc.a
HIST=../../generate/brat.o ../../generate/util.o


all: genfunc-2d totient_ord_phase oned find-zero slice

genfunc-2d: genfunc-2d.o genfunc.o $(HIST) $(FUNC)
oned: genfunc.o $(FUNC)
find-zero: genfunc.o $(FUNC)
slice: genfunc.o $(FUNC)

totient_ord_phase: genfunc-2d.o genfunc.o
	ln -f genfunc-2d totient_ord_phase
	ln -f genfunc-2d totient_exp_phase
	ln -f genfunc-2d totient_exp_mag
	ln -f genfunc-2d totient_big
	ln -f genfunc-2d carmichael
	ln -f genfunc-2d carmichael_big
	ln -f genfunc-2d mobius_exp_mag
	ln -f genfunc-2d mobius_big
	ln -f genfunc-2d divisor_exp_mag
	ln -f genfunc-2d divisor_big
	ln -f genfunc-2d divisor_twist
	ln -f genfunc-2d divisor_uncircle
	ln -f genfunc-2d sigma_one
	ln -f genfunc-2d sigma_two
	ln -f genfunc-2d little_omega
	ln -f genfunc-2d big_omega_exp_mag
	ln -f genfunc-2d big_omega_big
	ln -f genfunc-2d unitary
	ln -f genfunc-2d liouv_lambda
	ln -f genfunc-2d mertens_m
	ln -f genfunc-2d mangoldt_lambda
	ln -f genfunc-2d mangoldt_lambda_big
	ln -f genfunc-2d exp_mangoldt_lambda
	ln -f genfunc-2d exp_mango_big
	ln -f genfunc-2d thue_morse
	ln -f genfunc-2d thue_morse_big
	ln -f genfunc-2d isqrt_big
	ln -f genfunc-2d partition_big
	ln -f genfunc-2d xperiment

genfunc.o: genfunc.h
genfunc-2d.o: genfunc.h

oned.o: genfunc.h
find-zero.o: genfunc.h
slice.o: genfunc.h

.o:
	cc -o $@ $^ -L$(ANANT) -lanant -lgmp -ldb -lpthread -lm -lstdc++

.c.o:
	cc -c -g $^  -std=gnu11 -Wall -O2 $(INCLUDES)

.C.o:
	cc -c $^  -std=gnu11 -Wall -O2 $(INCLUDES)
