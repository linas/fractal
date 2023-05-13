
INCLUDES =  -I../../tools/inc -I../../misc/bignum \
	-I ../../generate

FUNC=../../tools/lib/libfunc.a
HIST=../../generate/brat.o ../../generate/util.o


all: sum-1d sum-2d

sum-1d: sum-1d.o sum.o
sum-2d: sum-2d.o sum.o $(HIST)

.o:
	cc -o $@ $^ -lm -lstdc++


.C.o:
	cc -c $^  -Wall -O2 $(INCLUDES)

