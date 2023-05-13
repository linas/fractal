
ANANT=../anant-git/src

INCLUDES =  -I../../tools/inc -I$(ANANT) \
   -I ../../generate

FUNC=../../tools/lib/libfunc.a
HIST=../../generate/brat.o ../../generate/util.o

all: multi plic newton

plic: plic.o
newton: newton.o

multi: multi.o
	cc -o $@ $^ $(FUNC) -lpthread -lm -lstdc++

.o:
	cc -o $@ $^ -L$(ANANT) $(FUNC) -lanant -lgmp -ldb -lpthread -lm -lstdc++

.c.o:
	cc -c -g $^ -std=gnu11 -Wall -O2 $(INCLUDES)

.C.o:
	cc -c $^  -Wall -O2 $(INCLUDES)
