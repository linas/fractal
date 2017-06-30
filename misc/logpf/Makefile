
ANANT=../anant-git/src

INCLUDES =  -I../../tools/inc -I$(ANANT) \
   -I ../../generate

FUNC=../../tools/lib/libfunc.a
HIST=../../generate/brat.o ../../generate/util.o


all: xperiment genfunc-2d

genfunc-2d: genfunc-2d.o $(HIST) $(FUNC)

xperiment: genfunc-2d.o
	ln -f genfunc-2d disclog
	ln -f genfunc-2d disclog_big
	ln -f genfunc-2d loggpf_mag
	ln -f genfunc-2d loggpf_phase

.o:
	cc -o $@ $^ -L$(ANANT) -lanant -lgmp -ldb -lpthread -lm -lstdc++

.c.o:
	cc -c -g $^  -std=gnu11 -Wall -O2 $(INCLUDES)

.C.o:
	cc -c $^  -std=g++11 -Wall -O2 $(INCLUDES)
