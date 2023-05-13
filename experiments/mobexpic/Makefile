
ANANT=../anant-git/src

INCLUDES =  -I../../tools/inc -I$(ANANT) \
   -I ../../generate

FUNC=../../tools/lib/libfunc.a
HIST=../../generate/brat.o ../../generate/util.o


all: dirichlet genfunc-2d

genfunc-2d: genfunc-2d.o $(HIST) $(FUNC)
unity: $(FUNC)

dirichlet: genfunc-2d
	ln -f genfunc-2d dirichlet
	ln -f genfunc-2d dirichlet_big

.o:
	cc -o $@ $^ -L$(ANANT) -lanant -lgmp -ldb -lpthread -lm -lstdc++

.c.o:
	cc -c -g $^  -std=gnu11 -Wall -O2 $(INCLUDES)

.C.o:
	cc -c $^  -std=gnu++11 -Wall -O2 $(INCLUDES)
