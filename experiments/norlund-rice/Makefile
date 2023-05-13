
# CC=gcc-4.0
CC= cc
CFLAGS= -g -O2 -Wall

TOP = ../../tools
INC = $(TOP)/inc
LIB = $(TOP)/lib

EXES = baez dsubn-analytic zeta-imag

all: $(EXES)

# executables
baez.o: 
zeta-imag.o:

baez: baez.o 
dsubn-analytic:	dsubn-analytic.o
zeta-imag: zeta-imag.o

.o:
	$(CC) -o $@ $^ -L$(LIB) -lfunc -lgsl -lgslcblas -lm

.c.o:
	$(CC) -c $(CFLAGS) -I $(INC) $<

clean:
	rm -f *.o $(EXES) a.out x y j.c x.c x.dat *.png x.gplot

realclean: clean

dist:
	tar -zcvf x *.c *.h *.gplot README Makefile .cvsignore */*.c */*.h */*.gplot */Makefile */.cvsignore
