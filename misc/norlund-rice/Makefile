
# CC=gcc-4.0
CC= cc
CFLAGS= -g -O2 -Wall

TOP = ../../tools
INC = $(TOP)/inc
LIB = $(TOP)/lib

EXES = baez zeta-imag

all: $(EXES)

# executables
baez.o: 
zeta-imag.o:

baez: baez.o $(LIB)/binomial.o $(LIB)/harmonic.o
zeta-imag: zeta-imag.o $(LIB)/binomial.o $(LIB)/harmonic.o

.o:
	$(CC) -o $@ $^ -lgsl -lgslcblas -lm

.c.o:
	$(CC) -c $(CFLAGS) -I $(INC) $<

clean:
	rm -f *.o $(EXES) a.out x y j.c x.c x.dat *.png x.gplot

realclean: clean

dist:
	tar -zcvf x *.c *.h *.gplot README Makefile .cvsignore */*.c */*.h */*.gplot */Makefile */.cvsignore
