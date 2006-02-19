
CC=gcc-4.0
CFLAGS= -g -O2 -Wall

TOP = ../../tools
INC = $(TOP)/inc
LIB = $(TOP)/lib

EXES = baez

all: $(EXES)

# executables
baez.o: 

baez: baez.o

.o:
	$(CC) -o $@ $^ -lgsl -lgslcblas -lm

.c.o:
	$(CC) -c $(CFLAGS) -I $(INC) $<

clean:
	rm -f *.o $(EXES) a.out x y j.c x.c x.dat *.png x.gplot

realclean: clean

dist:
	tar -zcvf x *.c *.h *.gplot README Makefile .cvsignore */*.c */*.h */*.gplot */Makefile */.cvsignore
