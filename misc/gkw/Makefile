
TOP = ../../tools
INC = $(TOP)/inc
LIB = $(TOP)/lib

EXES= asymptote diag hurwitz similar

all: $(EXES)

asymptote:	asymptote.o 
diag: diag.o
hurwitz:	hurwitz.o 
similar: similar.o

.c.o:
	cc -c -g -O2 -I $(INC) $<

.o:
	cc -g -o $@ $^ -L$(LIB) -lgsl -lgslcblas -lfunc -lm -lstdc++


clean:
	rm -f *.o junk x

realclean: clean
	rm -f $(EXES)  *.png

