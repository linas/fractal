
TOP = ../../tools
INC = $(TOP)/inc
LIB = $(TOP)/lib

EXES= asymptote diag hurwitz 

all: $(EXES)

asymptote:	asymptote.o 
diag: diag.o
hurwitz:	hurwitz.o 

.c.o:
	cc -c -g -O2 -I $(INC) $<

.o:
	cc -g -o $@ $^ -L$(LIB) -lgsl -lgslcblas -lfunc -lm 


clean:
	rm -f *.o junk x

realclean: clean
	rm -f $(EXES)  *.png

