
TOP = ../../tools
INC = $(TOP)/inc
LIB = $(TOP)/lib

EXES= asymptote hurwitz 

all: $(EXES)

asymptote:	asymptote.o 
hurwitz:	hurwitz.o 

.c.o:
	cc -c -g -O2 -I $(INC) $<

.o:
	cc -g -o $@ $^ -L$(LIB) -lgsl -lgslcblas -lfunc -lm 


clean:
	rm -f *.o junk x

realclean: clean
	rm -f $(EXES)  *.png

