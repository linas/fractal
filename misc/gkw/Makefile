
TOP = ../../tools
INC = $(TOP)/inc
LIB = $(TOP)/lib

EXES= hurwitz polylog-test

all: $(EXES)

hurwitz:	hurwitz.o 
polylog-test:	polylog-test.o polylog.o

.c.o:
	cc -c -g -O2 -I $(INC) $<

.o:
	cc -g -o $@ $^ -L$(LIB) -lgsl -lgslcblas -lfunc -lm 


clean:
	rm -f *.o junk x

realclean: clean
	rm -f $(EXES)  *.png

