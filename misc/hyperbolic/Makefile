
TOP = ../../tools
INC = $(TOP)/inc
LIB = $(TOP)/lib

EXES= tree

all: $(EXES)

tree: tree.o eps.o

.c.o:
	cc -c -g -O2 -I $(INC) $<

.o:
	cc -g -o $@ $^ -L$(LIB) -lm 


clean:
	rm -f *.o junk x

realclean: clean
	rm -f $(EXES)  *.png


tree.o:	tree.c eps.h flt.h
