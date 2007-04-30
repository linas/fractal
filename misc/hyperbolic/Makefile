
TOP = ../../tools
INC = $(TOP)/inc
LIB = $(TOP)/lib

EXES= mod-tree tree

all: $(EXES)

mod-tree: mod-tree.o eps.o flt-eps.o
tree: tree.o eps.o flt-eps.o

.c.o:
	cc -c -g -O2 -I $(INC) $<

.o:
	cc -g -o $@ $^ -L$(LIB) -lm 


clean:
	rm -f *.o junk x

realclean: clean
	rm -f $(EXES)  *.png


tree.o:	tree.c eps.h flt.h flt-eps.h
