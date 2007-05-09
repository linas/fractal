
TOP = ../../tools
INC = $(TOP)/inc
LIB = $(TOP)/lib

EXES= mod-tree perimeter rotations tree

all: $(EXES)

mod-tree: mod-tree.o eps.o flt-eps.o
perimeter: perimeter.o 
rotations: rotations.o 
tree: tree.o eps.o flt-eps.o

.c.o:
	cc -c -g -O2 -I $(INC) $<

.o:
	cc -g -o $@ $^ -L$(LIB) -lfunc -lm 


clean:
	rm -f *.o junk x

realclean: clean
	rm -f $(EXES)  *.png


rotations.o:	rotations.c
perimeter.o:	perimeter.c 
mod-tree.o:	mod-tree.c eps.h flt-eps.h
tree.o:	tree.c eps.h flt-eps.h
