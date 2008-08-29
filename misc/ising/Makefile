
TOP = ../../tools
INC = $(TOP)/inc
LIB = $(TOP)/lib

EXES= ising kac lr-tree partition takagi-gral

all: $(EXES)


ising:	ising.o
kac:	kac.o 
lr-tree:	lr-tree.o
partition:	partition.o
takagi-gral:	takagi-gral.o

.c.o:
	cc -c -g -O2 -I $(INC) $<

.o:
	cc -g -o $@ $^ -L $(LIB) -lfunc -lstdc++ -lm 


clean:
	rm -f *.o junk x

realclean: clean
	rm -f $(EXES)  *.png

