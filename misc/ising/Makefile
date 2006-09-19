
TOP = ../../tools
INC = $(TOP)/inc
LIB = $(TOP)/lib

EXES= ising kac lr-tree partition takagi-gral

all: $(EXES)


ising:	ising.o $(LIB)/question.o $(LIB)/Farey.o $(LIB)/gcf.o
kac:	kac.o $(LIB)/question.o $(LIB)/Farey.o $(LIB)/gcf.o
lr-tree:	lr-tree.o
partition:	partition.o
takagi-gral:	takagi-gral.o

.c.o:
	cc -c -g -O2 -I $(INC) $<

.o:
	cc -g -o $@ $^ -lstdc++ -lm 


clean:
	rm -f *.o junk x

realclean: clean
	rm -f $(EXES)  *.png

