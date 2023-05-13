
TOP = ../../tools
INC = $(TOP)/inc
LIB = $(TOP)/lib

EXES= counting eigen eigen-extrap gral-matrix gral-explore sym

all: $(EXES)

counting: counting.o
eigen: eigen.o gral-simple.o
eigen-extrap: eigen-extrap.o gral-simple.o
gral-matrix: gral-matrix.o
gral-explore: gral-explore.o gral-simple.o
sym: sym.o

.c.o:
	cc -c -g -O2 $< -I $(INC)

.o:
	cc -g -o $@ $^ -L $(LIB) -lfunc -lstdc++ -lm
# 	cc -g -o $@ $^ -L $(LIB) -lfunc -llapack -lf77blas -latlas -lg2c -lm


clean:
	rm -f *.o junk x

realclean: clean
	rm -f $(EXES)  *.png

