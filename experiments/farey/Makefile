
TOP = ../../tools
INC = $(TOP)/inc
LIB = $(TOP)/lib


EXES= binom divisor-rep trinary prt-farey st-xform

all: $(EXES)

binom:	binom.o 
density:	density.o
divisor-rep:	divisor-rep.o 
prt-farey:	prt-farey.o 
st-xform:	st-xform.o
trinary:	trinary.o


.C.o:
	gcc -c $^ -I $(INC)

.o:
	gcc -o $* $^ -L $(LIB) -lfunc -lm -lstdc++


.c.o:
	gcc -c $^ -I $(INC)

clean:
	rm -f *.o a.out

realclean: clean
	rm -rf $(EXES)
