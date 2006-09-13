
TOP = ../../tools
INC = $(TOP)/inc
LIB = $(TOP)/lib


EXES= binom trinary prt-farey

all: $(EXES)

binom:	binom.o $(LIB)/gcf.o $(LIB)/binomial.o
density:	density.o
prt-farey:	prt-farey.o $(LIB)/Farey.o $(LIB)/gcf.o
trinary:	trinary.o


.C.o:
	gcc -c $^ -I $(INC)

.o:
	g++ -o $* $^ -L $(LIB) -lm


.c.o:
	gcc -c $^ -I $(INC)

clean:
	rm -f *.o a.out

realclean: clean
	rm -rf $(EXES)
