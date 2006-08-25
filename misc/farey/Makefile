

EXES= binom density trinary

all: $(EXES)

binom:	binom.o gcf.o
density:	density.o
trinary:	trinary.o


TOP = ../../tools
INC = $(TOP)/inc
LIB = $(TOP)/lib

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
