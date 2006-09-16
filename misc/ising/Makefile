
TOP = ../../tools
INC = $(TOP)/inc
LIB = $(TOP)/lib

EXES= ising takagi-gral

all: $(EXES)


ising:	ising.o $(LIB)/question.o $(LIB)/Farey.o $(LIB)/gcf.o
takagi-gral:	takagi-gral.o

.c.o:
	cc -c -O2 -I $(INC) $<

.o:
	cc -o $@ $^ -lstdc++ -lm 


clean:
	rm -f *.o junk x

realclean: clean
	rm -f $(EXES)  *.png

