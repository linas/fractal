
TOP = ../../tools
INC = $(TOP)/inc
LIB = $(TOP)/lib

EXES= hurwitz polylog

all: $(EXES)

hurwitz:	hurwitz.o 
polylog:	polylog.o 

.c.o:
	cc -c -g -O2 -I $(INC) $<

.o:
	cc -g -o $@ $^ -L$(LIB) -lfunc -lm 


clean:
	rm -f *.o junk x

realclean: clean
	rm -f $(EXES)  *.png

