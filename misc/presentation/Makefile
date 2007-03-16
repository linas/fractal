#
# Makefile
#
TOP = ../../tools
INC = $(TOP)/inc
LIB = $(TOP)/lib

EXES = present

all: $(EXES)

present:	present.o 

present.o:	present.c

.c.o:
	cc -c -g -pg -O2 -I $(INC) $<

.o:
	cc -pg -o $@ $^ -lm 


clean:
	rm -f *.o a.out junk j.c

realclean: clean
	rm -f $(EXES)
