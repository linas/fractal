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
	cc -c -g -O2 -I $(INC) $<

.o:
	cc -o $@ $^ -lm


clean:
	rm -f *.o a.out junk j.c

realclean: clean
	rm -f $(EXES)
