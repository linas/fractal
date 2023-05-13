#
# Makefile
#
TOP = ../../tools
INC = $(TOP)/inc
LIB = $(TOP)/lib

EXES = treemap

all: $(EXES)

treemap:	treemap.o 

treemap.o:	treemap.c

.c.o:
	cc -c -g -O2 -I $(INC) $<

.o:
	cc -o $@ $^ -lm


clean:
	rm -f *.o a.out junk j.c

realclean: clean
	rm -f $(EXES)
