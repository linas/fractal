#
# Makefile
#
TOP = ../../tools
INC = $(TOP)/inc
LIB = $(TOP)/lib

EXES = plouffe

all: $(EXES)

plouffe:	plouffe.o $(LIB)/bernoulli.o $(LIB)/binomial.o $(LIB)/harmonic.o

plouffe.o:	plouffe.c

.c.o:
	cc -c -O2 -I $(INC) $<

.o:
	cc -o $@ $^ -lm


clean:
	rm -f *.o a.out junk j.c

realclean: clean
	rm -f $(EXES)
