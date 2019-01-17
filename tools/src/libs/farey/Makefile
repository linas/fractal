#
#
# HISTORY:
# created Linas Vepstas January 1994

TOP=../../..

BIN = $(TOP)/bin
LIB = $(TOP)/lib
INC = $(TOP)/inc

# CC = /bin/cc -O2 -O -I $(INC) -DAIXV3
# CC = /bin/fcc -f2 -O -I $(INC)
# CC = cc -Wall -g -O2 -I $(INC) -DANSI_C -DLINUX
CC = cc -Wall -g -O3 -std=gnu11 -I $(INC)
CCP = c++ -Wall -g -O3 -std=c++14 -I $(INC)

RM= rm -f

OBJS=bernoulli.o binomial.o bitops.o \
     cache.o dirichlet.o Farey.o FareyTree.o gcf.o gpf.o \
     harmonic.o isqrt.o modular.o moebius.o necklace.o \
     prime.o question.o question-new.o stern-brocot.o \
     stirling.o totient.o

INCS=bernoulli.h binomial.h bitops.h cache.h cplex.h Farey.h FareyTree.h \
     flt.h gcf.h gpf.h harmonic.h isqrt.h \
     modular.h moebius.h necklace.h prime.h question.h stirling.h totient.h

all: inc $(OBJS)

inc:	$(INCS)
	touch inc
	echo -n "$? " | xargs -d " " -I '{}'  echo cp '{}' $(INC)/'{}'
	echo -n "$? " | xargs -d " " -I '{}'  cp '{}' $(INC)/'{}'

bernoulli.o:	bernoulli.c bernoulli.h
binomial.o:	binomial.c binomial.h
bitops.o: bitops.c bitops.h
cache.o:	cache.c cache.h
dirichlet.o: dirichlet.c
Farey.o:	Farey.C Farey.h
FareyTree.o:	FareyTree.C FareyTree.h
gcf.o:	gcf.c gcf.h
gpf.o:	gpf.c gpf.h
harmonic.o:	harmonic.c harmonic.h
isqrt.o: isqrt.c isqrt.h
modular.o:	modular.c modular.h
moebius.o:	moebius.c moebius.h cache.h
necklace.o: necklace.c necklace.h
prime.o:	prime.c prime.h
question.o:	question.C question.h Farey.h
question-new.o:	question-new.c flt.h question.h
stirling.o: stirling.c stirling.h
totient.o:	totient.c totient.h

.c.o:
	$(CC) -c $<
	ar cruU $(LIB)/libfunc.a $@

.C.o:
	$(CCP) -c $<
	ar cruU $(LIB)/libfunc.a $@

clean:
	cd test; make clean
	cd maps; make clean
	cd prods; make clean
	cd loops; make clean
	$(RM) core a.out *.o tmp junk glop inc

realclean: clean
	cd test; make realclean
	cd maps; make realclean
	cd prods; make realclean
	cd loops; make realclean
	$(RM) inc

dist:
	tar -zcvf x *.C *.h *.c Makefile .gitgnore README
