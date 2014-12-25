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

RM= rm -f

OBJS=bernoulli.o binomial.o cache.o Farey.o FareyTree.o gcf.o harmonic.o \
	  modular.o moebius.o Prime.o question.o question-new.o stirling.o \
     totient.o

INCS=bernoulli.h binomial.h cache.h cplex.h Farey.h FareyTree.h \
     flt.h gcf.h harmonic.h \
     modular.h moebius.h Prime.h question.h stirling.h totient.h

all: inc $(OBJS)

inc:	$(INCS)
	touch inc
	echo -n "$? " | xargs -d " " -I '{}'  echo cp '{}' $(INC)/'{}'
	echo -n "$? " | xargs -d " " -I '{}'  cp '{}' $(INC)/'{}'

bernoulli.o:	bernoulli.c bernoulli.h
binomial.o:	binomial.c binomial.h
cache.o:	cache.c cache.h
Farey.o:	Farey.C Farey.h
FareyTree.o:	FareyTree.C FareyTree.h
gcf.o:	gcf.c gcf.h
harmonic.o:	harmonic.c harmonic.h
modular.o:	modular.c modular.h
moebius.o:	moebius.c moebius.h cache.h
Prime.o:	Prime.c Prime.h
question.o:	question.C question.h Farey.h
question-new.o:	question-new.c flt.h question.h
stirling.o: stirling.c stirling.h
totient.o:	totient.c totient.h

.c.o:
	$(CC) -c $<
	ar cru $(LIB)/libfunc.a $@

.C.o:
	$(CC) -c $<
	ar cru $(LIB)/libfunc.a $@

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
	$(RM) Prime inc

dist:
	tar -zcvf x *.C *.h *.c Makefile .gitgnore README
