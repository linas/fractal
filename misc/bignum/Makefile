
# build Gnu MP code
#
# CC = cc -pg
CC = cc


EXES= bernoulli bound bs bsubn cache-fill cont_frac db-merge db-prt \
      dsubn exact hurwitz minima mp_zeta_test picturez \
      poly-perf sho stieltjes trunc zero-cross

MPLIB=libzmp.a

all:  $(MPLIB) $(EXES)

MPOBJS= db-cache.o mp-binomial.o mp-cache.o mp-consts.o \
	mp-gamma.o mp-hyper.o mp-misc.o mp-polylog.o mp-trig.o mp-zeta.o

bernoulli:	bernoulli.o $(MPLIB)
bound:	bound.o
cache-fill:	cache-fill.o $(MPLIB)
db-merge:	db-merge.o $(MPLIB)
db-prt:	db-prt.o $(MPLIB)
cont_frac:	cont_frac.o
dsubn:	dsubn.o $(MPLIB)
exact:	exact.o $(MPLIB)
hurwitz:	hurwitz.o $(MPLIB)
mp_zeta_test:	mp_zeta_test.o $(MPLIB)
picturez:	picturez.o $(MPLIB)
poly-perf:	poly-perf.o $(MPLIB)
sho:	sho.o $(MPLIB)
stieltjes:	stieltjes.o $(MPLIB)
trunc:	trunc.o
zero-cross:	zero-cross.o

# libarary objects
db-cache.o: db-cache.h
mp-binomial.o: mp-binomial.h mp-cache.h mp-misc.h mp-trig.h
mp-cache.o: mp-cache.h mp-complex.h
mp-consts.o: mp-consts.h mp-binomial.h mp-complex.h mp-trig.h mp-zeta.h
mp-gamma.o: mp-gamma.h mp-binomial.h mp-complex.h mp-consts.h mp-misc.h mp-trig.h mp-zeta.h
mp-hyper.o: mp-hyper.h mp-complex.h mp-misc.h
mp-misc.o: mp-misc.h mp-complex.h
mp-polylog.o: mp-polylog.h mp-binomial.h mp-cache.h mp-complex.h mp-consts.h mp-gamma.h mp-misc.h mp-trig.h mp-zeta.h
mp-trig.o: mp-trig.h mp-binomial.h mp-cache.h mp-complex.h mp-misc.h
mp-zeta.o: mp-zeta.h db-cache.h mp-binomial.h mp-cache.h mp-complex.h mp-consts.h mp-trig.h

# binaries
bernoulli.o: mp-polylog.h mp-misc.h
bs.o: mp-zeta.h mp-misc.h
bsubn.o: mp-zeta.h mp-binomial.h
cache-fill.o: db-cache.h mp-zeta.h mp-misc.h
db-merge.o: db-cache.h mp-misc.h
dsubn.o: mp-zeta.h mp-binomial.h mp-misc.h
exact.o: mp-complex.h mp-misc.h mp-polylog.h
hurwitz.o: mp-complex.h mp-misc.h
mp_zeta_test.o: mp-zeta.h mp-binomial.h mp-complex.h mp-consts.h mp-gamma.h mp-misc.h mp-polylog.h mp-trig.h
picturez.o: mp-zeta.h
sho.o: mp-complex.h mp-hyper.h mp-trig.h
stieltjes.o: mp-zeta.h mp-binomial.h mp-misc.h

bs:	bs.o $(MPLIB)
	$(CC) -o bs $^ -lgmp -lgsl -lcblas -ldb -lm -lc

bsubn:	bsubn.o $(MPLIB)
	$(CC) -o bsubn $^ -lgmp -lgsl -lcblas -ldb -lm -lc

minima: minima.o
	$(CC) -o minima $< -lgsl -lcblas -lm -lc

mp_zeta_test:	mp_zeta_test.o $(MPLIB)
	$(CC) -o mp_zeta_test $^ -lgmp -lgsl -lcblas -ldb -lm -lc

$(MPLIB): $(MPOBJS)
	ar cru $(MPLIB) $^

clean:
	rm -f core tmp junk glop a.out *.o

realclean:  clean
	rm -f $(EXES) *.png *.jpeg *.tex *.dat *.ps *.eps

dist:
	tar -zcvf mp.tgz *.c *.h *.gplot *.txt bzero.dat Makefile .cvsignore

.c.o:
	$(CC) -c -g -O2 -Wall $<

.o:
	$(CC) -g -o $* $^ -lgmp -ldb -lm
