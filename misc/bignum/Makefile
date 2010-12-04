
# build Gnu MP code
#
# CC = cc -pg
CC = cc

ANANT = anant

EXES= bernoulli bound bs bsubn cache-fill cont_frac \
		count-prime db-merge db-prt \
      dsubn entropy exact gkw-diag hurwitz hurwitz-euler \
		minima misc-test mp_quest_test mp_zeta_test picturez \
      poly-perf qint sawtooth sho stieltjes trunc \
		walsh yarh zero-cross zeta-zeros

MPLIB=anant/libanant.a

all:  $(MPLIB) $(EXES)
	cd anant
	make

bernoulli:	bernoulli.o $(MPLIB)
bound:	bound.o
cache-fill:	cache-fill.o $(MPLIB)
count-prime:	count-prime.o $(MPLIB)
db-merge:	db-merge.o $(MPLIB)
db-prt:	db-prt.o $(MPLIB)
cont_frac:	cont_frac.o
dsubn:	dsubn.o $(MPLIB)
entropy:	entropy.o $(MPLIB)
exact:	exact.o $(MPLIB)
gkw-diag:	gkw-diag.o $(MPLIB)
hurwitz:	hurwitz.o $(MPLIB)
hurwitz-euler:	hurwitz-euler.o $(MPLIB)
misc-test:	misc-test.o $(MPLIB)
mp_quest_test:	mp_quest_test.o $(MPLIB)
mp_zeta_test:	mp_zeta_test.o $(MPLIB)
picturez:	picturez.o $(MPLIB)
poly-perf:	poly-perf.o $(MPLIB)
qint:	qint.o $(MPLIB)
sawtooth:	sawtooth.o $(MPLIB)
sho:	sho.o $(MPLIB)
stieltjes:	stieltjes.o $(MPLIB)
trunc:	trunc.o
walsh:	walsh.o $(MPLIB)
yarh:	yarh.o $(MPLIB)
zero-cross:	zero-cross.o
zeta-zeros:	zeta-zeros.o $(MPLIB)


# binaries
bernoulli.o: $(ANANT)/mp-complex.h $(ANANT)/mp-polylog.h $(ANANT)/mp-misc.h
bs.o: $(ANANT)/mp-zeta.h $(ANANT)/mp-misc.h
bsubn.o: $(ANANT)/mp-zeta.h $(ANANT)/mp-binomial.h
cache-fill.o: $(ANANT)/db-cache.h $(ANANT)/mp-zeta.h $(ANANT)/mp-misc.h
count-prime.o: $(ANANT)/mp-quest.h
db-merge.o: $(ANANT)/db-cache.h $(ANANT)/mp-misc.h
dsubn.o: $(ANANT)/mp-zeta.h $(ANANT)/mp-binomial.h $(ANANT)/mp-misc.h
exact.o: $(ANANT)/mp-complex.h $(ANANT)/mp-misc.h $(ANANT)/mp-polylog.h
gkw-diag.o: $(ANANT)/mp-gkw.h
hurwitz.o: $(ANANT)/mp-complex.h $(ANANT)/mp-misc.h
hurwitz-euler.o: $(ANANT)/mp-complex.h $(ANANT)/mp-misc.h $(ANANT)/mp-trig.h
mp_zeta_test.o: $(ANANT)/mp-zeta.h $(ANANT)/mp-binomial.h $(ANANT)/mp-complex.h $(ANANT)/mp-consts.h $(ANANT)/mp-gamma.h $(ANANT)/mp-misc.h $(ANANT)/mp-polylog.h $(ANANT)/mp-trig.h
picturez.o: $(ANANT)/mp-zeta.h
qint.o: $(ANANT)/mp-quest.h
sho.o: $(ANANT)/mp-complex.h $(ANANT)/mp-hyper.h $(ANANT)/mp-trig.h
stieltjes.o: $(ANANT)/mp-zeta.h $(ANANT)/mp-binomial.h $(ANANT)/mp-misc.h

bs:	bs.o $(MPLIB)
	$(CC) -o bs $^ -lgmp -lgsl -lgslcblas -ldb -lm -lc

bsubn:	bsubn.o $(MPLIB)
	$(CC) -o bsubn $^ -lgmp -lgsl -lgslcblas -ldb -lm -lc

minima: minima.o
	$(CC) -o minima $< -lgsl -lgslcblas -lm -lc

mp_zeta_test:	mp_zeta_test.o $(MPLIB)
	$(CC) -o mp_zeta_test $^ -lgmp -lgsl -lgslcblas -ldb -lm -lc

clean:
	rm -f core tmp junk glop a.out *.o

realclean:  clean
	rm -f $(EXES) *.png *.jpeg *.tex *.dat *.ps *.eps

dist:
	tar -zcvf mp.tgz *.c *.h *.gplot *.txt bzero.dat Makefile .cvsignore

.c.o:
	$(CC) -c -g -O2 -Wall $< -I $(ANANT)

.o:
	$(CC) -g -o $* $^ -lgmp -ldb -lm
