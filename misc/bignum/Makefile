
# build Gnu MP code
#


EXES= bound bs bsubn cache-fill cont_frac db-merge db-prt \
      dsubn minima mp_zeta_test picturez \
      stieltjes trunc zero-cross

all:  $(EXES)

MPOBJS= db-cache.o mp-binomial.o mp-cache.o mp-consts.o mp-misc.o mp-trig.o

bound:	bound.o
cache-fill:	cache-fill.o mp_zeta.o $(MPOBJS)
db-merge:	db-merge.o $(MPOBJS)
db-prt:	db-prt.o $(MPOBJS)
cont_frac:	cont_frac.o
dsubn:	dsubn.o mp_zeta.o $(MPOBJS)
mp_zeta_test:	mp_zeta_test.o mp_zeta.o $(MPOBJS)
picturez:	picturez.o mp_zeta.o $(MPOBJS)
stieltjes:	stieltjes.o mp_zeta.o $(MPOBJS)
trunc:	trunc.o
zero-cross:	zero-cross.o

# libarary objects
db-cache.o: db-cache.h
mp-binomial.o: mp-binomial.h mp-cache.h mp-trig.h
mp-cache.o: mp-cache.h
mp-consts.o: mp-consts.h mp-binomial.h mp-trig.h
mp-misc.o: mp-misc.h
mp-trig.o: mp-trig.h mp-cache.h mp-binomial.h
mp_zeta.o: mp_zeta.h db-cache.h mp-cache.h

# binaries
bs.o: mp_zeta.h mp-misc.h
bsubn.o: mp_zeta.h mp-binomial.h
cache-fill.o: db-cache.h mp_zeta.h mp-misc.h
db-merge.o: db-cache.h mp-misc.h
dsubn.o: mp_zeta.h mp-binomial.h mp-misc.h
mp_zeta_test.o: mp_zeta.h mp-consts.h
stieltjes.o: mp_zeta.h mp-binomial.h mp-misc.h

bs:	bs.o mp_zeta.o $(MPOBJS)
	$(CC) -o bs $^ -lgmp -lgsl -lcblas -ldb -lm -lc

bsubn:	bsubn.o mp_zeta.o $(MPOBJS)
	$(CC) -o bsubn $^ -lgmp -lgsl -lcblas -ldb -lm -lc

minima: minima.o
	$(CC) -o minima $< -lgsl -lcblas -lm -lc

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
