############### Symbol Definitions ###################### 

INCLUDES =	-I../tools/inc
LIBS =		/usr/lib

CC = cc -Wall -g -O2 $(INCLUDES)

OLD_EXES =  brat-gapper stalk stalkmov mstop migrate measure offset \
        lyapunov next squige age phase orig whack \
        circout circin circmov brat \
        abs clamp takelog reclamp dump renorm \
        draw dx image movie paste stats sum \
        angle curl diff div card mandel morph \
        moire automatic wind winds ray \
        manvert cutoff brat draw dx image

NEW_EXES = affine alpha baker circle cutoff divisor elliptic \
        erdos euler-q fdist gamma hermite hurwitz ising ising-moment mand-flow \
        mp_zeta mobius plouffe polylog q-exp sho swap takagi totient zeta \
        ray automatic radius

all:    $(NEW_EXES)

dist:
	tar -zcvf gen.tgz *.c *.C *.h Makefile .cvsignore cvt

clean:
	rm -f x j.c a.out  *.o glop* junk* *.flo  *.txt 

realclean:
	rm -f *.gif *.jpeg *.png $(OLD_EXES) $(NEW_EXES)

.c.o:
	$(CC) -c $<

.C.o:
	$(CC) -c $<

.o:
	$(CC) -o $@ $^ -lm -lstdc++

coord-xforms.o:	coord-xforms.h
image.o: image.h util.h
gap-tongue.o: image.h util.h
man.o: image.h util.h
opers.o: opers.h
renorm.o: opers.h
util.o:	util.h

affine.o: affine.C
alpha.o: alpha.C
baker.o: baker.C
brat-beigen.o: brat-beigen.C
brat-gap-hair.o: brat-gap-hair.C
brat-gapper.o: brat-gapper.C
circle.o: circle.C brat.h
cutoff.o: cutoff.C
divisor.o: divisor.C brat.h
euler-q.o: euler-q.C
fdist.o: fdist.C
gamma.o: gamma.C
hermite.o: hermite.C
hurwitz.o: hurwitz.C brat.h
ising.o: ising.C
mobius.o: mobius.C
mp_zeta.o: mp_zeta.C
plouffe.o: plouffe.C
polylog.o: polylog.C brat.h
q-exp.o: q-exp.C brat.h
swap.o: swap.C
sho.o: sho.C brat.h coord-xforms.h
takagi.o: takagi.C brat.h
totient.o: totient.C brat.h
zeta.o: zeta.C brat.h

FUNC=../tools/lib/libfunc.a
MP=../misc/bignum/libzmp.a
GMP=-lgmp
GSL=-lgsl -lcblas

################## Dependencies and compile/link #########
#
affine: brat.o affine.o util.o
alpha: brat.o alpha.o util.o $(GSL)
baker: brat.o baker.o util.o
brat-beigen: brat.o brat-beigen.o util.o
brat-gap-hair: brat.o brat-gap-hair.o util.o $(FUNC)
brat-gapper: brat.o brat-gapper.o util.o $(FUNC)
circle: brat.o circle.o util.o
cutoff: brat.o cutoff.o coord-xforms.o util.o 
divisor: brat.o divisor.o util.o $(FUNC)
elliptic: brat.o elliptic.o util.o
erdos: brat.o erdos.o coord-xforms.o util.o $(FUNC)
euler-q: brat.o euler-q.o coord-xforms.o util.o $(FUNC)
fdist: brat.o fdist.o util.o $(FUNC)
gamma: brat.o gamma.o util.o coord-xforms.o $(MP) $(GMP) -ldb
hermite: brat.o hermite.o util.o 
hurwitz: brat.o hurwitz.o util.o $(MP) $(GMP) -ldb
ising: brat.o ising.o util.o $(FUNC)
ising-moment: brat.o ising-moment.o util.o $(FUNC)
mand-flow: brat.o mand-flow.o util.o
mobius: brat.o mobius.o util.o $(FUNC)
mp_zeta: brat.o mp_zeta.o util.o $(MP) $(GMP) -ldb
plouffe: brat.o plouffe.o util.o
polylog: brat.o polylog.o util.o $(MP) $(GMP) -ldb
q-exp: brat.o q-exp.o util.o
swap: brat.o swap.o util.o
sho: brat.o sho.o util.o coord-xforms.o $(MP) $(GMP) -ldb
takagi: brat.o takagi.o util.o
totient: brat.o totient.o util.o $(FUNC)
zeta: brat.o zeta.o util.o $(GSL)


manvert stalk stalkmov mstop migrate measure offset lyapunov next squige age phase orig whack circout circin circmov brat: brat.o  util.o $(FAREY)
	$(CC) -o brat brat.o util.o -lm -lstdc++
	ln -f brat manvert
	ln -f brat stalk
	ln -f brat mstop
	ln -f brat measure
	ln -f brat offset
	ln -f brat lyapunov
	ln -f brat next
	ln -f brat squige
	ln -f brat age
	ln -f brat migrate
	ln -f brat phase
	ln -f brat orig
	ln -f brat whack
	ln -f brat circout
	ln -f brat circin
	ln -f brat circmov
	ln -f brat stalkmov
	ln -f brat wind
	ln -f brat winds
	ln -f brat morph

card: card.c  util.o
	$(CC) -o card card.c util.o -lm

draw: draw.c  util.o
	$(CC) -o draw -I /home/linas/linas/fractal-tools/inc draw.c util.o /home/linas/linas/fractal-tools/lib/Farey.o -lm
#	$(CC) -o draw -I /src/tools/inc draw.c util.o -lm
#	$(CC) -o draw -I /src/tools/inc draw.c util.o /src/tools/lib/Farey.o -lm

sum diff paste angle curl div: sum.o  util.o opers.o
	$(CC) -o sum sum.o util.o opers.o -lm
	ln -f sum diff
	ln -f sum paste
	ln -f sum angle
	ln -f sum curl
	ln -f sum div

dx:	dx.c  util.o opers.o
	$(CC)  -o dx dx.c util.o opers.o -lm

abs clamp takelog reclamp dump renorm:	renorm.o  util.o opers.o
	$(CC)  -o renorm renorm.o opers.o util.o -lm
	ln -f renorm abs
	ln -f renorm clamp
	ln -f renorm takelog
	ln -f renorm reclamp
	ln -f renorm dump

stats:	stats.c  util.o opers.o
	$(CC)  -o stats stats.c opers.o util.o -lm

image: image.o  gap-tongue.o man.o util.o ../tools/lib/Farey.o
	$(CC) -o image $^ -lm -lstdc++
	ln -f image mandel
	ln -f image moire

movie: movie.o  man.o util.o
	$(CC) -o movie movie.o man.o util.o -lm

radius: radius.o 

automatic: radius
	ln -f radius automatic

ray:	ray.o
