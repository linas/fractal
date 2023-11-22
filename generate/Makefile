############### Symbol Definitions ######################

INCLUDES = -I../tools/inc \
           -I /usr/local/include/anant
LIBS =		/usr/lib

CC = cc -Wall -g -O2 $(INCLUDES)
CXX = cc -std=gnu++11 -Wall -g -O2 $(INCLUDES)

OLD_EXES =  brat-gapper stalk stalkmov mstop migrate measure offset \
        lyapunov next squige age phase orig whack \
        circout circin circmov brat \
        draw dx image movie paste stats sum \
        angle curl diff div card mandel morph \
        moire automatic wind winds ray \
        manvert cutoff brat draw dx image

NEW_EXES = affine alpha baker chirikov \
        circle-mom cutoff divisor elliptic \
        erdos euler-q fdist gamma gauss-red gkd gkw gkw-integrand \
        haar hardy hermite hurwitz ising \
        ising-moment mand-flow \
        mp_zeta mobius plouffe q-exp sho swap takagi totient zeta \
        ray automatic radius

TOOLS = renorm

all:    $(NEW_EXES) $(TOOLS)

dist:
	tar -zcvf gen.tgz *.c *.C *.h Makefile .cvsignore cvt

clean:
	rm -f x j.c a.out  *.o glop* junk* *.flo  *.txt

realclean:
	rm -f *.gif *.jpeg *.png $(OLD_EXES) $(NEW_EXES)

.c.o:
	$(CC) -c $<

.C.o:
	$(CXX) -c $<

.o:
	$(CC) -o $@ $^ -lpthread -lm -lstdc++

coord-xforms.o:	coord-xforms.h
image.o: image.h util.h
gap-tongue.o: image.h util.h
man.o: image.h util.h
opers.o: opers.h
renorm.o: opers.h
util.o:	util.h

affine.o: affine.C brat.h
alpha.o: alpha.C brat.h
baker.o: baker.C brat.h
brat-beigen.o: brat-beigen.C
brat-gap-hair.o: brat-gap-hair.C
brat-gapper.o: brat-gapper.C
chirikov.o: chirikov.C brat.h
circle-mom.o: circle-mom.C brat.h
cutoff.o: cutoff.C
divisor.o: divisor.C brat.h
elliptic.o: elliptic.C brat.h
erdos.o: erdos.C brat.h
euler-q.o: euler-q.C brat.h
fdist.o: fdist.C brat.h
gamma.o: gamma.C brat.h
gauss-red.o: gauss-red.C brat.h
gkd.o: gkd.C brat.h
gkw.o: gkw.C brat.h
gkw-integrand.o: gkw-integrand.C brat.h
haar.o: haar.C brat.h
hardy.o: hardy.C brat.h
hermite.o: hermite.C brat.h
hurwitz.o: hurwitz.C brat.h
ising.o: ising.C brat.h
mobius.o: mobius.C brat.h
mp_zeta.o: mp_zeta.C
plouffe.o: plouffe.C brat.h
q-exp.o: q-exp.C brat.h
swap.o: swap.C brat.h
sho.o: sho.C brat.h coord-xforms.h
takagi.o: takagi.C brat.h
totient.o: totient.C brat.h
zeta.o: zeta.C brat.h

FUNC=../tools/lib/libfunc.a -lpthread
# MP=../experiments/anant-git/src/libanant.a -ldb -lpthread
MP=-lanant -ldb -lpthread
GMP=-lgmp
GSL=-lgsl -lgslcblas

################## Dependencies and compile/link #########
#
affine: brat.o affine.o util.o
alpha: brat.o alpha.o util.o $(GSL)
baker: brat.o baker.o util.o
brat-beigen: brat.o brat-beigen.o util.o
brat-gap-hair: brat.o brat-gap-hair.o util.o $(FUNC)
brat-gapper: brat.o brat-gapper.o util.o $(FUNC)
chirikov: brat.o chirikov.o util.o
circle-mom: brat.o circle-mom.o util.o
cutoff: brat.o cutoff.o coord-xforms.o util.o
divisor: brat.o divisor.o util.o $(FUNC)
elliptic: brat.o elliptic.o util.o
erdos: brat.o erdos.o coord-xforms.o util.o $(FUNC)
euler-q: brat.o euler-q.o coord-xforms.o util.o $(FUNC)
fdist: brat.o fdist.o util.o $(FUNC)
gamma: brat.o gamma.o util.o coord-xforms.o $(MP) $(GMP)
gauss-red: brat.o gauss-red.o util.o
gkd: brat.o gkd.o util.o
gkw: brat.o gkw.o util.o $(FUNC) $(MP) $(GMP)
gkw-integrand: brat.o gkw-integrand.o util.o $(FUNC) $(MP) $(GMP)
haar: brat.o haar.o util.o
hardy: brat.o hardy.o util.o $(FUNC)
hermite: brat.o hermite.o util.o
hurwitz: brat.o hurwitz.o util.o $(MP) $(GMP)
ising: brat.o ising.o util.o $(FUNC)
ising-moment: brat.o ising-moment.o util.o $(FUNC)
mand-flow: brat.o mand-flow.o util.o
mobius: brat.o mobius.o util.o $(FUNC)
mp_zeta: brat.o mp_zeta.o util.o $(MP) $(GMP)
plouffe: brat.o plouffe.o util.o
q-exp: brat.o q-exp.o util.o
swap: brat.o swap.o util.o
sho: brat.o sho.o util.o coord-xforms.o $(MP) $(GMP)
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

sum diff paste mask angle curl div: sum.o  util.o opers.o
	$(CC) -o sum sum.o util.o opers.o -lm
	ln -f sum diff
	ln -f sum paste
	ln -f sum mask
	ln -f sum angle
	ln -f sum curl
	ln -f sum div

dx:	dx.c  util.o opers.o
	$(CC)  -o dx dx.c util.o opers.o -lm

abs clamp takelog takeroot recip reclamp dump renorm:	renorm.o  util.o opers.o
	$(CC)  -o renorm renorm.o opers.o util.o -lm
	ln -f renorm abs
	ln -f renorm clamp
	ln -f renorm takelog
	ln -f renorm takeroot
	ln -f renorm recip
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
