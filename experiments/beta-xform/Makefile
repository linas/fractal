
CC = cc -Wall -g -O2 $(INCLUDES)

GENDIR = ../../generate
BRAT = $(GENDIR)/brat.o $(GENDIR)/util.o

FUNCDIR = ../../tools/inc
FUNC=../../tools/lib/libfunc.a
LAPACK=/usr/lib/x86_64-linux-gnu/liblapacke.a \
       /usr/lib/x86_64-linux-gnu/liblapack.a \
       /usr/lib/x86_64-linux-gnu/blas/libblas.a \
       /usr/lib/gcc/x86_64-linux-gnu/12/libgfortran.a \
       /usr/lib/x86_64-linux-gnu/libquadmath.so.0

all: bernie bigbern remap unbernie reigen skew fp-diagram \
	carry-sum rexfer treefn treestep sampler julia julie psi \
	matrix psieigen psiblas nilpot bitshow complex psiunitary \
	bergman psifp midpo psidist psidelt psileft moments psimom \
	psivec psimid psimult irred zero tribo midmap fiter expan \
	orbit fpcoef poincare subdiag beta witt deguess almost \
	alldisk betadisk momenta zetadisk normie irred-viz irred-fraction \
	irred-gf irred-disk irred-midpoint irred-allowed irred-tacky \
	irred-height irred-tree irred-seek eventually

.C.o:
	$(CC) -c $< -std=gnu++11 -I $(GENDIR) -I $(FUNCDIR)

.c.o:
	$(CC) -c $< -I $(FUNCDIR)

.o:
	$(CC) -o $@ $^ -lgmp -lpthread -lm -lstdc++

bernie.o: bernie.C $(FUNCDIR)/bitops.h
bigbern.o: bigbern.C
remap.o: remap.c
reigen.o: reigen.c
skew.o: skew.c
unbernie.o: unbernie.C
fp-diagram.o: fp-diagram.C
carry-sum.o: carry-sum.C
rexfer.o: rexfer.c
treefn.o: treefn.c
treestep.o: treestep.C
sampler.o: sampler.c
julia.o: julia.c
julie.o: julie.C
psi.o: psi.c
matrix.o: matrix.C psi.c
psieigen.o: psieigen.c psi.c psibig.c
psiblas.o: psiblas.c psi.c psibig.c
nilpot.o: nilpot.c
bitshow.o: bitshow.c $(FUNCDIR)/bitops.h
complex.o: complex.c
psiunitary.o: psiunitary.c psi.c psibig.c psifp.c
bergman.o: bergman.C psi.c psibig.c psifp.c
psifp.o: psifp.c psi.c psibig.c
midpo.o: midpo.c psi.c psibig.c
psidist.o: psidist.c psi.c psibig.c
psidelt.o: psidelt.C psi.c psibig.c
psileft.o: psileft.c psi.c psibig.c
moments.o: moments.c
psimom.o: psimom.c psileft.c psi.c psibig.c
psivec.o: psivec.c psileft.c psi.c psibig.c
psimid.o: psimid.c psi.c psibig.c
psimult.o: psimult.c psi.c psibig.c
zero.o: zero.C
tribo.o: tribo.c
midmap.o: midmap.c
fiter.o: fiter.c psi.c
expan.o: expan.c psi.c psibig.c
orbit.o: orbit.c psi.c psibig.c
fpcoef.o: fpcoef.c psi.c psibig.c expan.c
	$(CC) -c $< -I $(FUNCDIR) -fopenmp

poincare.o: poincare.C
subdiag.o: subdiag.C psi.c psibig.c
beta.o: beta.c
witt.o: witt.c
betadisk.o: betadisk.C psibig.c
zetadisk.o: zetadisk.C psibig.c
normie.o: normie.C psibig.c
irred.o: irred.c irred-gold.c
irred-viz.o: irred-viz.c irred-gold.c
irred-fraction.o: irred-fraction.c irred-gold.c
irred-gf.o: irred-gf.C irred-gold.c
irred-disk.o: irred-disk.C irred-gold.c
irred-midpoint.o: irred-midpoint.c irred-gold.c
irred-allowed.o: irred-allowed.c irred-gold.c
irred-tacky.o: irred-tacky.c irred-gold.c
irred-height.o: irred-height.c irred-gold.c
irred-tree.o: irred-tree.c irred-gold.c
irred-seek.o: irred-seek.c irred-gold.c
eventually.o: eventually.c

bernie: bernie.o $(BRAT) $(FUNC)
bigbern: bigbern.o $(BRAT)
remap: remap.o
reigen: reigen.o
skew: skew.o
unbernie: unbernie.o $(BRAT)
fp-diagram: fp-diagram.o $(BRAT)
carry-sum: carry-sum.o $(BRAT)
rexfer: rexfer.o
treefn: treefn.o
treestep: treestep.o  $(BRAT)
sampler: sampler.o
julia: julia.o
julie: julie.o $(BRAT)
psi: psi.o
matrix: matrix.o $(BRAT)
psieigen: psieigen.o -lgsl -lgslcblas
psiblas: psiblas.o $(LAPACK)
nilpot: nilpot.o
bitshow: bitshow.o $(FUNC)
complex: complex.o
psiunitary: psiunitary.o -lgsl -lgslcblas
bergman: bergman.o -lgsl -lgslcblas # $(BRAT)
psifp: psifp.o -lgsl -lgslcblas
midpo: midpo.o -lgsl -lgslcblas
psidist: psidist.o
psidist: psidist.o
psidelt: psidelt.o $(BRAT)
psileft: psileft.o
moments: moments.o
psimom: psimom.o
psivec: psivec.o
psimid: psimid.o
psimult: psimult.o
irred: irred.o
zero: zero.o $(BRAT)
tribo: tribo.o
midmap: midmap.o
fiter: fiter.o
expan: expan.o
orbit: orbit.o
fpcoef: fpcoef.o
	$(CC) -o $@ $^ -lgmp -lgomp -lm

poincare: poincare.o $(BRAT)
subdiag: subdiag.o $(BRAT)
beta: beta.o
witt: witt.o
deguess: deguess.o
almost: almost.o
alldisk: alldisk.o $(BRAT)
betadisk: betadisk.o $(BRAT)
zetadisk: zetadisk.o
momenta: momenta.o $(BRAT)
irred-gf: irred-gf.o $(BRAT)
irred-disk: irred-disk.o $(BRAT)
