
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

all: sampler psi \
	matrix psieigen psiblas psiunitary \
	bergman psifp midpo psidist psidelt psileft moments psimom \
	psivec psimid psimult midmap fiter expan \
	orbit fpcoef almost deguess subdiag \
	alldisk betadisk momenta normie

.C.o:
	$(CC) -c $< -std=gnu++11 -I $(GENDIR) -I $(FUNCDIR)

.c.o:
	$(CC) -c $< -I $(FUNCDIR)

.o:
	$(CC) -o $@ $^ -lgmp -lpthread -lm -lstdc++

sampler.o: sampler.c
psi.o: psi.c
matrix.o: matrix.C psi.c
psieigen.o: psieigen.c psi.c psibig.c
psiblas.o: psiblas.c psi.c psibig.c
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
midmap.o: midmap.c
fiter.o: fiter.c psi.c
expan.o: expan.c psi.c psibig.c
orbit.o: orbit.c psi.c psibig.c
fpcoef.o: fpcoef.c psi.c psibig.c expan.c
	$(CC) -c $< -I $(FUNCDIR) -fopenmp

betadisk.o: betadisk.C psibig.c
subdiag.o: subdiag.C psi.c psibig.c
normie.o: normie.C psibig.c

sampler: sampler.o
psi: psi.o
matrix: matrix.o $(BRAT)
psieigen: psieigen.o -lgsl -lgslcblas
psiblas: psiblas.o $(LAPACK)
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
midmap: midmap.o
fiter: fiter.o
expan: expan.o
orbit: orbit.o
fpcoef: fpcoef.o
	$(CC) -o $@ $^ -lgmp -lgomp -lm

deguess: deguess.o
almost: almost.o
alldisk: alldisk.o $(BRAT)
betadisk: betadisk.o $(BRAT)
subdiag: subdiag.o $(BRAT)
momenta: momenta.o $(BRAT)
