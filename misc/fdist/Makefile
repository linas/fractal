############### Symbol Definitions ###################### 

TOP = ../../tools
INC = $(TOP)/inc
LIB = $(TOP)/lib

CC = /usr/bin/cc -O2 -I $(INC) 

EXES= akk f-wavelet fanal fdist fourier ftag janal jdist \
      mink-meas rdist shift takagi-dist

all: 	$(EXES)

clean:
	rm -f core tmp junk glop a.out *.o

realclean:	clean
	rm -f $(EXES) *.png *.jpeg *.tex *.dat *.ps *.eps

dist:
	tar -zcvf fdist.tgz *.c *.C *.gplot Makefile .cvsignore

.C.o:
	$(CC) -c -g -O2 $<

.o:
	$(CC) -o $* $< -L$(LIB) -lfunc -lgsl -lcblas -lgmpxx -lgmp -lm -lstdc++ 


akk.o:  akk.C 
f-wavelet.o:  f-wavelet.C $(LIB)/libfunc.a $(INC)/FareyTree.h
fanal.o:  fanal.C $(LIB)/libfunc.a $(INC)/FareyTree.h
fdist.o:  fdist.C $(LIB)/libfunc.a $(INC)/FareyTree.h
fourier.o:  fourier.C $(LIB)/libfunc.a $(INC)/FareyTree.h
ftag.o:  ftag.C $(LIB)/libfunc.a
janal.o:  janal.C 
jdist.o:  jdist.C $(LIB)/libfunc.a $(INC)/FareyTree.h
mink-meas.o:  mink-meas.C $(LIB)/libfunc.a
rdist.o:  rdist.C $(LIB)/libfunc.a $(INC)/FareyTree.h
shift.o:  shift.C $(LIB)/libfunc.a
takagi-dist.o:  takagi-dist.C $(LIB)/libfunc.a
