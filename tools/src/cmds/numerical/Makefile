
#
# Makefile
#
# created by Linas Vepstas December 1992
#

BIN = $(TOP)/bin
LIB = $(TOP)/lib
INC = $(TOP)/inc

# CC = /bin/fcc -f2 -O -I $(INC) -DAIX221
# CC = /bin/cc -O -I $(INC) -DAIX315
CC = /usr/bin/cc -O2 -I $(INC) -DANSI_C

RM= rm -f

EXES = aucovar autocorr correlation distrib dvar filter spectrum station

all: $(EXES)

clean:
	rm -f *.o core tmp junk a.out

realclean:	clean
	rm -f $(EXES)

# --------------------------------------------------------------------
aucovar:	aucovar.c $(LIB)/fileio.o $(LIB)/var.o $(INC)/fileio.h $(INC)/var.h
	$(CC) -O -o aucovar aucovar.c $(LIB)/fileio.o $(LIB)/var.o -lm
	cp aucovar $(BIN)/aucovar
	rm -f aucovar.o aucovar
	touch aucovar
  
autocorr:	autocorr.c $(LIB)/fileio.o $(LIB)/var.o $(INC)/fileio.h $(INC)/var.h
	$(CC) -O -o autocorr autocorr.c $(LIB)/fileio.o $(LIB)/var.o -lm
	cp autocorr $(BIN)/autocorr
	rm -f autocorr.o autocorr
	touch autocorr
  
correlation:	correlation.c $(LIB)/fileio.o $(LIB)/var.o $(INC)/fileio.h $(INC)/var.h
	$(CC) -O -o correlation correlation.c $(LIB)/fileio.o $(LIB)/var.o -lm
	cp correlation $(BIN)/correlation
	rm -f correlation.o correlation
	touch correlation
  
dvar:	dvar.c $(LIB)/fileio.o $(LIB)/scale.o $(INC)/fileio.h $(INC)/scale.h
	$(CC) -O -o dvar dvar.c $(LIB)/fileio.o $(LIB)/scale.o -lm
	cp dvar $(BIN)/dvar
	rm -f dvar.o dvar
	touch dvar
  

distrib:	distrib.c
	$(CC) -o distrib distrib.c -lm
	cp distrib $(BIN)/distrib
	rm -f distrib
	touch distrib
  

filter:	filter.c $(LIB)/fileio.o $(INC)/fileio.h
	$(CC) -O -o filter filter.c $(LIB)/fileio.o
	cp filter $(BIN)/filter
	rm -f filter.o filter
	touch filter

spectrum:	spectrum.c $(LIB)/fileio.o $(LIB)/spectra.o $(INC)/fileio.h $(INC)/spectra.h
	$(CC) -O -o spectrum spectrum.c $(LIB)/fileio.o $(LIB)/spectra.o -lm
	cp spectrum $(BIN)/spectrum
	rm -f spectrum.o spectrum
	touch spectrum
  
station:	station.c $(LIB)/fileio.o $(LIB)/var.o $(INC)/fileio.h $(INC)/var.h
	$(CC) -O -o station station.c $(LIB)/fileio.o $(LIB)/var.o -lm
	cp station $(BIN)/station
	rm -f station.o station
	touch station
  
