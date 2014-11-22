

#
# HISTORY:
# Linas Vepstas November 1993

TOP = /usr/home/src/tools

BIN = $(TOP)/bin
LIB = $(TOP)/lib
INC = $(TOP)/inc

# CC = /bin/fcc -f2 -O -I $(INC) -DAIX221
# CC = /bin/cc -O -I $(INC) -DAIX315
CC = /usr/bin/cc -O2 -I $(INC) -DANSI_C

OBJS=  scale.o spectra.o var.o 

all: inc $(OBJS)

clean:
	rm -f core tmp junk a.out *.o

realclean:	clean
	rm -f inc


# ------------------------------------------------------------
inc:	scale.h spectra.h var.h
	touch inc
	cp scale.h $(INC)/scale.h
	cp spectra.h $(INC)/spectra.h
	cp var.h $(INC)/var.h

scale.o:	scale.c scale.h
	$(CC) -c scale.c
	cp scale.o $(LIB)/scale.o
	rm scale.o
	touch scale.o

spectra.o:	spectra.c spectra.h
	$(CC) -c spectra.c
	cp spectra.o $(LIB)/spectra.o
	rm spectra.o
	touch spectra.o

var.o:	var.c var.h
	$(CC) -c var.c
	cp var.o $(LIB)/var.o
	rm var.o
	touch var.o

