
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
CC = /usr/bin/cc -O2 -I $(INC)

all: chop concat
	cd image; make all

chop:	chop.c 
	$(CC) -o chop chop.c 
	cp chop $(BIN)/chop
	rm chop
	touch chop

concat:	concat.c 
	$(CC) -o concat concat.c 
	cp concat $(BIN)/concat
	rm concat
	touch concat

clean:
	cd numerical; make clean
	rm -f junk tmp core a.out *.o glop

realclean: clean
	cd numerical; make realclean
	rm -f concat chop
	
