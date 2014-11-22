############### Symbol Definitions ###################### 

TOP = /src/tools
INC = $(TOP)/inc
LIB = $(TOP)/lib

CC = /bin/fcc -f2 -O -I $(INC)

all: testspec

################## Dependencies and compile/link #########

testspec: testspec.c $(LIB)/stoch.o $(LIB)/spectra.o $(INC)/stoch.h $(INC)/spectra.h
	$(CC) -a -o testspec testspec.c $(LIB)/stoch.o $(LIB)/spectra.o
	rm testspec.o

