############### Symbol Definitions ###################### 

INCLUDES =	/usr/include
LIBS =		/usr/lib
XLIB =		$(LIBS)/libX11.a 

# CC = /bin/cc -O -O2 -I$(INCLUDES)
CC = /usr/bin/cc -O -O2 -I$(INCLUDES)

all:	flo2mtv

################## Dependencies and compile/link #########

# blit
#	(Feb 1989) like moire, but an attempt to put image into 
#		root window.  Works but crashes. (Bug in server)
#	(Mar 1989) fixed with update 1713
#	(April 1989) converted to a parallel algorithm and renamed from
#		background to parallel.  Links to parablit.o instead of
#		bitblit.o
#	(March 1990) cast into four general purpose tools

testblit: testblit.c blit.o $(XLIB)
	$(CC) -a -o testblit testblit.c blit.o $(XLIB)

blit.o: blit.c
	$(CC) -c blit.c 


# addblit
#	(March 1990) gutted xpic to create this module.  It basically
#	just opens a window and adds stuff to it.

testadd: testadd.c globals.o addblit.o blit.o $(XLIB)
	$(CC) -a -o testadd testadd.c globals.o addblit.o blit.o $(XLIB)

addblit.o: addblit.c
	$(CC) -c addblit.c 

globals.o: globals.c
	$(CC) -c globals.c 


flo2mtv:	flo2mtv.c
	$(CC) flo2mtv.c -o flo2mtv

clean: 
	rm -f *.o

realclean: clean
	rm -f flo2mtv
