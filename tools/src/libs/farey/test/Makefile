############### Symbol Definitions ###################### 

TOP = ../../../..
INC = $(TOP)/inc
LIB = $(TOP)/lib

# fcc for IBM RT
# CC = /bin/fcc -f2 -O -I $(INC)

CC = /usr/bin/cc -O2 -I $(INC) 

EXES= zgap zgap-check zlist fairy xlist elist flist gap ginv glist sym imap \
	terms x zscale zinv newgap ftree

LISTS= fnums enums

all: 	$(EXES)

clean:
	rm -f core tmp junk glop a.out *.o

realclean:	clean
	rm -f $(EXES) $(ENUMS)

.C.o:
	$(CC) -c $<

.o:
	$(CC) -o $* $< $(LIB)/Farey.o $(LIB)/FareyTree.o -lm -lstdc++

fairy.o:  fairy.C $(LIB)/Farey.o $(INC)/Farey.h
ftree.o:  ftree.C $(LIB)/FareyTree.o $(INC)/FareyTree.h
newgap.o:  newgap.C $(LIB)/Farey.o $(INC)/Farey.h
################## Dependencies and compile/link #########
# elist
# 	(Jan 1994) excercise Farey Number converter
#
elist: elist.c $(LIB)/Farey.o $(INC)/Farey.h
	$(CC) -o elist elist.c $(LIB)/Farey.o -lm

# xlist
# 	(Jan 1994) excercise Farey Number converter
#
xlist: xlist.c $(LIB)/Farey.o $(INC)/Farey.h
	$(CC) -o xlist xlist.c $(LIB)/Farey.o -lm
	$(RM) xlist.o

integ: integ.c $(LIB)/Farey.o $(INC)/Farey.h
	$(CC) -o integ integ.c $(LIB)/Farey.o -lm
	$(RM) integ.o

# flist
# 	(Jan 1994) excercise Farey Number converter
#
flist: flist.c $(LIB)/Farey.o $(INC)/Farey.h
	$(CC) -o flist flist.c $(LIB)/Farey.o -lm

# glist
# 	(Jan 1994) excercise Farey Number converter
#
glist: glist.c $(LIB)/Farey.o $(INC)/Farey.h
	$(CC) -o glist glist.c $(LIB)/Farey.o -lm

tlist: tlist.c $(LIB)/Farey.o $(INC)/Farey.h
	$(CC) -o tlist tlist.c $(LIB)/Farey.o -lm

# terms
# 	(Jan 1994) excercise Farey Number converter
#
terms: terms.c $(LIB)/Farey.o $(INC)/Farey.h
	$(CC) -o terms terms.c $(LIB)/Farey.o -lm

# zinv
# 	(Jan 1994) excercise Farey Number converter
#
zinv: zinv.c $(LIB)/Farey.o $(INC)/Farey.h
	$(CC) -o zinv zinv.c $(LIB)/Farey.o -lm

# zscale
# 	(Jan 1994) excercise Farey Number converter
#
zscale: zscale.c $(LIB)/Farey.o $(INC)/Farey.h
	$(CC) -o zscale zscale.c $(LIB)/Farey.o -lm
	rm zscale.o

# ginv
# 	(Jan 1994) excercise Farey Number converter
#
ginv: ginv.c $(LIB)/Farey.o $(INC)/Farey.h
	$(CC) -o ginv ginv.c $(LIB)/Farey.o -lm


# sym
# 	(Jan 1994) excercise Farey Number converter
#
sym: sym.c $(LIB)/Farey.o $(INC)/Farey.h
	$(CC) -o sym sym.c $(LIB)/Farey.o -lm


# gap
# 	(Jan 1994) excercise Farey Number converter
#
gap: gap.c $(LIB)/Farey.o $(LIB)/Prime.o $(INC)/Farey.h $(INC)/Prime.h
	$(CC) -o gap gap.c $(LIB)/Farey.o $(LIB)/Prime.o -lm


# x
# 	(Jan 1994) excercise Farey Number converter
#
x: x.c $(LIB)/Farey.o $(INC)/Farey.h
	$(CC) -o x x.c $(LIB)/Farey.o -lm -lX11


############################################################

# ------------------------------------------------------------
gral:	integ
	integ 5000 > gral

# xnums.ps:
ps: gral
	filter gral x In | xgraph -t "Integral of Farey Map" 
	filter gral x IIn | xgraph -t "Double Integral of Farey Map" 

# ------------------------------------------------------------
xnums:	xlist
	xlist 2000 > xnums

# xnums.ps:
xnums.ps: xnums
	filter xnums i f | xgraph -t "Exp'd Continued Fractions"

# ------------------------------------------------------------
enums:
	elist 5000 > enums

enums.ps:
	filter enums i f | xgraph -t "Exp'd Continued Fractions"

# ------------------------------------------------------------
fnums:
	flist 5000 > fnums

fnums.ps:
	filter fnums i f | xgraph -t "Farey Numbers"

# ------------------------------------------------------------

tnums:
#	tlist 3000 0.01 > tnum0.01 
	tlist 1000 0.1 > tnum0.1 
#	tlist 1000 0.4 > tnum0.4 
#	tlist 1000 0.6 > tnum0.6 
#	tlist 1000 0.8 > tnum0.8
	tlist 1000 1.0 > tnum1.0

tnums.ps:
#	filter tnum0.01 i f | xgraph -t "Farey Number with t=0.01" &
	filter tnum0.1 i f | xgraph -t "Farey Number with t=0.1" &
#	filter tnum0.4 i f | xgraph -t "Farey Number with t=0.4" &
#	filter tnum0.6 i f | xgraph -t "Farey Number with t=0.6" & 
#	filter tnum0.8 i f | xgraph -t "Farey Number with t=0.8" &
	filter tnum1.0 i f | xgraph -t "Farey Number with t=1.0" &

# ------------------------------------------------------------

gnums:
	glist 3000 0.01 > gnum0.01 
	glist 1000 0.1 > gnum0.1 
	glist 1000 0.4 > gnum0.4 
	glist 1000 0.6 > gnum0.6 
	glist 1000 0.69314718 > gnum.ln2 
	glist 1000 0.8 > gnum0.8
	glist 1000 1.0 > gnum1.0

gnums.ps:
	filter gnum0.01 i f | xgraph -t "Farey Number with t=0.01" &
	filter gnum0.1 i f | xgraph -t "Farey Number with t=0.1" &
	filter gnum0.4 i f | xgraph -t "Farey Number with t=0.4" &
	filter gnum0.6 i f | xgraph -t "Farey Number with t=0.6" & 
	filter gnum.ln2 i f | xgraph -t "Farey Number with t=ln2" &
	filter gnum0.8 i f | xgraph -t "Farey Number with t=0.8" &
	filter gnum1.0 i f | xgraph -t "Farey Number with t=1.0" &

# ------------------------------------------------------------

znums:
#	zlist 3000 -1.2 > znum-1.2
#	zlist 3000 -1.0 > znum-1.0
#	zlist 3000 -0.8 > znum-0.8
#	zlist 3000 -0.5 > znum-0.5
	zlist 3000 -0.3 > znum-0.3
	zlist 3000 -0.25 > znum-0.25
#	zlist 3000 -0.2 > znum-0.2
#	zlist 3000 -0.03 > znum-0.03
#	zlist 3000 0.03 > znum0.03
#	zlist 3000 0.2 > znum0.2
#	zlist 3000 0.5 > znum0.5
#	zlist 3000 0.8 > znum0.8
#	zlist 3000 1.0 > znum1.0
#	zlist 3000 1.2 > znum1.2
#	zlist 3000 1.6 > znum1.6
#	zlist 3000 2.0 > znum2.0
#	zlist 3000 4.0 > znum4.0
#	zlist 3000 10.0 > znum10.0

znums.ps:
#	filter znum-1.2 i f | xgraph -t "Real Number with z= -1.2" &
#	filter znum-1.0 i f | xgraph -t "Real Number with z= -1.0" &
#	filter znum-0.8 i f | xgraph -t "Real Number with z= -0.8" &
#	filter znum-0.5 i f | xgraph -t "Real Number with z= -0.5" &
	filter znum-0.3 i f | xgraph -t "Real Number with z= -0.3" &
	filter znum-0.25 i f | xgraph -t "Real Number with z= -0.25" &
#	filter znum-0.2 i f | xgraph -t "Real Number with z= -0.2" &
#	filter znum-0.03 i f | xgraph -t "Real Number with z= -0.03" &
#	filter znum0.03 i f | xgraph -t "Real Number with z=0.03" &
#	filter znum0.2 i f | xgraph -t "Real Number with z=0.2" &
#	filter znum0.5 i f | xgraph -t "Real Number with z=0.5" &
#	filter znum0.8 i f | xgraph -t "Real Number with z=0.8" &
#	filter znum1.0 i f | xgraph -t "Real Number with z=1.0" &
#	filter znum1.2 i f | xgraph -t "Real Number with z=1.2" &
#	filter znum1.6 i f | xgraph -t "Real Number with z=1.6" &
#	filter znum2.0 i f | xgraph -t "Real Number with z=2.0" &
#	filter znum4.0 i f | xgraph -t "Real Number with z=4.0" &
#	filter znum10.0 i f | xgraph -t "Real Number with z=10.0" &

# ------------------------------------------------------------

zinvs:
	zinv 3000 0.03 > zinv0.03
	zinv 3000 0.2 > zinv0.2
	zinv 3000 0.5 > zinv0.5
	zinv 3000 2.0 > zinv2.0
	zinv 3000 4.0 > zinv4.0
	zinv 3000 10.0 > zinv10.0

zinvs.ps:
	filter zinv0.03 i f | xgraph -t "Reciprocal Real with z=0.03" &
	filter zinv0.2 i f | xgraph -t "Reciprocal Real with z=0.2" &
	filter zinv0.5 i f | xgraph -t "Reciprocal Real with z=0.5" &
	filter zinv2.0 i f | xgraph -t "Reciprocal Real with z=2.0" &
	filter zinv4.0 i f | xgraph -t "Reciprocal Real with z=4.0" &
	filter zinv10.0 i f | xgraph -t "Reciprocal Real with z=10.0" &

# ------------------------------------------------------------

station.ps:
	cat station.z | xgraph -t "Pure Gaussian Noise Stationarity" -x "loz10 (t2-t1)" -y "loz10 Var (X(t2)-X(t1))"
	cat station.b | xgraph -t "Stationarity of Brownian Motion" -x "loz10 (t2-t1)" -y "loz10 Var (X(t2) - X(t1))"

# ------------------------------------------------------------

smooth:
	spectrum zauss.z i z > smooth.z
	spectrum brown.b i z > smooth.b


smooth.ps:
	cat smooth.z | xgraph -t "Smooth Loz Power Spectrum of Gaussian Noise" -x "loz10 (omeza)" -y "loz10 Power"
	cat smooth.b | xgraph -t "Smooth Loz Power Spectrum of Brownian Motion" -x "loz10 (omeza)" -y "loz10 Power"

