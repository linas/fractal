############### Symbol Definitions ###################### 

TOP = /src/tools
INC = $(TOP)/inc
LIB = $(TOP)/lib

CC = /bin/fcc -f2 -O -I $(INC)

all: testrand teststoch testgauss

# ------------------------------------------------------------

clean:
	rm -f core a.out *.o tmp junk glop

realclean:	clean
	rm -f testrand teststoch testgauss *.b *.g *.ps


################## Dependencies and compile/link #########
#
# testrand
#         (Dec 1988) excercises the random number generator
testrand: testrand.c 
	$(CC) -o testrand testrand.c 

# teststoch
# 	(Dec 1992) excercise gaussian noise generator
#
teststoch: teststoch.c $(LIB)/stoch.o $(INC)/stoch.h
	$(CC) -o teststoch teststoch.c $(LIB)/stoch.o -lm
	rm teststoch.o


# testgauss
# 	(Dec 1992) excercise gaussian gaussal noise generator
#
testgauss: testgauss.c $(LIB)/stoch.o $(INC)/stoch.h
	$(CC) -o testgauss testgauss.c $(LIB)/stoch.o
	rm testgauss.o

############################################################
# The following generate the various spectral analysis of the noise data
# Dec 1993

# ------------------------------------------------------------
dvar:
	dvar gauss.g i g > dvar.g
	dvar brown.b i g > dvar.b

dvar.ps:
	cat dvar.g | xgraph -t "Pure Gaussian Noise Scaling" -x "log10 sigma" -y "log10 RMS Derivative"
	cat dvar.b | xgraph -t "Scaling of Brownian Motion" -x "log10 sigma" -y "log10 RMS Derivative"

# ------------------------------------------------------------

covar:
	aucovar gauss.g i g > covar.g
	aucovar brown.b i g > covar.b


covar.ps:
	cat covar.g | xgraph -t "Pure Gaussian Noise AutoCovariance" -x "(t2-t1)" -y "CoV (X(t1) , X(t2))"
	cat covar.b | xgraph -t "AutoCovariance of Brownian Motion" -x "(t2-t1)" -y "CoV (X(t1) , X(t2))"

# ------------------------------------------------------------

corr:
	autocorr gauss.g i g > corr.g
	autocorr brown.b i g > corr.b


corr.ps:
	cat corr.g | xgraph -t "Pure Gaussian Noise AutoCorrelation" -x "(t2-t1)" -y "Corr (X(t1) , X(t2))"
	cat corr.b | xgraph -t "AutoCorrelation of Brownian Motion" -x "(t2-t1)" -y "Corr (X(t1) , X(t2))"

# ------------------------------------------------------------

station:
	station gauss.g i g > station.g
	station brown.b i g > station.b


station.ps:
	cat station.g | xgraph -t "Pure Gaussian Noise Stationarity" -x "log10 (t2-t1)" -y "log10 Var (X(t2)-X(t1))"
	cat station.b | xgraph -t "Stationarity of Brownian Motion" -x "log10 (t2-t1)" -y "log10 Var (X(t2) - X(t1))"

# ------------------------------------------------------------

smooth:
	spectrum gauss.g i g > smooth.g
	spectrum brown.b i g > smooth.b


smooth.ps:
	cat smooth.g | xgraph -t "Smooth Log Power Spectrum of Gaussian Noise" -x "log10 (omega)" -y "log10 Power"
	cat smooth.b | xgraph -t "Smooth Log Power Spectrum of Brownian Motion" -x "log10 (omega)" -y "log10 Power"

