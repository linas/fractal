o

############### Symbol Definitions ###################### 

INC = $(TOP)/inc
LIB = $(TOP)/lib

CC = /bin/fcc -f2 -O -I $(INC)

all: testfract fbm

################## Dependencies and compile/link #########
#

# testgauss
# 	(Dec 1992) excercise gaussian noise generator
#
testgauss: testgauss.c $(LIB)/stoch.o $(INC)/stoch.h
	$(CC) -o testgauss testgauss.c $(LIB)/stoch.o
	rm testgauss.o


# testfract
# 	(Dec 1992) excercise gaussian fractal noise generator
#
testfract: testfract.c $(LIB)/stoch.o $(INC)/stoch.h
	$(CC) -o testfract testfract.c $(LIB)/stoch.o
	rm testfract.o

# fbm
# 	(Dec 1992) excercise gaussian fractal noise generator
#
fbm: fbm.c $(LIB)/stoch.o $(INC)/stoch.h
	$(CC) -o fbm fbm.c $(LIB)/stoch.o
	rm fbm.o

# lacune
# 	(Dec 1992) excercise gaussian fractal noise generator
#
lacune: lacune.c $(LIB)/stoch.o $(INC)/stoch.h
	$(CC) -o lacune lacune.c $(LIB)/stoch.o
	rm lacune.o


##########################################################
#

# -------------------------------------------------------
# fract builds the raw data fractal noise files
fract:
	fbm 0.1 > fbm.1
	fbm 0.3 > fbm.3
	fbm 0.5 > fbm.5
	fbm 0.7 > fbm.7
	fbm 0.9 > fbm.9

# -------------------------------------------------------
#fract.ps sets up graphs of the raw data
fract.ps:
	filter fbm.1 i f |xgraph -t "H=0.1 Fractal Noise" -x "Sample Points"
	filter fbm.3 i f |xgraph -t "H=0.3 Fractal Noise" -x "Sample Points"
	filter fbm.5 i f |xgraph -t "H=0.5 Fractal Noise" -x "Sample Points"
	filter fbm.7 i f |xgraph -t "H=0.7 Fractal Noise" -x "Sample Points"
	filter fbm.9 i f |xgraph -t "H=0.9 Fractal Noise" -x "Sample Points"

scale:
	dvar fbm.1 i f > dvar.1
	dvar fbm.3 i f > dvar.3
	dvar fbm.5 i f > dvar.5
	dvar fbm.7 i f > dvar.7
	dvar fbm.9 i f > dvar.9

scale.ps:
	cat dvar.1 |xgraph -t "H=0.1 Fractal Noise Scaling" -x "log10 Sigma" -y "log10 RMS Derivative"
	cat dvar.3 |xgraph -t "H=0.3 Fractal Noise Scaling" -x "log10 Sigma" -y "log10 RMS Derivative"
	cat dvar.5 |xgraph -t "H=0.5 Fractal Noise Scaling" -x "log10 Sigma" -y "log10 RMS Derivative"
	cat dvar.7 |xgraph -t "H=0.7 Fractal Noise Scaling" -x "log10 Sigma" -y "log10 RMS Derivative"
	cat dvar.9 |xgraph -t "H=0.9 Fractal Noise Scaling" -x "log10 Sigma" -y "log10 RMS Derivative"


# -------------------------------------------------------
# compute & graph autocorrelations
auto:
	autocorr fbm.1 i f > auto.1
	autocorr fbm.3 i f > auto.3
	autocorr fbm.5 i f > auto.5
	autocorr fbm.7 i f > auto.7
	autocorr fbm.9 i f > auto.9

auto.ps:
	cat auto.1 |xgraph -t "H=0.1 AutoCorrelation" -x "(t1-t2)" -y "Corr (X(t1), X(t2))"
	cat auto.3 |xgraph -t "H=0.3 AutoCorrelation" -x "(t1-t2)" -y "Corr (X(t1), X(t2))"
	cat auto.5 |xgraph -t "H=0.5 AutoCorrelation" -x "(t1-t2)" -y "Corr (X(t1), X(t2))"
	cat auto.7 |xgraph -t "H=0.7 AutoCorrelation" -x "(t1-t2)" -y "Corr (X(t1), X(t2))"
	cat auto.9 |xgraph -t "H=0.9 AutoCorrelation" -x "(t1-t2)" -y "Corr (X(t1), X(t2))"

# -------------------------------------------------------
# compute & graph autocovariance

covar:
	aucovar fbm.1 i f > covar.1
	aucovar fbm.3 i f > covar.3
	aucovar fbm.5 i f > covar.5
	aucovar fbm.7 i f > covar.7
	aucovar fbm.9 i f > covar.9

covar.ps:
	cat covar.1 |xgraph -t "H=0.1 AutoCovariance" -x "(t1-t2)" -y "CoV (X(t1), X(t2))"
	cat covar.3 |xgraph -t "H=0.3 AutoCovariance" -x "(t1-t2)" -y "CoV (X(t1), X(t2))"
	cat covar.5 |xgraph -t "H=0.5 AutoCovariance" -x "(t1-t2)" -y "CoV (X(t1), X(t2))"
	cat covar.7 |xgraph -t "H=0.7 AutoCovariance" -x "(t1-t2)" -y "CoV (X(t1), X(t2))"
	cat covar.9 |xgraph -t "H=0.9 AutoCovariance" -x "(t1-t2)" -y "CoV (X(t1), X(t2))"


grspec:
	spectrum fbm.1 i f | xgraph -x "log10 freq" -y "log10 Power" -t "H=0.1 Fractal Motion Power Spectrum" &
	spectrum fbm.5 i f | xgraph -x "log10 freq" -y "log10 Power" -t "H=0.5 Fractal Motion Power Spectrum" &
	spectrum fbm.9 i f | xgraph -x "log10 freq" -y "log10 Power" -t "H=0.9 Fractal Motion Power Spectrum" &


grstn:
	station fbm.1 i f |xgraph -x "log10 (t2-t1)" -y "log10 (Var [X(t2) - X(t1)])" -t "H=0.1 Fractal Motion Stationarity" &
	station fbm.3 i f |xgraph -x "log10 (t2-t1)" -y "log10 (Var [X(t2) - X(t1)])" -t "H=0.3 Fractal Motion Stationarity" &
	station fbm.5 i f |xgraph -x "log10 (t2-t1)" -y "log10 (Var [X(t2) - X(t1)])" -t "H=0.5 Fractal Motion Stationarity" &
	station fbm.7 i f |xgraph -x "log10 (t2-t1)" -y "log10 (Var [X(t2) - X(t1)])" -t "H=0.7 Fractal Motion Stationarity" &
	station fbm.9 i f |xgraph -x "log10 (t2-t1)" -y "log10 (Var [X(t2) - X(t1)])" -t "H=0.9 Fractal Motion Stationarity" &

