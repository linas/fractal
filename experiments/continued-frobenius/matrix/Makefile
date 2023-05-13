
#
# Makefile that uses lapack for matrix operations
#

TOP = ../../../tools
INC = $(TOP)/inc
LIB = $(TOP)/lib


all: hyper eigen comb hyper

# OBJS=../ache.o

comb: comb.o $(OBJS)
eigen: eigen.o $(OBJS)

hyper: hyper.o 

.o:
	cc -o $@ $^ -L $(LIB) -lfunc -llapack -lf77blas -latlas -lg2c -lm

.c.o:
	cc -g -O2 -c $< -I.. -I$(INC)

.C.o:
	gcc -g -O2  -DANSI_C -c $< -I.. -I$(INC)


dist:
	tar -zcvf x *.c *.h *.gplot .cvsignore Makefile
