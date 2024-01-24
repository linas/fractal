
CC = cc -Wall -g -O2 $(INCLUDES)

GENDIR = ../../../generate
BRAT = $(GENDIR)/brat.o $(GENDIR)/fileio.o $(GENDIR)/util.o

PERIODIC = ../periodic

all: recheck force orbit visitation borel reforce tsum unrolling \
	unwrap unstack alpha gamma

.C.o:
	$(CC) -c $< -std=gnu++11 -I $(GENDIR)

.c.o:
	$(CC) -c $< -I $(PERIODIC)

.o:
	$(CC) -o $@ $^ -lm -lstdc++

visitation.o: visitation.c visits.c
unrolling.o: unrolling.c unref.c unutil.c
unwrap.o: unwrap.c unref.c unutil.c
unstack.o: unstack.c unref.c unutil.c
gamma.o: gamma.c unutil.c

recheck: recheck.o
force: force.o
orbit: orbit.o
visitation: visitation.o
borel: borel.o $(BRAT)
reforce: reforce.o
tsum: tsum.o
unrolling: unrolling.o
unwrap: unwrap.o
unstack: unstack.o
alpha: alpha.o
gamma: gamma.o
