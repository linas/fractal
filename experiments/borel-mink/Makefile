
TOP = ../../tools
FUNCDIR = $(TOP)/inc

FUNC= $(TOP)/lib/libfunc.a

GENDIR = ../../generate
BRAT = $(GENDIR)/brat.o $(GENDIR)/util.o

all: borel borel-dbg

borel: borel.o $(FUNC) $(BRAT)
borel-dbg: borel-dbg.o $(FUNC)

borel-dbg.o: borel.cc
	cc -c $^ -o $@ -std=gnu++11 -Wall -O2 -g -I $(FUNCDIR) -DDEBUG

.o:
	cc -o $@ $^ -lpthread -lm -lstdc++ 

.cc.o:
	cc -c $^ -std=gnu++11 -Wall -O2 -g -I $(FUNCDIR) -I $(GENDIR)
