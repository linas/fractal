
TOP = ../../tools
FUNCDIR = $(TOP)/inc

FUNC= $(TOP)/lib/libfunc.a

GENDIR = ../../generate
BRAT = $(GENDIR)/brat.o $(GENDIR)/util.o

all: borel

borel: borel.o $(FUNC)

.o:
	cc -o $@ $^ -lm -lstdc++ 

.cc.o:
	cc -c $^ -std=gnu++11 -Wall -O2 -I $(FUNCDIR) -I $(GENDIR)
