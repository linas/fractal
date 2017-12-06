
CC = cc -std=gnu++11 -Wall -g -O2 $(INCLUDES)

GENDIR = ../../generate
BRAT = $(GENDIR)/brat.o $(GENDIR)/util.o

FUNCDIR = ../../tools/inc
FUNC=../../tools/lib/libfunc.a

all: bernie

.C.o:
	$(CC) -c $< -I $(GENDIR) -I $(FUNCDIR)

.o:
	$(CC) -o $@ $^ $(FUNC) -lpthread -lm -lstdc++

bernie.o: bernie.C

bernie: bernie.o $(BRAT)
