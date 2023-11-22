
CC = cc -std=gnu++11 -Wall -g -O2 $(INCLUDES)

GENDIR = ../../generate
BRAT = $(GENDIR)/brat.o $(GENDIR)/util.o

all: polylog

.C.o:
	$(CC) -c $< -I $(GENDIR)

.o:
	$(CC) -o $@ $^ -lpthread -lm -lstdc++

polylog.o: polylog.C

polylog: polylog.o $(BRAT)
