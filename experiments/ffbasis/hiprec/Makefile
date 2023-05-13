
CC = cc -std=gnu11

EXES = topsin
MPLIB=libanant.a

all: $(EXES)

topsin: topsin.o

.c.o:
	$(CC) -c -g -O2 -Wall -I. $<

.o:
	$(CC) -g -o $* $^ $(MPLIB) -lgmp -ldb -lm

