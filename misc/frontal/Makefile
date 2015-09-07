
INCLUDES =  -I .

all: lytic-1d lytic-parts

lytic-1d: lytic-1d.o lytic.o ../sum/sum.o
lytic-parts: lytic-parts.o lytic.o ../sum/sum.o

.o:
	cc -o $@ $^ -lm -lstdc++


.C.o:
	cc -c $^  -Wall -O2 $(INCLUDES)
