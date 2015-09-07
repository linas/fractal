
INCLUDES =  -I .

all: lytic

lytic: lytic.o ../sum/sum.o

.o:
	cc -o $@ $^ -lm -lstdc++


.C.o:
	cc -c $^  -Wall -O2 $(INCLUDES)
