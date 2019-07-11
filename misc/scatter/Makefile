
INCLUDES = -I ../../generate

HIST=../../generate/brat.o ../../generate/util.o


all: scatter

scatter: scatter.o
scatter.o: ../../generate/brat.h

.o:
	cc -o $@ $^ -lm -lstdc++

.C.o:
	c++ -c $^  -std=gnu++11 -Wall -g -O2 $(INCLUDES)
