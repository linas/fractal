
CC = cc -Wall -g -O2 $(INCLUDES)

all: eigen

.c.o:
   $(CC) -c $< -I $(FUNCDIR)

.o:
   $(CC) -o $@ $^ -lm -lstdc++

eigen: eigen.o -lgsl -lgslcblas
