
all: integ frobenius

integ:  integ.o
frobenius:	frobenius.o

.o:
	cc -o $@ $< -lm

.c.o:
	cc -g -c -O3 $^


clean:
	rm -f *.o frobenius 
