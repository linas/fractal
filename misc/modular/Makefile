

CC = /usr/bin/cc -O2  

EXES= orbit


all: 	$(EXES)

clean:
	rm -f core tmp junk glop a.out *.o

realclean:	clean
	rm -f $(EXES) 

.c.o:
	$(CC) -c -O2 $<

.o:
	$(CC) -o $* $^ -lglut -lGLU -lGL -lm -lstdc++


orbit.o: orbit.c
three-rep.o: three-rep.c

orbit: orbit.o three-rep.o

