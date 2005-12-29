

EXES= density trinary

all: $(EXES)

density:	density.o

trinary:	trinary.o


.C.o:
	gcc -c $^ -I ../../tools/inc

.o:
	g++ -o $* $^ ../../tools/lib/Farey.o -lm


.c.o:
	gcc -c $^ 

clean:
	rm -f *.o a.out

realclean: clean
	rm -rf $(EXES)
