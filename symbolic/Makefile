
CC = /usr/bin/cc -O2  

EXES= tensor


all: 	$(EXES)

clean:
	rm -f core tmp junk glop a.out *.o

realclean:	clean
	rm -f $(EXES) 

.c.o:
	$(CC) -c -O2 $<

.o:
	$(CC) -o $* $^ -lm


tensor:	tensor.o

