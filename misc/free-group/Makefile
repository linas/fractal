
CC = /usr/bin/cc -O2  

EXES= free


all: 	$(EXES)

clean:
	rm -f core tmp junk glop a.out *.o

realclean:	clean
	rm -f $(EXES) 

.c.o:
	$(CC) -c -O2 $<

.o:
	$(CC) -o $* $^ -lm


free.o: free.c
free: free.o 

