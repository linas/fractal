
CC = /usr/bin/cc -g -O2  

EXES= tensor


all: 	$(EXES)

clean:
	rm -f core tmp junk glop a.out *.o

realclean:	clean
	rm -f $(EXES) 

.c.o:
	$(CC) -c -O2 $<  `pkg-config --cflags glib`

.o:
	$(CC) -o $* $^ -lm `pkg-config --libs glib`


tensor:	tensor.o

