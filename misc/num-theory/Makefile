############### Symbol Definitions ###################### 


EXES= ndist series

all: 	$(EXES)

clean:
	rm -f core tmp junk glop a.out *.o

realclean:	clean
	rm -f $(EXES) 

dist:
	tar -zcvf num.tgz *.c *.h *.gplot Makefile .cvsignore

.c.o:
	$(CC) -c -O2 $<

.o:
	$(CC) -o $* $^ -lm


ndist:	ndist.o gcf.o totient.o
series:	series.o divisor.o gcf.o totient.o

ndist.o:	ndist.c gcf.h totient.h 

divisor.o:	divisor.c divisor.h
gcf.o:	gcf.c gcf.h
totient.o:	totient.c totient.h
