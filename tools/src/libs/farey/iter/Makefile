
SRCS= fiter.c

OBJS= $(SRCS:.c=.o)

EXES= $(SRCS:.c=)

LIBS= ../../../../lib/Farey.o

all: $(EXES)

# $(EXES):	$(OBJS)

.o:
	cc -o $* $< $(LIBS) -lm


$(OBJS):	$(SRCS)
	cc -c -I.. $<


clean:
	rm -f *.o junk core tmp a.out

realclean: clean
	rm -f $(EXES)
