#
# FILE:
# Makefile


LIBSRCS=			\
	pairhash.C		\
	wordhash.C

LIBOBJS=$(LIBSRCS:.C=.o)

.C.o:
	cc -c -g -O2 $< 

analyze: analyze.o $(LIBOBJS)
	cc analyze.o -o analyze $(LIBOBJS)

clean:
	rm -f junk glop tmp a.out core
	
realclean: clean
	rm -f $(LIBOBJS)
	makedepend
	rm -f Makefile.bak

depend:
	makedepend -- $(CLAGS) $(DEFINES) $(INCLUDES) \
	-I /usr/lib/g++-include			\
	 -I /usr/lib/gcc-lib/i486-linux/2.5.8/include -- \
	$(LIBSRCS) $(EXESRCS)

# ------------------------------------------------------------
# DO NOT DELETE THIS LINE -- makedepend depends on it.
