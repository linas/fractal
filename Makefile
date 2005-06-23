

all:

clean:
	cd generate; make clean
	cd image; make clean
	cd queue; make clean
	cd tools; make clean

realclean:
	cd generate; make realclean
	cd image; make realclean
	cd queue; make realclean
	cd tools; make realclean

