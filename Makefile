

all:

clean:
	cd generate; make clean
	cd misc; make clean
	cd queue; make clean
	cd symbolic; make clean
	cd tools; make clean

realclean:
	cd generate; make realclean
	cd misc; make realclean
	cd queue; make realclean
	cd symbolic; make realclean
	cd tools; make realclean

