
all: tarball

tarball:
	tar -zcvf paper.tgz Makefile *.lyx */*.lyx */*.png */*.jpg */*.gif */README */.cvsignore


clean:
	rm -f paper.tgz *.tex *.lyx~
	cd cantor; make clean
	cd de_Rham; make clean
	cd dedekind; make clean
	cd diff; make clean
	cd fdist; make clean
	cd frac-zeta; make clean
	cd gap; make clean
	cd gkw; make clean
	cd measure; make clean
	cd minkowski; make clean
	cd poch-zeta; make clean
	cd rat; make clean
	cd rotations; make clean
	cd scrapbook; make clean
	cd series; make clean
	cd takagi; make clean
	cd yarh; make clean

realclean: clean
