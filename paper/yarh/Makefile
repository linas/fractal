

# These are originally png's, but must be converted to ps for arxiv!
IMGS= rzero-15.ps rzero-151.ps rzero-1051.ps \
		rzero-11051.ps rzero-151051.ps r-sum-1051.ps r-sum-11051.ps \
		r-sum-121051.ps zero-s12-151051.ps zero-s13-151051.ps \
		zero-s23-151051.ps shadow-11051.ps

# These are originally eps drawings.
GRAPHS= swap-s12.ps swap-s13.ps swap-s23.ps \
		swap-s14.ps zero_im_z12_close.ps \
		zero_im_z12_full.ps zero_im_z12_mid.ps zero_re_z12_close.ps \
		zero_re_z12_full.ps zero_re_z12_mid.ps power-re.ps power-im.ps \

all: arxiv

yarh.tex: yarh.lyx
	lyx yarh.lyx -e latex
	/usr/bin/latex yarh.tex
	/usr/bin/bibtex yarh.aux
	/usr/bin/latex yarh.tex

arxiv: yarh.tex $(IMGS) $(GRAPHS)
	tar -zcvf yarh.tar.gz yarh.tex yarh.bbl $(IMGS) $(GRAPHS)
	touch arxiv

yarh.pdf: yarh.lyx
	lyx yarh.lyx -e pdf

clean:
	rm -f .log *.pdf *.ps *.dvi *.lyx~ 

realclean: clean

.SUFFIXES: .eps .ps .png .jpg

.eps.ps:
	cp $^ $*.ps

.png.jpg:
	cat $^ |  pngtopnm | pnmtojpeg > $*.jpg

.png.ps:
	# cat $^ |  pngtopnm | pnmscale 0.65 -nomix | pnmtopng -compression 9 > tmp.png
	convert $^ $*.ps
