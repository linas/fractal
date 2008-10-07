
# Make the submission to arxiv xxxxxxxxxxxmath.DS/xxxxx
#
# 
#
GRPHS= dyadic.eps \
       farey.eps \
       quest.ps \
       three-adic.ps \
       akk.ps

all: mink-exact.lyx $(GRPHS)
	lyx mink-exact.lyx -e latex
	/usr/bin/latex mink-exact.tex
	/usr/bin/bibtex mink-exact.aux

dist: mink-exact.tex mink-exact.aux mink-exact.bbl $(GRPHS) 
	tar -zcvf mink-exact.tar $^
	ls -la mink-exact.tar.gz

clean:
	rm -f *.log *.blg *.dvi

realclean: clean
   rm -f *.log *.tex *.pdf *.dvi *.lyx~

.SUFFIXES: .ps .png

.png.ps:
	cat $^ |  pngtopnm | pnmscale 0.65 -nomix | pnmtopng -compression 9 > tmp.png
	convert tmp.png $*.ps
