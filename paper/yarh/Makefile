all:

dist:
	cp swap-s12.eps swap-s12.ps
	cp swap-s23.eps swap-s23.ps
	cp swap-s13.eps swap-s13.ps
	cp swap-s14.eps swap-s14.ps
	cp zero_im_z12_close.eps zero_im_z12_close.ps
	cp zero_im_z12_full.eps zero_im_z12_full.ps
	cp zero_im_z12_mid.eps zero_im_z12_mid.ps
	cp zero_re_z12_close.eps zero_re_z12_close.ps
	cp zero_re_z12_full.eps zero_re_z12_full.ps
	cp zero_re_z12_mid.eps zero_re_z12_mid.ps
	cp power-re.eps power-re.ps
	cp power-im.eps power-im.ps
	tar -zcvf yarh.tar.gz yarh.tex swap-s12.ps swap-s13.ps swap-s23.ps \
		swap-s14.ps rzero-15.png rzero-151.png rzero-1051.png \
		rzero-11051.png rzero-151051.png r-sum-1051.png r-sum-11051.png \
		r-sum-121051.png zero-s12-151051.png zero-s13-151051.png \
		zero-s23-151051.png shadow-11051.png zero_im_z12_close.ps \
		zero_im_z12_full.ps zero_im_z12_mid.ps zero_re_z12_close.ps \
		zero_re_z12_full.ps zero_re_z12_mid.ps power-re.ps power-im.ps

clean:
	rm -f .log *.pdf *.dvi *.lyx~ 

realclean: clean
