
#
# Linas Vepstas dec 1994

TOP = /usr/home/src/tools

all: 
	cd $(TOP)/src/libs/file; make
	cd $(TOP)/src/libs/stoch; make
	cd $(TOP)/src/libs/spectra; make
	cd $(TOP)/src/cmds; make 


clean: 
	cd $(TOP)/src/libs/file; make clean
	cd $(TOP)/src/libs/stoch; make clean
	cd $(TOP)/src/libs/spectra; make clean
	cd $(TOP)/src/cmds; make clean 


realclean: 
	cd $(TOP)/src/libs/file; make realclean
	cd $(TOP)/src/libs/stoch; make realclean
	cd $(TOP)/src/libs/spectra; make realclean
	cd $(TOP)/src/cmds; make realclean 
	cd $(TOP)/lib; rm -f *
	cd $(TOP)/bin; rm -f *
	cd $(TOP)/inc; rm -f *


