
REALTOP = /clone/tools

all: 
	cd $(REALTOP)/src/libs/file; make TOP=$(REALTOP)
	cd $(REALTOP)/src/libs/stoch; make TOP=$(REALTOP)
	cd $(REALTOP)/src/libs/spectra; make TOP=$(REALTOP)
	cd $(REALTOP)/src/cmds; make  TOP=$(REALTOP)

clean: 
	cd $(REALTOP)/src/libs/file; make clean TOP=$(REALTOP)
	cd $(REALTOP)/src/libs/stoch; make clean TOP=$(REALTOP)
	cd $(REALTOP)/src/libs/spectra; make clean TOP=$(REALTOP)
	cd $(REALTOP)/src/cmds; make clean  TOP=$(REALTOP)
	cd $(REALTOP)/lib; make clean  TOP=$(REALTOP)
	cd $(REALTOP)/bin; make clean  TOP=$(REALTOP)

realclean: 
	cd $(REALTOP)/src/libs/file; make realclean TOP=$(REALTOP)
	cd $(REALTOP)/src/libs/stoch; make realclean TOP=$(REALTOP)
	cd $(REALTOP)/src/libs/spectra; make realclean TOP=$(REALTOP)
	cd $(REALTOP)/src/cmds; make realclean  TOP=$(REALTOP)
	cd $(REALTOP)/lib; make realclean  TOP=$(REALTOP)
	cd $(REALTOP)/bin; make realclean  TOP=$(REALTOP)
