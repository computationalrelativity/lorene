SHELL=/bin/sh

install: features.dvi refguide

features.dvi:
	latex features
	rm -f *.aux *.log

refguide:
	cd $(HOME_LORENE)/C++/Doc++; $(MAKE)

refguide.dvi:
	cd $(HOME_LORENE)/C++/Doc++; $(MAKE) refguide.tex
	latex refguide
	latex refguide
	rm -f *.aux *.log
	xdvi refguide &

uninstall:
	cd $(HOME_LORENE)/C++/Doc++; $(MAKE) -i uninstall
	rm -f *.aux *.log *.dvi
	rm -fr refguide
	rm -f refguide.tex



