SHELL=/bin/sh

SRC = lorene.C $(HOME_LORENE)/C++/Include/*.h 

install: features.dvi refguide

features.dvi: features.tex
	latex features
	rm -f *.aux *.log

refguide: $(SRC)
	mv refguide html
	doxygen lorene.dox
	mv html refguide

uninstall:
	cd $(HOME_LORENE)/C++/Doc++; $(MAKE) -i uninstall
	rm -f *.aux *.log *.dvi
	rm -fr refguide
	rm -f refguide.tex



