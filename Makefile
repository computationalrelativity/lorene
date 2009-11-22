SHELL=/bin/sh

all: doc cpp fortran export

install: all

doc:
	cd Doc; $(MAKE) -i

cpp:
	cd C++; $(MAKE)

fortran:
	cd F77; $(MAKE)

export:
	cd Export/C++; $(MAKE)

uninstall:
	cd Doc; $(MAKE) -i uninstall
	cd C++; $(MAKE) -i uninstall
	cd F77; $(MAKE) -i uninstall
	cd Export/C++; $(MAKE) -i uninstall
	rm -fr bin
	rm -fr Lib
	rm -f Test/*.o
	rm -f Test/test_fft Test/test_lapack Test/test_pgplot

test: 
	cd Test; $(MAKE)


.PHONY: all install doc cpp fortran export uninstall test
