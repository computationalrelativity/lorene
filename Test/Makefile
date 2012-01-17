SHELL=/bin/sh

include $(HOME_LORENE)/local_settings		  # defines which compiler,...

LIB	  = -L$(HOME_LORENE)/Lib -llorene_g -llorenef77_g \
		$(LIB_GSL) \
		$(LIB_LAPACK) \
		$(LIB_PGPLOT) \
		$(LIB_CXX) 

.SUFFIXES : .o .C

install: test_fft test_lapack # test_pgplot

test_fft: test_fft.o
	$(CXX) -o $@ $(CXXFLAGS_G) test_fft.o $(LIB)

test_lapack: test_lapack.o
	$(CXX) -o $@ $(CXXFLAGS_G) test_lapack.o $(LIB)

test_pgplot: test_pgplot.o
	$(CXX) -o $@ $(CXXFLAGS_G) test_pgplot.o $(LIB)


.C.o:
	$(CXX)  -c $(CXXFLAGS_G) $(INC) $<

