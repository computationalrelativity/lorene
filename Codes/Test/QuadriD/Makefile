SHELL=/bin/sh

include $(HOME_LORENE)/local_settings		  # defines which compiler,...

LIB_G = -L$(HOME_LORENE)/Lib -llorene_g -llorenef77_g 
LIB = -L$(HOME_LORENE)/Lib -llorene -llorenef77 

.SUFFIXES : .o .C

EXE = test_qmetrique
SRCC = test_qmetrique.C 
SRCF = 
OBJ = $(SRCC:.C=.o) $(SRCF:.f=.o)

$(EXE): $(OBJ)
	$(CXX) -o $@ $(CXXFLAGS_G) $(OBJ) $(LIB_G) $(LIB_LAPACK) \
					$(LIB_GSL) $(LIB_CXX)

.C.o:
	$(CXX)  -c $(CXXFLAGS_G) $(INC) $<
.f.o:
	$(F77)  -c $(F77FLAGS_G) $<

uninstall:
	rm -f $(OBJ)
	rm -f $(EXE)
	rm -f *.d

