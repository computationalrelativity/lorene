SHELL=/bin/sh

include $(HOME_LORENE)/local_settings		  # defines which compiler,...

LIB_G = -L$(HOME_LORENE)/Lib -llorene_g -llorenef77_g 

.SUFFIXES : .o .C .f

EXE = test_tensor
SRCC = test_tensor.C 
SRCF = 
OBJ = $(SRCC:.C=.o) $(SRCF:.f=.o)

$(EXE): $(OBJ)
	$(CXX) -o $@ $(CXXFLAGS_G) $(OBJ) $(LIB_G) $(LIB_LAPACK) \
					$(LIB_PGPLOT) $(LIB_CXX)

.C.o:
	$(CXX)  -c $(CXXFLAGS_G) $(INC) $<

.f.o:
	$(F77)  -c $(F77FLAGS_G) $<

clean:
	rm -f $(OBJ)
	rm -f $(EXE)


