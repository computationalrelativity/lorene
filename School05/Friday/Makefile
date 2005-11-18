SHELL=/bin/sh

include $(HOME_LORENE)/local_settings		  # defines which compiler,...

LIB_G = -L$(HOME_LORENE)/Lib -llorene_g -llorenef77_g 
LIB = -L$(HOME_LORENE)/Lib -llorene -llorenef77 

.SUFFIXES : .o .C

EXE1 = $(SRC1:.C=)
SRC1 = wave_expl.C
OBJ1 = $(SRC1:.C=.o)

EXE2 = $(SRC2:.C=)
SRC2 = wave_impl.C
OBJ2 = $(SRC2:.C=.o)

EXE3 = $(SRC3:.C=)
SRC3 = wave_impl_bc3.C
OBJ3 = $(SRC3:.C=.o)


all: $(EXE1) $(EXE2) $(EXE3)

$(EXE1): $(OBJ1)
	$(CXX) -o $@ $(CXXFLAGS_G) $(OBJ1) $(LIB_G) $(LIB_LAPACK) \
					$(LIB_GSL) $(LIB_PGPLOT) $(LIB_CXX)

$(EXE2): $(OBJ2)
	$(CXX) -o $@ $(CXXFLAGS_G) $(OBJ2) $(LIB_G) $(LIB_LAPACK) \
					$(LIB_GSL) $(LIB_PGPLOT) $(LIB_CXX)

$(EXE3): $(OBJ3)
	$(CXX) -o $@ $(CXXFLAGS_G) $(OBJ3) $(LIB_G) $(LIB_LAPACK) \
	$(LIB_PGPLOT) $(LIB_GSL) $(LIB_CXX) 


.C.o:
	$(CXX)  -c $(CXXFLAGS_G) $(INC) $<

clean:
	rm -f $(OBJ1) $(OBJ2) $(OBJ3)
	rm -f $(EXE1) $(EXE2) $(EXE3)


