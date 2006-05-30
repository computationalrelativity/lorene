SHELL=/bin/sh

include $(HOME_LORENE)/local_settings		  # defines which compiler,...

LIB_G = -L$(HOME_LORENE)/Lib -llorene_g -llorenef77_g 
LIB = -L$(HOME_LORENE)/Lib -llorene -llorenef77 

.SUFFIXES : .o .C

EXE1 = $(SRC1:.C=)
SRC1 = simple_fields.C
OBJ1 = $(SRC1:.C=.o)

EXE2 = $(SRC2:.C=)
SRC2 = kerr_test.C
OBJ2 = $(SRC2:.C=.o)

EXE3 = $(SRC3:.C=)
SRC3 = demo_tensor.C
OBJ3 = $(SRC3:.C=.o)

EXE4 = $(SRC4:.C=)
SRC4 = sol_fields.C
OBJ4 = $(SRC4:.C=.o)

EXE5 = $(SRC5:.C=)
SRC5 = sol_kerr.C
OBJ5 = $(SRC5:.C=.o)

exo: $(EXE1) $(EXE2) $(EXE3)
sol: $(EXE4) $(EXE5)
all: exo sol

$(EXE1): $(OBJ1)
	$(CXX) -o $@ $(CXXFLAGS_G) $(OBJ1) $(LIB_G) $(LIB_LAPACK) \
					$(LIB_GSL) $(LIB_PGPLOT) $(LIB_CXX)

$(EXE2): $(OBJ2)
	$(CXX) -o $@ $(CXXFLAGS_G) $(OBJ2) $(LIB_G) $(LIB_LAPACK) \
					$(LIB_GSL) $(LIB_PGPLOT) $(LIB_CXX)

$(EXE3): $(OBJ3)
	$(CXX) -o $@ $(CXXFLAGS_G) $(OBJ3) $(LIB_G) $(LIB_LAPACK) \
	$(LIB_PGPLOT) $(LIB_GSL) $(LIB_CXX) 
	cp $(HOME_LORENE)/Visu/OpenDX/visu_section.net .

$(EXE4): $(OBJ4)
	$(CXX) -o $@ $(CXXFLAGS_G) $(OBJ4) $(LIB_G) $(LIB_LAPACK) \
					$(LIB_GSL) $(LIB_PGPLOT) $(LIB_CXX)

$(EXE5): $(OBJ5)
	$(CXX) -o $@ $(CXXFLAGS_G) $(OBJ5) $(LIB_G) $(LIB_LAPACK) \
					$(LIB_GSL) $(LIB_PGPLOT) $(LIB_CXX)

.C.o:
	$(CXX)  -c $(CXXFLAGS_G) $(INC) $<

clean:
	rm -f $(OBJ1) $(OBJ2) $(OBJ3) $(OBJ4) $(OBJ5)
	rm -f $(EXE1) $(EXE2) $(EXE3) $(EXE4) $(EXE5)


