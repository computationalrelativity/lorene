SHELL=/bin/sh

include $(HOME_LORENE)/local_settings		  # defines which compiler,...

LIB_G = -L$(HOME_LORENE)/Lib -llorene_g -llorenef77_g 

.SUFFIXES : .o .C

SRC = cheby.C leg.C

OBJ = $(SRC:.C=.o)

SRC1 = solver.C
OBJ1 = $(SRC1:.C=.o)
EXE1 = $(SRC1:.C=)

SRC2 = sol_tau.C
OBJ2 = $(SRC2:.C=.o)
EXE2 = $(SRC2:.C=)

SRC3 = sol_coloc.C
OBJ3 = $(SRC3:.C=.o)
EXE3 = $(SRC3:.C=)

SRC4 = sol_galerkin.C
OBJ4 = $(SRC4:.C=.o)
EXE4 = $(SRC4:.C=)

SRC5 = sol_multi_tau.C
OBJ5 = $(SRC5:.C=.o)
EXE5 = $(SRC5:.C=)

SRC6 = sol_multi_homogeneous.C
OBJ6 = $(SRC6:.C=.o)
EXE6 = $(SRC6:.C=)

SRC7 = sol_multi_variation.C
OBJ7 = $(SRC7:.C=.o)
EXE7 = $(SRC7:.C=)

$(EXE1): $(OBJ1) $(OBJ)
	$(CXX) -o $@ $(CXXFLAGS_G) $(OBJ1) $(OBJ) $(LIB_G) $(LIB_LAPACK) \
					$(LIB_GSL) $(LIB_PGPLOT) $(LIB_CXX)

$(EXE2): $(OBJ2) $(OBJ)
	$(CXX) -o $@ $(CXXFLAGS_G) $(OBJ2) $(OBJ) $(LIB_G) $(LIB_LAPACK) \
					$(LIB_GSL) $(LIB_PGPLOT) $(LIB_CXX)
					
$(EXE3): $(OBJ3) $(OBJ)
	$(CXX) -o $@ $(CXXFLAGS_G) $(OBJ3) $(OBJ) $(LIB_G) $(LIB_LAPACK) \
					$(LIB_GSL) $(LIB_PGPLOT) $(LIB_CXX)

$(EXE4): $(OBJ4) $(OBJ)
	$(CXX) -o $@ $(CXXFLAGS_G) $(OBJ4) $(OBJ) $(LIB_G) $(LIB_LAPACK) \
					$(LIB_GSL) $(LIB_PGPLOT) $(LIB_CXX)

$(EXE5): $(OBJ5) $(OBJ)
	$(CXX) -o $@ $(CXXFLAGS_G) $(OBJ5) $(OBJ) $(LIB_G) $(LIB_LAPACK) \
					$(LIB_GSL) $(LIB_PGPLOT) $(LIB_CXX)

$(EXE6): $(OBJ6) $(OBJ)
	$(CXX) -o $@ $(CXXFLAGS_G) $(OBJ6) $(OBJ) $(LIB_G) $(LIB_LAPACK) \
					$(LIB_GSL) $(LIB_PGPLOT) $(LIB_CXX)

$(EXE7): $(OBJ7) $(OBJ)
	$(CXX) -o $@ $(CXXFLAGS_G) $(OBJ7) $(OBJ) $(LIB_G) $(LIB_LAPACK) \
					$(LIB_GSL) $(LIB_PGPLOT) $(LIB_CXX)

.C.o:
	$(CXX)  -c $(CXXFLAGS_G) $(INC) $<
