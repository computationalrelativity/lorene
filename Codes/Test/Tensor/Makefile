SHELL=/bin/sh

include $(HOME_LORENE)/local_settings		  # defines which compiler,...

LIB_G = -L$(HOME_LORENE)/Lib -llorene_g -llorenef77_g 
LIBLOR = $(HOME_LORENE)/Lib/liblorene_g.a $(HOME_LORENE)/Lib/liblorenef77_g.a 

.SUFFIXES : .o .C

SRC0 = test_Amu.C
OBJ0 = $(SRC0:.C=.o)
EXE0 = $(SRC0:.C=)

SRC1 = test_connect.C
OBJ1 = $(SRC1:.C=.o)
EXE1 = $(SRC1:.C=)

SRC2 = test_tensor.C
OBJ2 = $(SRC2:.C=.o)
EXE2 = $(SRC2:.C=)

SRC3 = test_poisson_angu.C
OBJ3 = $(SRC3:.C=.o)
EXE3 = $(SRC3:.C=)

SRC4 = test_vector_df.C
OBJ4 = $(SRC4:.C=.o)
EXE4 = $(SRC4:.C=)

SRC5 = test_sym_tensor_tt.C
OBJ5 = $(SRC5:.C=.o)
EXE5 = $(SRC5:.C=)

SRC6 = test_tens_poisson.C
OBJ6 = $(SRC6:.C=.o)
EXE6 = $(SRC6:.C=)

SRC7 = test_sym_tensor.C
OBJ7 = $(SRC7:.C=.o)
EXE7 = $(SRC7:.C=)

SRC8 = test_sym_tensor_decomp.C
OBJ8 = $(SRC8:.C=.o)
EXE8 = $(SRC8:.C=)

SRC9 = test_scalar_mtbl.C
OBJ9 = $(SRC9:.C=.o)
EXE9 = $(SRC9:.C=)

$(EXE0): $(OBJ0) $(LIBLOR) 
	$(CXX) -o $@ $(CXXFLAGS_G) $(OBJ0) $(LIB_G) $(LIB_LAPACK) \
					$(LIB_GSL) $(LIB_PGPLOT) $(LIB_CXX)

$(EXE1): $(OBJ1) $(LIBLOR) 
	$(CXX) -o $@ $(CXXFLAGS_G) $(OBJ1) $(LIB_G) $(LIB_LAPACK) \
					$(LIB_GSL) $(LIB_PGPLOT) $(LIB_CXX)
$(EXE2): $(OBJ2) $(LIBLOR) 
	$(CXX) -o $@ $(CXXFLAGS_G) $(OBJ2) $(LIB_G) $(LIB_LAPACK) \
					$(LIB_GSL) $(LIB_PGPLOT) $(LIB_CXX)
$(EXE3): $(OBJ3) $(LIBLOR) 
	$(CXX) -o $@ $(CXXFLAGS_G) $(OBJ3) $(LIB_G) $(LIB_LAPACK) \
					$(LIB_GSL) $(LIB_PGPLOT) $(LIB_CXX)
$(EXE4): $(OBJ4) $(LIBLOR) 
	$(CXX) -o $@ $(CXXFLAGS_G) $(OBJ4) $(LIB_G) $(LIB_LAPACK) \
					$(LIB_GSL) $(LIB_PGPLOT) $(LIB_CXX)

$(EXE5): $(OBJ5) $(LIBLOR) 
	$(CXX) -o $@ $(CXXFLAGS_G) $(OBJ5) $(LIB_G) $(LIB_LAPACK) \
					$(LIB_GSL) $(LIB_PGPLOT) $(LIB_CXX)

$(EXE6): $(OBJ6) $(LIBLOR) 
	$(CXX) -o $@ $(CXXFLAGS_G) $(OBJ6) $(LIB_G) $(LIB_LAPACK) \
					$(LIB_GSL) $(LIB_PGPLOT) $(LIB_CXX)

$(EXE7): $(OBJ7) $(LIBLOR) 
	$(CXX) -o $@ $(CXXFLAGS_G) $(OBJ7) $(LIB_G) $(LIB_LAPACK) \
					$(LIB_GSL) $(LIB_PGPLOT) $(LIB_CXX)

$(EXE8): $(OBJ8) $(LIBLOR) 
	$(CXX) -o $@ $(CXXFLAGS_G) $(OBJ8) $(LIB_G) $(LIB_LAPACK) \
					$(LIB_GSL) $(LIB_PGPLOT) $(LIB_CXX)

$(EXE9): $(OBJ9) $(LIBLOR) 
	$(CXX) -o $@ $(CXXFLAGS_G) $(OBJ9) $(LIB_G) $(LIB_LAPACK) \
					$(LIB_GSL) $(LIB_PGPLOT) $(LIB_CXX)

all: $(EXE1) $(EXE2) $(EXE3) $(EXE4) $(EXE5) $(EXE6) $(EXE7) $(EXE8) $(EXE9)


.C.o:
	$(CXX)  -c $(CXXFLAGS_G) $(INC) $<

uninstall:
	rm -f $(OBJ0) $(OBJ1) $(OBJ2) $(OBJ3) $(OBJ4) $(OBJ5) $(OBJ6) $(OBJ7) $(OBJ8) $(OBJ9)
	rm -f $(EXE0) $(EXE1) $(EXE2) $(EXE3) $(EXE4) $(EXE5) $(EXE6) $(EXE7) $(EXE8) $(EXE9)

