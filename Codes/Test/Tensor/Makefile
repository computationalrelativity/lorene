SHELL=/bin/sh

include $(HOME_LORENE)/local_settings		  # defines which compiler,...

LIB_G = -L$(HOME_LORENE)/Lib -llorene_g -llorenef77_g 

.SUFFIXES : .o .C

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


$(EXE1): $(OBJ1)
	$(CXX) -o $@ $(CXXFLAGS_G) $(OBJ1) $(LIB_G) $(LIB_LAPACK) \
					$(LIB_PGPLOT) $(LIB_CXX)

$(EXE2): $(OBJ2)
	$(CXX) -o $@ $(CXXFLAGS_G) $(OBJ2) $(LIB_G) $(LIB_LAPACK) \
					$(LIB_PGPLOT) $(LIB_CXX)

$(EXE3): $(OBJ3)
	$(CXX) -o $@ $(CXXFLAGS_G) $(OBJ3) $(LIB_G) $(LIB_LAPACK) \
					$(LIB_PGPLOT) $(LIB_CXX)
$(EXE4): $(OBJ4)
	$(CXX) -o $@ $(CXXFLAGS_G) $(OBJ4) $(LIB_G) $(LIB_LAPACK) \
					$(LIB_PGPLOT) $(LIB_CXX)


.C.o:
	$(CXX)  -c $(CXXFLAGS_G) $(INC) $<

uninstall:
	rm -f $(OBJ1) $(OBJ2) $(OBJ3) $(OBJ4) 
	rm -f $(EXE1) $(EXE2) $(EXE3) $(EXE4)

