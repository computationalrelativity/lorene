SHELL=/bin/sh

include $(HOME_LORENE)/local_settings		  # defines which compiler,...

LIB_G = -L$(HOME_LORENE)/Lib -llorene_g -llorenef77_g 

.SUFFIXES : .o .C

SRC1 = test_vdf_poisson.C
OBJ1 = $(SRC1:.C=.o)
EXE1 = $(SRC1:.C=)

SRC2 = test_poisson_vect_spher.C
OBJ2 = $(SRC2:.C=.o)
EXE2 = $(SRC2:.C=)

SRC3 = test_poisson_div_free.C
OBJ3 = $(SRC3:.C=.o)
EXE3 = $(SRC3:.C=)

$(EXE1): $(OBJ1) $(HOME_LORENE)/Lib/liblorene_g.a
	$(CXX) -o $@ $(CXXFLAGS_G) $(OBJ1) $(LIB_G) $(LIB_LAPACK) \
					$(LIB_PGPLOT) $(LIB_CXX)
$(EXE2): $(OBJ2) $(HOME_LORENE)/Lib/liblorene_g.a
	$(CXX) -o $@ $(CXXFLAGS_G) $(OBJ2) $(LIB_G) $(LIB_LAPACK) \
					$(LIB_PGPLOT) $(LIB_CXX)
$(EXE3): $(OBJ3) $(HOME_LORENE)/Lib/liblorene_g.a
	$(CXX) -o $@ $(CXXFLAGS_G) $(OBJ3) $(LIB_G) $(LIB_LAPACK) \
					$(LIB_PGPLOT) $(LIB_CXX)

all: $(EXE1) $(EXE2) $(EXE3)


.C.o:
	$(CXX)  -c $(CXXFLAGS_G) $(INC) $<

uninstall:
	rm -f $(OBJ1)  
	rm -f $(EXE1) 
	rm -f $(OBJ2)  
	rm -f $(EXE2) 

