SHELL=/bin/sh

include $(HOME_LORENE)/local_settings		  # defines which compiler,...

LIB_G = -L$(HOME_LORENE)/Lib -llorene_g -llorenef77_g 
LIBLOR = $(HOME_LORENE)/Lib/liblorene_g.a $(HOME_LORENE)/Lib/liblorenef77_g.a 

.SUFFIXES : .o .C

SRC1 = wave_evol.C
OBJ1 = $(SRC1:.C=.o)
EXE1 = $(SRC1:.C=)

SRC2 = einstein.C
OBJ2 = $(SRC2:.C=.o)
EXE2 = $(SRC2:.C=)

SRC3 = read_tslice.C
OBJ3 = $(SRC3:.C=.o)
EXE3 = $(SRC3:.C=)

SRC4 = visu_evol.C
OBJ4 = $(SRC4:.C=.o)
EXE4 = $(SRC4:.C=)

SRC5 = analyse_evol.C
OBJ5 = $(SRC5:.C=.o)
EXE5 = $(SRC5:.C=)

$(EXE1): $(OBJ1) $(LIBLOR) 
	$(CXX) -o $@ $(CXXFLAGS_G) $(OBJ1) $(LIB_G) $(LIB_LAPACK) \
	$(LIB_PGPLOT) $(LIB_GSL) $(LIB_CXX) 

$(EXE2): $(OBJ2) $(LIBLOR) 
	$(CXX) -o $@ $(CXXFLAGS_G) $(OBJ2) $(LIB_G) $(LIB_LAPACK) \
	$(LIB_PGPLOT) $(LIB_GSL) $(LIB_CXX)

$(EXE3): $(OBJ3) $(LIBLOR) 
	$(CXX) -o $@ $(CXXFLAGS_G) $(OBJ3) $(LIB_G) $(LIB_LAPACK) \
	$(LIB_PGPLOT) $(LIB_GSL) $(LIB_CXX)

$(EXE4): $(OBJ4) $(LIBLOR) 
	$(CXX) -o $@ $(CXXFLAGS_G) $(OBJ4) $(LIB_G) $(LIB_LAPACK) \
	$(LIB_PGPLOT) $(LIB_GSL) $(LIB_CXX)

$(EXE5): $(OBJ5) $(LIBLOR) 
	$(CXX) -o $@ $(CXXFLAGS_G) $(OBJ5) $(LIB_G) $(LIB_LAPACK) \
	$(LIB_PGPLOT) $(LIB_GSL) $(LIB_CXX)

all: $(EXE1) $(EXE2) $(EXE3) $(EXE4) $(EXE5)


.C.o:
	$(CXX)  -c $(CXXFLAGS_G) $(INC) $<

uninstall:
	rm -f $(OBJ1) $(OBJ2) $(OBJ3) $(OBJ4) $(OBJ5)
	rm -f $(EXE1) $(EXE2) $(EXE3) $(EXE4) $(EXE5)

