SHELL=/bin/sh

include $(HOME_LORENE)/local_settings		  # defines which compiler,...

LIB_G = -L$(HOME_LORENE)/Lib -llorene_g -llorenef77_g
LIBLOR = $(HOME_LORENE)/Lib/liblorene_g.a $(HOME_LORENE)/Lib/liblorenef77_g.a 

.SUFFIXES : .o .C

SRC1 = GraVaStar.C
OBJ1 = $(SRC1:.C=.o)
EXE1 = $(SRC1:.C=)
SRC2 = read_star.C
OBJ2 = $(SRC2:.C=.o)
EXE2 = $(SRC2:.C=)

$(EXE1): $(OBJ1) $(LIBLOR) 
	$(CXX) -o $@ $(CXXFLAGS_G) $(OBJ1) $(LIB_G) $(LIB_LAPACK) \
				        $(LIB_GSL) $(LIB_PGPLOT) $(LIB_CXX)

$(EXE2): $(OBJ2) $(LIBLOR) 
	$(CXX) -o $@ $(CXXFLAGS_G) $(OBJ2) $(LIB_G) $(LIB_LAPACK) \
				        $(LIB_GSL) $(LIB_PGPLOT) $(LIB_CXX)
.C.o:
	$(CXX)  -c $(CXXFLAGS_G) $(INC) $<

all: $(EXE1) $(EXE2)

uninstall:
	rm -f $(OBJ1) 
	rm -f $(EXE1)
	rm -f $(OBJ2) 
	rm -f $(EXE2)
