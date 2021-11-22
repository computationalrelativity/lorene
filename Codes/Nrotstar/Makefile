SHELL=/bin/sh

include $(HOME_LORENE)/local_settings		  # defines which compiler,...

LIB_G = -L$(HOME_LORENE)/Lib -llorene_g -llorenef77_g
LIB = -L$(HOME_LORENE)/Lib -llorene -llorenef77
LIBLOR = $(HOME_LORENE)/Lib/liblorene_g.a $(HOME_LORENE)/Lib/liblorenef77_g.a 

.SUFFIXES : .o .C

SRC1 = nrotstar.C
OBJ1 = $(SRC1:.C=.o)
EXE1 = $(SRC1:.C=)

SRC2 = nrotseq.C
OBJ2 = $(SRC2:.C=.o)
EXE2 = $(SRC2:.C=)

SRC3 = read_star.C
OBJ3 = $(SRC3:.C=.o)
EXE3 = $(SRC3:.C=)

$(EXE1): $(OBJ1) $(LIBLOR) 
	$(CXX) -o $@ $(CXXFLAGS_G) $(OBJ1) $(LIB_G) $(LIB_LAPACK) \
				        $(LIB_GSL) $(LIB_PGPLOT) $(LIB_CXX)
$(EXE2): $(OBJ2) $(LIBLOR)
	$(CXX) -o $@ $(CXXFLAGS_G) $(OBJ2) $(LIB_G) $(LIB_LAPACK) \
					$(LIB_GSL) $(LIB_PGPLOT) $(LIB_CXX)
$(EXE3): $(OBJ3) $(LIBLOR)
	$(CXX) -o $@ $(CXXFLAGS_G) $(OBJ3) $(LIB_G) $(LIB_LAPACK) \
					$(LIB_GSL) $(LIB_PGPLOT) $(LIB_CXX)
.C.o:
	$(CXX)  -c $(CXXFLAGS_G) $(INC) $<

all: $(EXE1) $(EXE2) $(EXE3)

uninstall clean:
	rm -f $(OBJ1) $(OBJ2) $(OBJ3)
	rm -f $(EXE1) $(EXE2) $(EXE3)
