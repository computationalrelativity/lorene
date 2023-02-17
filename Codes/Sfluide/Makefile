SHELL=/bin/sh

include $(HOME_LORENE)/local_settings		  # defines which compiler,...

LIB_G = -L$(HOME_LORENE)/Lib -llorene_g -llorenef77_g -lgsl -lgslcblas

.SUFFIXES : .o .C

EXE1 = sfstar
SRC1 = sfstar.C 
OBJ1 = $(SRC1:.C=.o)

EXE2 = glitch_rise
SRC2 = glitch_rise.C 
OBJ2 = $(SRC2:.C=.o)

all : $(EXE1) $(EXE2) 

$(EXE1): $(OBJ1) 
	$(CXX) -o $@ $(CXXFLAGS_G) $(OBJ1) $(LIB_G) $(LIB_LAPACK) \
					$(LIB_GSL) $(LIB_PGPLOT) $(LIB_CXX)

$(EXE2): $(OBJ2) 
	$(CXX) -o $@ $(CXXFLAGS_G) $(OBJ2) $(LIB_G) $(LIB_LAPACK) \
					$(LIB_GSL) $(LIB_PGPLOT) $(LIB_CXX)

.C.o:
	$(CXX)  -c $(CXXFLAGS_G) $(INC) $<

uninstall:
	rm -f $(OBJ1)
	rm -f $(EXE1)
	rm -f $(OBJ2)
	rm -f $(EXE2)
	rm -f *.d


clean: 
	rm *.o
