SHELL=/bin/sh

include $(HOME_LORENE)/local_settings		  # defines which compiler,...

LIB_G = -L$(HOME_LORENE)/Lib -llorene_g -llorenef77_g -lgsl -lgslcblas

.SUFFIXES : .o .C

EXE = sfstar
SRC = sfstar.C
OBJ = $(SRC:.C=.o)

CXXFLAGS_G = -g -Wall

$(EXE): $(OBJ) 
	$(CXX) -o $@ $(CXXFLAGS_G) $(OBJ) $(LIB_G) $(LIB_LAPACK) \
					$(LIB_GSL) $(LIB_PGPLOT) $(LIB_CXX)

.C.o:
	$(CXX)  -c $(CXXFLAGS_G) $(INC) $<

uninstall:
	rm -f $(OBJ)
	rm -f $(EXE)
	rm -f *.d


clean: 
	rm *.o
