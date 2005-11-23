SHELL=/bin/sh

include $(HOME_LORENE)/local_settings		  # defines which compiler,...

LIB_G = -L$(HOME_LORENE)/Lib -llorene_g -llorenef77_g 

.SUFFIXES : .o .C

SRC = sphere.C monopole.C monopole_functions.C monopole_config.C
OBJ = $(SRC:.C=.o)
EXE = $(SRC:.C=)

$(EXE): $(OBJ)
	$(CXX) -o $@ $(CXXFLAGS_G) $(OBJ) $(LIB_G) $(LIB_LAPACK) \
					$(LIB_GSL) $(LIB_PGPLOT) $(LIB_CXX)


.C.o:
	$(CXX)  -c $(CXXFLAGS_G) $(INC) $<
