SHELL=/bin/sh

include $(HOME_LORENE)/local_settings		  # defines which compiler,...

LIB = $(HOME_LORENE)/Lib
SRC = *.f 

install: $(LIB)/liblorenef77_g.a $(LIB)/liblorenef77.a

$(LIB)/liblorenef77_g.a : $(SRC)
	$(F77) -c $(F77FLAGS_G) $?
	
$(LIB)/liblorenef77.a : $(SRC)
	$(F77) -c $(F77FLAGS) $?
	
