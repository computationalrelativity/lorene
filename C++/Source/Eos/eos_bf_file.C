/*
 * Methods for Eos_bifluid and file manipulation
 *
 * (see file eos_bifluid.h for documentation)
 */

/*
 *   Copyright (c) 2001 Jerome Novak
 *
 *   This file is part of LORENE.
 *
 *   LORENE is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   LORENE is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with LORENE; if not, write to the Free Software
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */


char eos_bf_file_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.5  2003/12/05 15:09:47  r_prix
 * adapted Eos_bifluid class and subclasses to use read_variable() for
 * (formatted) file-reading.
 *
 * Revision 1.4  2002/10/16 14:36:34  j_novak
 * Reorganization of #include instructions of standard C++, in order to
 * use experimental version 3 of gcc.
 *
 * Revision 1.3  2002/01/11 14:09:34  j_novak
 * Added newtonian version for 2-fluid stars
 *
 * Revision 1.2  2001/12/04 21:27:53  e_gourgoulhon
 *
 * All writing/reading to a binary file are now performed according to
 * the big endian convention, whatever the system is big endian or
 * small endian, thanks to the functions fwrite_be and fread_be
 *
 * Revision 1.1.1.1  2001/11/20 15:19:27  e_gourgoulhon
 * LORENE
 *
 * Revision 1.1  2001/06/21  15:22:15  novak
 * Initial revision
 *
 *
 * $Header$
 *
 */
 
// Headers C
#include <stdlib.h>

// Header Lorene
#include "headcpp.h"
#include "eos_bifluid.h"
#include "utilitaires.h"

		//--------------------------------------//
		//  Identification virtual functions	//
		//--------------------------------------//


int Eos_bf_poly::identify() const	{ return 1; }

int Eos_bf_poly_newt::identify() const	{ return 2; }


		//---------------------------------------------//
		//    EOS construction from a binary file      //
		//---------------------------------------------//

Eos_bifluid* Eos_bifluid::eos_from_file(FILE* fich) {
    
    Eos_bifluid* p_eos ; 
    
    // Type (class) of EOS :
    int identificator ;     
    fread_be(&identificator, sizeof(int), 1, fich) ;		

    switch(identificator) {
	
	case 1 : {
	    p_eos = new Eos_bf_poly(fich) ; 
	    break ; 
	}
	
	default : {
	    cout << "Eos_bifluid::eos_from_file : unknown type of EOS !" << endl ; 
	    cout << " identificator = " << identificator << endl ; 
	    abort() ; 
	    break ; 
	}
	
    }
    
    return p_eos ; 
    
}

		//----------------------------------------------//
		//    EOS construction from a formatted file    //
		//----------------------------------------------//

Eos_bifluid* Eos_bifluid::eos_from_file(char *fname) {
    
    int identificator ; 

    // EOS identificator : 
    if (read_variable (fname, "ident", identificator) != 0)
      {
	cerr << "ERROR: Could not read the required variable 'ident' in " << fname << endl;
	exit (-1);
      }

    Eos_bifluid* p_eos ; 
    
    switch(identificator) {
	
	case 1 : {
	    p_eos = new Eos_bf_poly(fname) ; 
	    break ; 
	}
	
	case 2 : {
	    p_eos = new Eos_bf_poly_newt(fname) ; 
	    break ; 
	}
	
	default : {
	    cout << "Eos_bifluid::eos_from_file : unknown type of EOS !" << endl ; 
	    cout << " identificator = " << identificator << endl ; 
	    abort() ; 
	    break ; 
	}
	
    }
    
    return p_eos ; 
    
}






