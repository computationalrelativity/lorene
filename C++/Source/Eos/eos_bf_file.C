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
 * Revision 1.10  2015/06/10 14:39:17  a_sourie
 * New class Eos_bf_tabul for tabulated 2-fluid EoSs and associated functions for the computation of rotating stars with such EoSs.
 *
 * Revision 1.9  2014/10/13 08:52:52  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.8  2014/10/06 15:13:06  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.7  2014/04/25 10:43:51  j_novak
 * The member 'name' is of type string now. Correction of a few const-related issues.
 *
 * Revision 1.6  2008/08/19 06:42:00  j_novak
 * Minor modifications to avoid warnings with gcc 4.3. Most of them concern
 * cast-type operations, and constant strings that must be defined as const char*
 *
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
#include <cstdlib>

// Header Lorene
#include "headcpp.h"
#include "eos_bifluid.h"
#include "utilitaires.h"

		//--------------------------------------//
		//  Identification virtual functions	//
		//--------------------------------------//


namespace Lorene {
int Eos_bf_poly::identify() const	{ return 1; }

int Eos_bf_poly_newt::identify() const	{ return 2; }

int Eos_bf_tabul::identify() const	{ return 3; }

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
	
        case 3 : {
	    p_eos = new Eos_bf_tabul(fich) ; 
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

Eos_bifluid* Eos_bifluid::eos_from_file(const char *fname) {
    
    int identificator ; 

    // EOS identificator : 
    if (read_variable (fname, const_cast<char*>("ident"), identificator) != 0)
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

		//----------------------------------------------//
		//    EOS construction from a formatted file    //
		//----------------------------------------------//

Eos_bifluid* Eos_bifluid::eos_from_file(ifstream& fich) {
    
    int identificator ; 
    char blabla[80] ;

    // EOS identificator : 
    fich >> identificator ; fich.getline(blabla, 80) ;

    Eos_bifluid* p_eos ; 
    
    switch(identificator) {
	
	case 3 : {
	    p_eos = new Eos_bf_tabul(fich) ; 
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






}
