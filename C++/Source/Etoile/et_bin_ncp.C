/*
 *  Methods for class Et_bin_ncp
 *
 */

/*
 *   Copyright (c) 2002  Francois Limousin
 *
 *   This file is part of LORENE.
 *
 *   LORENE is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License version 2
 *   as published by the Free Software Foundation.
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


char et_bin_ncp_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2002/12/09 10:46:50  f_limousin
 * Methods for class Et_bin_ncp.
 *
 *
 *
 *
 * $Header$
 *
 */


// Lorene headers
#include "et_bin_ncp.h"


			    //--------------//
			    // Constructors //
			    //--------------//

// Standard constructor
// --------------------

Et_bin_ncp::Et_bin_ncp(Map& mp_i, int nzet_i, bool relat, const Eos& eos_i,
		       bool irrot, const Base_vect& ref_triad_i, const Metrique& flat0) 
             : Etoile_bin(mp_i, nzet_i, relat, eos_i, irrot, ref_triad_i),
               gamma(mp_i),
	       flat(flat0),
               gamma_tilde(mp_i, gamma, flat0) {

}

// Copy constructor
// ----------------

Et_bin_ncp::Et_bin_ncp(const Et_bin_ncp& et)
      	   : Etoile_bin(et),
	     gamma(et.gamma),
	     flat(et.flat),
             gamma_tilde(et.gamma_tilde){}


			    //------------//
			    // Destructor //
			    //------------//

Et_bin_ncp::~Et_bin_ncp(){}



			    //--------------//
			    //  Assignment  //
			    //--------------//

// Assignment to another Et_bin_ncp
// --------------------------------
void Et_bin_ncp::operator=(const Et_bin_ncp& et) {

// Assignement of proper quantities of class Etoile_bin
gamma = et.gamma ; 
flat = et.flat ;
gamma_tilde = et.gamma_tilde ;

}



			    //--------------//
			    //	  Outputs   //
			    //--------------//

// Save in a file
// --------------
void Et_bin_ncp::sauve(FILE* fich) const {
    
    Etoile_bin::sauve(fich) ; 
    
    gamma.sauve(fich) ; 
    flat.sauve(fich) ; 
    gamma_tilde.sauve(fich) ; 
        
}

// Printing
// --------

ostream& Et_bin_ncp::operator>>(ostream& ost) const {
    
    Etoile_bin::operator>>(ost) ; 
    
    ost << endl ; 
    ost << "Star in a binary system with non conformally flat metric" 
	<< endl ; 
    ost << "--------------------------------------------------------" 
	<< endl ; 

    ost << "Gamma : " << gamma << endl ; 
    ost << "Flat : " << flat << endl ; 
    ost << "Gamma_tilde : " << gamma_tilde << endl ; 

    return ost ; 

}
