/*
 *  Member functions of the Cmp class for the computation of integrals.
 */

/*
 *   Copyright (c) 1999-2001 Eric Gourgoulhon
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


char cmp_integ_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2001/11/20 15:19:27  e_gourgoulhon
 * Initial revision
 *
 * Revision 1.1  1999/12/09  10:50:21  eric
 * Initial revision
 *
 *
 * $Header$
 *
 */

// Headers Lorene
#include "map.h"
#include "cmp.h"

		    //-----------------------------------//
		    //	   Integral over all space	 //
		    //-----------------------------------//

double Cmp::integrale() const {
    
    const Tbl& integ = integrale_domains() ; 
    
    int nz = mp->get_mg()->get_nzone() ; 
    
    double resu = integ(0) ; 
    for (int l=1; l<nz; l++) {
	resu += integ(l) ; 
    }
    
    return resu ; 
}

		    //-----------------------------------//
		    //	   Integrals in each domain	 //
		    //-----------------------------------//

const Tbl& Cmp::integrale_domains() const {
    
    // Protection
    assert(etat != ETATNONDEF) ;

    // If the integrals have not been previously computed, the 
    //  computation must be done by the appropriate routine of the mapping : 

    if (p_integ == 0x0) {
        p_integ = mp->integrale(*this) ;
    }
    
    return *p_integ ;

}

