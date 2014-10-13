/*
 *  Member functions of the Scalar class for the computation of integrals.
 */

/*
 *   Copyright (c) 2003 Eric Gourgoulhon & Jerome Novak
 *   
 *   Copyright (c) 1999-2001 Eric Gourgoulhon (Cmp version)
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


char scalar_integ_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.4  2014/10/13 08:53:46  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.3  2004/07/26 16:02:23  j_novak
 * Added a flag to specify whether the primitive should be zero either at r=0
 * or at r going to infinity.
 *
 * Revision 1.2  2004/06/14 15:28:17  e_gourgoulhon
 * Added method primr().
 *
 * Revision 1.1  2003/09/25 09:33:36  j_novak
 * Added methods for integral calculation and various manipulations
 *
 *
 * $Header$
 *
 */

// Headers Lorene
#include "tensor.h"
#include "cmp.h"

		    //-----------------------------------//
		    //	   Integral over all space	 //
		    //-----------------------------------//

namespace Lorene {
double Scalar::integrale() const {
    
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

const Tbl& Scalar::integrale_domains() const {
    
    // Protection
    assert(etat != ETATNONDEF) ;

    // If the integrals have not been previously computed, the 
    //  computation must be done by the appropriate routine of the mapping : 

    if (p_integ == 0x0) {
      Cmp orig(*this) ;
      p_integ = mp->integrale(orig) ;
    }
    
    return *p_integ ;

}


            //----------------------//
            //  Radial primitive    // 
            //----------------------//

Scalar Scalar::primr(bool null_infty) const {

    Scalar resu(*mp) ; 

    mp->primr(*this, resu, null_infty) ; 
    
    return resu ; 
}




}
