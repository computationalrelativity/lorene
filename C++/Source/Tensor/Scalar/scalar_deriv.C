/*
 * Computations of partial derivatives d/dx, d/dy and d/dz of a Scalar.
 */

/*
 *   Copyright (c) 1999-2001 Eric Gourgoulhon
 *   Copyright (c) 1999-2001 Philippe Grandclement
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


char scalar_deriv_C[] = "$Header$" ;



/*
 * $Id$
 * $Log$
 * Revision 1.1  2003/09/25 08:13:52  j_novak
 * Added method for calculating derivatives
 *
 *
 * $Header$
 *
 */
 
// Headers C
#include <stdlib.h>

// Headers Lorene
#include "tensor.h"
#include "cmp.h"

			//---------------------//
			//	d/dr	       //
			//---------------------//

const Scalar& Scalar::dsdr() const {

    // Protection
    assert(etat != ETATNONDEF) ;

    // If the derivative has not been previously computed, the 
    //  computation must be done by the appropriate routine of the mapping : 

    if (p_dsdr == 0x0) {
      Cmp orig(*this) ;
      Cmp deriv(mp) ;
      mp->dsdr(orig, deriv) ;
      p_dsdr = new Scalar(deriv) ; 
    }
    
    return *p_dsdr ;

}

			//------------------------//
			//	1/r d/dtheta      //
			//------------------------//

const Scalar& Scalar::srdsdt() const {

    // Protection
    assert(etat != ETATNONDEF) ;

    // If the derivative has not been previously computed, the 
    //  computation must be done by the appropriate routine of the mapping : 

    if (p_srdsdt == 0x0) {
      Cmp orig(*this) ;
      Cmp deriv(mp) ;
      mp->srdsdt(orig, deriv) ;
      p_srdsdt = new Scalar(deriv) ;
    }
    
    return *p_srdsdt ;

}


			//----------------------------------//
			//	1/(r sin(theta) d/dphi	    //
			//----------------------------------//

const Scalar& Scalar::srstdsdp() const {

    // Protection
    assert(etat != ETATNONDEF) ;

    // If the derivative has not been previously computed, the 
    //  computation must be done by the appropriate routine of the mapping : 

    if (p_srstdsdp == 0x0) {
      Cmp orig(*this) ;
      Cmp deriv(mp) ;
      mp->srstdsdp(orig, deriv) ;
      p_srstdsdp = new Scalar(deriv) ; 
    }
    
    return *p_srstdsdp ;

}

			//---------------------//
			//	d/dx	       //
			//---------------------//

const Scalar& Scalar::dsdx() const {

    // Protection
    assert(etat != ETATNONDEF) ;

    // If the derivative has not been previously computed, the 
    //  computation must be done by the appropriate routine of the mapping : 

    if (p_dsdx == 0x0) {
      Cmp orig(*this) ;
      Cmp deriv_r(dsdr()) ;
      Cmp deriv_t(srdsdt()) ;
      Cmp deriv_p(srstdsdp()) ;
      Cmp deriv_x(mp) ;
      mp->comp_x_from_spherical(deriv_r, deriv_t, deriv_p, deriv_x) ;
      p_dsdx = new Scalar(deriv_x) ; 
    }
    
    return *p_dsdx ;

}

			//---------------------//
			//	d/dy	       //
			//---------------------//

const Scalar& Scalar::dsdy() const {

    // Protection
    assert(etat != ETATNONDEF) ;

    // If the derivative has not been previously computed, the 
    //  computation must be done by the appropriate routine of the mapping : 

    if (p_dsdy == 0x0) {
      Cmp orig(*this) ;
      Cmp deriv_r(dsdr()) ;
      Cmp deriv_t(srdsdt()) ;
      Cmp deriv_p(srstdsdp()) ;
      Cmp deriv_y(mp) ;
      mp->comp_y_from_spherical(deriv_r, deriv_t, deriv_p, deriv_y) ;
      p_dsdy = new Scalar(deriv_y) ; 
    }
    
    return *p_dsdy ;

}

			//---------------------//
			//	d/dz	       //
			//---------------------//

const Scalar& Scalar::dsdz() const {

    // Protection
    assert(etat != ETATNONDEF) ;

    // If the derivative has not been previously computed, the 
    //  computation must be done by the appropriate routine of the mapping : 

    if (p_dsdz == 0x0) {
      Cmp orig(*this) ;
      Cmp deriv_r(dsdr()) ;
      Cmp deriv_t(srdsdt()) ;
      Cmp deriv_p(srstdsdp()) ;
      Cmp deriv_z(mp) ;
      mp->comp_z_from_spherical(deriv_r, deriv_t, deriv_z) ;
      p_dsdz = new Scalar(deriv_z) ; 
    }
    
    return *p_dsdz ;

}

			//---------------------//
			//	d/dx^i	       //
			//---------------------//

const Scalar& Scalar::deriv(int i) const {
    
    switch (i) {
	
	case 0 : {
	    return dsdx() ; 
	}
	
	case 1 : {
	    return dsdy() ; 
	}
	
	case 2 : {
	    return dsdz() ; 
	}
	
	default : {
	    cout << "Scalar::deriv : index i out of range !" << endl ; 
	    cout << "  i = " << i << endl ; 
	    abort() ; 
	    return dsdx() ;  // Pour satisfaire le compilateur !
	}
	
    }
    
}

			//---------------------//
			//     Laplacian       //
			//---------------------//

const Scalar& Scalar::laplacien(int zec_mult_r) const {

    // Protection
    assert(etat != ETATNONDEF) ;

    // If the Laplacian has not been previously computed, the 
    //  computation must be done by the appropriate routine of the mapping : 
    if ( (p_lap == 0x0) || (zec_mult_r != ind_lap) ) {
	if (p_lap != 0x0) {
	    delete p_lap ;  // the Laplacian had been computed but with
			    //  a different value of zec_mult_r
	}
	Cmp orig(*this) ;
	Cmp laplace(mp) ;
	mp->laplacien(orig, zec_mult_r, laplace) ;
	p_lap = new Scalar(laplace) ;
	ind_lap = zec_mult_r ;
    }
    
    return *p_lap ;
    
}
    
   
    
    
