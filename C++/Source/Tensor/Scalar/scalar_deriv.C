/*
 * Computations of partial derivatives d/dx, d/dy and d/dz of a Scalar.
 */

/*
 *   Copyright (c) 2003 Eric Gourgoulhon & Jerome Novak
 *
 *   Copyright (c) 1999-2001 Eric Gourgoulhon (for a preceding Cmp version)
 *   Copyright (c) 1999-2001 Philippe Grandclement (for a preceding Cmp version)
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
 * Revision 1.7  2003/10/29 13:14:03  e_gourgoulhon
 * Added integer argument to derivative functions dsdr, etc...
 * so that one can choose the dzpuis of the result (default=2).
 *
 * Revision 1.6  2003/10/17 13:46:15  j_novak
 * The argument is now between 1 and 3 (instead of 0->2)
 *
 * Revision 1.5  2003/10/15 16:03:38  j_novak
 * Added the angular Laplace operator for Scalar.
 *
 * Revision 1.4  2003/10/15 10:43:58  e_gourgoulhon
 * Added new methods dsdt and stdsdp.
 *
 * Revision 1.3  2003/10/11 14:43:29  e_gourgoulhon
 * Changed name of local Cmp "deriv" to "derivee" (in order not
 * to shadow the member deriv).
 *
 * Revision 1.2  2003/10/10 15:57:29  j_novak
 * Added the state one (ETATUN) to the class Scalar
 *
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
			//         d/dr        //
			//---------------------//

const Scalar& Scalar::dsdr(int ced_mult_r) const {

    // Protection
    assert(etat != ETATNONDEF) ;

    // If the derivative has not been previously computed, the 
    //  computation must be done by the appropriate routine of the mapping : 

	if (p_dsdr == 0x0) {
		if (etat == ETATUN) {
	p_dsdr = new Scalar(*mp) ;
	p_dsdr->set_etat_zero() ;
      }
      else {
	Cmp orig(*this) ;
	Cmp derivee(mp) ;
	mp->dsdr(orig, derivee) ;
	p_dsdr = new Scalar(derivee) ; 

		switch (ced_mult_r) {
			case 0 : {
				p_dsdr->dec_dzpuis(2) ; 
				break ;
			}

			case 2 : break ; 
			
			default : {
				cout << "Scalar::dsdr : unexpected value of ced_mult_r !"
					 << endl << "  ced_mult_r = " << ced_mult_r << endl ; 
				abort() ;
				break ; 
			}

		}	  

      }	  
	  
    }
    
    return *p_dsdr ;

}

			//--------------------//
			//    1/r d/dtheta    //
			//--------------------//

const Scalar& Scalar::srdsdt(int ced_mult_r) const {

    // Protection
    assert(etat != ETATNONDEF) ;

    // If the derivative has not been previously computed, the 
    //  computation must be done by the appropriate routine of the mapping : 

    if (p_srdsdt == 0x0) {
      if (etat == ETATUN) {
	p_srdsdt = new Scalar(*mp) ;
	p_srdsdt->set_etat_zero() ;
      }
      else {
	Cmp orig(*this) ;
	Cmp derivee(mp) ;
	mp->srdsdt(orig, derivee) ;
	p_srdsdt = new Scalar(derivee) ;
		switch (ced_mult_r) {
			case 0 : {
				p_srdsdt->dec_dzpuis(2) ; 
				break ;
			}

			case 2 : break ; 
			
			default : {
				cout << "Scalar::srdsdt : unexpected value of ced_mult_r !"
					 << endl << "  ced_mult_r = " << ced_mult_r << endl ; 
				abort() ;
				break ; 
			}

		}	  

      }
	  
    }
    return *p_srdsdt ;

}


			//------------------------------//
			//    1/(r sin(theta) d/dphi    //
			//------------------------------//

const Scalar& Scalar::srstdsdp(int ced_mult_r) const {

    // Protection
    assert(etat != ETATNONDEF) ;

    // If the derivative has not been previously computed, the 
    //  computation must be done by the appropriate routine of the mapping : 

    if (p_srstdsdp == 0x0) {
      if (etat == ETATUN) {
	p_srstdsdp = new Scalar(*mp) ;
	p_srstdsdp->set_etat_zero() ;
      }
      else {
	Cmp orig(*this) ;
	Cmp derivee(mp) ;
	mp->srstdsdp(orig, derivee) ;
	p_srstdsdp = new Scalar(derivee) ; 
		switch (ced_mult_r) {
			case 0 : {
				p_srstdsdp->dec_dzpuis(2) ; 
				break ;
			}

			case 2 : break ; 
			
			default : {
				cout << "Scalar::srstdsdp : unexpected value of ced_mult_r !"
					 << endl << "  ced_mult_r = " << ced_mult_r << endl ; 
				abort() ;
				break ; 
			}

		}	  

      }
	  
    }
    return *p_srstdsdp ;

}

			//--------------------//
			//      d/dtheta      //
			//--------------------//

const Scalar& Scalar::dsdt() const {
    
    assert(etat != ETATNONDEF) ;	// Protection

    // If the derivative has not been previously computed, the 
    //  computation must be done by the appropriate routine of the mapping : 

	if (p_dsdt == 0x0) {
	
		p_dsdt = new Scalar(*mp) ;

		if (etat == ETATUN) {
			p_dsdt->set_etat_zero() ;	
		}
		else {
			mp->dsdt(*this, *p_dsdt) ;
		}
    }
	
    return *p_dsdt ;

}

			//------------------------------//
			//      1/sin(theta) d/dphi     //
			//------------------------------//

const Scalar& Scalar::stdsdp() const {
    
    assert(etat != ETATNONDEF) ;	// Protection

    // If the derivative has not been previously computed, the 
    //  computation must be done by the appropriate routine of the mapping : 

	if (p_stdsdp == 0x0) {
	
		p_stdsdp = new Scalar(*mp) ;

		if (etat == ETATUN) {
			p_stdsdp->set_etat_zero() ;	
		}
		else {
			mp->stdsdp(*this, *p_stdsdp) ;
		}
    }
	
    return *p_stdsdp ;

}

			//-----------------//
			//      d/dx       //
			//-----------------//

const Scalar& Scalar::dsdx(int ced_mult_r) const {

    // Protection
    assert(etat != ETATNONDEF) ;

    // If the derivative has not been previously computed, the 
    //  computation must be done by the appropriate routine of the mapping : 

    if (p_dsdx == 0x0) {
      if (etat == ETATUN) {
	p_dsdx = new Scalar(*mp) ;
	p_dsdx->set_etat_zero() ;
      }
      else {
	Cmp orig(*this) ;
	Cmp deriv_r(dsdr()) ;
	Cmp deriv_t(srdsdt()) ;
	Cmp deriv_p(srstdsdp()) ;
	Cmp deriv_x(mp) ;
	mp->comp_x_from_spherical(deriv_r, deriv_t, deriv_p, deriv_x) ;
	p_dsdx = new Scalar(deriv_x) ; 

		switch (ced_mult_r) {
			case 0 : {
				p_dsdx->dec_dzpuis(2) ; 
				break ;
			}

			case 2 : break ; 
			
			default : {
				cout << "Scalar::dsdx : unexpected value of ced_mult_r !"
					 << endl << "  ced_mult_r = " << ced_mult_r << endl ; 
				abort() ;
				break ; 
			}

		}	  
      }
    }
    
    return *p_dsdx ;

}

			//-----------------//
			//      d/dy       //
			//-----------------//

const Scalar& Scalar::dsdy(int ced_mult_r) const {

    // Protection
    assert(etat != ETATNONDEF) ;

    // If the derivative has not been previously computed, the 
    //  computation must be done by the appropriate routine of the mapping : 

    if (p_dsdy == 0x0) {
      if (etat == ETATUN) {
	p_dsdy = new Scalar(*mp) ;
	p_dsdy->set_etat_zero() ;
      }
      else {
	Cmp orig(*this) ;
	Cmp deriv_r(dsdr()) ;
	Cmp deriv_t(srdsdt()) ;
	Cmp deriv_p(srstdsdp()) ;
	Cmp deriv_y(mp) ;
	mp->comp_y_from_spherical(deriv_r, deriv_t, deriv_p, deriv_y) ;
	p_dsdy = new Scalar(deriv_y) ; 

		switch (ced_mult_r) {
			case 0 : {
				p_dsdy->dec_dzpuis(2) ; 
				break ;
			}

			case 2 : break ; 
			
			default : {
				cout << "Scalar::dsdy : unexpected value of ced_mult_r !"
					 << endl << "  ced_mult_r = " << ced_mult_r << endl ; 
				abort() ;
				break ; 
			}

		}	  

      }
    }
    return *p_dsdy ;

}

			//-----------------//
			//      d/dz       //
			//-----------------//

const Scalar& Scalar::dsdz(int ced_mult_r) const {

    // Protection
    assert(etat != ETATNONDEF) ;

    // If the derivative has not been previously computed, the 
    //  computation must be done by the appropriate routine of the mapping : 

    if (p_dsdz == 0x0) {     
      if (etat == ETATUN) {
	p_dsdz = new Scalar(*mp) ;
	p_dsdz->set_etat_zero() ;
      }
      else {
	Cmp orig(*this) ;
	Cmp deriv_r(dsdr()) ;
	Cmp deriv_t(srdsdt()) ;
	Cmp deriv_p(srstdsdp()) ;
	Cmp deriv_z(mp) ;
	mp->comp_z_from_spherical(deriv_r, deriv_t, deriv_z) ;
	p_dsdz = new Scalar(deriv_z) ; 

		switch (ced_mult_r) {
			case 0 : {
				p_dsdz->dec_dzpuis(2) ; 
				break ;
			}

			case 2 : break ; 
			
			default : {
				cout << "Scalar::dsdz : unexpected value of ced_mult_r !"
					 << endl << "  ced_mult_r = " << ced_mult_r << endl ; 
				abort() ;
				break ; 
			}

		}	  
      }
    }
    return *p_dsdz ;

}

			//-----------------//
			//      d/dx^i     //
			//-----------------//

const Scalar& Scalar::deriv(int i, int ced_mult_r) const {
    
    switch (i) {
	
	case 1 : {
	    return dsdx(ced_mult_r) ; 
	}
	
	case 2 : {
	    return dsdy(ced_mult_r) ; 
	}
	
	case 3 : {
	    return dsdz(ced_mult_r) ; 
	}
	
	default : {
	    cout << "Scalar::deriv : index i out of range !" << endl ; 
	    cout << "  i = " << i << endl ; 
	    abort() ; 
	    return dsdx(ced_mult_r) ;  // Pour satisfaire le compilateur !
	}
	
    }
    
}

			//---------------------//
			//     Laplacian       //
			//---------------------//

const Scalar& Scalar::laplacian(int zec_mult_r) const {

    // Protection
    assert(etat != ETATNONDEF) ;

    // If the Laplacian has not been previously computed, the 
    //  computation must be done by the appropriate routine of the mapping : 
    if ( (p_lap == 0x0) || (zec_mult_r != ind_lap) ) {
	if (p_lap != 0x0) {
	    delete p_lap ;  // the Laplacian had been computed but with
			    //  a different value of zec_mult_r
	}
	if (etat == ETATUN) {
	  p_lap = new Scalar(*mp) ;
	  p_lap->set_etat_zero() ;
	  ind_lap = zec_mult_r ;
	}
	else {
	  Cmp orig(*this) ;
	  Cmp laplace(mp) ;
	  mp->laplacien(orig, zec_mult_r, laplace) ;
	  p_lap = new Scalar(laplace) ;
	  ind_lap = zec_mult_r ;
	}
    }
    
    return *p_lap ;
    
}
    
			//-----------------------------//
			//     Angular Laplacian       //
			//-----------------------------//

const Scalar& Scalar::lapang() const {

    // Protection
    assert(etat != ETATNONDEF) ;

    // If the Laplacian has not been previously computed, the 
    //  computation must be done by the appropriate routine of the mapping : 
    if ( p_lapang == 0x0 ) {
      if (etat == ETATUN) {
	p_lapang = new Scalar(*mp) ;
	p_lapang->set_etat_zero() ;
      }
      else {
	p_lapang = new Scalar(*mp) ;
	mp->lapang(*this, *p_lapang) ;
      }
    }
    
    return *p_lapang ;
    
}
   
    
    
