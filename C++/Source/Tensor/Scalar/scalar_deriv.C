/*
 * Computations of partial derivatives d/dx, d/dy and d/dz of a Scalar.
 */

/*
 *   Copyright (c) 2003-2004 Eric Gourgoulhon & Jerome Novak
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
 * Revision 1.17  2005/11/24 09:25:08  j_novak
 * Added the Scalar version for the Laplacian
 *
 * Revision 1.16  2005/09/15 15:51:27  j_novak
 * The "rotation" (change of triad) methods take now Scalars as default
 * arguments.
 *
 * Revision 1.15  2005/05/25 16:11:05  j_novak
 * Better handling of the case with no compactified domain.
 *
 * Revision 1.14  2004/08/24 09:14:52  p_grandclement
 * Addition of some new operators, like Poisson in 2d... It now requieres the
 * GSL library to work.
 *
 * Also, the way a variable change is stored by a Param_elliptic is changed and
 * no longer uses Change_var but rather 2 Scalars. The codes using that feature
 * will requiere some modification. (It should concern only the ones about monopoles)
 *
 * Revision 1.13  2004/06/22 08:50:00  p_grandclement
 * Addition of everything needed for using the logarithmic mapping
 *
 * Revision 1.12  2004/02/26 22:51:34  e_gourgoulhon
 * Added methods derive_cov, derive_con and derive_lie.
 *
 * Revision 1.11  2004/01/28 14:02:06  j_novak
 * Suppressed base handling.
 *
 * Revision 1.10  2004/01/28 13:25:42  j_novak
 * The ced_mult_r arguments have been suppressed from the Scalar::*dsd* methods.
 * In the div/mult _r_dzpuis, there is no more default value.
 *
 * Revision 1.9  2004/01/27 15:10:02  j_novak
 * New methods Scalar::div_r_dzpuis(int) and Scalar_mult_r_dzpuis(int)
 * which replace div_r_inc*. Tried to clean the dzpuis handling.
 * WARNING: no testing at this point!!
 *
 * Revision 1.8  2003/11/03 13:37:59  j_novak
 * Still dzpuis...
 *
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
 
// Headers Lorene
#include "scalar.h"
#include "map.h"
#include "tensor.h"
#include "cmp.h"

			//---------------------//
			//         d/dr        //
			//---------------------//

namespace Lorene {
const Scalar& Scalar::dsdr() const {

    // Protection
    assert(etat != ETATNONDEF) ;

    // If the derivative has not been previously computed, the 
    //  computation must be done by the appropriate routine of the mapping : 

    if (p_dsdr == 0x0) {
      p_dsdr = new Scalar(*mp) ;
      if (etat == ETATUN) {
	p_dsdr->set_etat_zero() ;
      }
      else {
	mp->dsdr(*this, *p_dsdr) ;

      }
    }
     
    int dzp = (dzpuis == 0) ? 2 : dzpuis+1 ;
    if (mp->get_mg()->get_type_r(mp->get_mg()->get_nzone() - 1) != UNSURR)
	dzp = 0 ;
    p_dsdr->set_dzpuis(dzp) ;

    return *p_dsdr ;

}

			//--------------------//
			//    1/r d/dtheta    //
			//--------------------//

const Scalar& Scalar::srdsdt() const {

    // Protection
    assert(etat != ETATNONDEF) ;

    // If the derivative has not been previously computed, the 
    //  computation must be done by the appropriate routine of the mapping : 

    if (p_srdsdt == 0x0) {
	p_srdsdt = new Scalar(*mp) ;
	if (etat == ETATUN) {
	  p_srdsdt->set_etat_zero() ;
	}
	else {
	  mp->srdsdt(*this, *p_srdsdt) ;
	}
    }

    int dzp = (dzpuis == 0) ? 2 : dzpuis+1 ;
    if (mp->get_mg()->get_type_r(mp->get_mg()->get_nzone() - 1) != UNSURR)
	dzp = 0 ;
    p_srdsdt->set_dzpuis(dzp) ;

    return *p_srdsdt ;

}


			//------------------------------//
			//    1/(r sin(theta)) d/dphi    //
			//------------------------------//

const Scalar& Scalar::srstdsdp() const {

    // Protection
    assert(etat != ETATNONDEF) ;

    // If the derivative has not been previously computed, the 
    //  computation must be done by the appropriate routine of the mapping : 

    if (p_srstdsdp == 0x0) {
      p_srstdsdp = new Scalar(*mp) ;
      if (etat == ETATUN) {
	p_srstdsdp->set_etat_zero() ;
      }
      else {
	mp->srstdsdp(*this, *p_srstdsdp) ;
      }
    }

    int dzp = (dzpuis == 0) ? 2 : dzpuis+1 ;
    if (mp->get_mg()->get_type_r(mp->get_mg()->get_nzone() - 1) != UNSURR)
	dzp = 0 ;
    p_srstdsdp->set_dzpuis(dzp) ;

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
	
	
	p_dsdt->set_dzpuis(dzpuis) ;

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
	p_stdsdp->set_dzpuis(dzpuis) ;

    return *p_stdsdp ;

}

			//-----------------//
			//      d/dx       //
			//-----------------//

const Scalar& Scalar::dsdx() const {

    // Protection
    assert(etat != ETATNONDEF) ;

    // If the derivative has not been previously computed, the 
    //  computation must be done by the appropriate routine of the mapping : 

    if (p_dsdx == 0x0) {
      p_dsdx = new Scalar(*mp) ;
      if (etat == ETATUN) {
	p_dsdx->set_etat_zero() ;
      }
      else {
	mp->comp_x_from_spherical(dsdr(), srdsdt(), srstdsdp(), *p_dsdx) ;
      }
    }	

    int dzp = (dzpuis == 0) ? 2 : dzpuis+1 ;
    if (mp->get_mg()->get_type_r(mp->get_mg()->get_nzone() - 1) != UNSURR)
	dzp = 0 ;
    p_dsdx->set_dzpuis(dzp) ;

    return *p_dsdx ;

}

			//-----------------//
			//      d/dy       //
			//-----------------//

const Scalar& Scalar::dsdy() const {

    // Protection
    assert(etat != ETATNONDEF) ;

    // If the derivative has not been previously computed, the 
    //  computation must be done by the appropriate routine of the mapping : 

    if (p_dsdy == 0x0) {
      p_dsdy = new Scalar(*mp) ;
      if (etat == ETATUN) {
	p_dsdy->set_etat_zero() ;
      }
      else {
	mp->comp_y_from_spherical(dsdr(), srdsdt(), srstdsdp(), *p_dsdy) ;
      }
    }

    int dzp = (dzpuis == 0) ? 2 : dzpuis+1 ;
    if (mp->get_mg()->get_type_r(mp->get_mg()->get_nzone() - 1) != UNSURR)
	dzp = 0 ;
    p_dsdy->set_dzpuis(dzp) ;

    return *p_dsdy ;

}

			//-----------------//
			//      d/dz       //
			//-----------------//

const Scalar& Scalar::dsdz() const {

    // Protection
    assert(etat != ETATNONDEF) ;

    // If the derivative has not been previously computed, the 
    //  computation must be done by the appropriate routine of the mapping : 

    if (p_dsdz == 0x0) {     
      p_dsdz = new Scalar(*mp) ;
      if (etat == ETATUN) {
	p_dsdz->set_etat_zero() ;
      }
      else {
	mp->comp_z_from_spherical(dsdr(), srdsdt(), *p_dsdz) ;
      }
    }

    int dzp = (dzpuis == 0) ? 2 : dzpuis+1 ;
    if (mp->get_mg()->get_type_r(mp->get_mg()->get_nzone() - 1) != UNSURR)
	dzp = 0 ;
    p_dsdz->set_dzpuis(dzp) ;

    return *p_dsdz ;

}

			//-----------------//
			//      d/dx^i     //
			//-----------------//

const Scalar& Scalar::deriv(int i) const {
    
    switch (i) {
	
	case 1 : {
	    return dsdx() ; 
	}
	
	case 2 : {
	    return dsdy() ; 
	}
	
	case 3 : {
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

                    //--------------------------//
                    //  Covariant derivatives   //
                    //--------------------------//

const Vector& Scalar::derive_cov(const Metric& gam) const {
  
    const Vector* p_resu = 
        dynamic_cast<const Vector*>( &(Tensor::derive_cov(gam)) ) ;

    assert(p_resu != 0x0) ;

    return *p_resu ;
  
}


const Vector& Scalar::derive_con(const Metric& gam) const {
  
    const Vector* p_resu = 
        dynamic_cast<const Vector*>( &(Tensor::derive_con(gam)) ) ;

    assert(p_resu != 0x0) ;

    return *p_resu ;
  
}


                    //--------------------------//
                    //       Lie derivative     //
                    //--------------------------//


Scalar Scalar::derive_lie(const Vector& vv) const {

    Scalar resu(*mp) ; 
    
    compute_derive_lie(vv, resu) ;
    
    return resu ; 
    
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
	p_lap = new Scalar(*mp) ;
	mp->laplacien(*this, zec_mult_r, *p_lap) ;
	ind_lap = zec_mult_r ;
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

    p_lapang->set_dzpuis(dzpuis) ;
    
    return *p_lapang ;
    
}
   
    
    
			//---------------------//
			//         d/dradial   //
			//---------------------//

const Scalar& Scalar::dsdradial() const {

    // Protection
    assert(etat != ETATNONDEF) ;

    // If the derivative has not been previously computed, the 
    //  computation must be done by the appropriate routine of the mapping : 

    if (p_dsdradial == 0x0) {
      p_dsdradial = new Scalar(*mp) ;
      if (etat == ETATUN) {
	p_dsdradial->set_etat_zero() ;
      }
      else {
	mp->dsdradial(*this, *p_dsdradial) ;
      }
    }

    int dzp = (dzpuis == 0) ? 2 : dzpuis+1 ;
    if (mp->get_mg()->get_type_r(mp->get_mg()->get_nzone() - 1) != UNSURR)
	dzp = 0 ;
    p_dsdradial->set_dzpuis(dzp) ;

    return *p_dsdradial ;

}

                         //-----------------//
			//      d/drho   //
			//-----------------//

const Scalar& Scalar::dsdrho() const {

    // Protection
    assert(etat != ETATNONDEF) ;

    // If the derivative has not been previously computed, the 
    //  computation must be done by the appropriate routine of the mapping : 

    if (p_dsdrho == 0x0) {     
      p_dsdrho = new Scalar(*mp) ;
      if (etat == ETATUN) {
	p_dsdrho->set_etat_zero() ;
      }
      else {
	Scalar der_r (dsdr()) ;
	Scalar der_t (srdsdt()) ;
	Valeur val (der_r.get_spectral_va().mult_st() + 
		    der_t.get_spectral_va().mult_ct()) ;

	Scalar res (*mp) ;
	res = val ;
 
	*p_dsdrho = res ;	
      }
    }

    int dzp = (dzpuis == 0) ? 2 : dzpuis+1 ;
    if (mp->get_mg()->get_type_r(mp->get_mg()->get_nzone() - 1) != UNSURR)
	dzp = 0 ;
    p_dsdrho->set_dzpuis(dzp) ;

    return *p_dsdrho ;

}
}
