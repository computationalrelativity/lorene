/*
 *  Methods of class Vector related to eta and mu
 *
 *   (see file vector.h for documentation)
 *
 */

/*
 *   Copyright (c) 2005 Eric Gourgoulhon & Jerome Novak
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


 

/*
 * $Id$
 * $Log$
 * Revision 1.6  2024/01/26 17:44:25  g_servignat
 * Updated the Pseudopolytrope_1D class to be consistent with the paper (i.e. with a GPP in the middle)
 *
 * Revision 1.5  2016/12/05 16:18:18  j_novak
 * Suppression of some global variables (file names, loch, ...) to prevent redefinitions
 *
 * Revision 1.4  2014/10/13 08:53:45  j_novak
 * Lorene classes and functions now belong to the namespace Lorene.
 *
 * Revision 1.3  2014/10/06 15:13:21  j_novak
 * Modified #include directives to use c++ syntax.
 *
 * Revision 1.2  2008/08/27 08:52:23  jl_cornou
 * Added fonctions for angular potential A
 *
 * Revision 1.1  2005/02/14 13:01:50  j_novak
 * p_eta and p_mu are members of the class Vector. Most of associated functions
 * have been moved from the class Vector_divfree to the class Vector.
 *
 *
 * $Header$
 *
 */

// Headers C
#include <cstdlib>
#include <cassert>

// Headers Lorene
#include "tensor.h"

			//--------------//
			//     eta      //
			//--------------//
			
			
namespace Lorene {
const Scalar& Vector::eta() const {


    if (p_eta == 0x0) {   // a new computation is necessary
	
	// All this has a meaning only for spherical components:
#ifndef NDEBUG 
	const Base_vect_spher* bvs = dynamic_cast<const Base_vect_spher*>(triad) ;
	assert(bvs != 0x0) ; 
#endif
	
	// eta is computed from its definition:
	Scalar sou_eta = *cmp[1] ;  //V^th
	sou_eta.div_tant() ;
	sou_eta += cmp[1]->dsdt() + cmp[2]->stdsdp();
		// cout << "sou_eta: " << *sou_eta.get_spectral_va().c << endl;
		// sou_eta.set_spectral_va().ylm() ;
		// sou_eta.spectral_display(" ", 1e-10) ;
		// abort() ;
	// Resolution of the angular Poisson equation for eta
	// --------------------------------------------------
	p_eta = new Scalar( sou_eta.poisson_angu() ) ; 
	
    }

    return *p_eta ; 

}


			//--------------//
			//     mu       //
			//--------------//
			
			
const Scalar& Vector::mu() const {


    if (p_mu == 0x0) {   // a new computation is necessary
	
	// All this has a meaning only for spherical components:
#ifndef NDEBUG 
	const Base_vect_spher* bvs = dynamic_cast<const Base_vect_spher*>(triad) ;
	assert(bvs != 0x0) ; 
#endif
	
	Scalar tmp = *cmp[2] ; 	// V^ph
	tmp.div_tant() ; 		// V^ph / tan(th)
	
	// dV^ph/dth + V^ph/tan(th) - 1/sin(th) dV^th/dphi 
	tmp += cmp[2]->dsdt() - cmp[1]->stdsdp() ; 
	
	// Resolution of the angular Poisson equation for mu
	// --------------------------------------------------
	p_mu = new Scalar( tmp.poisson_angu() ) ;  
	
    }
    
    return *p_mu ; 

}

			//-----------//
			//    A      //
			//-----------//

const Scalar& Vector::A() const {
	
	
	if (p_A == 0x0) {   // A new computation is necessary
	
		// All this has a meaning only for spherical components :
#ifndef NDEBUG 
	const Base_vect_spher* bvs = dynamic_cast<const Base_vect_spher*>(triad) ;
	assert(bvs != 0x0) ; 
#endif
	
	// p_eta doit être calculé
	if (p_eta == 0x0) { Scalar etatmp = this->eta(); }

	Scalar tmp = -*cmp[0] ;   // -V^r
	tmp.div_r_dzpuis(2);		 // -V^r/r

	Scalar eta_tilde = *p_eta ;
	Scalar etad = eta_tilde.dsdr() ; 
	eta_tilde.div_r_dzpuis(2);
	etad.set_dzpuis(2);
	tmp += etad + eta_tilde ; // d eta / dr + eta/r

	p_A = new Scalar (tmp) ;
	}
	
	return *p_A ;
}
	




			//----------------//
			//  update_vtvp   //
			//----------------//
			

void Vector::update_vtvp() {

    assert( (p_eta != 0x0) && (p_mu != 0x0) ) ; 

    // V^theta :
    *cmp[1] = p_eta->dsdt() - p_mu->stdsdp() ; 

    // V^phi : 
    *cmp[2] = p_eta->stdsdp() + p_mu->dsdt() ; 
    
    Scalar* p_eta_tmp = p_eta ;  //## in order not to delete p_eta and p_mu
    p_eta = 0x0 ;
    Scalar* p_mu_tmp = p_mu ;
    p_mu = 0x0 ;	
    Vector::del_deriv() ;
    
    p_eta = p_eta_tmp ;
    p_mu = p_mu_tmp ;
    
}			


void Vector::set_vr_eta_mu(const Scalar& vr_i, const Scalar& eta_i,
			   const Scalar& mu_i) {
    
    // All this has a meaning only for spherical components:
    assert( dynamic_cast<const Base_vect_spher*>(triad) != 0x0 ) ; 
    assert(&vr_i.get_mp() == &eta_i.get_mp()) ; 

    del_deriv() ; 
    
    // V^r
    *cmp[0] = vr_i ; 

    p_eta = new Scalar( eta_i ) ; 	// eta
    
    p_mu = new Scalar( mu_i ) ; 	// mu 
		
    update_vtvp() ;

    return ;
}





}
