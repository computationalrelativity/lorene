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


char vector_divfree_aux[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2008/08/27 09:01:27  jl_cornou
 * Methods for solving Dirac systems for divergence free vectors
 *
 * Revision 1.1  2005/02/14 13:01:50  j_novak
 * p_eta and p_mu are members of the class Vector. Most of associated functions
 * have been moved from the class Vector_divfree to the class Vector.
 *
 *
 * $Header$
 *
 */

// C++ headers
#include "headcpp.h"

// C headers
#include <stdlib.h>

// Lorene headers
#include "metric.h"
#include "nbr_spx.h"
#include "utilitaires.h"

// Headers C
#include <stdlib.h>
#include <assert.h>

// Headers Lorene
#include "tensor.h"

			//----------------//
			//  update_etavr  //
			//----------------//

void Vector_divfree::update_etavr() {

	assert(p_A != 0x0) ;

	Scalar eta_tilde(*mp) ;
	Scalar vr(*mp) ;
	sol_Dirac_A(*p_A, eta_tilde, vr, 0x0) ;

	*cmp[0] = vr ;
	p_eta = &eta_tilde ;

	Scalar* p_eta_tmp = p_eta ;  //## in order not to delete p_eta
    	p_eta = 0x0 ;
    	Vector::del_deriv() ;
    
    	p_eta = p_eta_tmp ;
	
}

			//---------------//
			// set_vr_eta_mu //
			//---------------//

 void Vector_divfree::set_vr_eta_mu(const Scalar& vr_i, const Scalar& eta_i,
 			   const Scalar& mu_i) {
     
     // All this has a meaning only for spherical components:
     assert( dynamic_cast<const Base_vect_spher*>(triad) != 0x0 ) ; 
     assert(&vr_i.get_mp() == &eta_i.get_mp()) ; 
 
     // V^r
     *cmp[0] = vr_i ; 
     
     p_eta = new Scalar( eta_i ) ; 	// eta
     
     p_mu = new Scalar( mu_i ) ; 	// mu 
 		
     update_vtvp() ;
 
     return ;
 }


			//--------------//
			//   set_A_mu   //
			//--------------//

void Vector_divfree::set_A_mu(const Scalar& A_i, const  Scalar& mu_i, const Param* par_bc) {

    // All this has a meaning only for spherical components:
    assert( dynamic_cast<const Base_vect_spher*>(triad) != 0x0 ) ; 
    assert(&A_i.get_mp() == &mu_i.get_mp()) ; 

    del_deriv() ; 
    

	p_A = new Scalar (A_i) ;
	p_mu = new Scalar (mu_i) ;
	
	Scalar eta_tilde(*mp) ;
	Scalar vr(*mp) ;
	sol_Dirac_A(*p_A, vr, eta_tilde, par_bc) ;
	*cmp[0] = vr ;
	p_eta = new Scalar(eta_tilde) ;
	p_eta->set_dzpuis(2);
	
//update_etavr();
//set_vr_eta_mu(*cmp[0], *p_eta, *p_mu) ;
	
	update_vtvp();	

	return ;
} 

