/*
 *  Methods of class Vector_divfree related to eta and mu
 *
 *   (see file vector.h for documentation)
 *
 */

/*
 *   Copyright (c) 2003 Eric Gourgoulhon & Jerome Novak
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


char vector_df_etamu_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.6  2004/12/28 10:37:24  j_novak
 * Better way of enforcing zero divergence.
 *
 * Revision 1.5  2004/03/29 07:13:39  e_gourgoulhon
 * Treatment of dvr.etat = ETATZERO in eta().
 *
 * Revision 1.4  2003/10/20 10:13:28  e_gourgoulhon
 * Corrected error in the computation of mu.
 *
 * Revision 1.3  2003/10/17 16:36:05  e_gourgoulhon
 * In method update_vtvp(): replaced del_deriv() by
 *   Vector::del_deriv().
 *
 * Revision 1.2  2003/10/16 21:38:02  e_gourgoulhon
 * Added methods mu() and update_vtvp().
 *
 * Revision 1.1  2003/10/16 15:28:25  e_gourgoulhon
 * First version
 *
 *
 * $Header$
 *
 */

// Headers C
#include <stdlib.h>
#include <assert.h>

// Headers Lorene
#include "tensor.h"

			//--------------//
			//     eta      //
			//--------------//
			
			
const Scalar& Vector_divfree::eta() const {


	if (p_eta == 0x0) {   // a new computation is necessary
		
		// All this has a meaning only for spherical components:
		#ifndef NDEBUG 
		const Base_vect_spher* bvs = dynamic_cast<const Base_vect_spher*>(triad) ;
		assert(bvs != 0x0) ; 
		#endif

		// eta is computed from the divergence-free condition:

		Scalar dvr = - cmp[0]->dsdr() ; 	// - dV^r/dr  ( -r^2 dV^r/dr in the CED)
		dvr.mult_r_dzpuis(0) ;

		// Final result for the V^r source for eta:
		dvr -= 2. * (*cmp[0]) ; 
		
		// Resolution of the angular Poisson equation for eta
		// --------------------------------------------------
		p_eta = new Scalar( dvr.poisson_angu() ) ; 
	
	}

	return *p_eta ; 

}


			//--------------//
			//     mu       //
			//--------------//
			
			
const Scalar& Vector_divfree::mu() const {


	if (p_mu == 0x0) {   // a new computation is necessary
		
		// All this has a meaning only for spherical components:
		#ifndef NDEBUG 
		const Base_vect_spher* bvs = dynamic_cast<const Base_vect_spher*>(triad) ;
		assert(bvs != 0x0) ; 
		#endif

		Scalar tmp = *cmp[2] ; 	// V^ph
		tmp.div_tant() ; 		// V^ph / tan(th)
		
		// dV^ph/dth + V^ph/tan(th) - 1/sin(th) dV^th/dphi 
		tmp = cmp[2]->dsdt() + tmp - cmp[1]->stdsdp() ; 
		
		// Resolution of the angular Poisson equation for mu
		// --------------------------------------------------
		p_mu = new Scalar( tmp.poisson_angu() ) ;  

	}

	return *p_mu ; 

}


			//----------------//
			//  update_vtvp   //
			//----------------//
			

void Vector_divfree::update_vtvp() {

	assert( (p_eta != 0x0) && (p_mu != 0x0) ) ; 
	
	// V^theta :
	*cmp[1] = p_eta->dsdt() - p_mu->stdsdp() ; 
	
	// V^phi : 
	*cmp[2] = p_eta->stdsdp() + p_mu->dsdt() ; 
	
	Vector::del_deriv() ; //## in order not to delete p_eta and p_mu
	
}			









