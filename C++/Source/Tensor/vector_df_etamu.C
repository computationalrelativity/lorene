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
		
		// treatment of the CED : 
		int nzm1 = mp->get_mg()->get_nzone() - 1 ; // index of the CED
		Scalar dvr_ext(*mp) ;
		dvr_ext.allocate_all() ; 
		dvr_ext.annule(0,nzm1-1) ; // zero in all domains but the CED
		dvr_ext.set_domain(nzm1) = dvr.domain(nzm1) ; // equal to dvr in the CED
		dvr_ext.set_spectral_base( (dvr.get_spectral_va()).get_base() ) ; 
		dvr_ext.div_r() ; // division by r in the CED
		
		// Multiplication by r in all domains but the CED:
		dvr.annule_domain(nzm1) ; 
		dvr.mult_r() ; 
		
		// Adding the CED part
		dvr.set_domain(nzm1) = dvr_ext.domain(nzm1) ; 

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

		abort() ; 

	}

	return *p_mu ; 

}











