/*
 *  Resolution of the divergence-free vector Poisson equation
 *
 *    (see file vector.h for documentation).
 *
 */

/*
 *   Copyright (c) 2003 Eric Gourgoulhon & Jerome Novak
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

char vector_df_poisson_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2003/10/27 09:03:01  j_novak
 * Added an assert on the triad
 *
 * Revision 1.2  2003/10/20 19:44:43  e_gourgoulhon
 * First full version
 *
 * Revision 1.1  2003/10/20 14:45:27  e_gourgoulhon
 * Not ready yet...
 *
 *
 * $Header$
 *
 */

// C headers
//#include <>

// Lorene headers
#include "tensor.h"


Vector_divfree Vector_divfree::poisson() const {

  // All this has a meaning only for spherical components:
#ifndef NDEBUG 
  const Base_vect_spher* bvs = dynamic_cast<const Base_vect_spher*>(triad) ;
  assert(bvs != 0x0) ; 
#endif
        // Solution for the r-component
	// ----------------------------
	
	Scalar source_r = *(cmp[0]) ; 
	source_r.mult_r() ; 
	
	Scalar khi = source_r.poisson() ; 
	khi.div_r() ;   // khi now contains V^r
	
	// Solution for mu
	// ---------------
	
	Scalar mu_resu = mu().poisson() ;
	
	// Final result
	// ------------
	
	Vector_divfree resu(*mp, *triad, *met_div) ; 

	resu.set_vr_mu(khi, mu_resu) ; 
	
	return resu ;

}
