/*
 *  Resolution of the TT tensor Poisson equation
 *
 *    (see file sym_tensor.h for documentation).
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

char sym_ttt_poisson_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.1  2003/11/07 16:53:52  e_gourgoulhon
 * First version
 *
 *
 * $Header$
 *
 */

// C headers
//#include <>

// Lorene headers
#include "tensor.h"


Sym_tensor_tt Sym_tensor_tt::poisson() const {

	// All this has a meaning only for spherical components:
	assert(dynamic_cast<const Base_vect_spher*>(triad) != 0x0) ; 

    // Solution for the rr-component
	// ----------------------------
	
	Scalar source_rr = operator()(1,1) ; 
	
	// Multiplication by r^2 :
	source_rr.mult_r() ; 	
	source_rr.set_dzpuis( source_rr.get_dzpuis() + 1 ) ;
	source_rr.mult_r_ced() ; 
	source_rr.mult_r() ; 	
	source_rr.set_dzpuis( source_rr.get_dzpuis() + 1 ) ;
	source_rr.mult_r_ced() ; 
	
	Scalar khi = source_rr.poisson() ; 
	
	khi.div_r() ;   		// division 
	khi.div_r() ;   		// by r^2		--> khi now contains h^{rr}
	
	
	// Solution for mu
	// ---------------
	
	Scalar mu_resu = mu().poisson() ;
	
	// Final result
	// ------------
	
	Sym_tensor_tt resu(*mp, *triad, *met_div) ; 

	resu.set_rr_mu(khi, mu_resu) ; 
	
	return resu ;

}
