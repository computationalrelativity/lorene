/*
 * Method of class Star_bin to compute the extrinsic curvature tensor
 *
 */

/*
 *   Copyright (c) 2004 Francois Limousin
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


char star_bin_extr_curv_C[] = "$Header$" ;

/*
 * $Id$
 * $Log$
 * Revision 1.3  2004/02/27 09:52:41  f_limousin
 * Correction of an error on the computation of kcar_auto.
 *
 * Revision 1.2  2004/01/20 15:18:00  f_limousin
 * First version
 *
 *
 * $Header$
 *
 */

// Headers Lorene
#include "star.h"

void Star_bin::extrinsic_curvature(){
    
    // Gradient tilde (with respect to the spherical coordinates
    //           of the mapping)
    // D~_j beta^i 
    
    const Tensor& dshift = shift_auto.derive_con(flat) ; 
    
     // Trace of D~_j beta^i : 
    Scalar div_shift = shift_auto.divergence(flat) ; 

    // Computation of K^{ij}
    // See Eq (49) from Gourgoulhon et al. (2001)
    // ------------------------------------------

    for (int i=1; i<=3; i++) 
	for (int j=1; j<=i; j++) {
	    tkij_auto.set(i, j) = dshift(i, j) + dshift(j, i) - 
		double(2) /double(3) * div_shift * (gtilde.con())(i,j) ; 
	}
    
    
    tkij_auto = 0.5 * tkij_auto / nnn ;   
     
    // Computation of K_{ij} K^{ij}
    // ----------------------------
    
    Tensor tkij_auto_cov = tkij_auto.down(0, gtilde).down(1, gtilde) ;
  
    kcar_comp = contract(tkij_auto_cov, 0, 1, tkij_auto, 0, 1, true) ; 
      
    kcar_auto.std_spectral_base() ; 
    
}
