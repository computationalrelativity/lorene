/*
 * Method of class Et_bin_ncp to compute the extrinsic curvature tensor
 *
 */

/*
 *   Copyright (c) 2003 Francois Limousin
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


char et_bin_ncp_extr_curv_C[] = "$Header$" ;

/*
 * $Header$
 *
 */

// Headers Lorene
#include "et_bin_ncp.h"

void Et_bin_ncp::extrinsic_curvature(){
    
    // Components of shift_auto with respect to the Cartesian triad
    //  (d/dx, d/dy, d/dz) of the mapping : 
    Tenseur shift_auto_local = shift_auto ; 
    shift_auto_local.change_triad( mp.get_bvect_cart() ) ; 
    
    // Gradient tilde (with respect to the Cartesian coordinates
    //           of the mapping)
    // D~_j beta^i 
    
    Tenseur dshift =  shift_auto_local.derive_con(gtilde) ; 
    
    // Return to the absolute reference frame
    dshift.change_triad(ref_triad) ; 
    
    // Trace of D~_j beta^i : 
    Tenseur div_shift = contract(shift_auto_local.derive_cov(gtilde), 0, 1) ; 
    
    // Computation of K^{ij}
    // See Eq (49) from Gourgoulhon et al. (2001)
    // -----------------------------------------
    tkij_auto.set_etat_qcq() ; 
    for (int i=0; i<3; i++) {
	for (int j=i; j<3; j++) {
	    tkij_auto.set(i, j) = dshift(i, j) + dshift(j, i) - double(2) /double(3) * div_shift() * (gtilde.con())(i,j) ; 
	}
    }
    
    tkij_auto =  0.5 * tkij_auto / nnn ; 
    
    tkij_auto.set_std_base() ;

    // Computation of K_{ij} K^{ij}
    // --------------------------------
    
    kcar_auto.set_etat_qcq() ; 
    
    kcar_auto.set() = 0 ; 
    
    kcar_auto.set_std_base() ;

    for (int i=0; i<3; i++) {
	for (int j=0; j<3; j++) {
	
	    kcar_auto.set() +=  contract(operator*(met_gamma.cov(), contract(operator*(met_gamma.cov(), tkij_auto), 1, 3)), 1, 3)() ; 
	
	}
    }
    
    kcar_auto.set_std_base() ; 
    
}
